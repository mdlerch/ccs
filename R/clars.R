# Try to add tree methodology
clarsO <- function(x, y, cost, maxk = 50, eps = 1e-6, trace = FALSE, costfunc = NULL)
{
    # default costfunc is just sum of used variables
    if (is.null(costfunc))
    {
        costfunc <- function(Active)
        {
            sum(cost[Active])
        }
    }
    x <- scale(x)
    # variable setup
    n <- nrow(x); p <- ncol(x)

    if (missing(cost)) cost <- rep(1, p)

    # sets
    Inactive <- rep(TRUE, p)
    Active <- rep(FALSE, p)

    activeMatrix <- matrix(FALSE, nrow = maxk + 1, ncol = p)
    skipMatrix <- matrix(FALSE, nrow = maxk + 1, ncol = p)
    treeMatrix <- matrix(0, nrow = maxk + 1, ncol = maxk + 1)
    treenum <- 1
    tree <- NULL
    maxdepth <- 1

    Gram <- t(x) %*% x

    # number of steps
    k <- 1
    # number of variables currently in model
    nv <- 0

    mu <- matrix(0, nrow = maxk + 1, ncol = n)
    muC <- rep(0, n)
    beta <- matrix(0, nrow = maxk + 1, ncol = p)
    betaC <- rep(0, p)

    # Equation 2.8
    e <- y - muC
    cvec <- t(x) %*% e
    cmax <- max(abs(cvec))
    svec <- abs(cvec) / cost
    smax <- max(svec)
    j <- svec >= smax - eps
    Active <- j
    newj <- which(j)

    tree <- 1

    gmax <- 0

    trace.out <- data.frame(var = colnames(x), cost = cost)

    if (trace)
    {
        cat("\nIteration: ")
        cat(k)
        cat("\nCurrent number of variables: 0")
        cat("\n")
        cat("\nTry to put in next ")
        cat(as.character(trace.out$var[j]))
        cat("\n")
        trace.out$In <- ifelse(Active, "In", "Out")
        trace.out$In[!Active & j] <- "next"
        # trace.out$cont <- 0
        trace.out$cvec <- cvec
        trace.out$score <- cvec / cost
        trace.out$gammaP <- 0
        trace.out$gammaN <- 0
        trace.out$gamtilde <- 0
        trace.out$beta <- 0
        print(trace.out)
        cat("\n\n")
    }

    skip <- rep(FALSE, p)

    while (nv < p & k < maxk)
    {

        cat("tree: ")
        cat(tree)
        cat("\n")
        k <- k + 1
        activeMatrix[k, ] <- Active

        # 1. We have a new variable entering. Update residuals
        Inactive <- !Active
        nv <- sum(Active)
        e <- y - muC
        cvec <- t(x) %*% e
        cmax <- max(abs(cvec[Active]))

        # 2. Find unit-vector of equal projection.
        Signs <- sign(cvec[Active])
        XA <- x[ , Active] * rep(1, n) %*% t(Signs)
        gA <- t(XA) %*% XA
        one <- rep(1, sum(Active))
        AA <- 1 / sqrt(one %*% solve(gA) %*% one)
        w <- AA %*% t(solve(gA) %*% one)
        u <- XA %*% t(w)

        # 3. Increment model fit in the direction of u.

        a <- t(x) %*% u

        ## to OLS solution
        gmax <- drop(cmax / AA)

        # cat("\n")
        # cat("Iteration: ")
        # cat(k - 1)
        # cat("\n")
        # cat("Tree: ")
        # cat(tree)
        # cat("\n")
        # cat("activeMatrix: ")
        # cat(activeMatrix[tree[1], ])
        # cat("\n")
        # cat("activeMatrix: ")
        # cat(activeMatrix[tree[2], ])
        # cat("\n")
        # cat("beta: ")
        # cat(betaC)
        # cat("\n")
        # cat("gmax: ")
        # cat(gmax)
        # cat("\n")
        # cat("Skips: ")
        # cat(skipMatrix[tree[1], ])
        # cat("\n")


        if (nv == p)
        {
            gamma <- gmax
        } else {

            # calculate each variables next price
            # TODO: can this be made more efficient?
            price <- rep(0, p)
            for (i in 1:p)
            {
                ActiveN <- Active
                ActiveN[i] <- TRUE
                price[i] <- costfunc(ActiveN)
            }
            # Two gamma calculations
            gamP <- rep(0, p)
            gamN <- rep(0, p)
            gamP[Inactive] <- (cmax - cvec[Inactive]) / (AA - a[Inactive])
            gamN[Inactive] <- (cmax + cvec[Inactive]) / (AA + a[Inactive])

            # First choice, gamma between 0 and cmax/A
            gamvec <- apply(cbind(gamP, gamN), 1, mingt0)
            # TODO: how big can gamma be? can it be negative
            OK <- Inactive & (gamvec < 2 * gmax) &! skipMatrix[tree[1], ]

            if (any(OK))
            {
                cvecP <- cmax - gamvec * AA
                score <- abs((cvecP) / price)
                # TODO: this includes all "acceptable"
                best <- max(abs(score[OK]))
                newj <- which(best == score)
                gamma <- gamvec[newj]
                # actual updates later in code
                tree <- c(k, tree)
            } else {
                # complete the ols
                gamma <- gmax
                muC <- muC + drop(gamma) * u
                betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
                mu[k, ] <- muC
                beta[k, ] <- betaC
                Active[newj] <- TRUE
                activeMatrix[k, ] <- Active

                treeMatrix[ , treenum] <- c(tree, rep(0, maxk + 1 - length(tree)))
                treenum <- treenum + 1
                maxdepth <- max(maxdepth, length(tree))

                # revert
                if (length(tree) == 1)
                {
                    # TODO: might need to set up a skipMatrix for tree == 0
                    skipper <- which(activeMatrix[tree[1], ])
                    betaC <- rep(0, p)
                    muC <- rep(0, n)
                    Active <- rep(FALSE, p)
                } else {
                    toskip <- which(activeMatrix[tree[1], ] &! activeMatrix[tree[2], ])
                    skipMatrix[tree[2], toskip] <- TRUE
                    betaC <- beta[tree[2], ]
                    muC <- mu[tree[2], ]
                    Active <- activeMatrix[tree[2], ]
                    tree <- tree[-1]
                    newj <- NA
                }

                # TODO: do i need this?
                next
            }

            direction <- sign(gamma)

            j[newj] <- TRUE
        }

        skip <- rep(FALSE, p)
        # find gamma tilde for each j
        gammaj <- rep(0, p)
        # TODO: should this be betaC?
        gammaj[Active] <- -beta[k, Active] / (w * Signs)
        # If there are any gammaj that will eventually cross
        flag.cross <- FALSE
        # TODO: gamma.tilde could also be negative if we choose a negative
        # gamma...
        if (any(direction * gammaj > 0))
        {
            if (nv == p)
            {
                gamma <- gmax
                temp <- gammaj
                temp[temp <= 0] <- max(temp) + 1
                outj <- which.min(temp)
                gamma.tilde <- gammaj[outj]
                if (gamma.tilde < gamma)
                {
                    gamma <- gamma.tilde
                    j[outj] <- FALSE
                    nv <- nv - 2
                }
            }
            # get the first to cross
            temp <- gammaj
            temp[temp <= 0] <- max(temp) + 1
            outj <- which.min(temp)
            gamma.tilde <- gammaj[outj]
            if (gamma.tilde < gamma)
            {
                flag.cross <- TRUE
                gamma <- gamma.tilde
                j[outj] <- FALSE
                j[newj] <- FALSE
                nv <- nv - 2
                skip[outj] <- TRUE
            }
        }

        # If gamma is longer than cmax/AA make a pit stop along the way
        if (drop(gamma) > gmax)
        {
            # short gamma step
            shortgamma <- gmax
            mu[k, ] <- muC + drop(shortgamma) * u
            beta[k, Active] <- betaC[Active] + drop(shortgamma) * w * Signs
            # whole gamma step
            muC <- muC + drop(gamma) * u
            betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
            mu[k + 1, ] <- muC
            beta[k + 1, Active] <- betaC[Active]
            # TODO: do i want this?
            activeMatrix[k + 1, ] <- Active
            k <- k + 1
        } else {
            muC <- muC + drop(gamma) * u
            betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
            mu[k, ] <- muC
            beta[k, Active] <- betaC[Active]
        }


        if (any(skip))
        {
            if (trace)
            {
                if (beta[k, skip] != 0)
                {
                    cat("Should have been zero.  Was: ")
                    cat(beta[k, skip])
                }
            }
            beta[k, skip] <- 0
        }
        if (trace)
        {
            cat("\nIteration: ")
            cat(k)
            cat("\nCurrent number of variables: ")
            cat(nv)
            # if (flag.contender)
            # {
            #     cat("\nNext in via contender: ")
            #     cat(as.character(trace.out$var[newj]))
            # } else {
                cat("\nNext in via lars: ")
                cat(as.character(trace.out$var[newj]))
            # }
            if (flag.cross)
            {
                cat("\nBut variable ")
                cat(as.character(trace.out$var[outj]))
                cat(" is crossing zero")
            }
            cat("\n")
            cat("Selected gamma: ")
            cat(drop(gamma))
            cat("\n")
            cat("cmax / A: ")
            cat(gmax)
            cat("\n")
            cat("Current SS: ")
            cat(var(e))
            cat("\n")

            trace.out$In <- ifelse(Inactive, "Out", "In")
            trace.out$In[!Active & j] <- "next"
            trace.out$In[skip] <- "skip"
            # trace.out$cont <- ifelse(contenders, "Y", "N")
            trace.out$cvec <- cvec
            trace.out$score <- score
            trace.out$gammaP <- gamP
            trace.out$gammaN <- gamN
            trace.out$gamtilde <- gammaj
            ww <- rep(0, p)
            ww[Active] <- w
            trace.out$w <- ww
            trace.out$beta <- beta[k, ]
            print(trace.out)
            cat("\n\n")
        }
        Active <- j
        activeMatrix[k, ] <- Active
    }

    treeMatrix[ , treenum] <- c(tree, rep(0, maxk + 1 - length(tree)))
    maxdepth <- max(maxdepth, length(tree))

    list(mu = mu[1:k, ], beta = beta[1:k, ], tree = treeMatrix[1:maxdepth, 1:treenum])
}
