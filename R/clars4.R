# Try to add tree methodology
clars4 <- function(x, y, cost, maxk = 50, eps = 1e-6, trace = FALSE, costfunc = NULL)
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


    Gram <- t(x) %*% x

    # number of steps
    k <- 0
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

    activeMatrix <- matrix(FALSE, nrow = maxk, ncol = p)
    skipMatrix <- matrix(FALSE, nrow = maxk, ncol = p)
    tree <- NULL

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
        trace.out$cvecP <- 0
        trace.out$score <- cvec /cost
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
        k <- k + 1

        # 1. We have a new variable entering. Update residuals
        Inactive <- !Active
        nv <- nv + 1
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
            legal <- gamvec < 2 * gmax

            # Are there any variables that meet conditions?
            if (length(tree))
            {
                if (tree[1] == 0)
                {
                    OK <- Inactive & legal &! skip
                } else {
                    OK <- Inactive & legal &! skip &! skipMatrix[tree[1], ]
                }
            } else {
                OK <- Inactive & legal &! skip
            }
            if (any(OK))
            {
                tree <- c(k, tree)
                cvecP <- cmax - gamvec * AA
                score <- abs((cvecP) / price)
                best <- max(abs(score[OK]))
                newj <- which(best == score)
                gamma <- gamvec[newj]

            } else {
                # Need to revert along the tree branch
                # TODO: left off here
                # 1. Make a fake step here
                beta[k, ] <- 0
                # 2. Reset conditions to previous
                betaC <- beta[tree[1], ]
                muC <- mu[tree[1], ]
                # 3. Record variable we must skip
                skipper <- which(activeMatrix[tree[2], ] &! activeMatrix[tree[3], ])

                skipMatrix[tree[1], skipper] <- TRUE

                Active <- activeMatrix[tree[1], ]

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
                    cat("\n")
                    cat("None available!!!")
                    cat("\n")
                    cat("Remove variable: ")
                    cat(colnames(x)[skipper])
                    cat("\n")
                    cat("cmax / A: ")
                    cat(cmax / AA)
                    cat("\n")
                    cat("Current SS: ")
                    cat(var(e))
                    cat("\n")

                    trace.out$In <- ifelse(Inactive, "Out", "In")
                    trace.out$In[!Active & j] <- "next"
                    trace.out$In[skip] <- "skip"
                    # trace.out$cont <- ifelse(contenders, "Y", "N")
                    trace.out$cvecP <- cvecP
                    trace.out$score <- score
                    trace.out$gammaP <- gamP
                    trace.out$gammaN <- gamN
                    trace.out$gamtilde <- gammaj
                    ww <- rep(0, p)
                    ww[Active] <- w
                    trace.out$w <- ww
                    trace.out$beta <- beta[k + 1, ]
                    print(trace.out)
                    cat("\n\n")
                }
                cat("gmax: ")
                cat(gmax)
                cat("\n")
                cat("gamvec: ")
                cat(gamvec)
                cat("\n")
                cat(colnames(x)[skipper])
                cat(" was a bad idea")
                cat("\n")
                cat("tree: ")
                cat(tree)
                cat("\n")
                cat("Active:" )
                cat(Active)
                cat("\n")
                cat("actmat:" )
                cat(activeMatrix[tree[1], ])
                cat("\n")
                cat("Base off of ")
                cat(tree[1])
                cat("\n")
                cat("skippers are: ")
                cat(skipMatrix[tree[1], ])
                cat("\n")
                cat("updated beta is: ")
                cat(betaC)
                cat("\n")


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
                gamma <- cmax / AA
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
        if (drop(gamma) > cmax / AA)
        {
            muC <- muC + drop(gamma) * u
            betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
            # short gamma step
            shortgamma <- cmax / AA
            mu[k + 1, ] <- mu[k, ] + drop(shortgamma) * u
            beta[k + 1, Active] <- betaC[Active] + drop(shortgamma) * w * Signs
            # whole gamma step
            mu[k + 2, ] <- muC
            beta[k + 2, Active] <- betaC[Active] + drop(gamma) * w * Signs
            activeMatrix[k + 1, ] <- Active

            k <- k + 1
        } else {
            muC <- muC + drop(gamma) * u
            betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
            mu[k + 1, ] <- muC
            beta[k + 1, Active] <- betaC[Active] + drop(gamma) * w * Signs
        }


        if (any(skip))
        {
            if (trace)
            {
                if (beta[k + 1, skip] != 0)
                {
                    cat("Should have been zero.  Was: ")
                    cat(beta[k + 1, skip])
                }
            }
            beta[k + 1, skip] <- 0
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
            cat(cmax / AA)
            cat("\n")
            cat("Current SS: ")
            cat(var(e))
            cat("\n")

            trace.out$In <- ifelse(Inactive, "Out", "In")
            trace.out$In[!Active & j] <- "next"
            trace.out$In[skip] <- "skip"
            # trace.out$cont <- ifelse(contenders, "Y", "N")
            trace.out$cvecP <- cvecP
            trace.out$score <- score
            trace.out$gammaP <- gamP
            trace.out$gammaN <- gamN
            trace.out$gamtilde <- gammaj
            ww <- rep(0, p)
            ww[Active] <- w
            trace.out$w <- ww
            trace.out$beta <- beta[k + 1, ]
            print(trace.out)
            cat("\n\n")
        }
        Active <- j
        activeMatrix[k, ] <- Active
    }

    list(beta = beta[1:(k + 1), ])
}
