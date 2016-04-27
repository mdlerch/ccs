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

    muC <- rep(0, n)
    beta <- matrix(0, nrow = maxk + 1, ncol = p)
    mu <- matrix(0, nrow = maxk + 1, ncol = n)

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
            legal <- gamvec < gmax
            nonzero <- gamvec != 0

            # Are there any variables that meet conditions?
            OK <- Inactive & legal & nonzero &! skip
            if (any(OK))
            {
                tree <- c(k - 1, tree)
                base <- k + 1
                cvecP <- cmax - gamvec * AA
                score <- abs((cvecP) / price)
                best <- max(abs(score[OK]))
                newj <- which(best == score)
                gamma <- gamvec[newj]

            } else {
                # Need to revert along the tree branch
                # TODO: left off here
                skipper <- which(activeMatrix[k, ] &! activeMatrix[tree[1], ])
                base <- base - 1
                cat(colnames(x)[skipper])
                cat(" was a bad idea")
                cat("\n")
                cat("Base off of ")
                cat(base - 1)
                cat("\n")
                cat("\n")
                print(beta[1:(k + 1), ])
                print(beta[base, ])
                cat("\n")
                gamma <- gmax
                muC <- mu[base, ]
                Active <- activeMatrix[base, ]
                skip <- rep(FALSE, p)
                skip[skipper] <- TRUE
                next
            }

            direction <- sign(gamma)

            j[newj] <- TRUE
        }

        skip <- rep(FALSE, p)
        # find gamma tilde for each j
        gammaj <- rep(0, p)
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
            # short gamma step
            shortgamma <- cmax / AA
            mu[k + 1, ] <- mu[k, ] + drop(shortgamma) * u
            beta[k + 1, Active] <- beta[k, Active] + drop(shortgamma) * w * Signs
            # whole gamma step
            mu[k + 2, ] <- muC
            beta[k + 2, Active] <- beta[k, Active] + drop(gamma) * w * Signs
            activeMatrix[k + 1, ] <- Active

            k <- k + 1
        } else {
            muC <- muC + drop(gamma) * u
            mu[k + 1, ] <- muC
            beta[k + 1, Active] <- beta[k, Active] + drop(gamma) * w * Signs
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
        activeMatrix[k + 1, ] <- Active
    }

    list(beta = beta[1:(k + 1), ])
}
