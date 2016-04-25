clars <- function(x, y, cost, maxk = 50, eps = 1e-6, trace = FALSE, costfunc = NULL)
{
    # default costfunc is just sum of used variables
    if (is.null(costfunc))
    {
        costfunc <- function(Active)
        {
            sum(cost[Active])
        }
    }
    lastgam <- 10000
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

    mu <- rep(0, n)
    beta <- matrix(0, nrow = maxk + 1, ncol = p)
    mul <- matrix(0, nrow = maxk + 1, ncol = n)

    # Equation 2.8
    e <- y - mu
    cvec <- t(x) %*% e
    cmax <- max(abs(cvec))
    svec <- abs(cvec) / cost
    smax <- max(svec)
    j <- svec >= smax - eps
    Active <- j

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

        # 1. Find the next variable to add.
        Inactive <- !Active
        nv <- nv + 1
        e <- y - mu
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

        if (nv == p)
        {
            gamma <- cmax / AA
        } else {
            gmax <- cmax / AA
            price1 <- sum(cost[Active])
            # TODO: this can be made more efficient
            price <- rep(0, p)
            for (i in 1:p)
            {
                ActiveN <- Active
                ActiveN[i] <- TRUE
                price[i] <- costfunc(ActiveN)
            }
            gamP <- rep(0, p)
            gamN <- rep(0, p)
            gamP[Inactive] <- (cmax - cvec[Inactive]) / (AA - a[Inactive])
            gamN[Inactive] <- (cmax + cvec[Inactive]) / (AA + a[Inactive])

            # First choice, gamma between 0 and cmax/A
            gamvec <- apply(cbind(gamP, gamN), 1, mingt0)
            legal <- gamvec < drop(cmax / AA)


            # if first choice
            if (any(Inactive &! skip & legal))
            {
                cvecP <- cmax - gamvec * AA
                score <- abs((cvecP) / price)
                best <- max(abs(score[Inactive &! skip & legal]))
                newj <- which(best == score)
                gamma <- gamvec[newj]
                lastgam <- gamma
            } else {
                # second choice |gammas| smaller than cmax / A
                # chosen by min c * price
                gamvecN <- apply(cbind(gamP, gamN), 1, maxlt0)
                legal <- abs(gamvec) < drop(cmax / AA) & gamvec != 0
                if (any(Inactive &! skip & legal))
                {
                    cvecP <- cmax - gamvec * AA
                    score <- abs((cvecP) * price)
                    best <- min(abs(score[Inactive &! skip & legal]))
                    newj <- which(best == score)
                    gamma <- gamvec[newj]
                } else {
                    # third choice any gammas
                    # chosen by min c * price
                    gamvec <- apply(cbind(gamP, gamN), 1, minabs)
                    cvecP <- cmax - gamvec * AA
                    score <- abs((cvecP) * price)
                    best <- min(abs(score[Inactive &! skip]))
                    newj <- which(best == score)
                    gamma <- gamvec[newj]
                }
            }

            cvecP <- cmax - gamvec * AA
            score <- abs((cvecP) / price)
            best <- max(abs(score[Inactive &! skip]))
            newj <- which(best == score)
            gamma <- gamvec[newj]

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
            mu <- mu + drop(gamma) * u
            # short gamma step
            shortgamma <- cmax / AA
            mul[k + 1, ] <- mul[k, ] + drop(shortgamma) * u
            beta[k + 1, Active] <- beta[k, Active] + drop(shortgamma) * w * Signs
            # whole gamma step
            mul[k + 2, ] <- mu
            beta[k + 2, Active] <- beta[k, Active] + drop(gamma) * w * Signs
            k <- k + 1
        } else {
            mu <- mu + drop(gamma) * u
            mul[k + 1, ] <- mu
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
    }

    list(beta = beta[1:(k + 1), ])
}
