clars <- function(x, y, cost, maxk = 50, eps = 1e-6, trace = FALSE)
{
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
        trace.out$cvec <- cvec
        trace.out$score <- cvec / cost
        trace.out$scoreP <- 0
        trace.out$scoreN <- 0
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
            ### METHOD 1
            # opts <- Inactive & !skip
            # gamP <- rep(0, p)
            # gamN <- rep(0, p)
            # gamP[opts] <- (cmax - cvec[opts]) / (AA - a[opts])
            # gamN[opts] <- (cmax + cvec[opts]) / (AA + a[opts])

            # gammas <- apply(cbind(gamP, gamN), 1, mindist)


            # idx <- which.max((cvec[opts] - gammas[opts] * a[opts]) / cost[opts])
            # newj <- which(cumsum(opts) == idx)[1]

            # gamma <- gammas[newj]
            # direction <- sign(gamma)
            ### METHOD 1

            #### METHOD 2
            # gamP <- rep(0, p)
            # gamN <- rep(0, p)
            # gamP[Inactive] <- (cmax - cvec[Inactive]) / (AA - a[Inactive])
            # gamN[Inactive] <- (cmax + cvec[Inactive]) / (AA + a[Inactive])

            # rP <- apply(cbind(gamP), 1, function(g) e - g * u)
            # rN <- apply(cbind(gamN), 1, function(g) e - g * u)

            # scoreP <- apply(rP, 2, var) * cost
            # scoreN <- apply(rN, 2, var) * cost

            # bestP <- min(abs(scoreP[Inactive &! skip]))
            # bestN <- min(abs(scoreN[Inactive &! skip]))

            # if (bestP < bestN)
            # {
            #     newj <- which(bestP == abs(scoreP))
            #     gamma <- gamP[newj]
            # } else {
            #     newj <- which(bestN == abs(scoreN))
            #     gamma <- gamN[newj]
            # }
            # direction <- sign(gamma)
            ### METHOD 2

            ### METHOD 3
            # scoreP <- cvec / cost
            # scoreN <- cvec / cost

            # best <- max(cvec[Inactive & !skip] / cost[Inactive & !skip])
            # newj <- which(cvec / cost == best)

            # gamP <- rep(0, p)
            # gamN <- rep(0, p)
            # gamP[Inactive] <- (cmax - cvec[Inactive]) / (AA - a[Inactive])
            # gamN[Inactive] <- (cmax + cvec[Inactive]) / (AA + a[Inactive])

            # gamvec <- c(gamP[newj], gamN[newj])
            # temp <- min(abs(gamvec))
            # direction <- sign(gamvec[which(temp == abs(gamvec))])
            # gamma <- direction * temp
            ### METHOD 3

            ### METHOD 4
            # price <- sum(cost[Active])

            # gamP <- rep(0, p)
            # gamN <- rep(0, p)
            # gamP[Inactive] <- (cmax - cvec[Inactive]) / (AA - a[Inactive])
            # gamN[Inactive] <- (cmax + cvec[Inactive]) / (AA + a[Inactive])
            # scoreP <- rep(0, p)
            # scoreN <- rep(0, p)
            # scoreP[Inactive] <- (cmax - gamP[Inactive] * a[Inactive]) / (price + cost[Inactive])
            # scoreN[Inactive] <- (cmax - gamN[Inactive] * a[Inactive]) / (price + cost[Inactive])
            # nextmax <- max(abs(scoreP[Inactive]), abs(scoreN[Inactive]))
            # newj <- which(nextmax == c(abs(scoreP), abs(scoreN)))
            # if (newj > p)
            # {
            #     newj <- newj - p
            #     gamma <- gamN[newj]
            # } else {
            #     gamma <- gamP[newj]
            # }
            # direction <- sign(gamma)
            ### METHOD 4

            ### METHOD 5
            price <- sum(cost[Active])

            gamP <- rep(0, p)
            gamN <- rep(0, p)
            gamP[Inactive] <- (cmax - cvec[Inactive]) / (AA - a[Inactive])
            gamN[Inactive] <- (cmax + cvec[Inactive]) / (AA + a[Inactive])
            gamvec <- apply(cbind(gamP, gamN), 1, mingt0)
            scoreP <- cvec - gamvec * a
            scoreN <- (cvec - gamvec * a) / (price + cost)
            # if (exists("LASSO") & LASSO)
            # {
            #     best <- max(abs(scoreP[Inactive &! skip]))
            #     newj <- which(best == abs(scoreP))
            # } else {
                best <- max(abs(scoreN[Inactive &! skip]))
                newj <- which(best == abs(scoreN))
            # }
            gamma <- gamvec[newj]
            direction <- sign(gamma)
            ### METHOD 5

            j[newj] <- TRUE
        }

        # find gamma tilde for each j
        gammaj <- rep(0, p)
        gammaj[Active] <- -beta[k, Active] / (w * Signs)
        # If there are any gammaj that will eventually cross
        flag.cross <- FALSE
        skip <- rep(FALSE, p)
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
        mu <- mu + drop(gamma) * u
        beta[k + 1, Active] <- beta[k, Active] + drop(gamma) * w * Signs
        mul[k + 1, ] <- mu
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
            cat("Current SS: ")
            cat(var(e))
            cat("\n")

            trace.out$In <- ifelse(Inactive, "Out", "In")
            trace.out$In[!Active & j] <- "next"
            trace.out$In[skip] <- "skip"
            # trace.out$cont <- ifelse(contenders, "Y", "N")
            trace.out$cvec <- cvec
            trace.out$score <- cvec / cost
            trace.out$scoreP <- scoreP
            trace.out$scoreN <- scoreN
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
