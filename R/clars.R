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
    r <- y - mu
    cvec <- t(x) %*% r
    cmax <- max(abs(cvec))
    svec <- abs(cvec) / cost
    smax <- max(svec)
    j <- svec >= smax - eps
    Active <- j

    trace.out <- data.frame(var = colnames(x))

    if (trace)
    {
        cat("\nIteration: ")
        cat(k)
        cat("\nCurrent number of variables: 0")
        cat("\n")
        cat("\nTry to put in next ")
        cat(as.character(trace.out$var[j]))
        cat("\n")
        trace.out$Active <- ifelse(Active, "In", "Out")
        trace.out$Active[!Active & j] <- "next"
        trace.out$cont <- 0
        trace.out$cvec <- cvec
        trace.out$cost <- cost
        trace.out$score <- cvec / cost
        trace.out$gammaP <- 0
        trace.out$gammaN <- 0
        trace.out$gamtilde <- 0
        trace.out$gamma_sel <- 0
        trace.out$beta <- 0
        print(trace.out)
        cat("\n\n")
    }

    while (nv < p & k < maxk)
    {
        k <- k + 1

        # 1. Find the next variable to add.
        Inactive <- !Active
        nv <- nv + 1
        # TODO: don't really like this
        r <- y - mu
        cvec <- t(x) %*% r
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

        # 3 options for how far to go along u (gamma).  Flow chart:
        #  a) Are there any variables more correlated than what we have in model
        #  already: TODO: should this be scaled by price?
        #    Set gamma so that these variables enter and are now equal to the
        #    already in the model variables as far as correlation.
        #  b) Find gamma until next variables are tied.
        #    lars selects smallest gamma, we select smallest gamma * cost
        #  Given an above gamma, does any variable cross zero?
        #    This suggests an option to remove a variable (i.e. a cost) so go
        #    for it. TODO: stop at that point?  keep going?

        a <- t(x) %*% u

        gamP <- rep(0, p)
        gamN <- rep(0, p)
        gamP[Inactive] <- (cmax - cvec[Inactive]) / (AA - a[Inactive])
        gamN[Inactive] <- (cmax + cvec[Inactive]) / (AA + a[Inactive])

        gammaP <- rep(0, p)
        gammaN <- rep(0, p)
        gammaP[Inactive] <- apply(cbind(gamP[Inactive], gamN[Inactive]),
                                  1, mingt0)
        gammaN[Inactive] <- apply(cbind(gamP[Inactive], gamN[Inactive]),
                                  1, minlt0)




        contenders <- rep(FALSE, p)
        contenders[which(abs(cvec) > cmax)] <- TRUE
        temp <- gammaP * cost
        temp[temp == 0] <- max(temp) + 1
        newj <- which.min(temp)
        gamma <- gammaP[newj]

        flag.contender <- FALSE
        flag.cross <- FALSE

        direction <- 1
        if (any(contenders))
        {
            tempp <- gammaP * cost
            tempp[!contenders] <- 0
            tempn <- gammaN * cost
            tempn[!contenders] <- 0
            temp <- c(tempp, abs(tempn))
            gamma <- min(temp[temp != 0])
            newj <- which(gamma == c(tempp, abs(tempn)))
            if (newj > p)
            {
                gamma <- gammaN[newj - p]
                newj <- newj - p
            } else {
                gamma <- gammaP[newj]
            }
            flag.contender <- TRUE
        }
        j[newj] <- TRUE

        # find gamma tilde for each j
        # populate with large values (for non-active variables)
        gammaj <- rep(0, p)
        gammaj[Active] <- -beta[k, Active] / (w * Signs)
        # If there are any gammaj that will eventually cross
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
            }
        }
        if (trace)
        {
            cat("\nIteration: ")
            cat(k)
            cat("\nCurrent number of variables: ")
            cat(nv)
            if (flag.contender)
            {
                cat("\nTry to put in contender ")
                cat(as.character(trace.out$var[newj]))
            } else {
                cat("\nTry to put in next ")
                cat(as.character(trace.out$var[newj]))
            }
            if (flag.cross)
            {
                cat("\nBut variable ")
                cat(as.character(trace.out$var[outj]))
                cat(" is crossing zero")
            }
            cat("\n")
            trace.out$Active <- ifelse(Active, "In", "Out")
            trace.out$Active[!Active & j] <- "next"
            trace.out$cont <- ifelse(contenders, "Y", "N")
            trace.out$cvec <- cvec
            trace.out$cost <- cost
            trace.out$score <- gammaP * cost
            trace.out$gammaP <- gammaP
            trace.out$gammaN <- gammaN
            trace.out$gamtilde <- gammaj
            trace.out$gamma_sel <- drop(gamma)
            trace.out$beta <- beta[k, ]
            print(trace.out)
            cat("\n\n")
        }
        mu <- mu + drop(gamma) * u
        beta[k + 1, Active] <- beta[k, Active] + drop(gamma) * w * Signs
        mul[k + 1, ] <- mu
        Active <- j

    }

    list(beta = beta[1:(k + 1), ])
}
