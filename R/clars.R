clars <- function(x, y, cost, maxk = 50, eps = 1e-6, trace = FALSE, type = 1)
{
    x <- scale(x)
    # variable setup
    n <- nrow(x)
    p <- ncol(x)

    if (missing(cost))
    {
        cost <- rep(1, p)
    }

    # current prediction
    mu <- rep(0, n)

    # sets
    Inactive <- rep(TRUE, p)
    Active <- rep(FALSE, p)

    Gram <- t(x) %*% x

    # number lars iterations.  lars will do p steps.  number of lasso steps is
    # not a fixed quantity
    k <- 0

    # number of variables currently in model (for lars, equivalent to step
    # number)
    nv <- 0

    beta <- matrix(0, nrow = maxk + 1, ncol = p)
    mul <- matrix(0, nrow = maxk + 1, ncol = n)

    # Equation 2.8
    # get first entry outside of loop
    r <- y - mu
    cvec <- t(x) %*% r
    cmax <- max(abs(cvec))
    svec <- abs(cvec) / cost
    smax <- max(svec)
    j <- svec >= smax - eps
    if (type == 3)
    {
        j <- cvec >= cmax - eps
    }
    Active <- j

    trace.out <- data.frame(var = colnames(x))

    while (nv < p & k < maxk)
    {
        k <- k + 1

        # 1. Find the next variable to add.

        # Equation 2.9
        Inactive <- !Active
        nv <- nv + 1
        # TODO: don't really like this
        r <- y - mu
        cvec <- t(x) %*% r
        cmax <- max(abs(cvec[j]))

        # 2. Find unit-vector of equal projection.
        # Following equations 2.4 through 2.6

        Signs <- sign(cvec[Active])
        # Equation 2.4
        XA <- x[ , Active] * rep(1, n) %*% t(Signs)
        # Equation 2.5
        gA <- t(XA) %*% XA
        one <- rep(1, sum(Active))
        # Equation 2.5
        AA <- 1 / sqrt(one %*% solve(gA) %*% one)
        # Equation 2.6 NOTE add the Signs to match package
        w <- AA %*% t(solve(gA) %*% one)
        # Equation 2.6
        u <- XA %*% t(w)

        # 3. Increment model fit in the direction of u.
        #  New estimate will be mu + \gamma * u where \gamma is large enough
        #  such that the next input variable will now be equally correlated.

        if (nv == p)
        {
            # cheat and just use OLS
            beta[k + 1, ] <- coef(lm(y ~ x - 1))
        } else
        {
            # Equation 2.11
            a <- t(x) %*% u
            # Equation 2.13
            gammas <- rep(0, p)
            gammas[Inactive] <- apply(
                cbind((cmax - cvec[Inactive]) / (AA - a[Inactive]),
                      (cmax + cvec[Inactive]) / (AA + a[Inactive])), 1, mingt0)

            # TODO: this makes no sense
            # Type 1 find smallest gamma / cost
            if (type == 1)
            {
                temp <- gammas / cost
                temp[temp == 0] <- max(temp) + 1
                newj <- which.min(temp)
                gamma <- gammas[newj]
                j[newj] <- TRUE
            }
            # Type 2 drop variable if crossing 0
            if (type == 2)
            {
                # start off like type 1
                temp <- gammas * cost
                temp[temp == 0] <- max(temp) + 1
                newj <- which.min(temp)
                gamma <- gammas[newj]
                # find gamma tilde for each j
                # populate with large values (for non-active variables)
                gammaj <- rep(0, p)
                gammaj[Active] <- -beta[k, Active] / (w * Signs)
                j[newj] <- TRUE
                # If there are any gammaj that will eventually cross
                if (any(gammaj > 0))
                {
                    # get the first to cross
                    temp <- gammaj
                    temp[temp <= 0] <- max(temp) + 1
                    outj <- which.min(temp)
                    gamma.tilde <- gammaj[outj]
                    if (gamma.tilde < gamma)
                    {
                        gamma <- gamma.tilde
                        j[outj] <- FALSE
                        j[newj] <- FALSE
                        if (trace) {
                            cat("\n\nVariable crossing 0, ")
                            cat(as.character(trace.out$var[outj]))
                            cat("\n")
                        }
                        nv <- nv - 2
                    }
                }
            }
            if (type == 3)
            {
                # start off like type 1
                temp <- gammas
                temp[temp == 0] <- max(temp) + 1
                newj <- which.min(temp)
                gamma <- gammas[newj]
                # find gamma tilde for each j
                # populate with large values (for non-active variables)
                gammaj <- rep(0, p)
                gammaj[Active] <- -beta[k, Active] / (w * Signs)
                j[newj] <- TRUE
                # If there are any gammaj that will eventually cross
                if (any(gammaj > 0))
                {
                    # get the first to cross
                    temp <- gammaj
                    temp[temp <= 0] <- max(temp) + 1
                    outj <- which.min(temp)
                    gamma.tilde <- gammaj[outj]
                    if (gamma.tilde < gamma)
                    {
                        gamma <- gamma.tilde
                        j[outj] <- FALSE
                        j[newj] <- FALSE
                        if (trace) {
                            cat("\n\nVariable crossing 0, ")
                            cat(as.character(trace.out$var[outj]))
                            cat("\n")
                        }
                        nv <- nv - 2
                    }
                }
            }
            if (trace)
            {
                cat("Iteration: ")
                cat(k)
                cat("\nCurrent number of variables: ")
                cat(nv)
                cat("\nTry to put in: ")
                cat(as.character(trace.out$var[newj]))
                cat("\n")
                trace.out$Active <- ifelse(Active, "In", "Out")
                # trace.out$j <- ifelse(j, "In", "Out")
                trace.out$cvec <- cvec
                trace.out$cost <- cost
                trace.out$score <- (cmax - gammas * AA) / cost
                trace.out$gamma <- gammas
                trace.out$gamrate <- gammas * cost
                trace.out$gamtilde <- gammaj
                trace.out$gamma_sel <- gamma
                trace.out$beta <- beta[k, ]
                print(trace.out)
                cat("\n\n")
                # cat("gamma: "); cat(gamma)
                # cat("cmax/A: "); cat(cmax/AA)
                # cat("\nNew score: "); cat(smax - gamma * AA); cat("\nsmax: ")
                # cat(smax); cat("\nnextmax: "); cat(max(svec[Inactive])); cat("\n")
            }
            mu <- mu + gamma * u
            beta[k + 1, Active] <- beta[k, Active] + gamma * w * Signs
            mul[k + 1, ] <- mu
            Active <- j
        }

    }

    list(beta = beta[1:(k + 1), ])
}
