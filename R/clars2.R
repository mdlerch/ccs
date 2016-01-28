clars2 <- function(x, y, cost, maxk = 1000, eps = 1e-6, trace = FALSE)
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
    svec <- abs(cvec) / cost
    smax <- max(svec)
    j <- svec >= smax - eps


    while (nv < p & k < maxk)
    {
        k <- k + 1

        # 1. Find the next variable to add.

        if (trace) { cat("\nIteration: "); cat(k); cat("\n"); }
        #print(cbind(cvec, svec)) }

        # Equation 2.9
        Active <- Active | j
        Inactive <- !Active
        nv <- nv + 1
        # TODO: don't really like this
        cmax <- max(abs(cvec[j]))

        if (trace) { cat("selected: "); cat(colnames(x)[Active]); cat("\n") }
        if (trace) { print(cbind(cvec[Active], svec[Active])) }
        if (trace) { print(cbind(cvec[Inactive], svec[Inactive])) }

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
            gammas <- apply(
                cbind((cmax - cvec[Inactive]) / (AA - a[Inactive]),
                      (cmax + cvec[Inactive]) / (AA + a[Inactive])), 1, mingt0)
            fgammas <- gammas[which(gammas > 0)]
            best <- which.min(fgammas / cost[Inactive])
            gamma <- fgammas[best]

            # Find gamma tilde for each j (the point at which the variable will
            # cross 0)
            gammaj <- rep(1000000, p)
            gammaj[Active] <- - beta[k, Active] / (w * Signs)

            js <- which(gammaj > 0)

            lass <- FALSE
            if (length(js))
            {
                gamma.tilde <- min(gammaj[js])
                if (gamma.tilde < gamma)
                {
                    gamma <- gamma.tilde
                    out <- which(gamma == gammaj)
                    nv <- nv - 1
                    lass <- TRUE
                    if (trace)
                    {
                        cat("Variable crossing 0")
                    }
                }
            }

            j <- rep(FALSE, p)
            j[which(cumsum(Inactive) == best)[1]] <- TRUE
            if (trace)
            {
                cat("gamma: "); cat(gamma)
                cat("cmax/A: "); cat(cmax/AA)
                cat("\nNew score: "); cat(smax - gamma * AA); cat("\nsmax: ")
                cat(smax); cat("\nnextmax: "); cat(max(svec[Inactive])); cat("\n")
            }
            mu <- mu + gamma * u
            beta[k + 1, Active] <- beta[k, Active] + gamma * w * Signs
            mul[k + 1, ] <- mu
        }

    }

    list(beta = beta[1:(k + 1), ])
}
