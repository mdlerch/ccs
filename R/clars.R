clars <- function(x, y, cost, maxk = 1000, eps = 1e-6, trace = TRUE)
{
    x <- scale(x)
    # variable setup
    n <- nrow(x)
    p <- ncol(x)

    if (missing(cost))
    {
        cost <- rep(1, p)
    }

    cost <- cost / sum(cost)


    # current prediction
    mu <- rep(0, n)

    r <- y - mu
    C <- t(x) %*% r / sd(r) / (n - 1)

    CM <- max(abs(C))
    price <- cost * CM

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

    while (nv < p & k < maxk)
    {
        k <- k + 1

        # 1. Find the next variable to add.
        # Calculate the projections (correlations) along the residuals (y - mu)
        # on each x. Find the largest of these projections and add that variable
        # to the to the active list.
        # TODO: is it possible that more than one variable could be added?

        # Equation 2.8
        r <- y - mu
        # cvec <- t(x) %*% r / sd(r) / (n - 1)
        cvec <- t(x) %*% r
        price <- cost * max(cvec)


        # score vector
        svec <- abs(cvec) - price

        if (trace)
        {
            cat("\nIteration: ")
            cat(k)
            cat("\n")
            print(cbind(cvec, svec))
        }

        # Equation 2.9
        smax <- max(svec)
        j <- svec >= smax - eps
        Active <- Active | j
        Inactive <- !Active
        nv <- nv + 1

        if (trace)
        {
            cat("selected: ")
            cat(colnames(x)[j])
            cat("\n")
        }

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

            temp <- c((smax - svec[Inactive]) / (AA - a[Inactive]),
                    (smax + cvec[Inactive] + price[Inactive]) / (AA - a[Inactive]))
            # TODO what is cmax / AA option?

            gamma <- min(temp[temp > eps])
            if (trace)
            {
                cat("gamma: ")
                cat(which(temp == gamma))
                cat(" ")
                cat(gamma)
                cat("\nNew score: ")
                cat(smax - gamma * AA)
                cat("\nsmax: ")
                cat(smax)
                cat("\nnextmax: ")
                cat(max(svec[Inactive]))
                cat("\n")
            }
            mu <- mu + gamma * u
            beta[k + 1, Active] <- beta[k, Active] + gamma * w * Signs
            mul[k + 1, ] <- mu
        }

    }

    list(beta = beta[1:(k + 1), ])
}
