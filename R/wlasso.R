mwlasso <- function(x, y, w, maxk = 1000, eps = 1e-6)
{
    if (missing(w))
    {
        w <- 1 / abs(coef(lm(y ~ x - 1)))
    }
    we <- w

    x <- scale(x, scale = w)
    # variable setup
    n <- nrow(x)
    p <- ncol(x)

    # current prediction
    mu <- rep(0, n)

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

    # Doing a lasso step?
    lass <- FALSE

    while (nv < p & k < maxk)
    {
        k <- k + 1

        # 1. Find the next variable to add.
        # Calculate the projections (correlations) along the residuals (y - mu)
        # on each x. Find the largest of these projections and add that variable
        # to the to the active list.
        # TODO: is it possible that more than one variable could be added?

        # Equation 2.8
        cvec <- t(x) %*% (y - mu)
        if (lass)
        {
            cmax <- max(abs(cvec[-out]))
            Active[out] <- FALSE
            Inactive <- !Active
        } else {
            # Equation 2.9
            cmax <- max(abs(cvec))
            j <- abs(cvec) >= cmax - eps
            Active <- Active | j
            Inactive <- !Active
            nv <- nv + 1
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
        AA <- 1/sqrt(one %*% solve(gA) %*% one)
        # Equation 2.6 NOTE add the Signs to match package
        w <- AA %*% t(solve(gA) %*% one)
        # Equation 2.6
        u <- XA %*% t(w)

        # 2. Find unit-vector of equal projection.

        # This is another method to do step 2, may be more efficient than above.
        # R <- chol(Gram[Active, Active])
        # GA1 <- backsolve(R, backsolvet(R, Signs))
        # AA2 <- 1 / sqrt(sum(GA1 * Signs))
        # w2 <- AA2 %*% GA1
        # u2 <- x[ , Active] %*% t(w2)

        # 3. Increment model fit in the direction of u.
        #  New estimate will be mu + \gamma * u where \gamma is large enough
        #  such that the next input variable will now be equally correlated.

        # Equation 2.11
        a <- t(x) %*% u
        # Equation 2.13
        temp <- c((cmax - cvec[Inactive]) / (AA - a[Inactive]),
                  (cmax + cvec[Inactive]) / (AA + a[Inactive]))
        gamma <- min(temp[temp > eps], cmax / AA)

        gammaj <- rep(1000000, p)
        gammaj[Active] <- - beta[k, Active] / (w * Signs)

        js <- which(gammaj > 0)

        lass <- FALSE
        if (length(js) > 0)
        {
            gamma.tilde <- min(gammaj[js])
            if (gamma.tilde < gamma)
            {
                gamma <- gamma.tilde
                out <- which(gamma == gammaj)
                nv <- nv - 1
                lass <- TRUE
            }
        }

        mu <- mu + gamma * u
        beta[k + 1, Active] <- beta[k, Active] + gamma * w * Signs
        mul[k + 1, ] <- mu

    }

    return(beta[1:(k + 1), ])
}
