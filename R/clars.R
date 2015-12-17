# x <- as.matrix(read.csv("./X.csv"))[ , -1]
# y <- as.matrix(read.csv("./y.csv"))[ , -1]

# source("./R/util.R")

library(lars)
data(diabetes)
x <- scale(diabetes$x)
y <- diabetes$y
maxk <- 15
eps <- 1e-9



clars <- function(x, y, cost, maxk = 1000, eps = 1e-6)
{
    x <- scale(x)
    # variable setup
    n <- nrow(x)
    p <- ncol(x)

    if missing(cost)
    {
        cost <- rep(1, p)
    }


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

        # convert to cost per
        cvec <- cvec / cost * Inactive

        # Equation 2.9
        cmax <- max(abs(cvec))
        j <- abs(cvec) >= cmax - eps
        Active <- Active | j
        Inactive <- !Active
        nv <- nv + 1

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

        if (nv == p)
        {
            # cheat and just use OLS
            beta[k + 1, Active] <- coef(lm(y ~ x - 1))
        } else
        {
            # Equation 2.11
            a <- t(x) %*% u
            # Equation 2.13
            temp <- c((cmax - cvec[Inactive]) / (AA - a[Inactive]),
                      (cmax + cvec[Inactive]) / (AA + a[Inactive]))
            gamma <- min(temp[temp > eps], cmax / AA)
            mu <- mu + gamma * u
            beta[k + 1, Active] <- beta[k, Active] + gamma * w * Signs
            mul[k + 1, ] <- mu
        }

    }

    list(beta = beta[1:(k + 1), ])
}
