x <- as.matrix(read.csv("./X.csv"))[ , -1]
y <- as.matrix(read.csv("./y.csv"))[ , -1]
source("./R/util.R")

mlars <- function(x, y, maxk = 1000, eps = 1e-6)
{
    # variable setup
    n <- nrow(x)
    p <- ncol(x)

    # current prediction
    mu <- rep(0, n)

    Inactive <- rep(TRUE, p)
    Active <- rep(FALSE, p)

    # USE GRAM
    # x <- x / sqrt(n - 1)
    # lars packages converts x to this norm then undo's when returning.
    # Not sure if there is some advantage to this, but I won't bother

    Gram <- t(x) %*% x

    # not used anywhere
    # lassocond <- FALSE

    # number lars iterations.  lars will do p steps.  number of lasso steps is
    # not a fixed quantity
    k <- 0

    # number of variables currently in model (for lars, equivalent to step
    # number)
    nv <- 0

    # x <- scale(x)

    beta <- matrix(0, nrow = maxk + 1, ncol = p)

    while (nv < p & k < maxk)
    {
        k <- k + 1

        # find the correlations with residuals
        cvec <- t(x) %*% (y - mu)
        # strongest correlation
        cmax <- max(abs(cvec))

        # which variable(s) are most correlated?
        j <- abs(cvec) >= cmax - eps

        # update Active and Inactive vectors
        Active <- Active | j
        Inactive <- !Active

        # TODO: could the increment ever be greater than 1?
        nv <- nv + 1

        # get correlation directions
        Signs <- sign(cvec[Active])

        # Equation 2.4
        XA <- x[ , Active] * rep(1, n) %*% t.default(Signs)
        # Equation 2.5
        gA <- t(XA) %*% XA
        one <- rep(1, sum(Active))
        # Equation 2.5
        AA <- 1/sqrt(one %*% solve(gA) %*% one)
        # Equation 2.6
        w <- AA %*% t(solve(gA) %*% one)
        # Equation 2.6
        u <- cbind(XA %*% t(w2)


        # old version?
        # R <- chol(Gram[Active, Active])
        # R
        # GA1 <- backsolve(R, backsolvet(R, Signs))
        # AA <- 1 / sqrt(sum(GA1 * Signs))
        # w <- AA %*% GA1
        # u <- x[ , Active] %*% t(w)

        # If all variables active go to lsq solution
        if (nv == p)
        {
            # gamma <- (cvec / AA)[1]
            beta[k + 1, Active] <- coef(lm(y ~ x - 1))
        } else
        {
            # eqn 2.11
            a <- t(x) %*% u
            # lars calc
            # u %*% x
            # eqn 2.13
            temp <- c( (cmax - cvec[Inactive]) / (AA - a[Inactive]),
                      (cmax + cvec[Inactive]) / (AA + a[Inactive]) )
            gamma <- min(temp[temp > eps], cmax / AA)
            mu <- mu + gamma * u
            beta[k + 1, Active] <- beta[k, Active] + gamma * w
        }

    }

    beta[1:(k+1), ]
}
