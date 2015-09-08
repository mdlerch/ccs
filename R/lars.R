# x <- as.matrix(read.csv("../X.csv"))[ , -1]
# y <- as.matrix(read.csv("../y.csv"))[ , -1]
# source("./util.R")

mlars <- function(x, y, maxk = 1000, eps = 1e-6)
{
    # variable setup
    n <- nrow(x)
    p <- ncol(x)

    # current position
    mu <- rep(0, n)
    Inactive <- rep(TRUE, p)
    Active <- rep(FALSE, p)

    # USE GRAM
    # x <- x / sqrt(n - 1)
    # lars packages converts x to this norm then undo's when returning.
    # Not sure if there is some advantage to this, but I won't bother

    Gram <- t(x) %*% x
    lassocond <- FALSE

    k <- 0
    nv <- 0
    loop <- TRUE

    # x <- scale(x)

    beta <- matrix(0, nrow = maxk + 1, ncol = p)

    while (nv < p & k < maxk & loop)
    {
        k <- k + 1
        cvec <- t(x) %*% (y - mu)
        cmax <- max(abs(cvec))

        j <- abs(cvec) >= cmax - eps

        Active <- Active | j
        Inactive <- !Active
        nv <- nv + 1

        Signs <- sign(cvec[Active])

        # TODO: watch out for really big stuff
        R <- chol(Gram[Active, Active])

        # TODO: figure out what is going on here!!!
        GA1 <- backsolve(R, backsolvet(R, Signs))
        AA <- 1 / sqrt(sum(GA1 * Signs))
        # eqn 2.6 b
        w <- AA %*% GA1

        # Equiangular vector
        u <- x[ , Active] %*% t(w)

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
