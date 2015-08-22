
backsolvet <- function(r, x, k)
{
    if(missing(k))
        k <- nrow(r)
    r <- t(r)[k:1, k:1, drop = F]
    x <- as.matrix(x)[k:1,  , drop = F]
    x <- backsolve(r, x)
    drop(x[k:1,  ])
}



lars <- function(X, y)
{

    set.seed(3)
    X <- matrix(rnorm(500), 50, 10)
    b <- rnorm(10)
    y <- X %*% b + rnorm(50)

    X <- scale(X)
    y <- y - mean(y)

    eps <- 1e-6

    # variable setup
    n <- nrow(X)
    p <- ncol(X)

    maxk <- 512 * p

    # current position
    mu <- rep(0, n)
    Inactive <- 1:p
    Active <- NULL

    # USE GRAM
    Gram <- t(X) %*% X
    lassocond <- FALSE

    k <- 0
    nv <- 0
    loop <- TRUE

    beta <- matrix(0, nrow = p, ncol = maxk)

    while (nv < p & k < maxk & loop)
    {
        k <- k + 1
        cvec <- t(X) %*% (y - mu)
        cmax <- max(abs(cvec))

        j <- abs(cvec) >= cmax - eps

        if (!lasssocond)
        {
            Active <- c(Active, j)
            Inactive <- Inactive[-which(Inactive == j)]
            nv <- nv + 1
        }

        Signs <- sign(cvec[Active])

        # TODO: watch out for really big stuff
        R <- chol(Gram[Active, Active])

        # TODO: figure out what is going on here!!!
        GA1 <- backsolve(R, backsolvet(R, Signs))
        AA <- 1 / sqrt(sum(GA1 * Signs))
        w <- AA %*% GA1

        # Equiangular vector
        u <- X[ , Active] %*% w

        # If all variables active go to lsq solution
        if (nv == p)
        {
            gamma <- cvec / AA
        } else
        {
            # eqn 2.11
            a <- t(X) %*% u
            # eqn 2.13
            temp <- c( (cmax - cvec[Inactive]) / (AA - a[Inactive]),
                       (cmax + cvec[Inactive]) / (AA + a[Inactive]) )
            gamma <- min(temp[temp > eps], Cmax / AA)
        }

        # TODO: there is a lasso step here

        # eqn 2.12
        mu = mu + gamma  %*% u

        # TODO: what is Cardi?  assume "empty" Assume stop == 0

        beta <- beta[Active, k + 1] <- beta[Active, k] + gamma %*% w




        # TODO: there is a lasso step here








    }



}


