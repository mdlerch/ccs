lars <- function(X, y)
{

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
        R <- choleski(Gram[A, A])


        u <- X[ , A] %*% w





    }



}


