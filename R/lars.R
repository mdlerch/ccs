set.seed(19)
x <- matrix(rnorm(200), 20, 10)
y <- rnorm(20)
eps <- 0.00001

# TODO: rewrite this function
# Get updated R from choleski decomposition
updateR <- function(xnew, R = NULL, xold, Gram = FALSE, eps = 0.00001)
{
    ###Gram argument determines the nature of xnew and xold
    xtx <- if(Gram) xnew else sum(xnew^2)
    norm.xnew <- sqrt(xtx)
    if(is.null(R))
    {
        R <- matrix(norm.xnew, 1, 1)
        attr(R, "rank") <- 1
        return(R)
    }
    Xtx <- if(Gram) xold else drop(t(xnew) %*% xold)
    r <- backsolvet(R, Xtx)
    rpp <- norm.xnew^2 - sum(r^2)
    rank <- attr(R, "rank")	### check if R is machine singular
    if(rpp <= eps)
    {
        rpp <- eps
    } else {
        rpp <- sqrt(rpp)
        rank <- rank + 1
    }
    R <- cbind(rbind(R, 0), c(r, rpp))
    attr(R, "rank") <- rank
    R
}



lars <- function(x, y)
{
    # algorithm params
    maxK <- 10000

    # get quantities
    n <- length(y)
    m <- ncol(x)
    # intercept vector?
    one <- rep(1, n)

    # standardize x
    meanx <- drop( one %*% x ) / n
    normx <- sqrt(drop( one %*% x ^ 2 ))
    x <- scale(x)

    # center y
    mu <- mean(y)
    y <- scale(y, center = TRUE, scale = FALSE)

    # x not in model
    inactive <- seq(m)
    # im?
    im <- inactive

    # initials
    Cvec <- drop( t(y) %*% x )
    ssy <- sum( y ^ 2 )
    resids <- y

    # GRAM is X'X
    Gram <- t(x) %*% x

    beta <- matrix(0, maxK + 1, m)

    # active set of x
    active <- NULL

    drops <- FALSE
    # sign of correlation
    Sign <- NULL

    done <- FALSE
    k <- 0
    while (k < maxK & !done)
    {
        k  <- k + 1

        # remaining projections
        C <- Cvec[inactive]

        # largest inactive projection
        Cmax <- max(abs(C))

        # TODO: assuming drops is FALSE

        new <- abs(C) >= Cmax - eps
        C <- C[!new]
        new <- inactive[new]

        for (inew in new)
        {
            cholR <- updateR(Gram[inew, inew], R, drop(Gram[inew, active]),
                             Gram = TRUE, eps = eps)
        }





    }

}


