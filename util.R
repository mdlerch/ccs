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
    } else
    {
        rpp <- sqrt(rpp)
        rank <- rank + 1
    }
    R <- cbind(rbind(R, 0), c(r, rpp))
    attr(R, "rank") <- rank
    R
}

# TODO: rewrite this function
backsolvet <- function(r, x, k)
{
    if(missing(k))
        k <- nrow(r)
    r <- t(r)[k:1, k:1, drop = F]
    x <- as.matrix(x)[k:1,  , drop = F]
    x <- backsolve(r, x)
    drop(x[k:1,  ])
}
