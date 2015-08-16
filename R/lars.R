# TODO: Removed action/actions variables

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

    first.in <- integer(m)

    # initials
    Cvec <- drop( t(y) %*% x )
    ssy <- sum( y ^ 2 )
    resids <- y

    # GRAM is X'X
    Gram <- t(x) %*% x
    R <- NULL

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
            # TODO: assuming use.Gram
            R <- updateR(Gram[inew, inew], R, drop(Gram[inew, active]),
                             Gram = TRUE, eps = eps)
        }

        if (attr(R, "rank") == length(active))
        {
            # if we've reached a singularity remove that variable
            nR <- seq(length(active))
            R <- R[nR, nR, drop = FALSE]
            attr(R, "rank") <- length(active)
            ignores <- c(ignores, inew)
        } else
        {
            # TODO: could probably figure out a way to remove this part
            #  if statements are slow and this is probably going to be off
            #  for a long time
            if (first.in[inew] == 0)
            {
                # First time the variable goes in
                first.in[inew] <- k
            }
            active <- c(active, inew)
            Sign <- c(Sign, sign(Cvec[inew]))
        }

        # TODO:  what is Gi1?
        Gi1 <- backsolve(R, backsolvet(R, Sign))

        # TODO: Ignore forward.stagewise method

        A <- 1 / sqrt(sum( Gi1 * Sign))
        w <- A * Gi1





    }

}


