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
    inActive <- seq(m)
    # im?
    im <- inActive

    first.in <- integer(m)

    # initials
    Cvec <- drop( t(y) %*% x )
    ssy <- sum( y ^ 2 )
    resids <- y

    # GRAM is X'X
    Gram <- t(x) %*% x
    R <- NULL

    beta <- matrix(0, maxK + 1, m)

    # active set -- script{A}
    Active <- NULL

    # ignores
    ignores <- NULL

    drops <- FALSE
    # sign of correlation
    Sign <- NULL

    done <- FALSE
    k <- 0
    while (k < maxK & !done)
    {
        k  <- k + 1

        # remaining projections
        C <- Cvec[inActive]
        # largest inActive projection
        Cmax <- max(abs(C))

        # new index to add
        new <- abs(C) >= Cmax - eps

        # new is now Active
        C <- C[!new]
        # get index of new
        new <- inActive[new]
        Sign <- c(Sign, sign(Cvec[new]))

        Active <- c(Active, new)
        gActive_inv <- solve(t(x[ , Active]) %*% x[ , Active])


        A <- 1 / sqrt( sum( gActive_inv * Sign ) )
        w <- drop(A * gActive_inv)
        u <- x[, Active] %*% w

        # TODO: not working
        a <- drop(t(x) %*% t(u))




        # TODO: assuming use.gram

        if (length(Active) >= min(n - 1, m - length(ignores)))
        {
            gamhat <- Cmax / A
        } else
        {
            a <- drop(w %*% Gram[Active, -c(Active, ignores), drop = FALSE])
            gam <- c( (Cmax - C) / (A - a), (Cmax + C) / (A + a))
            gamhat <- min(gam[gam > eps], Cmax / A)
        }

        # TODO: assuming type == "lasso"
        dropid <- NULL
        b1 <- beta[k, Active]

        z1 <- -b1 / w
        zmin <- min(z1[z1 > eps], gamhat)

        if (zmin < gamhat)
        {
            gamhat <- zmin
            drops <- z1 == zmin
        } else
        {
            drops <- FALSE
        }

        beta[k + 1, ] <- beta[k, ]
        beta[k + 1, Active] <- beta[k + 1, Active] + gamhat * w

        #TODO: assume use.Gam
        Cvec <- Cvec - gamhat * Gram[ , Active, drop = FALSE] %*% w

        Gamrat <- 




    }

}


        # Maybe more efficient?
        # for (inew in new)
        # {
        #     # TODO: assuming use.Gram
        #     R <- updateR(Gram[inew, inew], R, drop(Gram[inew, Active]),
        #                      Gram = TRUE, eps = eps)
        # }

        # if (attr(R, "rank") == length(Active))
        # {
        #     # if we've reached a singularity remove that variable
        #     nR <- seq(length(Active))
        #     R <- R[nR, nR, drop = FALSE]
        #     attr(R, "rank") <- length(Active)
        #     ignores <- c(ignores, inew)
        # } else
        # {
        #     # TODO: could probably figure out a way to remove this part
        #     #  if statements are slow and this is probably going to be off
        #     #  for a long time
        #     if (first.in[inew] == 0)
        #     {
        #         # First time the variable goes in
        #         first.in[inew] <- k
        #     }
        #     Active <- c(Active, inew)
        #     Sign <- c(Sign, sign(Cvec[inew]))
        # }

        # # TODO: understand Gi1 (I think it is the same as script{g}_A^{-1}
        # Gi1 <- backsolve(R, backsolvet(R, Sign))
        # OR:
        # Gi1 <- solve(t(x[ , new]) %*% x[ , new])

