lars <- function(x, y)
{
    # algorithm params
    maxK <- 10000

    # get quantities
    n <- length(y)
    p <- ncol(x)

    # standardize x
    x <- scale(x)

    # center y
    mu <- mean(y)
    y <- scale(y, center = TRUE, scale = FALSE)

    done <- FALSE
    k <- 0
    while (k < maxK & !done)
    {
        k  <- k + 1

        C <- t(y) %*% x



    }

}


