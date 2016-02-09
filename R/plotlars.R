plotlars <- function(lars.object)
{
    betas <- lars.object$beta
    k <- nrow(betas)
    p <- ncol(betas)

    plot(0, 0, type = 'n', xlim = c(0, k), ylim = c(min(betas), max(betas)))
    for (v in 1:p)
    {
        lines(1:k, betas[ , v], col = v)
        lines(1:k, betas[ , v], col = v, pch = 16)
    }
}


plotlars2 <- function(lars.object, ylim = c(-30, 30))
{
    betas <- lars.object$beta
    k <- nrow(betas)
    p <- ncol(betas)

    plot(0, 0, type = 'n', xlim = c(0, k), ylim = ylim)
    for (v in 1:p)
    {
        lines(1:k, betas[ , v], col = v)
        lines(1:k, betas[ , v], col = v, pch = 16)
    }
}
