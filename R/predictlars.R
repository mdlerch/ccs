predictlars <- function(object, x, s)
{
    betas <- object$beta

    k <- nrow(betas)

    if (s < 1 || s > k)
    {
        stop("s is not valid")
    }

    i <- floor(s)
    f <- ceiling(s)

    if (i != f)
    {
        beta.i <- betas[i, ]
        beta.f <- betas[f, ]

        frac <- s - i

        m <- beta.f - beta.i

        # y = m * x + b
        beta.out <- m * frac + beta.i

    } else {
        beta.out <- betas[i, ]
    }

    x %*% beta.out
}
