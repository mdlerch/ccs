predictlarsbeta <- function(object, x, s)
{
    betas <- object$betas

    p <- ncol(betas)
    k <- nrow(betas)

    if (s < 0 || s > k)
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

    beta.out <- m * frac + beta.i

    } else {
        beta.out <- betas[i, ]
    }

    beta.out
}

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

predictlasso <- function(object, x, s)
{
    betas <- object$beta

    k <- nrow(betas)

    norms <- apply(abs(betas), 1, sum)

    if (s < norms[1] || s > norms[k])
    {
        stop("s is not valid")
    }


    if (!(s %in% norms))
    {
        sorted <- sort(c(s, norms))
        i <- which(s == sorted) - 1
        f <- i + 1
        beta.i <- betas[i, ]
        beta.f <- betas[f, ]

        m <- (beta.f - beta.i) / (norms[f] - norms[i])

        frac <- s - norms[i]

        # y = m * x + b

        # dx1 <- s - norms[i]
        # dx2 <- norms[f] - s
        # dx <- norms[f] - norms[i]
        # beta.out <- dx2 / dx * beta.i + dx1 / dx * beta.f

        beta.out <- m * frac + beta.i

    } else {
        beta.out <- betas[s == norms, ]
    }

    beta.out
}
