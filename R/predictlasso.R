predictlasso <- function(object, x, s)
{
    betas <- object$beta

    k <- nrow(betas)

    norms <- apply(betas, 1, function(x) sum(abs(x)))

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

        frac <- s - i

        m <- beta.f - beta.i

        # y = m * x + b
        beta.out <- m * frac + beta.i

    } else {
        beta.out <- betas[s == norms, ]
    }

    beta.out
}
