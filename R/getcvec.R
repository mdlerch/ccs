getcvec <- function(y, mu, x)
{
    e <- y - muC
    t(x) %*% e
}

getu <- function(x, Active, cvec)
{
    Signs <- sign(cvec[Active])
    XA <- x[ , Active] * rep(1, n) %*% t(Signs)
    gA <- t(XA) %*% XA
    one <- rep(1, sum(Active))
    AA <- 1 / sqrt(one %*% solve(gA) %*% one)
    w <- AA %*% t(solve(gA) %*% one)
    XA %*% t(w)
}

getgamma <- function(cmax, cvec, AA, a, Active, Inactive, skip, price)
{
            gamP <- rep(0, p)
            gamN <- rep(0, p)
            gamP[Inactive] <- (cmax - cvec[Inactive]) / (AA - a[Inactive])
            gamN[Inactive] <- (cmax + cvec[Inactive]) / (AA + a[Inactive])

            # First choice, gamma between 0 and cmax/A
            gamvec <- apply(cbind(gamP, gamN), 1, mingt0)
            legal <- gamvec < gmax
            nonzero <- gamvec != 0

            # Are there any variables that meet conditions?
            OK <- Inactive & legal & nonzero &! skip

            if (any(OK))
            {
                tree <- c(k - 1, tree)
                cvecP <- cmax - gamvec * AA
                score <- abs((cvecP) / price)
                best <- max(abs(score[OK]))
                newj <- which(best == score)
                gamma <- gamvec[newj]
                return(gamma)

            } else {
                return(FALSE)
            }
}
