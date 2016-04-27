cls <- function(x, y, cost, maxk = 50, eps = 1e-6, trace = FALSE, costfunc = NULL)
{
    # default costfunc is just sum of used variables
    if (is.null(costfunc))
    {
        costfunc <- function(Active)
        {
            sum(cost[Active])
        }
    }
    x <- scale(x)
    # variable setup
    n <- nrow(x); p <- ncol(x)

    if (missing(cost)) cost <- rep(1, p)

    # sets
    Inactive <- rep(TRUE, p)
    Active <- rep(FALSE, p)

    activeMatrix <- matrix(FALSE, nrow = maxk, ncol = p)
    skipMatrix <- matrix(FALSE, nrow = maxk, ncol = p)
    tree <- NULL

    Gram <- t(x) %*% x

    # number of steps
    k <- 1
    # number of variables currently in model
    nv <- 0

    mu <- matrix(0, nrow = maxk + 1, ncol = n)
    muC <- rep(0, n)
    beta <- matrix(0, nrow = maxk + 1, ncol = p)
    betaC <- rep(0, p)

    # Equation 2.8
    e <- y - muC
    cvec <- t(x) %*% e
    cmax <- max(abs(cvec))
    svec <- abs(cvec) / cost
    smax <- max(svec)
    j <- svec >= smax - eps
    Active <- j

    tree <- 1

    while (nv < p & k < maxk)
    {
        k <- k + 1
        activeMatrix[k, ] <- Active

        Inactive <- !Active
        nv <- nv + 1
        e <- y - muC
        cvec <- t(x) %*% e
        cmax <- max(abs(cvec[Active]))

        Signs <- sign(cvec[Active])
        XA <- x[ , Active] * rep(1, n) %*% t(Signs)
        gA <- t(XA) %*% XA
        one <- rep(1, sum(Active))
        AA <- 1 / sqrt(one %*% solve(gA) %*% one)
        w <- AA %*% t(solve(gA) %*% one)
        u <- XA %*% t(w)

        a <- t(x) %*% u

        ## to OLS solution
        gmax <- drop(cmax / AA)

        price <- rep(0, p)
        for (i in 1:p)
        {
            ActiveN <- Active
            ActiveN[i] <- TRUE
            price[i] <- costfunc(ActiveN)
        }
        # Two gamma calculations
        gamP <- rep(0, p)
        gamN <- rep(0, p)
        gamP[Inactive] <- (cmax - cvec[Inactive]) / (AA - a[Inactive])
        gamN[Inactive] <- (cmax + cvec[Inactive]) / (AA + a[Inactive])

        # First choice, gamma between 0 and cmax/A
        gamvec <- apply(cbind(gamP, gamN), 1, mingt0)






    }
