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

    activeMatrix <- matrix(FALSE, nrow = maxk + 1, ncol = p)
    skipMatrix <- matrix(FALSE, nrow = maxk + 1, ncol = p)
    treeMatrix <- matrix(0, nrow = maxk + 1, ncol = maxk + 1)
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
    newj <- which(j)

    tree <- 1

    gmax <- 0

    while (nv < p & k < maxk)
    {
        cat("tree: ")
        cat(tree)
        cat("\n")
        k <- k + 1
        activeMatrix[k, ] <- Active

        Inactive <- !Active
        nv <- sum(Active)
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
        if (nv == p)
        {
            gamma <- gmax
            muC <- muC + drop(gamma) * u
            betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
            mu[k, ] <- muC
            beta[k, ] <- betaC
            Active[newj] <- TRUE
            activeMatrix[k, ] <- Active
            tree <- c(k, tree)
        } else {
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

            # cat("Skip these: ")
            # cat(colnames(x)[skipMatrix[tree[1], ]])
            # cat("\n")
            OK <- Inactive & (gamvec < 2 * gmax) &! skipMatrix[tree[1], ]

            if (any(OK))
            {
                # TODO: At this point we know that we have a gamma less than 2 *
                # gmax, so first try the ones less than gmax
                under <- OK & (gamvec < gmax)
                cvecP <- cmax - gamvec * AA
                # score <- abs((cvecP) / price)
                # best <- max(abs(score[OK]))
                score <- abs((cvecP) / price)
                if (any(under))
                {
                    best <- max(abs(score[under]))
                } else {
                    best <- min(abs(score[OK]))
                }
                newj <- which(best == score)
                gamma <- gamvec[newj]
                # TODO: add lasso part later

                muC <- muC + drop(gamma) * u
                betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
                mu[k, ] <- muC
                beta[k, ] <- betaC
                Active[newj] <- TRUE
                activeMatrix[k, ] <- Active
                tree <- c(k, tree)
            } else {
                # complete the ols
                gamma <- gmax
                muC <- muC + drop(gamma) * u
                betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
                mu[k, ] <- muC
                beta[k, ] <- betaC
                Active[newj] <- TRUE
                activeMatrix[k, ] <- Active

                # revert
                if (length(tree) == 1)
                {
                    betaC <- rep(0, p)
                    muC <- rep(0, n)
                    Active <- rep(FALSE, p)
                } else {
                    toskip <- which(activeMatrix[tree[1], ] &! activeMatrix[tree[2], ])
                    skipMatrix[tree[2], toskip] <- TRUE
                    betaC <- beta[tree[2], ]
                    muC <- mu[tree[2], ]
                    Active <- activeMatrix[tree[2], ]
                    tree <- tree[-1]
                    newj <- NA
                }
            }
        }
    }

    return(list(beta = beta[1:k, ]))
}
