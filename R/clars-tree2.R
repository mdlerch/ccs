clarstree <- function(x, y, cost, maxk = 50, eps = 1e-6, trace = FALSE, costfunc = NULL)
{
    # default costfunc is just sum of used variables
    if (is.null(costfunc))
    {
        costfunc <- function(Active)
        {
            sum(cost[Active])
        }
    }

    # setup
    n <- nrow(x); p <- ncol(x)
    if (missing(cost)) cost <- rep(1, p)

    Inactive <- rep(TRUE, p)
    Active <- rep(FALSE, p)
    activeMatrix <- matrix(FALSE, nrow = maxk + 1, ncol = p)
    skipMatrix <- matrix(FALSE, nrow = maxk + 1, ncol = p)
    # how many trees have we had so far
    newtree <- 0
    newtreenum <- 0
    newtreestart <- FALSE
    newskipMatrix <- matrix(FALSE, nrow = maxk + 1, ncol = p)
    newtreeMatrix <- matrix(FALSE, nrow = maxk + 1, ncol = p)

    # distance to OLS solution (set to 0 to initialize)
    gmax <- 0
    # how deep did the deepest tree go?  For outputting treeMatrix
    maxdepth <- 1
    # setup mu and beta matrix and current values
    mu <- matrix(0, nrow = maxk + 1, ncol = n)
    muC <- rep(0, n)
    beta <- matrix(0, nrow = maxk + 1, ncol = p)
    betaC <- rep(0, p)
    # number of steps
    k <- 1
    # number of variables currently in model
    nv <- 0
    # pre-compute the gram matrix, makes things more efficient
    Gram <- t(x) %*% x

    ### STEP 1
    ##########
    x <- scale(x)
    y <- y - mean(y)
    r <- y - muC

    ### STEP 2
    ##########
    cvec <- t(x) %*% r
    cmax <- max(abs(cvec))
    svec <- abs(cvec) / cost
    smax <- max(svec)
    j <- svec >= smax - eps
    larsj <- cvec == cmax
    Active <- j
    newj <- which(j)
    r2 <- var(r) * (n - 1)

    if (which(larsj) != which(j))
    {
        cat("Not a lars selection\n")
        newtree <- 1
        newtreestart <- TRUE
        newtreenum <- 1
        newskipMatrix[newtreenum, j] <- TRUE
    } else {
        cat("clars == lars\n")
    }

    trace.out <- data.frame(var = colnames(x), cost = cost)
    if (trace)
    {
        cat("\nIteration: ")
        cat(k)
        cat("\nCurrent number of variables: 0")
        cat("\n")
        cat("\nTry to put in next ")
        cat(as.character(trace.out$var[j]))
        cat("\n")
        trace.out$In <- ifelse(Active, "In", "Out")
        trace.out$In[!Active & j] <- "next"
        # trace.out$cont <- 0
        trace.out$cvec <- cvec
        trace.out$score <- cvec / cost
        trace.out$gammaP <- 0
        trace.out$gammaN <- 0
        trace.out$gamtilde <- 0
        trace.out$beta <- 0
        print(trace.out)
        cat("\n\n")
    }

    skip <- rep(FALSE, p)
    while (nv < p & k < maxk)
    {

        k <- k + 1
        activeMatrix[k, ] <- Active
        Inactive <- !Active
        nv <- sum(Active)
        r <- y - muC
        cvec <- t(x) %*% r
        cmax <- max(abs(cvec[Active]))

        ### STEP 3
        ##########
        Signs <- sign(cvec[Active])
        XA <- x[ , Active] * rep(1, n) %*% t(Signs)
        gA <- t(XA) %*% XA
        one <- rep(1, sum(Active))
        AA <- 1 / sqrt(one %*% solve(gA) %*% one)
        w <- AA %*% t(solve(gA) %*% one)
        u <- XA %*% t(w)
        a <- t(x) %*% u
        # to OLS solution
        gmax <- drop(cmax / AA)

        if (nv == p)
        {
            gamma <- gmax
        } else {
            # calculate each variables next price
            # TODO: can this be made more efficient?
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

            temp <- c(gamP[Inactive], gamN[Inactive])

            larsgamm <- min(temp[temp > eps], cmax / AA)

            # TODO: TODO: TODO:
            # TODO: TODO: TODO:
            # SOMETHING SEEMS WRONG WITH SCORE
            # ALGORITHM IS ALLOWING gamvec > gmax to have greater preference
            # than gamvec < gmax

            # First choice, gamma between 0 and cmax/A
            gamvec <- apply(cbind(gamP, gamN), 1, mingt0)

            if (newtreestart)
            {
                OK <- Inactive & (gamvec < gmax) &! newskipMatrix[newtreenum, ] &! (gamvec == 0)
            } else {
                OK <- Inactive & (gamvec < gmax) &! (gamvec == 0)
            }

            if (any(OK))
            {
                # STEP 4
                theta <- acos(t(r) %*% u / sqrt(r2))
                r2j <- r2 + gamvec^2 - 2 * sqrt(r2) * gamvec * cos(theta)
                cvecP <- (cmax - gamvec * AA) * sqrt(r2j)
                score <- abs((cvecP) / price)
                # TODO: this includes all "acceptable"
                best <- max(abs(score[OK]))
                newj <- which(best == score)
                gamma <- gamvec[newj]

                newtreestart <- FALSE
                # let me know if I'm not making the lars pick
                if (gamma != larsgamm)
                {
                    cat("Not a lars selection")
                    if (newtree == 0)
                    {
                        cat(".  Starting a new tree at ")
                        cat(k)
                        cat(".\n")
                        newtree <- k
                        newtreestart <- TRUE
                    } else {
                        cat(", but already in a tree.\n")
                    }
                }
            } else {
                # need to revert the tree.

                if (newtreestart == TRUE)
                {
                    cat("Major problems\n")
                }

                # complete the ols
                cat("\n")
                cat("Completing the OLS")
                cat("\n")
                gamma <- gmax
                muC <- muC + drop(gamma) * u
                betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
                mu[k, ] <- muC
                beta[k, ] <- betaC
                activeMatrix[k, ] <- Active

                # revert
                cat("\n")
                cat("Reverting the tree")
                cat("\n")
                print(newskipMatrix[newtreenum, ])

                if (newtree == 1)
                {
                    # if we go back to step 1
                    cat("Back to the start.\n")
                    betaC <- beta[newtree, ]
                    muC <- mu[newtree, ]
                    Active <- activeMatrix[newtree, ]
                    # repeat steps 1 and 2
                    ### STEP 1
                    ##########
                    r <- y - muC

                    ### STEP 2
                    ##########
                    cvec <- t(x) %*% r
                    cmax <- max(abs(cvec))
                    svec <- abs(cvec) / cost
                    smax <- max(svec[!newskipMatrix[newtreenum, ]])
                    j <- svec == smax
                    larsj <- cvec == cmax
                    Active <- j
                    newj <- which(j)
                    r2 <- var(r) * (n - 1)

                    cat("LARS selection: ")
                    cat(which(larsj))
                    cat("\nclars selection: ")
                    cat(which(j))
                    cat("\n")

                    if (which(larsj) != which(j))
                    {
                        cat("Not a lars selection\n")
                        newtreenum <- newtreenum + 1
                        newskipMatrix[newtreenum, j] <- TRUE
                        cat("No longer allowed at this node: ")
                        newtreestart <- TRUE
                    } else {
                        cat("clars == lars\n")
                        # we do start with a lars step
                        newtree <- 0
                        newtreestart <- FALSE
                    }

                } else if (newtree != 0) {
                    cat("Back to step ")
                    cat(newtree)
                    cat(".\n")
                    betaC <- beta[newtree - 1, ]
                    muC <- mu[newtree - 1, ]
                    Active <- activeMatrix[newtree - 1, ]
                    # TODOTODO might need to actually go back to k - 1
                    newtreestart <- TRUE
                }

                # TODO: do i need this?
                next
            }

            direction <- sign(gamma)

            j[newj] <- TRUE
        }

        skip <- rep(FALSE, p)
        # find gamma tilde for each j
        gammaj <- rep(0, p)
        # TODO: should this be betaC?
        gammaj[Active] <- -betaC[Active] / (w * Signs)
        # If there are any gammaj that will eventually cross
        flag.cross <- FALSE
        # TODO: gamma.tilde could also be negative if we choose a negative
        # gamma...
        if (any(direction * gammaj > 0))
        {
            if (nv == p)
            {
                gamma <- gmax
                temp <- gammaj
                temp[temp <= 0] <- max(temp) + 1
                outj <- which.min(temp)
                gamma.tilde <- gammaj[outj]
                if (gamma.tilde < gamma)
                {
                    gamma <- gamma.tilde
                    j[outj] <- FALSE
                    nv <- nv - 2
                }
            }
            # get the first to cross
            temp <- gammaj
            temp[temp <= 0] <- max(temp) + 1
            outj <- which.min(temp)
            gamma.tilde <- gammaj[outj]
            if (gamma.tilde < gamma)
            {
                flag.cross <- TRUE
                gamma <- gamma.tilde
                j[outj] <- FALSE
                j[newj] <- FALSE
                nv <- nv - 2
                skip[outj] <- TRUE
            }
        }

        # If gamma is longer than cmax/AA make a pit stop along the way
        if (drop(gamma) > gmax)
        {
            # short gamma step
            shortgamma <- gmax
            mu[k, ] <- muC + drop(shortgamma) * u
            beta[k, Active] <- betaC[Active] + drop(shortgamma) * w * Signs
            # whole gamma step
            muC <- muC + drop(gamma) * u
            betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
            mu[k + 1, ] <- muC
            beta[k + 1, Active] <- betaC[Active]
            # TODO: do i want this?
            activeMatrix[k + 1, ] <- Active
            k <- k + 1
        } else {
            muC <- muC + drop(gamma) * u
            betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
            mu[k, ] <- muC
            beta[k, Active] <- betaC[Active]
        }


        if (any(skip))
        {
            if (trace)
            {
                if (beta[k, skip] != 0)
                {
                    cat("Should have been zero.  Was: ")
                    cat(beta[k, skip])
                }
            }
            beta[k, skip] <- 0
        }
        if (trace)
        {
            cat("\nIteration: ")
            cat(k)
            cat("\nLARS eq: ")
            cat(newtree == 0)
            cat("\nAt a node: ")
            cat(newtreestart)
            cat("\nCurrent number of variables: ")
            cat(nv)
            # if (flag.contender)
            # {
            #     cat("\nNext in via contender: ")
            #     cat(as.character(trace.out$var[newj]))
            # } else {
                cat("\nNext in via lars: ")
                cat(as.character(trace.out$var[newj]))
            # }
            if (flag.cross)
            {
                cat("\nBut variable ")
                cat(as.character(trace.out$var[outj]))
                cat(" is crossing zero")
            }
            cat("\n")
            cat("Selected gamma: ")
            cat(drop(gamma))
            cat("\n")
            cat("cmax / A: ")
            cat(gmax)
            cat("\n")
            cat("Current SS: ")
            cat(var(r))
            cat("\n")

            trace.out$In <- ifelse(Inactive, "Out", "In")
            trace.out$In[!Active & j] <- "next"
            trace.out$In[skip] <- "skip"
            # trace.out$cont <- ifelse(contenders, "Y", "N")
            trace.out$cvec <- cvec
            trace.out$score <- score
            trace.out$gammaP <- gamP
            trace.out$gammaN <- gamN
            trace.out$gamtilde <- gammaj
            ww <- rep(0, p)
            ww[Active] <- w
            trace.out$w <- ww
            trace.out$beta <- beta[k, ]
            print(trace.out)
            cat("\n\n")
        }
        Active <- j
        activeMatrix[k, ] <- Active
    }

    list(mu = mu[1:k, ], beta = beta[1:k, ],
         activeMatrix = activeMatrix[1:k, ])
}
