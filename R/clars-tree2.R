clarstree <- function(x, y, cost, maxk = 50, eps = 1e-6, trace = FALSE, costfunc = NULL, ofgmax = 1)
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
    newtree <- 0
    # how many trees have we had so far
    newtreenum <- 0
    newtreestart <- FALSE
    newskipMatrix <- matrix(FALSE, nrow = maxk + 1, ncol = p)
    newtreeMatrix <- matrix(FALSE, nrow = maxk + 1, ncol = p)
    newtreebeta <- matrix(FALSE, nrow = maxk + 1, ncol = p)
    newtreemu <- matrix(FALSE, nrow = maxk + 1, ncol = n)


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

    cat("lars choice: ")
    cat(which(larsj))
    cat("\nclars choice: ")
    cat(which(j))
    cat("\n")

    if (which(larsj) != which(j))
    {
        newtree <- 1
        # newtreestart <- TRUE
        newtreenum <- 1
        newskipMatrix[newtreenum, j] <- TRUE
        newtreebeta[newtreenum, ] <- betaC
        newtreemu[newtreenum, ] <- muC
        cat("Recording node at ")
        cat(k)
        cat("\nclars skips for this node:\n")
        print(newskipMatrix[newtreenum, ])
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

        cat("\n")

        cat("Iteration: ")
        cat(k)
        cat("\n")
        if (newtree != 0 & newtreestart)
        {
            cat("At the node of a clars tree\n")
            cat("with beta =\n")
            print(betaC)
        } else if (newtree != 0 & !newtreestart) {
            cat("In a clars tree\n")
        } else {
            cat("Following lars\n")
        }

        activeMatrix[k, ] <- Active
        Inactive <- !Active
        nv <- sum(Active)
        r <- y - muC
        cvec <- t(x) %*% r
        cmax <- max(abs(cvec[Active]))

        # cat("Active\n")
        # print(Active)
        # print(activeMatrix)
        # cat("beta\n")
        # print(beta)

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

            larsgamm <- min(temp[temp > eps], cmax / AA)
            # First choice, gamma between 0 and cmax/A
            gamvec <- apply(cbind(gamP, gamN), 1, mingt0)


            larsj <- which(larsgamm == gamvec)
            cat("lars choice: ")
            cat(larsj)
            cat("\n")
            # cat("newtreestart\n")
            # print(newtreestart)
            # cat("Inactive\n")
            # print(Inactive)
            # cat("gamvec < gmax\n")
            # print(gamvec < gmax)
            # cat("newskipMatrix[newtreenum, ]\n")
            # print(newskipMatrix[newtreenum, ])
            # cat("gamvec\n")
            # print(gamvec)
            if (newtreestart)
            {
                OK <- Inactive & (gamvec < ofgmax*gmax) &! newskipMatrix[newtreenum, ] &! (gamvec == 0)
            } else {
                OK <- Inactive & (gamvec < ofgmax*gmax) &! (gamvec == 0)
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

                cat("clars choice: ")
                cat(newj)
                cat("\n")

                # let me know if I'm not making the lars pick
                if (gamma != larsgamm)
                {
                    if (newtreestart)
                    {
                        newskipMatrix[newtreenum, newj] <- TRUE
                        cat("Updating skips to\n")
                        print(newskipMatrix[newtreenum, ])
                        newtreestart <- FALSE
                    }
                    if (newtree == 0)
                    {
                        cat("Starting a new tree at ")
                        cat(k)
                        cat(".\n")
                        newtree <- k
                        newtreenum <- newtreenum + 1
                        newskipMatrix[newtreenum, newj] <- TRUE
                        newtreebeta[newtreenum, ] <- betaC
                        newtreemu[newtreenum, ] <- muC

                        cat("clars skips for this node are\n")
                        print(newskipMatrix[newtreenum, ])
                        newtreestart <- FALSE
                    } else {
                        cat("Already in a tree based on node ")
                        cat(newtree)
                        cat("\n")
                        newtreestart <- FALSE
                    }
                } else {
                    if (newtreestart)
                    {
                        cat("Now aligned with lars with beta = \n")
                        print(betaC)
                        newtree <- 0
                        newtreestart <- FALSE
                    }
                }
            } else {
                # need to revert the tree.
                cat("no clars choice available\n")

                if (newtreestart == TRUE)
                {
                    stop("MAJOR PROBLEMS 1\n")
                }

                # complete the ols
                # TODOTODO: need to k+=1
                cat("Completing the OLS")
                cat("\n")
                gamma <- gmax
                muC <- muC + drop(gamma) * u
                betaC[Active] <- betaC[Active] + drop(gamma) * w * Signs
                mu[k, ] <- muC
                beta[k, ] <- betaC
                activeMatrix[k, ] <- Active

                cat("Back to step ")
                cat(newtree)
                cat("\n")

                # revert
                if (newtree == 1)
                {
                    # if we go back to step 1
                    betaC <- beta[newtree, ]
                    muC <- mu[newtree, ]
                    Active <- activeMatrix[newtree, ]
                    j <- Active
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

                    cat("lars choice: ")
                    cat(which(larsj))
                    cat("\nclars choice: ")
                    cat(which(j))
                    cat("\n")

                    if (which(larsj) != which(j))
                    {
                        newskipMatrix[newtreenum, j] <- TRUE
                        cat("clars skips for this node:\n")
                        print(newskipMatrix[newtreenum, ])
                        newtreestart <- FALSE
                    } else {
                        cat("Hey, Now aligned with lars with beta = \n")
                        print(betaC)
                        cat("\n")
                        # we do start with a lars step
                        newtree <- 0
                        newtreestart <- FALSE
                    }

                } else if (newtree != 0) {
                    cat("Here's where problem lies\n")
                    print(newskipMatrix[newtreenum, ])


                    # not newtree - 1 but previous newtree?
                    betaC <- newtreebeta[newtreenum, ]
                    muC <- newtreemu[newtreenum, ]
                    Active <- activeMatrix[newtree, ]
                    j <- Active
                    # TODOTODO might need to actually go back to k - 1
                    newtreestart <- TRUE
                    cat("The beta get's set to:\n")
                    print(betaC)
                    cat("Should it be:\n")
                    print(newtreebeta[newtreenum, ])

                }

                # TODO: do i need this? I think so
                next
            }

            direction <- sign(gamma)

            j[newj] <- TRUE
        }


        print(newj)
        print(gamma)


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

        # TODOTODO this should happen any more
        # If gamma is longer than cmax/AA make a pit stop along the way
        if (drop(gamma) > gmax)
        {
            cat("THIS SHOULD NOT BE HAPPENING\n")
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
    }

    list(mu = mu[1:k, ], beta = beta[1:k, ],
         activeMatrix = activeMatrix[1:k, ])
}
