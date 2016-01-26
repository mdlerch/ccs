library(lars)
data(diabetes)
x <- scale(diabetes$x)
y <- diabetes$y
maxk <- 15
eps <- 1e-9
x <- scale(x)
n <- nrow(x)
p <- ncol(x)
set.seed(42)
cost <- round(runif(p, 10, 100)) / 10
y <- y - mean(y)
mu <- rep(0, n)
beta <- rep(0, p)
r <- y - mu
cvec <- t(x) %*% r
per <- abs(cvec / cost)
which.max(per)
trace <- TRUE
x <- scale(x)
# variable setup
n <- nrow(x)
p <- ncol(x)
mu <- rep(0, n)
r <- y - mu
C <- t(x) %*% r / sd(r) / (n - 1)
Inactive <- rep(TRUE, p)
Active <- rep(FALSE, p)
Gram <- t(x) %*% x
# number lars iterations.  lars will do p steps.  number of lasso steps is
# not a fixed quantity
k <- 0
# number of variables currently in model (for lars, equivalent to step
# number)
nv <- 0

beta <- matrix(0, nrow = maxk + 1, ncol = p)
mul <- matrix(0, nrow = maxk + 1, ncol = n)

while (nv < p & k < maxk)
{
    k <- k + 1

    # 1. Find the next variable to add.
    # Calculate the projections (correlations) along the residuals (y - mu)
    # on each x. Find the largest of these projections and add that variable
    # to the to the active list.
    # TODO: is it possible that more than one variable could be added?

    # Equation 2.8
    r <- y - mu
    # cvec <- t(x) %*% r / sd(r) / (n - 1)
    cvec <- t(x) %*% r

    # score vector
    svec <- abs(cvec) / cost

    if (trace)
    {
        cat("\nIteration: ")
        cat(k)
        cat("\n")
        print(cbind(cvec, svec))
    }

    # TODO: TODO: WORK HERE
    # Equation 2.9
    # TODO: This part may need to change
    # TODO: this is important!
    # TODO: probably just "manually" add new to active
    # Remember (from selecting gamma) what the new variable will be
    smax <- max(svec)
    j <- svec >= smax - eps
    Active <- Active | j
    Inactive <- !Active
    nv <- nv + 1
    cmax <- abs(cvec[j])

    if (trace)
    {
        cat("selected: ")
        cat(colnames(x)[Active])
        cat("\n")
    }

    # 2. Find unit-vector of equal projection.
    # Following equations 2.4 through 2.6

    Signs <- sign(cvec[Active])
    # Equation 2.4
    XA <- x[ , Active] * rep(1, n) %*% t(Signs)
    # Equation 2.5
    gA <- t(XA) %*% XA
    one <- rep(1, sum(Active))
    # Equation 2.5
    AA <- 1 / sqrt(one %*% solve(gA) %*% one)
    # Equation 2.6 NOTE add the Signs to match package
    w <- AA %*% t(solve(gA) %*% one)
    # Equation 2.6
    u <- XA %*% t(w)

    # 3. Increment model fit in the direction of u.
    #  New estimate will be mu + \gamma * u where \gamma is large enough
    #  such that the next input variable will now be equally correlated.

    # TODO: What is the best gamma to choose?
    # TODO: Take the largest C which implies that the next step will cover
    # the most distance?  Then the most cost efficient is maybe the model
    # without the next predictor?

    if (nv == p)
    {
        # cheat and just use OLS
        beta[k + 1, ] <- coef(lm(y ~ x - 1))
    } else
    {
        # Equation 2.11
        a <- t(x) %*% u
        # Equation 2.13

        # 3. a. Each variable has "two" gammas, select shortest positive
        mingt0 <- function(x)
        {
            if (sum(x > 0) > 0)
            {
                return(min(x[x > 0]))
            }
            return(0)
        }

        # TODO: add in the "zero" distance?

        gammas <- apply(
            cbind((cmax - cvec[Inactive]) / (AA - a[Inactive]),
                  (cmax + cvec[Inactive]) / (AA + a[Inactive])), 1, mingt0)

        gammas <- gammas[which(gammas > 0)]

        # 3. b. Find what the cvecs would be after gamma
        # TODO: Need to double check that cmax - gammas * AA > 0?  Where do
        # "flips" happen?  Look into this.

        best <- which.max((cmax - gammas * AA) / cost[Inactive])

        # 3. c. Select the "best" gamma

        # TODO: WORK ON HERE DOWN
        # TODO: SHOULD BE SIMPLE CALC TO FIND CORRELATIONS
        # TODO: VERIFY THAT THIS WILL ACTUALLY WORK...
        # TODO:    ...as in that by skipping a "better" predictor we can still
        # add it in later without screwing ourselves over.

        gamma <- gammas[best]

        if (trace)
        {
            cat("gamma: ")
            cat(best)
            cat(" ")
            cat(gamma)
            cat("\nNew score: ")
            cat(cmax - gamma * AA)
            cat("\ncmax: ")
            cat(cmax)
            cat("\nnextmax: ")
            cat(max(cvec[Inactive]))
            cat("\n")
        }
        mu <- mu + gamma * u
        beta[k + 1, Active] <- beta[k, Active] + gamma * w * Signs
        mul[k + 1, ] <- mu
    }

}
