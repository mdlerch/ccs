source("util.R")

X <- as.matrix(read.csv("./X.csv")[-1])
y <- read.csv("./y.csv")$V1

eps <- 1e-6

# variable setup
n <- nrow(X)
p <- ncol(X)

maxk <- 512 * p

# current position
mu <- rep(0, n)
Inactive <- rep(TRUE, p)
Active <- rep(FALSE, p)

# USE GRAM
X <- X / sqrt(n - 1) # in order to match lars
Gram <- t(X) %*% X
lassocond <- FALSE

k <- 0
nv <- 0
loop <- TRUE

# X <- scale(X)

beta <- matrix(0, nrow = maxk, ncol = p)

while (nv < p & k < maxk & loop)
{
    k <- k + 1
    cvec <- t(X) %*% (y - mu)
    cmax <- max(abs(cvec))

    j <- abs(cvec) >= cmax - eps

    if (!lassocond)
    {
        Active <- Active | j
        Inactive <- !Active
        nv <- nv + 1
    }

    Signs <- sign(cvec[Active])

    # TODO: watch out for really big stuff
    R <- chol(Gram[Active, Active])

    # TODO: figure out what is going on here!!!
    GA1 <- backsolve(R, backsolvet(R, Signs))
    AA <- 1 / sqrt(sum(GA1 * Signs))
    w <- drop(AA %*% GA1)

    # Equiangular vector
    u <- X[ , Active] %*% w

    # If all variables active go to lsq solution
    if (nv == p)
    {
        gamma <- cvec / AA
    } else
    {
        # eqn 2.11
        a <- t(X) %*% u
        # eqn 2.13
        temp <- c( (cmax - cvec[Inactive]) / (AA - a[Inactive]),
                  (cmax + cvec[Inactive]) / (AA + a[Inactive]) )
        gamma <- min(temp[temp > eps], cmax / AA)
    }

    # TODO: there is a lasso step here

    # eqn 2.12
    mu = mu + gamma * u

    # TODO: what is Cardi?  assume "empty" Assume stop == 0

    beta[k + 1, Active] <- beta[k, Active] + gamma * w


    # TODO: there is a lasso step here
}

