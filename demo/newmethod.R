set.seed(42)
library(lars)
data(diabetes)
x <- scale(diabetes$x)
y <- diabetes$y
cost <- round(runif(p, 10, 100)) / 10

maxk <- 15
eps <- 1e-9
trace <- TRUE

mingt0 <- function(x)
{
    if (sum(x > 0) > 0)
    {
        return(min(x[x > 0]))
    }
    return(0)
}

n <- nrow(x)
p <- ncol(x)
y <- y - mean(y)
mu <- rep(0, n)
beta <- rep(0, p)
Inactive <- rep(TRUE, p)
Active <- rep(FALSE, p)
Gram <- t(x) %*% x
k <- 0
nv <- 0
beta <- matrix(0, nrow = maxk + 1, ncol = p)
mul <- matrix(0, nrow = maxk + 1, ncol = n)

r <- y - mu
cvec <- t(x) %*% r
svec <- abs(cvec) / cost
smax <- max(svec)
j <- svec >= smax - eps

while (nv < p & k < maxk)
{
    k <- k + 1
    r <- y - mu
    cvec <- t(x) %*% r
    svec <- abs(cvec) / cost

    if (trace) { cat("\nIteration: "); cat(k); cat("\n"); print(cbind(cvec, svec)) }

    Active <- Active | j
    Inactive <- !Active
    nv <- nv + 1
    cmax <- abs(cvec[j])

    if (trace) { cat("selected: "); cat(colnames(x)[Active]); cat("\n") }

    # 2. Find unit-vector of equal projection.
    Signs <- sign(cvec[Active])
    XA <- x[ , Active] * rep(1, n) %*% t(Signs)
    gA <- t(XA) %*% XA
    one <- rep(1, sum(Active))
    AA <- 1 / sqrt(one %*% solve(gA) %*% one)
    w <- AA %*% t(solve(gA) %*% one)
    u <- XA %*% t(w)

    # 3. Increment model fit in the direction of u.
    if (nv == p)
    {
        # cheat and just use OLS
        beta[k + 1, ] <- coef(lm(y ~ x - 1))
    } else
    {
        a <- t(x) %*% u
        # TODO: add in the "zero" distance?
        gammas <- apply(
            cbind((cmax - cvec[Inactive]) / (AA - a[Inactive]),
                  (cmax + cvec[Inactive]) / (AA + a[Inactive])), 1, mingt0)
        fgammas <- gammas[which(gammas > 0)]
        best <- which.max((cmax - fgammas * AA) / cost[Inactive])
        gamma <- fgammas[best]
        j <- rep(FALSE, p)
        j[which(cumsum(Inactive) == best)[1]] <- TRUE
        if (trace)
        {
            cat("gamma: "); cat(best); cat(" "); cat(gamma); cat("\nNew score: ")
            cat(cmax - gamma * AA); cat("\ncmax: "); cat(cmax); cat("\nnextmax: ")
            cat(max(cvec[Inactive])); cat("\n")
        }
        mu <- mu + gamma * u
        beta[k + 1, Active] <- beta[k, Active] + gamma * w * Signs
        mul[k + 1, ] <- mu
    }

}
