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
cost <- round(runif(p, 1, 10))
price <- cost / max(cost)

y <- y - mean(y)

mu <- rep(0, n)
beta <- rep(0, p)


r <- y - mu
cvec <- t(x) %*% r
per <- abs(cvec / price)
which.max(per)



