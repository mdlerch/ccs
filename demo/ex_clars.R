library(lars)
data(diabetes)
x <- scale(diabetes$x)
y <- diabetes$y
maxk <- 15
eps <- 1e-9
x <- scale(x)
n <- nrow(x)
p <- ncol(x)

cost <- round(runif(p, 1, 10))

mu <- rep(0, n)

r <- y - mu

C <- t(x) %*% r / sd(r) / (n - 1)

CM <- max(abs(C))

S <- abs(C) - cost * CM / max(cost)
cbind(S, C, cost)

idx <- which.max(S)
SM <- S[idx]

u <- x[ , idx]

