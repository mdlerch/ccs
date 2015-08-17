n <- 20
m <- 10

set.seed(19)
X <- matrix(rnorm(n * m), m, n)
coeff <- runif(n, -5, 5)
y <- X %*% coeff + rnorm(m, 3)
costvec <- runif(n, 1, 10)

source("../R/brute.R")

out <- brute_cosso(y, X, costvec)

plot(out$scores ~ out$costs)



