library(lars)
library(cosso)
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


clars(x, y, trace = TRUE)
mlars(x, y, trace = TRUE)

-
out <- mlars(x, y)$beta

mu <- rep(0, n)

r <- y - mu

r_lm <- y - x %*% solve(t(x) %*% x) %*% t(x) %*% y

C <- t(x) %*% r / sd(r) / (n - 1)

CM <- max(abs(C))
p <- price * CM

S <- abs(C) - p
cbind(S, C, p, price, cost)

idx <- which.max(S)
SM <- S[idx]

u <- x[ , idx]

for (gam in seq(0, 10, 1))
{
    rp <- r - u * gam
    cp <- t(x) %*% rp / sd(rp) / (n - 1)
    Sp <- cp - p
    print(which.max(Sp))
}





