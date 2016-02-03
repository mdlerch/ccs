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

cout2 <- clars(x, y, cost = cost, trace = TRUE, type = 2, maxk = 50)
plotlars(cout2)

cout3 <- clars(x, y, cost = cost, trace = TRUE, type = 3, maxk = 50)

cout1 <- clars(x, y, cost = cost, trace = TRUE, type = 1)

cout2 <- clars2(x, y, cost = cost, trace = TRUE)
cout3 <- clars2(x, y, trace = TRUE)
mout <- mlars(x, y, trace = TRUE)
lout <- lars(x, y, type = "lar")
louts <- lars(x, y, type = "lasso")

idx <- 4; sd(y - predictlars(cout, x, idx)); sd(y - predictlars(mout, x, idx))

