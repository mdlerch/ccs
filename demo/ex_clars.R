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


cout <- clars(x, y)
cout2 <- clars2(x, y)
mout <- mlars(x, y)

idx <- 4; sd(y - predictlars(cout, x, idx)); sd(y - predictlars(mout, x, idx))

