library(lars)
data(diabetes)
x <- scale(diabetes$x)
n <- nrow(x)
p <- ncol(x)
y <- diabetes$y
maxk <- 50
eps <- 1e-9
trace <- TRUE

costfunc <- function(A) { sum(cost[A]) }

set.seed(95)
cost <- rep(1, p)
