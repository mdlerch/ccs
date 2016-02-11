library(lars)
library(cosso)
data(diabetes)
x <- scale(diabetes$x)
n <- nrow(x)
p <- ncol(x)
y <- diabetes$y
maxk <- 50
eps <- 1e-9

set.seed(95)

# 1. All variables cost same 
cost <- rep(1, p)
canalyze(x, y, cost, hold = 0.2)


