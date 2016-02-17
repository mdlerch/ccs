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

cout <- clars(x, y, cost = cost, trace = TRUE, maxk = 100)
mout <- mlars(x, y)
plotlars(cout2)

clarsscore <- evalclars(cout2, x, y, cost)

mlarsscore <- evalclars(mout, x, y, cost)

plot(clarsscore$score ~ clarsscore$modelcost, type = "p", ylim = c(2500, 3600), col = "blue", pch = 16, cex = 2)
points(mlarsscore$score ~ mlarsscore$modelcost, col = "red", pch = 16, cex = 2)

