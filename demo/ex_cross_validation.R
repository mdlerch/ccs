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

ntest <- round(0.10 * n)
test <- sample(1:n, ntest, replace = FALSE)
train <- setdiff(1:n, test)

xtrain <- x[train, ]
ytrain <- y[train]
xtest <- x[test, ]
ytest <- y[test]

cout <- clars(xtrain, ytrain, cost = cost, trace = TRUE, type = 3, maxk = 100)
mout <- mlars(xtrain, ytrain)

clarsscore <- evalclars(cout, xtest, ytest, cost)
mlarsscore <- evalclars(mout, xtest, ytest, cost)

plot(clarsscore$score ~ clarsscore$modelcost, type = "p", ylim = c(2500, 3600), col = "blue", pch = 16, cex = 2)
points(mlarsscore$score ~ mlarsscore$modelcost, col = "red", pch = 16, cex = 2)

