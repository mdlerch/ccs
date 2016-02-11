library(lars)
library(cosso)
data(diabetes)
x <- scale(diabetes$x)
y <- diabetes$y
maxk <- 50
eps <- 1e-9
x <- scale(x)
n <- nrow(x)
p <- ncol(x)
set.seed(42)
cost <- round(runif(p, 10, 100)) / 10
ntest <- round(0.30 * n)
test <- sample(1:n, ntest, replace = FALSE)
train <- setdiff(1:n, test)
xtrain <- x[train, ]
ytrain <- y[train]
xtest <- x[test, ]
ytest <- y[test]

LASSO <- FALSE
cout <- clars(xtrain, ytrain, cost = cost, trace = TRUE, maxk = 20)
mout <- mlars(xtrain, ytrain, trace = TRUE, maxk = 20)
nout <- clars(xtrain, ytrain, trace = TRUE, maxk = 20)
lout <- mlasso(xtrain, ytrain)

clarsscore <- evalclars(cout, xtrain, ytrain, cost)
mlarsscore <- evalclars(mout, xtrain, ytrain, cost)
lassscore <- evalclars(lout, xtrain, ytrain, cost)
ylims <- c(2500, 4000)
xlims <- c(min(clarsscore$modelcost, mlarsscore$modelcost),
           max(clarsscore$modelcost, mlarsscore$modelcost))
plot(clarsscore$score ~ clarsscore$modelcost, type = "p", ylim = ylims, xlim = xlims, col = "blue", pch = 16, cex = 2)
points(mlarsscore$score ~ mlarsscore$modelcost, col = "red", pch = 16, cex = 2)
points(lassscore$score ~ lassscore$modelcost, col = "purple", pch = 16, cex = 2)


clarsscore <- evalclars(cout, xtest, ytest, cost)
mlarsscore <- evalclars(mout, xtest, ytest, cost)
lassscore <- evalclars(lout, xtest, ytest, cost)
ylims <- c(2500, 4000)
xlims <- c(min(clarsscore$modelcost, mlarsscore$modelcost),
           max(clarsscore$modelcost, mlarsscore$modelcost))
plot(clarsscore$score ~ clarsscore$modelcost, type = "p", ylim = ylims, xlim = xlims, col = "blue", pch = 16, cex = 2)
points(mlarsscore$score ~ mlarsscore$modelcost, col = "red", pch = 16, cex = 2)
points(lassscore$score ~ lassscore$modelcost, col = "purple", pch = 16, cex = 2)


cout <- clars(x, y, trace = TRUE, maxk = 30)
mout <- mlars(x, y)
lout <- mlasso(x, y)
clarsscore <- evalclars(cout, x, y, cost)
mlarsscore <- evalclars(mout, x, y, cost)
lassscore <- evalclars(lout, x, y, cost)
ylims <- c(2500, 4000)
xlims <- c(min(clarsscore$modelcost, mlarsscore$modelcost),
           max(clarsscore$modelcost, mlarsscore$modelcost))
plot(clarsscore$score ~ clarsscore$modelcost, type = "p", ylim = ylims, xlim = xlims, col = "blue", pch = 16, cex = 2)
points(mlarsscore$score ~ mlarsscore$modelcost, col = "red", pch = 16, cex = 2)
points(lassscore$score ~ lassscore$modelcost, col = "purple", pch = 16, cex = 2)

plotlars(cout)
