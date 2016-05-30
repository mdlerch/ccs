# TODO: better strategy: normalize by saturated model cost?
library(lars)
library(cosso)
data(diabetes)
x <- scale(diabetes$x)
n <- nrow(x)
p <- ncol(x)
x2e <- matrix(rnorm(n * p, 0, 0.5), nrow = n, ncol = p)
x2 <- x + x2e
set.seed(95)
idx <- sample(1:n, round(n * 0.2), replace = FALSE)
x2[ , 2] <- x[ , 2]
x2[idx, 2] <- ifelse(x[idx, 2] == unique(x[ , 2])[1], unique(x[ , 2])[2],
                     unique(x[ , 2])[1])
x <- cbind(x, x2)
x <- scale(x)
p <- ncol(x)

y <- diabetes$y
maxk <- 50
eps <- 1e-9
best <- order(abs(lm(y ~ x - 1)$coef))
worst <- order(-abs(lm(y ~ x - 1)$coef))

set.seed(95)
train <- sample(1:n, round(n * 0.8), replace = FALSE)
test <- setdiff(1:n, train)
xtrain <- x[train, ]
ytrain <- y[train]
xtrain <- scale(xtrain)
ytrain <- ytrain - mean(ytrain)
xtest <- x[test, ]
ytest <- y[test]


###########################################################################
##                         half cost duplicates                          ##
###########################################################################


set.seed(95)
cost <- rep(1, p / 2)
cost <- c(cost, cost / 2)

clarsP <- clars(xtrain, ytrain, cost, maxk = 50, trace = FALSE)
larsP <- lars(xtrain, ytrain, type = "lar")
lassP <- lars(xtrain, ytrain, type = "lasso")

clarsP.eval <- evalclars(clarsP, xtrain, ytrain, cost)
larsP.eval <- evalclars(list(beta = larsP$beta), xtrain, ytrain, cost)
lassP.eval <- evalclars(list(beta = lassP$beta), xtrain, ytrain, cost)

xlim <- c(min(clarsP.eval$modelcost, larsP.eval$modelcost, lassP.eval$modelcost),
          max(clarsP.eval$modelcost, larsP.eval$modelcost, lassP.eval$modelcost))
ylim <- c(min(clarsP.eval$score, larsP.eval$score, lassP.eval$score),
          max(clarsP.eval$score, larsP.eval$score, lassP.eval$score))

plot(clarsP.eval$score ~ clarsP.eval$modelcost, type = "p", xlim = xlim, ylim = ylim,
     col = "blue", pch = 16, cex = 4, xlab = "Model Cost", ylab = "Score")

points(larsP.eval$score ~ larsP.eval$modelcost, col = "red", pch = 16, cex = 3)
points(lassP.eval$score ~ lassP.eval$modelcost, col = "purple", pch = 16, cex = 2)

legend('topright', c("clars", "lars", "lasso"), pch = 16, col = c("blue", "red", "purple"))

###########################################################################
##                        double cost duplicates                         ##
###########################################################################


set.seed(95)
cost <- rep(1, p / 2)
cost <- c(cost, cost * 2)

clarsP <- clars(xtrain, ytrain, cost, maxk = 50, trace = FALSE)
larsP <- lars(xtrain, ytrain, type = "lar")
lassP <- lars(xtrain, ytrain, type = "lasso")

clarsP.eval <- evalclars(clarsP, xtrain, ytrain, cost)
larsP.eval <- evalclars(list(beta = larsP$beta), xtrain, ytrain, cost)
lassP.eval <- evalclars(list(beta = lassP$beta), xtrain, ytrain, cost)

xlim <- c(min(clarsP.eval$modelcost, larsP.eval$modelcost, lassP.eval$modelcost),
          max(clarsP.eval$modelcost, larsP.eval$modelcost, lassP.eval$modelcost))
ylim <- c(min(clarsP.eval$score, larsP.eval$score, lassP.eval$score),
          max(clarsP.eval$score, larsP.eval$score, lassP.eval$score))

plot(clarsP.eval$score ~ clarsP.eval$modelcost, type = "p", xlim = xlim, ylim = ylim,
     col = "blue", pch = 16, cex = 4, xlab = "Model Cost", ylab = "Score")

points(larsP.eval$score ~ larsP.eval$modelcost, col = "red", pch = 16, cex = 3)
points(lassP.eval$score ~ lassP.eval$modelcost, col = "purple", pch = 16, cex = 2)

legend('topright', c("clars", "lars", "lasso"), pch = 16, col = c("blue", "red", "purple"))

