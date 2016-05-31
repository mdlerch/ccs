# TODO: better strategy: normalize by saturated model cost?
library(lars)
library(cosso)
data(diabetes)
x <- scale(diabetes$x)
n <- nrow(x)
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
##                      all variables cost the same                      ##
###########################################################################

set.seed(95)
cost <- rep(1, p)

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

clarsT.eval <- evalclars(clarsP, xtest, ytest, cost)
larsT.eval <- evalclars(list(beta = larsP$beta), xtest, ytest, cost)
lassT.eval <- evalclars(list(beta = lassP$beta), xtest, ytest, cost)

xlim <- c(min(clarsT.eval$modelcost, larsT.eval$modelcost, lassT.eval$modelcost),
          max(clarsT.eval$modelcost, larsT.eval$modelcost, lassT.eval$modelcost))
ylim <- c(min(clarsT.eval$score, larsT.eval$score, lassT.eval$score),
          max(clarsT.eval$score, larsT.eval$score, lassT.eval$score))

plot(clarsT.eval$score ~ clarsT.eval$modelcost, type = "p", xlim = xlim, ylim = ylim,
     col = "blue", pch = 16, cex = 4, xlab = "Model Cost", ylab = "Score")

points(larsT.eval$score ~ larsT.eval$modelcost, col = "red", pch = 16, cex = 3)
points(lassT.eval$score ~ lassT.eval$modelcost, col = "purple", pch = 16, cex = 2)

legend('topright', c("clars", "lars", "lasso"), pch = 16, col = c("blue", "red", "purple"))

###########################################################################
##                  random cost same order of magnitude                  ##
###########################################################################

set.seed(95)
cost <- runif(p, 1, 10)

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


clarsT.eval <- evalclars(clarsP, xtest, ytest, cost)
larsT.eval <- evalclars(list(beta = larsP$beta), xtest, ytest, cost)
lassT.eval <- evalclars(list(beta = lassP$beta), xtest, ytest, cost)

xlim <- c(min(clarsT.eval$modelcost, larsT.eval$modelcost, lassT.eval$modelcost),
          max(clarsT.eval$modelcost, larsT.eval$modelcost, lassT.eval$modelcost))
ylim <- c(min(clarsT.eval$score, larsT.eval$score, lassT.eval$score),
          max(clarsT.eval$score, larsT.eval$score, lassT.eval$score))

plot(clarsT.eval$score ~ clarsT.eval$modelcost, type = "p", xlim = xlim, ylim = ylim,
     col = "blue", pch = 16, cex = 4, xlab = "Model Cost", ylab = "Score")

points(larsT.eval$score ~ larsT.eval$modelcost, col = "red", pch = 16, cex = 3)
points(lassT.eval$score ~ lassT.eval$modelcost, col = "purple", pch = 16, cex = 2)

legend('topright', c("clars", "lars", "lasso"), pch = 16, col = c("blue", "red", "purple"))

###########################################################################
##                random cost different order of magnitude                ##
###########################################################################

set.seed(95)
cost <- exp(runif(p, 1, 10))

clarsP <- clars(xtrain, ytrain, cost, maxk = 50, trace = TRUE)
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


clarsT.eval <- evalclars(clarsP, xtest, ytest, cost)
larsT.eval <- evalclars(list(beta = larsP$beta), xtest, ytest, cost)
lassT.eval <- evalclars(list(beta = lassP$beta), xtest, ytest, cost)

xlim <- c(min(clarsT.eval$modelcost, larsT.eval$modelcost, lassT.eval$modelcost),
          max(clarsT.eval$modelcost, larsT.eval$modelcost, lassT.eval$modelcost))
ylim <- c(min(clarsT.eval$score, larsT.eval$score, lassT.eval$score),
          max(clarsT.eval$score, larsT.eval$score, lassT.eval$score))

plot(clarsT.eval$score ~ clarsT.eval$modelcost, type = "p", xlim = xlim, ylim = ylim,
     col = "blue", pch = 16, cex = 4, xlab = "Model Cost", ylab = "Score")

points(larsT.eval$score ~ larsT.eval$modelcost, col = "red", pch = 16, cex = 3)
points(lassT.eval$score ~ lassT.eval$modelcost, col = "purple", pch = 16, cex = 2)

legend('topright', c("clars", "lars", "lasso"), pch = 16, col = c("blue", "red", "purple"))

###########################################################################
##                best cost more (best = largest ols coef                ##
###########################################################################

set.seed(95)
cost <- sort(runif(p, 1, 10))[best]

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


clarsT.eval <- evalclars(clarsP, xtest, ytest, cost)
larsT.eval <- evalclars(list(beta = larsP$beta), xtest, ytest, cost)
lassT.eval <- evalclars(list(beta = lassP$beta), xtest, ytest, cost)

xlim <- c(min(clarsT.eval$modelcost, larsT.eval$modelcost, lassT.eval$modelcost),
          max(clarsT.eval$modelcost, larsT.eval$modelcost, lassT.eval$modelcost))
ylim <- c(min(clarsT.eval$score, larsT.eval$score, lassT.eval$score),
          max(clarsT.eval$score, larsT.eval$score, lassT.eval$score))

plot(clarsT.eval$score ~ clarsT.eval$modelcost, type = "p", xlim = xlim, ylim = ylim,
     col = "blue", pch = 16, cex = 4, xlab = "Model Cost", ylab = "Score")

points(larsT.eval$score ~ larsT.eval$modelcost, col = "red", pch = 16, cex = 3)
points(lassT.eval$score ~ lassT.eval$modelcost, col = "purple", pch = 16, cex = 2)

legend('topright', c("clars", "lars", "lasso"), pch = 16, col = c("blue", "red", "purple"))

###########################################################################
##               best cost less (best = largest ols coef)                ##
###########################################################################

set.seed(95)
cost <- sort(runif(p, 1, 10))[worst]

clarsP <- clars(xtrain, ytrain, cost, maxk = 50, trace = FALSE)
mlasso(xtrain ,ytrain)
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


clarsT.eval <- evalclars(clarsP, xtest, ytest, cost)
larsT.eval <- evalclars(list(beta = larsP$beta), xtest, ytest, cost)
lassT.eval <- evalclars(list(beta = lassP$beta), xtest, ytest, cost)

xlim <- c(min(clarsT.eval$modelcost, larsT.eval$modelcost, lassT.eval$modelcost),
          max(clarsT.eval$modelcost, larsT.eval$modelcost, lassT.eval$modelcost))
ylim <- c(min(clarsT.eval$score, larsT.eval$score, lassT.eval$score),
          max(clarsT.eval$score, larsT.eval$score, lassT.eval$score))

plot(clarsT.eval$score ~ clarsT.eval$modelcost, type = "p", xlim = xlim, ylim = ylim,
     col = "blue", pch = 16, cex = 4, xlab = "Model Cost", ylab = "Score")

points(larsT.eval$score ~ larsT.eval$modelcost, col = "red", pch = 16, cex = 3)
points(lassT.eval$score ~ lassT.eval$modelcost, col = "purple", pch = 16, cex = 2)

legend('topright', c("clars", "lars", "lasso"), pch = 16, col = c("blue", "red", "purple"))
