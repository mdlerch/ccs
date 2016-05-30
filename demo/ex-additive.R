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

clars1 <- clars(xtrain, ytrain, cost, maxk = 50, trace = TRUE)
mlasso(xtrain ,ytrain)
lars1 <- lars(xtrain, ytrain, type = "lar")
lass1 <- lars(xtrain, ytrain, type = "lasso")

clars1.eval <- evalclars(clars1, xtrain, ytrain, cost)
lars1.eval <- evalclars(list(beta = lars1$beta), xtrain, ytrain, cost)
lass1.eval <- evalclars(list(beta = lass1$beta), xtrain, ytrain, cost)

xlim <- c(min(clars1.eval$modelcost, lars1.eval$modelcost, lass1.eval$modelcost),
          max(clars1.eval$modelcost, lars1.eval$modelcost, lass1.eval$modelcost))
ylim <- c(min(clars1.eval$score, lars1.eval$score, lass1.eval$score),
          max(clars1.eval$score, lars1.eval$score, lass1.eval$score))

plot(clars1.eval$score ~ clars1.eval$modelcost, type = "p", xlim = xlim, ylim = ylim,
     col = "blue", pch = 16, cex = 4, xlab = "Model Cost", ylab = "Score")

points(lars1.eval$score ~ lars1.eval$modelcost, col = "red", pch = 16, cex = 3)
points(lass1.eval$score ~ lass1.eval$modelcost, col = "purple", pch = 16, cex = 2)

legend('topright', c("clars", "lars", "lasso"), pch = 16, col = c("blue", "red", "purple"))

###########################################################################
##                  random cost same order of magnitude                  ##
###########################################################################

set.seed(95)
cost <- runif(p, 1, 10)
out <- canalyze(x, y, cost, hold = 0.2)

cout <- clars(x, y, cost, maxk = 100, trace = TRUE)
c4 <- clars4(x, y, cost, maxk = 100, trace = TRUE)

plotcanalyze(out, "full")
plotcanalyze(out, "train")
plotcanalyze(out, "test")

###########################################################################
##                random cost different order of magnitude                ##
###########################################################################

cost <- exp(runif(p, 1, 10))
out <- canalyze(x, y, cost, hold = 0.2)

plotcanalyze(out, "full")
plotcanalyze(out, "train")
plotcanalyze(out, "test")

###########################################################################
##                best cost more (best = largest ols coef                ##
###########################################################################

cost <- sort(runif(p, 1, 10))[best]
out <- canalyze(x, y, cost, hold = 0.2)

plotcanalyze(out, "full")
plotcanalyze(out, "train")
plotcanalyze(out, "test")

###########################################################################
##               best cost less (best = largest ols coef)                ##
###########################################################################

cost <- sort(runif(p, 1, 10))[worst]
out <- canalyze(x, y, cost, hold = 0.2)

plotcanalyze(out, "full")
plotcanalyze(out, "train")
plotcanalyze(out, "test")
