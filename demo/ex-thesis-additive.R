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

set.seed(98)
cost <- runif(p, 1, 10)
c4 <- clars4(x, y, cost, maxk = 9, trace = TRUE)

cls(x, y, cost, maxk = 50)
cls(x, y, maxk = 12)
mlars(x, y)$beta

set.seed(95)

cost <- rep(1, p)
out <- canalyze(x, y, cost, hold = 0.2)
cost <- runif(p, 1, 10)
out <- canalyze(x, y, cost, hold = 0.2)
c4 <- clars4(x, y, cost, maxk = 9, trace = TRUE)




# 1. All variables cost same
cost <- rep(1, p)
out <- canalyze(x, y, cost, hold = 0.2)


    n <- nrow(x)
    ntest <- round(hold * n)
    test <- sample(1:n, ntest, replace = FALSE)
    train <- setdiff(1:n, test)
    x <- x[train, ]
    y <- y[train]
cout <- clars(x, y, cost, maxk = 22, trace = TRUE)

c3po <- clarsB(x, y, maxk = 22, trace = TRUE)
mout <- mlars(x, y, maxk = 3, trace = TRUE)


plotcanalyze(out, "full")
plotcanalyze(out, "train")
plotcanalyze(out, "test")

# 2. Random cost same order of magnitude
cost <- runif(p, 1, 10)
out <- canalyze(x, y, cost, hold = 0.2)

cout <- clars(x, y, cost, maxk = 100, trace = TRUE)
c4 <- clars4(x, y, cost, maxk = 100, trace = TRUE)

plotcanalyze(out, "full")
plotcanalyze(out, "train")
plotcanalyze(out, "test")

# 2. Random cost differing of magnitude
cost <- exp(runif(p, 1, 10))
out <- canalyze(x, y, cost, hold = 0.2)

plotcanalyze(out, "full")
plotcanalyze(out, "train")
plotcanalyze(out, "test")

# 3. Best cost more. Best = largest lm coef
cost <- sort(runif(p, 1, 10))[best]
out <- canalyze(x, y, cost, hold = 0.2)

plotcanalyze(out, "full")
plotcanalyze(out, "train")
plotcanalyze(out, "test")

# 3. Best cost less. Best = largest lm coef
cost <- sort(runif(p, 1, 10))[worst]
out <- canalyze(x, y, cost, hold = 0.2)

plotcanalyze(out, "full")
plotcanalyze(out, "train")
plotcanalyze(out, "test")
