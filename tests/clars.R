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


###########################################################################
##              clars == LARS == mlars for cost = constant               ##
###########################################################################

a <- clars(x, y, cost = rep(1, p))
b <- lars(x, y, type = "lar")
d <- mlars(x, y)

all.equal(a$beta, b$beta)


eps <- 1e-6
if (any(!(abs(a$beta - b$beta) < eps)))
{
    stop("LARS and clars don't match")
}

if (any(!(abs(a$beta - d$beta) < eps)))
{
    stop("mlars and clars don't match")
}


###########################################################################
##                    clars with same order magnitude                    ##
###########################################################################

set.seed(42)
cost <- round(runif(p, 10, 100)) / 10

cout <- clars(x, y, cost = cost, trace = TRUE, maxk = 50)
mout <- clars(x, y, trace = TRUE, maxk = 50)
cout2 <- list(beta = cout$beta2)

plotlars(cout)
plotlars(cout2)
plotlars(mout)


clarsscore <- evalclars(cout, x, y, cost)
clar2score <- evalclars(cout2, x, y, cost)
mlarsscore <- evalclars(mout, x, y, cost)

plot(clarsscore$score ~ clarsscore$modelcost, type = "p", ylim = c(2500, 3600), col = "blue", pch = 16, cex = 2)
points(mlarsscore$score ~ mlarsscore$modelcost, col = "red", pch = 16, cex = 2)
points(clar2score$score ~ clar2score$modelcost, col = "green", pch = 16, cex = 2)

