library(cosso)
library(lars)

data(diabetes)

y <- diabetes$y
x <- scale(diabetes$x)

lfit <- lars(x = x, y = y, type = "lar")
lfit$beta
lfit$mu

lfit.mu <- predict.lars(lfit, newx = x, type = "fit")$fit

mfit <- mlars(x = x, y = y)
mfit[[2]][ , 1]


sum(abs(lfit$beta - mfit))
apply(lfit$beta - mfit, 1, sum)


all.equal(lfit$beta, mfit)
