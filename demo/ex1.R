library(lqa)
library(parcor)
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
mfit <- mfit[[1]]


zapsmall(lfit$beta - mfit)


sum(abs(lfit$beta - mfit))
apply(lfit$beta - mfit, 1, sum)


all.equal(lfit$beta, mfit)

lafit <- lars(x = x, y = y, type = "lasso", trace = TRUE)
lafit$beta

mafit <- mlasso(x = x, y = y)
mafit

lafit$beta - mafit

mwlasso(x = x, y = y)

adafit <- adalasso(X = x, y = y, intercept = FALSE)
