library(lqa)
library(parcor)
library(cosso)
library(lars)
data(diabetes)
y <- diabetes$y
x <- scale(diabetes$x)

lfit <- lars(x = x, y = y, type = "lar")
lafit <- lars(x = x, y = y, type = "lasso")

malarsfit <- mlars(x = x, y = y)
malassofit <- mlasso(x = x, y = y)

objlars <- list(betas = malarsfit)
objlasso <- list(betas = malassofit)

s <- 20
predictlars(objlasso, x, s)
predict.lars(lfit, type = "coefficients", mode = "norm", s = s)$coefficients

lafit$beta - mafit

mwlasso(x = x, y = y)

adafit <- adalasso(X = x, y = y, intercept = FALSE)
