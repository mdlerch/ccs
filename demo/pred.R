library(cosso)
library(lars)
data(diabetes)
y <- diabetes$y
x <- scale(diabetes$x)


larsfit <- lars(x = x, y = y, type = "lar")
lassfit <- lars(x = x, y = y, type = "lasso")
mlarsfit <- mlars(x = x, y = y)
mlassfit <- mlasso(x = x, y = y)

predictlars(mlarsfit, x, 1.4)
predict.lars(larsfit, x, s = 1.4, type = "coef", mode = "step")$coef

# NOTE: Original lars worries about scaling
predictlasso(mlassfit, x, 1500 / 21)
predict.lars(lassfit, x, s = 1500, type = "coef", mode = "norm")$coef
