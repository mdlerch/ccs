library(cosso)
library(lars)
data(diabetes)
y <- diabetes$y
x <- scale(diabetes$x)

larsfit <- lars(x = x, y = y, type = "lar")
lassfit <- lars(x = x, y = y, type = "lasso")
mlarsfit <- mlars(x = x, y = y)
mlassfit <- mlasso(x = x, y = y)

sum(abs(larsfit$beta - mlarsfit$beta) > 1e-10)
sum(abs(lassfit$beta - mlassfit$beta) > 1e-10)
