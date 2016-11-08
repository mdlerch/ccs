library(lars)
library(ccs)
data(diabetes)
x <- scale(diabetes$x)
y <- diabetes$y
maxk <- 15
eps <- 1e-9
x <- scale(x)
n <- nrow(x)
p <- ncol(x)
beta = matrix(0, nrow = maxk, ncol = p)

a <- lassocpp(x, beta, y, maxk = 20)
larscpp(x, beta, y, maxk = 15)
b <- mlasso(x, y, trace = TRUE, maxk=13)

i <- 10
rbind(a[i, ], b$beta[i, ])
i <- 11
rbind(a[i, ], b$beta[i, ])
i <- 12
rbind(a[i, ], b$beta[i, ])
i <- 13
rbind(a[i, ], b$beta[i, ])

60.11927
513.2237
175.5532
259.3675
88.65915
43.67793
135.9841
54.0156
5.567232
41.99966
7.270701
27.97002



