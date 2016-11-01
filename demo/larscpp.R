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

a <- larscpp(x, beta, y, maxk = 12)
out <- mlars(x, y, trace = TRUE, maxk=13)

AA: 17.85715
 a: 3.886778 1.851389 21 8.303722 5.245326 5.484568 -7.703031 8.689939 9.369332 8.16228
AA: 16.36749
 a: 5.628973 2.939797 17.85715 9.74125 9.449649 7.155949 -9.451004 12.73898 17.85715 10.
53716

1  age Out   6387.845  791.80273 1057.8302 0.000000
2  sex Out   1464.022  964.77590  936.5804 0.000000
3  bmi  In  19938.140    0.00000    0.0000 2.862822
4  map Out  15009.575  388.18982 1192.6033 0.000000
5   tc Out   7208.343  808.00129 1034.3359 0.000000
6  ldl Out   5917.476  903.65928  976.2522 0.000000
7  hdl Out -13422.051 1162.25328  490.0432 0.000000
8  tch Out  14634.544  430.83431 1164.4579 0.000000
9  ltg Out  19238.913   60.11927 1290.0203 0.000000
10 glu Out  13003.679  540.16299 1129.6037 0.000000

1  age Out   6154.174 1023.9839 1057.2117  0.00000
2  sex Out   1352.718 1161.2591  963.0428  0.00000
3  bmi  In  18675.636    0.0000    0.0000 17.23308
4  map Out  14510.361  513.2237 1202.4607  0.00000
5   tc Out   6892.998 1401.4426  936.3467  0.00000
6  ldl Out   5587.748 1223.0292  970.0269  0.00000
7  hdl Out -12958.950 1158.4298  680.0599  0.00000
8  tch Out  14112.111  891.6310 1071.6304  0.00000
9  ltg  In  18675.636    0.0000    0.0000 14.37025
10 glu Out  12512.969  841.8955 1098.4100  0.00000

1  age Out  3265.2517  678.4428 542.9948  0.000000
2  sex Out  -156.0552  805.3807 451.2286  0.000000
3  bmi  In  9510.9203    0.0000   0.0000 20.702760
4  map  In  9510.9203    0.0000   0.0000  3.773164
5   tc Out  2043.2145 1011.6144 455.7304  0.000000
6  ldl Out  1915.1455  795.6201 492.7580  0.000000
7  hdl Out -8108.4707  712.0030 175.5532  0.000000
8  tch Out  7574.1650  394.2242 614.0822  0.000000
9  ltg  In  9510.9203    0.0000   0.0000 17.853135
10 glu Out  7105.0463  461.4621 603.7475  0.000000

cvec <- t(x) %*% y
    Inactive <- rep(TRUE, p)
    Active <- rep(FALSE, p)
        cmax <- max(abs(cvec))
        j <- abs(cvec) >= cmax - eps
        Active <- Active | j
        Inactive <- !Active
        Signs <- sign(cvec[Active])
        # Equation 2.4
        XA <- x[ , Active] * rep(1, n) %*% t(Signs)
        # Equation 2.5
        gA <- t(XA) %*% XA
        one <- rep(1, sum(Active))
        # Equation 2.5
        AA <- 1/sqrt(one %*% solve(gA) %*% one)
        # Equation 2.6 NOTE add the Signs to match package
        w <- AA %*% t(solve(gA) %*% one)
        # Equation 2.6
        u <- XA %*% t(w)


n <- nrow(x)
p <- ncol(x)

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

