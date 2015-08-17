brute_cosso <- function(y, x, cost, savefunc = function(x) {summary(x)[["sigma"]]})
{
    m <- ncol(x)
    n <- nrow(x)

    IN <- matrix(0, 2 ^ m - 1, m)
    modelcost <- numeric(2 ^ m - 1)
    modelscore <- numeric(2 ^ m - 1)
    for (i in 1:(2 ^ m - 1))
    {
        contr <- as.integer(intToBits(i))[1:m]
        IN[i, ] <- contr
        modelcost[i] <- contr %*% cost
        X <- x[ , contr == 1]
        modelscore[i] <- savefunc(lm(y ~ X))
    }

    list(vars = IN, scores = modelscore, costs = modelcost)
}


