canalyze <- function(x, y, cost, hold = 0.2, maxk = 500, eps = 1e-5)
{
    result <- list()

    x <- scale(x)

    n <- nrow(x)
    ntest <- round(hold * n)
    test <- sample(1:n, ntest, replace = FALSE)
    train <- setdiff(1:n, test)

    xtrain <- x[train, ]
    ytrain <- y[train]

    xtest <- x[test, ]
    ytest <- y[test]

    cout <- clars(x, y, cost, maxk = maxk)
    cou2 <- clars2(x, y, cost, maxk = maxk)
    mout <- mlars(x, y, maxk = maxk)
    lout <- mlasso(x, y)

    result$clars_full <- evalclars(cout, x, y, cost)
    result$clar2_full <- evalclars(cou2, x, y, cost)
    result$lars_full <- evalclars(mout, x, y, cost)
    result$lass_full <- evalclars(lout, x, y, cost)

    cout <- clars(xtrain, ytrain, cost, maxk = maxk)
    cou2 <- clars2(xtrain, ytrain, cost, maxk = maxk)
    mout <- mlars(xtrain, ytrain, maxk = maxk)
    lout <- mlasso(xtrain, ytrain)

    result$clars_train <- evalclars(cout, xtrain, ytrain, cost)
    result$clar2_train <- evalclars(cou2, xtrain, ytrain, cost)
    result$lars_train <- evalclars(mout, xtrain, ytrain, cost)
    result$lass_train <- evalclars(lout, xtrain, ytrain, cost)

    result$clars_test <- evalclars(cout, xtest, ytest, cost)
    result$clar2_test <- evalclars(cou2, xtest, ytest, cost)
    result$lars_test <- evalclars(mout, xtest, ytest, cost)
    result$lass_test <- evalclars(lout, xtest, ytest, cost)

    return(result)
}

