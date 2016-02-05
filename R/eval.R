evalclars <- function(obj, x, y, cost)
{
    betas <- obj$beta

    k <- nrow(betas)
    p <- ncol(betas)

    score <- numeric(k)
    modelcost <- numeric(k)
    for (i in 1:nrow(betas))
    {
        used <- rep(0, p)
        used[betas[i, ] != 0] <- 1
        ypred <- x %*% betas[i, ]
        score[i] <- var(y - ypred)
        modelcost[i] <- sum(used * cost)
    }

    return(data.frame(modelcost, score))
}
