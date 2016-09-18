evalclars <- function(obj, x, y, cost, costfunc = NULL)
{
    if (is.null(costfunc))
    {
        costfunc <- function(Active)
        {
            sum(cost[Active])
        }
    }

    betas <- obj$beta

    k <- nrow(betas)
    p <- ncol(betas)
    n <- length(y)

    score <- numeric(k)
    modelcost <- numeric(k)
    for (i in 1:nrow(betas))
    {
        used <- betas[i, ] != 0
        ypred <- x %*% betas[i, ]
        score[i] <- var(y - ypred) * (n - 1) / n
        modelcost[i] <- costfunc(used)
    }

    return(data.frame(modelcost, score))
}
