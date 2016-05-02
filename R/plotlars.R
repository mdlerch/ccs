plotlars <- function(lars.object)
{
    betas <- lars.object$beta
    k <- nrow(betas)
    p <- ncol(betas)

    plot(0, 0, type = 'n', xlim = c(0, k), ylim = c(min(betas), max(betas)))
    for (v in 1:p)
    {
        lines(1:k, betas[ , v], col = v)
        lines(1:k, betas[ , v], col = v, pch = 16)
    }
}


plotlars2 <- function(lars.object, ylim = c(-30, 30))
{
    betas <- lars.object$beta
    k <- nrow(betas)
    p <- ncol(betas)

    plot(0, 0, type = 'n', xlim = c(0, k), ylim = ylim)
    for (v in 1:p)
    {
        lines(1:k, betas[ , v], col = v)
        lines(1:k, betas[ , v], col = v, pch = 16)
    }
}

plotclars <- function(lars.object, cost)
{
    eps <- 1e-5
    betas <- lars.object$beta
    k <- nrow(betas)
    p <- ncol(betas)

    plot(0, 0, type = 'n', xlim = c(0, k), ylim = c(min(betas), max(betas)))
    for (v in 1:p)
    {
        lines(1:k, betas[ , v], col = v)
        lines(1:k, betas[ , v], col = v, pch = 16)
    }
    par(new = TRUE)
    plot(0, 0, type = 'n', xlim = c(0, k), ylim = c(0, sum(cost)), xaxt = 'n', yaxt = 'n')
    axis(4)
    mtext("y2",side=4,line=3)
    for (v in 1:k)
    {
        points(v, sum(cost[abs(betas[v, ]) > eps]))
    }
}

plotcanalyze <- function(obj, type = "full", bonus = FALSE)
{
    out <- paste(c("clars", "lars", "lass", "clar2"), "_", type, sep = "")
    xlims <- c(min(obj[[out[1]]]$modelcost, obj[[out[2]]]$modelcost,
                   obj[[out[3]]]$modelcost, obj[[out[4]]]$modelcost),
               max(obj[[out[1]]]$modelcost, obj[[out[2]]]$modelcost,
                   obj[[out[3]]]$modelcost, obj[[out[4]]]$modelcostt))
    ylims <- c(min(obj[[out[1]]]$score, obj[[out[2]]]$score,
                   obj[[out[3]]]$score, obj[[out[4]]]$score),
               max(obj[[out[1]]]$score, obj[[out[2]]]$score,
                   obj[[out[3]]]$score, obj[[out[4]]]$score))
    plot(obj[[out[1]]]$score ~ obj[[out[1]]]$modelcost, type = "p",
         ylim = ylims, xlim = xlims, col = "blue", pch = 16, cex = 4,
         xlab = "Model Cost", ylab = "Score")
    points(obj[[out[2]]]$score ~ obj[[out[2]]]$modelcost, col = "red", pch = 16, cex = 3)
    points(obj[[out[3]]]$score ~ obj[[out[3]]]$modelcost, col = "purple", pch = 16, cex = 2)
    if (bonus)
    {
        points(obj[[out[4]]]$score ~ obj[[out[4]]]$modelcost, col = "green", pch = 16, cex = 1)
    }

    if (!bonus)
    {
        legend('topright', c("clars", "lars", "lasso"), pch = 16, col = c("blue", "red", "purple"))
    } else {
        legend('topright', c("clars", "lars", "lasso", "clar2"), pch = 16, col = c("blue", "red", "purple", "green"))
    }
}

