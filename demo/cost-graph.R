cost <- function(x, a = 100)
{
    if(x > a)
    {
        return(a * exp((x - a) / a))
    }
    return(x)
}

top <- 200
plot(0:top, sapply(0:top, cost), type = 'l', xlab = "True cost", ylab = "Effective cost")
