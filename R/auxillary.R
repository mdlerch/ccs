
mingt0 <- function(x)
{
    if (sum(x > 0) > 0)
    {
        return(min(x[x > 0]))
    }
    return(0)
}

minlt0 <- function(x)
{
    if (sum(x < 0) > 0)
    {
        return(min(x[x < 0]))
    }
    return(0)
}



