
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


mindist <- function(x)
{
    if (sum(x > 0) != 0)
    {
        temp <- x[x != 0]
        return(temp[which(abs(temp) == min(abs(temp)))])
    }
    return(0)
}
