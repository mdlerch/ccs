source("util.R")

X <- as.matrix(read.csv("./X.csv")[-1])
y <- read.csv("./y.csv")$V1
