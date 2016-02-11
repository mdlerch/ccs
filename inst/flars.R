# LARS (slightly adapted) from the R lars package

# x <- as.matrix(read.csv("./X.csv")[-1])
# y <- read.csv("./y.csv")$V1

# x <- scale(diabetes$x)
# y <- diabetes$y

flars <- function(x, y, eps = .Machine$double.eps, max.steps = 100)
  {
      nm <- dim(x)
      n <- nm[1]
      m <- nm[2]
      im <- inactive <- seq(m)
      one <- rep(1, n)
      vn <- dimnames(x)[[2]]
      ### Center x and y, and scale x, and save the means and sds
      meanx <- drop(one %*% x)/n
      x <- scale(x, meanx, FALSE)	# centers x
      normx <- sqrt(drop(one %*% (x^2)))
      nosignal<-normx/sqrt(n) < eps
      if(any(nosignal))# ignore variables with too small a variance
      {
          ignores<-im[nosignal]
          inactive<-im[-ignores]
          normx[nosignal]<-eps*sqrt(n)
      } else
      {
          ignores <- NULL #singularities; augmented later as well
      }
      names(normx) <- NULL
      x <- scale(x, FALSE, normx)	# scales x
          Gram <- t(x) %*% x	#Time saving
      mu <- mean(y)
      y <- drop(y - mu)
      Cvec <- drop(t(y) %*% x)
      ssy <- sum(y^2)	### Some initializations
      residuals <- y
      beta <- matrix(0, max.steps + 1, m)	# beta starts at 0
      Gamrat <- NULL
      arc.length <- NULL
      R2 <- 1
      RSS <- ssy
      first.in <- integer(m)
      active <- NULL	# maintains active set
      actions <- as.list(seq(max.steps))
      Sign <- NULL	# Keeps the sign of the terms in the model
      R <- NULL	###
      ### Now the main loop over moves
      k <- 0
      while((k < max.steps) & (length(active) < min(m - length(ignores),n-1)) )
      {
          action <- NULL
          k <- k + 1
          C <- Cvec[inactive]	#
          ### identify the largest nonactive gradient
          Cmax <- max(abs(C))	### Check if we are in a DROP situation
          if(!any(drops)) {
              new <- abs(C) >= Cmax - eps
              C <- C[!new]	# for later
              new <- inactive[new]	# Get index numbers
              ### We keep the choleski R  of X[,active] (in the order they enter)
              for(inew in new) {
                  if(use.Gram) {
                      R <- updateR(Gram[inew, inew], R, drop(Gram[
                                                             inew, active]), Gram = TRUE,eps=eps)
                  }
                  else {
                      R <- updateR(x[, inew], R, x[, active], Gram
                                   = FALSE,eps=eps)
                  }
                  if(attr(R, "rank") == length(active)) {
                      ##singularity; back out
                      nR <- seq(length(active))
                      R <- R[nR, nR, drop = FALSE]
                      attr(R, "rank") <- length(active)
                      ignores <- c(ignores, inew)
                      action <- c(action,  - inew)
                      if(trace)
                          cat("LARS Step", k, ":\t Variable", inew, 
                              "\tcollinear; dropped for good\n")	#
                  }
                  else {
                      if(first.in[inew] == 0)
                          first.in[inew] <- k
                      active <- c(active, inew)
                      Sign <- c(Sign, sign(Cvec[inew]))
                      action <- c(action, inew)
                      if(trace)
                          cat("LARS Step", k, ":\t Variable", inew, 
                              "\tadded\n")	#
                  }
              }
          }
          Gi1 <- backsolve(R, backsolvet(R, Sign))
          ### Now we have to do the forward.stagewise dance
          ### This is equivalent to NNLS
          dropouts<-NULL
          A <- 1/sqrt(sum(Gi1 * Sign))
          w <- A * Gi1	# note that w has the right signs
          ### Now we see how far we go along this direction before the
          ### next competitor arrives. There are several cases
          ###
          ### If the active set is all of x, go all the way
          if(length(active) >=  min(n-1, m - length(ignores) ) ) {
              gamhat <- Cmax/A
          }
          else {
              a <- drop(w %*% Gram[active,  - c(active,ignores), drop = FALSE])
              gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
              ### Any dropouts will have gam=0, which are ignored here
              gamhat <- min(gam[gam > eps], Cmax/A)
          }
          beta[k + 1,  ] <- beta[k,  ]
          beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
          Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% w
          Gamrat <- c(Gamrat, gamhat/(Cmax/A))
          arc.length <- c(arc.length, gamhat)
          ### Check if we have to drop any guys
          if(!is.null(vn))
              names(action) <- vn[abs(action)]
          actions[[k]] <- action
          inactive <- im[ - c(active, ignores)]
      }
      beta <- beta[seq(k + 1),  ]	#
      dimnames(beta) <- list(paste(0:k), vn)	### Now compute RSS and R2
      residuals <- y - x %*% t(beta)
      beta <- scale(beta, FALSE, normx)
      RSS <- apply(residuals^2, 2, sum)
      R2 <- 1 - RSS/RSS[1]
      Cp <- ((n - k - 1) * RSS)/rev(RSS)[1] - n + 2 * seq(k + 1)
      object <- list(call = call, type = TYPE, R2 = R2, RSS = RSS, Cp = Cp,
                     actions = actions[seq(k)], entry = first.in, Gamrat = Gamrat,
                     arc.length = arc.length, Gram = if(use.Gram) Gram else NULL,
                     beta = beta, mu = mu, normx = normx, meanx = meanx)
      class(object) <- "lars"
      object
  }
