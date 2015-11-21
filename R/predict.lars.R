predict.lars <- function(object, newx, s, type = c("fit", "coefficients"),
                         mode = c("step", "fraction", "norm","lambda"), ...)
{
    mode <- match.arg(mode)
    type <- match.arg(type)
    if(missing(newx) & type == "fit") {
        warning("Type=fit with no newx argument; type switched to coefficients")
        type <- "coefficients"
    }
    betas <- object$beta
    if(object$type!="LASSO"&&mode%in%c("fraction", "norm"))#detects discontinuities in norm
        betas=betabreaker(object)
    dimnames(betas) <- list(NULL,dimnames(betas)[[2]])
    sbetas <- scale(betas, FALSE, 1/object$normx)	#scaled for unit norm x
    kp <- dim(betas)
    k <- kp[1]
    p <- kp[2]
    steps <- seq(k)
    if(missing(s)) {
        s <- steps
        mode <- "step"
    }
    sbeta <- switch(mode,
                    step = {
                        if(any(s < 0) | any(s > k))
                            stop("Argument x out of range")
                        steps
                    }
                    ,
                    fraction = {
                        if(any(s > 1) | any(s < 0))
                            stop("Argument s out of range")
                        nbeta <- drop(abs(sbetas) %*% rep(1, p))
                        nbeta/nbeta[k]
                    }
                    ,
                    norm = {
                        nbeta <- drop(abs(sbetas) %*% rep(1, p))
                        if(any(s > nbeta[k]) | any(s < 0))
                            stop("Argument s out of range")
                        nbeta
                    }
                    ,
                    lambda={
                        lambdas=object$lambda
                        s[s>max(lambdas)]=max(lambdas)
                        s[s<0]=0
                        c(lambdas,0)
                    }
                    )

    sfrac <- (s - sbeta[1])/(sbeta[k] - sbeta[1])
    sbeta <- (sbeta - sbeta[1])/(sbeta[k] - sbeta[1])
    usbeta<-unique(sbeta)
    useq<-match(usbeta,sbeta)
    sbeta<-sbeta[useq]
    betas<-betas[useq,,drop=FALSE]
    coord <- approx(sbeta, seq(sbeta), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newbetas <- ((sbeta[right] - sfrac) * betas[left,  , drop = FALSE] +
                 (sfrac - sbeta[left]) * betas[right,  , drop = FALSE]) /
                 (sbeta[right] - sbeta[left])
    newbetas[left == right,  ] <- betas[left[left == right],  ]
    robject <- switch(type,
                      coefficients = list(s = s, fraction = sfrac, mode = mode, 
                                          coefficients = drop(newbetas)),
                      fit = list(s = s, fraction = sfrac, mode = mode,
                                 fit = drop(scale(newx, object$meanx, FALSE) %*% t(newbetas)) + object$mu))
    robject
}

# If a beta switches sign (ie goes from in and positive to in and negative)
# do something fancy, otherwise just return beta
betabreaker <- function(object)
{
    frac.step <- function(x){
        ##  step=x[1]-1
        ##  start=x[2]
        ##  end=x[3]
        ##  step -start/(end-start)
        x[1]-1 -x[2]/(x[3]-x[2])
    }
    beta <- object$beta
    sbeta <- sign(beta)
    kp <- dim(beta)
    k <- kp[1];p <- kp[2]
    dsbeta <- abs(sbeta[-1,]-sbeta[-k,])
    if(any(dsbeta==2)){
        bbeta <- matrix(cbind(
                              step.end=rep(1:(k-1),p),
                              beta.start=as.vector(beta[-k,]),
                              beta.end=as.vector(beta[-1,])
                              )[dsbeta==2],ncol=3)
        fsteps <- apply(bbeta,1,frac.step)
        new.beta <- predict(object,type="coefficient",s=fsteps+1,mode="step")$coef
        new.beta <- rbind(beta,new.beta)
        fo <- c(seq(k)-1,fsteps)
        beta <-  new.beta[order(fo),]
        dimnames(beta)[[1]] <- format(round(sort(fo),2))
    }
    beta

}
