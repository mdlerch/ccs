lasso.adapt.bic2<-function(x,y){

    # adaptive lasso from lars with BIC stopping rule 
    # this one uses the "known variance" version of BIC with RSS/(full model mse)
    # must use a recent version of R so that normalize=FALSE can be used in lars

    require(lars)
    ok <- complete.cases(x,y)
    x <- x[ok,]                            # get rid of na's
    y <- y[ok]                             # since regsubsets can't handle na's
    m <- ncol(x)
    n <- nrow(x)
    x <- as.matrix(x)                      # in case x is not a matrix

    #  standardize variables like lars does 
    one <- rep(1, n)
    meanx <- drop(one %*% x)/n
    xc <- scale(x, meanx, FALSE)         # first subtracts mean
    normx <- sqrt(drop(one %*% (xc^2)))
    names(normx) <- NULL
    xs <- scale(xc, FALSE, normx)        # now rescales with norm (not sd)

    out.ls <- lm(y~xs)                      # ols fit on standardized
    beta.ols <- out.ls$coeff[2:(m+1)]       # ols except for intercept
    w <- abs(beta.ols)                      # weights for adaptive lasso
    xs <- scale(xs,center=FALSE,scale=1/w)  # xs times the weights
    object <- lars(xs,y,type="lasso",normalize=FALSE)

    # get min BIC
    # bic=log(n)*object$df+n*log(as.vector(object$RSS)/n)   # rss/n version
    sig2f <- summary(out.ls)$sigma^2        # full model mse
    bic2 <- log(n)*object$df+as.vector(object$RSS)/sig2f       # Cp version
    step.bic2 <- which.min(bic2)            # step with min BIC

    fit <- predict.lars(object,xs,s=step.bic2,type="fit",mode="step")$fit
    coeff <- predict.lars(object,xs,s=step.bic2,type="coef",mode="step")$coefficients
    coeff <- coeff*w/normx                  # get back in right scale
    st <- sum(coeff !=0)                    # number nonzero
    mse <- sum((y-fit)^2)/(n-st-1)          # 1 for the intercept

    # this next line just finds the variable id of coeff. not equal 0
    if(st>0) x.ind<-as.vector(which(coeff !=0)) else x.ind<-0
    intercept <- as.numeric(mean(y)-meanx%*%coeff)
    return(list(fit=fit,st=st,mse=mse,x.ind=x.ind,coeff=coeff,intercept=intercept,object=object,
                bic2=bic2,step.bic2=step.bic2))
}
