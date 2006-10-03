"check.model.weed" <-
function(
                             ## data
                             x,xy,y,wx,wxy,i,
                             ## output of MCMC run
                             res,                           
                             ## options
                             nit,thin,burnin,bin,nqqplot,nresamp=100)
  {
    if(burnin>0){subburn <- -(1:burnin)}else{subburn <- 1:(nit/thin)}
    if(!is.null(xy)){nxy <- nrow(xy)}
    if(!is.null(y)){ny <- nrow(y)}else{ny <- 0}

    
    ## estimated parameters
    alpha <- mean(res$alpha.MC[subburn])
    beta <- mean(res$beta.MC[subburn])
    kappa <-  mean(res$kappa.MC[subburn])
    lambda <-  mean(res$lambda.MC[subburn])
    tau <-  mean(res$tau.MC[subburn])

    ## residuals
    if(nxy>0){epsxy <- i[1:nxy]/(lambda*wxy)}
    if(ny>0)
      {
        wy.pred <- apply(res$wy.MC[,subburn],1,mean)
        epsy <- i[(nxy+1):(nxy+ny)]/(lambda*wy.pred)
      }
    eps <- epsxy

    ############
    ## residuals 
    vario.resid <- EmpiricalVariogram(x=xy,
                                    data=eps,
                                    bin=bin,
                                    grid=FALSE)

#    get(getOption("device"))()
    par(mfrow=c(2,3))
   
    vario.resamp <- matrix(nr=length(bin)+1,nc=nresamp)
    for(iresamp in 1:nresamp)
      {
        resamp <- rgamma(shape=1/tau,rate=1/(lambda*tau),n=nxy)
        vario.resamp[,iresamp] <- EmpiricalVariogram(x=xy,
                                                     data=resamp,
                                                     bin=bin,
                                                     grid=FALSE)$emp.vario
      }

    vario.envelope <- matrix(nr=length(bin)+1,nc=4)
    vario.envelope[,1] <- apply(vario.resamp,1,quantile,prob=0.025)
    vario.envelope[,2] <- apply(vario.resamp,1,quantile,prob=0.5)
    vario.envelope[,3] <- apply(vario.resamp,1,quantile,prob=0.95)
    vario.envelope[,4] <- apply(vario.resamp,1,quantile,prob=0.975)
    
    
    plot(vario.resid$centers,vario.resid$emp.vario,type="b",
         ylim=c(0,max(vario.envelope)),
         col=2,sub="Residuals",
         xlab="Spatial lag",ylab="Semi variance",xlim=c(0,max(bin)))
    lines(vario.resid$centers,vario.envelope[,1],lty=2)
    lines(vario.resid$centers,vario.envelope[,4],lty=2)
    lines(vario.resid$centers,vario.envelope[,2],lty=3)
    lines(vario.resid$centers,vario.envelope[,3],lty=3)
    
    hist(eps,prob=TRUE,breaks=50,main="",sub="Residuals",xlab="")
    lines(seq(0,max(eps),0.1),
          dgamma(x=seq(0,max(eps),0.1),shape=1/tau,rate=1/(lambda*tau)),col=2)
    qqplot(rgamma(shape=1/tau,rate=1/(lambda*tau),n=nqqplot),eps,
           ylab="Residuals",xlab="Quantiles of theor. Gamma dist.");
    abline(0,1,lty=2,col=2)

    ## h values at x and xy sites 
    lag <- as.matrix(dist(rbind(x,xy),diag=T))
    U <- chol(exp(-lag/kappa))
    invU <- solve( U)
    g.est <- qnorm(pgamma(c(wx,wxy),shape=alpha,rate=beta,log.p=TRUE),log.p=TRUE)
    h.est <- g.est %*% invU

    vario.g.est<- EmpiricalVariogram(x=rbind(x,xy),
                                 data=g.est,
                                 bin=bin,
                                 grid=FALSE)
 
    g.resamp <- GaussRF(x=rbind(x,xy),
                 grid=FALSE,
                 model="exponential",
                 param=c(0,variance=1,nugget=0,scale=kappa),
                 n=nresamp)

    for(iresamp in 1:nresamp)
      {
        vario.resamp[,iresamp] <- EmpiricalVariogram(x=rbind(x,xy),
                                           data=g.resamp[,iresamp],
                                           bin=bin,
                                           grid=FALSE)$emp.vario
      }
    vario.envelope <- matrix(nr=length(bin)+1,nc=4)
      vario.envelope[,1] <- apply(vario.resamp,1,quantile,prob=0.025)
    vario.envelope[,2] <- apply(vario.resamp,1,quantile,prob=0.5)
    vario.envelope[,3] <- apply(vario.resamp,1,quantile,prob=0.95)
    vario.envelope[,4] <- apply(vario.resamp,1,quantile,prob=0.975)

    plot(vario.g.est$centers,vario.g.est$emp.vario,type="b",ylim=c(0,max(vario.envelope)),
         col=2,xlab="Spatial lag",ylab="Semi variance",
         sub="Variogram of transformed Gaussian values",xlim=c(0,max(bin)))
    lines(vario.resid$centers,vario.envelope[,1],lty=2)
    lines(vario.resid$centers,vario.envelope[,4],lty=2)
    lines(vario.resid$centers,vario.envelope[,2],lty=3)
    lines(vario.resid$centers,vario.envelope[,3],lty=3)
  

    hist(h.est,breaks=50,prob=TRUE,main="",sub="Transformed Gaussian values",xlab="")
    lines(seq(min(h.est),max(h.est),0.1),
          dnorm(seq(min(h.est),max(h.est),0.1)),col=2)
    qqnorm(h.est,main="",sub="Transformed Gaussian values");abline(0,1,col=2)
    
    
  }

