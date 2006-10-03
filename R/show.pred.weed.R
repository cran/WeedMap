"show.pred.weed" <-
function(sim=NULL,res,param=TRUE,pairs=TRUE,wy=TRUE,wz=TRUE,nit,thin,burnin)
  {
    if(burnin>0){subburn <- -(1:burnin)}else{subburn <- 1:(nit/thin)}
    if(!is.null(res$wy.MC)){wy.pred <- apply(res$wy.MC[,subburn],1,mean)}
    if(nrow(res$wz.MC)==1)
      {
       wz.pred <- mean(res$wz.MC[,subburn])}else
    {
      wz.pred <- apply(res$wz.MC[,subburn],1,mean)
    }
    if(param)
      {
        get(getOption("device"))()
        par(mfrow=c(4,2))

        ##
        plot(res$alpha.MC/res$beta.MC,type="l",
             sub=paste(mean(res$alpha.MC[subburn]/res$beta.MC[subburn])),
             main=paste(mean(res$alpha.MC[subburn])/mean(res$beta.MC[subburn])),ylab="mu")
         abline(h=mean(res$alpha.MC[subburn]/res$beta.MC[subburn]),lty=2)
        if(!is.null(sim))
          {abline(h=mean(res$alpha.MC/res$beta.MC),lty=2);abline(h=sim$mu,lty=2,col=2)}
        if(!is.null(sim$wx))
          {abline(h=mean(sim$wx),col=3,lty=3)}

        ##
        plot(res$alpha.MC/res$beta.MC^2,type="l",
             sub=paste(mean(res$alpha.MC[subburn]/res$beta.MC[subburn]^2)),
             main=(paste(mean(res$alpha.MC[subburn])/mean(res$beta.MC[subburn]^2))),ylab="sigma^2")
        abline(h=mean(res$alpha.MC[subburn]/res$beta.MC[subburn]^2),lty=2)
          if(!is.null(sim)){abline(h=sim$sigma^2,lty=2,col=2)}
        if(!is.null(sim$wx))
          {abline(h=var(sim$wx),col=3,lty=3)}

        ##
        plot(res$alpha.MC,type="l",
             main=paste(mean(res$alpha.MC[subburn])),ylab="alpha")
        abline(h=mean(res$alpha.MC[subburn]),lty=2);
        if(!is.null(sim))
          {abline(h=sim$alpha,lty=2,col=2)}

        ##
        plot(res$beta.MC,type="l",
             main=paste(mean(res$beta.MC[subburn])),ylab="beta")
        abline(h=mean(res$beta.MC[subburn]),lty=2);
        if(!is.null(sim))
          {abline(h=sim$beta,lty=2,col=2)}

        if(!is.null(res$i))
          {
            plot(res$tau.MC,type="l",
                 main=paste(mean(res$tau.MC[subburn])))
            abline(h=mean(res$tau.MC[subburn]),lty=2);
            if(!is.null(sim))
              {abline(h=sim$tau,lty=2,col=2)
               #abline(h=sqrt(var(sim$i/(sim$lambda*sim$wy))),lty=3,col=3)
             }
            
            plot(res$lambda.MC,type="l",main=paste(mean(res$lambda.MC)))
            abline(h=mean(res$lambda.MC[subburn]),lty=2);
            if(!is.null(sim))
              {
                abline(h=sim$lambda,lty=2,col=2)
#                abline(h=mean(sim$i/sim$wy),lty=3,col=3)
              }
          }

        plot(res$kappa.MC,type="l",ylim=c(0,res$kappa.max),
             main=paste(mean(res$kappa.MC[subburn])))
        abline(h=mean(res$kappa.MC),lty=2);
        if(!is.null(sim)){abline(h=sim$param.cov[4],lty=2,col=2)}

        plot(1:10,1:10,xlab="",ylab="",type="n",axes=FALSE)
        text(5,5,"MCMC simulation of parameters",cex=1.8)
      }
    
    if(pairs)
      {
        get(getOption("device"))()
        par(mfrow=c(1,2))
  
        if(!is.null(res$y))
          {
            mae=floor(mean(abs(sim$wy-wy.pred))*1000)/1000
            plot(sim$wy,wy.pred,
                 xlab="Simulated w data at sites y",
                 ylab="Prediction of w data at sites y",
                 sub=paste("MAE=",mae))
            abline(0,1)
            
      ##       plot(wy.pred,sim$i,
##                  ylab="Simulated image index at sites y",
##                  xlab="Prediction of w data at sites y")
##             abline(0,mean(res$lambda.MC))
##             abline(0,sim$lambda,col=2,lty=2)
          }
        
        mae=floor(mean(abs(sim$wz-wz.pred))*1000)/1000
            plot(sim$wz,wz.pred,
                 xlab="Simulated data at sites z",
                 ylab="Prediction at sites z",
                 sub=paste("MAE=",mae))
        abline(0,1)
          
      }

    
    if(wy)
      {
        if(!is.null(res$wy.MC))
          {
            get(getOption("device"))()
            par(mfrow=c(4,4))
            for(ky in 1:16)
              {
   ##              plot(1:(nit/thin),1:(nit/thin),
##                      ylim=c(min(c(sim$wy,res$wy.MC)),max(c(sim$wy,res$wy.MC))),
##                          type="n")
##                 points(1:(nit/thin),res$wy.MC[ky,],type="l")
                plot(1:(nit/thin),res$wy.MC[ky,],type="l",ylab="")
                abline(h=mean(res$wy.MC[ky,subburn]),col=1,lty=2,lwd=1.5)
                if(!is.null(sim)){abline(h=sim$wy[ky],col=2,lty=2,lwd=1.5)}
              }
          }
      }

    
    if(wz)
      {
        get(getOption("device"))()
        par(mfrow=c(4,4))
        for(kz in 1:16)
          {
            ##     plot(1:(nit/thin),1:(nit/thin),
            ##       ylim=c(min(c(sim$wz,res$wz.MC,res$wz.MC)),max(c(sim$wz,res$wz.MC,res$wz.MC))),
            ##          type="n")
            plot(1:(nit/thin),res$wz.MC[kz,],type="l",ylab="")
            abline(h=mean(res$wz.MC[kz,]),col=1,lty=2,lwd=1.5)
            if(!is.null(sim)){abline(h=sim$wz[kz],col=2,lty=2)}
          }
      }
    
    
    
  }

