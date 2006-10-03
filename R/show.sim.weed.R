"show.sim.weed" <-
function(sim)
{
  get(getOption("device"))()
  par(mfrow=c(3,3))
  
  #setplot(seq(0,1,length=sim$npix[1]),seq(0,1,length=sim$npix[2]))

    if(!is.null(sim$x))
      {
        plot(sim$x,main="Weed measurement sites only",
             xlab="",ylab="",xlim=c(-0.02,1.02),ylim=c(-0.02,1.02))
      }
  
  if(!is.null(sim$xy))
    {
      plot(sim$xy,main="Homotopy sites",
           xlab="",ylab="",xlim=c(-0.02,1.02),ylim=c(-0.02,1.02))
    }

  if(!is.null(sim$wy))
    {
      plot(sim$y,main="Image index sites only",
           xlab="",ylab="",xlim=c(-0.02,1.02),ylim=c(-0.02,1.02))
    }

  if(!is.null(sim$npix))
    {
      image(seq(0,1,length=sim$npix[1]),
            seq(0,1,length=sim$npix[2]),
            matrix(nr=sim$npix[1],nc=sim$npix[2],
                   sim$wgrid,byrow=FALSE),
            xlab="",ylab="",col=terrain.colors(25),main="True weed field",
            xlim=c(-0.02,1.02),ylim=c(-0.02,1.02))
    }
  
  if(!is.null(sim$x) | !is.null(sim$xy))
    {
      image.plot(as.image(Z=c(sim$wx,sim$wxy),x=rbind(sim$x,sim$xy)),
                 main="Exact count data",xlim=c(-0.02,1.02),ylim=c(-0.02,1.02),
                 col=terrain.colors(25))
    }

  if(!is.null(sim$i))
    {
      image.plot(as.image(Z=sim$i,x=rbind(sim$xy,sim$y)),
                 main="Image index data",xlim=c(-0.02,1.02),ylim=c(-0.02,1.02),
                 col=terrain.colors(25))
    }
  
  plot(sim$z,xlab="",ylab="",
       main="Prediction sites",xlim=c(-0.02,1.02),ylim=c(-0.02,1.02),pch=16,cex=.5)
  
  if(!is.null(sim$wxy) | !is.null(sim$wy))
    {
      plot(c(sim$wxy,sim$wy),sim$i,
           xlab="Simulated exact counts",ylab="Simulated image index")
      abline(0,sim$lambda,col=2)
    }
  
}

