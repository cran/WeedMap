"sim.weed" <-
function(nx,ny,nxy,nz,param.cov,mu,sigma,lambda,tau,nbin,
                          true.field=FALSE,npix=NULL,z.on.grid=TRUE)
  {
    alpha <- mu^2/sigma^2
    beta <- mu/sigma^2
    
    ## coord sampled in [0,1]x[0,1]
    if(nx>0){x <- matrix(nr=nx,nc=2,runif(2*nx))}else{x <- NULL}
    if(ny>0){y <- matrix(nr=ny,nc=2,runif(2*ny))}else{y <- NULL}
    if(nxy>0){xy <- matrix(nr=nxy,nc=2,runif(2*nxy))}else{xy <- NULL}
    if(z.on.grid)
      {
        z <- expand.grid(seq(0,1,length=sqrt(nz)),seq(0,1,length=sqrt(nz)))
        z <- as.matrix(z)
        colnames(z) <- NULL
                                        #grid(0,1,0,1,sqrt(nz),sqrt(nz))+.01*runif(1)
      }else
    {
      z <- cbind(runif(nz),runif(nz))
    }
    
    
    
    coord <- rbind(x,xy,y,z)
    if(true.field == TRUE)
      {
        coord.grid <- matrix(nc=2,nr=npix[1]*npix[2],NA)
        coord.grid[,1] <- rep(seq(0,1,length=npix[1]),npix[2])
        coord.grid[,2] <- as.vector(matrix(nr=npix[1],nc=npix[2],byr=TRUE,
                                           rep(seq(0,1,length=npix[2]),npix[1])))
      }
    if(true.field==TRUE){coord.all <- rbind(coord,coord.grid)}else
    {coord.all <- coord}

    ## sampling Gaussian random field
    g <- GaussRF(x=coord.all,
                 grid=FALSE,
                 model="exponential",
                 param=param.cov,
                 n=1)
    
    ## transforming into weed field
    w = qgamma(pnorm(g),shape=alpha,rate=beta)

    ## binned weed data
    bin <- seq(from=floor(min(w)),to=ceiling(max(w)),length=nbin)
    v <- numeric(nx+nxy)
    if((nx+nxy) > 0)
      {
        for(k in 1:(nx+nxy))
          {
            for(kk in 1:length(bin))
              {
                if(w[k] > bin[kk]) {v[k] <- kk}
              }
          }
      }
    
    ## add noise to get image values
    i = lambda*w * rgamma(n=length(w),shape=1/tau,rate=1/tau)
    
    wgrid=w[-(1:(nx+nxy+ny+nz))]
    igrid=i[-(1:(nx+nxy+ny+nz))]
    if(nx>0)
      {
        wx <- w[1:nx]
        vx <- v[1:nx]
      }else
    {
      wx <- NULL
      vx <- NULL
    }
    if(nxy>0)
      {
        wxy <- w[(nx+1):(nx+nxy)]
        vxy <- v[(nx+1):(nx+nxy)]
      }else
    {
      wxy <-NULL
      vxy <- NULL
    }
    if(ny>0){wy <- w[(nx+nxy+1):(nx+nxy+ny)]}else{wy <-NULL}
    if((nxy+ny)>0){i <- i[(nx+1):(nx+nxy+ny)]}else{i <- NULL}
    wz <- w[(nx+nxy+ny+1):(nx+nxy+ny+nz)]

    if(true.field==FALSE){coord.grid <-NULL} 

    list(x=x,
         xy=xy,
         y=y,
         z=z,
         coord.grid=coord.grid,
         wx=wx,
         wxy=wxy,
         wy=wy,
         wz=wz,
         vx=vx,
         vxy=vxy,
         i=i,
         wgrid=wgrid,
         igrid=igrid,
         bin=bin,
         param.cov=param.cov,
         mu=mu,
         sigma=sigma,
         alpha=alpha,
         beta=beta,
         lambda=lambda,
         tau=tau,
         npix=npix)
  }

