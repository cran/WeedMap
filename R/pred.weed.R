"pred.weed" <-
function(nit,
                      thin=1,
                      ## data
                      x=NULL,xy=NULL,y=NULL,z,
                      wx=NULL,wxy=NULL,i=NULL,
                      ## init
                      alpha=NULL,beta=NULL,lambda=NULL,tau=NULL,kappa=NULL,
                      ## proposals 
                      sd.prop.h=0.1,
                      sd.prop.alpha=0.5,
                      sd.prop.beta=0.5,
                      sd.prop.lambda=0.1,
                      sd.prop.tau=0.1,
                      delta.prop.kappa=2,
                      # priors
                      mprior.alpha,
                      vprior.alpha,
                      mprior.beta,
                      vprior.beta,
                      mprior.kappa,
                      vprior.kappa=999,
                      mprior.tau=1,
                      vprior.tau=999,
                      mprior.lambda=1,
                      vprior.lambda=999,
                      n.kappa=1,
                      kappa.max=0)
  {
    ## x:  sites where only w is observed
    ## y:  sites where only i is observed
    ## xy: sites where w and i are observed 
    ## z:  sites where w will be predicted

    ## variables
    ## hx,hxy,hy,hz
    ## gx,gxy,gy,gz
    ## wx,wxy,wy,wz
    ## i
    ## alpha,beta,lambda,tau,kappa
    ## wy.MC,wz.MC
    ## alpha.MC,beta.MC,lambda.MC,tau.MC,kappa.MC
    
    
    ####################
    # Dimension output arrays
    nx <- ifelse(is.null(x),0,nrow(x))
    nxy <- ifelse(is.null(xy),0,nrow(xy))
    ny <- ifelse(is.null(y),0,nrow(y))
    nz <- nrow(z)
    if(ny>0){wy.MC <- matrix(nr=ny,nc=nit/thin)}else{wy.MC <- NULL}
    wz.MC <- matrix(nr=nz,nc=nit/thin)
    alpha.MC <- rep(NA,nit/thin)
    beta.MC <- rep(NA,nit/thin)
    lambda.MC <- rep(NA,nit/thin)
    tau.MC <- rep(NA,nit/thin)
    kappa.MC <- rep(NA,nit/thin)
    seq.kappa <- seq(kappa.max/n.kappa,kappa.max,length=n.kappa)
    if(is.null(kappa))
      {
        ikappa <- sample(x=1:n.kappa,size=1)
      }else{
        ikappa <- max(1,floor(n.kappa*kappa/kappa.max))
      }
    kappa <- seq.kappa[ikappa]
    

    ##################################
    ## Compute Choleski decompositions
    print(" Computing Choleski decompositions ")
    coord <- rbind(x,xy,y,z)
    lag <- as.matrix(dist(coord,diag=T))
    U <- array(dim=c(nrow(coord),nrow(coord),n.kappa))
    invU <- array(dim=c(nrow(coord),nrow(coord),n.kappa))
    jac <- rep(NA,n.kappa)
    for(ikap in 1:n.kappa)
      {
        kap <- seq.kappa[ikap]
        U[,,ikap] <- chol(exp(-lag/kap))
        invU[,,ikap] <- solve( U[,,ikap])
        jac[ikap] <- prod(diag(invU[,,ikap]))
      }
    
    ###############################
    ## Init variables

    ## alpha, beta
    if((nx+nxy)>0)
      {
        if(is.null(alpha) & is.null(beta))
          {
            alpha <- mean(c(wx,wxy))^2/var(c(wx,wxy))
            beta <- mean(c(wx,wxy))/var(c(wx,wxy))
          }
        if(is.null(alpha) & !is.null(beta))
          {
            alpha <- mean(c(wx,wxy))*beta
          } 
        if(!is.null(alpha) & is.null(beta))
          {
            beta <- alpha/mean(c(wx,wxy))
          }
      }else
    {
      if(is.null(alpha) & is.null(beta))
        {
          alpha <- mean(i)^2/var(i)
          beta <- mean(i)/var(i)
        }
      if(is.null(alpha) & !is.null(beta))
        {
          alpha <- mean(i)*beta
        } 
      if(!is.null(alpha) & is.null(beta))
        {
          beta <- alpha/mean(i)
        }
    }
    ## Init h,g values at x 
    if(is.null(wx))
      {
        hx <- gx <- NULL
      }else
      {
        gx <- qnorm(pgamma(wx,shape=alpha,rate=beta,log.p=TRUE),log.p=TRUE)
        hx <- gx %*% invU[1:nx,1:nx,ikappa]
      }
    ## Init h,g values at xy
    if(is.null(wxy))
      {
        hxy <- gxy <- NULL
      }else
      {
        gxy <- qnorm(pgamma(wxy,shape=alpha,rate=beta,log.p=TRUE),log.p=TRUE)
        hxy <- (c(gx,gxy) %*% invU[1:(nx+nxy),1:(nx+nxy),ikappa])[(nx+1):(nx+nxy)]
        #hxy <- gxy %*% invU[(nx+1):(nx+nxy),(nx+1):(nx+nxy),ikappa]
      }

    ## Init h,g,w values at y 
    if(is.null(i) | ny==0)
      {
        hy <- gy <- wy <- NULL
      }else
    {
      if(ny>0)
        {
          hy <- rnorm(n=ny)
          gy <- (c(hx,hxy,hy) %*% U[1:(nx+nxy+ny),1:(nx+nxy+ny),ikappa])[(nx+nxy+1):(nx+nxy+ny)]
          wy <- qgamma(pnorm(gy,log.p=TRUE),shape=alpha,rate=beta,log.p=TRUE)
        }
    }

    ## lambda    
    if(is.null(lambda))
      {
        lambda <- mean(i/c(wxy,wy))
      }
    ## eps    
    if(!is.null(i))
      {
        eps <- i / (lambda*c(wxy,wy))
      }else {eps <- NULL}
    ## tau
    if(is.null(tau) & !is.null(i))
      {
        tau=sqrt(var(eps))
      }
    ## Init h,g,w values at z
    hz <- rnorm(n=nz)
    gz <- (c(hx,hxy,hy,hz) %*% U[,,ikappa])[(nx+nxy+ny+1):(nx+nxy+ny+nz)]
    wz <- qgamma(pnorm(gz,log.p=TRUE),shape=alpha,rate=beta,log.p=TRUE)
    ## Store intial state
    if(ny>0){wy.MC[,1] <- wy}
    wz.MC[,1] <- wz
    alpha.MC[1] <- alpha
    beta.MC[1] <- beta
    lambda.MC[1] <- lambda
    if(!is.null(i)){tau.MC[1] <- tau}
    kappa.MC[1] <- kappa

    
    ####################
    ## Start Markov chain
    for(kit in 2:nit)
      {
        if(kit%%thin == 0){print(kit/nit*100)}

        ############
        ## Update wy
        if(ny > 0 & sd.prop.h != 0) 
          {
            ## pick and update a component of wy
            ky <- sample(1:ny,1)

            ## propose new hy value
            hy.star <- hy
            hy.star[ky] <- hy.star[ky] + rnorm(n=1,sd=sd.prop.h)
            ## compute proposed gy value
            gy <-  (c(hx,hxy,hy) %*% U[1:(nx+nxy+ny),1:(nx+nxy+ny),ikappa])[(nx+nxy+1):(nx+nxy+ny)]
            gy.star <- (c(hx,hxy,hy.star) %*% U[1:(nx+nxy+ny),1:(nx+nxy+ny),ikappa])[(nx+nxy+1):(nx+nxy+ny)]
            wy <- qgamma(pnorm(gy,log.p=TRUE),
                               shape=alpha,rate=beta,log.p=TRUE)
            wy.star <-  qgamma(pnorm(gy.star,log.p=TRUE),
                               shape=alpha,rate=beta,log.p=TRUE)
            eps.star <- i / (lambda*c(wxy,wy.star))
            ## compute MH ratio
            lR <- -0.5*(hy.star[ky]^2    - hy[ky]^2) +
              ## jacobian of h->w
              sum(dnorm(gy,log=TRUE)) - sum(dnorm(gy.star,log=TRUE))  +
                sum(dgamma(wy.star,shape=alpha,rate=beta,log=TRUE)) -
                  sum(dgamma(wy,shape=alpha,rate=beta,log=TRUE)) +
                    ## density of eps
                    sum(dgamma(x=eps.star,shape=1/tau,rate=1/tau,log=TRUE)) - 
                      sum(dgamma(x=eps,shape=1/tau,rate=1/tau,log=TRUE)) +
                       ## jacobian of w->i
                        sum(log(wy)-log(wy.star))                  
            R <- exp(lR)
            b <- rbinom(n=1,size=1,prob=min(1,R))
            if(b ==1)
              {
                hy <- hy.star
                gy <- gy.star
                wy <- wy.star
                eps <- eps.star
              }
          }
        if(kit%%thin == 0)
          {
            wy.MC[,kit/thin] <- wy
 ##            gz <- (c(hx,hxy,hy,hz) %*% U[,,ikappa])[(nx+nxy+ny+1):(nx+nxy+ny+nz)]
##             wz <- qgamma(pnorm(gz,log.p=TRUE),
##                          shape=alpha,rate=beta,log.p=TRUE)
##             wz.MC[,kit/thin] <- wz
          }


        #######################################
        ## Update wz by a Gibbs step at sites z
        hz <- rnorm(nz)
        gz <- (c(hx,hxy,hy,hz) %*% U[,,ikappa])[(nx+nxy+ny+1):(nx+nxy+ny+nz)]
        wz <- qgamma(pnorm(gz,log.p=TRUE),
                     shape=alpha,rate=beta,log.p=TRUE)
        if(kit%%thin == 0)
          {
            wz.MC[,kit/thin] <- wz
          }

        ###############
        ## Update alpha
        #print("Update alpha")
        alpha.star <- alpha + rnorm(1)*sd.prop.alpha
        if((alpha.star >0) & (sd.prop.alpha !=0))
          {
            w <- c(wx,wxy,wy,wz)
            g <- c(gx,gxy,gy,gz)
            h <- c(hx,hxy,hy,hz)
            g.star <- qnorm(pgamma(w,shape=alpha.star,rate=beta,log.p=TRUE),
                            log.p=TRUE)
            h.star <- g.star %*% invU[,,ikappa]
            lR <- -0.5*(sum(h.star^2)-sum(h^2)) +
              sum(dnorm(g,log=TRUE)) - sum(dnorm(g.star,log=TRUE)) +
                sum(dgamma(w,shape=alpha.star,rate=beta,log=TRUE)) -
                  sum(dgamma(w,shape=alpha,rate=beta,log=TRUE)) +
                    dgamma(x=alpha.star,
                           shape=mprior.alpha^2/vprior.alpha,
                           rate=mprior.alpha/vprior.alpha,log=TRUE) -
                             dgamma(x=alpha,
                                    shape=mprior.alpha^2/vprior.alpha,
                                    rate=mprior.alpha/vprior.alpha,log=TRUE)
            R <- exp(lR)
            b <- rbinom(n=1,size=1,prob=min(1,R))
            if(b == 1)
              {
                alpha <- alpha.star
                if(!is.null(wx))
                   {
                     hx <- h.star[1:nx]
                     gx <- g.star[1:nx]
                   }
                if(nxy>0)
                  {
                    hxy <- h.star[(nx+1):(nx+nxy)]
                    gxy <- g.star[(nx+1):(nx+nxy)]
                  }
                if(ny>0)
                  {
                    hy <- h.star[(nx+nxy+1):(nx+nxy+ny)]
                    gy <- g.star[(nx+nxy+1):(nx+nxy+ny)]
                  }
                hz <- h.star[(nx+nxy+ny+1):(nx+nxy+ny+nz)]
                gz <- g.star[(nx+nxy+ny+1):(nx+nxy+ny+nz)]
              }
          }
        if(kit%%thin == 0)
          {
            alpha.MC[kit/thin] <- alpha
          }

        ################################
        ## Update jointly beta and kappa
        #print("Update jointly beta and kappa")
        beta.star <- beta + rnorm(n=1)*sd.prop.beta
        ikappa.star <- ikappa +
          floor((2*delta.prop.kappa+1)*runif(1))-delta.prop.kappa
        if(((beta.star > 0) & ((ikappa.star >= 1)  & (ikappa.star <= n.kappa)))
           & (sd.prop.beta !=0 | delta.prop.kappa !=0))
          {
            kappa.star <- seq.kappa[ikappa.star]
            h <- c(hx,hxy,hy,hz)
            g <- c(gx,gxy,gy,gz)
            w <- c(wx,wxy,wy,wz)
            g.star <- qnorm(pgamma(w,shape=alpha,rate=beta.star,log.p=TRUE),
                            log.p=TRUE)
            h.star <- g.star %*% invU[,,ikappa.star]
            lR <- sum(-0.5*(h.star^2-h^2)) +
              log(jac[ikappa.star]/jac[ikappa]) +
                sum(dnorm(g,log=TRUE))-
                  sum(dnorm(g.star,log=TRUE))+
                    sum(dgamma(w,shape=alpha,rate=beta.star,log=TRUE))-
                      sum(dgamma(w,shape=alpha,rate=beta,log=TRUE)) +
                        dgamma(x=beta.star,
                               shape=mprior.beta^2/vprior.beta,
                               rate=mprior.beta/vprior.beta,log=TRUE) - 
                                 dgamma(x=beta,
                                        shape=mprior.beta^2/vprior.beta,
                                        rate=mprior.beta/vprior.beta,log=TRUE) +
                                          dgamma(x=kappa.star,log=TRUE,
                                                 shape=mprior.kappa^2/vprior.kappa,
                                                 rate=mprior.kappa/vprior.kappa) -
                                            dgamma(x=kappa,log=TRUE,
                                                 shape=mprior.kappa^2/vprior.kappa,
                                                 rate=mprior.kappa/vprior.kappa)
            R <- exp(lR)
            b <- rbinom(n=1,size=1,prob=min(1,R))
            if(b == 1)
              {
                kappa <- kappa.star
                ikappa <- ikappa.star
                beta <- beta.star
                if(!is.null(wx))
                   {
                     hx <- h.star[1:nx]
                     gx <- g.star[1:nx]
                   }
                if(nxy>0)
                  {
                    hxy <- h.star[(nx+1):(nx+nxy)]
                    gxy <- g.star[(nx+1):(nx+nxy)]
                  }
                if(ny>0)
                  {
                    hy <- h.star[(nx+nxy+1):(nx+nxy+ny)]
                    gy <- g.star[(nx+nxy+1):(nx+nxy+ny)]
                  }
                hz <- h.star[(nx+nxy+ny+1):(nx+nxy+ny+nz)]
                gz <- g.star[(nx+nxy+ny+1):(nx+nxy+ny+nz)]
              } 
          }
        if(kit%%thin == 0)
          {
            beta.MC[kit/thin] <- beta
            kappa.MC[kit/thin] <- kappa
          }

        ## joint update (lambda,tau)
     # print("joint update (lambda,tau)")
        if(!is.null(i) & (sd.prop.lambda !=0 | sd.prop.tau !=0))
          {
            lambda.star <- lambda + rnorm(n=1)*sd.prop.lambda
            eps.star <- i / (lambda.star*c(wxy,wy))
            tau.star <- tau + rnorm(n=1)*sd.prop.tau
            if((lambda.star > 0) & (tau.star > 0))
              {
                lR <- sum(dgamma(x=eps.star,shape=1/tau.star,rate=1/tau.star,log=TRUE)) -
                  sum(dgamma(x=eps,shape=1/tau,rate=1/tau,log=TRUE)) +
                    (nxy+ny)*(log(lambda)-log(lambda.star)) +
                    dgamma(x=lambda.star,
                           shape=mprior.lambda^2/vprior.lambda,
                           rate=mprior.lambda/vprior.lambda,log=TRUE) - 
                             dgamma(x=lambda,
                                    shape=mprior.lambda^2/vprior.lambda,
                                    rate=mprior.lambda/vprior.lambda,log=TRUE) +
                                      dgamma(x=tau.star,
                                             shape=mprior.tau^2/vprior.tau,
                                             rate=mprior.tau/vprior.tau,log=TRUE) - 
                                               dgamma(x=tau,
                                                      shape=mprior.tau^2/vprior.tau,
                                                      rate=mprior.tau/vprior.tau,log=TRUE)
                R <- exp(lR)
                b <- rbinom(n=1,size=1,prob=min(1,R))
                if(b == 1)
                  {
                    lambda <- lambda.star
                    eps <- eps.star
                    tau <- tau.star
                  }
              }
            if(kit%%thin == 0)
              {
                lambda.MC[kit/thin] <- lambda
                tau.MC[kit/thin] <- tau
                                        #if(!is.null(i)){tau.MC[kit/thin] <- tau}
              }
          }
            
      }
    

    list(x=x,
         y=y,
         z=z,
         wx=wx,
         i=i,
         nit=nit,
         thin=thin,
         wy.MC=wy.MC,
         wz.MC=wz.MC,
         alpha.MC=alpha.MC,
         beta.MC=beta.MC,
         lambda.MC=lambda.MC,
         tau.MC=tau.MC,
         kappa.MC=kappa.MC,
         n.kappa=n.kappa,
         kappa.max=kappa.max)
  }

