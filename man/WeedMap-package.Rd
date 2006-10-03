\name{WeedMap-package}
\alias{WeedMap-package}
\alias{WeedMap}
\docType{package}
\title{Spatial prediction of weed intensities}
\description{Simulation, inference and prediction for a Bayesian spatial
  statistical model for weed intensities and a covariate}
\details{
\tabular{ll}{
Package: \tab WeedMap\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2006-09-28\cr
License: \tab GPL\cr
}
The function \code{\link{sim.weed}} makes simulation from the model. The
simulated dataset can be graphically displayed with \code{\link{show.sim.weed}}. The
function \code{\link{pred.weed}} makes inference and prediction whose results
can be graphically displayed with \code{\link{show.pred.weed}} (monitoring of Markov
chain simulation) and with \code{\link{check.model.weed}} (goodness of fit assessment).
}
\author{Gilles Guillot
   \url{www.inapg.inra.fr/ens_rech/mathinfo/personnel/guillot/welcome.html}
}
\references{ G. Guillot, N. Loren, M. Rudemo, Bayesian spatial prediction of weed intensities from  exact count 
  data  and picture based indexes, 2006, submitted}

\keyword{ package }

\examples{


## Simulate a data set
sim <- sim.weed(nx=30, ny=20, nxy=20, nz=49,
                param.cov=c(mean=0,variance=1,nugget=0,scale=.1),
                mu=80, sigma=70, lambda=1, tau=0.2, nbin=10,
                true.field = TRUE, npix =  c(100,100), z.on.grid = TRUE)



## show the graphics
show.sim.weed(sim)

\dontrun{

## make joint inference and prediction
res <- pred.weed(nit=10000,
                 thin=10,
                 ## data
                 x=sim$x,
                 xy=sim$xy,
                 y=sim$y,
                 z=sim$z,
                 wx=sim$wx,
                 wxy=sim$wxy,
                 i=sim$i,
                 ## init
                 #alpha=alpha,
                 #beta=beta,
                 lambda=1,
                 #tau=tau,
                 #kappa=kappa,
                 ## proposals
                 sd.prop.h=0.1,
                 sd.prop.alpha=0.1,
                 sd.prop.beta=0.01,
                 sd.prop.lambda=0.,
                 sd.prop.tau=0.5,
                 delta.prop.kappa=2,
                 ## priors
                 mprior.alpha=0.625,
                 vprior.alpha=1,
                 mprior.beta=0.0125,
                 vprior.beta=1,
                 mprior.kappa=.5,
                 vprior.kappa=999,
                 mprior.lambda=1,
                 vprior.lambda=1,
                 mprior.tau=0.1,
                 vprior.tau=10,
                 n.kappa=30,
                 kappa.max=5*sim$param.cov[4])

show.pred.weed(sim=sim,
               res=res,
               param=TRUE,
               pairs=TRUE,
               wy=FALSE,
               wz=FALSE,
               nit=res$nit,
               thin=res$thin,
               burnin=500)

check.model.weed(x=sim$x,
                 xy=sim$xy,
                 y=sim$y,
                 wx=sim$wx,
                 wxy=sim$wxy,
                 i=sim$i,
                 ## output of MCMC run
                 res=res,                           
                 ## options
                 nit=res$nit,
                 thin=res$thin,
                 burnin=500,
                 bin=seq(.1,.5,.05),
                 nqqplot=500,
                 nresamp=200)

}
}
