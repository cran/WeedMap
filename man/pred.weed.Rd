\name{pred.weed}
\alias{pred.weed}

\title{Spatial prediction of  weed counts}
\description{Makes joint Bayesian inference and spatial prediction
  within a model for weed count data and a covariate}

\usage{
pred.weed(nit, thin = 1, x = NULL, xy = NULL, y = NULL, z, wx = NULL,
wxy = NULL,
i = NULL, alpha = NULL, beta = NULL, lambda = NULL, tau = NULL, kappa =
NULL,
sd.prop.h = 0.1, sd.prop.alpha = 0.5, sd.prop.beta = 0.5, sd.prop.lambda
= 0.1, sd.prop.tau = 0.1,
delta.prop.kappa = 2, mprior.alpha, vprior.alpha, mprior.beta, vprior.beta, mprior.kappa, vprior.kappa = 999, mprior.tau = 1, vprior.tau = 999, mprior.lambda = 1, vprior.lambda = 999, n.kappa = 1, kappa.max = 0)
}

\arguments{
  \item{nit}{Number of iterations}
  \item{thin}{Thinning of the Markov chain}
  \item{x}{Coordinates of sites where weed counts only are given }
  \item{xy}{Coordinates of sites where weeds and images are given}
  \item{y}{Coordinates of sites where image indexes only are given }
  \item{z}{Coordinates of sites where predcited values are sought}
  \item{wx}{Weed values at sites \code{x}}
  \item{wxy}{Weed values at sites \code{xy}} 
  \item{i}{Image values at sites \code{xy} and \code{y}. Values have to
    be     concatenated in the same vector in this order}
  \item{alpha}{Init value for the shape parameter of \code{w}}
  \item{beta}{Init value for the scale  parameter of \code{w}}
  \item{lambda}{Init value for the scaling  parameter relating  \code{w} and \code{i} as
    \eqn{i = \lambda w \epsilon}}
  \item{tau}{Init value for the variance of noise \eqn{epsilon}}
  \item{kappa}{Init value for the spatial scale parameter of the Gausian random field}
  \item{sd.prop.h}{Standard deviation of the Gaussian increment in the
    proposal of Gaussian components \eqn{h}}
  \item{sd.prop.alpha}{Standard deviation of the Gaussian increment in the
    proposal for \code{alpha}}
  \item{sd.prop.beta}{Standard deviation of the Gaussian increment in the
    proposal for \code{beta}}
  \item{sd.prop.lambda}{Standard deviation of the Gaussian increment in the
    proposal for \code{lambda}}
  \item{sd.prop.tau}{Standard deviation of the Gaussian increment in the
    proposal for \code{tau}}
  \item{delta.prop.kappa}{Maximum number of increments allowed (the
    amplitude of an increment being \code{kappa.max}/\code{n.kappa}) for a
    move in the proposal of \code{kappa}}
  \item{mprior.alpha}{A priori mean of \code{alpha} }
  \item{vprior.alpha}{A priori variance of \code{alpha} }
  \item{mprior.beta}{ A priori mean of \code{beta} }
  \item{vprior.beta}{ A priori variance of \code{beta} }
  \item{mprior.kappa}{A priori mean of \code{kappa} }
  \item{vprior.kappa}{A priori variance of \code{kappa} }
  \item{mprior.tau}{ A priori mean of \code{tau} }
  \item{vprior.tau}{A priori variance of \code{tau} }
  \item{mprior.lambda}{  A priori mean of \code{lambda} }
  \item{vprior.lambda}{A priori variance of \code{lambda} }
  \item{n.kappa}{Number of steps in the discretisation of the support of
    \code{kappa}}
  \item{kappa.max}{Maximum value in the truncation of the support of \code{kappa}}
}
\details{If  standard deviation of the Gaussian increment in the update
  of alpha, beta, lambda or tau, or if the step in the
  increment of kappa is equal to 0, then this variable is not
processed in the MCMC run and stays at its initial value. This is the
way to specify that inference should not be made on a one or several
variable. See examples in \code{\link{WeedMap}} where lambda is
initialised at 1 and not updated (\code{sd.prop.lambda=0}).

If init values are not given, the corresponding parameters are 
initialised from the prior.  
}

\value{A list whose elements are:
  \item{x}{Coordinates of sites where weed counts only are given }
  \item{y}{Coordinates of sites where image indexes only are given }
  \item{z}{Coordinates of sites where predcited values are sought}
  \item{wx}{Weed values at sites \code{x}}
  \item{i}{Image values at sites \code{xy} and \code{y}. Values are
    concatenated in the same vector in this order}  
  \item{nit}{Number of iterations}
  \item{thin}{Thinning of the Markov chain}
  \item{wy.MC}{A matrix with \code{ny} rows and \code{nit/thin}
    columns containing sampled values of \code{wy}}
  \item{wz.MC}{A matrix with \code{nz} rows and \code{nit/thin}
    columns containing sampled values of \code{wz}}
  \item{alpha.MC}{A vector of length \code{nit/thin} of simulated
    \code{alpha} values}
  \item{beta.MC}{A vector of length \code{nit/thin} of simulated
    \code{beta} values}
  \item{lambda.MC}{A vector of length \code{nit/thin} of simulated
    \code{lambda} values}
  \item{tau.MC}{A vector of length \code{nit/thin} of simulated
    \code{tau} values}
  \item{kappa.MC}{A vector of length \code{nit/thin} of simulated
    \code{kappa} values}
  \item{n.kappa}{Number of steps in the discretisation of the support of
    \code{kappa}}
  \item{kappa.max}{Maximum value in the truncation of the support of \code{kappa}}
}
\references{ G. Guillot, N. Loren, M. Rudemo, Bayesian spatial prediction of weed intensities from  exact count 
  data  and picture based indexes, 2006, submitted}
\author{Gilles Guillot}



\keyword{}

    