\name{sim.weed}
\alias{sim.weed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate weed exact count data and a covariate}
\description{
  The data simulated are of two kinds: exact weed counts and image
  derived indexes. The counts \eqn{w} arise from a transformed Gaussian model
  with Gamma marginal. The image indexes \eqn{i} relate to \eqn{w} as
  \eqn{i(s) = \lambda w(s) * \epsilon(s)}, where \eqn{\epsilon} is a spatially
  uncorrelated noise. The simulation is carried out on the unit square. 
}
\usage{
sim.weed(nx, ny, nxy, nz, param.cov, mu, sigma, lambda, tau, nbin, true.field = FALSE, npix = NULL, z.on.grid = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nx}{ number of sites where only \eqn{w} is observed }
  \item{ny}{ number of sites where only \eqn{i} is observed }
  \item{nxy}{  number of sites where \eqn{w} and \eqn{i} are observed }
  \item{nz}{  number of sites where \eqn{w} will be predicted. If
    \code{z.on.grid=TRUE}, \code{nz} should be a square.}
  \item{param.cov}{ Vector of parameters of the underlying Gaussian
    randomFields. The simulation of the Gaussian random fields performed
    assuming an exponential covariance function. The vector of
    parameters should have four components, namely: mea, variance,
    nugget, scale and should be given as
    e.g.:\code{c(mean=0,variance=1,nugget=0,scale=.1)}. See the
    documentation of \code{\link[RandomFields]{GaussRF}} in package
    \code{\link[RandomFields]{RandomFields}} for details.}
  \item{mu}{ Mean of \eqn{w} }
  \item{sigma}{ Variance  of \eqn{w} }
  \item{lambda}{ Scaling factor relating \eqn{w} and \eqn{i}}
  \item{tau}{ Variance of \eqn{\epsilon}}
  \item{nbin}{Number of bins for the binned data}
  \item{true.field}{Logical: shoud values on a dense grid be given}
  \item{npix}{ A vector giving the number of pixels horizontally and
    vertically of the grid if \code{true.field=TRUE}}
  \item{z.on.grid}{Logical: set TRUE if sites of prediction \eqn{z} are
    required on a grid; then\code{nz} should be a square.}
}
\details{}
\value{
  \item{x}{Coordinates of sites where weed counts only are given }
  \item{xy}{Coordinates of sites where weeds and images are given}
  \item{y}{Coordinates of sites where image indexes only are given }
  \item{z}{Coordinates of sites where predcited values are sought}
  \item{coord.grid}{Logical telling whether all the realisation of the
    weed random field is also given on a grid}
  \item{wx}{Weed values at sites \code{x}}
  \item{wy}{Weed values at sites \code{y}}
  \item{wxy}{Weed values at sites \code{xy}}
  \item{wz}{Weed values at sites \code{z}}
  \item{vx}{Binned weed values at sites \code{x}}
  \item{vxy}{Binned weed values at sites \code{xy}}
  \item{i}{Image values at sites \code{xy} and \code{y}. Values are
    concatenated in the same vector in this order}
  \item{wgrid}{Weed counts values at the node of a grid}
  \item{igrid}{Image index values at the node of a grid}
  \item{bin}{Binning of the weed counts}
  \item{param.cov}{Parameters of the covariance function of the
    underlying Gaussian random field}
  \item{mu}{Mean of \code{w}}
  \item{sigma}{Variance of \code{w}}
  \item{alpha}{Shape parameter of \code{w}}
  \item{beta}{Scale  parameter of \code{w}}
  \item{lambda}{Scaling  parameter relating  \code{w} and \code{i} as
    \eqn{i = \lambda w \epsilon}}
  \item{tau}{Variance of noise \eqn{epsilon}}
  \item{npix}{number of pixel in the horizontal direction (same value is
    assumed for the vertical direction)}
}
\references{ G. Guillot, N. Loren, M. Rudemo, Bayesian spatial prediction of weed intensities from  exact count 
  data  and picture based indexes, 2006, submitted}
\author{Gilles Guillot}



\keyword{}

