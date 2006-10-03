\name{check.model.weed}
\alias{check.model.weed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Check goodness of fit of the model}
\description{Plot various graphics to assess goddness of fit of the
  transformed Gaussian model}
\usage{
check.model.weed(x, xy, y, wx, wxy, i, res, nit, thin, burnin, bin, nqqplot, nresamp = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Coordinates of sites where weed counts only are given }
  \item{xy}{Coordinates of sites where weeds and images are given}
  \item{y}{Coordinates of sites where image indexes only are given }
  \item{wx}{Weed values at sites \code{x}}
  \item{wxy}{Weed values at sites \code{xy}}
   \item{i}{Image values at sites \code{xy} and \code{y}. Values are
    concatenated in the same vector in this order}
  \item{res}{Output of \code{pred.weed}}
   \item{nit}{Number of iterations}
   \item{thin}{Thinning of the Markov chain}
    \item{burnin}{Number of  values to discard before computing statistics
    about the empirical posterior}
   \item{bin}{Binning for the empirical variograms}
  \item{nqqplot}{Number of observation in the qqplot}
  \item{nresamp}{Number of realisation resampled to get enveloppe estimates}
}
\details{}
\value{None}

\references{ G. Guillot, N. Loren, M. Rudemo, Bayesian spatial prediction of weed intensities from  exact count 
  data  and picture based indexes, 2006, submitted}
\author{Gilles Guillot}

\examples{}

\keyword{}
