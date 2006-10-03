\name{show.pred.weed}
\alias{show.pred.weed}

\title{Plot results of pred.weed}
\description{Plot traces of Markov chains simulated by
  \code{pred.weed}, compute statistics on the posterior distribution,
  compare to known true values. (if any)}

\usage{
show.pred.weed(sim = NULL, res, param = TRUE, pairs = TRUE, wy = TRUE, wz = TRUE, nit, thin, burnin)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim}{A list, typically the output of function \code{sim.weed}}
  \item{res}{A list, typically the output of function \code{pred.weed}}
  \item{param}{Logical: if \code{TRUE} graphs of simulated parameters contained in
    \code{res} are displayed}
  \item{pairs}{Logical: if \code{TRUE} pair plots of predicted weed
    values are ploted against known true values provided as outputof
    function \code{sim} are displayed (for debuging and convergence checking)}
  \item{wy}{Logical: if \code{TRUE} chains corresponding to the nine
    first sites in \code{y} are displayed}
  \item{wz}{Logical: if \code{TRUE} chains corresponding to the nine
    first sites in \code{z} are displayed}
  \item{nit}{Number of Markov iterations}
  \item{thin}{Thinning of the chain}
  \item{burnin}{Number of  values to discard before computing statistics
    about the empirical posterior. The number of iterations discarded is
    the number of stored iterations, e.g. if \code{nit=10000} and 
    \code{thin=10}, then only 1000 iterations are stored and
    \code{burnin=100} will remove the 100 first which correspond to
    actually  thin*burnin=1000 Markov iterations.}
}

\value{No value returned}

\references{G. Guillot, N. Loren, M. Rudemo, Bayesian spatial prediction of weed intensities from  exact count 
  data  and picture based indexes, 2006, submitted}
\author{Gilles Guillot}



\keyword{}

