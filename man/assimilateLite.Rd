\name{assimilateLite}
\alias{assimilateLite}
\title{ Ensemble Assimilation Lite version}
\description{
 Main assimilation function called by \code{analyseUG}
}
\usage{
assimilateLite(prm, A, HA, dy, R, debugmode, Kfname)
}

\arguments{
  \item{prm}{LIST, assimilation parameters}
  \item{A}{REAL, [n,m], 'dgCMatrix' or 'dgeMatrix', perturbation matrix,
  where \code{n}: augmented state vector length, and \code{m}: ensemble members}
  \item{HA}{REAL, [p,m], perturbations mapped into the dual of the pbservation
  space, where \code{p}: number of observations}
  \item{dy}{REAL, [p], innovation vector}
  \item{R}{REAL, [p,p], 'ddiMatrix', 'dsyMatrix', or 'dtCMatrix': observation error covariance}
  \item{debugmode}{LOGICAL, whether to activate verbose output}
  \item{Kfname}{ CHARACTER, filename, including path, of the output file
  to write the Kalman gain to}
}
\details{
This function would most commonly be called by \code{\link{analyseUG}}
}
\value{
A LIST with the components
  \item{dx}{vector of increments}
  \item{A}{matrix of updated perturbations}
}
\author{
Javier Garcia-Pintado
}
