\name{assimilateLite}
\alias{assimilateLite}
\title{ Ensemble Assimilation Lite version}
\description{
 Main assimilation function called by \code{analyseUG}
}
\usage{
assimilateLite(prm, A, HA, dy, R, debugmode, Kfname,
               xpos = NULL, xllen = NULL, ypos = NULL, mpi = NULL, sglst = NULL)
}

\arguments{
  \item{prm}{LIST, assimilation parameters}
  \item{A}{REAL, [n,m], 'dgCMatrix' or 'dgeMatrix', perturbation matrix,
  where \code{n}: augmented state vector length, and \code{m}: ensemble members}
  \item{HA}{REAL, [p,m], perturbations mapped into the dual of the observation
  space, where \code{p} is the number of observations}
  \item{dy}{REAL, [p], innovation vector}
  \item{R}{REAL, [p,p], 'ddiMatrix', 'dsyMatrix', or 'dtCMatrix': observation error covariance}
  \item{debugmode}{LOGICAL, whether to activate verbose output}
  \item{Kfname}{ CHARACTER, filename, including path, of the output file
    to write the Kalman gain to}
  \item{xpos}{REAL, [n,2], augmented state geographical positions [0,0] for global variables}
  \item{xllen}{REAL, [n], localisation length}
  \item{ypos}{REAL, [p,2], observation coordinates}
  \item{mpi}{LOGICAL}
  \item{sglst}{LIST. If not NULL, SpatialGraph-based distance is used for localisation}
}
\details{
This function would most commonly be called by
  \code{\link{analyseUG}}. The arguments \code{xpos}, \code{xllen},
  \code{ypos}, \code{mpi}, and \code{sglst} are only proided for input
  compatibility with \code{assimilate()}, which is under
  development. These are oriented to conduct observation localization, and
  are unused and set to NULL within this function. 
}
\value{
A LIST with the components
  \item{dx}{vector of increments}
  \item{A}{matrix of updated perturbations}
}
\author{
Javier Garcia-Pintado
}
