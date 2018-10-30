\name{analyseUG}
\alias{analyseUG}
\title{
 ensemble contructor as input to assimilation
}
\description{
 Construct the matrices and objects needed as input to
 \code{assimilate}. Structured grid input metadata are given via the
 \code{G} argument. Unstructured grid input metadata are given via
 \code{gLON} and \code{gLAT}.
}
\usage{
analyseUG(G = NULL, gLON = NULL, gLAT = NULL, fit, prm, X, yls = NULL,
gauDA = NULL, gauTS = NULL, aug = NULL, ana = NULL, dsn = ".", debugmode = TRUE,
retdy = FALSE, mpi = FALSE, analysis_scn = "a0", theta = NULL, lite=TRUE)
}
\arguments{
  \item{G}{OPT, LIST, grid metadata for structured grids with a
    format as \code{gmeta6} objects in package spgrass6} 
  \item{gLON}{OPT, LIST, geographical coordinates - unstructured grids [nx,ny] 'matrix' }
  \item{gLAT}{OPT, LIST, geographical coordinates - unstructured grids [nx,ny] 'matrix' }
  \item{fit}{INTEGER, forecast time index for saved output}
  \item{prm}{LIST, global assimilation parameters}
  \item{X}{LIST, state vector, structured as X[[xKIND]][[timelab]][[im]]}
  \item{yls}{OPT, LIST,    observation list, structured as yls[[yKIND]][[yTYPE]]$s lists.}
  \item{gauDA}{OPT, LIST,    y <-> HE matches, to include observations by state-vector augmentation. Augments x and y}
  \item{gauTS}{OPT, LIST, as gauDA lists, without the $y slot. Augments x, and can be used to obtain a smooth estimate (ETKS) at specific locations/variables in the domain. }
  \item{aug}{OPT, LIST,    augmentation blocks (e.g.; physical parameters) }
  \item{ana}{LIST,    analysis parameter for each variable type}
  \item{dsn}{CHARACTER != '', output files are stored in the path $dsn/results}
  \item{debugmode}{OPT, LOGICAL, TRUE}
  \item{retdy}{OPT, LOGICAL, default FALSE. TRUE for no assimilation, and instead return the innovation vector}
  \item{mpi}{OPT, LOGICAL, (forced to FALSE in rDAFLite) whether to parallelize the assimilation via MPI}
  \item{analysis_scn}{CHARACTER, label added to output file for eventual
    identification of assimilation scenarios}
  \item{theta}{OPT, LIST, required for finite difference sensitivity
    Kalman smoothers [pIKS,pMKS]}
  \item{lite}{LOGICAL, only TRUE in this version}
}
\details{
A detailed documentation of the arguments to this function is under
preparation as a vignette. Specifically, this involves \code{gauDA},
\code{prm}, \code{aug} and \code{ana} lists. Please, ask the
author/maintainer for documentation requests in the meantime.

\code}{analyseUG} conducts the preparation of the matrices needed for
the assimilation. The actual assimilation is conducted by the function
\code{assimilate}, which includes Euclidean based localization and
localization based on spatial networks, via the use of the SpatialGraph
package, as well as MPI parallelization.

For localization based on \code{SpatialGraph} objects, it is assumed that the
input \code{SpatialGraph} is stored as \code{file.path(dsn,prm$sgf)} and
precalculated along-network distances for nodes in the
\code{SpatialGraph} are stored in \code{file.path(dsn,prm$gridSGf)}.

The function \code{assimilate} is planned to be publicly released along
with the vignette documentation. In the interim time, the version
\code{assimilateLite}, which conducts global assimilation via the EnKF
or the ETKF, serves for the purpose of package
testing. \code}{analyseUG} also includes the Finite Difference
Sensitivity (FDS) smoothers described in García-Pintado and Paul (2018).
}
\value{
A LIST with the components
  \item{E}{Updated ensemble matrix}
  \item{xdf}{Data frame with state vector information}
  \item{dydf}{Data frame with observation and innovation information}
}
\author{
Javier Garcia-Pintado
}
\references{
García-Pintado J. and Paul, A., 2018. Evaluation of iterative Kalman
smoother schemes for multi-decadal past climate analysis with
comprehensive Earth system models. Geoscientific Model Development. 
}
\keyword{data assimilation}% use one of  RShowDoc("KEYWORDS")
\keyword{spatial}
