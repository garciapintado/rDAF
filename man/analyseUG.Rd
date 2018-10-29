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
gauDA = NULL, gauTS = NULL, aug = NULL, ana = NULL, dsn = "", debugmode = TRUE,
retdy = FALSE, mpi = FALSE, analysis_scn = "a0", theta = NULL, lite=TRUE)
}
\arguments{
  \item{G}{OPT, LIST, grid metadata with a structure similar to \code{meta6}
    objects in package \code{spgrass6}} 
  \item{gLON}{OPT, LIST,    Geographical coordinates - unstructured grids [nx,ny] 'matrix' }
  \item{gLAT}{OPT, LIST,    Geographical coordinates - unstructured grids [nx,ny] 'matrix' }
  \item{fit}{INTEGER, forecast time index for saved output}
  \item{prm}{LIST,    DA-specific data}
  \item{X}{LIST,    state vector, structured as X[[xKIND]][[timelab]][[im]]}
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
  
For localization based on \code{SpatialGraph} objects, it is assumed that the
input \code{SpatialGraph} is stored as \code{file.path(dsn,prm$sgf)} and
precalculated along-network distances for nodes in the
\code{SpatialGraph} are stored in \code{file.path(dsn,prm$gridSGf)}.

However, for the current version, only code{assimilateLite} is
available, with no localisation or MPI parallelization available
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
\keyword{data assimilation}% use one of  RShowDoc("KEYWORDS")
\keyword{spatial}% __ONLY ONE__ keyword per line
