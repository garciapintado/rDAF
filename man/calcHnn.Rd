\name{calcHnn}
\alias{calcHnn}
\title{ grid to point nearest neighbours }
\description{
 Fill in a LOGICAL matrix indicating for each location in a [p] set the
 corresponding nearest neighbour in a [ng] grid.
}
\usage{
calcHnn(G, gLON, gLAT, pos)
}

\arguments{
  \item{G}{LIST, a code{gmeta6} object as defined by the package
  spgrass6 with grid metadata for regular grids}
  \item{gLON}{REAL, [ng] coordinates for irregular grids, longitude}
  \item{gLAT}{REAL, [ng] coordinates for irregular grids, latitude}
  \item{pos}{REAL, [p,2], point coordinates}
}
\details{
  This function would most commonly be called by \code{\link{analyseUG}}.
  For irregular grids, nearest grid neighbours to the point dataset is
  calculated via \code{nn2dist} in the package splancs. 
}
\value{
A [p,ng] \code{dgCMatrix} as provided by the package Matrix, where
  \code{p} is the number of point locations, and \code{ng} is the number
  of cells in the grid. The grid is assumed to be sortened starting fro
  the upper left and advancing by row. 
}
