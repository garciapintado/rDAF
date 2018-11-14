\name{interppEarth}
\alias{interppEarth}
\title{Interpolates spatial locations in the globe}
\description{
  Interpolates spatial locations from irregular points into irregular points in
  the globe by use of akima interpp
}
\usage{
interppEarth(LON, LAT, z, xo, yo)
}
\arguments{
  \item{LON}{vector of input longitudes in [0,360]}
  \item{LAT}{vector of input latitudes in [-90,90]}
  \item{z}{data values}
  \item{xo}{vector of output longitudes}
  \item{yo}{vector of output latitudes}
}
\details{
  \code{interppEarth} stacks longitudes bands in input data as approximate continuous
  longitudinal domain}
\value{
A REAL vector of interpolated values at requested locations
}
\keyword{ utilities }
\keyword{ spatial }
