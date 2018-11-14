\name{approxM}
\alias{approxM}
\title{ Multivariate linear interpolation}
\description{ Conduct a multivariate interpolation by calling approx() }
\usage{ approxM(trmL, E, forward=TRUE) }
\arguments{
  \item{trmL}{[[n]] LIST - transformation matrices : one transformation matrix / row in E}
  \item{E}{[n,m] matrix or [m] vector}
  \item{forward}{LOGICAL, default to TRUE}
}
\details{
 NA values are returned for xout out of input bounds
}
\value{
A vector or matrix of interpolated values
}
