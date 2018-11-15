\name{deg2gms}
\alias{deg2gms}
\title{Convert decimal degrees to degrees, minutes and seconds}
\description{
Convert decimal degrees to degrees, minutes and seconds
}
\usage{
deg2gms(deg)
}
\arguments{
  \item{deg}{a scalar or vector of decimal degrees}
}
\value{
if input is a scalar the funtions returns a numeric vector, if input is
a vector of decimal degrees, the function returns a 3-columns matrix,
with names "g","m","s".
}
\examples{
deg <- c(12.3,124.5,358.1)
gms <- deg2gms(deg)
}
\keyword{ utilities }
\keyword{ spatial }
