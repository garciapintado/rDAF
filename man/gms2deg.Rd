\name{gms2deg}
\alias{gms2deg}
\title{Convert degrees, minutes and seconds to decimal degrees}
\description{
Convert degrees, minutes and seconds to decimal degrees
}
\usage{
gms2deg(gms)
}
\arguments{
  \item{gms}{a length-3 vector or a 3-columns matrix expressing sexagesimal degrees, minutes
    and seconds}
}
\value{
The function returns a numeric vector of decimal degrees
}
\examples{
deg <- c(12.3,124.5,358.1)
gms <- deg2gms(deg)
deg2 <- gms2deg(gms)
}
\keyword{ utilities }
\keyword{ spatial }
