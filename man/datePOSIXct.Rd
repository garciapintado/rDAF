\name{datePOSIXct}
\alias{datePOSIXct}
\title{
 Get date from POSIXct with lags
}
\description{
 Get dates from a POSIXct time vector with a possible added lag
}
\usage{
datePOSIXct(x, dlag = NULL, mlag = NULL, ylag = NULL)
}
\arguments{
  \item{x}{A POSIXct time vector}
  \item{dlag}{lag: days}
  \item{mlag}{lag: months}
  \item{ylag}{lag: years}
}
\details{
 Any combination of days months and years is possible as input
}
\value{
A CHARACTER vector of dates in "%Y-%m-%d" format
}
\seealso{
See Also \code{\link{tlag}}
}
