\name{tlag}
\alias{tlag}
\title{ Add a lag to POSIXct time vector }
\description{
Add a positive or negative lag to a POSIXct timestamp vector
}
\usage{
tlag(x, dlag = NULL, mlag = NULL, ylag = NULL)
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
A POSIXct time vector
}
