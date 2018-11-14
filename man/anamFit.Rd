\name{anamFit}
\alias{anamFit}
\title{ Fit an empirical Gaussian Anamorphosis}
\description{Fit an empirical Gaussian Anamorphosis}
\usage{
anamFit(x, extrapolate=TRUE, ndisc=101, useRGeostats=FALSE,
        xlim=c(-Inf,Inf), keepMoments=TRUE)
}
\arguments{
  \item{x}{REAL, input vector}
  \item{extrapolate}{LOGICAL, whether to go beyond input limits by
    multiplying the original range by two}
  \item{ndisc}{INTEGER, number of steps used to discretize the range of raw
    values, This parameter is used only for empirical anamorphosis and
    if the anamorphosis is done via \code{anam.fit} in package RGeostats }
  \item{useRGeostats}{LOGICAL, where to call \code{RGeostats::anam.fit}}
  \item{xlim}{REAL, bounds for the raw variables. For example, \code{xlim[1]==0}
    prevents negative values in the back-transformed raw values}
  \item{keepMoments}{LOGICAL, whether to preserve the mean and standard deviation
    in the transformed variable vector with respect to the raw}
}
\examples{
 xraw <- rlnorm(100)
 xGA <- anamFit(xraw)
 xgau <- approx(xGA[,'r'], xGA[,'g'], xraw)$y
 #hist(xraw, nclass=25)
 #hist(xgau, nclass=25)
}