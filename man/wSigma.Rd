\name{wSigma}
\alias{wSigma}
\title{Inflate covariance matrix}
\description{inflate covariance matrix keeping correlations}
\usage{
wSigma(sigma, w)
}
\arguments{
  \item{sigma}{REAL, [n,n], covariance matrix}
  \item{w}{REAL, [n], inflation weights}
}
\examples{
P <- matrix(c(1.0, 0.5, 0.2,
              0.5, 2.0, 0.3,
              0.2, 0.3, 1.5), 3)
w <- c(1,2,3)
wP    <- wSigma(P,w)

Dsqi  <- diag(1/sqrt(diag(P)))
rho   <- Dsqi \%*\% P \%*\% Dsqi       # rho
wDsqi <- diag(1/sqrt(diag(wP)))
wrho  <- wDsqi \%*\% wP \%*\% wDsqi   # rho
all.equal(rho, wrho)
}
