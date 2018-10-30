\name{getThetaG}
\alias{getThetaG}
\title{
 Wrapper around getEsensitivity for theta LISTS
}
\description{
 A theta LIST has a specific structure described in the vignette in
 this package. Two of the components of theta are \code{theta$MCgpar},
 which is a matrix of random model inputs, and \code{theta$HE}, which is
 the model background mapped into the dual of the observation space. The
 function gets the ensemble sensitivity of \code{theta$HE} to \code{theta$MCgpar}.
}
\usage{
getThetaG(theta, useEigen = NULL)
}
\arguments{
  \item{theta}{A list with at least the components \code{theta$b} , which is a
    named vector of perturbed model inputs, \code{theta$cost}, which serves
    just to screen out failed ensemble simulations (for which the cost
    function is NA), \code{theta$MCgpar}, or ensemble of model inputs,
    and \code{theta$HE}, or dual of the observation space}
  \item{useEigen}{argument passed to \code{\link{getEsensitivity}}
}
}
\details{
 The ensemble sensitivity is obtained
}
\value{
  A [p,ntheta] sensitivity matrix
}
\author{
 Javier Garcia-Pintado
}
