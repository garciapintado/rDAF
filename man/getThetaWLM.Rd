\name{getThetaWLM}
\alias{getThetaWLM}
\title{ Wrapper around univariate weighted linear regression for a theta
  LIST }
\description{ A theta LIST has a specific structure described in the vignette in
 this package. This function is a simple wrapper to automate weighted
 linear regressions for input and output organised in a \code{theta} LIST.
}
\usage{
getThetaWLM(theta, draw = FALSE, weighted = TRUE)
}
\arguments{
  \item{theta}{A list with at least the components \code{theta$b} , which is a
    named vector of perturbed model inputs, \code{theta$cost}, which serves
    to screen out failed ensemble simulations (for which the cost
    function is NA) and as possible weights, \code{theta$MCgpar}, or ensemble of model inputs,
    and \code{theta$HE}, or dual of the observation space}
  \item{draw}{LOGICAL, draw linear regression on the screen as they are calculated}
  \item{weighted}{LOGICAL, if TRUE, weights are inversely proportional
    to the value of the cost function provided in \code{theta$cost}}
}
\details{
 Univariate regression betwwen the observations and the individual model inputs
}
