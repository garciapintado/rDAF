\name{orthoEnsemble}
\alias{orthoEnsemble}
\title{ Conditional-sampling ensemble subject to constraints}
\description{
 Given a background mean and standard deviation for a number of variables, the function returns random samples
 from the background conditional distribution, with a given inflation
 factor and possibly subject to fixed-value and cross-variable constraints. See details below.
}
\usage{
orthoEnsemble(p, MC, m, MCcons = NULL, sdfac = 1)
}
\arguments{
  \item{p}{LIST with named input variables as components. Each component
    is then a vector of length two, containing its background mean and standard deviation
  }
  \item{MC}{Dataframe with metadata, where the \code{MC[,'ispar']} is a LOGICAL
    column indicating which variables are to be included in the output
    matrix, out of all those available in \code{p}. Output variables
    are those with both  \code{MC[,'ispar']} and  \code{MC[,'flag']} set
    to \code{TRUE}, with length \code{ntheta}. \code{MC} has the additional 'max' and 'min' columns
    indicating fix-valued contraints on output samples.
    \code{MC[,'dis']} indicates the random distribution type. Available values are 'rnorm',
   'truncnorm' and 'rlnorm'.}
  \item{m}{number of members in the ensemble, m has to be an exact
    multiplier of ntheta}
  \item{MCcons}{CHARACTER, 3-column matrix of cross-variable constraints. For each row, the first and third
  columns are variable names matching those in \code{p}, and the second
  column is one character in ['<','=','>']}
  \item{sdfac}{REAL, scalar or [ntheta] vector. Multiplicative factor of the standard deviation for the perturbations}
}
\details{
  Makes a conditional-sampling parameter Ensemble given the parameters
  in p and MC. This ensemble includes the control run [or background]
  as the component $b. Then one list member for each perturbed
  parameter.
  MCcons gives a number of equality and inequality constrains with one row per constraint. As it might happen that a background value violates a constrain, these are first evaluated on the background, then on the generated ensemble
}
\value{
A LIST with the components
  \item{b}{[ntheta] background vector}
  \item{thetaNam}{[ntheta] name of the output variables}
  \item{dtheta}{[ntheta,mperpar] matrix of perturbations, where mperpar
    is the number of perturbations per thetaNam variable}
  Then for each thetaNam, there is an additional component in the LIST
  with a vector or matrix (if more than one sample is requested for each
  thetaNam variable) where only thetaNam is perturbed
}
\author{
Javier Garcia-Pintado
}
