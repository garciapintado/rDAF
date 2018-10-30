\name{mcEnsemble}
\alias{mcEnsemble}
\title{
 Random ensemble subject to constraints
}
\description{
  Given a background mean and standard deviation for a number of
  variables, the function returns a random ensemble possibly subject to
  fixed-value and cross-variable constraints. See details below.
}
\usage{
mcEnsemble(p, MC, m, MCcons = NULL, Sigma = NULL, empiricalNorm = TRUE)
}
\arguments{
  \item{p}{LIST with named input variables as components. Each component
    is then a vector of length two, containing its background mean and standard deviation
  }
  \item{MC}{Dataframe with metadata, where the \code{MC[,'ispar']} is a LOGICAL
    column indicating which variables are to be included in the output
    matrix, out of all those available in \code{p}}. The number of
  output variables in the ensemble is \code{ntheta==sum(CF$MC$flag)}
  \code{MC} has the additional 'max' and 'min' columns
  indicating fix-valued contraints on output samples.
  \code{MC[,'dis']} indicates the random distribution type. Available values are 'rnorm',
  'truncnorm' and 'rlnorm'.
  \item{m}{number of members in the ensemble, m has to be an exact
    multiplier of ntheta}
  \item{MCcons}{CHARACTER, 3-column matrix of cross-variable constraints. For each row, the first and third
  columns are variable names matching those in \code{p}, and the second
  column is one character in ['<','=','>']}
  \item{Sigma}{REAL, OPT, [ntheta,ntheta], full covariance matrix. If
    Sigma is provided, only those element in Sigma with normal distribution will be considered for covariances}
  \item{empiricalNorm}{whether to use empircal norm for those variables
    with an 'rnorm' distribution}
}
\details{
 Makes a random sample subject to the abovementioned fixed-value
 constraints and cross-variable constraints.
}
\value{
An [ntheta,m] matrix
}
