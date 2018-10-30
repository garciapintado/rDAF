\name{getEsensitivity}
\alias{getEsensitivity}
\title{ Get mean ensemble sensitivities }
\description{
 Get mean ensemble sensitivites of selected output in a model to the
 model inputs by multivariate regression
}
\usage{
getEsensitivity(X, Y, useEigen = FALSE, eigop = TRUE)
}
\arguments{
  \item{X}{REAL, [n,m] input matrix, one member in each column}
  \item{Y}{REAL, [p,m] output matrix, one member in each column}
  \item{useEigen}{LOGICAL, whether to froce eigenvalue decomposition to
    solve the system}
  \item{eigop}{LOGICAL, FALSE is not practical, only available to emulate  (A4) in Ancell and Hakim (2007)}
}
\details{
 Return a [p,n] sensitivity matrix \code{G} representing the mean ensemble
 sensitivity of the output to the input variables. The function uses either singular value
 decomposition or eigenvalue decomposition.
}
\value{
 A [p,n] sensitivity matrix
}
