\name{sqrtm_rs}
\alias{sqrtm_rs}
\title{ square root of a square matrix }
\description{
 Simple square root of a real symmetric matrix via eigendecomposition
}
\usage{
sqrtm_rs(x)
}
\arguments{
  \item{x}{REAL, [n,n], symmetric matrix}
}
\details{
  Simple square root of a real symmetric matrix via eigendecomposition
}
\value{
  A symmetric matrix
}
\examples{
 A <- matrix(c(4,1,1,4),nrow=2)
 Asq <- sqrtm_rs(A)
}
