\name{nnGpos}
\alias{nnGpos}
\title{nearest neighbour index from point matrix to grid}
\description{
 Find the nearest neighbour index to a grid from a matrix
 of point coordinates
}
\usage{
nnGpos(G, pos)
}
\arguments{
  \item{G}{a "gmeta6" object with grid metadata}
  \item{pos}{a two-column matrix with point coordinates}
}
\details{
The function assumes pos indicates [east,north] positions. The returned
value is a vector indicating, for each row in pos, the nearest index in
the grid refers to indexes sorted byrow=TRUE (C-like) starting by the
upper-left corner and advancing row by row, from North to South.
}
\value{
  A vector
}
\author{Javier Garcia-Pintado, e-mail:
  \email{j.garcia-pintado@reading.ac.uk}}
\keyword{ utilities }
\keyword{ spatial }
