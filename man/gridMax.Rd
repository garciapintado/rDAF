\name{gridMax}
\alias{gridMax}
\title{Maximum value of a grid block}
\description{ Get the maximum (or minimum) value and corresponding
  coordinates from a grid block, defined from logical masks for the grid dimensions}
\usage{ gridMax(x, y, z, xmsk, ymsk, sFUN=max) }
\arguments{
  \item{x}{coordinates of first dimension}
  \item{y}{coordinates of second dimension}
  \item{z}{numeric matrix, or vector matching x,y dimensions, advancing first by x indices}
  \item{xmsk}{LOGICAL vector}
  \item{ymsk}{LOGICAL vector}
  \item{sFUN}{either 'max' (default), or 'min'}
}
\details{
 If more than a maximum/minimum value is found, the first one is returned
}
\value{
A named vector of length 3 with x,y,z values.
}
