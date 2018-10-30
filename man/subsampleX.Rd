\name{subsampleX}
\alias{subsampleX}
\title{Get a member from an X list of ensemble variables}
\description{
 Extract one member for all variable kinds and times in a X object. See
 details below.  
}
\usage{
subsampleX(X, ims)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A X object as described in the details below}
  \item{ims}{member index}
}
\details{
 X is list, where each component represents a state vector kind of
 variable (an 'xKIND'). Each X[[xKIND]] has the slot '$val', which is
 itself a list, with each component representing a time. These times
 match the slot 'X[[xKIND]]$timelab', and each 'it' time contains a
 matrix representing an ensemble state vector, each member in a column.
 The function cycles about variable kinds and times and extract one
 single member, preserving the X structure.
}
\value{
 A single-member X object
}
\author{
Javier Garcia-Pintado
}
