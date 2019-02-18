\name{getObjects}
\alias{getObjects}
\title{Get all objects in the calling function as a list}
\description{
 This function is houthgt to be called from another function, from which
 all object in the local environent are returned as a named list
}
\usage{
getObjects(except)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{except}{CHARACTER vector indicating objects to exclude}
}
\value{
 A list
}
\examples{
 # not run
 if (1 > 2) {
   anyfun <- function(model) {
     arg <- ls()
     EXEROOT <- 'loquesea'
     FOO     <- 1:10
     getObjects(c('arg',arg))
   }
   xx <- anyfun('model1')
 } # end not run
}

