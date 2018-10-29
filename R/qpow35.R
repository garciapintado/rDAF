qpow35 <- function(x, inverse=FALSE) {
  x[x < 0.0] <- 0.0
  if (!inverse) {
    x <- x^(3.0/5.0)
  } else {
    x <- x^(5.0/3.0)
  }
  return(x)
}
