linearScale <- function(x, a=0.0, b=1.0, inverse=FALSE) {
  # a, b :: scalar, or length(ncol(x)) vectors or dimensions matching x 
  # This is so, to allow for time augmentation of the parameter in analysis()
  
  if (length(a) != 1 && length(a) != ncol(x) && length(a) != nrow(x)*ncol(x))
    stop('linearScale:: a length not exact divisor of x length')
  if (length(a) == ncol(x))
    a <- rep(a, each=nrow(x))          # expand [m -> m*n]

  if (length(b) != 1 && length(b) != ncol(x) && length(b) != nrow(x)*ncol(x))
    stop('linearScale:: b length not exact divisor of x length')
  if (length(b) == ncol(x))
    b <- rep(b, each=nrow(x))          # expand [m -> m*n]
  
  if (!inverse) {         # forward
    x <- a + b * x
  } else {                # inverse
    x <- (x - a) / b
  }
  return(x)
}
