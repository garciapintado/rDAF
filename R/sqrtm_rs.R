sqrtm_rs <- function(x) {
  # square root of a real symmetric matric by diagonalisation
  # warning: no check is done for symmetry. It is up to the user to call this function appropriately
  CL <- eigen(x, symmetric=TRUE) # C=CL$vectors will be then an orthogonal matrix, so C^{-1} == t(C) 
  #C  <- CL$vectors
  #L  <- CL$values
  #xsq <- C %*% diag(sqrt(L)) %*% t(C)
  xsq <- tcrossprod(CL$vectors %*% diag(sqrt(CL$values)), CL$vectors)
  return(xsq)
}
