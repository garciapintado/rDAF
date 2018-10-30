getThetaG <- function(theta, useEigen=NULL) {
  thetaNam <- names(theta$b)
  sm <- !is.na(theta$cost)
  X <- theta$MCgpar[thetaNam,sm]
  Y <- theta$HE[,sm]
  getEsensitivity(X, Y, useEigen=useEigen)
}
