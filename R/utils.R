blkdiag <- function(...) {
  # +++ purpose +++
  # build a matrix by diagonally stacked block matrices
  dots <- list(...)
	if (! all(sapply(dots, is.matrix)) ||
	    ! all(sapply(dots, is.numeric)) )
		stop("All input arguments in '...' must be numeric matrices")

	nrows <- sapply(dots, nrow)
	ncols <- sapply(dots, ncol)
	if (any(nrows == 0) || any(ncols == 0))
		stop("All input matrices '...' must be non-empty.")

	n <- sum(nrows)
	N <- c(0, cumsum(nrows))
	m <- sum(ncols)
	M <- c(0, cumsum(ncols))
	A <- matrix(0, nrow = n, ncol = m)
	k <- length(dots)
	for (i in 1:k) {
		A[(N[i]+1):N[i+1], (M[i]+1):M[i+1]] <- dots[[i]]
	}
	return(A)
} # end function blkdiag

deg2gms <- function(deg) {
 # +++ purpose +++
 # convert decimal degrees to degrees minutes and seconds
 g <- trunc(deg)
 minrest <- (deg - g) * 60
 m <- trunc(minrest)
 s <- (minrest - m) * 60
 if (length(g) == 1){
   ans <- c(g,m,s)
   names(ans) <- c('g','m','s')
 } else {
   ans <- matrix(c(g,m,s),ncol=3)
   colnames(ans) <- c('g','m','s')
 }
 return(ans)
} # end function deg2gms

gms2deg <- function(gms) {
 # +++ purpose +++
 # convert degrees minutes and seconds to decimal degrees
 if (length(gms) == 3)
   gms <- matrix(gms,nrow=1)
 ans <- gms[,1] + gms[,2]/60 + gms[,3]/3600
 return(ans)
} # end function gms2deg

hypot <- function(a,b){
  sqrt(a^2 + b^2)
} # end function hypot

gridMax <- function(x, y, z, xmsk=rep(TRUE,length(x)), ymsk=rep(TRUE,length(y)), sFUN=max) {
  ng <- length(x) * length(y)
  if (length(z) != ng)
    stop('z length does not match given dimensions')
  if (length(xmsk) != length(x) && length(ymsk) != length(y))
    stop('Logical masks do not match given dimensions')
  if (!is.matrix(z))
    z <- matrix(z,length(x),length(y))
  x <- x[xmsk]                                              # block
  y <- y[ymsk]
  z <- z[xmsk,ymsk]
  ind <- which(z == sFUN(z), arr.ind=TRUE)
  ans <- c(x[ind[1]], y[ind[2]], z[ind[1], ind[2]])
  names(ans) <- c('x','y','z')
  return(ans)
} # end function gridMax()

#fooMax <- matStat(ol[[1]][['MOC']][['1850-01']][,1], lat_aux_grid, moc_z, moc_latmsk, moc_zmsk, sFUN=max)


wSigma <- function(sigma, w) {
  # inflates a covariance keeping correlation
  diag(sqrt(w)) %*% sigma %*%  diag(sqrt(w))   
} # end function wSigma()

