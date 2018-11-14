approxM <- function(trmL, E, forward=TRUE) {
  # multivariate version of approx() linear interpolation
  #
  # trmL :: [[n]] list - transformation matrices : one transformation matrix / row in E
  # E    :: [n,m] matrix or [m] vector

  # notes:  
  # NA values are returned for xout out of input bounds

  if (class(E) != 'numeric')
    n <- nrow(E)
  else
    n <- 1

  if (length(trmL) != n)
    stop('ganam :: --ERR 001--')

  if (forward)
    trid <- 1:2
  else                                        # inverse
    trid <- 2:1
 
  tE <- E
  if (n == 1) { # scalar input
    tE <- approx(trmL[[1]][,trid[1]],y=trmL[[1]][,trid[2]],xout=E)$y 
    return(tE)
  }
  # else: vector-valued ensemble    
  for (i in 1:n) {
    tE[i,] <- approx(trmL[[i]][,trid[1]],y=trmL[[i]][,trid[2]],xout=E[i,])$y 
  }
  return(tE)
} # end function approxM
