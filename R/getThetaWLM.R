getThetaWLM <- function(theta, draw=FALSE, weighted=TRUE) {
  # this function is a simple wrapper around a weighted linear regression
  
  # theta$thetaNam : input names   [ntheta]
  # theta$b        : background    [>= ntheta]
  # theta$d        : perturbations [ntheta,mperpar]
  # theta$cost     : weights^{-1}  [m = 1 + mperpar*ntheta]
  # theta$HE       : dual of the observation space [p,m]

  ntheta   <- nrow(theta$d)
  mperpar  <- ncol(theta$d)
  thetaNam <- rownames(theta$d)

  m <- 1 + mperpar * ntheta
  p <- nrow(theta$HE)

  if (weighted)
    rw <- 1/theta$cost                                                # weights for linear regression
  else
    rw <- rep(1,m)

  G <- matrix(NA, p, ntheta)
  for (i in 1:ntheta) {
    vname <- thetaNam[i]
    cat('getThetaWLM: ',vname,'\n')
    immin <- (i-1)*mperpar + 1 + 1                                    # xtra 1 for the background
    immax <- immin + mperpar - 1
    ims <- c(1,immin:immax)                                           # [1+mperpar]
    x <- c(theta$b[vname],  theta$b[vname] + theta$d[i,])             # [1+mperpar]
    for (iy in 1:p) {
      HE <- theta$HE[iy,ims]
      lm0 <- lm(HE~x,weights=rw[ims])
      G[iy,i] <- coef(lm0)[2]
      if (draw) {
        cat('getThetaWLM:: drawing is deactivated in the packaged version')  
        #dev.new()
        #plot(x, HE, xlim=range(x,na.rm=TRUE), ylim=range(HE, na.rm=TRUE))
        #abline(lm0, col='blue')
      }
    }
  }
  return(G)
}

