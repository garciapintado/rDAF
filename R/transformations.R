anamFit <- function(x, extrapolate=TRUE, ndisc=101, useRGeostats=FALSE, xlim=c(-Inf,Inf), keepMoments=TRUE) {
  # fit an empirical Gaussian Anamorphosis transformation
  # This can done through RGeostats::anam.fit with empirical Gaussian Anamorphosis,
  # but an option is available if this package is no available
  #  
  # return matrix with two colums for eventual linear interpolation with stats:approx()
  # 
  # xlim are soft limits. If these are exceeded by bound in the estimated transform function, they are disregarded  
  # and a warning is issued
 
  x   <- x[!is.na(x)]
  n   <- length(x)
  xsde <- sd(x)                                                    # prior moments
  xbar <- mean(x)

  if (xsde == 0)
    return(NULL)
  if (!('RGeostats' %in% installed.packages())) {
    if (useRGeostats) {
      cat("anamFit: RGeostats is not installed. useRGeostats forced to FALSE\n")
      useRGeostats <- FALSE
    }
    db.create <- function(x){x}
    anam.fit <- function(x,ndisc) {
      usup <- seq(0,1,length=ndisc)
      x2U  <- cbind(quantile(x, probs=usup, na.rm=TRUE, names=FALSE, type=8),
                    usup)                                                       # P_x(x)
      U2N  <- cbind(qnorm(usup),usup)                          
      rg   <- cbind(x2U[,1],U2N[,1])                                            # P^{-1}_{\tilde{x}}\(  P_x(x) 
      if (rg[1,2] == -Inf) {
        slp <- diff(rg[2:3,2])/diff(rg[2:3,1])
        rg[1,2] <- rg[2,2] - slp * diff(rg[1:2,1])
      }
      if (rg[ndisc,2] == Inf) {
        slp <- diff(rg[ndisc - c(2,1), 2])/diff(rg[ndisc - c(2,1),1])
        rg[ndisc,2] <- rg[ndisc-1,2] + slp * diff(tail(rg[,1],2))
      }
      return(rg) 
    } # end function anam.fit
  } # end local db.create() and anam.fit()
  
  if (useRGeostats) {                                     # importFrom("RGeostats","db.create","anam.fit")  
    #require('RGeostats')
    x0 <- min(0, x, na.rm=TRUE)
    db   <- db.create(list(x1=1:n))                                # RGeostat::db.create
    db@items$x1 <- x - x0                                          # empirical-type anamorphosis in RGeostats does not allow for negative numbers
    anam <-  anam.fit(db, name='x1', type='emp', draw=FALSE, ndisc=ndisc) # RGeostat::anam.fit
    rg   <-  matrix(anam@tdisc,nrow=anam@ndisc)                    # anam@tdisc: (empirical) discretization array
    rg   <- rg[,c(2,1)]                                            # [raw,Gaussian] matrix | plot(rg[,1],rg[,2])
    if (x0 < 0)                                                    # add back x0 if needed
      rg[,1] <- rg[,1] + x0
  } else {
    rg <- anam.fit(x,ndisc)
  }

  if (xlim[1] > rg[1,1]) {
    cat("anamFit:: warning! xlim[1] > minimum raw value in the transformation. Setting the latter as limit\n")
    xlim[1] <- rg[1,1]
  }
  if (xlim[2] < rg[nrow(rg),1]) {
    cat("anamFit:: warning! xlim[2] < maximum raw value in the transformation. Setting the latter as limit\n")
    xlim[2] <- rg[nrow(rg),1]
  }

  if (extrapolate) {                                                  # extent in both directions
    incg <- diff(range(rg[,2]))*2                               
    slp <- diff(rg[1:2,2])/diff(rg[1:2,1])                            # first segment slope              
    rg <- rbind(rg[1,]-c(incg/slp,incg),rg)                
    if (rg[1,1] < xlim[1])                                            # e.g. xlim[1]==0 to prevent negative values in original variable
      rg[1,1] <- xlim[1]
    slp <- diff(tail(rg[,2],2))/diff(tail(rg[,1],2))                  # last segment slope              
    rg <- rbind(rg, tail(rg,1)+c(incg/slp,incg))                
    if (rg[nrow(rg),1] > xlim[2])
      rg[nrow(rg),1] <- xlim[2]
    #  incr <- diff(range(rg[,1])) 
    #  incg <- diff(range(rg[,2]))
    #  rg <- rbind(rg[1,]-c(incr,incg),rg)
    #  rg <- rbind(rg, tail(rg,1)+c(incr,incg))
  }
  colnames(rg) <- c('r','g')                                          # raw, Gaussian
  rownames(rg) <- NULL 

  # preserve two first statistical moments
  if (keepMoments) {
    xg <- approxM(list(rg),x)
    rg[,2] <- (rg[,2] - mean(xg)) / sd(xg) * xsde + xbar
    xg <- approxM(list(rg),x)
    if (!all.equal(c(mean(xg),sd(xg)),c(xbar,xsde)) )
      stop('anamFit: moment preservation failed')
  }
  return(rg)
} # end function anamFit()

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
} # end function linearScale()

qpow35 <- function(x, inverse=FALSE) {
  x[x < 0.0] <- 0.0
  if (!inverse) {
    x <- x^(3.0/5.0)
  } else {
    x <- x^(5.0/3.0)
  }
  return(x)
} # end function qpow35
