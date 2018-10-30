orthoEnsemble <- function(p, MC, m, MCcons=NULL, sdfac=1.0) {
  # makes a conditional-sampling parameter Ensemble given the parameters in p & MC
  # this ensemble includes the control run [or background] as $b,
  # then one list member for each perturbed parameter
  
  # ensemble includes all MC$ispar=TRUE parameters in each parameter-space state vector
  # 
  # MCcons gives a number of equality and inequality constrains in a 3-column matrix, with one row per constraint
  # as it might happen that a background value violates a constrain, these are first evaluated on the background, then on the generated ensemble

  if (!all.equal(names(p),rownames(MC)))
    stop('orthoEnsemble:: mismatch between p and MC names')

  pboo <- MC$ispar                               # LOGICALq
  pnames   <- rownames(MC)[pboo] 
  thetaNam <- rownames(MC)[pboo & MC$flag]       # names of just uncertain parameters
  ntheta <- length(thetaNam)                     # number of uncertain parameters

  if (length(sdfac) == 1)
   sdfac <- rep(sdfac,ntheta)
  if (length(sdfac) != ntheta)
   stop('orthoEnsemble:: ---ERR001---')

  if (m == 1) {
    mperpar <- 0
  } else {
    if (m %% ntheta != 0 && m > 1)
      stop('orthoEnsemble:: ---ERR002---')
    mperpar <- m/ntheta
  }

  gpar     <- list()
  gpar$b   <- sapply(p[pboo],function(x){x[1]})

  # check b is within constraints
  for (i in 1:length(thetaNam)) {
    vnam <- thetaNam[i]
    gpar$b[vnam] <- max(gpar$b[vnam], MC[vnam,'min'])
    gpar$b[vnam] <- min(gpar$b[vnam], MC[vnam,'max'])
  }

  gpar$thetaNam <- thetaNam

  if (m ==1)
    return(gpar)

  dtheta <- matrix(NA,ntheta,mperpar)
  rownames(dtheta) <- thetaNam

  for (i in 1:length(thetaNam)) {
    vnam <- thetaNam[i]
    cat('drawing sample from: ',vnam,'\n')
    gpar[[vnam]] <- matrix(gpar$b, nrow=sum(pboo), ncol=mperpar) # recycle background in each column
    rownames(gpar[[vnam]]) <- pnames

    if (MC[vnam,'dis'] %in% c('rnorm','truncnorm'))
      sdev <- p[[vnam]][2]
    else if (MC[vnam,'dis'] == 'rlnorm')
      sdev <- p[[vnam]][1]*0.5                    # assume cv=0.5 for lognormal variables

    
    if (mperpar == 1) {                           # finite difference approximation
       dtheta[i,] <- sdev*sdfac[i]                # sd(c(b,dtheta)) == sdev*sqrt(2)/2
    } else {
       dtheta[i,] <- rnorm(mperpar, mean=0, sd=sdev*sdfac[i])       
    }
    gpar[[vnam]][vnam,] <- gpar[[vnam]][vnam,] + dtheta[i,]  # perturbed ensemble 
  }


  # add fix constraints
  for (i in 1:length(thetaNam)) {
    vnam <- thetaNam[i]
    if (!is.na(MC[vnam,'min']))
     gpar[[vnam]][vnam,] <- pmax(gpar[[vnam]][vnam,], MC[vnam,'min'])
    if (!is.na(MC[vnam,'max']))
     gpar[[vnam]][vnam,] <- pmin(gpar[[vnam]][vnam,], MC[vnam,'max'])
  }

  # add dependence constrains
  if (!is.null(MCcons)) { # the 2nd term constraints the 1st 
    for (ic in 1:nrow(MCcons)) {
      vnam1 <- MCcons[ic,1]
      vnam2 <- MCcons[ic,3]
      # $b
      if (MCcons[ic,2] == '<')
        gpar$b[vnam1] <- min(gpar$b[vnam1],gpar$b[vnam2])
      if (MCcons[ic,2] == '>')
        gpar$b[vnam1] <- max(gpar$b[vnam1],gpar$b[vnam2])
      if (MCcons[ic,2] == '=')
        gpar$b[vnam1] <- gpar$b[vnam2]
      # [[thetaNam]]
      for (i in 1:length(thetaNam)) {
        lnam <- thetaNam[i]
        if (MCcons[ic,2] == '<')
          gpar[[lnam]][vnam1,] <- pmin(gpar[[lnam]][vnam1,], gpar[[lnam]][vnam2,])
        if (MCcons[ic,2] == '>')
          gpar[[lnam]][vnam1,] <- pmax(gpar[[lnam]][vnam1,], gpar[[lnam]][vnam2,])
        if (MCcons[ic,2] == '=')
          gpar[[lnam]][vnam1,] <- gpar[[lnam]][vnam2,]
      }
    }
  }
  gpar$dtheta   <- dtheta

  for (i in 1:length(thetaNam)) {
    vnam <- thetaNam[i]
    gpar$dtheta[vnam,] <- gpar[[vnam]][vnam,] - gpar$b[vnam]
  }

  return(gpar)
}
