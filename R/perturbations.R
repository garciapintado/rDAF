mcEnsemble <- function(p, MC, m, MCcons=NULL, Sigma=NULL, empiricalNorm=TRUE) {
  # make a Monte Carlo Ensemble given the parameters in p & MC

  # if Sigma is provided, it has to be [ntheta,ntheta] dimensions, where ntheta==sum(CF$MC$flag)
  # only those element in Ptheta with normal distribution will be considered for covariances

  #if (empiricalNorm)
  #  require(MASS)
  #if (any(MC[,'dis'] == 'truncnorm'))
  #  require(truncnorm)
  
  if (!all.equal(names(p),rownames(MC)))
    stop('mcEnsemble:: missmatch of names between p and MC')

  pboo <- MC$ispar                              # LOGICAL
  if (!all.equal(pboo,MC$flag))
    stop('mcEnsemble:: MC$ispar != MC$flag')
  
  gpar <-  matrix(NA, nrow=sum(pboo), ncol=m)   # standard DA ensemble deployment
  rownames(gpar) <- rownames(MC)[pboo]
  for (i in 1:nrow(gpar)) {
    vnam <- rownames(gpar)[i]
    if (is.character(p[[vnam]])) {
      vvec <- rep(NA, m)  
    } else {
      if (length(p[[vnam]]) > 2)
        stop('randomisation just available for numeric scalars')
      if (MC[vnam,'flag'] && m > 1) { # numeric random
        if (is.logical(p[[vnam]]))
          vvec <- sample(c(1,0),m,TRUE)
        else if (MC[vnam,'dis'] == 'rnorm')
          vvec <- rnorm(m, mean=p[[vnam]][1],sd=p[[vnam]][2])   # ~N(u,sd)
        else if (MC[vnam,'dis'] == 'rlnorm')
          vvec <- exp(rnorm(m, mean=log(p[[vnam]][1]),sd=log(p[[vnam]][2])))   # exp(rnorm()) == rlnorm()  ~lognormal
        else if (MC[vnam,'dis'] == 'truncnorm') {
          #tmean <- p[[vnam]][1]
          #tsd   <- p[[vnam]][2]
          vvec <- rtruncnorm(m, a=MC[vnam,'min'], mean=p[[vnam]][1],sd=p[[vnam]][2])   # exp(rnorm()) == rlnorm()  ~lognormal
          #for (il in 1:5) {  
          #  vvec <- rtruncnorm(m, a=MC[vnam,'min'], mean=tmean,sd=tsd)  
          #  tmean <- tmean - (mean(vvec) - p[[vnam]][1])
          #  tsd   <- tsd *  p[[vnam]][2] / sd(vvec)
          #  cat('vvec[mu,sd]=[',mean(vvec),',',sd(vvec),']\n')
          #}
        } else  if (MC[vnam,'dis'] == 'runif')
         vvec <- runif(m, min(p[[vnam]]),max(p[[vnam]]))       # ~U(min,max)
        else
          stop('main:: parameter distribution not properly specified in MC')
      } else
        vvec <- rep(p[[vnam]][1], m)
    }
    gpar[i,] <- vvec
  }

  if (empiricalNorm && m > 1) {
    rnormid  <- which(MC$flag == TRUE & MC$dis == 'rnorm')
    if (length(rnormid) > 0) {
      if (!is.null(Sigma)) {
        if (sum(MC$flag) != ncol(Sigma))
          stop('mcEnsemble: sum(MC$flag) != ncol(Sigma)')
        Sid <- MC$dis[MC$flag] == 'rnorm'
        Sigma <- Sigma[Sid,Sid]
      } else {
        Sigma <- diag(sapply(p[rnormid],function(x){x[2]})^2)
      }
      gparNorm <- mvrnorm(n=m, 
                          mu=sapply(p[rnormid],function(x){x[1]}),
                          Sigma=Sigma,
                          empirical=TRUE) # transposed with respect DA standard
      gpar[rnormid,] <- t(gparNorm)
    } # end if rnormid
  }
  
  # add constraints
  for (i in 1:nrow(gpar)) {
    vnam <- rownames(gpar)[i]
    vvec <- gpar[i,]
    if (!is.na(MC[vnam,'min']))
      vvec <- pmax(vvec,MC[vnam,'min'])
    if (!is.na(MC[vnam,'max']))
      vvec <- pmin(vvec,MC[vnam,'max'])
    gpar[i,] <- vvec
  }
  
  if (!is.null(MCcons)) { # the 2nd term constraints the 1st 
   for (ic in 1:nrow(MCcons)) {
     if (MCcons[ic,2] == '<')
       gpar[MCcons[ic,1],] <- pmin(gpar[MCcons[ic,1],],gpar[MCcons[ic,3],])
     if (MCcons[ic,2] == '>')
       gpar[MCcons[ic,1],] <- pmax(gpar[MCcons[ic,1],],gpar[MCcons[ic,3],])
     if (MCcons[ic,2] == '=')
       gpar[MCcons[ic,1],] <- gpar[MCcons[ic,3],]
   }
  }
  
  return(gpar)
} # end function mcEnsemble()

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
} # end function orthoEnsemble()
