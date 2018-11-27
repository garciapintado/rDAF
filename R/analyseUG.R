analyseUG <- function(G=NULL, gLON=NULL, gLAT=NULL, fit, prm, X, yls=NULL, gauDA=NULL, gauTS=NULL, aug=NULL, ana=NULL,
                      dsn='.', debugmode=TRUE, retdy=FALSE, mpi=FALSE, analysis_scn='a0', theta=NULL, lite=TRUE) {

  # ensemble contructor as input to assimilation - analyseUG stands for unstructured grids
  #
  #
  # G            :: OPT, LIST,    gmeta6 class object      - structured grids. [$cells, $xseq (W-E longitudes),
  #                               $ryseq (N-S latitudes), $cols (grid dimension: longitude), $rows (grid dimension: latitude)
  # gLON,gLAT    :: OPT, LIST,    Geographical coordinates - unstructured grids [nx,ny] 'matrix' 
  # fit          ::      INTEGER, forecast time index [for output data]
  # prm          ::      LIST,    DA-specific data
  # X            ::      LIST,    state vector, structured as X[[xKIND]][[timelab]][[im]]
  # yls          :: OPT, LIST,    observation list, structured as yls[[yKIND]][[yTYPE]]$s lists. 
  # gauDA        :: OPT, LIST,    y <-> HE matches, to include observations by state-vector augmentation. Augments x and y
  # gauTS        :: OPT, LIST,    as gauDA lists, without the $y slot. Augments x, and can be used to obtain a smooth estimate (ETKS)
  #                               at specific locations/variables in the domain. 
  # aug          :: OPT, LIST,    augmentation blocks (parameters)
  # ana          ::      LIST,    analysis parameter for each variable type
  # DEM          :: OPT, byrow=TRUE vector representation. Digital Terrain Model - just required for water level variables
  # dsn          :: CHARACTER != '' 
  # retdy        :: OPT, TRUE for no assimilation, and instead return the innovation vector [just overpass observations]
  # mpi          :: whether to parallel DA in assimilate()
  # analysis_scn :: CHARACTER - appended to background and analysis files 
  # theta        :: OPT, LIST - required for parameter space KF [prm$method]. List components: 
  #                             $d  :: [ntheta] vector of parameter perturbations - has to comply with X, so the ensemble for X components
  #                                    match first the background, then one column for each $d element
  #                             $b  :: [ntheta] vector of background parameter values
  #                             $Pb :: [ntheta,ntheta] parameters error background covariance matrix
  # lite         :: LOGICAL, wheter to use the lite version: no MPI, no localization, no inflation, no rotation for ensemble filters   
  # revisions:
  # 2016-01-16 JGP
  #   re-engineering from former version
  #   now flexible assimilation building blocks for any number of variables to assimilate & analyse
  #   forward operator is only by 0-order nearest neighbourhood [no scaling, no representation, no interpolation]
  # 2017-03-15
  # For water level observations, the distributed observations have to be passed now as gauDA structures.
  # This makes this function more standard at the cost of increasing memory requirements.
  #
  # Notes: either G, or [gLON,gLAT] input assumes a common grid for all gridded variables, which must match the given dimensions. 
  #        State input via X[[vKIND]]$val[[im]], im=1,...m must then match the given metadata.
  #        That is a vector of length G$cells indicating coordinates as 1st by nx [longitude, in general],
  #        starting from the grid ul corner.   
  #        Thus, for staggered grids, variables on each grid should be analysed independently
  # 2017-10-18
  #   include pIKS, pMKS schemes. These are not an ensemble schemes, in the sense that covariances are explicit. 
  #   The method uses finite differences though (or conditional sampling), which are provided as a list of ensemble matrices,
  #   as input estimates of the sensitivites of the observations to the model inputs, to construct a Kalman update for the inputs
  # 2018-11-25
  #   Modifications to make it amenable to R-CRAN 
    
  options(warn=2)
  if (lite)
    mpi <- FALSE
  
  if (!is.null(G) && !is.null(gLON))
    stop('analyseUG: use just either G [regular grids] or gLON,gLAT [irregular grids] for horizontal grid definition')

  if (is.null(G))
   structG <- FALSE                                                             # irregular grid - although it also may be used for regular grids
  else
   structG <- TRUE                                                              # structured regular grid

  if (!file.exists(file.path(dsn,'results')))
    dir.create(file.path(dsn,'results'), recursive=TRUE)
  
  if (prm$method %in%  c('pIKS','pMKS')) {
    if (is.null(theta))
      stop('analyseUG:: ---ERR: theta list not provided for parameter-space method---')
    m <- length(theta$d)+1                                                      # extra member for the background in X
  } else 
    m <- prm$m                                                                  # ensemble size
  if (!is.null(G))
    ng <- G$cells                                                               # cells in each horizontal grid
  else
    ng <- length(gLON)

  x     <- NULL                                                                 # [n]   - state vector
  xpos  <- NULL                                                                 # [n,2] - geographical positions - [NA,NA] for global
  xz    <- NULL                                                                 # [n]   - vertical position      - [NA] for global 
  xKIND <- NULL                                                                 # [n]   - variable KIND
  xtime <- NULL                                                                 # [n] CHAR, state variable timestamp
  xgrd  <- NULL                                                                 # [n] LOGICAL, TRUE for gridded data | FALSE for sparse in space

  p     <- 0                                                                    # [1]   - number of observations
  y     <- NULL                                                                 # [p]   - observation vector
  ypos  <- NULL                                                                 # [p,2] - geographical positions - [NA,NA] global obs
  yz    <- NULL                                                                 # [p]   - vertical positions     - [NA]   global obs
  yKIND <- NULL                                                                 # [p]   - variable KIND
  yTYPE <- NULL                                                                 # [p]   - observation TYPE
  ytime <- NULL                                                                 # [p] CHAR, observation timestamp
  ygrd  <- NULL                                                                 # [p] LOGICAL, TRUE for variables mapped from gridded variables QC

  E     <- NULL                                                                 # [n,m]
  HE    <- NULL                                                                 # [p,m]
  R     <- NULL                                                                 # [p,p]
  H     <- NULL                                                                 # [p,n]

  getCooG <- function(G) {
    # get coordinates in a regular lattice as a 2-col matrix
    cbind(rep(G$xseq,G$rows),rep(G$ryseq,each=G$cols))}

  wSigma <- function(sigma, w) {
    # inflates a covariance keeping correlation
    diag(sqrt(w)) %*% sigma %*%  diag(sqrt(w))   
  }

  ## make y, ypos, HE
  # 1 :: gridded variables <-> satellite observations

  if (structG)
    gpos <- getCooG(G)                                                          # horizontal grid center positions
  else
    gpos <- cbind(as.numeric(gLON),as.numeric(gLAT))                            # horizontal grid center positions  

  # just for SpatialGraph-based localisation
  if (!is.null(prm$loc_boo)) {
    if (prm$loc_boo && !is.null(prm$loc_distype)) {
      if (prm$loc_distype == 'SG') {                                            # along-network distance calculations
        prm$sglst    <- list()                                                  # information for SpatialGraph-based distance
        prm$sglst$sg <- readRDS(file.path(dsn,prm$sgf))                         # 'SpatialGraph'
        gposSG       <- readRDS(file.path(dsn,prm$gridSGf))                     # precalculated along-network distances for grid nodes
      }
    }
  }

  if (is.list(X[[1]]$val[[1]])) {
   ismatX <- FALSE
  } else if (is.matrix(X[[1]]$val[[1]])) {
   ismatX <- TRUE
  } else {
   stop('data in X could not be identified either as list or matrix')
  }

  getXtimes <- function(x){names(x$val)[which(sapply(x$val,length) != 0)] }
  vts  <- lapply(X, FUN=getXtimes)                                       # list : CHARACTER timelab / gridded vKIND
  nvts <- sapply(vts,length)                                             # vector: number of timesteps per gridded x vKIND
  if (!ismatX)
    vlen <- sapply(X, function(x){length(x$val[[vts[[1]]]][[1]])})         # vector: state-vector length for this variable
  else
    vlen <- sapply(X,function(x){nrow(x$val[[1]])})
  nzs  <- vlen/ng                                                        # vector: number of vertical levels / variable
  xzl  <- lapply(X,FUN=function(x){x$z})                                 # list: vertical levels / variable
  if( !identical(as.numeric(nzs), as.numeric(sapply(xzl,'length'))) ) 
    stop('analyseUG: xzl vector does not match given grid layers')

  # stack X as E
  E     <- matrix(NA, nrow=sum(nvts*vlen), ncol=m)                             # just gridded data
  xpos  <- matrix(NA, nrow=sum(nvts*vlen), ncol=2)                             # init for gridded state vector time-stacks
  xpos[,1] <- gpos[,1]                                                         # recycle [~F95 spread(gpos,1,sum(nvts))]
  xpos[,2] <- gpos[,2] 
  xKIND <- rep(names(X),      times=nvts*vlen)
  xtime <- rep(do.call(c,vts),times=rep(vlen,nvts))                            # CHARACTER
  xgrd  <- rep(TRUE,nrow(xpos))                                                # spatially distributed [i.e. either structured or unstructured grids]
  xzlg  <- lapply(xzl, function(x,ng){rep(x,each=ng)},ng=ng)

  for (iv in 1:length(X)) {
    xz    <- c(xz,rep(xzlg[[iv]],nvts[iv]))
  }; rm(xzlg,xzl)

  FUN <- function(x,vlen) {if (length(x)==0) 
                             x <- rep(NA,vlen)
                           else 
                             x <- as.vector(x)
                           return(x)}
  for (iv in 1:length(X)) {                                                   # iterate through variables
    vKIND <- names(X)[iv]
    for (vit in 1:length(X[[iv]]$val)) {                                      # through times
      vtime <- names(X[[iv]]$val)[vit]                                        # CHAR
      xboo  <- as.logical(xKIND == vKIND & xtime == vtime)
      if (sum(xboo) == 0)
        next
      if (is.list(X[[iv]]$val[[vit]]) && length(X[[iv]]$val[[vit]]) == m) {
        E[xboo,] <- vapply(X[[iv]]$val[[vit]],FUN=FUN, FUN.VALUE=rep(0,sum(xboo)), vlen=vlen[iv])
      } else if (is.matrix(X[[iv]]$val[[vit]]) && ncol(X[[iv]]$val[[vit]]) == m) {
        E[xboo,] <- X[[iv]]$val[[vit]]
      } else
        stop('X[[iv]]$val[[vit]] not identified as an ensemble matrix or list')
      if (!is.null(prm$sglst)) {
        if (is.null(prm$sglst$xposSG))
          prm$sglst$xposSG <- gposSG                                            # init
        else
          prm$sglst$xposSG <- rbind(prm$sglst$xposSG,gposSG)                    # augment
      }
    }
  }

 # y, ypos, yKIND, yTYPE, ytime, HE
 if (!is.null(yls)) {
  for (iv in 1:length(yls)) {                                                   # iterate through variable KIND
    vKIND <- names(yls)[iv]
    if (all(sapply(yls[[vKIND]], function(x){is.null(x$s)})))                  # restricted to yFMT=='s'
      next

    cat('mapping X into',vKIND,'observations\n')
    if (!(vKIND %in% names(X)))
      stop('yls$',vKIND,' not a valid vKIND in X')
         #      cat('analyse:: stacking satellite observations|',vKIND,'|vit:',vit,'\n')
    for (j in 1:length(yls[[vKIND]])) {                                        # iterate through yTYPE observations
      if (is.null(yls[[vKIND]][[j]]$s$data) || !yls[[vKIND]][[j]]$s$use)                                  # matching observation type exists
        next
      yTYPEo <- names(yls[[vKIND]])[j] 
      vytimes <- yls[[vKIND]][[yTYPEo]]$s$data$timelab                                 # CHARACTER     
      yfield  <- yls[[vKIND]][[yTYPEo]]$s$data$yfield
      if (is.null(yfield))
        yfield <- 'y'

        #if (!inherits(vytimes,'POSIXct'))
        #  stop('analyse:: yls:',vKIND,' times not POSIXct')
      cat('vytimes:',vytimes,'\n')
      y_usew <- yls[[vKIND]][[yTYPEo]]$s$w
      if (is.null(y_usew))
        ylsusew <- FALSE
      
      for (yit in 1:length(vytimes)) {
        cat('OP:',vytimes[yit],'\n')
        if (!(vytimes[yit] %in% xtime))
          next
        y_op    <- yls[[vKIND]][[yTYPEo]]$s$data[[vytimes[yit]]][[yfield]]                        # [p_op]
        ypos_op <- as.matrix(yls[[vKIND]][[yTYPEo]]$s$data[[vytimes[yit]]][['pos']])              # [p_op x 2] coordinate matrix
        yz_op   <- yls[[vKIND]][[yTYPEo]]$s$data[[vytimes[yit]]][['z']]                           # [p_op]
        qced_op <- yls[[vKIND]][[yTYPEo]]$s$data[[vytimes[yit]]][['qced']]                        # [p_op] flag to assimilate just quality-controlled observations
        if (is.null(qced_op))
          qced_op <- rep(TRUE,length(y_op))
        yposOK <- ypos_op[,1] >= min(gpos[,1]) & ypos_op[,1] <= max(gpos[,1]) &                  # indices of observations within the domain for the analysis
                  ypos_op[,2] >= min(gpos[,2]) & ypos_op[,2] <= max(gpos[,2])
        qced_op <-  qced_op & yposOK
        if (y_usew) {                                                                           # apply observation weights -> inverse multiplier of R
          yw_op <- yls[[vKIND]][[yTYPEo]]$s$data[[vytimes[yit]]][['w']]
          if (is.null(yw_op))
            stop('analyse:: --observation weights requested but not available--')
        } else {
          yw_op <- rep(1.0, length(y_op))
        }

        y_op    <- y_op[qced_op]
        ypos_op <- ypos_op[qced_op,]
        yz_op   <- yz_op[qced_op]
        yw_op   <- yw_op[qced_op]
        p_op    <- length(y_op)

        if (!is.null(yls[[vKIND]][[yTYPEo]]$s$data[[vytimes[yit]]][['R']])) {                                   # 1st try: R matrix input a list with as many R-matrices as overpasses
          R_op <- yls[[vKIND]][[yTYPEo]]$s$data[[vytimes[yit]]][['R']][qced_op,qced_op]
          if (y_usew) {
            if (all(R_op[!diag(nrow(R_op))] == 0)) {    # observation weight scaling
              R_op <- diag(diag(R_op)/yw_op)            # if it is diagonal
            } else { # full covariance scaling required
              R_op <- wSigma(R_op,1/yw_op)
            }
          }; rm(y_usew)
        } else if (!is.null(yls[[vKIND]][[yTYPEo]]$s$data[[vytimes[yit]]][['r']])) {                            # 2nd try: vector variance, specific to observations - uncorrelated errors
          R_op <- Diagonal(x=yls[[vKIND]][[yTYPEo]]$s$data[[vytimes[yit]]][['r']][qced_op]/yw_op)
        } else if (is.numeric(yls[[vKIND]][[yTYPEo]]$s$r)) {                                                    # 3rd try: common scalar variance for all this [[yTYPEo]][[yFMT]]
          if (length(yls[[vKIND]]$s$r) != 1)
            stop('analyse:: length of scalar observation error variance != 1')
          R_op <- Diagonal(x=rep(yls[[vKIND]][[yTYPEo]]$s$r, p_op)/yw_op)
        } else {
          stop('analyse:: observation error variance not found for vKIND:',vKIND)
        }
        if (!all(dim(R_op) == p_op))
          stop('analyse:: yls$',vKIND,'[[',yTYPEo,']] R dimensions do not match p_op at ',vytimes[yit])
        if (nrow(ypos_op) != p_op || length(yz_op) != p_op)
          stop('analyse:: yls$',vKIND,'[[',yTYPEo,']] position dimensions do not match p_op at ',vytimes[yit])

        H_op  <- calcHnn(G, gLON, gLAT, ypos_op)                              # [p_op,ng] dgCMatrix
        xboo  <- as.logical(xKIND == vKIND & xtime == vytimes[yit])
        if (sum(xboo) != ng)
          stop('analyseUG: full grid not identified for ',vKIND,'.',yTYPEo,'|time:',vytimes[yit])
        #HE_op <- H_op %*% E[xboo,]                                            # [p_op x ng].[ng x m] dgeMatrix E-NA columns maps into HE_op NA columns
        #HE_op <- as.matrix(HE_op)                                             # [p_op, m] matrix
        Ho <- Matrix(0,nrow(H_op),nrow(E))                                    # [p_op,n]
        Ho[,xboo] <- H_op

        p     <- p + p_op
        y     <- c(y, y_op)
        ypos  <- rbind(ypos, ypos_op)
        yz    <- c(yz,   yz_op)
        yKIND <- c(yKIND, rep(vKIND,p_op))
        yTYPE <- c(yTYPE, rep(yTYPEo,p_op))
        ytime <- c(ytime, rep(vytimes[yit],p_op))
        ygrd  <- c(ygrd,  rep(TRUE,  p_op))
        #HE   <- rbind(HE, HE_op)       not needed -- recalculated later in the transformed space
        if (is.null(R)) {                                                  # 1st observation dataset
          R <- R_op
          H <- Ho                                                          # dgCMatrix
        } else {                                                           # with possible intermediate x augmentation
          R <- bdiag(R,R_op)
          H <- rbind(H,Ho)
        }
        rm(p_op, y_op, ypos_op, yz_op, qced_op, yposOK, yw_op, R_op, H_op, Ho, xboo)
      } # end for overpass
    } # end for yTYPE
  } # end for vKIND
 } # end if not NULL yls

  # 2 :: gauDA [y <-> HE] state-vector and observation vector augmentation
  cat('analyse:: time series observations\n')
  p_gau <- 0
  if (length(gauDA) > 0) {
    E_gau     <- NULL
    y_gau     <- NULL
    ypos_gau  <- NULL
    yz_gau    <- NULL
    yKIND_gau <- NULL                                                      # generic observation KIND
    yTYPE_gau <- NULL                                                      # specific observation TYPE
    ytime_gau <- NULL
    Rblk      <- NULL
    for (iv in 1:length(gauDA)) {                                          # variable types
      vKIND   <- names(gauDA)[iv]      
      p_this  <- length(gauDA[[iv]]$y)                                     # available observations
      if (p_this == 0)
        next
      E_gau     <- rbind(E_gau,gauDA[[iv]]$HE)
      y_gau     <- c(y_gau, gauDA[[iv]]$y)
      ypos_gau  <- rbind(ypos_gau, gauDA[[iv]]$pos )
      yz_gau    <- c(yz_gau, gauDA[[iv]]$z)
      yKIND_gau <- c(yKIND_gau, rep(vKIND, p_this))
      yTYPE_gau <- c(yTYPE_gau, gauDA[[iv]]$yTYPE)
      ytime_gau <- c(ytime_gau, gauDA[[iv]]$timelab)
      
      if (!is.null(gauDA[[iv]]$R))                   # full R matrix input [e.g. preprocessed online for temporal correlation model]
        Rgau <- gauDA[[iv]]$R
      else
        Rgau <- Diagonal(x=gauDA[[iv]]$r)            # ddiMatrix
      if (is.null(Rblk))
        Rblk <- Rgau
      else
        Rblk <- bdiag(Rblk,Rgau)
    } # end for iv TS assim
    p_gau <- length(y_gau)
    ygrd <- c(ygrd,rep(FALSE,p_gau))
  } # end if(length(gauDA))

  if (p_gau > 0) {
    E     <- rbind(E, E_gau)
    xpos  <- rbind(xpos, ypos_gau)
    xz    <- c(xz, yz_gau)
    if (!is.null(prm$sglst))
      prm$sglst$xposSG <- rbind(prm$sglst$xposSG,
                                pointsToLines(ypos_gau, prm$sglst$sg@e))
    xKIND <- c(xKIND, yKIND_gau)
    if (is.null(xtime)) {
      xtime <- ytime_gau
    } else {
      xtime <- c(xtime, ytime_gau)
    }
    xgrd <- c(xgrd, rep(FALSE,p_gau))
    y     <- c(y, y_gau)
    ypos  <- rbind(ypos, ypos_gau)
    yz    <- c(yz, yz_gau)
    yKIND <- c(yKIND, yKIND_gau)
    yTYPE <- c(yTYPE, yTYPE_gau)
    ytime <- c(ytime, ytime_gau)


    #HE  <- rbind(HE, E_gau)                                           # as  HE_gau == E_gau -- not needed as recalculated later in transformed space
    if (is.null(R)) {
      H <- Matrix(0,p_gau,nrow(E))
      H[,(nrow(E)-p_gau+1):nrow(E)] <- Diagonal(p_gau)
      R <- Rblk
    } else {
      H <- bdiag(H,Diagonal(p_gau))
      R <- bdiag(R,Rblk)
    }
    rm(p_gau, E_gau, y_gau, ypos_gau, yKIND_gau, yTYPE_gau, ytime_gau)
  } # end if (p_gau > 0)

  if (prm$rfactor != 1) {
    cat('analyseUG:: inflating R by rfactor:',prm$rfactor,'\n')
    R <- prm$rfactor * R                                                           # e.g. for inflating R for fractional assimilation
  }
  # observations done!

  # 3 :: augment state-vector with gauTS to get smoothed time series
  cat('analyse:: gauTS augmentation',as.character(Sys.time()),'\n')
  if (!is.null(gauTS))
    stop('analyse :: code not ready as TS smoother - pending augmentation by gauTS')

  # 4 :: augment state-vector with parameters [boundary conditions + model parameters]

  if (!is.null(aug)) {
    for (iv in 1:length(aug)) {
      vKIND <- names(aug)[iv]
      E <- rbind(E,aug[[iv]]$E)
      xpos  <- rbind(xpos, aug[[iv]]$pos)
      xz    <- c(xz, aug[[iv]]$z)
      if (!is.null(prm$sglst))
        prm$sglst$xposSG <- rbind(prm$sglst$xposSG,
                                  pointsToLines(matrix(aug[[iv]]$pos,ncol=2), prm$sglst$sg@e)) # note this is meaningless for global parameters
      xKIND <- c(xKIND, rep(vKIND,length(aug[[iv]]$timelab)))
      xtime <- c(xtime, aug[[iv]]$timelab)
      xgrd  <- c(xgrd,  rep(FALSE,  length(aug[[iv]]$timelab)))
    }
  }

  x <- rowMeans(E, na.rm=TRUE)                                             # [n,1] ensemble mean (augmented state-vector) - physical space
  if (debugmode) {
    cat('analyse:: writing debugmode=TRUE files',as.character(Sys.time()),'\n')
    #xvar <- varSparse(E)
    xvar <- apply(E,1,var,na.rm=TRUE)
    saveRDS(E,       file=file.path(dsn,"results",paste("Eaug_",gsub('a','b',analysis_scn),"_",fit,".rds",sep="")))
    saveRDS(xvar,    file=file.path(dsn,"results",paste("xvar_",gsub('a','b',analysis_scn),"_",fit,".rds",sep="")))

  }

  if (length(y) == 0 && !retdy)
    return(list(E=E, xdf=data.frame(xKIND=xKIND, xgrd=xgrd, xtime=xtime, xpos=xpos)))

  if (nrow(E) > ncol(H))                                                   # if 'gauTS' or 'aug' augmentation exists - append columns with 0s
    H <- cbind2(H,Matrix(0,nrow(H),nrow(E)-ncol(H)))                       # [p,n] dgCMatrix :: H does not include transformations here

  if (prm$method %in% c('pIKS','pMKS')) {
    if (is.null(theta))
      stop('theta list not provided for pKs parameter estimation')
    theta$R <- R

    if (is.null(theta$G)) {                                                # sensitivity matrix: if not externally estimated
      if (ncol(E) != length(theta$b) + 1)                                  #                     get here by simple forward FD
        stop('analyseUG:: dim(E) not compliant for finite-difference sensitivity estimation for pKs')
      Gx  <- (E[,-1] -  E[,1]) / rep(theta$d, each=nrow(E))                # [n,ntheta] finite difference approximation to the Jacobian of the observations with respect to the model parameters
      theta$G  <- H %*% Gx                                                 # [p,ntheta] = [p,n][n,nbeta]
    }
    if (prm$method == 'pMKS') {
      dyb <- y - H %*% E[,1]                                               # for a fully nonlinear H, this should be calculated out of analyseUG
      theta$PHT  <- theta$Pthetab %*% t(theta$G)
    } else {         # pIKS
      dyb <- y - H %*% E[,1] - theta$G %*% (theta$b0 - theta$b)            #  "            " | theta$b0 is theta_b, theta$b is theta_l
      theta$PHT  <- theta$Pthetab0 %*% t(theta$G)
    }
    dyb <- as.matrix(dyb)                                                  # 'dgeMatrix' from the above operation does not allow for data.frame storage
    theta$HPHT <- theta$G %*% theta$PHT
    theta$kappa <- kappa(theta$HPHT+R)                                     # condition number
    if (theta$kappa > 1.0E06)
      stop('analyseUG:: theta$kappa > 1.0E02 [',theta$kappa,'] -- consider prescaling to improve covariance matrix condition')
    theta$K     <- theta$PHT %*% solve(theta$HPHT + R)                     # [ntheta,p]
    I_KG       <- diag(rep(1,length(theta$b))) - theta$K %*% theta$G
    if (prm$method == 'pMKS') {    
      theta$a     <- theta$b + theta$K %*% dyb
      theta$Pthetaa <- I_KG %*% theta$Pthetab %*% t(I_KG) + theta$K %*% R %*% t(theta$K)                          # Jazwinski (1970), p. 270 after Aoki0 (1967)
    } else {         # pIKS
      theta$a     <- theta$b0 + theta$K %*% dyb
      theta$Pthetaa <- I_KG %*% theta$Pthetab0 %*% t(I_KG) + theta$K %*% R %*% t(theta$K)                         # Jazwinski (1970), p. 270 after Aoki0 (1967)      
    }

    if (debugmode) {
      dydf <- data.frame(y=y, dyb=dyb, yKIND=yKIND, yTYPE=yTYPE, ygrd=ygrd, ytime=ytime, ypos.lon=ypos[,1], ypos.lat=ypos[,2],
                         yz=yz, stringsAsFactors=FALSE)
      xdf <- data.frame(xa=NA, xb=E[,1], dx=NA, xKIND=xKIND, xgrd=xgrd, xtime=xtime, xpos.lon=xpos[,1],xpos.lat=xpos[,2],
                        xz=xz, stringsAsFactors=FALSE)                                                            # physical space analysis, increments & metadata
      saveRDS(dydf,file=file.path(dsn,"results", paste("dydf_",analysis_scn,"_",fit,".rds",sep="")))
      saveRDS(xdf, file=file.path(dsn,"results", paste("xdf_", analysis_scn,"_",fit,".rds",sep="")))
    }

    #theta$Pthetaa <- I_KG %*% theta$Pthetab # standard
    return(theta)
  }

  #sdB <- sdSparse(E)                                                        # pre-transform for eventual inflation
  sdB <- apply(E,1,sd,na.rm=TRUE)

  # 6 :: adapt data to assimilate()

  # get x-based localisation lengths
  cat('analyse:: getting localisation lengths',as.character(Sys.time()),'\n')
  xKINDs <- unique(xKIND)
  if (!all(xKINDs %in% names(ana)))
   stop('analyse:: analysis parameters not provided for all xKINDs')
  lls <- sapply(ana, function(x){x$ll})
  xllen  <- rep(0,nrow(E))
  for (iv in 1:length(xKINDs)) {
    xllen[xKIND == xKINDs[iv]] <- lls[xKINDs[iv]]
  }; rm(lls)

  #browser()

  # forward transform :: individualised functions for specific kind of variables in the ensemble
  # note : observations are not transformed - thus for observed vKIND variables, if observations need to be transformed,
  #        it has to be done prior to call analyse(), and the forward transform function of HE here is into
  #        the forward transform space of y. 
  # Also, ana[[vKIND]]$trf should be set to 'none' for this function for vKIND variables in gauDA, and transformations done externally if needed
  cat('analyse:: getting forward transforms',as.character(Sys.time()),'\n')
  for (iv in 1:length(xKINDs)) {
    vKIND <- xKINDs[iv]
    if (ana[[vKIND]]$trf == 'none')
      next
    if (ana[[vKIND]]$trf %in% c('autoGA','GA'))  # has to be done externally
      next
    else {
      xids <- which(xKIND == vKIND)
      #yids <- which(yKIND == vKIND)
    }
    if (ana[[vKIND]]$trf == 'qpow35') {
      E[xids,] <- qpow35(E[xids,])
      #y[yids]  <- qpow35(y[yids])     not ready because R would then have to be transformed as well
    } else if (ana[[vKIND]]$trf == 'linearScale') {
      E[xids,] <- linearScale(E[xids,], b=ana[[vKIND]]$linmul) # warning forced similar for all yTYPES 
    } else if (ana[[vKIND]]$trf == 'log') { # natural logarithm  
      E[xids,] <- log(E[xids,])
    } else
      stop('analyse :: transformation function ',ana[[vKIND]]$trf,' for ',vKIND,' not available')
  }
  cat('analyse:: adapting data to assimilate',as.character(Sys.time()),'\n')

  # HE <- calcHE_Hm(H,E,Hm=NULL)                                        # not necessarily linear - get HE - transformed [observation] space
  HE <-  as.matrix(H %*% E)

  if (length(y) > 0) {
    yb <- rowMeans(HE, na.rm=TRUE)                                    # Hunt_al2007 notation \overline{y^b}
    dyb <- y - yb                                                     # [p]   average prior innovation vector - observation space
    dydf <- data.frame(y=y, dyb=dyb, yKIND=yKIND, yTYPE=yTYPE, ygrd=ygrd, ytime=ytime, ypos.lon=ypos[,1], ypos.lat=ypos[,2],
                       yz=yz, stringsAsFactors=FALSE) 
    if (retdy)                                                        #       to QC
      return(dydf)

    if (!is.null(prm$sglst))
     prm$sglst$yposSG <- pointsToLines(ypos, prm$sglst$sg@e)
  } else {
   if (retdy)
     return(NULL)
  }

  HA <- Matrix(HE - rowMeans(HE, na.rm=TRUE))      # [p,m] perturbations of the transformed HE matrix (recycled Hx)
  rm(HE)                                           # warning: transformation of observed variables not ready
                                                   #          this would involve trf(R),trf(y)
                                                   #          and re-calculate HA, dyb in the transformed space to assimilate()
  xf <- rowMeans(E, na.rm=TRUE)                    # transformed space       - state mean
  A <- E - xf                                      # [n,m] transformed space - ensemble perturbations
  A <- Matrix(A)                                   #  'dgCMatrix', else 'dgeMatrix' 
  #R <- Matrix(R)                                  # R is 'dtCMatrix'[triangular sparse matrix],  or 'ddiMatrix' [diagonal]
  sm <-  !is.na(E[which(!is.na(sdB))[1],])                              # LOGICAL
  prm$n <- nrow(E)
  rm(E) # release memory for assimilation - recreate later

  Kfname <- file.path(dsn,"results", paste("K_",analysis_scn,"_",fit,".rds",sep="")) # filename to save Kalman gain matrix

  cat('analyse:: start assimilation |',as.character(Sys.time()),'| p =',length(dyb),'\n')
  if (lite) 
    dxA <- assimilateLite(prm, A[,sm], HA[,sm], dyb, R, debugmode, Kfname)
  else {
    if (!exists('assimilate')) {
      cat('analyseUG: WARNING:: switch to assimilateLite(), as assimilate() not available in this version')
      assimilate <- assimilateLite
    }
    dxA <- assimilate(prm, A[,sm], HA[,sm], dyb, R, debugmode, Kfname,
                      xpos, xllen, ypos, mpi=mpi, sglst=prm$sglst)
  }
  cat('analyse:: end assimilation   |',as.character(Sys.time()),'\n')

                                                                      # transformed space - analysis mean
  xa <- xf + as.numeric(dxA$dx)                                       # 'dgeMatrix' to standard vector
  
  if (prm$method != "EnOI") {
    E <- matrix(NA,length(xa),m)
    E[,sm] <- as.matrix(dxA$A + xa)                                   # transformed space - analysis ensemble  [cast from dgeMatrix]
  }

  # Hxa <- rowMeans(calcHE_Hm(H,E,Hm), na.rm=TRUE)
  Hxa <-  rowMeans(H %*% E, na.rm=TRUE)
  dydf$dya <- y - Hxa                                                 # [p]   average posterior innovation vector - observation space

  # inverse transform into model space
  xKINDs <- unique(xKIND)
  for (iv in 1:length(xKINDs)) {
    vKIND <- xKINDs[iv]
    if (ana[[vKIND]]$trf == 'none')
      next
    if (ana[[vKIND]]$trf %in% c('autoGA','GA'))  # has to be done externally
      next
    else
      xids <- which(xKIND == vKIND)
    if (ana[[vKIND]]$trf == 'qpow35') {
      E[xids,] <- qpow35(E[xids,], inverse=TRUE)
    } else if (ana[[vKIND]]$trf == 'linearScale') {
      E[xids,] <- linearScale(E[xids,], b=ana[[vKIND]]$linmul, inverse=TRUE)
    } else if (ana[[vKIND]]$trf == 'log') {
      E[xids,] <- exp(E[xids,])
    } else
      stop('analyse:: transformation function not available')
  }

  xa <- rowMeans(E, na.rm=TRUE)                                       # physical space - analysis mean
  A  <- E - xa                                                        # physical space - analysis perturbations

  # inflate in physical ensemble space
  # note: variables with forward/inverse transform external
  #       to analyseUG are however inflated in the transformed space!
  # inflation and rotation
  if (prm$method != 'EnOI') {

    for (iv in 1:length(xKINDs)) {
      vKIND <- xKINDs[iv]
      infac <- ana[[vKIND]]$infac
      #cat('analyse::vKIND:',vKIND,'|infac:',infac,'\n')
      if (infac == 1.0)
        next
      else
        xids <- which(xKIND == vKIND)
      if (infac < 1.0) {                                   # pseudo-automatic
        if (length(xids) > 1) {
          # sdA <- sdSparse(A[xids,])
          sdA <- apply(A[xids,],1,sd,na.rm=TRUE)
        } else
          sdA <- sd(A[xids,])
        infav <- (1.0 - infac) * sdB[xids] / sdA + infac   
        infav <- pmin(infav,10.0)                           # hard threshold. Consider improving
        infav[is.na(infav)] <- 0.0
      } else {                                             # fix infac
        infav  <- rep(infac, length(xids))
      }
      A[xids,] <- A[xids,] * infav                         # recycle infav for each col
    } # end for inflation

    # TODO: randomly rotate the ensemble every prm$rotate steps
    #if (prm$rotate && prm$rotate_ampl != 0) {
    #  if (prm$rotate_ampl == 1)
    #    A[xids,] <- A[xids,] * genU(m)
    #  else
    #    A[xids,] <- A[xids,] * genUeps(m, prm$rotate_ampl)
    #}
    E <- A + xa                                                       # physical space - analysis ensemble
  } # end inflation/rotation

  xdf <- data.frame(xa=xa, xb=x, dx=xa-x, xKIND=xKIND, xgrd=xgrd, xtime=xtime, xpos.lon=xpos[,1],xpos.lat=xpos[,2],
                    xz=xz, stringsAsFactors=FALSE)  # physical space analysis, increments & metadata
   
  if (debugmode) {
    xvar = apply(E,1,'var', na.rm=TRUE)
    saveRDS(E,   file=file.path(dsn,"results", paste("Eaug_",analysis_scn,"_",fit,".rds",sep="")))
    saveRDS(xvar,file=file.path(dsn,"results", paste("xvar_",analysis_scn,"_",fit,".rds",sep="")))
    saveRDS(dydf,file=file.path(dsn,"results", paste("dydf_",analysis_scn,"_",fit,".rds",sep="")))
    saveRDS(xdf, file=file.path(dsn,"results", paste("xdf_", analysis_scn,"_",fit,".rds",sep="")))
  }
  
  return(list(E=E, xdf=xdf, dydf=dydf))

} # end function analyseUG
