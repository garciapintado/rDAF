# function [dx, A] = assimilate(prm, A, xpos, xllen, HA, dy, ypos, R, debugmode, Kfname, mpi=TRUE, sglst=NULL)
#
# Calculate correction of the ensemble mean and update ensemble perturbations
#
#
# TYPE(T_prm)             :: prm       - system parameters
# REAL, DIMENSION(:,:)    :: A         - [n,m] ensemble perturbations
# REAL, DIMENSION(:,:)    :: HA        - [p,m] ensemble observations
# REAL, DIMENSION(:)      :: dy        - [p]   vector of innovations, dy = y - Hx
# REAL, DIMENSION(:,:)    :: R         - [p,p] observation error covariance matrix
# LOGICAL                 :: debugmode -
# CHARACTER(len=*)        :: Kfname    -  complete path/filename to save the K matrix
# LOGICAL, OPTIONAL       :: mpi
# TYPE(T_sglst)           :: sglst     - SpatialGraph data
#
# return list with two components:
#        $dx     - correction of the mean dx = K dy (numeric vector)
#        $A      - updated ensemble perturbations       [n x m] matrix
#
# File:           assimilate.R
#
# Created:        2012-03-23 JGP
#
# Last modified:  2016-01-05 JGP
#
# Author:         Javier Garcia-Pintado
#                 MARUM, University of Bremen
#
# Purpose:        Core code for EnKF-type DA - Lite version [no localization, no MPI]
#
# Description:    This procedure calculates the analysis correction and updates
#                 the ensemble by applying scheme specified by prm$method.
#                 It assumes that the system runs in the asynchronous regime,
#                 i.e. that both increments `dy' and ensemble observations `HA'
#                 contain values recorded at the time of observations.
#                     Not all schemes are available for every localisation
#                     method. Following are the lists of available schemes for
#                     each method.
#                     * With no localisation:
#                       EnKF
#                       ETKF    ::     This is the "symmetric" solution of
#                                      the ETKF (an implementation of the
#                                      Ensemble Square Root Filter, ESRF).
#                                      See eq(13) in Sakov_Oke2008a
#                  The function returns $dx$, and $A^a$, according with
#                  eqs. (3) and (4) in Sakov et al., 2010.
##
## This file is part of the free-software rDAFlite
## See LICENSE for details.

assimilateLite <- function(prm, A, HA, dy, R, debugmode, Kfname,
                           xpos = NULL, xllen = NULL, ypos = NULL, mpi = NULL, sglst = NULL) {
  # prm       : list with assimilation parameters
  # A         : [n,m] 'dgCMatrix' or 'dgeMatrix', n: state vector length, m: ensemble members
  # HA        : [p,m] REAL, where p is the number of observations
  # dy        : [p]   innovation vector
  # ypos      : [p,2] matrix, observation coordinates
  # R         : [p,p] 'ddiMatrix', 'dsyMatrix', or 'dtCMatrix': observation error covariance
  # debugmode : LOGICAL, to trigger output writing
  # Kfname    : CHARACTER, path/filename to write the Kalman gain, if debugmode=TRUE

  # arguments to provide compatibility with assimilate()  
  if (!is.null(xpos))
    cat('assimilateLite:: xpos forced to NULL\n')
  if (!is.null(xllen))
    cat('assimilateLite:: xllen forced to NULL\n')
  if (!is.null(ypos))
    cat('assimilateLite:: ypos forced to NULL\n')
  if (!is.null(mpi))
    cat('assimilateLite:: mpi forced to NULL\n')
  if (!is.null(sglst))
    cat('assimilateLite:: sglst forced to NULL\n')
  xpos  <- NULL
  xllen <- NULL
  ypos  <- NULL
  mpi   <- NULL
  sglst <- NULL

  stopifnot(prm$method %in% c('EnKF','ETKF'))                         # just EnKF, ETKF defined by now

  n <- nrow(A)                                                          # includes any possible augmentation
  m <- ncol(A)
  p <- length(dy)                                                     # number of observations == nrow(HA)
  
  dx <- Matrix(0,n,1, sparse=FALSE)                                  # init dx
  if (p < 1) {
    A <- A
    return(list(dx=dx,A=A))                                           # no state correction, and return background anomalies
  }

  if (any(dim(R) != p))
    stop('assimilate :: dim(R) != [p,p]')

  #   Sakov, P., G. Evensen, and L. Bertino (2010): Asynchronous data
  #   assimilation with the EnKF. Tellus 62A, 24-29.

  #Rsq    <- Matrix(expm::sqrtm(R))                                    # [p,p] R^{1/2} 'ddiMatrix' or 'dsyMatrix' 
  Ri     <- Matrix::solve(R)                                          # [p,p] R^{-1}  'dtCMatrix' or 'dsyMatrix' method. Call LAPACK 
  if (isDiagonal(Ri))
    Risq   <- sqrtm_rs(Ri)
  else
    Risq   <- Matrix(expm::sqrtm(Ri))                                   # [p,p] R^{-1/2}'dtCMatrix' or 'ddiMatrix'

  Rdiag  <- diag(R)                                                   # [p] just for output. Warning! neglects non-diagonal R values
  s <- Risq %*% dy / sqrt(m-1)                                        # [p,1] 'dgeMatrix', eq. (5) being R a full measurement error covariance matrix 
  S <- Risq %*% HA / sqrt(m-1)                                        # [p,m] 'dgeMatrix', eq. (6) Being R as above. S = R^{-1/2}HX = R^{-1/2}Y
                                                                      #       also eq. (13) in Sakov_Bertino (2010)  
  # EnKF global perturbation
  if (prm$method == 'EnKF') {
    #Dtilde <- Rsq %*% matrix(rnorm(p*m),p,m)                                     # [p,m] sample from R
    #D      <- Risq %*% Dtilde / sqrt(m-1) == Risq %*% Rsq %*% Matrix(...) # [p,m] = [p,p][p,m]
    D <- Matrix(rnorm(p*m),p,m) / sqrt(m - 1)
    D <- (D - rowMeans(D))
  }


  # no localisation
    if (prm$method == 'ETKF') {                                       # [crossprod(S) == t(S) %*% S]
      M <- solve(Diagonal(m) + crossprod(S))                          # Eq.(7a) [m,m] = [m,m] + [m,p][p,m] | 'dsyMatrix' <- class(Diagonal(m))='ddiMatrix' | class(M)='dsyMatrix' [Symmetric real matrix in non-packed storage]
      G <- tcrossprod(M,S)                                            #          [m,p] = [m,m][m,p]         | 'dgeMatrix'
    } else if (m <= p) { # not ETKF (do not allocate M)
      G <- tcrossprod(solve(Diagonal(m) + crossprod(S)),S)            # Eq.(7a) 
    } else {
      G <- crossprod(S, solve(Diagonal(p) + tcrossprod(S)))           # Eq.(7b) [m,p] = [m,p]([p,p] + [p,m][m,p])
    }
    dx <- A %*% G %*% s                                               # Eq.(3)  [n,1] = [n,m][m,p][p,1]    | 'dgeMatrix' [n,1]

    if (debugmode)
      K <- (A / sqrt(m-1)) %*% G %*% Risq                             # 'dgeMatrix' [n,p] = [n,m][m,p][p,p] Just for post-analysis | K = X G R^{-1/2}
                                                                      # e.g. [eq. (14) in Sakov & Bertino, 2010)
    if (prm$method == 'EnKF') {
      Tr <- G %*% (D - S)                                             # Eq.(8) [m,m] = [m,p][p,m]
      A  <- A %*% (Diagonal(m) + Tr)                                  # Eq.(4) [n,m] = [n,m][m,m]
    } else if (prm$method == 'DEnKF') {
      Tr <- -0.5 * G %*% S                                            # Eq.(10) [m,n]
      A <- A %*% (Diagonal(m) + Tr)                                   # Eq.(4) [n,m] = [n,m][mm]
    } else if (prm$method == 'ETKF') {

     # if (isDiagonal(M))
        A <- A %*% sqrtm_rs(M)
      #else
      #  A <- A %*% Matrix(expm::sqrtm(M))                               # Eq.(9) -> Eq.(4) ; A = A(I+T) = A(I + [I+S'S]^{-1/2}-I) =
    } else if (prm$method == 'EnOI') {                                #                      = A([I+S'S]^{-1/2}) = AM^{2}
      # do nothing
    } else {
      stop("EnKF: error: assimilate(): method",prm$method," is not defined for the assyncronous EnKF\n")
    }

  if (debugmode) {
    saveRDS(K, file=Kfname)
    rm(K)
  }
  return(list(dx=dx, A=A))
} # end function assimilateLite()


