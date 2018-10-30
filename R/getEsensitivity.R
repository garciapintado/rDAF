getEsensitivity <- function(X, Y, useEigen=FALSE, eigop=TRUE) {
  # get ensemble sensitivity of the observation space to control variables
  # the definition of ensemble sentivity here matches that of
  # Hakim & Torn (2008)
  # Chen & Oliver (2011)
  # Ancell & Hakim (2007) Appendix (what they call "ensemble estimation of adjoint sensitivity")
  # This differs from the simplified linear-regression sensitivities of Ancell & Hakim (2007), where covariances in the
  # input variables are included in the sensitivity estimates.
  #require(Matrix)

  dX <- X - rowMeans(X)                                                  # input perturbations : [d,m]
  dY <- Y - rowMeans(Y)                                                  # output perturbations: [p,m]
  d  <- nrow(dX)
  p  <- nrow(dY)
  m  <- ncol(dX)

  if (is.null(useEigen)) {
  if (m >= d)
      useEigen <- FALSE
  else
      useEigen <- TRUE
  }
  if (!useEigen) {                                                       # singular value decomposition
    #dX.svd <- svd(dX, nu=nrow(dX),nv=ncol(dX))
    dX.svd <- svd(dX)                                                    # U S Vt
    #cat('dim(U):',dim(dX.svd$u),'\n')
    #dZ <- dX.svd$u %*% diag(dX.svd$d) %*% t(dX.svd$v)                   # dZ == dX
    svI <- 1/dX.svd$d
    if (m <= d)
      svI[m:length(svI)] <- 0
    svMI <- matrix(0,ncol(dX.svd$v),ncol(dX.svd$u))
    diag(svMI) <- svI
    dXI <- dX.svd$v %*% svMI %*% t(dX.svd$u)                             # [m,m][m,u][u,u]  pseudo-inverse
    G  <- dY %*% dXI                                                     # [p,d]
  } else {                                                               # eigenvalue decomposition
    crossX <- crossprod(dX)                                              # [m,m]
    CL <- eigen(crossX, symmetric=TRUE)
    #crossX.kappa <- kappa(crossX, exact=TRUE)
    crossX.rcond <- rcond(crossX)
    crossX.rank  <- rankMatrix(crossX)
    #cat('crossX.kappa:',crossX.kappa,'\n')
    #cat('crossX.rcond:',crossX.rcond,'\n')                               # -> 1 fine | -> 0 bad
    #cat('crossX.rank:',crossX.rank,'\n')

    # crossZ <- CL$vectors %*% diag(CL$values) %*% t(CL$vectors)         # == crossX
    eigvI <- 1/CL$values                                                 # [d]
    if (m > d)
      eigvI[(d+1):m] <- 0
    else
      eigvI[m] <- 0
    eigvI[abs(CL$values[1]/CL$values) > 1.0E10] <- 0

    if (eigop) {  # default for low-rank: after (corrected) Hakim & Torn (2008) Eq(13)
      crossXI <- CL$vectors %*% diag(eigvI) %*% t(CL$vectors)                          # [m,m][m,m][m,m]
      if (p < d)
        G <- tcrossprod( dY %*% crossXI , dX)                                          # [p,d] = ([p,m][m,m])[m,d]
      else
        G <- dY %*% tcrossprod(crossXI , dX)                                           # [p,d] = [p,m] ([m,m][m,d])
    } else {                                                      # eigop != 1 Pxx explicitly calculated
      cat('warning: eigop=',eigop,' | full Pxx explicitly calculated\n')
      XV <- dX %*% CL$vectors                                     # not practical, just to follow:
      PexxI <- tcrossprod(XV %*% diag(eigvI^2), XV)               # [d,d] (A4) in Ancell and Hakim (2007)
      G <- dY %*% t(dX) %*% PexxI                                 # m > d => all.equal(GeSV,Gt) == TRUE
    }
  }
  return(G)
}  # end function getEsensitivity()

