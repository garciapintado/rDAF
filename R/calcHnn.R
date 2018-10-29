calcHnn <- function(G=NULL, gLON=NULL, gLAT=NULL, pos) {
  #
  # simple nearest neighbour association from points to a grid
  #
  # G         :: gmeta6 S3 - faster than giving gLON, gLAT for regular grids
  # gLON,gLAT :: coordinates for irregular grids
  # pos  :: 2-col matrix with locations
  #
  #
  #
  # File:           calcHnn.R
  #
  # Created:        15/01/2016
  #
  # Last modified:  15/01/2016
  #
  #
  # Author:         Javier Garcia-Pintado
  #                 ESSC [NCEO]
  #
  # get H [as a dgCMatrix] from single pixels nearest to observations
  #
  #require(Matrix)

  p  <- nrow(pos)

  if (!is.null(G) && !is.null(gLON))
    stop('calcHnn: use just either G [regular grids] or gLON,gLAT [irregular grids] for horizontal grid definition')

  if (!is.null(G)) {
    ng  <- G$cells                                                              # cells in each horizontal grid
    Gid <- nnGpos(G,pos)                                                        # in hydrosim
  } else {
    ng  <- length(gLON)
    Gid <- n2dist(cbind(as.numeric(gLON),as.numeric(gLAT)),pos)$neighs         # in splancs 
  }

  dgTH <- spMatrix(p,ng,i=1:p,j=Gid,x=rep(1.0,p))           # dgTMatrix
  dgCH <- as(dgTH,'dgCMatrix')                              # coerce to dgCMatrix
  return(dgCH)
} # function calcHnn
