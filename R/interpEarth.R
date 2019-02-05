interpEarth <- function(LON, LAT, z, xo, yo, output = 'grid', retnum = TRUE, patchwin=90) {
  # irregular points onto regular grid around the Earth by use of akima interpolation
  # require(interp)

  if (any(LON < -180) || any(LON > 360))
    stop('interpEarth: input longitudes out of [-180,360]')
  if (any(xo < -180) || any(xo > 360))
    stop('interpEarth: output longitudes out of [-180,360]')
  if (diff(range(LON)) > 360)
    stop('interpEarth: input longitudes expand more than 360 degrees')
  if (diff(range(xo)) > 360)
    stop('interpEarth: output longitudes expand more than 360 degrees')
  if (!(output %in% c('grid','points')))
    stop('interpEarth: output must be grid or points')

  if (any(diff(xo) <= 0))
    stop('xo not monotonically increasing')
  if (any(diff(yo) <= 0))
    stop('yo not monotonically increasing')
  
  # express input and output longitudes centered at 180 (in [0,360])
  LON360 <- LON
  if (any(LON < 0))                                     # input in [-180,0)
    LON360[LON < 0] <- LON[LON < 0] + 360               # in [0,360]  

  xo360 <- xo
  if (any(xo < 0))                                     # output in [-180,0)
    xo360[xo < 0] <- xo[xo < 0] + 360                  # in [0,360]  
  
  # recycle patchwin at each side
  lid <- which(LON360 > (360 - patchwin))  # patch to the left
  rid <- which(LON360 < patchwin)          # patch to the right
  LONaug <- c(LON360[lid]-360, LON360, LON360[rid]+360)
  LATaug <- c(LAT[lid],        LAT,    LAT[rid]) 
  zaug   <- c(z[lid],          z,      z[rid])

  if (output != 'grid')
    zo <- interp::interp(LONaug, LATaug, zaug, xo=xo360, yo=yo, output=output)$z
  else {
    zo <- matrix(NA,length(xo),length(yo))
    if (any(xo < 0))
      zo[xo <  0,] <-  akima::interp(LONaug, LATaug, zaug, xo=xo[xo < 0]+360, yo=yo)$z
    if (any(xo >= 0))
      zo[xo >= 0,] <-  akima::interp(LONaug, LATaug, zaug, xo=xo[xo >= 0], yo=yo)$z
  }  
  if (output == 'grid' && retnum)                            # write as numeric from ul corner and advancing by row=TRUE
    zo <- as.numeric(zo[,ncol(zo):1])
  return(zo)
} # end function interpEarth()
