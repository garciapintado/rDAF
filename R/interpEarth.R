interpEarth <- function(LON, LAT, z, xo, yo, output = 'grid', retnum = TRUE) {
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
  
  # express input and output longitudes in [0,360]
  LON360 <- LON
  if (any(LON360 < 0))                                     # input in [-180,0)
    LON360[LON360 < 0] <- LON360[LON360 < 0] + 360         # in [0,360]  

  xo360 <- xo
  if (any(xo360 < 0))                                      # output in [-180,0)
    xo360[xo360 < 0]   <- xo360[xo360 < 0] + 360           # in [0,360]  

  
  # augment source data by duplicating points for ranges [-90,0] & [360,450]  
  LONaug <- LON
  LATaug <- LAT
  zaug   <- z
  
  ldid   <- which(LON360 > 180)                            # left duplicated indices
  LONaug <- c(LONaug,LON360[ldid]-360)
  LATaug <- c(LATaug,LAT[ldid])
  zaug   <- c(zaug,zaug[ldid])

  rdid   <- which(LON360 < 180)                            # right duplicated indices
  LONaug <- c(LONaug,LON360[rdid]+360)
  LATaug <- c(LATaug,LAT[rdid])
  zaug   <- c(zaug,zaug[rdid])
  
  if (output != 'grid')
    zo <- interp::interp(LONaug, LATaug, zaug, xo=xo360, yo=yo, output=output)$z
  else {
   # if (requireNamespace("akima")) 
      zo <- akima::interp(LONaug, LATaug, zaug, xo=xo360, yo=yo)$z  # much faster
   # else
   #   stop('akima::interp not available')
  }  
  if (output == 'grid' && retnum)                            # write as numeric from ul corner and advancing by row=TRUE
    zo <- as.numeric(zo[,ncol(zo):1])
  return(zo)
} # end function interpEarth()
