interppEarth <- function(LON, LAT, z, xo, yo) {
  # irregular onto irregular points around the Earth by use of akima
  # all input/output longitude in [0,360] 
  # interpolation is divided into 3 chunks:
  #  left     center    right
  # [0,90]   [90,270] [270,360]

  #require(akima)
  LON180 <- LON 
  LON180[LON180 > 180] <- LON180[LON180 > 180] - 360           # in [-180,180]

  no <- length(xo)
  if (length(yo) != no)
   stop('interppEarth --ERR001--')
  if (any(xo < 0) || any(xo > 360)) 
   stop('interppEarth --ERR002--')

  zo <- rep(NA,no)

  oboo <- xo < 90                                                               # [0,90]    LON
  if (sum(oboo) > 0)
    zo[oboo] <- interpp(LON180, LAT, z, xo=xo[oboo],       yo=yo[oboo])$z

  oboo <- xo>=90 & xo <= 270                                                    # [90,270]  LON
  if (sum(oboo) > 0) 
    zo[oboo] <- interpp(LON,    LAT, z, xo=xo[oboo],       yo=yo[oboo])$z    

  oboo <- xo > 270                                                              # [270,360] LON
  if (sum(oboo) > 0) 
    zo[oboo] <- interpp(LON180, LAT, z, xo=xo[oboo] - 360, yo=yo[oboo])$z    
 
  return(zo)
}
