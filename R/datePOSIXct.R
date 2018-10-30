datePOSIXct <- function(x, dlag=NULL, mlag=NULL, ylag=NULL) {
  # day of year CESM envrun$RUN_STARTDATE format
  # x    :: POSIXct
  # dlag :: extra lag in days
  # mlag :: extra lag in months
  # ylag :: extra lag in years

  x <- tlag(x, dlag=dlag, mlag=mlag, ylag=ylag)
  x <- format(x, format='%Y-%m-%d')
  return(x)
}
