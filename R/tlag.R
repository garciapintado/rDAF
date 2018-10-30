tlag <- function(x,  dlag = NULL, mlag = NULL, ylag = NULL) {
  # add a positive or negative lag to a POSIXct timestamp
  # x    :: POSIXct
  # dlag :: days lag (integer)
  # mlag :: months lag (integer)
  # ylag :: years lag (integer)

  if (!is.null(dlag)) {
    x <- seq(from=x, by=sign(dlag)*24*3600, length.out=abs(as.integer(dlag))+1)
    x <- x[length(x)]
  }
  if (!is.null(mlag)) {
    if (mlag > 0) {
      x <- seq(from=x, by='month', length.out=abs(as.integer(mlag))+1)
      x <- x[length(x)]
    } else if (mlag < 0) {
      xs <- x - (abs(mlag)+1)*31*24*3600
      dayxtr <- abs(as.numeric(format(xs, format='%d')) - as.numeric(format(x, format='%d')))
      xs <- xs - dayxtr*24*3600
      xseq <- seq(from=xs, by='month', length.out=abs(as.integer(mlag))+5) # 5 to assure series covers end
      xseq <- rev(xseq[xseq <= x])
      x  <- xseq[abs(mlag)+1]
    }
  }
  if (!is.null(ylag)) {
    if (ylag > 0) {
      x <- seq(from=x, by='year', length.out=abs(as.integer(ylag))+1)
      x <- x[length(x)]
    } else if (ylag < 0) {
      xs <- x - (abs(ylag)+1)*366*24*3600
      monxtr <- NULL
      while(is.null(monxtr) || monxtr != 0) {
       monxtr <-  abs(as.numeric(format(xs, format='%m')) - as.numeric(format(x, format='%m')))
       xs <- xs - monxtr*30*24*3600
      }
      dayxtr <- as.numeric(format(xs, format='%d')) - as.numeric(format(x, format='%d'))
      xs <- xs - dayxtr*24*3600
      xseq <- seq(from=xs, by='year', length.out=abs(as.integer(ylag))+50)
      xseq <- rev(xseq[xseq <= x])
      x  <- xseq[abs(ylag)+1]
    }
  }
  return(x)
} # end function tlag
