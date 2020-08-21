#' @title Rebin light curve data frame
#' @description Takes the input data frame light curve and rebins in time accordingly.
#' @author Adam Gonzalez
#' @param lightcurve data frame of light curve
#' @param dt time bin size in seconds
#' @return Data frame of the rebinned light curve.
#' @examples binned.lightcurve <- lc.rebin(lightcurve, 100)
#' @export


lc.rebin <- function(lightcurve, dt){
  lightcurve <- subset(lightcurve, is.nan(lightcurve$RATE) == FALSE)
  x <- lightcurve$TIME
  y.vals <- lightcurve$RATE
  y.errs <- lightcurve$RATE.ERR
  b.vals <- lightcurve$BKG
  b.errs <- lightcurve$BKG.ERR
  t0 <- x[1]
  x <- x - t0
  new.x <- seq(from=min(x), to=max(x), by=dt)
  y.vals.bin <- rep(0, length(new.x))
  y.errs.bin <- rep(0, length(new.x))
  b.vals.bin <- rep(0, length(new.x))
  b.errs.bin <- rep(0, length(new.x))

  for (j in 1:(length(new.x)-1)){
    y.count <- 0
    b.count <- 0
    y.vals.sum <- 0
    y.errs.sum <- 0
    b.vals.sum <- 0
    b.errs.sum <- 0
    for (i in 1:length(x)){
      if (new.x[j] <= x[i]){
        if (x[i] < new.x[j+1]){
          y.count <- y.count + 1
          b.count <- b.count + 1
          y.vals.sum <- y.vals.sum + y.vals[i]
          y.errs.sum <- y.errs.sum + y.errs[i]^2
          b.vals.sum <- b.vals.sum + b.vals[i]
          b.errs.sum <- b.errs.sum + b.errs[i]^2
        }
      }
    }
    y.vals.bin[j] <- y.vals.sum / y.count
    y.errs.bin[j] <- sqrt(y.errs.sum) / y.count
    b.vals.bin[j] <- b.vals.sum / b.count
    b.errs.bin[j] <- sqrt(b.errs.sum) / b.count
  }

  return(data.frame(TIME=new.x[-length(new.x)]+t0, TIME.ERR=dt/2, RATE=y.vals.bin[-length(new.x)], RATE.ERR=y.errs.bin[-length(new.x)], BKG=b.vals.bin[-length(new.x)], BKG.ERR=b.errs.bin[-length(new.x)]))
}
