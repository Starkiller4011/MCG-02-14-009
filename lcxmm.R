#' @title Convert epiclcorr .fits file into a data frame
#' @description Reads in the corrected light curve file from XMM SAS epiclccorr task.
#' @param filename .fits file from epiclcorr
#' @param tzero if TRUE will set the output light curve to start at 0 seconds
#' @return Data frame of the light curve with structure: time, dt, source, error, background, error.
#' @examples lightcurve <- lc.xmm(lccor.fits)
#' @export


lc.xmm <- function(filename, tzero = FALSE){
  lc.fits <- FITSio::readFITS(filename)
  lc <- data.frame(lc.fits$col[1], lc.fits$col[1], lc.fits$col[2], lc.fits$col[3], lc.fits$col[5], lc.fits$col[6])
  colnames(lc) <- c("TIME", "TIME.ERR", "RATE", "RATE.ERR", "BKG", "BKG.ERR")
  if (tzero == TRUE){
    lc$TIME <- lc$TIME - lc$TIME[1]
  }
  lc$TIME.ERR <- (lc$TIME[2]-lc$TIME[1])/2
  return(lc)
}
