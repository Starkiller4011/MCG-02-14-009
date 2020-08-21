#' @title Compute the cross spectrum
#' @description Takes the input data frames of two light curves and computes the cross spectrum between them.
#' @param lightcurve1 data frame of the first light curve
#' @param lightcurve2 data frame of the second light curve
#' @return Data frame of the cross spectrum.
#' @examples raw.cross <- crossspec(lightcurve1, lightcurve2)
#' @export
crossspec <- function(lightcurve1, lightcurve2){
  if (lightcurve1$TIME[1] != 0){
    lightcurve1$TIME <- lightcurve1$TIME - lightcurve1$TIME[1]
  }
  if (lightcurve2$TIME[1] != 0){
    lightcurve2$TIME <- lightcurve2$TIME - lightcurve2$TIME[1]
  }
  dt <- lightcurve1$TIME[2] - lightcurve1$TIME[1]
  N <- length(lightcurve1$TIME)
  src1.mean <- mean(lightcurve1$RATE)
  bkg1.mean <- mean(lightcurve1$BACKV)
  noise1 <- 2 * (src1.mean + bkg1.mean) / (src1.mean)^2
  dft1 <- fft(lightcurve1$RATE)
  dft1 <- dft1[-1]
  dft1 <- dft1[1:(N/2)]
  power1 <- (abs(dft1)^2) * (2 * dt) / (N * src1.mean^2)
  src2.mean <- mean(lightcurve2$RATE)
  bkg2.mean <- mean(lightcurve2$BACKV)
  noise2 <- 2 * (src2.mean + bkg2.mean) / (src2.mean)^2
  dft2 <- fft(lightcurve2$RATE)
  dft2 <- dft2[-1]
  dft2 <- dft2[1:(N/2)]
  power2 <- (abs(dft2)^2) * (2 * dt) / (N * src2.mean^2)
  freq <- (1:length(lightcurve1$TIME))/(N*dt)
  freq <- freq[1:(N/2)]
  cross <- Conj(dft2) * dft1 * (2*dt) / (N*src1.mean*src2.mean)
  return(data.frame(FRQ=freq, POW.1=power1, PN.1=noise1, POW.2 = power2, PN.2=noise2, CROSS=cross))
}

#' @title Compute the lag-frequency spectrum
#' @description Takes the input data frames of a cross spectrum and computes the lag-frequency spectrum.
#' @param lightcurve1 data frame of the first light curve
#' @param lightcurve2 data frame of the second light curve
#' @param gfac geometric binning factor (default 1.5)
#' @return Data frame of the lag-frequency spectrum.
#' @examples spec <- lagfrq(lightcurve1, lightcurve2, gfac = 1.5)
#' @export
lagfrq2 <- function(cross.data, gfac = 1.4){
  # cross.data <- crossspec(lightcurve1, lightcurve2)
  frq.bins <- min(cross.data$FRQ)
  while (max(frq.bins)*gfac <= max(cross.data$FRQ)){
    frq.bins <- c(frq.bins, max(frq.bins)*gfac)
  }
  frq.bins <- subset(frq.bins, frq.bins >= 3*min(cross.data$FRQ))
  mid.frqs <- 0
  for (i in 1:(length(frq.bins)-1)){
    mid.frqs <- c(mid.frqs, 10^( mean( c(log10(frq.bins[i]),log10(frq.bins[i+1])) ) ))
  }
  mid.frqs <- mid.frqs[-1]
  neg.frqs <- mid.frqs - frq.bins[1:(length(frq.bins)-1)]
  pos.frqs <- frq.bins[2:length(frq.bins)] - mid.frqs
  cross.data.bin <- data.frame(FRQ = mid.frqs,
                               FRQ.NEG = neg.frqs,
                               FRQ.POS = pos.frqs,
                               POW.1 = rep(0, length.out = length(mid.frqs)),
                               PN.1 = rep(cross.data$PN.1, length.out = length(mid.frqs)),
                               POW.2 = rep(0, length.out = length(mid.frqs)),
                               PN.2 = rep(cross.data$PN.2, length.out = length(mid.frqs)),
                               CROSS = rep(0, length.out = length(mid.frqs)),
                               COUNT = rep(0, length.out = length(mid.frqs)))
  for (k in 1:length(mid.frqs)){
    count <- 0
    psum.1 <- 0.0
    psum.2 <- 0.0
    csum <- 0.0
    for (m in 1:length(cross.data$FRQ)){
      if (frq.bins[k] <= cross.data$FRQ[m] && cross.data$FRQ[m] < frq.bins[k+1]){
        count <- count + 1
        psum.1 <- psum.1 + cross.data$POW.1[m]
        psum.2 <- psum.2 + cross.data$POW.2[m]
        csum <- csum + cross.data$CROSS[m]
      }
    }
    cross.data.bin$POW.1[k] <- psum.1/count
    cross.data.bin$POW.2[k] <- psum.2/count
    cross.data.bin$CROSS[k] <- csum/count
    cross.data.bin$COUNT[k] <- count
  }
  nsq <- 0
  coh <- rep(0, length.out = length(mid.frqs))
  lag <- rep(0, length.out = length(mid.frqs))
  lag.err <- rep(0, length.out = length(mid.frqs))
  for (i in 1:length(mid.frqs)){
    lag[i] <- Arg(cross.data.bin$CROSS[i]) / (2*pi*mid.frqs[i])
    nsq <- ((cross.data.bin$POW.1[i]-cross.data.bin$PN.1[i])*cross.data.bin$PN.2[i] + (cross.data.bin$POW.2[i]-cross.data.bin$PN.2[i])*cross.data.bin$PN.1[i] + cross.data.bin$PN.1[i]*cross.data.bin$PN.2[i]) / cross.data.bin$COUNT[i]
    if (nsq < abs(cross.data.bin$CROSS[i]^2)){
      coh[i] <- (abs(cross.data.bin$CROSS[i])^2 - nsq) / (cross.data.bin$POW.1[i]*cross.data.bin$POW.2[i])
      lag.err[i] <- sqrt( (1-coh[i]) / (2*coh[i]*cross.data.bin$COUNT[i]) ) / (2*pi*mid.frqs[i])
    } else {# (nsq >= abs(cross.data.bin$CROSS[i])^2){
      coh[i] <- NA_real_
      lag.err[i] <- NA_real_
    }
  }
  return(data.frame(FRQ = mid.frqs, FRQ.NEG = neg.frqs, FRQ.POS = pos.frqs, COH = coh, LAG = lag, LAG.ERR = lag.err))
}