#!/usr/bin/env Rscript

########## Setup ##########
options(scipen=999)

########## Load required libraries ##########
require(FITSio)
require(sour)
require(ggplot2)
require(AstroR)

########## XMM Data ##########
bin.width = 1000

########## XMM Data ##########
band.all <- xmm.pn.lc("data/lightcurves/pn/old/0550640101pn_lccor_300-10000.fits")
band.soft <- xmm.pn.lc("data/lightcurves/pn/old/0550640101pn_lccor_300-1000.fits")
band.hard <- xmm.pn.lc("data/lightcurves/pn/old/0550640101pn_lccor_2000-10000.fits")
band.uv <- read.csv("data/lightcurves/uvw1/0550640101_uvw1.lc")
colnames(band.uv) <- c("TIME", "RATE", "ERROR")
band.uv$TIMED <- rep((band.uv$TIME[2] - band.uv$TIME[1])/2, length(band.uv$TIME))
t <- band.all$TIME[1]
band.all$TIME <- band.all$TIME - t
t <- band.soft$TIME[1]
band.soft$TIME <- band.soft$TIME - t
t <- band.hard$TIME[1]
band.hard$TIME <- band.hard$TIME - t

band.all <- subset(band.all, band.all$TIME < 100000)
band.soft <- subset(band.soft, band.soft$TIME < 100000)
band.hard <- subset(band.hard, band.hard$TIME < 100000)
band.uv <- subset(band.uv, band.uv$TIME < 100000)

band.all.binned <- bin.lc(band.all, bin.width)
band.soft.binned <- bin.lc(band.soft, bin.width)
band.hard.binned <- bin.lc(band.hard, bin.width)

band.all.binned <- band.all.binned[complete.cases(band.all.binned),]
band.soft.binned <- band.soft.binned[complete.cases(band.soft.binned),]
band.hard.binned <- band.hard.binned[complete.cases(band.hard.binned),]

band.all.binned <- subset(band.all.binned, band.all.binned$ERROR < as.numeric(quantile(band.all.binned$ERROR, 0.95)))
band.soft.binned <- subset(band.soft.binned, band.soft.binned$ERROR < as.numeric(quantile(band.soft.binned$ERROR, 0.95)))
band.hard.binned <- subset(band.hard.binned, band.hard.binned$ERROR < as.numeric(quantile(band.hard.binned$ERROR, 0.95)))

plot(lc.plot.error(band.all.binned, bin.width = 0.005))
plot(lc.plot.error(band.soft.binned, bin.width = 0.005))
plot(lc.plot.error(band.hard.binned, bin.width = 0.005))

band.uv.shifted1 <- data.frame(band.uv)
band.uv.shifted1$TIME <- band.uv.shifted1$TIME + (-10000)
band.uv.shifted2 <- data.frame(band.uv)
band.uv.shifted2$TIME <- band.uv.shifted2$TIME + (-10000)

plot(lc.plot(band.soft.binned))
plot(lc.plot(band.hard.binned))
plot(lc.plot(band.uv))

plot(lc.overplot(band.all.binned, band.uv, names = c("0.3-10 keV", "UVW1")))
plot(lc.overplot(band.all.binned, band.uv.shifted1, names = c("0.3-10 keV", "UVW1")))
plot(lc.overplot(band.all.binned, band.uv.shifted2, names = c("0.3-10 keV", "UVW1")))

plot(lc.overplot(band.hard.binned, band.uv, names = c("2-10 keV", "UVW1")))
plot(lc.overplot(band.hard.binned, band.uv.shifted1, names = c("2-10 keV", "UVW1")))
plot(lc.overplot(band.hard.binned, band.uv.shifted2, names = c("2-10 keV", "UVW1")))

plot(lc.overplot(band.soft.binned, band.uv, names = c("0.3-1 keV", "UVW1")))
plot(lc.overplot(band.soft.binned, band.uv.shifted1, names = c("0.3-1 keV", "UVW1")))
plot(lc.overplot(band.soft.binned, band.uv.shifted2, names = c("0.3-1 keV", "UVW1")))

plot(lc.overplot(band.soft.binned, band.hard.binned, names = c("0.3-1 keV", "2-10 keV")))

dcf.lc1 <- data.frame(t = band.soft.binned$TIME, y = band.soft.binned$RATE, dy = band.soft.binned$ERROR)
dcf.lc2 <- data.frame(t = band.hard.binned$TIME, y = band.hard.binned$RATE, dy = band.hard.binned$ERROR)
result <- sour::cross_correlate(dcf.lc1, dcf.lc2, dtau = bin.width)
plot(result$tau, result$ccf, type = "l", bty = "n", xlab = "time delay [s]", ylab = "CCF")



dcf.lc1 <- data.frame(t = band.hard.binned$TIME, y = band.hard.binned$RATE, dy = band.hard.binned$ERROR)
dcf.lc1$t <- dcf.lc1$t / 1000
dcf.lc2 <- data.frame(t = band.uv$TIME, y = band.uv$RATE, dy = band.uv$ERROR)
dcf.lc2$t <- dcf.lc2$t / 1000
ccf.result <- sour::cross_correlate(dcf.lc1, dcf.lc2, method="dcf", dtau = bin.width/1000)
ccf <- data.frame(ccf.result$tau, ccf.result$ccf)


plot(ccf.result$tau, ccf.result$ccf, type = "l", bty = "n", xlab = "time delay [ks]", ylab = "CCF")

# result <- sour::cross_correlate(dcf.lc1, dcf.lc2, method = "iccf", local.est = TRUE, dtau = bin.width, nsim = 10000)
# hist(result$cent.dist, breaks = 50, col = "steelblue1", main = "", border = NA, prob = TRUE, xlab = "centroid delay (day)")
