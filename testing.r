#!/usr/bin/env Rscript

########## Setup ##########
options(scipen=999)

########## Load required libraries ##########
require(FITSio)
require(ggplot2)
require(AstroR)

########## Load XMM Data ##########
bin.width = 1500
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

########## Clean Data ##########
band.all <- subset(band.all, band.all$TIME < 100000)
band.soft <- subset(band.soft, band.soft$TIME < 100000)
band.hard <- subset(band.hard, band.hard$TIME < 100000)
band.uv <- subset(band.uv, band.uv$TIME < 100000)

########## Bin Data ##########
band.all.binned <- bin.lc(band.all, bin.width)
band.soft.binned <- bin.lc(band.soft, bin.width)
band.hard.binned <- bin.lc(band.hard, bin.width)
band.all.binned <- band.all.binned[complete.cases(band.all.binned),]
band.soft.binned <- band.soft.binned[complete.cases(band.soft.binned),]
band.hard.binned <- band.hard.binned[complete.cases(band.hard.binned),]
band.all.binned <- subset(band.all.binned, band.all.binned$ERROR < as.numeric(quantile(band.all.binned$ERROR, 0.95)))
band.soft.binned <- subset(band.soft.binned, band.soft.binned$ERROR < as.numeric(quantile(band.soft.binned$ERROR, 0.95)))
band.hard.binned <- subset(band.hard.binned, band.hard.binned$ERROR < as.numeric(quantile(band.hard.binned$ERROR, 0.95)))

########## Plot Error distributions ##########
plot(lc.plot.error(band.all.binned, bin.width = 0.001))
plot(lc.plot.error(band.soft.binned, bin.width = 0.001))
plot(lc.plot.error(band.hard.binned, bin.width = 0.001))

########## Plot light curves ##########
plot(lc.plot(band.all.binned))
plot(lc.plot(band.soft.binned))
plot(lc.plot(band.hard.binned))
plot(lc.plot(band.uv))

