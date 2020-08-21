#!/usr/bin/env Rscript

########## Setup ##########
options(scipen=999)

########## Load required libraries ##########
require(FITSio)
require(sour)
require(ggplot2)
require(AstroR)

########## XMM Data ##########
broadband <- xmm.pn.lc("data/lightcurves/pn/old/0550640101pn_lccor_300-10000.fits")
sband <- xmm.pn.lc("data/lightcurves/pn/old/0550640101pn_lccor_300-1000.fits")
hband <- xmm.pn.lc("data/lightcurves/pn/old/0550640101pn_lccor_2000-10000.fits")

plot(lc.plot(broadband))

broadband <- subset(broadband, broadband$TIME < 352230000)
sband <- subset(sband, sband$TIME < 352230000)
hband <- subset(hband, hband$TIME < 352230000)

broadband <- subset(broadband, broadband$ERROR < 1)
sband <- subset(sband, sband$ERROR < 1)
hband <- subset(hband, hband$ERROR < 1)

sband <- sband[which(sband$TIME %in% hband$TIME),]
hband <- hband[which(hband$TIME %in% sband$TIME),]
broadband <- broadband[which(broadband$TIME %in% sband$TIME),]
plot(lc.plot.error(broadband))
plot(lc.plot.error(sband))
plot(lc.plot.error(hband))

broadband.binned <- bin.lc(broadband, 100)
sband.binned <- bin.lc(sband, 100)
hband.binned <- bin.lc(hband, 100)

broadband.binned <- subset(broadband.binned, broadband.binned$ERROR < 0.3)
sband.binned <- subset(sband.binned, sband.binned$ERROR < 0.3)
hband.binned <- subset(hband.binned, hband.binned$ERROR < 0.3)
plot(lc.plot.error(broadband.binned))
plot(lc.plot.error(sband.binned))
plot(lc.plot.error(hband.binned))

sband.binned <- sband.binned[which(sband.binned$TIME %in% hband.binned$TIME),]
hband.binned <- hband.binned[which(hband.binned$TIME %in% sband.binned$TIME),]
broadband.binned <- broadband.binned[which(broadband.binned$TIME %in% sband.binned$TIME),]

plot(lc.plot(broadband.binned))

plot(hs.plot(sband.binned,hband.binned))





# ccf
a <- bin.lc(broadband, 1000)
a <- data.frame(t = a$TIME, y = a$RATE, dy = a$ERROR, dt = a$TIMED)
s <- a$t[1]
a$t <- a$t - s
a <- subset(a, a$t < 100000)
a <- a[complete.cases(a),]
aav <- mean(a$y)
am <- max(a$y)
a$y <- (a$y/aav)
xr <- ggplot() +
  geom_linerange(data = a, aes(x = t, ymin = y-dy, ymax = y+dy), color = "black") +
  geom_linerange(data = a, aes(xmin = t-dt, xmax = t+dt, y = y), color = "black") +
  labs(x = "Time [s]", y = "Rate [count/s]") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.text.x = element_text(margin=margin(6,6,10,6,"pt")),
        axis.text.y = element_text(margin=margin(6,6,10,6,"pt")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(xr)


bt <- read.csv("data/lightcurves/uvw1/0550640101_uvw1.lc")
colnames(bt) <- c("t", "y", "dy")
bt <- subset(bt, bt$t < 100000)
bt <- bt[complete.cases(bt),]
# bmin <- min(bt$y)
# print(bmin)
# bt$y <- bt$y - bmin
bav <- mean(bt$y)
print(bav)
# bm <- max(bt$y)
bt$y <- (bt$y - bav) + 1
# bt$y <- bt$y/bav
bdt <- (bt$t[2] - bt$t[1])/2
bdt <- rep(bdt, length(bt$t))
b <- data.frame(t = bt$t, y = bt$y, dy = bt$dy, dt = bdt)
uv <- ggplot() +
  geom_linerange(data = b, aes(x = t, ymin = y-dy, ymax = y+dy), color = "black") +
  geom_linerange(data = b, aes(xmin = t-dt, xmax = t+dt, y = y), color = "black") +
  labs(x = "Time [s]", y = "Rate [count/s]") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.text.x = element_text(margin=margin(6,6,10,6,"pt")),
        axis.text.y = element_text(margin=margin(6,6,10,6,"pt")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(uv)

comp <- ggplot() +
  geom_linerange(data = a, aes(x = t, ymin = y-dy, ymax = y+dy), color = "black") +
  geom_linerange(data = a, aes(xmin = t-dt, xmax = t+dt, y = y), color = "black") +
  geom_linerange(data = b, aes(x = t, ymin = y-dy, ymax = y+dy), color = "red") +
  geom_linerange(data = b, aes(xmin = t-dt, xmax = t+dt, y = y), color = "red") +
  labs(x = "Time [s]", y = "Rate [count/s]") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.text.x = element_text(margin=margin(6,6,10,6,"pt")),
        axis.text.y = element_text(margin=margin(6,6,10,6,"pt")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(comp)

a$t <- a$t - 14000
comp <- ggplot() +
  geom_linerange(data = a, aes(x = t, ymin = y-dy, ymax = y+dy), color = "black") +
  geom_linerange(data = a, aes(xmin = t-dt, xmax = t+dt, y = y), color = "black") +
  geom_linerange(data = b, aes(x = t, ymin = y-dy, ymax = y+dy), color = "red") +
  geom_linerange(data = b, aes(xmin = t-dt, xmax = t+dt, y = y), color = "red") +
  labs(x = "Time [s]", y = "Rate [count/s]") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.text.x = element_text(margin=margin(6,6,10,6,"pt")),
        axis.text.y = element_text(margin=margin(6,6,10,6,"pt")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(comp)

a$t <- a$t + 14000

result <- sour::cross_correlate(a, b, dtau = 100)
plot(result$tau, result$ccf, type = "l", bty = "n", xlab = "time delay [s]", ylab = "CCF")

result <- sour::cross_correlate(a, b, method = "iccf", local.est = TRUE, dtau = 100, nsim = 200000)
hist(result$cent.dist, breaks = 50, col = "steelblue1", main = "", border = NA, prob = TRUE, xlab = "centroid delay (day)")





hard.rat <- hard.ratio(sband.binned, hband.binned)
hard.rat <- hard.rat[complete.cases(hard.rat),]
plot(hr.plot(hard.rat))

broadband.binned <- bin.lc(broadband, 2000)
sband.binned <- bin.lc(sband, 2000)
hband.binned <- bin.lc(hband, 2000)

broadband.binned <- subset(broadband.binned, broadband.binned$ERROR < 0.2)
sband.binned <- subset(sband.binned, sband.binned$ERROR < 0.2)
hband.binned <- subset(hband.binned, hband.binned$ERROR < 0.2)

sband.binned <- sband.binned[which(sband.binned$TIME %in% hband.binned$TIME),]
hband.binned <- hband.binned[which(hband.binned$TIME %in% sband.binned$TIME),]
broadband.binned <- broadband.binned[which(broadband.binned$TIME %in% sband.binned$TIME),]

hard.rat <- hard.ratio(sband.binned, hband.binned)
hard.rat <- hard.rat[complete.cases(hard.rat),]
plot(hr.plot(hard.rat))

ff.df <- prep.ff(sband.binned, hband.binned)
plot(ff.plot(ff.df))

# bins <- unlist(list.files("data/lightcurves/pn/new", pattern = ".*lccor.*"))
# energy <- rep(0, length(bins))
# energy.error <- rep(0, length(bins))
# p.fvar <- rep(0, length(bins))
# p.fvar.error <- rep(0, length(bins))
# v.fvar <- rep(0, length(bins))
# v.fvar.error <- rep(0, length(bins))
# e.fvar <- rep(0, length(bins))
# e.fvar.error <- rep(0, length(bins))
# 
# tlc <- xmm.lc("data/lightcurves/pn/new/0550640101pn_lccor_1024-1220.fits")
# tlc <- subset(tlc, tlc$TIME < 352230000)
# tlc <- tlc[complete.cases(tlc),]
# tlc <- bin.lc(tlc, 200)
# tlc <- subset(tlc, tlc$ERROR < 0.15)
# plot(lc.plot(tlc))
# test <- fvar.edelson(tlc)
# print(mean(tlc$ERROR^2))
# print((var(tlc$RATE)-mean(tlc$ERROR^2))/(mean(tlc$RATE)^2))
# 
# i <- 1
# for (bin in bins) {
#   lc.path <- paste(list("data/lightcurves/pn/new/",bin), collapse = "")
#   print(lc.path)
#   lc <- xmm.lc(lc.path)
#   lc <- lc[complete.cases(lc),]
#   lc <- bin.lc(lc, 1000)
#   lc <- lc[complete.cases(lc),]
#   fvar.date <- fvar.edelson(lc)
#   fvar.datv <- fvar.vaughan(lc)
#   fvar.datp <- fvar.ponti(lc)
#   r <- sub("0550640101pn_lccor_", "", bin)
#   r <- sub(".fits", "", r)
#   e <- as.numeric(unlist(strsplit(r, split = "-", fixed = TRUE)))
#   energy[i] <- mean(e)
#   energy.error <- ((e[2]-e[1])/2)
#   p.fvar[i] <- fvar.datp[1]
#   p.fvar.error[i] <- fvar.datp[2]
#   v.fvar[i] <- fvar.datv[1]
#   v.fvar.error[i] <- fvar.datv[2]
#   e.fvar[i] <- fvar.date[1]
#   e.fvar.error[i] <- fvar.date[2]
#   i <- i + 1
# }
# fvar.df <- data.frame(energy, energy.error, e.fvar, e.fvar.error, v.fvar, v.fvar.error, p.fvar, p.fvar.error)
# colnames(fvar.df) <- c("ENERGY", "ENERGYD", "EFV", "EFE", "VFV", "VFE", "PFV", "PFE")


########## Suzaku XIS ##########
xi03.broadband <- suzaku.lc("data/lightcurves/xis/xi03_500-10000_sc.fits", "data/lightcurves/xis/xi03_500-10000_bg.fits")
xi1.broadband <- suzaku.lc("data/lightcurves/xis/xi1_500-10000_sc.fits", "data/lightcurves/xis/xi1_500-10000_bg.fits")
xi03.sband <- suzaku.lc("data/lightcurves/xis/xi03_500-1000_sc.fits", "data/lightcurves/xis/xi03_500-1000_bg.fits")
xi03.hband <- suzaku.lc("data/lightcurves/xis/xi03_2000-10000_sc.fits", "data/lightcurves/xis/xi03_2000-10000_bg.fits")

xi03.sband <- xi03.sband[complete.cases(xi03.sband),]
xi03.hband <- xi03.hband[complete.cases(xi03.hband),]
xi03.broadband <- xi03.broadband[complete.cases(xi03.broadband),]
xi1.broadband <- xi1.broadband[complete.cases(xi1.broadband),]

xi03.sband <- xi03.sband[which(xi03.sband$TIME %in% xi03.hband$TIME),]
xi03.hband <- xi03.hband[which(xi03.hband$TIME %in% xi03.sband$TIME),]
xi03.broadband <- xi03.broadband[which(xi03.broadband$TIME %in% xi03.sband$TIME),]
xi1.broadband <- xi1.broadband[which(xi1.broadband$TIME %in% xi03.broadband$TIME),]

xi03.broadband.binned <- bin.lc(xi03.broadband, 5760)
xi1.broadband.binned <- bin.lc(xi1.broadband, 5760)
xi03.sband.binned <- bin.lc(xi03.sband, 5760)
xi03.hband.binned <- bin.lc(xi03.hband, 5760)

plot(lc.plot.error(xi03.broadband.binned, 0.001))
plot(lc.plot.error(xi1.broadband.binned, 0.001))
plot(lc.plot.error(xi03.sband.binned, 0.001))
plot(lc.plot.error(xi03.hband.binned, 0.001))

plot(lc.plot(xi03.broadband.binned))
plot(lc.plot(xi1.broadband.binned))
plot(hs.plot(xi03.sband.binned, xi03.hband.binned))

xi03.hard.rat <- hard.ratio(xi03.sband.binned, xi03.hband.binned)
xi03.hard.rat <- xi03.hard.rat[complete.cases(xi03.hard.rat),]
plot(hr.plot(xi03.hard.rat))

xi03.ff.df <- prep.ff(xi03.sband.binned, xi03.hband.binned)
plot(ff.plot(xi03.ff.df))
