#!/usr/bin/env Rscript

########## Setup ##########
options(scipen=999)

########## Load required libraries ##########
require(FITSio)
require(sour)
require(ggplot2)
require(AstroR)

########## Set Global Constants ##########
xmm.bin.width = 1500
suz.bin.width = 5000
nsims <- 1000
beta <- 1.87

########## Load XMM Data ##########
xmm.band.all <- xmm.pn.lc("data/lightcurves/pn/old/0550640101pn_lccor_300-10000.fits")
xmm.band.soft <- xmm.pn.lc("data/lightcurves/pn/old/0550640101pn_lccor_300-1000.fits")
xmm.band.hard <- xmm.pn.lc("data/lightcurves/pn/old/0550640101pn_lccor_2000-10000.fits")
xmm.band.uv <- read.csv("data/lightcurves/uvw1/0550640101_uvw1.lc")
colnames(xmm.band.uv) <- c("TIME", "RATE", "ERROR")
xmm.band.uv$BACKV <- 0
xmm.band.uv$BACKE <- 0
xmm.band.uv$TIMED <- rep((xmm.band.uv$TIME[2] - xmm.band.uv$TIME[1])/2, length(xmm.band.uv$TIME))
t <- xmm.band.all$TIME[1]
xmm.band.all$TIME <- xmm.band.all$TIME - t
t <- xmm.band.soft$TIME[1]
xmm.band.soft$TIME <- xmm.band.soft$TIME - t
t <- xmm.band.hard$TIME[1]
xmm.band.hard$TIME <- xmm.band.hard$TIME - t

########## Load Suzaku Data ##########
suz.band.all <- suzaku.lc("data/lightcurves/xis/xi03_500-10000_sc.fits", "data/lightcurves/xis/xi03_500-10000_bg.fits")
suz.band.soft <- suzaku.lc("data/lightcurves/xis/xi03_500-1000_sc.fits", "data/lightcurves/xis/xi03_500-1000_bg.fits")
suz.band.hard <- suzaku.lc("data/lightcurves/xis/xi03_2000-10000_sc.fits", "data/lightcurves/xis/xi03_2000-10000_bg.fits")

########## Clean Data ##########
xmm.band.all <- subset(xmm.band.all, xmm.band.all$TIME < 100000)
xmm.band.soft <- subset(xmm.band.soft, xmm.band.soft$TIME < 100000)
xmm.band.hard <- subset(xmm.band.hard, xmm.band.hard$TIME < 100000)
xmm.band.uv <- subset(xmm.band.uv, xmm.band.uv$TIME < 100000)
suz.band.all <- subset(suz.band.all, suz.band.all$TIME < 100000)
suz.band.soft <- subset(suz.band.soft, suz.band.soft$TIME < 100000)
suz.band.hard <- subset(suz.band.hard, suz.band.hard$TIME < 100000)

########## Bin Data ##########
xmm.band.uv.binned <- xmm.band.uv
xmm.band.all.binned <- bin.lc(xmm.band.all, xmm.bin.width)
xmm.band.soft.binned <- bin.lc(xmm.band.soft, xmm.bin.width)
xmm.band.hard.binned <- bin.lc(xmm.band.hard, xmm.bin.width)
suz.band.all.binned <- bin.lc(suz.band.all, suz.bin.width)
suz.band.soft.binned <- bin.lc(suz.band.soft, suz.bin.width)
suz.band.hard.binned <- bin.lc(suz.band.hard, suz.bin.width)
xmm.band.all.binned <- xmm.band.all.binned[complete.cases(xmm.band.all.binned),]
xmm.band.soft.binned <- xmm.band.soft.binned[complete.cases(xmm.band.soft.binned),]
xmm.band.hard.binned <- xmm.band.hard.binned[complete.cases(xmm.band.hard.binned),]
suz.band.all.binned <- suz.band.all.binned[complete.cases(suz.band.all.binned),]
suz.band.soft.binned <- suz.band.soft.binned[complete.cases(suz.band.soft.binned),]
suz.band.hard.binned <- suz.band.hard.binned[complete.cases(suz.band.hard.binned),]
xmm.band.all.binned <- subset(xmm.band.all.binned, xmm.band.all.binned$ERROR < as.numeric(quantile(xmm.band.all.binned$ERROR, 0.95)))
xmm.band.soft.binned <- subset(xmm.band.soft.binned, xmm.band.soft.binned$ERROR < as.numeric(quantile(xmm.band.soft.binned$ERROR, 0.95)))
xmm.band.hard.binned <- subset(xmm.band.hard.binned, xmm.band.hard.binned$ERROR < as.numeric(quantile(xmm.band.hard.binned$ERROR, 0.95)))
suz.band.all.binned <- subset(suz.band.all.binned, suz.band.all.binned$ERROR < as.numeric(quantile(suz.band.all.binned$ERROR, 0.95)))
suz.band.soft.binned <- subset(suz.band.soft.binned, suz.band.soft.binned$ERROR < as.numeric(quantile(suz.band.soft.binned$ERROR, 0.95)))
suz.band.hard.binned <- subset(suz.band.hard.binned, suz.band.hard.binned$ERROR < as.numeric(quantile(suz.band.hard.binned$ERROR, 0.95)))

########## Plot Error distributions ##########
plot(lc.plot.error(xmm.band.all.binned, bin.width = 0.001))
plot(lc.plot.error(xmm.band.soft.binned, bin.width = 0.001))
plot(lc.plot.error(xmm.band.hard.binned, bin.width = 0.001))
plot(lc.plot.error(xmm.band.uv.binned, bin.width = 0.0005))
plot(lc.plot.error(suz.band.all.binned, bin.width = 0.0005))
plot(lc.plot.error(suz.band.soft.binned, bin.width = 0.0005))
plot(lc.plot.error(suz.band.hard.binned, bin.width = 0.0005))

########## Plot light curves ##########
plot(lc.plot(xmm.band.all.binned, background = TRUE))
plot(lc.plot(xmm.band.soft.binned, background = TRUE))
plot(lc.plot(xmm.band.hard.binned, background = TRUE))
plot(lc.plot(xmm.band.uv.binned, background = TRUE))
plot(lc.plot(suz.band.all.binned, background = TRUE))
plot(lc.plot(suz.band.soft.binned, background = TRUE))
plot(lc.plot(suz.band.hard.binned, background = TRUE))

########## Hardness Ratios ##########
xmm.hr <- hard.ratio(xmm.band.soft.binned, xmm.band.hard.binned)
suz.hr <- hard.ratio(suz.band.soft.binned, suz.band.hard.binned)
plot(hr.plot(xmm.hr))
plot(hr.plot(suz.hr))

########## Flux Flux Plots ##########
xmm.ff <- prep.ff(xmm.band.soft.binned, xmm.band.hard.binned)
suz.ff <- prep.ff(suz.band.soft.binned, suz.band.hard.binned)
plot(ff.plot(xmm.ff))
plot(ff.plot(suz.ff))

########## DCF Contours Helper ##########
dcf.contours <- function(target, nsims, beta, scale.factor, shift.factor) {
  simmed.lc <- sim.lc(beta, bins = length(target$t), length = target$t[length(target$t)], scale.factor =scale.factor, shift.factor = shift.factor)
  colnames(simmed.lc) <- c("t", "dt", "y", "dy")
  sdcf <- cross_correlate(target, simmed.lc, method = "dcf", dtau = xmm.bin.width, use.errors = TRUE)
  sim.results <- data.frame(numeric(0))
  for (i in 2:length(sdcf$tau)) {
    sim.results <- cbind(sim.results, data.frame(numeric(0)))
  }
  colnames(sim.results) <- sdcf$tau
  run <- 0
  while (run <= nsims) {
    simmed.lc <- sim.lc(1.87, bins = length(target$t), length = target$t[length(target$t)], scale.factor = scale.factor, shift.factor = shift.factor)
    colnames(simmed.lc) <- c("t", "dt", "y", "dy")
    sdcf <- cross_correlate(target, simmed.lc, method = "dcf", dtau = xmm.bin.width, use.errors = TRUE)
    cdf <- data.frame(sdcf$ccf[1])
    for (i in 2:length(sdcf$tau)) {
      cdf <- cbind(cdf, data.frame(sdcf$ccf[i]))
    }
    colnames(cdf) <- sdcf$tau
    sim.results <- rbind(sim.results, cdf)
    run <- run + 1
  }
  contours <- data.frame(del=c(0,0,0))
  for (i in 1:length(sdcf$tau)) {
    t <- sim.results[[i]]
    p1 <- quantile(t, 0.90, na.rm = TRUE)
    p2 <- quantile(t, 0.95, na.rm = TRUE)
    p3 <- quantile(t, 0.9999, na.rm = TRUE)
    r <- data.frame(c(p1[[1]],p2[[1]],p3[[1]]))
    colnames(r) <- c(sdcf$tau[i])
    contours <- cbind(contours, r)
  }
  contours <- contours[-1]
  t <- as.data.frame(t(data.frame(sdcf$tau)))
  colnames(t) <- sdcf$tau
  contours <- rbind(contours, t)
  contours <- as.data.frame(t(contours))
  colnames(contours) <- c("p90","p95","p99","tau")
  contours $n90 <- (-1)*contours$p90
  contours $n95 <- (-1)*contours$p95
  contours $n99 <- (-1)*contours$p99
  return(contours)
}

########## Hard Soft DCF ##########
xmm.bu.dcf.lc1 <- data.frame(t = xmm.band.uv$TIME, y = xmm.band.uv$RATE, dy = xmm.band.uv$ERROR)
xmm.bu.dcf.lc2 <- data.frame(t = xmm.band.all.binned$TIME, y = xmm.band.all.binned$RATE, dy = xmm.band.all.binned$ERROR)
xmm.bu.dcf <- cross_correlate(xmm.bu.dcf.lc1, xmm.bu.dcf.lc2, method = "dcf", dtau = xmm.bin.width, use.errors = TRUE)
xmm.bu.contours <- dcf.contours(target = xmm.bu.dcf.lc1, nsims = nsims, beta = beta, scale.factor = 1.5, shift.factor = 2.1)
xmm.bu.dcf.df <- data.frame(tau = xmm.bu.dcf$tau, dtau = xmm.bu.dcf$tau, ccf = xmm.bu.dcf$ccf, lower = xmm.bu.dcf$lower, upper = xmm.bu.dcf$upper)
xmm.bu.dcf.df$dtau <- (xmm.bu.dcf.df$tau[2]-xmm.bu.dcf.df$tau[1])/2
xmm.bu.dcf.plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black") +
  geom_line(data = xmm.bu.dcf.df, aes(x = tau, y = ccf), color = "black") +
  # geom_linerange(data = xmm.bu.dcf.df, aes(x = tau, ymin = lower, ymax = upper), color = "black") +
  # geom_linerange(data = xmm.bu.dcf.df, aes(xmin = tau-dtau, xmax = tau+dtau, y = ccf), color = "black") +
  geom_line(data = xmm.bu.contours, aes(x = tau, y = p99), color = "red") +
  geom_line(data = xmm.bu.contours, aes(x = tau, y = n99), color = "red") +
  geom_line(data = xmm.bu.contours, aes(x = tau, y = p95), color = "blue") +
  geom_line(data = xmm.bu.contours, aes(x = tau, y = n95), color = "blue") +
  geom_line(data = xmm.bu.contours, aes(x = tau, y = p90), color = "green") +
  geom_line(data = xmm.bu.contours, aes(x = tau, y = n90), color = "green") +
  labs(x = "tau", y = "ccf") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.text.x = element_text(margin=margin(6,6,10,6,"pt")),
        axis.text.y = element_text(margin=margin(6,6,10,6,"pt")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(xmm.bu.dcf.plot)

xmm.su.dcf.lc1 <- data.frame(t = xmm.band.uv$TIME, y = xmm.band.uv$RATE, dy = xmm.band.uv$ERROR)
xmm.su.dcf.lc2 <- data.frame(t = xmm.band.soft.binned$TIME, y = xmm.band.soft.binned$RATE, dy = xmm.band.soft.binned$ERROR)
xmm.su.dcf <- cross_correlate(xmm.su.dcf.lc1, xmm.su.dcf.lc2, method = "dcf", dtau = xmm.bin.width)
xmm.su.contours <- dcf.contours(target = xmm.su.dcf.lc1, nsims = nsims, beta = beta, scale.factor = 1.0, shift.factor = 0.8)
xmm.su.dcf.df <- data.frame(tau = xmm.su.dcf$tau, dtau = xmm.su.dcf$tau, ccf = xmm.su.dcf$ccf, lower = xmm.su.dcf$lower, upper = xmm.su.dcf$upper)
xmm.su.dcf.df$dtau <- (xmm.su.dcf.df$tau[2]-xmm.su.dcf.df$tau[1])/2
xmm.su.dcf.plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black") +
  geom_linerange(data = xmm.su.dcf.df, aes(x = tau, ymin = lower, ymax = upper), color = "black") +
  geom_linerange(data = xmm.su.dcf.df, aes(xmin = tau-dtau, xmax = tau+dtau, y = ccf), color = "black") +
  geom_line(data = xmm.su.contours, aes(x = tau, y = p99), color = "red") +
  geom_line(data = xmm.su.contours, aes(x = tau, y = n99), color = "red") +
  geom_line(data = xmm.su.contours, aes(x = tau, y = p95), color = "blue") +
  geom_line(data = xmm.su.contours, aes(x = tau, y = n95), color = "blue") +
  geom_line(data = xmm.su.contours, aes(x = tau, y = p90), color = "green") +
  geom_line(data = xmm.su.contours, aes(x = tau, y = n90), color = "green") +
  labs(x = "tau", y = "ccf") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.text.x = element_text(margin=margin(6,6,10,6,"pt")),
        axis.text.y = element_text(margin=margin(6,6,10,6,"pt")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(xmm.su.dcf.plot)

xmm.hu.dcf.lc1 <- data.frame(t = xmm.band.uv$TIME, y = xmm.band.uv$RATE, dy = xmm.band.uv$ERROR)
xmm.hu.dcf.lc2 <- data.frame(t = xmm.band.hard.binned$TIME, y = xmm.band.hard.binned$RATE, dy = xmm.band.hard.binned$ERROR)
xmm.hu.dcf <- cross_correlate(xmm.hu.dcf.lc1, xmm.hu.dcf.lc2, method = "dcf", dtau = xmm.bin.width)
xmm.hu.contours <- dcf.contours(target = xmm.hu.dcf.lc1, nsims = nsims, beta = beta, scale.factor = 0.5, shift.factor = 0.4)
xmm.hu.dcf.df <- data.frame(tau = xmm.hu.dcf$tau, dtau = xmm.hu.dcf$tau, ccf = xmm.hu.dcf$ccf, lower = xmm.hu.dcf$lower, upper = xmm.hu.dcf$upper)
xmm.hu.dcf.df$dtau <- (xmm.hu.dcf.df$tau[2]-xmm.hu.dcf.df$tau[1])/2
xmm.hu.dcf.plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black") +
  geom_linerange(data = xmm.hu.dcf.df, aes(x = tau, ymin = lower, ymax = upper), color = "black") +
  geom_linerange(data = xmm.hu.dcf.df, aes(xmin = tau-dtau, xmax = tau+dtau, y = ccf), color = "black") +
  geom_line(data = xmm.hu.contours, aes(x = tau, y = p99), color = "red") +
  geom_line(data = xmm.hu.contours, aes(x = tau, y = n99), color = "red") +
  geom_line(data = xmm.hu.contours, aes(x = tau, y = p95), color = "blue") +
  geom_line(data = xmm.hu.contours, aes(x = tau, y = n95), color = "blue") +
  geom_line(data = xmm.hu.contours, aes(x = tau, y = p90), color = "green") +
  geom_line(data = xmm.hu.contours, aes(x = tau, y = n90), color = "green") +
  labs(x = "tau", y = "ccf") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.text.x = element_text(margin=margin(6,6,10,6,"pt")),
        axis.text.y = element_text(margin=margin(6,6,10,6,"pt")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(xmm.hu.dcf.plot)

xmm.hs.dcf.lc1 <- data.frame(t = xmm.band.soft.binned$TIME, y = xmm.band.soft.binned$RATE, dy = xmm.band.soft.binned$ERROR)
xmm.hs.dcf.lc2 <- data.frame(t = xmm.band.hard.binned$TIME, y = xmm.band.hard.binned$RATE, dy = xmm.band.hard.binned$ERROR)
xmm.hs.dcf <- cross_correlate(xmm.hs.dcf.lc1, xmm.hs.dcf.lc2, method = "dcf", dtau = xmm.bin.width)
xmm.hs.contours <- dcf.contours(target = xmm.hs.dcf.lc1, nsims = nsims, beta = beta, scale.factor = 0.5, shift.factor = 0.4)
xmm.hs.dcf.df <- data.frame(tau = xmm.hs.dcf$tau, dtau = xmm.hs.dcf$tau, ccf = xmm.hs.dcf$ccf, lower = xmm.hs.dcf$lower, upper = xmm.hs.dcf$upper)
xmm.hs.dcf.df$dtau <- (xmm.hs.dcf.df$tau[2]-xmm.hs.dcf.df$tau[1])/2
xmm.hs.dcf.plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black") +
  geom_linerange(data = xmm.hs.dcf.df, aes(x = tau, ymin = lower, ymax = upper), color = "black") +
  geom_linerange(data = xmm.hs.dcf.df, aes(xmin = tau-dtau, xmax = tau+dtau, y = ccf), color = "black") +
  geom_line(data = xmm.hs.contours, aes(x = tau, y = p99), color = "red") +
  geom_line(data = xmm.hs.contours, aes(x = tau, y = n99), color = "red") +
  geom_line(data = xmm.hs.contours, aes(x = tau, y = p95), color = "blue") +
  geom_line(data = xmm.hs.contours, aes(x = tau, y = n95), color = "blue") +
  geom_line(data = xmm.hs.contours, aes(x = tau, y = p90), color = "green") +
  geom_line(data = xmm.hs.contours, aes(x = tau, y = n90), color = "green") +
  labs(x = "tau", y = "ccf") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.text.x = element_text(margin=margin(6,6,10,6,"pt")),
        axis.text.y = element_text(margin=margin(6,6,10,6,"pt")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(xmm.hs.dcf.plot)

suz.hs.dcf.lc1 <- data.frame(t = suz.band.soft.binned$TIME, y = suz.band.soft.binned$RATE, dy = suz.band.soft.binned$ERROR)
suz.hs.dcf.lc2 <- data.frame(t = suz.band.hard.binned$TIME, y = suz.band.hard.binned$RATE, dy = suz.band.hard.binned$ERROR)
suz.hs.dcf <- cross_correlate(suz.hs.dcf.lc1, suz.hs.dcf.lc2, method = "dcf", dtau = suz.bin.width)
suz.hs.contours <- dcf.contours(target = suz.hs.dcf.lc1, nsims = nsims, beta = beta, scale.factor = 0.5, shift.factor = 0.4)
suz.hs.dcf.df <- data.frame(tau = suz.hs.dcf$tau, dtau = suz.hs.dcf$tau, ccf = suz.hs.dcf$ccf, lower = suz.hs.dcf$lower, upper = suz.hs.dcf$upper)
suz.hs.dcf.df$dtau <- (suz.hs.dcf.df$tau[2]-suz.hs.dcf.df$tau[1])/2
suz.hs.dcf.plot <- ggplot() +
  geom_hline(yintercept = 0, color = "black") +
  geom_linerange(data = suz.hs.dcf.df, aes(x = tau, ymin = lower, ymax = upper), color = "black") +
  geom_linerange(data = suz.hs.dcf.df, aes(xmin = tau-dtau, xmax = tau+dtau, y = ccf), color = "black") +
  geom_line(data = suz.hs.contours, aes(x = tau, y = p99), color = "red") +
  geom_line(data = suz.hs.contours, aes(x = tau, y = n99), color = "red") +
  geom_line(data = suz.hs.contours, aes(x = tau, y = p95), color = "blue") +
  geom_line(data = suz.hs.contours, aes(x = tau, y = n95), color = "blue") +
  geom_line(data = suz.hs.contours, aes(x = tau, y = p90), color = "green") +
  geom_line(data = suz.hs.contours, aes(x = tau, y = n90), color = "green") +
  labs(x = "tau", y = "ccf") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.text.x = element_text(margin=margin(6,6,10,6,"pt")),
        axis.text.y = element_text(margin=margin(6,6,10,6,"pt")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(suz.hs.dcf.plot)
