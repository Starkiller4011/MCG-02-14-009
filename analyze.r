#!/usr/bin/env Rscript

########## Setup ##########
options(scipen=999)

########## Load required libraries ##########
require(FITSio)
require(sour)
require(ggplot2)
require(AstroR)

########## GLOBAL CONSTANTS ##########
xmm.bin.width <- 1000
suzaku.bin.width <- 5760
image.width <- 7.24
image.height <- 3.46
nsims <- 10000
beta <- 1.87
dcf.bin.width <- xmm.bin.width

########## GLOBAL PATH CONSTANTS ##########
cwd <- normalizePath(file.path('.'), winslash = "/")
lc.path <- normalizePath(file.path(cwd, 'data', 'lightcurves'), winslash = "/")
pn.lc.path <- normalizePath(file.path(lc.path, 'pn', 'lc'), winslash = "/")
uvw1.lc.path <- normalizePath(file.path(lc.path, 'uvw1'), winslash = "/")
xis.lc.path <- normalizePath(file.path(lc.path, 'xis', 'new', 'lc'), winslash = "/")
pn.fvar.path <- normalizePath(file.path(lc.path, 'pn', 'fvar'), winslash = "/")
xis.fvar.path <- normalizePath(file.path(lc.path, 'xis', 'new', 'fvar'), winslash = "/")
plots.path <- normalizePath(file.path(cwd, "plots"), winslash = "/")
pn.plots.path <- normalizePath(file.path(plots.path, "PN_Spectrum"), winslash = "/")
xis.plots.path <- normalizePath(file.path(plots.path, "XIS_Spectrum"), winslash = "/")
rgs.plots.path <- normalizePath(file.path(plots.path, "RGS_Spectrum"), winslash = "/")
timing.plots.path <- normalizePath(file.path(plots.path, "Timing"), winslash = "/")

########## Set File Paths ##########
xmm.broadband.path <- normalizePath(file.path(pn.lc.path, '0550640101pn_lccor_300-10000.fits'), winslash = "/")
xmm.soft.path <- normalizePath(file.path(pn.lc.path, '0550640101pn_lccor_300-2000.fits'), winslash = "/")
xmm.hard.path <- normalizePath(file.path(pn.lc.path, '0550640101pn_lccor_2000-10000.fits'), winslash = "/")
xmm.uvw1.path <- normalizePath(file.path(uvw1.lc.path, "0550640101_uvw1.lc"), winslash = "/")
suzaku.broadband.path <- c(normalizePath(file.path(xis.lc.path, '05_10_src.lc'), winslash = "/"), normalizePath(file.path(xis.lc.path, '05_10_bkg.lc'), winslash = "/"))
suzaku.soft.path <- c(normalizePath(file.path(xis.lc.path, '05_2_src.lc'), winslash = "/"), normalizePath(file.path(xis.lc.path, '05_2_bkg.lc'), winslash = "/"))
suzaku.hard.path <- c(normalizePath(file.path(xis.lc.path, '2_10_src.lc'), winslash = "/"), normalizePath(file.path(xis.lc.path, '2_10_bkg.lc'), winslash = "/"))

########## Load Data ##########
xmm.broadband <- xmm.pn.lc(xmm.broadband.path)
xmm.soft <- xmm.pn.lc(xmm.soft.path)
xmm.hard <- xmm.pn.lc(xmm.hard.path)
xmm.uvw1 <- read.csv(xmm.uvw1.path)
colnames(xmm.uvw1) <- c("TIME", "RATE", "ERROR")
xmm.uvw1$TIMED <- 0
for (i in 2:length(xmm.uvw1$TIME)) {
  xmm.uvw1$TIMED[i] <- (xmm.uvw1$TIME[i] - xmm.uvw1$TIME[i-1])/2
}
xmm.uvw1$BACKV <- 0.1
xmm.uvw1$BACKE <- 0.1
suzaku.broadband <- suzaku.lc(suzaku.broadband.path[1], suzaku.broadband.path[2])
suzaku.soft <- suzaku.lc(suzaku.soft.path[1], suzaku.soft.path[2])
suzaku.hard <- suzaku.lc(suzaku.hard.path[1], suzaku.hard.path[2])

########## Set Origin ##########
xmm.broadband <- lc.setOrigin(xmm.broadband, 0)
xmm.soft <- lc.setOrigin(xmm.soft, 0)
xmm.hard <- lc.setOrigin(xmm.hard, 0)
suzaku.broadband <- lc.setOrigin(suzaku.broadband, 0)
suzaku.soft <- lc.setOrigin(suzaku.soft, 0)
suzaku.hard <- lc.setOrigin(suzaku.hard, 0)

########## Remove Background Flaring ##########
xmm.broadband <- subset(xmm.broadband, xmm.broadband$TIME < 90000)
xmm.soft <- subset(xmm.soft, xmm.soft$TIME < 90000)
xmm.hard <- subset(xmm.hard, xmm.hard$TIME < 90000)

########## Bin Data ##########
xmm.broadband <- bin.lc(xmm.broadband, xmm.bin.width)
xmm.soft <- bin.lc(xmm.soft, xmm.bin.width)
xmm.hard <- bin.lc(xmm.hard, xmm.bin.width)
suzaku.broadband <- bin.lc(suzaku.broadband, suzaku.bin.width)
suzaku.soft <- bin.lc(suzaku.soft, suzaku.bin.width)
suzaku.hard <- bin.lc(suzaku.hard, suzaku.bin.width)

########## Get Hardness Data ##########
xmm.hrdf <- hard.ratio(xmm.soft, xmm.hard)
suzaku.hrdf <- hard.ratio(suzaku.soft, suzaku.hard)

########## Get Flux-Flux Data ##########
xmm.ffdf <- prep.ff(xmm.soft, xmm.hard)
suzaku.ffdf <- prep.ff(suzaku.soft, suzaku.hard)

########## Create Plots ##########
## XMM
plot(lc.plot(xmm.uvw1, background = TRUE, plt.title = "XMM UVW1 Light Curve"))
ggsave(normalizePath(file.path(pn.plots.path, "xmm_uvw1_bg.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)
plot(lc.plot(xmm.broadband, background = TRUE, plt.title = "XMM Broadband Light Curve"))
ggsave(normalizePath(file.path(pn.plots.path, "xmm_broadband.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)
plot(hs.plot(xmm.soft, xmm.hard, plt.title = "XMM Hard and Soft Light Curves"))
ggsave(normalizePath(file.path(pn.plots.path, "xmm_hs.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)
plot(hr.plot(xmm.hrdf, plt.title = "XMM Hardness Ratio"))
ggsave(normalizePath(file.path(pn.plots.path, "xmm_hr.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)
plot(ff.plot(xmm.ffdf, plt.title = "XMM Flux-Flux Plot"))
ggsave(normalizePath(file.path(pn.plots.path, "xmm_ff.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)

## Suzaku
plot(lc.plot(suzaku.broadband, background = TRUE, plt.title = "Suzaku Broadband Light Curve"))
ggsave(normalizePath(file.path(xis.plots.path, "suzaku_broadband.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)
plot(hs.plot(suzaku.soft, suzaku.hard, plt.title = "Suzaku Hard and Soft Light Curves"))
ggsave(normalizePath(file.path(xis.plots.path, "suzaku_hs.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)
plot(hr.plot(suzaku.hrdf, plt.title = "Suzaku Hardness Ratio"))
ggsave(normalizePath(file.path(xis.plots.path, "suzaku_hr.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)
plot(ff.plot(suzaku.ffdf, plt.title = "Suzaku Flux-Flux Plot"))
ggsave(normalizePath(file.path(xis.plots.path, "suzaku_ff.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)

## XMM UVW1
plot(lc.plot(xmm.uvw1, plt.title = "XMM UVW1 Light Curve", background = FALSE))
ggsave(normalizePath(file.path(pn.plots.path, "xmm_uvw1.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)

## XMM UVW1 vs Broadband
# Normalize Light Curves
xmm.uvw1.norm <- normalize.lc(xmm.uvw1)
xmm.broadband.norm <- normalize.lc(xmm.broadband)
# Scale UVW1 - ** Interpretation is a dampening of variability when the X-Rays impact the disc, possibly due to a cooler corona? 
xmm.uvw1.norm$RATE <- (xmm.uvw1.norm$RATE - 1)
xmm.broadband.norm$RATE <- (xmm.broadband.norm$RATE - 1)
rate.scale.factor <- (max(xmm.broadband.norm$RATE)/max(xmm.uvw1.norm$RATE))
error.scale.factor <- (max(xmm.broadband.norm$ERROR)/max(xmm.uvw1.norm$ERROR))
xmm.uvw1.norm$RATE <- (xmm.uvw1.norm$RATE * rate.scale.factor)
xmm.uvw1.norm$ERROR <- (xmm.uvw1.norm$ERROR * error.scale.factor)
xmm.uvw1.norm$RATE <- (xmm.uvw1.norm$RATE + 1)
xmm.broadband.norm$RATE <- (xmm.broadband.norm$RATE + 1)

########## DCF ###########

#### Helper Functions ####
dcf.contours <- function(target, nsims, bin.width, beta, scale.factor, shift.factor) {
  simmed.lc <- sim.lc(beta, bins = length(target$t), length = target$t[length(target$t)], scale.factor =scale.factor, shift.factor = shift.factor)
  colnames(simmed.lc) <- c("t", "dt", "y", "dy")
  sdcf <- cross_correlate(target, simmed.lc, method = "dcf", dtau = bin.width, use.errors = TRUE)
  sim.results <- data.frame(numeric(0))
  for (i in 2:length(sdcf$tau)) {
    sim.results <- cbind(sim.results, data.frame(numeric(0)))
  }
  colnames(sim.results) <- sdcf$tau
  run <- 0
  while (run <= nsims) {
    simmed.lc <- sim.lc(1.87, bins = length(target$t), length = target$t[length(target$t)], scale.factor = scale.factor, shift.factor = shift.factor)
    colnames(simmed.lc) <- c("t", "dt", "y", "dy")
    sdcf <- cross_correlate(target, simmed.lc, method = "dcf", dtau = bin.width, use.errors = TRUE)
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
    p1 <- quantile(t, 0.90)
    p2 <- quantile(t, 0.95)
    p3 <- quantile(t, 0.9999)
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

## UVW1 vs Broadband DCF ##
dcflc.xmm.uvw1 <- data.frame(t = xmm.uvw1$TIME, y = xmm.uvw1$RATE, dy = xmm.uvw1$ERROR)
dcflc.xmm.bb <- data.frame(t = xmm.broadband$TIME, y = xmm.broadband$RATE, dy = xmm.broadband$ERROR)

dcf.uvw1.bb <- cross_correlate(dcflc.xmm.uvw1, dcflc.xmm.bb, method = "dcf", dtau = dcf.bin.width, use.errors = TRUE)
dcf.uvw1.bb.contours <- dcf.contours(target = dcflc.xmm.uvw1, nsims = nsims, bin.width = dcf.bin.width, beta = beta, scale.factor = (max(xmm.broadband.norm$RATE)/mean(xmm.broadband$RATE)), shift.factor = mean(xmm.broadband$RATE))

dcf.uvw1.bb.df <- data.frame(tau = dcf.uvw1.bb$tau, dtau = dcf.uvw1.bb$tau, dcf = dcf.uvw1.bb$ccf)
dcf.uvw1.bb.df$dtau <- (dcf.uvw1.bb.df$tau[2]-dcf.uvw1.bb.df$tau[1])/2

contour.99 <- data.frame(tau = dcf.uvw1.bb.contours$tau, upper = dcf.uvw1.bb.contours$p99, lower = dcf.uvw1.bb.contours$n99)
contour.95 <- data.frame(tau = dcf.uvw1.bb.contours$tau, upper = dcf.uvw1.bb.contours$p95, lower = dcf.uvw1.bb.contours$n95)
contour.90 <- data.frame(tau = dcf.uvw1.bb.contours$tau, upper = dcf.uvw1.bb.contours$p90, lower = dcf.uvw1.bb.contours$n90)

colors <- c("DCF" = "black","99.99%" = "blue", "95%" = "orange", "90%" = "red")
dcf.uvw1.bb.plot <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_line(data = dcf.uvw1.bb.df, aes(x = tau, y = dcf), size = 1, color = "black") +
  geom_line(data = contour.99, aes(x = tau, y = upper, color = "99.99%"), linetype = "dashed", size = 1) +
  geom_line(data = contour.99, aes(x = tau, y = lower, color = "99.99%"), linetype = "dashed", size = 1) +
  geom_line(data = contour.95, aes(x = tau, y = upper, color = "95%"), linetype = "dotdash", size = 1) +
  geom_line(data = contour.95, aes(x = tau, y = lower, color = "95%"), linetype = "dotdash", size = 1) +
  geom_line(data = contour.90, aes(x = tau, y = upper, color = "90%"), linetype = "dotted", size = 1) +
  geom_line(data = contour.90, aes(x = tau, y = lower, color = "90%"), linetype = "dotted", size = 1) +
  scale_x_continuous(breaks = round(seq(-30000, 30000, by = 10000), 1)) +
  scale_y_continuous(breaks = round(seq(-1, 1, by = 0.25), 1)) + ylim(-1, 1) +
  scale_color_manual(name = "Confidence", values = colors) +
  labs(x = "tau [s]", y = "dcf", title = "Broadband UV DCF") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.ticks.length = unit(-3, "pt"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(margin=margin(6,6,10,6,"pt"), size=14),
        axis.text.y = element_text(margin=margin(6,6,10,6,"pt"), size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(dcf.uvw1.bb.plot)
ggsave(normalizePath(file.path(pn.plots.path, "xmm_uvw1_dcf.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)


## Plot Normalized and scaled UVW1 vs Broadband Unshifted

plot(lc.overplot(xmm.broadband.norm, xmm.uvw1.norm, rate.scale.factor, plt.title = "Normalized Unshifted UVW1 and PN Broadband", y1lab = "[normalized count/s]", y2lab = "[normalized count/s]", names = c("Broadband", "UVW1")))

plot(lc.plot(xmm.uvw1.norm, plt.title = "XMM UVW1 Light Curve", background = FALSE))
plot(lc.overplot(xmm.broadband.norm, xmm.uvw1.norm, plt.title = "Normalized Unshifted UVW1 and PN Broadband", names = c("Broadband", "UVW1")))
ggsave(normalizePath(file.path(pn.plots.path, "xmm_uvw1_normed_unshifted.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)

## Apply approximate shift
SHIFT <- -12000
xmm.uvw1.norm.shifted <- xmm.uvw1.norm
xmm.uvw1.norm.shifted$TIME <- (xmm.uvw1.norm.shifted$TIME + SHIFT)

plot(lc.overplot(xmm.broadband.norm, xmm.uvw1.norm.shifted, plt.title = "Normalized Shifted UVW1 and PN Broadband", names = c("Broadband", "UVW1")))
ggsave(normalizePath(file.path(pn.plots.path, "xmm_uvw1_normed_shifted.png"), winslash = "/", mustWork = FALSE), width = image.width, height = image.height)






