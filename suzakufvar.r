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
suz.bin.width = 5760
fvar.path <- "data/lightcurves/xis/new/fvar/"

bb.src.path <- "data/lightcurves/xis/lc/2_10_src.lc"
bb.bak.path <- "data/lightcurves/xis/lc/2_10_bkg.lc"

lc <- suzaku.lc(bb.src.path, bb.bak.path)
lc <- bin.lc(lc, suz.bin.width)
plot(lc.plot(lc, background = TRUE))


low <- list(400,1000,1500,2000,4000)
high <- list(1000,1500,2000,4000,10000)

e <- list()
e.err <- list()
fvar <- list()
fvar.err <- list()

for (i in 1:length(low)) {
  l <- (low[[i]]/1000)
  h <- (high[[i]]/1000)
  diff <- ((h-l)/2)
  mid <- l+diff
  e[i] <- mid
  e.err[i] <- diff
  src <- paste(fvar.path,paste(paste(toString(low[i]),toString(high[i]),sep="_"),"_src.lc",sep=""),sep="")
  bkg <- paste(fvar.path,paste(paste(toString(low[i]),toString(high[i]),sep="_"),"_bkg.lc",sep=""),sep="")
  lc <- suzaku.lc(src, bkg)
  lc <- bin.lc(lc,suz.bin.width)
  fv <- fvar.ponti(lc)
  fvar[i] <- fv[[1]]
  fvar.err <- fv[[2]]
  plot(lc.plot(lc,background = FALSE))
}

e <- unlist(e)
e.err <- unlist(e.err)
fvar <- unlist(fvar)
fvar.err <- unlist(fvar.err)

suzaku.fvar <- data.frame(e, e.err, fvar, fvar.err)
colnames(suzaku.fvar) <- c("e", "e.err", "fvar", "fvar.err")

plt <- ggplot(data=suzaku.fvar) +
  geom_linerange(aes(x = e, ymin = fvar-fvar.err, ymax = fvar+fvar.err), color = "black") +
  geom_linerange(aes(xmin = e-e.err, xmax = e+e.err, y=fvar), color = "black") +
  geom_hline(yintercept = 0, color = "green") +
  labs(x = "Energy (keV)", y = "FVAR", title = "Suzaku Fractional Variability") +
  scale_x_log10(breaks = c(0.5,1.0,2.0,5.0,10.0)) +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.ticks.length = unit(-3, "pt"),
        axis.text.x = element_text(margin = margin(6,6,10,6,"pt"), size = 14),
        axis.text.y = element_text(margin = margin(6,6,10,6,"pt"), size = 14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot(plt)
ggsave("plots/XIS_Spectrum/suzaku_fvar.png")




