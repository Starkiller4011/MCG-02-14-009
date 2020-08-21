#!/usr/bin/env Rscript

########## Setup ##########
options(scipen=999)

########## Load required libraries ##########
require(FITSio)
require(ggplot2)
require(AstroR)

########## Set defaults ##########
cwd <- normalizePath('.')
lc_path <- normalizePath(file.path(cwd, 'data', 'lightcurves'))
pn_lc_path <- normalizePath(file.path(lc_path, 'pn/new'))
uvw1_lc_path <- normalizePath(file.path(lc_path, 'uvw1'))
broadband <- list("300-10000")
hr_bins <- list("300-1000", "2000-10000")
lag_f_bins <- list("300-1000", "1000-4000")
#lag_e_bins <- list("300-425", "425-600", "600-850", "850-1200", "1200-1750", "1750-2500", "2500-3500", "3500-5000", "5000-7000", "7000-10000")
lag_e_bins <- list("300-358", "358-426", "426-508", "508-605", "605-721", "721-859", "859-1024", "1024-1220", "1220-1454", "1454-1733", "1733-2064", "2064-2460", "2460-2931", "2931-3493", "3493-4162", "4162-4960", "4960-5910", "5910-7043", "7043-8392", "8392-10000")

########## User set variables ##########
obsid <- "0550640101"
bg.cap <- TRUE
lc.time.cap = 90000
lc.bgr.cap = 0.1

########## Helper Functions ##########
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

get_lc_data <- function(e_bins) {
  lc.data <- list()
  i <- 1
  for (bin in e_bins) {
    lc_filename <- paste(list(obsid,"pn_lccor_",bin,".fits"), collapse = "")
    lc_path <- normalizePath(file.path(pn_lc_path,lc_filename))
    lc <- readFITS(lc_path)
    time <- lc$col[[1]]-lc$col[[1]][[1]]
    sc.rate <- lc$col[[2]]
    sc.error <- lc$col[[3]]
    bg.rate <- lc$col[[5]]
    bg.error <- lc$col[[6]]
    lc.data[[i]] <- data.frame(
      time,
      sc.rate,
      sc.error,
      bg.rate,
      bg.error,
      stringsAsFactors = FALSE
    )
    i = i+1
  }
  return(lc.data)
}

bin_lc <- function(lc, bin.width) {
  nbins <- floor((lc$time[[length(lc$time)]]-lc$time[[1]])/bin.width)
  time <- seq(from = (bin.width/2), to = ((bin.width*nbins)-(bin.width/2)), by = bin.width)
  sc.rate <- list()
  sc.error <- list()
  bg.rate <- list()
  bg.error <- list()
  i <- 1
  for (bin in time) {
    sub <- lc[which(lc$time >= bin-(bin.width/2)),]
    sub <- sub[which(sub$time < bin+(bin.width/2)),]
    sc.rate[[i]] <- mean(sub$sc.rate)
    serr <- list()
    j <- 1
    for (err in sub$sc.error) {
      serr[[j]] <- err^2
      j <- j+1
    }
    if (length(serr) > 0) {
      sc.error[[i]] <- sqrt(Reduce("+",serr))/length(serr)
    } else {
      sc.error[[i]] <- NaN
    }
    bg.rate[[i]] <- mean(sub$bg.rate)
    berr <- list()
    j <- 1
    for (err in sub$bg.error) {
      berr[[j]] <- err^2
      j <- j+1
    }
    if (length(berr) > 0) {
      bg.error[[i]] <- sqrt(Reduce("+",berr))/length(berr)
    } else {
      bg.error[[i]] <- NaN
    }
    i <- i+1
  }
  sc.rate <- unlist(sc.rate, use.names = FALSE)
  sc.error <- unlist(sc.error, use.names = FALSE)
  bg.rate <- unlist(bg.rate, use.names = FALSE)
  bg.error <- unlist(bg.error, use.names = FALSE)
  lc.new <- data.frame(
    time,
    sc.rate,
    sc.error,
    bg.rate,
    bg.error,
    stringsAsFactors = FALSE
  )
  return(lc.new)
}

lc_fvar <- function(lc) {
  rate <- lc$sc.rate
  error <- lc$sc.error
  X <- mean(rate)
  s2 <- var(rate)
  s2err <- mean(error^2)
  fvar <- sqrt((s2-s2err)/(X^2))
  fvar.err <- ((1/fvar)*sqrt(1/(2*length(rate)))*(s2/(X^2)))
  return(c(fvar,fvar.err))
}

power.trend <- function(repetitions = 5, x = c(0, 1), sd = 1, slope = 1, alpha = 0.05) {
  X <- rep(x, repetitions)
  ncp <- slope ^ 2 * sum((X - mean(X))^2) / sd ^ 2
  return(1 - pf(qf(1 - alpha, 1, length(X) - 2), 1, length(X) - 2, ncp = ncp))
}

########## Load Data ##########
broadband.data <- get_lc_data(broadband)
hardness.ratio.data <- get_lc_data(hr_bins)
lag.frequency.data <- get_lc_data(lag_f_bins)
lag.energy.data <- get_lc_data(lag_e_bins)

########## Broadband ##########
b.bin_size = 100
b.init.data <- broadband.data[[1]]
b.data <- bin_lc(b.init.data, b.bin_size)
if (bg.cap) {
  b.data <- b.data[which(b.data$bg.rate < lc.bgr.cap),]
}
b.data <- b.data[which(b.data$time < lc.time.cap),]
b.data <- b.data[complete.cases(b.data),]
bplt <- ggplot() +
  geom_linerange(data = b.data, aes(x = time, ymin = sc.rate-sc.error, ymax = sc.rate+sc.error), color = "black") +
  geom_linerange(data = b.data, aes(x = time, ymin = bg.rate-bg.error, ymax = bg.rate+bg.error), color = "grey") +
  #labs(x = "Time [s]", y = "Rate [count/s]", title = "0.3-10 keV Light Curve") +
  labs(x = "Time [s]", y = "Rate [count/s]") +
  theme_bw() + theme(plot.title = element_text(size = 18, hjust = 0.5),
                     plot.subtitle = element_text(size = 12, hjust = 0.5),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
plot(bplt)

########## Hardness ratio ##########
hr.bin_size = 100
s.data <- hardness.ratio.data[[1]]
h.data <- hardness.ratio.data[[2]]
s.data <- bin_lc(s.data, hr.bin_size)
h.data <- bin_lc(h.data, hr.bin_size)
if (bg.cap) {
  s.data <- s.data[which(s.data$bg.rate < lc.bgr.cap),]
  h.data <- h.data[which(h.data$bg.rate < lc.bgr.cap),]
}
s.data <- s.data[which(s.data$time %in% b.data$time),]
h.data <- h.data[which(h.data$time %in% b.data$time),]
s.data <- s.data[complete.cases(s.data),]
h.data <- h.data[complete.cases(h.data),]
hr.data <- data.frame(time = s.data$time, s.rate = s.data$sc.rate, s.error = s.data$sc.error, h.rate = h.data$sc.rate, h.error = h.data$sc.error, stringsAsFactors = FALSE)
hr.data$hr <- (hr.data$h.rate-hr.data$s.rate)/(hr.data$h.rate+hr.data$s.rate)
hr.data$hre <- sqrt(((((2*hr.data$s.rate)/(hr.data$h.rate+hr.data$s.rate)^2)*hr.data$h.error)^2)+((((2*hr.data$h.rate)/(hr.data$h.rate+hr.data$s.rate)^2)*hr.data$h.error)^2))
lm.fit <- lm(formula = hr~time, data = hr.data)
pval <- formatC(anova(lm.fit)$'Pr(>F)'[1], format = "e", digits = 2)
a <- as.character(format(unname(coef(lm.fit)[2]), digits = 2))
b <- as.numeric(format(unname(coef(lm.fit)[1]), digits = 2))
r <- as.character(format(summary(lm.fit)$adj.r.squared, digits = 3))
if (b < 0) {
  s <- "-"
} else {
  s <- "+"
}
b <- as.character(abs(b))
lm.eq <- paste(list("y=",a,"*x",s,b,", r-sqr=",r,", p-value=",pval), sep = "", collapse = "")
hr.predicted <- data.frame(hr.pred = predict(lm.fit, hr.data), hr.time = hr.data$time)
hsplt <- ggplot() +
  geom_linerange(data = s.data, aes(x = time, ymin = sc.rate-sc.error, ymax = sc.rate+sc.error), color = "black") +
  geom_linerange(data = s.data, aes(x = time, ymin = bg.rate-bg.error, ymax = bg.rate+bg.error), color = "red") +
  geom_linerange(data = h.data, aes(x = time, ymin = sc.rate-sc.error, ymax = sc.rate+sc.error), color = "blue") +
  geom_linerange(data = h.data, aes(x = time, ymin = bg.rate-bg.error, ymax = bg.rate+bg.error), color = "green") +
  #labs(x = "Time [s]", y = "Rate [count/s]", title = "Hard(2-10 keV) and Soft(0.3-1 keV) Light Curves") +
  labs(x = "Time [s]", y = "Rate [count/s]") +
  theme_bw() + theme(plot.title = element_text(size = 18, hjust = 0.5),
                     plot.subtitle = element_text(size = 12, hjust = 0.5),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
plot(hsplt)
hrplt <- ggplot() +
  # geom_point(data = hr.data, aes(x = time, y = hr), size = 1, color = "black") +
  geom_linerange(data = hr.data, aes(x = time, ymin = hr-hre, ymax = hr+hre), color = "black") +
  geom_line(data = hr.predicted, aes(x = hr.time, y = hr.pred), size = 1.2, color = "red") +
  geom_segment(aes(x = hr.data$time, xend = hr.data$time[[length(hr.data$time)]], y = mean(hr.data$hr), yend = mean(hr.data$hr)), size = 1.2, color = "blue") +
  coord_cartesian(ylim = c(-1, 0)) +
  #labs(x = "Time [s]", y = "Ratio", title = "Hardness Ratio [0.3-1 keV vs 2-10 keV]", subtitle = lm.eq) +
  labs(x = "Time [s]", y = "Hardness Ratio") +
  theme_bw() + theme(plot.title = element_text(size = 18, hjust = 0.5),
                     plot.subtitle = element_text(size = 12, hjust = 0.5),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
plot(hrplt)
# print(power.trend(repetitions = 1, x = hr.data$time, sd = sd(hr.data$hr), slope = coef(lm.fit)[2]))
# print(power.trend(repetitions = 1, x = hr.data$time, sd = sd(hr.data$hr), slope = 0.00000035))
# print(power.trend(repetitions = 1, x = hr.data$time, sd = sd(hr.data$hr), slope = 0.0000008))

########## Flux Flux Plot ##########
ff.bin_size = 500
ff.b.data <- broadband.data[[1]]
ff.b.data <- bin_lc(ff.b.data, ff.bin_size)
if (bg.cap) {
  ff.b.data <- ff.b.data[which(ff.b.data$bg.rate < lc.bgr.cap),]
}
ff.b.data <- ff.b.data[which(ff.b.data$time < lc.time.cap),]
ff.b.data <- ff.b.data[complete.cases(ff.b.data),]
ff.s.data <- hardness.ratio.data[[1]]
ff.h.data <- hardness.ratio.data[[2]]
ff.s.data <- bin_lc(ff.s.data, ff.bin_size)
ff.h.data <- bin_lc(ff.h.data, ff.bin_size)
if (bg.cap) {
  ff.s.data <- ff.s.data[which(ff.s.data$bg.rate < lc.bgr.cap),]
  ff.h.data <- ff.h.data[which(ff.h.data$bg.rate < lc.bgr.cap),] 
}
ff.s.data <- ff.s.data[which(ff.s.data$time %in% ff.b.data$time),]
ff.h.data <- ff.h.data[which(ff.h.data$time %in% ff.b.data$time),]
ff.s.data <- ff.s.data[complete.cases(ff.s.data),]
ff.h.data <- ff.h.data[complete.cases(ff.h.data),]
ff.data <- data.frame(time = ff.s.data$time, sr=ff.s.data$sc.rate, se=ff.s.data$sc.error, hr=ff.h.data$sc.rate, he=ff.h.data$sc.error, stringsAsFactors = FALSE)
lin.fit <- lm(formula = hr~sr, data = ff.data)
lpval <- formatC(anova(lin.fit)$'Pr(>F)'[1], format = "e", digits = 2)
la <- as.character(format(unname(coef(lin.fit)[2]), digits = 2))
lb <- as.numeric(format(unname(coef(lin.fit)[1]), digits = 2))
lr <- as.character(format(summary(lin.fit)$adj.r.squared, digits = 3))
if (lb < 0) {
  ls <- "-"
} else {
  ls <- "+"
}
lb <- as.character(abs(lb))
pl.fit <- lm(formula = log10(hr)~log10(sr), data = ff.data)
plpval <- formatC(anova(pl.fit)$'Pr(>F)'[1], format = "e", digits = 2)
pla <- as.character(format(unname(coef(pl.fit)[2]), digits = 2))
plb <- as.numeric(format(unname(coef(pl.fit)[1]), digits = 2))
plr <- as.character(format(summary(pl.fit)$adj.r.squared, digits = 3))
if (plb < 0) {
  pls <- "-"
} else {
  pls <- "+"
}
plb <- as.character(abs(plb))
lm.eq <- paste(list("linear: y=",la,"*x",ls,lb,", r-sqr=",lr,", p-value=",lpval,
                    "\npower law: log(y)=",pla,"*log(x)",pls,plb,", r-sqr=",plr,", p-value=",plpval),
               sep = "", collapse = "")
lin.pred <- data.frame(sr = ff.data$sr, hr = predict(lin.fit, ff.data))
pl.pred <- data.frame(sr = ff.data$sr, hr = 10^predict(pl.fit, ff.data))
ffplt <- ggplot() +
  geom_linerange(data = ff.data, aes(xmin = sr-se, xmax = sr+se, y = hr)) +
  geom_linerange(data = ff.data, aes(x = sr, ymin = hr-he, ymax = hr+he)) +
  geom_point(data = ff.data, aes(x = sr, y = hr), size = 1) +
  geom_line(data = lin.pred, aes(x = sr, y = hr), size = 1.2, color = "red") +
  geom_line(data = pl.pred, aes(x = sr, y = hr), size = 1.2, color = "blue") +
  scale_y_log10() + scale_x_log10() +
  #labs(x = "0.3-1 keV Rate [count/s]", y = "2-10 keV Rate [count/s]", title = "Flux Flux [0.3-1 keV vs 2-10 keV]", subtitle = lm.eq) +
  labs(x = "0.3-1 keV Rate [count/s]", y = "2-10 keV Rate [count/s]") +
  theme_bw() + theme(plot.title = element_text(size = 18, hjust = 0.5),
                     plot.subtitle = element_text(size = 12, hjust = 0.5),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
plot(ffplt)

########## Fvar ##########
fv.bin_size = 1000
fvars <- list()
fvar.errors <- list()
i <- 1
for (lc in lag.energy.data) {
  lc.data <- bin_lc(lc, fv.bin_size)
  if (bg.cap) {
    lc.data <- lc.data[which(lc.data$bg.rate < lc.bgr.cap),]
  }
  lc.data <- lc.data[which(lc.data$time < lc.time.cap),]
  lc.data <- lc.data[complete.cases(lc.data),]
  lc.fvar <- lc_fvar(lc.data)
  fvars[[i]] <- lc.fvar[[1]]
  fvar.errors[[i]] <- lc.fvar[[2]]
  i <- i + 1
}
fvars <- unlist(fvars, use.names = FALSE)
fvar.errors <- unlist(fvar.errors, use.names = FALSE)
es <- c(0.300, 0.358, 0.426, 0.508, 0.605, 0.721, 0.859, 1.024, 1.220, 1.454, 1.733, 2.064, 2.460, 2.931, 3.493, 4.162, 4.960, 5.910, 7.403, 8.392)
ee <- c(0.358, 0.426, 0.508, 0.605, 0.721, 0.859, 1.024, 1.220, 1.454, 1.733, 2.064, 2.460, 2.931, 3.493, 4.162, 4.960, 5.910, 7.403, 8.392, 10.000)

fvar.data <- data.frame(es, ee, fv = fvars, fve = fvar.errors, stringsAsFactors = FALSE)
fvar.data$em <- fvar.data$es + ((fvar.data$ee-fvar.data$es)/2)

fv.lin.fit <- lm(formula = fv~em, data = fvar.data)
fv.lpval <- formatC(anova(fv.lin.fit)$'Pr(>F)'[1], format = "e", digits = 2)
fv.la <- as.character(format(unname(coef(fv.lin.fit)[2]), digits = 2))
fv.lb <- as.numeric(format(unname(coef(fv.lin.fit)[1]), digits = 2))
fv.lr <- as.character(format(summary(fv.lin.fit)$adj.r.squared, digits = 3))
if (fv.lb < 0) {
  fv.ls <- "-"
} else {
  fv.ls <- "+"
}
fv.lb <- as.character(abs(fv.lb))
fv.pl.fit <- lm(formula = log10(fv)~log10(em), data = fvar.data)
fv.plpval <- formatC(anova(fv.pl.fit)$'Pr(>F)'[1], format = "e", digits = 2)
fv.pla <- as.character(format(unname(coef(fv.pl.fit)[2]), digits = 2))
fv.plb <- as.numeric(format(unname(coef(fv.pl.fit)[1]), digits = 2))
fv.plr <- as.character(format(summary(fv.pl.fit)$adj.r.squared, digits = 3))
if (fv.plb < 0) {
  fv.pls <- "-"
} else {
  fv.pls <- "+"
}
fv.plb <- as.character(abs(fv.plb))
fv.lm.eq <- paste(list("linear: y=",fv.la,"*x",fv.ls,fv.lb,", r-sqr=",fv.lr,", p-value=",fv.lpval,
                    "\npower law: log(y)=",fv.pla,"*log(x)",fv.pls,fv.plb,", r-sqr=",fv.plr,", p-value=",fv.plpval),
               sep = "", collapse = "")
fv.lin.pred <- data.frame(em = fvar.data$em, fv = predict(fv.lin.fit, fvar.data))
fv.pl.pred <- data.frame(em = fvar.data$em, fv = 10^predict(fv.pl.fit, fvar.data))

fvplt <- ggplot() +
  geom_linerange(data = fvar.data, aes(x = em, ymin = fv-fve, ymax = fv+fve)) +
  geom_linerange(data = fvar.data, aes(xmin = em-((ee-es)/2), xmax = em+((ee-es)/2), y = fv)) +
  geom_point(data = fvar.data, aes(x = em, y = fv), size = 1) +
  #geom_line(data = fv.lin.pred, aes(x = em, y = fv, group = 1), size = 1.2, color = "red") +
  geom_line(data = fv.pl.pred, aes(x = em, y = fv, group = 1), size = 1.2, color = "red") +
  scale_y_log10() +
  scale_x_log10(breaks = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0),
                labels = c("","","0.5","","","","","1.0","2.0","","","5.0","","","","","10.0")) +
  #scale_x_log10(breaks = c(0.3,0.3625,0.425,0.5125,0.6,0.725,0.85,1.025,1.2,1.475,1.75,2.125,2.5,3,3.5,4.25,5,6,7,8.5,10)) +
  #labs(x = "Energy (keV)", y = "Fractional Variability", title = "Fractional Variability", subtitle = fv.lm.eq) +
  labs(x = "Energy (keV)", y = "Fractional Variability") +
  theme_bw() + theme(plot.title = element_text(size = 18, hjust = 0.5),
                     #plot.subtitle = element_text(size = 12, hjust = 0.5),
                     axis.ticks.length = unit(-3, "pt"),
                     axis.text.x = element_text(margin=margin(6,6,10,6,"pt")),
                     axis.text.y = element_text(margin=margin(6,6,10,6,"pt")),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
plot(fvplt)
# print(power.trend(repetitions = 1, x = fvar.data$em, sd = sd(fvar.data$fv), slope = coef(fv.lin.fit)[2]))
# print(power.trend(repetitions = 1, x = fvar.data$em, sd = sd(fvar.data$fv), slope = coef(fv.pl.fit)[2]))
# print(power.trend(repetitions = 1, x = fvar.data$em, sd = sd(fvar.data$fv), slope = 0.015))
# 
# 
# print(fvar.data$es)
# print(fvar.data$em)
# print(fvar.data$ee)
# 

roundUp <- function(x,to=10)
{
  to*(x%/%to + as.logical(x%%to))
}

# print(roundUp(344,to=10))


nbins <- 20
spacing <- ((log10(10)-log10(0.3))/nbins)
# print(spacing)
bins <- seq(from = log10(0.3), to = log10(5), by = spacing)
bins <- 10^bins
bins <- bins*1000
# print(bins)
bins <- roundUp(bins, 1)
# print(bins)
s <- bins[-length(bins)]
e <- bins[-1]
# print(s)
# print(e)





