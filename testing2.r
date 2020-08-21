#!/usr/bin/env Rscript

################### Setup ###################
options(scipen=999)

########## Load required libraries ##########
require(FITSio)
require(sour)
require(ggplot2)
require(AstroR)
require(splines)

require(reticulate)
use_python("c:/Program Files/Python38/python.exe")

pyfits <- import("astropy.oi.fits")

arf.path <- "data/spectra/pn/101_pn.arf"
rmf.path <- "data/spectra/pn/101_pn.rmf"

arf <- xmm.arf(arf.path)
t <- makeFITSimHdr()
rmf <- readFITS(rmf.path, hdu = 2)

energies <- list()
areas <- list()

for(i in 1:length(arf$SPECRESP)) {
  energies[[i]] <- ((arf$ENERGY_LO[[i]]+arf$ENERGY_HI[[i]])/2)
  areas[[i]] <- arf$SPECRESP[[i]]
}

energies <- unlist(energies)
areas <- unlist(areas)

area.spline <- lm(areas ~ bs(energies))
print(summary(area.spline))

binned.areas <- list()

