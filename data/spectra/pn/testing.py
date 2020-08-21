#!/usr/bin/env python3
from __future__ import print_function,division
import astropy.io.fits as pyfits

rmf = pyfits.open("101_pn.rmf")

print(rmf['EBOUNDS'].data)


def read_rmf(rmf_path, energies):
	'''Finds channels corresponding to energy bins using response matrix'''
	rmf_file=pyfits.open(rmf_path)
	ebounds = rmf_file['EBOUNDS'].data
	emins=[]
	emaxs=[]
	channels=[]
	# Find channels corresponding to energy bins
	for row in ebounds:
		channels.append(row[0])
		emins.append(row[1])
		emaxs.append(row[2])
	channel_bins=[]
	for e in energies:
		for i in range(1,len(channels)-1):
			if emins[i]<=e<emins[i+1]:
				channel_bins.append(channels[i])
	
	return channel_bins
