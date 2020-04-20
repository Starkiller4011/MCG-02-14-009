#!/usr/bin/env python3
# %% General Imports
import os
import numpy as np
import pandas as pd
from astropy.io import fits

# %%
lc_path = os.path.abspath('../data/lightcurves')
lc_fits_path = os.path.join(lc_path, 'FITS')
lc_csv_path = os.path.join(lc_path, 'CSV')
for fits_file in os.listdir(lc_fits_path):
    path = os.path.join(lc_fits_path, fits_file)
    name = fits_file[:-5]
    lc_filename = name + '.lc'
    with fits.open(path, memmap=True) as hdu_list:
        time = hdu_list[1].data['TIME'].tolist().copy()
        rate = hdu_list[1].data['RATE'].tolist().copy()
        error = hdu_list[1].data['ERROR'].tolist().copy()
    outfile_path = os.path.join(lc_csv_path,lc_filename)
    with open(outfile_path, 'w') as outfile:
        for i in range(0,len(time)):
            outfile.write(str(time)+','+str(rate)+','+str(error)+'\n')
# %% testing
lc_name = '0550640101pn_bg_raw_300-425.fits'
path = os.path.join(lc_fits_path, lc_name)
name = lc_name[:-5]
lc_filename = name + '.lc'
with fits.open(path, memmap=True) as hdu_list:
    obs_info = hdu_list[0].header
    cols = hdu_list[1].columns
    print(obs_info)
    for key in obs_info:
        print(key + ': ' + str(obs_info[key]))
    for key in hdu_list[1].header:
        print(key + ': ' + str(hdu_list[1].header[key]))
    time = hdu_list[1].data['TIME'].tolist().copy()
    flux = hdu_list[1].data['RATE'].tolist().copy()
    error = hdu_list[1].data['ERROR'].tolist().copy()
    print(float(hdu_list[0].header['MJDREF']) * 86400)
    print(time[0])

# %%
print('** Done ** All lightcurves converted')