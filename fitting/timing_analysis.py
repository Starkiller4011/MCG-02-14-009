#!/usr/bin/env python3
# %% General Imports
import os
import numpy as np
import pandas as pd

# %%
energy_ranges = {'300-10000': {'sc':{'time':[],'rate':[],'error':[]},
                               'bg':{'time':[],'rate':[],'error':[]},
                               'cor':{'time':[],'rate':[],'error':[]}},
                 '300-1000': {'sc':{'time':[],'rate':[],'error':[]},
                              'bg':{'time':[],'rate':[],'error':[]},
                              'cor':{'time':[],'rate':[],'error':[]}},
                 '2000-10000': {'sc':{'time':[],'rate':[],'error':[]},
                                'bg':{'time':[],'rate':[],'error':[]},
                                'cor':{'time':[],'rate':[],'error':[]}},
                 '1000-4000': {'sc':{'time':[],'rate':[],'error':[]},
                               'bg':{'time':[],'rate':[],'error':[]},
                               'cor':{'time':[],'rate':[],'error':[]}},
                 '300-425': {'sc':{'time':[],'rate':[],'error':[]},
                             'bg':{'time':[],'rate':[],'error':[]},
                             'cor':{'time':[],'rate':[],'error':[]}},
                 '425-600': {'sc':{'time':[],'rate':[],'error':[]},
                             'bg':{'time':[],'rate':[],'error':[]},
                             'cor':{'time':[],'rate':[],'error':[]}},
                 '600-850': {'sc':{'time':[],'rate':[],'error':[]},
                             'bg':{'time':[],'rate':[],'error':[]},
                             'cor':{'time':[],'rate':[],'error':[]}},
                 '850-1200': {'sc':{'time':[],'rate':[],'error':[]},
                              'bg':{'time':[],'rate':[],'error':[]},
                              'cor':{'time':[],'rate':[],'error':[]}},
                 '1200-1750': {'sc':{'time':[],'rate':[],'error':[]},
                               'bg':{'time':[],'rate':[],'error':[]},
                               'cor':{'time':[],'rate':[],'error':[]}},
                 '1750-2500': {'sc':{'time':[],'rate':[],'error':[]},
                               'bg':{'time':[],'rate':[],'error':[]},
                               'cor':{'time':[],'rate':[],'error':[]}},
                 '2500-3500': {'sc':{'time':[],'rate':[],'error':[]},
                               'bg':{'time':[],'rate':[],'error':[]},
                               'cor':{'time':[],'rate':[],'error':[]}},
                 '3500-5000': {'sc':{'time':[],'rate':[],'error':[]},
                               'bg':{'time':[],'rate':[],'error':[]},
                               'cor':{'time':[],'rate':[],'error':[]}},
                 '5000-7000': {'sc':{'time':[],'rate':[],'error':[]},
                               'bg':{'time':[],'rate':[],'error':[]},
                               'cor':{'time':[],'rate':[],'error':[]}},
                 '7000-10000': {'sc':{'time':[],'rate':[],'error':[]},
                                'bg':{'time':[],'rate':[],'error':[]},
                                'cor':{'time':[],'rate':[],'error':[]}}}
lc_files_path = os.path.abspath('../data/lightcurves')
lc_path = os.path.join(lc_files_path, 'CSV')
obsid = '0550640101'

# %%
for energy_bin in energy_ranges:
    sc_lc_filename = obsid + 'pn_lc_raw_' + energy_bin + '.lc'
    sc_lc_path = os.path.join(lc_path, sc_lc_filename)
    bg_lc_filename = obsid + 'pn_bg_raw_' + energy_bin + '.lc'
    bg_lc_path = os.path.join(lc_path, bg_lc_filename)
    cor_lc_filename = obsid + 'pn_lccor_' + energy_bin + '.lc'
    cor_lc_path = os.path.join(lc_path, cor_lc_filename)
    print(energy_bin)
    print('  ' + str(sc_lc_filename))
    print('    ' + str(sc_lc_path))
    print('  ' + str(bg_lc_filename))
    print('    ' + str(bg_lc_path))
    print('  ' + str(cor_lc_filename))
    print('    ' + str(cor_lc_path))

# %%
