#!/usr/bin/env python3
# General Imports
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import astropy
from xspec import *

# Load local models and set defaults
AllModels.lmod('relxill')
AllModels.lmod('thcomp')
Xset.abund = 'wilm'
Fit.nIterations = 100000
Fit.statMethod = 'cstat'

# Data to load
sc_spectra = ['spectra/mos/101_m1_sc_opt.pha']
bg_spectra = ['spectra/mos/101_m1_bg_opt.pha']
responses = ['spectra/mos/101_m1.rmf']
arfs = ['spectra/mos/101_m1.arf']
# sc_spectra = ['spectra/pn/101_sc_opt.pha']
# bg_spectra = ['spectra/pn/101_bg_opt.pha']
# responses = ['spectra/pn/101_pn.rmf']
# arfs = ['spectra/pn/101_pn.arf']
# sc_spectra = ['spectra/pn/701_sc_opt.pha']
# bg_spectra = ['spectra/pn/701_bg_opt.pha']
# responses = ['spectra/pn/701_pn.rmf']
# arfs = ['spectra/pn/701_pn.arf']
# sc_spectra = ['spectra/xis/oxi03sc.pha']
# bg_spectra = ['spectra/xis/oxi03bg.pha']
# responses = ['spectra/xis/xi03.rsp']
# arfs = [None]
# sc_spectra = ['spectra/pn/101_sc_opt.pha', 'spectra/pn/701_sc_opt.pha', 'spectra/mos/101_m1_sc_opt.pha', 'spectra/mos/101_m2_sc_opt.pha', 'spectra/mos/701_m1_sc_opt.pha', 'spectra/mos/701_m2_sc_opt.pha']
# bg_spectra = ['spectra/pn/101_bg_opt.pha', 'spectra/pn/701_bg_opt.pha', 'spectra/mos/101_m1_bg_opt.pha', 'spectra/mos/101_m2_bg_opt.pha', 'spectra/mos/701_m1_bg_opt.pha', 'spectra/mos/701_m2_bg_opt.pha']
# responses = ['spectra/pn/101_pn.rmf', 'spectra/pn/701_pn.rmf', 'spectra/mos/101_m1.rmf', 'spectra/mos/101_m2.rmf', 'spectra/mos/701_m1.rmf', 'spectra/mos/701_m2.rmf']
# arfs = ['spectra/pn/101_pn.arf', 'spectra/pn/701_pn.arf', 'spectra/mos/101_m1.arf', 'spectra/mos/101_m2.arf', 'spectra/mos/701_m1.arf', 'spectra/mos/701_m2.arf']

def load_spectra(sc_spectra, bg_spectra, responses, arfs, sources):
    spectra = {}
    load_strings = []
    for i in range(0, sources):
        load_string = str(i+1) + ':' + str(i+1) + ' ' + str(sc_spectra[i])
        load_strings.append(load_string)
    for i in range(0, sources):
        load_string = str(sources+i+1) + ':' + str(sources+i+1) + ' ' + str(bg_spectra[i])
        load_strings.append(load_string)
    data_string = ''
    for load_string in load_strings:
        data_string += (load_string + ' ')
    AllData(data_string)
    for i in range(0, sources):
        spectra[i] = AllData(i+1)
        spectra[sources+i] = AllData(sources+i+1)
    for i in range(0, sources):
        for j in range(0, sources):
            if j == i:
                spectra[i].multiresponse[j] = responses[i]
                spectra[i].multiresponse[j].arf = arfs[i]
                spectra[sources+i].multiresponse[j] = responses[i]
                spectra[sources+i].multiresponse[j].arf = None
                spectra[i].multiresponse[sources+j] = responses[i]
                spectra[i].multiresponse[sources+j].arf = None
                spectra[sources+i].multiresponse[sources+j] = responses[i]
                spectra[sources+i].multiresponse[sources+j].arf = None
            else:
                spectra[i].multiresponse[j] = None
                spectra[sources+i].multiresponse[j] = None
                spectra[i].multiresponse[sources+j] = None
                spectra[sources+i].multiresponse[sources+j] = None
    return spectra

def plot_ratio(energy, ratio, energy_err, ratio_err, lims=True):
    avg_err = np.average(ratio_err)
    max_err = np.max(ratio_err)
    upper_slim = 1 + avg_err
    lower_slim = 1 - avg_err
    upper_lim = 1 + max_err
    lower_lim = 1 - max_err
    fig1, ax1 = plt.subplots()
    ax1.plot(energy, ratio, linestyle='-', color='r', linewidth=1)
    ax1.errorbar(energy, ratio, xerr=energy_err, yerr=ratio_err, color='r', fmt='none', elinewidth=1)
    ax1.axhline(y=1, color='lime', linestyle='-', linewidth=1)
    if lims:
        ax1.axhline(y=upper_slim, color='orange', linestyle='-', linewidth=1)
        ax1.axhline(y=lower_slim, color='orange', linestyle='-', linewidth=1)
        ax1.axhline(y=upper_lim, color='cyan', linestyle='-', linewidth=1)
        ax1.axhline(y=lower_lim, color='cyan', linestyle='-', linewidth=1)
    ax1.set_xscale('log')
    ax1.set_xticks([0.5,1,2,5,10])
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    return fig1

# Load data
print('#######################################')
sources = len(sc_spectra)
spectra = load_spectra(sc_spectra, bg_spectra, responses, arfs, sources)
print('#######################################')
AllData.show()
print('#######################################')

# Set source model
sm_s = Model('const*tbabs*po', 'source', 1)
sm_b = AllModels(2, 'source')
pars = {1: '1 -1', 2: '0.135 -1', 3: 1.89, 4: 1e-4}
sm_s.setPars(pars)
sm_b(1).values = 0.0
sm_b(1).frozen = True

# Generate background model
spectra[0].ignore('0.-2.0 5.0-**')
spectra[1].ignore('0.-0.3 10.-**')
params = 0
base_model = 'bknpo'
back_model = {}
bm = Model(base_model, 'back', 2)
bm.setPars({1:1.89, 2:5.0, 3:1.89, 4:1e-4})
for par in range(0, bm.nParameters):
    back_model[par+1] = bm(par+1).values[0]
    params += 1
print(back_model)
Fit.perform()
for par in range(0, bm.nParameters):
    back_model[par+1] = bm(par+1).values[0]
print(back_model)
Plot.device = '/xs'
Plot.xAxis = 'keV'
# Plot.add = True
Plot.setGroup('1 2')
Plot('ra')
energy = np.asarray(Plot.x(2))
energy_err = np.asarray(Plot.xErr(2))
ratio = np.asarray(Plot.y(2))
ratio_err = np.asarray(Plot.yErr(2))
avg_err = np.average(ratio_err)
max_err = np.max(ratio_err)
upper_lim = 1 + avg_err
lower_lim = 1 - avg_err
fig = plot_ratio(energy, ratio, energy_err, ratio_err)
Plot.device = '/xs'
Plot.xAxis = 'keV'
# Plot.add = True
Plot.setGroup('1 2')
Plot('ra')
plt.show()
Plot.device = '/null'
Plot.xAxis = 'keV'
# Plot.add = True
Plot.setGroup('1 2')
Plot('ra')
i = 0
while i < len(ratio):
    if ratio[i] > upper_lim:
        j = i
        r = ratio[j]
        e = []
        while r > upper_lim:
            e.append(energy[j])
            j += 1
            r = ratio[j]
        lineE = np.average(e)
        print(params)
        epsilon = 0.2
        lineE_str = str(lineE) + ' 0.01 ' + str(lineE-epsilon) + ' ' + str(lineE-epsilon) + ' ' + str(lineE+epsilon) + ' ' + str(lineE+epsilon)
        print(lineE_str)
        back_model[params+1] = lineE_str
        back_model[params+2] = '0.001 0.001 0 0 0.01 0.01'
        back_model[params+3] = 1e-4
        params += 3
        print(back_model)
        base_model += '+gauss'
        AllModels -= 'back'
        bm = Model(base_model, 'back', 2)
        bm.setPars(back_model)
        AllModels.show()
        Fit.perform()
        i = j + 1
    i += 1
for par in range(1, bm.nParameters+1):
    bm(par).frozen = True
Plot.device = '/xs'
Plot.xAxis = 'keV'
# Plot.add = True
Plot.setGroup('1 2')
Plot('ra')
energy = np.asarray(Plot.x(2))
energy_err = np.asarray(Plot.xErr(2))
ratio = np.asarray(Plot.y(2))
ratio_err = np.asarray(Plot.yErr(2))
avg_err = np.average(ratio_err)
max_err = np.max(ratio_err)
upper_lim = 1 + max_err
lower_lim = 1 - max_err
fig = plot_ratio(energy, ratio, energy_err, ratio_err)
Plot.device = '/xs'
Plot.xAxis = 'keV'
# Plot.add = True
Plot.setGroup('1 2')
Plot('ra')
plt.show()

models = {'const*tbabs*(po+bb)':{'pars': {1: '1 -1',
                                          2: '0.135 -1',
                                          3: '1.89', 4: 1e-4,
                                          5: 0.1, 6: 1e-4}},
        'const*tbabs*(po+bb+zgauss)':{'pars': {1: '1 -1',
                                               2: '0.135 -1',
                                               3: '1.89', 4: 1e-4,
                                               5: 0.1, 6: 1e-4,
                                               7:'6.4 -1', 8: '0.001 0.001 0 0 0.01 0.01', 9: '0.02845 -1', 10: 1e-4}},
        'const*tbabs*(po+bb+zgauss+zgauss)':{'pars': {1: '1 -1',
                                                      2: '0.135 -1',
                                                      3: '1.89', 4: 1e-4,
                                                      5: 0.1, 6: 1e-4,
                                                      7:'6.4 -1', 8: '0.001 0.001 0 0 0.01 0.01', 9: '0.02845 -1', 10: 1e-4,
                                                      11:'6.7 -1', 12: '0.001 0.001 0 0 0.01 0.01', 13: '0.02845 -1', 14: 1e-4}},
        'const*tbabs*(po+bb+zgauss+zgauss+zgauss)':{'pars': {1: '1 -1',
                                                             2: '0.135 -1',
                                                             3: '1.89', 4: 1e-4,
                                                             5: 0.1, 6: 1e-4,
                                                             7:'6.4 -1', 8: '0.001 0.001 0 0 0.01 0.01', 9: '0.02845 -1', 10: 1e-4,
                                                             11:'6.7 -1', 12: '0.001 0.001 0 0 0.01 0.01', 13: '0.02845 -1', 14: 1e-4,
                                                             15:'6.97 -1', 16: '0.001 0.001 0 0 0.01 0.01', 17: '0.02845 -1', 18: 1e-4}},
        'const*tbabs*(po+bb+bb+zgauss+zgauss+zgauss)':{'pars': {1: '1 -1',
                                                                2: '0.135 -1',
                                                                3: '1.89', 4: 1e-4,
                                                                5: 0.1, 6: 1e-4,
                                                                7: 0.1, 8: 1e-4,
                                                                9:'6.4 -1', 10: '0.001 0.001 0 0 0.01 0.01', 11: '0.02845 -1', 12: 1e-4,
                                                                13:'6.7 -1', 14: '0.001 0.001 0 0 0.01 0.01', 15: '0.02845 -1', 16: 1e-4,
                                                                17:'6.97 -1', 18: '0.001 0.001 0 0 0.01 0.01', 19: '0.02845 -1', 20: 1e-4}},
        'const*tbabs*(thcomp*diskbb)':{'pars': {1: '1 -1',
                                                2: '0.135 -1',
                                                3: '1.89', 4: 100, 5: 1, 6: '0.02845 -1',
                                                7: 0.1, 8: 1}},
        'const*tbabs*((thcomp*diskbb)+zgauss)':{'pars': {1: '1 -1',
                                                         2: '0.135 -1',
                                                         3: '1.89', 4: 100, 5: 1, 6: '0.02845 -1',
                                                         7: 0.1, 8: 1,
                                                         9:'6.4 -1', 10: '0.001 0.001 0 0 0.01 0.01', 11: '0.02845 -1', 12: 1e-4}},
        'const*tbabs*((thcomp*diskbb)+zgauss+zgauss)':{'pars': {1: '1 -1',
                                                                2: '0.135 -1',
                                                                3: '1.89', 4: 100, 5: 1, 6: '0.02845 -1',
                                                                7: 0.1, 8: 1,
                                                                9:'6.4 -1', 10: '0.001 0.001 0 0 0.01 0.01', 11: '0.02845 -1', 12: 1e-4,
                                                                13:'6.7 -1', 14: '0.001 0.001 0 0 0.01 0.01', 15: '0.02845 -1', 16: 1e-4}},
        'const*tbabs*((thcomp*diskbb)+zgauss+zgauss+zgauss)':{'pars': {1: '1 -1',
                                                                       2: '0.135 -1',
                                                                       3: '1.89', 4: 100, 5: 1, 6: '0.02845 -1',
                                                                       7: 0.1, 8: 1,
                                                                       9:'6.4 -1', 10: '0.001 0.001 0 0 0.01 0.01', 11: '0.02845 -1', 12: 1e-4,
                                                                       13:'6.7 -1', 14: '0.001 0.001 0 0 0.01 0.01', 15: '0.02845 -1', 16: 1e-4,
                                                                       17:'6.97 -1', 18: '0.001 0.001 0 0 0.01 0.01', 19: '0.02845 -1', 20: 1e-4}}}

spectra[0].notice('all')
spectra[1].notice('all')
spectra[0].ignore('0.-0.3 10.-**')
spectra[1].ignore('0.-0.3 10.-**')

for model in models:
    print(model)
    AllModels -= 'source'
    sm_s = Model(model, 'source', 1)
    sm_s.setPars(models[model]['pars'])
    sm_b = AllModels(2, 'source')
    sm_b(1).values = 0.0
    sm_b(1).frozen = True
    Fit.perform()
    for par in range(1, sm_s.nParameters + 1):
        models[model]['pars'][par] = sm_s(par).values
    models[model]['cstat'] = Fit.statistic
    models[model]['dof'] = Fit.dof
    Plot.device = '/xs'
    Plot.xAxis = 'keV'
    Plot.add = True
    Plot('ld ra')
    # Plot.device = '/null'
    # Plot.xAxis = 'keV'
    # Plot.setGroup('1 2')
    # Plot('ra')
    # energy = np.asarray(Plot.x(2))
    # energy_err = np.asarray(Plot.xErr(2))
    # ratio = np.asarray(Plot.y(2))
    # ratio_err = np.asarray(Plot.yErr(2))
    # fig = plot_ratio(energy, ratio, energy_err, ratio_err, lims=False)
    # plt.show()

print('\n##############################\nFits complete')
for model in models:
    print(str(model) + ': cstat: ' + ('%.2f' % models[model]['cstat']) + ', dof: ' + str(models[model]['dof']))