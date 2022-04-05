# -*- coding: cp1252 -*-
"""
Script for extracting and plotting the frequencies of vibration of a nanowire as a function of the SEM beam position
For a series of measurement studying the effect of defocus on the sensitivity of our setup.
Uses functions defined in map_frequence.py

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
import map_frequence as mf
from data_RSA_new import *
import os
# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

kwargs = dict(dat = '20211214',
              workdir="D:\\Documents\\Boulot\\Grenoble",
              datadir="D:\\Documents\\Boulot\\Grenoble\\Data",
              #batch = 'Fil 7_3',
              sens = '10',
              aller_retour = None, # None or 'both'
              savefigs = False,
              savefiles = False,
              plot_spectrums = False,
              old_position_correction = False,
              criteria_for_peak_selection = 'above_noise_max', # Optional, see possible values below
              threshold_peak = 3, # Optional, can be set in params file instead
              required_points_above_threshold = 3, # Optional, can be set in params file instead
              last_quarter_for_threshold = 'auto', # Optional, can be set in params file instead - Values : True, False or 'auto'
              )

batches = ['Fil 9_2', 'Fil 9_3', 'Fil 9_4', 'Fil 9_5', 'Fil 9_6',
           'Fil 9_7', 'Fil 9_8', 'Fil 9_9', 'Fil 9_10', 'Fil 9_11',
           'Fil 9_14', 'Fil 9_15']

on_focus = 5.070911 # Working distance when focused on the tip of the NW

WDs = [5.070911, 5.072, 5.073, 5.074, 5.075, 5.077, 5.079, 5.081, 5.083, 5.085,
       5.086, 5.078465] # Last: focused at clamping

freq_middle = 850e3 # Frequency between the two modes

res = []
for ii, batch in enumerate(batches):
    kwargs['batch'] = batch
    rectangle = '1'

    tu = mf.map_frequence(rectangle, **kwargs)
    close('all')
    res.append(tu)

fig_all1, ax_all1 = subplots()
fig_all2, ax_all2 = subplots()

f1s, f2s = [], []
df1, df2 = [], []
defoc = []
for ii, tu in enumerate(res):
    wd = WDs[ii]
    defocus = wd - on_focus

    f1, f2 = tu[5], tu[6]

    # Correct frequency inversion
    for jj in range(len(f1)):
        if not isnan(f1[jj]) and f1[jj] > freq_middle \
           or not isnan(f2[jj]) and f2[jj] < freq_middle \
           or not isnan(f1[jj]) and not isnan(f2[jj]) and f2[jj] < f1[jj]:
            temp = f1[jj]
            f1[jj] = f2[jj]
            f2[jj] = temp

    ff1 = f1[where(logical_not(isnan(f1)))[0]]
    ff2 = f2[where(logical_not(isnan(f2)))[0]]

    f1s.append(ff1.mean())
    f2s.append(ff2.mean())
    df1.append(ff1.var())
    df2.append(ff2.var())

    defoc.append(defocus/1e3)

    ax_all1.plot(ff1, 'o', label=str(wd))
    ax_all2.plot(ff2, 'o', label=str(wd))

ax_all1.legend()
ax_all2.legend()

""" Plot figures """
fig_freqs, ax_freqs = subplots()
xlabel('Defocalisation (mm)')
ylabel('Frequency (Hz)')

errorbar(defoc, f1s, yerr=df1, marker='o', linestyle='')
errorbar(defoc, f2s, yerr=df2, marker='o', linestyle='')

""" Save figures """
fig_freqs.savefig("D:\\Documents\\Boulot\\Grenoble\\Data\\%s\\defocus"%kwargs['dat'])
fig_all1.savefig("D:\\Documents\\Boulot\\Grenoble\\Data\\%s\\defocus_mode1"%kwargs['dat'])
fig_all2.savefig("D:\\Documents\\Boulot\\Grenoble\\Data\\%s\\defocus_mode2"%kwargs['dat'])
print("Figures saved")

show()
