# -*- coding: cp1252 -*-
"""
Script for fitting the spectrums of vibrations of a nanowire with lorentzian shape
Uses functions defined in fit_spectres.py

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
import fit_spectres as fs
from data_RSA_new import *

kwargs = dict(#dat = '20211217',
              #batch = 'Fil 23_2',
              #sens = '10', # '10' or '01'
              #aller_retour = None, # Can be 'A', 'R', or None
              #savefigs = True,
              #savefiles = True,
              plot_fits = True,
              old_position_correction = False,
              guess_gamma = 2*pi*1000,
              downsample = 1,
              )

results = fs.fit_spectrums_batch(**kwargs)
"""
for ii in range(1, 2):
    rectangle = str(ii)
    tu = fs.fit_spectrums(rectangle, **kwargs)
    f1s, f2s, gammas1, gammas2, ampls1, ampls2, offsets, popts, guesses, success, X_ini, Y_ini = tu
"""
show()
