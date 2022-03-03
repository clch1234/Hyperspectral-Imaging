# -*- coding: cp1252 -*-
"""
Script for extracting and plotting the frequencies of vibration of a nanowire as a function of the SEM beam position
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

#workdir = "C:\\Users\\cchardin\\Documents\\Grenoble new"
workdir = "D:\\Documents\\Boulot\\Grenoble"

kwargs = dict(#dat = '20211122',
              #batch = 'Fil 7_3',
              #sens = '01',
              #aller_retour = None, # None or 'both'
              #savefigs = False,
              #savefiles = False,
              plot_spectrums = False,
              old_position_correction = False,
              criteria_for_peak_selection = 'above_noise_max', # Optional, see possible values below
              threshold_peak = 3, # Optional, can be set in params file instead
              required_points_above_threshold = 3, # Optional, can be set in params file instead
              last_quarter_for_threshold = 'auto', # Optional, can be set in params file instead - Values : True, False or 'auto'
              )

"""
Possible values for criteria_for_peak_selection :
    - above_noise_max
    - above_noise_mean
"""
results = mf.map_frequence_batch(**kwargs)

"""
aller_retour_ini = kwargs['aller_retour']
sens_ini = kwargs['sens']
for ii in range(1, 2):
    rectangle = str(ii)
    if kwargs['aller_retour'] is None:
        tu = mf.map_frequence(rectangle, **kwargs)
        F, I, spectrogram, F1, F2, f1, f2, false_positives, X_ini, Y_ini, X, Y, bitmap = tu
    elif kwargs['aller_retour'] == 'both':
        kwargs['aller_retour'] = 'A'
        tu = mf.map_frequence(rectangle, **kwargs)
        FA, IA, spectrogramA, F1A, F2A, f1A, f2A, false_positivesA, X_iniA, Y_iniA, XA, YA, bitmapA = tu

        # Switch 'sens' and 'aller_retour'
        if kwargs['sens'] == '10':
            kwargs['sens'] = '01'
        elif kwargs['sens'] == '01':
            kwargs['sens'] = '10'
        kwargs['aller_retour'] = 'R'
        tu = mf.map_frequence(rectangle, **kwargs)
        FR, IR, spectrogramR, F1R, F2R, f1R, f2R, false_positivesR, X_iniR, Y_iniR, XR, YR, bitmapR = tu

        # Reinitialize for next rectangle
        kwargs['aller_retour'] = aller_retour_ini
        kwargs['sens'] = sens_ini

    else:
        raise ValueError('Wrong value of "aller_retour" : '+str(kwargs['aller_retour']))
"""
show()
