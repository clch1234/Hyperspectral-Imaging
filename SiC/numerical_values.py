# -*- coding: cp1252 -*-
"""
Numerical values for the SiC nanowires we measured

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""

# Length of the NW, in µm
lengths = {'Fil 1':400, # 405, Première mesure 12/2021
           'Fil 2':275,
           'Fil 3':346,
           'Fil 4':340,
           'Fil 5':430,
           'Fil 8':323,
           'Fil 9':367,
           'Fil 4_2':340, # Acq du 20/06/2022
           }

# Diameter of the NW, in nm
diameters = {'Fil 1':285, # 267, Première mesure 12/2021
             'Fil 2':450,
             'Fil 3':323,
             'Fil 4':550,
             'Fil 5':400,
             'Fil 8':380,
             'Fil 9':400, # Measured between 363 nm and 430 nm
             'Fil 4_2':550, # Acq du 20/06/2022
             }

# Step of the manual measurement, in µm
steps = {'Fil 1':40,
         'Fil 2':20,
         'Fil 4':20,
         'Fil 5':20,
         'Fil 8':20,
         'Fil 9':20,
         'Fil 4_2':20, # Acq du 20/06/2022
         }

# Position in µm of position 0
positions_0 = {'Fil 1':3,
               'Fil 5':1,
               'Fil 8':1,
               'Fil 9':3,
               'Fil 4_2':1.5, # Acq du 20/06/2022
               }

# Upper bound for frequency
freqs_max = {'Fil 1':3.5e3,
             'Fil 5':5e3,
             'Fil 8':6e3,
             'Fil 9':6e3,
             'Fil 4_2':8e3, # Acq du 20/06/2022
             }

# Lower bound for frequency
freqs_min = {
             }
