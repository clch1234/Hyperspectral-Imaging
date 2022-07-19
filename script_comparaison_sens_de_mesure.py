"""
Script for plotting speciific acquisition with oposite directions of measurement
Secifically : Fil 7 ; Fil 7_2 ; Fil 8 ; Fil 8_2 (14/12/2021)

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
import os
from functions import *

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

# Arguments that still need to be set
remove_aberrant = True # Remove aberrant frequency values
normalize_from_global_image = False # Use global image taken before acquisition for position normalization

# Select data
dat = '20211214'
rectangle = '1'
AR = None

data = [('Fil 7', '10'),
        ('Fil 7_2', '01'),
        ('Fil 8', '01'),
        ('Fil 8_2', '10')]

# Create figure
fig10, ax10 = subplots()
xlabel('Normalized position y/L')
ylabel('Normalised fresquency shift')# $\\frac{f- \left< f \right>}{\left< f \right>}$')

fig01, ax01 = subplots()
xlabel('Normalized position y/L')
ylabel('Normalised fresquency shift')# $\\frac{f- \left< f \right>}{\left< f \right>}$')

# Do the things
for ii in range(len(data)):
    batch, sens  = data[ii]

    datadir = 'D:\\Documents\\Boulot\\Grenoble\\Data'
    dirname = datadir+u'\\%s\\%s\\Rectangle %s\\'%(dat, batch, rectangle)

    """ Get the data """
    f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = popts_to_fit_data(datadir,
                                                                                       dat,
                                                                                       batch,
                                                                                       rectangle,
                                                                                       sens,
                                                                                       remove_aberrant=remove_aberrant,
                                                                                       AR=AR,
                                                                                       normalize_from_global_image=normalize_from_global_image,
                                                                                       only_two_peaks=False,
                                                                                       select_from_gradient=None,
                                                                                       )

    ff1s, ff2s = (f1s-nanmean(f1s))/nanmean(f1s), (f2s-nanmean(f2s))/nanmean(f2s)

    ax = ax10 if sens == '10' else ax01

    ax.plot(ys1, ff1s, 'o', label=batch+' Mode 1')
    if len(f2s) > 0:
        ax.plot(ys2, ff2s, 'o', label=batch+' Mode 2')

ax10.legend()
ax01.legend()

fig10.tight_layout()
fig01.tight_layout()

show()

figsdir = datadir+u'\\%s\\Comparaison Fils 7 et 8'%dat
if not os.path.isdir(figsdir):
    os.mkdir(figsdir)
fig10.savefig(figsdir+'\\Frequency 10')
fig01.savefig(figsdir+'\\Frequency 01')
