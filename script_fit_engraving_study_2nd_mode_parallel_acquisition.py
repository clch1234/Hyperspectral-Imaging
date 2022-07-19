# -*- coding: cp1252 -*-
"""
Script for fitting frequencies of vibration of a nanowire with the model proposed by Ludovic Bellon
Here on the 2nd mechanical mode, from the data obtained with the "old" continuous acquisition method

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
from scipy.optimize import curve_fit
import os
import json
import importlib.util
import functions_engraving_study as fes
import parameters_parallel_acquisition as ppa
from functions import *
from data_selector import DataSelector

# Data selector
ds = DataSelector(description="Select a single rectangle")
ds.select_directory()

rectangle = ds.directory.split('/')[-1].split(' ')[1]
batch = ds.directory.split('/')[-2]
dat = ds.directory.split('/')[-3]
datadir = '/'.join(ds.directory.split('/')[:-3])

sens = ds.sens_var.get()
if not sens in ('10', '01'):
    raise ValueError('Sens should be either "10" or "01", %s received.'%sens)
AR = ds.AR_var.get()
if AR == 'None':
    AR = None
if not AR in ('A', 'R', None):
    raise ValueError('AR should be either "A", "R" or None, %s received.'%AR)

savefigs = ds.savefigs_var.get()
savefiles = ds.savefiles_var.get()

##sens = '10' # '10' or '01'
##AR = None # Can be 'A', 'R' or None

do_fit21 = False
do_fit22 = True

remove_aberrant = True

""" Load data files """
for key, corres in ppa.correspondences.items():
    if corres == batch:
        break
else:
    raise ValueError('No correspondence found for batch = %s !'%batch)

dirname = '/'.join(ds.directory.split('/')[:-1])
filesdir = datadir + u'\\%s\\%s\\Data_files\\Parallel_acquisition'%(dat, batch)
figsdir = datadir + u'\\%s\\%s\\Figures\\Engraving_study\\Parallel_acquisition'%(dat, batch)
savedir = filesdir+u'\\Fits_%s_%s'%(batch, rectangle)

if AR == 'A':
    figsdir += '_Aller'
elif AR == 'R':
    figsdir += '_Retour'

if not os.path.isdir(figsdir):
    os.makedirs(figsdir)
if not os.path.isdir(filesdir):
    os.makedirs(filesdir)

f1s = load(filesdir+'\\freqs1.npy', allow_pickle=True)
ys1 = load(filesdir+'\\ys1.npy', allow_pickle=True)

f2s = load(filesdir+'\\freqs2.npy', allow_pickle=True)
ys2 = load(filesdir+'\\ys2.npy', allow_pickle=True)

print('Files loaded')

""" Normalize y by length """
# Load X, Y and bitmap
bitmap = loadtxt(dirname+u'\\Rectangle 1\\Bitmap.txt')
X_ini = loadtxt(dirname+u'\\Data_files\\%s_rectangle_1_raw_X.txt'%corres)
Y_ini = loadtxt(dirname+u'\\Data_files\\%s_rectangle_1_raw_Y.txt'%corres)

# Detection of the edges of the NW and left-right separation
centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(bitmap, X_ini)

# Deduce basis position and NW length
length = abs(Y_ini[first_line_found, 0] - Y_ini[last_line_found, 0])
basis = Y_ini[first_line_found, 0]

ys1 = (ys1 - basis) / length
ys2 = (ys2 - basis) / length

""" Define fit functions """
if sens == '10':
    dOmega2 = fes.dOmega2_10
elif sens == '01':
    dOmega2 = fes.dOmega2_01

def Omega(yy, Omega0, dR0, dT_max):
    return Omega0 * (1 + dOmega2(yy, dR0, dT_max))
vect_Omega = vectorize(Omega)

if do_fit21 and do_fit22:
    # Fit both modes with shared parameters
    def get_merged_y_and_f(ys1, ys2, f1s, f2s):
        y_merged = concatenate((ys1, ys2))
        f_merged = concatenate((f1s, f2s))
        return y_merged, f_merged

    def to_fit(yy, Omega10, Omega20, dR0, dT_max):
        y1 = yy[:len(ys1)]
        y2 = yy[len(ys1):]

        res1 = vect_Omega(y1, Omega10, dR0, dT_max)
        res2 = vect_Omega(y2, Omega20, dR0, dT_max)

        return concatenate((res1, res2))

elif do_fit21 or do_fit22:
    # Fit only one mode
    to_fit = Omega

vect_to_fit = vectorize(to_fit)

""" Do the fits """
if do_fit21 and do_fit22:
    y_to_fit, f_to_fit = get_merged_y_and_f(ys1, ys2, f1s, f2s)
    Omega_to_fit = 2*pi*f_to_fit
    ini_guess = [2*pi*f1s.mean(), 2*pi*f2s.mean(), .001, 500]
    bounds = [[2*pi*f1s.min()*.8, 2*pi*f2s.min()*.8, 0, 0],
              [2*pi*f1s.max()*1.2, 2*pi*f2s.max()*1.2, 1, inf]]

    print('Fitting both modes with shared parameters.')
    fit = curve_fit(to_fit, #vect_to_fit,
                    y_to_fit,
                    Omega_to_fit,
                    p0=ini_guess,
                    bounds=bounds)
    Omega10, Omega20, dR0, dT_max = fit[0]

elif do_fit21:
    y_to_fit = ys1
    Omega_to_fit = 2*pi*f1s
    ini_guess = [2*pi*f1s.mean(), .001, 500]
    bounds = [[2*pi*f1s.min()*.8, 0, 0],
              [2*pi*f1s.max()*1.2, 1, inf]]

    print('Fitting mode 2.1 only.')
    fit = curve_fit(vect_to_fit,
                    y_to_fit,
                    Omega_to_fit,
                    p0=ini_guess,
                    bounds=bounds)
    Omega10, dR0, dT_max = fit[0]

elif do_fit22:
    y_to_fit = ys2
    Omega_to_fit = 2*pi*f2s
    ini_guess = [2*pi*f2s.mean(), .001, 500]
    bounds = [[2*pi*f2s.min()*.8, 0, 0],
              [2*pi*f2s.max()*1.2, 1, inf]]

    print('Fitting mode 2.2 only.')
    fit = curve_fit(vect_to_fit,
                    y_to_fit,
                    Omega_to_fit,
                    p0=ini_guess,
                    bounds=bounds)
    Omega20, dR0, dT_max = fit[0]

print('Done !')

""" Plot results """
print('Plotting figures.')
yy = linspace(0, 1, 51)
if do_fit21:
    y1 = linspace(ys1.min(), ys1.max(), 51)
    val1 = vect_Omega(y1, Omega10, dR0, dT_max)/2/pi
    val1_full = vect_Omega(yy, Omega10, dR0, dT_max)/2/pi
if do_fit22:
    y2 = linspace(ys2.min(), ys2.max(), 51)
    val2 = vect_Omega(y2, Omega20, dR0, dT_max)/2/pi
    val2_full = vect_Omega(yy, Omega20, dR0, dT_max)/2/pi

fig_fit, ax_fit = subplots()
if do_fit21:
    plot(ys1, f1s, 'o', c='C0', label='Mode 2.1')
    plot(yy, val1_full, 'r', label='Fit 2.1')
    axhline(Omega10/2/pi, color='C0', label='$\Omega_{1, 0}/2\pi = %f$'%(Omega10/2/pi))
if do_fit22:
    plot(ys2, f2s, 'o', c='C1', label='Mode 2.2')
    plot(yy, val2_full, 'g', label='Fit 2.2')
    axhline(Omega20/2/pi, color='C1', label='$\Omega_{2, 0}/2\pi = %f$'%(Omega20/2/pi))
suptitle('%s (%s)'%(batch, rectangle))
title('$\delta_{R_0} = %.6f$ - $\Delta T_{max} = %i$'%(dR0, dT_max))
xlabel('Beam position (normalized)')
ylabel('Frequency (Hz)')
legend()
tight_layout()

if do_fit21:
    fig1, ax1 = subplots()
    plot(ys1, f1s, 'o', label='Mode 2.1')
    plot(y1, val1, 'r', label='Fit 2.1')
    suptitle('%s (%s)'%(batch, rectangle))
    title('$\Omega_{1, 0}/2\pi = %i$ - $\delta_{R_0} = %.6f$ - $\Delta T_{max} = %i$'%(Omega10/2/pi, dR0, dT_max))
    xlabel('Beam position (normalized)')
    ylabel('Frequency (Hz)')
    legend()
    tight_layout()

if do_fit22:
    fig2, ax2 = subplots()
    plot(ys2, f2s, 'o', label='Mode 2.2')
    plot(y2, val2, 'r', label='Fit 2.2')
    suptitle('%s (%s)'%(batch, rectangle))
    title('$\Omega_{2, 0}/2\pi = %i$ - $\delta_{R_0} = %.6f$ - $\Delta T_{max} = %i$'%(Omega20/2/pi, dR0, dT_max))
    xlabel('Beam position (normalized)')
    ylabel('Frequency (Hz)')
    legend()
    tight_layout()

fig_contributions, ax_contributions = subplots()
dOmega2_R = fes.dOmega2_R_10 if sens == '10' else fes.dOmega2_R_01
valR = array([dOmega2_R(y, dR0) for y in yy])
valE = array([fes.dOmega2_T(y, dT_max) for y in yy])
plot(yy, valR, label='$d\Omega_R$')
plot(yy, valE, label='$d\Omega_E$')
plot(yy, valR + valE, label='$d\Omega_{tot}$')
axhline(0, color='k')
suptitle('%s (%s)'%(batch, rectangle))
title('$\delta_{R_0} = %.6f$ - $\Delta T_{max} = %i$'%(dR0, dT_max))
xlabel('Beam position (normalized)')
ylabel('$d\Omega$ (relative)')
legend()
tight_layout()

show()

if savefigs:
    fig_fit.savefig(figsdir+'\\Fits.png')
    fig_contributions.savefig(figsdir+'\\Contributions.png')
    if do_fit21:
        fig1.savefig(figsdir+'\\Mode_1.png')
    if do_fit22:
        fig2.savefig(figsdir+'\\Mode_2.png')
    print('Figures saved')

if savefiles:
    dd = {'dR0':dR0, 'dT_max':dT_max}
    if do_fit21:
        dd['Omega10'] = Omega10
    else:
        dd['Omega10'] = nan
    if do_fit22:
        dd['Omega20'] = Omega20
    else:
        dd['Omega20'] = nan
    with open(filesdir+'\\fit_engraving_study.json', 'w') as ff:
        json.dump(dd, ff)
    print('Fit values saved')
