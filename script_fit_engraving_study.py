# -*- coding: cp1252 -*-
"""
Script for fitting frequencies of vibration of a nanowire with the model proposed by Ludovic Bellon

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
##spec = importlib.util.spec_from_file_location("functions_engraving_study", "D:/Documents/Boulot/Scripts/Grenoble/Engraving_study/functions_engraving_study.py")
##fes = importlib.util.module_from_spec(spec)
##spec.loader.exec_module(fes)
import functions_engraving_study as fes
from functions import *
from data_selector import DataSelector

##savefigs = True
##savefiles = True

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

do_fit11 = True
do_fit12 = False

remove_aberrant = True
normalize_from_global_image = False

""" Load or create data files """
dirname = '/'.join(ds.directory.split('/')[:-1])
filesdir = datadir + u'\\%s\\%s\\Data_files'%(dat, batch)
figsdir = datadir + u'\\%s\\%s\\Figures\\Engraving_study'%(dat, batch)
savedir = filesdir+u'\\Fits_%s_%s'%(batch, rectangle)

if normalize_from_global_image:
    print('Using global image for normalization.')
    figsdir += '\\Norm_from_global_image'
    filesdir += '\\Norm_from_global_image'
    savedir += '_Norm_from_global_image'
if AR == 'A':
    figsdir += '_Aller'
elif AR == 'R':
    figsdir += '_Retour'
##if len([p for p in os.path.listdir(datadir+u'\\%s\\%s'%(dat, batch)) if 'Rectangle' in p or 'Réctangle' in p]):
##    figsdir += '_Rectangle %s'%rectangle
if not os.path.isdir(figsdir):
    os.makedirs(figsdir)
if not os.path.isdir(filesdir):
    os.makedirs(filesdir)

prefixe_AR = AR+'_' if AR in ('A', 'R') else ''

if os.path.isfile(savedir+'\\'+prefixe_AR+batch+'_gammas1.txt'):
    f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = load_fit_data(savedir, batch, prefixe_AR)
else:
    save_fit_data(datadir, dat, batch, rectangle, sens, remove_aberrant, AR, normalize_from_global_image)
    f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = load_fit_data(savedir, batch, prefixe_AR)
print('Files loaded')

""" Define fit functions """
if sens == '10':
    dOmega1 = fes.dOmega1_10
elif sens == '01':
    dOmega1 = fes.dOmega1_01

def Omega(yy, Omega0, dR0, dT_max):
    return Omega0 * (1 + dOmega1(yy, dR0, dT_max))
vect_Omega = vectorize(Omega)

if do_fit11 and do_fit12:
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

elif do_fit11 or do_fit12:
    # Fit only one mode
    to_fit = Omega

vect_to_fit = vectorize(to_fit)

""" Do the fits """
if do_fit11 and do_fit12:
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

elif do_fit11:
    y_to_fit = ys1
    Omega_to_fit = 2*pi*f1s
    ini_guess = [2*pi*f1s.mean(), .001, 500]
    bounds = [[2*pi*f1s.min()*.8, 0, 0],
              [2*pi*f1s.max()*1.2, 1, inf]]

    print('Fitting mode 1 only.')
    fit = curve_fit(vect_to_fit,
                    y_to_fit,
                    Omega_to_fit,
                    p0=ini_guess,
                    bounds=bounds)
    Omega10, dR0, dT_max = fit[0]

elif do_fit12:
    y_to_fit = ys2
    Omega_to_fit = 2*pi*f2s
    ini_guess = [2*pi*f2s.mean(), .001, 500]
    bounds = [[2*pi*f2s.min()*.8, 0, 0],
              [2*pi*f2s.max()*1.2, 1, inf]]

    print('Fitting mode 2 only.')
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
if do_fit11:
    y1 = linspace(ys1.min(), ys1.max(), 51)
    val1 = vect_Omega(y1, Omega10, dR0, dT_max)/2/pi
    val1_full = vect_Omega(yy, Omega10, dR0, dT_max)/2/pi
if do_fit12:
    y2 = linspace(ys2.min(), ys2.max(), 51)
    val2 = vect_Omega(y2, Omega20, dR0, dT_max)/2/pi
    val2_full = vect_Omega(yy, Omega20, dR0, dT_max)/2/pi

fig_fit, ax_fit = subplots()
if do_fit11:
    plot(ys1, f1s, 'o', c='C0', label='Mode 1.1')
    plot(yy, val1_full, 'r', label='Fit 1.1')
    axhline(Omega10/2/pi, color='C0', label='$\Omega_{1, 0}/2\pi = %f$'%(Omega10/2/pi))
if do_fit12:
    plot(ys2, f2s, 'o', c='C1', label='Mode 1.2')
    plot(yy, val2_full, 'g', label='Fit 1.2')
    axhline(Omega20/2/pi, color='C1', label='$\Omega_{2, 0}/2\pi = %f$'%(Omega20/2/pi))
suptitle('%s (%s)'%(batch, rectangle))
title('$\delta_{R_0} = %.6f$ - $\Delta T_{max} = %i$'%(dR0, dT_max))
xlabel('Beam position (normalized)')
ylabel('Frequency (Hz)')
legend()
tight_layout()

if do_fit11:
    fig1, ax1 = subplots()
    plot(ys1, f1s, 'o', label='Mode 1.1')
    plot(y1, val1, 'r', label='Fit 1.1')
    suptitle('%s (%s)'%(batch, rectangle))
    title('$\Omega_{1, 0}/2\pi = %i$ - $\delta_{R_0} = %.6f$ - $\Delta T_{max} = %i$'%(Omega10/2/pi, dR0, dT_max))
    xlabel('Beam position (normalized)')
    ylabel('Frequency (Hz)')
    legend()
    tight_layout()

if do_fit12:
    fig2, ax2 = subplots()
    plot(ys2, f2s, 'o', label='Mode 1.2')
    plot(y2, val2, 'r', label='Fit 1.2')
    suptitle('%s (%s)'%(batch, rectangle))
    title('$\Omega_{2, 0}/2\pi = %i$ - $\delta_{R_0} = %.6f$ - $\Delta T_{max} = %i$'%(Omega20/2/pi, dR0, dT_max))
    xlabel('Beam position (normalized)')
    ylabel('Frequency (Hz)')
    legend()
    tight_layout()

fig_contributions, ax_contributions = subplots()
dOmega1_R = fes.dOmega1_R_10 if sens == '10' else fes.dOmega1_R_01
valR = array([dOmega1_R(y, dR0) for y in yy])
valE = array([fes.dOmega1_T(y, dT_max) for y in yy])
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
    if do_fit11:
        fig1.savefig(figsdir+'\\Mode_1.png')
    if do_fit12:
        fig2.savefig(figsdir+'\\Mode_2.png')
    print('Figures saved')

if savefiles:
    dd = {'dR0':dR0, 'dT_max':dT_max}
    if do_fit11:
        dd['Omega10'] = Omega10
    else:
        dd['Omega10'] = nan
    if do_fit12:
        dd['Omega20'] = Omega20
    else:
        dd['Omega20'] = nan
    with open(filesdir+'\\fit_engraving_study.json', 'w') as ff:
        json.dump(dd, ff)
    print('Fit values saved')
