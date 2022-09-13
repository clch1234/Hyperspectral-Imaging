# -*- coding: cp1252 -*-
"""
Script for plotting the lorentzian fits of all spectrums measured for a nanowire

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from data_RSA_new import *
from functions import *
from data_selector import DataSelector
from fit_spectres import multi_lorentz
import datetime

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

# Arguments that still need to be set
remove_aberrant = True # Remove aberrant frequency values
normalize_from_global_image = False # Use global image taken before acquisition for position normalization
old_position_correction = False # Normally never needs to be changed

# Data selector
ds = DataSelector(description="Select a single rectangle")
ds.savefiles_var.set(False)
ds.savefiles_check['state'] = 'disabled'
ds.select_directory()

rectangle = ds.directory.split('/')[-1].split(' ')[1]
batch = ds.directory.split('/')[-2]
dat = ds.directory.split('/')[-3]
datadir = '/'.join(ds.directory.split('/')[:-3])
dirname = datadir+u'\\%s\\%s\\Réctangle %s\\'%(dat, batch, rectangle) if dat < '20211212' else datadir+u'\\%s\\%s\\Rectangle %s\\'%(dat, batch, rectangle)

sens = ds.sens_var.get()
if not sens in ('10', '01'):
    raise ValueError('Sens should be either "10" or "01", %s received.'%sens)
AR = ds.AR_var.get()
if AR == 'None':
    AR = None
if not AR in ('A', 'R', None):
    raise ValueError('AR should be either "A", "R" or None, %s received.'%AR)

savefigs = ds.savefigs_var.get()

print('Starting plot fits script.')
print('Batch :', batch, '\nRectangle :', rectangle)

""" Load or create data directory """
filesdir = datadir+u'\\%s\\%s\\Data_files'%(dat, batch)
if not os.path.isdir(dirname+'\\Fits'):
    os.mkdir(dirname+'\\Fits')

time = datetime.datetime.now()
fitsdir = ds.directory+'\\Fits\\'+time.strftime("%d-%m-%Y-%H-%M-%S")

if AR in ('A', 'R'):
    fitsdir += ' - '+AR
if not os.path.isdir(fitsdir):
    os.mkdir(fitsdir)
if not os.path.isdir(fitsdir+'\\Zoom mode 1'):
    os.mkdir(fitsdir+'\\Zoom mode 1')
    os.mkdir(fitsdir+'\\Zoom mode 2')

""" Load data """
prefixe_AR = AR+'_' if AR in ('A', 'R') else ''
sufixe_AR = ' '+AR if AR in ('A', 'R') else ''
offsets = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_offsets.txt'%(batch, rectangle))
f1s = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_f1s.txt'%(batch, rectangle))
gammas1 = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_gammas1.txt'%(batch, rectangle))
ampls1 = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_ampls1.txt'%(batch, rectangle))
f2s = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_f2s.txt'%(batch, rectangle))
gammas2 = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_gammas2.txt'%(batch, rectangle))
ampls2 = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_ampls2.txt'%(batch, rectangle))

popts = [(offsets[ii], f1s[ii]*2*pi, gammas1[ii], ampls1[ii], f2s[ii]*2*pi, gammas2[ii], ampls2[ii]) for ii in range(len(f1s))]

if AR is None:
    spectrogram = loadtxt(dirname+'\\Spectres.txt')
    bitmap = loadtxt(dirname+u'\\Bitmap.txt')
    X_ini = loadtxt(dirname+u'\\Tensions X.txt')
    Y_ini = loadtxt(dirname+u'\\Tensions Y.txt')
else:
    if dat < '20211212':
        spectrogram = loadtxt(dirname+'\\Spectres %s.txt'%aller_retour)
        bitmap = loadtxt(dirname+'\\Bitmap %s.txt'%aller_retour)
        X_ini = loadtxt(dirname+'\\Tensions X %s.txt'%aller_retour)
        Y_ini = loadtxt(dirname+'\\Tensions Y %s.txt'%aller_retour)
    else:
        if AR == 'A':
            spectrogram = loadtxt(dirname+'\\Aller\\Spectres.txt')
            if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                bitmap = loadtxt(dirname+'\\Aller\\Bitmap.txt')
                X_ini = loadtxt(dirname+'\\Aller\\Tensions X.txt')
                Y_ini = loadtxt(dirname+'\\Aller\\Tensions Y.txt')
            elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                bitmap = loadtxt(dirname+'\\Aller\\Bitmap')
                X_ini = loadtxt(dirname+'\\Aller\\Tensions X')
                Y_ini = loadtxt(dirname+'\\Aller\\Tensions Y')
            else:
                raise ValueError('No bitmap file found')
        elif AR == 'R':
            spectrogram = loadtxt(dirname+'\\Retour\\Spectres.txt')
            if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                bitmap = loadtxt(dirname+'\\Retour\\Bitmap.txt')
                X_ini = loadtxt(dirname+'\\Retour\\Tensions X.txt')
                Y_ini = loadtxt(dirname+'\\Retour\\Tensions Y.txt')
            elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                bitmap = loadtxt(dirname+'\\Retour\\Bitmap')
                X_ini = loadtxt(dirname+'\\Retour\\Tensions X')
                Y_ini = loadtxt(dirname+'\\Retour\\Tensions Y')
            else:
                raise ValueError('No bitmap file found')
        else:
            raise ValueError('aller_retour should be either "A" or "R"')
bitmap_ravel = bitmap.ravel()
filename = u'Paramétres spéctres.txt' if dat < '20211212' else u'Parametres spectres.txt'
params = load_params(dirname, filename)
F = get_corrected_freqs(dirname, params, old_corrections=old_position_correction, AR=AR, sens=sens)
I = repeat(range(spectrogram.shape[0]), spectrogram.shape[1]).reshape(spectrogram.shape)
xshape, yshape = [int(float(st)) for st in params[u'Résolution (pixels)'].split('x')] if dat < '20211212' else [int(float(st)) for st in params[u'Resolution (pixels)'].split('x')]

print('Files loaded.')

dict_figs = {}
for ii, spectrum in enumerate(spectrogram):
    if ii%100 == 0:
        print('Spectrum nb %i / %i'%(ii, len(spectrogram)))
    popt = popts[ii]
    ff = F[ii]
    if not isnan(popt[1]) and not isnan(popt[4]):
        fit = array([multi_lorentz(2*pi*f, popt) for f in ff])

        if popt[1] < popt[4]:
            omega1 = popt[1]
            gamma1 = popt[2]
            A1 = popt[3]
            omega2 = popt[4]
            gamma2 = popt[5]
            A2 = popt[6]
        else:
            print(ii, 'inversion des modes')
            omega1 = popt[4]
            gamma1 = popt[5]
            A1 = popt[6]
            omega2 = popt[1]
            gamma2 = popt[2]
            A2 = popt[3]

        fig, ax = subplots()
        plot(ff/1e3, spectrum, label='Data')
        plot(ff/1e3, fit, label='$\Omega_1/2\pi \simeq %i$\n$\Gamma_1/2\pi \simeq %i$\n$A_1 \simeq %.4f$\n$\Omega_2/2\pi \simeq %i$\n$\Gamma_2/2\pi \simeq %i$\n$A_2 \simeq %.4f$'%(omega1/2/pi, gamma1/2/pi, A1, omega2/2/pi, gamma2/2/pi, A2))
        legend()
        xlabel('Frequency (kHz)')
        ylabel('dB')
        title('Spectrum %i'%ii)

        axx = inset_axes(ax, width="30%", height ="30%", loc=2)
        axx.xaxis.set_visible(False)
        axx.yaxis.set_visible(False)
        axx.pcolor(X_ini, Y_ini, bitmap, cmap='gray', shading='nearest')
        axx.invert_yaxis()
        yy, xx = ii//xshape, ii%xshape
        axx.plot(X_ini[yy, xx], Y_ini[yy, xx], 'xr')

        # Zoom mode 1
        fig1, ax1 = subplots()
        plot(ff, spectrum, 'o-', label='Data')
        plot(ff, fit, label='$\Omega_1/2\pi = %s$\n$\Gamma_1/2\pi = %s$\n$A_1 = %s$'%(omega1/2/pi, gamma1/2/pi, A1))
        legend()
        xlabel('Frequency (Hz)')
        ylabel('dB')
        title('Spectrum %i - Mode 1'%ii)
        xlim(omega1/2/pi - max(5*gamma1/2/pi, 500), omega1/2/pi + max(5*gamma1/2/pi, 500))

        axx = inset_axes(ax1, width="30%", height ="30%", loc=2)
        axx.xaxis.set_visible(False)
        axx.yaxis.set_visible(False)
        axx.pcolor(X_ini, Y_ini, bitmap, cmap='gray', shading='nearest')
        axx.invert_yaxis()
        yy, xx = ii//xshape, ii%xshape
        axx.plot(X_ini[yy, xx], Y_ini[yy, xx], 'xr')

        # Zoom mode 2
        fig2, ax2 = subplots()
        plot(ff, spectrum, 'o-', label='Data')
        plot(ff, fit, label='$\Omega_2/2\pi = %s$\n$\Gamma_2/2\pi = %s$\n$A_2 = %s$'%(omega2/2/pi, gamma2/2/pi, A2))
        legend()
        xlabel('Frequency (Hz)')
        ylabel('dB')
        title('Spectrum %i - Mode 2'%ii)
        xlim(omega2/2/pi - max(5*gamma2/2/pi, 500), omega2/2/pi + max(5*gamma2/2/pi, 500))

        axx = inset_axes(ax2, width="30%", height ="30%", loc=2)
        axx.xaxis.set_visible(False)
        axx.yaxis.set_visible(False)
        axx.pcolor(X_ini, Y_ini, bitmap, cmap='gray', shading='nearest')
        axx.invert_yaxis()
        yy, xx = ii//xshape, ii%xshape
        axx.plot(X_ini[yy, xx], Y_ini[yy, xx], 'xr')

        if savefigs:
            fig.savefig(fitsdir+'\\spectrum_'+str(ii).zfill(5))
            fig1.savefig(fitsdir+'\\Zoom mode 1\\spectrum_'+str(ii).zfill(5))
            fig2.savefig(fitsdir+'\\Zoom mode 2\\spectrum_'+str(ii).zfill(5))
        dict_figs['%i'%ii] = fig
        dict_figs['%i_zoom_1'%ii] = fig1
        dict_figs['%i_zoom_2'%ii] = fig2
        close(fig)
        close(fig1)
        close(fig2)
