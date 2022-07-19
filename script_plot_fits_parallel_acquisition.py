# -*- coding: cp1252 -*-
"""
Script for plotting the lorentzian fits of all spectrums measured for a nanowire, with parallel acquisition

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
import parameters_parallel_acquisition as ppa

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

print('Starting plot fits (parallel acquisition) script.')
print('Batch :', batch, '\nRectangle :', rectangle)

""" Load or create figures directory """
time = datetime.datetime.now()
figsdir = datadir+u'\\%s\\%s\\Figures\\Parallel_acquisition\\Fits\\'%(dat, batch)+time.strftime("%d-%m-%Y-%H-%M-%S")

if not os.path.isdir(figsdir):
    os.makedirs(figsdir)
if not os.path.isdir(figsdir+'\\Zoom mode 1'):
    os.mkdir(figsdir+'\\Zoom mode 1')
    os.mkdir(figsdir+'\\Zoom mode 2')

""" Load data """
# Get correspondence with spec nb
for key, corres in ppa.correspondences.items():
    if corres == batch:
        break
else:
    raise ValueError('No correspondence found for batch = %s !'%batch)

dummy1, nw, dummy2, fil, dummy3, spec = key.split('_')

# Spectrogram and frequencies
spectrogram_dir = datadir+u'\\%s\\%s'%(dat, batch)
filename_RSA = 'NW'+nw+'_fil_'+fil+'_spectrogram_'+spec.zfill(4)+'.npy'

spectrogram = load(spectrogram_dir+'\\'+filename_RSA.split('.')[0]+'_without_duplicates'+'.npy')

frequencies_dir = datadir+u'\\%s'%dat
filename_freqs_RSA = 'NW'+nw+'_fil_'+fil+'_frequencies_spectrogram_'+spec.zfill(4)+'.txt'

with open(frequencies_dir + '\\' + filename_freqs_RSA, 'r') as ff:
    lines = ff.readlines()
f_centre = float(lines[0].split(' = ')[1])
span = float(lines[1].split(' = ')[1])
ff = linspace(f_centre-span/2, f_centre+span/2, len(spectrogram[0]))

# Fit results
filesdir = datadir+u'\\%s\\%s\\Data_files\\Parallel_acquisition\\Fits'%(dat, batch)

offsets = loadtxt(filesdir+'\\Fit_offsets.txt')
f1s = loadtxt(filesdir+'\\Fit_f1s.txt')
gammas1 = loadtxt(filesdir+'\\Fit_gammas1.txt')
ampls1 = loadtxt(filesdir+'\\Fit_ampls1.txt')
f2s = loadtxt(filesdir+'\\Fit_f2s.txt')
gammas2 = loadtxt(filesdir+'\\Fit_gammas2.txt')
ampls2 = loadtxt(filesdir+'\\Fit_ampls2.txt')

popts = [(offsets[ii], f1s[ii]*2*pi, gammas1[ii], ampls1[ii], f2s[ii]*2*pi, gammas2[ii], ampls2[ii]) for ii in range(len(f1s))]

print('Files loaded.')

dict_figs = {}
for ii, spectrum in enumerate(spectrogram):
    if ii%100 == 0:
        print('Spectrum nb %i / %i'%(ii, len(spectrogram)))
    popt = popts[ii]
    if not isnan(popt[1]) or not isnan(popt[4]):
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

        if isnan(omega1):
            lab = '$\Omega_2/2\pi \simeq %i$\n$\Gamma_2/2\pi \simeq %i$\n$A_2 \simeq %.4f$'%(omega2/2/pi, gamma2/2/pi, A2)
        elif isnan(omega2):
            lab = '$\Omega_1/2\pi \simeq %i$\n$\Gamma_1/2\pi \simeq %i$\n$A_1 \simeq %.4f$\n$'%(omega1/2/pi, gamma1/2/pi, A1)
        else:
            lab = '$\Omega_1/2\pi \simeq %i$\n$\Gamma_1/2\pi \simeq %i$\n$A_1 \simeq %.4f$\n$\Omega_2/2\pi \simeq %i$\n$\Gamma_2/2\pi \simeq %i$\n$A_2 \simeq %.4f$'%(omega1/2/pi, gamma1/2/pi, A1, omega2/2/pi, gamma2/2/pi, A2)

        fig, ax = subplots()
        plot(ff/1e3, spectrum, label='Data')
        plot(ff/1e3, fit, label=lab)
        legend()
        xlabel('Frequency (kHz)')
        ylabel('dB')
        title('Spectrum %i'%ii)

##        axx = inset_axes(ax, width="30%", height ="30%", loc=2)
##        axx.xaxis.set_visible(False)
##        axx.yaxis.set_visible(False)
##        axx.pcolor(X_ini, Y_ini, bitmap, cmap='gray', shading='nearest')
##        axx.invert_yaxis()
##        yy, xx = ii//xshape, ii%xshape
##        axx.plot(X_ini[yy, xx], Y_ini[yy, xx], 'xr')

        if savefigs:
            fig.savefig(figsdir+'\\spectrum_'+str(ii).zfill(5))
        dict_figs['%i'%ii] = fig
        close(fig)

        # Zoom mode 1
        if not isnan(omega1):
            fig1, ax1 = subplots()
            plot(ff, spectrum, 'o-', label='Data')
            plot(ff, fit, label='$\Omega_1/2\pi = %s$\n$\Gamma_1/2\pi = %s$\n$A_1 = %s$'%(omega1/2/pi, gamma1/2/pi, A1))
            legend()
            xlabel('Frequency (Hz)')
            ylabel('dB')
            title('Spectrum %i - Mode 1'%ii)
            xlim(omega1/2/pi - max(5*gamma1/2/pi, 500), omega1/2/pi + max(5*gamma1/2/pi, 500))

##        axx = inset_axes(ax1, width="30%", height ="30%", loc=2)
##        axx.xaxis.set_visible(False)
##        axx.yaxis.set_visible(False)
##        axx.pcolor(X_ini, Y_ini, bitmap, cmap='gray', shading='nearest')
##        axx.invert_yaxis()
##        yy, xx = ii//xshape, ii%xshape
##        axx.plot(X_ini[yy, xx], Y_ini[yy, xx], 'xr')

        if savefigs:
            fig1.savefig(figsdir+'\\Zoom mode 1\\spectrum_'+str(ii).zfill(5))
        dict_figs['%i_zoom_1'%ii] = fig1
        close(fig1)

        # Zoom mode 2
        if not isnan(omega2):
            fig2, ax2 = subplots()
            plot(ff, spectrum, 'o-', label='Data')
            plot(ff, fit, label='$\Omega_2/2\pi = %s$\n$\Gamma_2/2\pi = %s$\n$A_2 = %s$'%(omega2/2/pi, gamma2/2/pi, A2))
            legend()
            xlabel('Frequency (Hz)')
            ylabel('dB')
            title('Spectrum %i - Mode 2'%ii)
            xlim(omega2/2/pi - max(5*gamma2/2/pi, 500), omega2/2/pi + max(5*gamma2/2/pi, 500))

##        axx = inset_axes(ax2, width="30%", height ="30%", loc=2)
##        axx.xaxis.set_visible(False)
##        axx.yaxis.set_visible(False)
##        axx.pcolor(X_ini, Y_ini, bitmap, cmap='gray', shading='nearest')
##        axx.invert_yaxis()
##        yy, xx = ii//xshape, ii%xshape
##        axx.plot(X_ini[yy, xx], Y_ini[yy, xx], 'xr')

        if savefigs:
            fig2.savefig(figsdir+'\\Zoom mode 2\\spectrum_'+str(ii).zfill(5))
        dict_figs['%i_zoom_2'%ii] = fig2
        close(fig2)
