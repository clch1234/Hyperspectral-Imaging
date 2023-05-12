# -*- coding: cp1252 -*-
"""
Script for plotting the following:
    - Bitmap image of acquisition
    - Gradient along X-axis of bitmap
    - Map of peak amplitude

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
import functions_engraving_study as fes
from scipy.optimize import curve_fit

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

# Arguments that still need to be set
remove_aberrant = True # Remove aberrant frequency values
normalize_from_global_image = False # Use global image taken before acquisition for position normalization
old_position_correction = False # Normally never needs to be changed

""" Data selector """
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

""" Load data """
filesdir = datadir+u'\\%s\\%s\\Data_files'%(dat, batch)
savedir = filesdir+u'\\Fits_%s_%s'%(batch, rectangle)

# Load fits
prefixe_AR = AR+'_' if AR in ('A', 'R') else ''
sufixe_AR = ' '+AR if AR in ('A', 'R') else ''
if os.path.isfile(savedir+'\\'+prefixe_AR+batch+'_gammas1.txt'):
    f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = load_fit_data(savedir, batch, prefixe_AR)
else:
    save_fit_data(datadir, dat, batch, rectangle, sens, remove_aberrant, AR, normalize_from_global_image)
    f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = load_fit_data(savedir, batch, prefixe_AR)

# Load bitmap
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

centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(bitmap, X_ini)
X_corrected = get_normalized_X(X_ini, centers, widths, first_line_found, last_line_found)

if normalize_from_global_image:
    print('Normalizing from global image')
    Y_corrected = get_normalized_Y(datadir, dat, batch, rectangle, sens)#, color_box='r')
else:
    print('Normalizing from acquisition image')
    length = abs(Y_ini[first_line_found, 0] - Y_ini[last_line_found, 0])
    basis = Y_ini[last_line_found, 0]
    Y_corrected = (basis - Y_ini) / length

print('Files loaded.')

""" Calculate gradient """
grad = gradient(bitmap, axis=1)
print('Gradient calculated')

""" Plot """
fig_bitmap, ax_bitmap = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Raw image of measurement')
xlim(-2, 2)
ax_bitmap.set_aspect(10)
tight_layout()

fig_grad, ax_grad = subplots()
pcolor(X_corrected, Y_corrected, grad, shading='nearest') #cmap='gray', 
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Gradient of the raw image image of measurement')
colorbar()
xlim(-2, 2)
ax_grad.set_aspect(10)
tight_layout()

## This is a log plot of abs(grad) ##
fig_abs_grad, ax_abs_grad = subplots()
# Remove 0 value from grad, replace by next smallest (abs) value
grad[where(grad==0)] = nan
m = nanmin(abs(grad))
grad[where(isnan(grad))] = m
# Normalize for log plot
norm = colors.LogNorm(vmin=nanmin(abs(grad)), vmax=nanmax(abs(grad)))
# Plot
pcolor(X_corrected, Y_corrected, abs(grad), norm=norm, shading='nearest') #cmap='gray', 
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Gradient of the raw image image of measurement')
colorbar(label='Absolute value of gradient (arb.)')
xlim(-2, 2)
ax_abs_grad.set_aspect(10)
tight_layout()

fig_ampl1, ax_ampl1 = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
norm = colors.LogNorm(vmin=nanmin(ampls1), vmax=nanmax(ampls1))
sc = scatter(xs1, ys1, c=ampls1, norm=norm)
colorbar(sc, label='Amplitude (arb.)')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Amplitude Mode 1')
xlim(-2, 2)
ax_ampl1.set_aspect(10)
tight_layout()

fig_ampl2, ax_ampl2 = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
norm = colors.LogNorm(vmin=nanmin(ampls2), vmax=nanmax(ampls2))
sc = scatter(xs2, ys2, c=ampls2, norm=norm)
colorbar(sc, label='Amplitude (arb.)')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Amplitude Mode 2')
xlim(-2, 2)
ax_ampl2.set_aspect(10)
tight_layout()

show()

""" Save """
if savefigs:
    figsdir = datadir+u'\\%s\\%s\\Figures\\Bitmap_Gradient_Amplitude'%(dat, batch)
    if not os.path.isdir(figsdir):
        os.mkdir(figsdir)

    fig_bitmap.savefig(figsdir+'\\Bitmap.png')
    fig_grad.savefig(figsdir+'\\Gradient.png')
    fig_abs_grad.savefig(figsdir+'\\Gradient_Abs.png')
    fig_ampl1.savefig(figsdir+'\\Amplitude1.png')
    fig_ampl2.savefig(figsdir+'\\Amplitude2.png')
    print('Figures saved')
