# -*- coding: cp1252 -*-
"""
Script for plotting fit results as 1D plot and as maps on the SEM image
Also saves the local centers, widths, and egdes of the NW.

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
from matplotlib import colors
import os
from data_RSA_new import *
from functions import *
from data_selector import DataSelector

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

# Arguments that still need to be set
remove_aberrant = True # Remove aberrant frequency values
normalize_from_global_image = False # Use global image taken before acquisition for position normalization

# Data selector
ds = DataSelector(description="Select a single rectangle")
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
savefiles = ds.savefiles_var.get()

print('Starting plot data and maps script.')
print('Batch :', batch, '\nRectangle :', rectangle)

""" Load or create data files """
filesdir = datadir+u'\\%s\\%s\\Data_files'%(dat, batch)
figsdir = datadir+u'\\%s\\%s\\Figures\\Normalized_coordinates'%(dat, batch)
savedir = filesdir+u'\\Fits_%s_%s'%(batch, rectangle)
if normalize_from_global_image:
    print('Using global image for normalization.')
    figsdir += '\\Norm_from_global_image'
    savedir += '_Norm_from_global_image'
if not os.path.isdir(figsdir):
    os.mkdir(figsdir)

prefixe_AR = AR+'_' if AR in ('A', 'R') else ''

if os.path.isfile(savedir+'\\'+prefixe_AR+batch+'_gammas1.txt'):
    f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = load_fit_data(savedir, batch, prefixe_AR)
else:
    save_fit_data(datadir, dat, batch, rectangle, sens, remove_aberrant, AR, normalize_from_global_image)
    f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = load_fit_data(savedir, batch, prefixe_AR)
print('Files loaded')

""" Load the images with normalized coordinates """
sufixe_AR = ' '+AR if AR in ('A', 'R') else ''
X_ini = loadtxt(filesdir+'\\'+prefixe_AR+batch+'_rectangle_'+rectangle+'_raw_X.txt')
Y_ini = loadtxt(filesdir+'\\'+prefixe_AR+batch+'_rectangle_'+rectangle+'_raw_Y.txt')
if AR is None or dat < '20211212':
    bitmap = loadtxt(dirname+u'\\Bitmap'+sufixe_AR+'.txt')
elif AR == 'A':
    if dat < '20211212':
        bitmap = loadtxt(dirname+u'\\Bitmap A.txt')
    else:
        if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
            bitmap = loadtxt(dirname+'\\Aller\\Bitmap.txt')
        elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
            bitmap = loadtxt(dirname+'\\Aller\\Bitmap')
        else:
            raise ValueError('No bitmap file found')
elif AR == 'R':
    if dat < '20211212':
        bitmap = loadtxt(dirname+u'\\Bitmap R.txt')
    else:
        if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
            bitmap = loadtxt(dirname+'\\Retour\\Bitmap.txt')
        elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
            bitmap = loadtxt(dirname+'\\Retour\\Bitmap')
        else:
            raise ValueError('No bitmap file found')

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

""" Plot figures """
xs = X_corrected.ravel()
ys = Y_corrected.ravel()

# Figure frequencies
fig_freqs, ax_freqs = subplots()
if len(ys1) > 0:
    plot(ys1, f1s/1e3, 'o', label='Mode 1')
if len(ys2) > 0:
    plot(ys2, f2s/1e3, 'o', label='Mode 2')
xlabel('y/L')
ylabel('Frequency (kHz)')
title('Frequencies')
legend()

# Figure gamma
fig_gammas, ax_gammas = subplots()
if len(ys1) > 0:
    plot(ys1, gammas1/2/pi, 'o', label='Mode 1')
if len(ys2) > 0:
    plot(ys2, gammas2/2/pi, 'o', label='Mode 2')
xlabel('y/L')
ylabel('$\Gamma/2\pi$ (Hz)')
title('Damping Rates')
legend()

# Figure amplitude 
fig_ampls, ax_ampls = subplots()
if len(ys1) > 0:
    semilogy(ys1, ampls1, 'o', label='Mode1')
if len(ys2) > 0:
    semilogy(ys2, ampls2, 'o', label='Mode2')
xlabel('y/L')
ylabel('Amplitude')
title('Amplitudes')
legend()

# Maps
fig_map_f1, ax_map_f1 = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Frequency Mode 1')
suptitle('Rectangle '+rectangle)
if len(ys1) > 0:
    norm = colors.Normalize(vmin=nanmin(f1s)/1e3, vmax=nanmax(f1s)/1e3)
    sc = scatter(xs1, ys1, c=f1s/1e3, norm=norm)
    colorbar(sc, label='Frequency (kHz)')
xlim(nanmin(xs), nanmax(xs))
ylim(nanmin(ys), nanmax(ys))

fig_map_f2, ax_map_f2 = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Frequency Mode 2')
suptitle('Rectangle '+rectangle)
if len(ys2) > 0:
    norm = colors.Normalize(vmin=nanmin(f2s)/1e3, vmax=nanmax(f2s)/1e3)
    sc = scatter(xs2, ys2, c=f2s/1e3, norm=norm)
    colorbar(sc, label='Frequency (kHz)')
xlim(nanmin(xs), nanmax(xs))
ylim(nanmin(ys), nanmax(ys))

fig_map_gammas1, ax_map_gammas1 = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Damping Rate Mode 1')
suptitle('Rectangle '+rectangle)
if len(ys1) > 0:
    norm = colors.Normalize(vmin=nanmin(gammas1)/2/pi, vmax=nanmax(gammas1)/2/pi)
    sc = scatter(xs1, ys1, c=gammas1/2/pi, norm=norm)
    colorbar(sc)
xlim(nanmin(xs), nanmax(xs))
ylim(nanmin(ys), nanmax(ys))

fig_map_gammas2, ax_map_gammas2 = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Damping Rate Mode 2')
suptitle('Rectangle '+rectangle)
if len(ys2) > 0:
    norm = colors.Normalize(vmin=nanmin(gammas2)/2/pi, vmax=nanmax(gammas2)/2/pi)
    sc = scatter(xs2, ys2, c=gammas2/2/pi, norm=norm)
    colorbar(sc)
xlim(nanmin(xs), nanmax(xs))
ylim(nanmin(ys), nanmax(ys))

fig_map_ampls1, ax_map_ampls1 = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Amplitude Mode 1')
suptitle('Rectangle '+rectangle)
if len(ys1) > 0:
    norm = colors.LogNorm(vmin=nanmin(ampls1), vmax=nanmax(ampls1))
    sc = scatter(xs1, ys1, c=ampls1, norm=norm)
    colorbar(sc)
xlim(nanmin(xs), nanmax(xs))
ylim(nanmin(ys), nanmax(ys))

fig_map_ampls2, ax_map_ampls2 = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Amplitude Mode 2')
suptitle('Rectangle '+rectangle)
if len(ys2) > 0:
    norm = colors.LogNorm(vmin=nanmin(ampls2), vmax=nanmax(ampls2))
    sc = scatter(xs2, ys2, c=ampls2, norm=norm)
    colorbar(sc)
xlim(nanmin(xs), nanmax(xs))
ylim(nanmin(ys), nanmax(ys))

if savefigs:
    fig_freqs.savefig(figsdir+'\\frequencies')
    fig_gammas.savefig(figsdir+'\\gammas')
    fig_ampls.savefig(figsdir+'\\amplitudes')

    fig_map_f1.savefig(figsdir+'\\map_frequencies1')
    fig_map_f2.savefig(figsdir+'\\map_frequencies2')
    fig_map_gammas1.savefig(figsdir+'\\map_gammas1')
    fig_map_gammas2.savefig(figsdir+'\\map_gammas2')
    fig_map_ampls1.savefig(figsdir+'\\map_amplitudes1')
    fig_map_ampls2.savefig(figsdir+'\\map_amplitudes2')
    print('Figures saved')
else:
    print('savefigs set to', savefigs)

if savefiles:
    savetxt(filesdir+'\\%s_rectangle_%s_Centers.txt'%(batch, rectangle), centers)
    savetxt(filesdir+'\\%s_rectangle_%s_Widths.txt'%(batch, rectangle), widths)
    savetxt(filesdir+'\\%s_rectangle_%s_Edges.txt'%(batch, rectangle), edges)

show()
