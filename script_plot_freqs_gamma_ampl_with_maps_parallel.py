# -*- coding: cp1252 -*-
"""
Script for plotting fit results as 1D plot and (soon) as maps on the SEM image
This one is for the mode acquired in parallel (with the old acquisition method)

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
import parameters_parallel_acquisition as ppa

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

""" Load data files """
for key, corres in ppa.correspondences.items():
    if corres == batch:
        break
else:
    raise ValueError('No correspondence found for batch = %s !'%batch)

dirname = '/'.join(ds.directory.split('/')[:-1])
filesdir = datadir + u'\\%s\\%s\\Data_files\\Parallel_acquisition'%(dat, batch)
figsdir = datadir + u'\\%s\\%s\\Figures\\Normalized_coordinates\\Parallel_acquisition'%(dat, batch)

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

X_corrected = get_normalized_X(X_ini, centers, widths, first_line_found, last_line_found)
Y_corrected = (Y_ini - basis) / length

""" Plot figures """
xs = X_corrected.ravel()
ys = Y_corrected.ravel()

# Figure frequencies
fig_freqs, ax_freqs = subplots()
if len(ys1) > 0:
    plot(ys1, f1s/1e3, 'o', label='Mode 2.1')
if len(ys2) > 0:
    plot(ys2, f2s/1e3, 'o', label='Mode 2.2')
xlabel('y/L')
ylabel('Frequency (kHz)')
title('Frequencies')
legend()
tight_layout()

if len(ys1) > 0:
    fig_freqs1, ax_freqs1 = subplots()
    plot(ys1, f1s/1e3, 'o', label='Mode 2.1')
    xlabel('y/L')
    ylabel('Frequency (kHz)')
    title('Frequencies Mode 2.1')
    legend()
    tight_layout()

if len(ys2) > 0:
    fig_freqs2, ax_freqs2 = subplots()
    plot(ys2, f2s/1e3, 'o', label='Mode 2.2')
    xlabel('y/L')
    ylabel('Frequency (kHz)')
    title('Frequencies Mode 2.2')
    legend()
    tight_layout()
"""
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
xlim(max(-5, nanmin(xs[xs != -inf])), min(5, nanmax(xs[xs != inf])))
ylim(nanmin(ys[ys != -inf]), nanmax(ys[ys != inf]))

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
xlim(max(-5, nanmin(xs[xs != -inf])), min(5, nanmax(xs[xs != inf])))
ylim(nanmin(ys[ys != -inf]), nanmax(ys[ys != inf]))

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
xlim(max(-5, nanmin(xs[xs != -inf])), min(5, nanmax(xs[xs != inf])))
ylim(nanmin(ys[ys != -inf]), nanmax(ys[ys != inf]))

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
xlim(max(-5, nanmin(xs[xs != -inf])), min(5, nanmax(xs[xs != inf])))
ylim(nanmin(ys[ys != -inf]), nanmax(ys[ys != inf]))

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
xlim(max(-5, nanmin(xs[xs != -inf])), min(5, nanmax(xs[xs != inf])))
ylim(nanmin(ys[ys != -inf]), nanmax(ys[ys != inf]))

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
xlim(max(-5, nanmin(xs[xs != -inf])), min(5, nanmax(xs[xs != inf])))
ylim(nanmin(ys[ys != -inf]), nanmax(ys[ys != inf]))
"""
if savefigs:
    fig_freqs.savefig(figsdir+'\\frequencies')
    if len(ys1) > 0:
        fig_freqs1.savefig(figsdir+'\\frequencies1')
    if len(ys2) > 0:
        fig_freqs2.savefig(figsdir+'\\frequencies2')
##    fig_gammas.savefig(figsdir+'\\gammas')
##    fig_ampls.savefig(figsdir+'\\amplitudes')

##    fig_map_f1.savefig(figsdir+'\\map_frequencies1')
##    fig_map_f2.savefig(figsdir+'\\map_frequencies2')
##    fig_map_gammas1.savefig(figsdir+'\\map_gammas1')
##    fig_map_gammas2.savefig(figsdir+'\\map_gammas2')
##    fig_map_ampls1.savefig(figsdir+'\\map_amplitudes1')
##    fig_map_ampls2.savefig(figsdir+'\\map_amplitudes2')
    print('Figures saved')
else:
    print('savefigs set to', savefigs)

show()
