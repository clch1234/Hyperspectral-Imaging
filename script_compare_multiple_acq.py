# -*- coding: cp1252 -*-
"""
Script for comparing the frequency evolutions for multiple acquisitions on a single NW

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
ds = DataSelector(description="Select a whole acquisition")
ds.savefiles_var.set(False)
ds.savefiles_check['state'] = 'disabled'
ds.sens_label['state'] = 'disabled'
ds.sens_entry['state'] = 'disabled'
ds.AR_label['state'] = 'disabled'
ds.AR_combo['state'] = 'disabled'
ds.select_directory()

batch = ds.directory.split('/')[-1]
dat = ds.directory.split('/')[-2]
workdir = '/'.join(ds.directory.split('/')[:-3])
datadir = '/'.join(ds.directory.split('/')[:-2])

sens = ds.sens_var.get()
if not sens in ('10', '01'):
    raise ValueError('Sens should be either "10" or "01", %s received.'%sens)
aller_retour = ds.AR_var.get()
if aller_retour == 'None':
    aller_retour = None
if not aller_retour in ('A', 'R', None):
    raise ValueError('AR should be either "A", "R" or None, %s received.'%aller_retour)

savefigs = ds.savefigs_var.get()
savefiles = ds.savefiles_var.get()

rectangle = '1'

print('Starting.')

# Works for batches like 'Fil 1_2'
def idx_batch(batch):
    if ' ' in batch:
        return batch.split(' ')[1].split('_')[0]
    else:
        return ''
batches = [dirname for dirname in os.listdir(datadir+'\\'+dat) if idx_batch(batch) == idx_batch(dirname) and os.path.isdir(datadir+'\\%s\\%s\\Data_files'%(dat, dirname))]
print('Detected batches :', batches)

""" Create figures """
figsdir = datadir + '\\%s\\Figures\\'%dat + ' + '.join(batches)
if not os.path.isdir(figsdir):
    os.makedirs(figsdir)

idx_color = 0

fig_freqs, ax_freqs = subplots()
ax_freqs.set_xlabel('y/L')
ax_freqs.set_ylabel('Frequency (kHz)')
ax_freqs.set_title('Frequencies')

fig_splitting, ax_splitting = subplots()
ax_splitting.set_xlabel('y/L')
ax_splitting.set_ylabel('Mode splitting (Hz)')
ax_splitting.set_title('Mode Splitting')

##fig_gammas, ax_gammas = subplots()
##ax_gammas.set_xlabel('y/L')
##ax_gammas.set_ylabel('$\Gamma/2\pi$ (Hz)')
##ax_gammas.set_title('Damping Rates')
##
##fig_ampls, ax_ampls = subplots()
##ax_ampls.set_xlabel('y/L')
##ax_ampls.set_ylabel('Amplitude')
##ax_ampls.set_title('Amplitudes')

for batch in batches:
    print('\n\n')
    print('Batch :', batch, '\nRectangle :', rectangle)

    sens = input('Sens ? ')
    AR = input('AR (y/[n]) ? ')
    if AR in ('n', 'N', ''):
        ARs = (None,)
    elif AR in ('y', 'Y'):
        ARs = ('A', 'R')
    dirname = workdir+u'\\Data\\%s\\%s\\Rectangle %s'%(dat, batch, rectangle) if dat >= '20211212' else workdir+u'\\Data\\%s\\%s\\Réctangle %s'%(dat, batch, rectangle)

    color = 'C%i'%idx_color
    idx_color += 1

    for ii, aller_retour in enumerate(ARs):
        # Flip 'sens' for the second round
        if ii == 1:
            sens = '10' if sens == '01' else '01' if sens == '10' else ''

        print()
        print('Treating aller_retour = %s'%aller_retour)
        print('sens =', sens)

        """ Load or create data files """
        filesdir = datadir+u'\\%s\\%s\\Data_files'%(dat, batch)
        savedir = filesdir+u'\\Fits_%s_%s'%(batch, rectangle)
        if normalize_from_global_image:
            print('Using global image for normalization.')
            savedir += '_Norm_from_global_image'

        prefixe_AR = aller_retour+'_' if aller_retour in ('A', 'R') else ''

        # Load all points
        print()
        print('Loading all points for figure frequency')
        if os.path.isfile(savedir+'\\'+prefixe_AR+batch+'_gammas1.txt'):
            f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = load_fit_data(savedir, batch, prefixe_AR)
        else:
            save_fit_data(datadir, dat, batch, rectangle, sens, remove_aberrant, aller_retour, normalize_from_global_image)
            f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = load_fit_data(savedir, batch, prefixe_AR)

        # Load only points with two peaks
        print()
        print('Loading points with two peaks only for figure splitting')
        savedir_otp = savedir+'_Only_two_peaks\\'
        if os.path.isfile(savedir_otp+prefixe_AR+batch+'_gammas1.txt'):
            f1s_otp, gammas1_otp, ampls1_otp, xs1_otp, ys1_otp, f2s_otp, gammas2_otp, ampls2_otp, xs2_otp, ys2_otp = load_fit_data(savedir_otp, batch, prefixe_AR)
        else:
            save_fit_data(datadir, dat, batch, rectangle, sens, remove_aberrant, aller_retour, normalize_from_global_image, only_two_peaks=True)
            f1s_otp, gammas1_otp, ampls1_otp, xs1_otp, ys1_otp, f2s_otp, gammas2_otp, ampls2_otp, xs2_otp, ys2_otp = load_fit_data(savedir_otp, batch, prefixe_AR)

        print()
        print('Files loaded')

        """ Load the images with normalized coordinates """
        sufixe_AR = ' '+aller_retour if aller_retour in ('A', 'R') else ''
        X_ini = loadtxt(filesdir+'\\'+prefixe_AR+batch+'_rectangle_'+rectangle+'_raw_X.txt')
        Y_ini = loadtxt(filesdir+'\\'+prefixe_AR+batch+'_rectangle_'+rectangle+'_raw_Y.txt')
        if aller_retour is None or dat < '20211212':
            bitmap = loadtxt(dirname+u'\\Bitmap'+sufixe_AR+'.txt')
        elif aller_retour == 'A':
            if dat < '20211212':
                bitmap = loadtxt(dirname+u'\\Bitmap A.txt')
            else:
                if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                    bitmap = loadtxt(dirname+'\\Aller\\Bitmap.txt')
                elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                    bitmap = loadtxt(dirname+'\\Aller\\Bitmap')
                else:
                    raise ValueError('No bitmap file found')
        elif aller_retour == 'R':
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

        marker = '>' if sens == '01' else '<' if sens == '10' else 'o'

        # Figure frequencies
        if len(ys1) > 0:
            ax_freqs.plot(ys1, f1s/1e3, marker=marker, linestyle='', label=batch, c=color)
        if len(ys2) > 0:
            ax_freqs.plot(ys2, f2s/1e3, marker=marker, linestyle='', c=color)
        ax_freqs.legend()

        # Figure splitting
        if len(ys1_otp) > 0:
            ax_splitting.plot(ys1_otp, abs(f1s_otp-f2s_otp), marker=marker, linestyle='', label=batch, c=color)
        ax_splitting.legend()

##    # Figure gamma
##    if len(ys1) > 0:
##        ax_gammas.plot(ys1, gammas1/2/pi, marker=marker, linestyle='', label='%s - Mode 1'%batch)
##    if len(ys2) > 0:
##        ax_gammas.plot(ys2, gammas2/2/pi, marker=marker, linestyle='', label='%s - Mode 2'%batch)
##    ax_gammas.legend()
##
##    # Figure amplitude 
##    if len(ys1) > 0:
##        ax_ampls.semilogy(ys1, ampls1, marker=marker, linestyle='', label='%s - Mode 1'%batch)
##    if len(ys2) > 0:
##        ax_ampls.semilogy(ys2, ampls2, marker=marker, linestyle='', label='%s - Mode 2'%batch)
##    ax_ampls.legend()

if savefigs:
    fig_freqs.savefig(figsdir+'\\frequencies')
    fig_splitting.savefig(figsdir+'\\splitting')
##    fig_gammas.savefig(figsdir+'\\gammas')
##    fig_ampls.savefig(figsdir+'\\amplitudes')
    print('Figures saved')
else:
    print('savefigs set to', savefigs)

show()
