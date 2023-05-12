# -*- coding: cp1252 -*-
"""
Script for comparing the frequency evolutions for multiple acquisitions on a single NW
Here with the 2nd mechanical mode measured in parallel

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
from data_selector import DataSelector, BatchesSelector
from PyQt5 import QtWidgets, QtGui, QtCore
import sys
import parameters_parallel_acquisition as ppa

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

dat = '20220622'
batch = 'Fil 6'
rectangle = '1'
inversion = '0016'

savefigs = False

dirname = "D:\\Documents\\Boulot\\Grenoble\\Data\\%s"%dat

""" Create figures """
figsdir = dirname + '\\Figures\\%s multiple acquisitions'%batch
if not os.path.isdir(figsdir):
    os.makedirs(figsdir)

idx_color = 0

fig11, ax11 = subplots()
ax11.set_xlabel('y/L')
ax11.set_ylabel('Frequency (kHz)')
ax11.set_title('1st mechanical mode, 1st peak')

fig12, ax12 = subplots()
ax12.set_xlabel('y/L')
ax12.set_ylabel('Frequency (kHz)')
ax12.set_title('1st mechanical mode, 2nd peak')

fig21, ax21 = subplots()
ax21.set_xlabel('y/L')
ax21.set_ylabel('Frequency (kHz)')
ax21.set_title('2nd mechanical mode, 1st peak')

fig22, ax22 = subplots()
ax22.set_xlabel('y/L')
ax22.set_ylabel('Frequency (kHz)')
ax22.set_title('2nd mechanical mode, 2nd peak')

fig_chrono, ax = subplots()
ax.set_xlabel('y/L')
ax.set_ylabel('Frequency (kHz)')
ax.set_title('Mode 1.1 in chronological order')

""" Load and plot """
for key, corres in ppa.correspondences.items():
    dummy1, nw, dummy2, fil, dummy3, spec = key.split('_')

    if batch in corres:
        dir_coordinates = dirname + u'\\%s\\Data_files'%corres
        dir_bitmap = dirname + u'\\%s\\Rectangle %s'%(corres, rectangle)

        if spec < inversion:
            # Before we switch the RSAs : mode 1 on ViBr, mode 2 on continuous acquisition
            dir_mode2 = dirname + u'\\%s\\Data_files\\Parallel_acquisition'%corres
            dir_mode1 = dirname + u'\\%s\\Data_files\\Fits_%s_%s'%(corres, corres, rectangle)

            """ Load mode 2 """
            f21s = load(dir_mode2+'\\freqs1.npy', allow_pickle=True)
            ys21 = load(dir_mode2+'\\ys1.npy', allow_pickle=True)

            f22s = load(dir_mode2+'\\freqs2.npy', allow_pickle=True)
            ys22 = load(dir_mode2+'\\ys2.npy', allow_pickle=True)


            # Normalize y by length
            # Load X, Y and bitmap
            bitmap = loadtxt(dir_bitmap+u'\\Bitmap.txt')
            X_ini = loadtxt(dir_coordinates+u'\\%s_rectangle_1_raw_X.txt'%corres)
            Y_ini = loadtxt(dir_coordinates+u'\\%s_rectangle_1_raw_Y.txt'%corres)

            # Detection of the edges of the NW and left-right separation
            centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(bitmap, X_ini)

            # Deduce basis position and NW length
            length = abs(Y_ini[first_line_found, 0] - Y_ini[last_line_found, 0])
            basis = Y_ini[first_line_found, 0]

            ys21 = (ys21 - basis) / length
            ys22 = (ys22 - basis) / length

            # Plot
            ax21.plot(ys21, f21s/1e3, 'o', label=corres, c='C%i'%idx_color)
            ax22.plot(ys22, f22s/1e3, 'o', label=corres, c='C%i'%idx_color)

            """ Load mode 1 """
            f11s, gammas11, ampls11, xs11, ys11, f12s, gammas12, ampls12, xs12, ys12 = load_fit_data(dir_mode1, corres)
            
            # Plot
            ax11.plot(ys11, f11s/1e3, 'o', label=corres, c='C%i'%idx_color)
            ax12.plot(ys12, f12s/1e3, 'o', label=corres, c='C%i'%idx_color)

            """ Test """
            yy = int(corres.split('_')[-1]) + ys11 if ppa.dd_sens[key] == '01' else int(corres.split('_')[-1]) + (1-ys11)
            ax.plot(yy, f11s/1e3, 'o', label=corres, c='C%i'%idx_color)

        else:
            # After we switch the RSAs : mode 1 on continuous acquisition, mode 2 on ViBr
            dir_mode1 = dirname + u'\\%s\\Data_files\\Parallel_acquisition'%corres
            dir_mode2 = dirname + u'\\%s\\Data_files\\Fits_%s_%s'%(corres, corres, rectangle)

            """ Load mode 1 """
            f11s = load(dir_mode1+'\\freqs1.npy', allow_pickle=True)
            ys11 = load(dir_mode1+'\\ys1.npy', allow_pickle=True)

            f12s = load(dir_mode1+'\\freqs2.npy', allow_pickle=True)
            ys12 = load(dir_mode1+'\\ys2.npy', allow_pickle=True)

            # Normalize y by length
            # Load X, Y and bitmap
            bitmap = loadtxt(dir_bitmap+u'\\Bitmap.txt')
            X_ini = loadtxt(dir_coordinates+u'\\%s_rectangle_1_raw_X.txt'%corres)
            Y_ini = loadtxt(dir_coordinates+u'\\%s_rectangle_1_raw_Y.txt'%corres)

            # Detection of the edges of the NW and left-right separation
            centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(bitmap, X_ini)

            # Deduce basis position and NW length
            length = abs(Y_ini[first_line_found, 0] - Y_ini[last_line_found, 0])
            basis = Y_ini[first_line_found, 0]

            ys11 = (ys11 - basis) / length
            ys12 = (ys12 - basis) / length

            # Plot
            ax11.plot(ys11, f11s/1e3, 'o', label=corres, c='C%i'%idx_color)
            ax12.plot(ys12, f12s/1e3, 'o', label=corres, c='C%i'%idx_color)

            """ Test """
            yy = int(corres.split('_')[-1]) + ys11 if ppa.dd_sens[key] == '01' else int(corres.split('_')[-1]) + (1-ys11)
            ax.plot(yy, f11s/1e3, 'o', label=corres, c='C%i'%idx_color)

            """ Load mode 2 """
            f21s, gammas21, ampls21, xs21, ys21, f22s, gammas22, ampls22, xs22, ys22 = load_fit_data(dir_mode2, corres)
            
            # Plot
            ax21.plot(ys21, f21s/1e3, 'o', label=corres, c='C%i'%idx_color)
            ax22.plot(ys22, f22s/1e3, 'o', label=corres, c='C%i'%idx_color)

        idx_color += 1

ax11.legend()
ax12.legend()
ax21.legend()
ax22.legend()
fig11.tight_layout()
fig12.tight_layout()
fig21.tight_layout()
fig22.tight_layout()
fig_chrono.tight_layout()

show()

if savefigs:
    fig11.savefig(figsdir+'\\Mode_11')
    fig12.savefig(figsdir+'\\Mode_12')
    fig21.savefig(figsdir+'\\Mode_21')
    fig22.savefig(figsdir+'\\Mode_22')
    fig_chrono.savefig(figsdir+'\\Mode_11_chrono')
    print('Figures saved')
else:
    print('savefigs set to', savefigs)

show()
