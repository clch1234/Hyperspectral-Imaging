# -*- coding: cp1252 -*-
"""
Script for comparing the frequency evolutions for multiple acquisitions on a single NW
Here only for the mode measured in ViBr

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
import functions_engraving_study as fes
from scipy.optimize import curve_fit

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

dat = '20211216'
batch_01 = 'Fil 19_2'
batch_10 = 'Fil 19'
first_direction = '10' # Only useful if batch_10 == batch_01, ignored otherwise
rectangle = '1'
mode = 1 # nb of mechanical mode

savefigs = False
do_fits = True

dirname = "D:\\Documents\\Boulot\\Grenoble\\Data\\%s"%dat

""" Create figures """
if savefigs:
    figsdir = dirname + '\\Figures\\Sum_diff_%s_%s'%(batch_01, batch_10)
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

fig_diff1, ax_diff1 = subplots()
ax_diff1.set_xlabel('y/L')
ax_diff1.set_ylabel('Frequency difference (kHz) Mode 1')
ax_diff1.set_title('Difference "10" - "01"')

fig_sum1, ax_sum1 = subplots()
ax_sum1.set_xlabel('y/L')
ax_sum1.set_ylabel('Frequency difference (kHz)')
ax_sum1.set_title('Average ("10" + "01") / 2 Mode 1')

""" Load mode 1 """
if batch_10 == batch_01:
    if first_direction == '10':
        prefixe_10 = 'A_'
        prefixe_01 = 'R_'
        suffixe_10 = '' if dat >= '20211212' else ' A'
        suffixe_01 = '' if dat >= '20211212' else ' R'
    elif first_direction == '01':
        prefixe_10 = 'R_'
        prefixe_01 = 'A_'
        suffixe_10 = '' if dat >= '20211212' else ' R'
        suffixe_01 = '' if dat >= '20211212' else ' A'
else:
    prefixe_10 = ''
    prefixe_01 = ''
    suffixe_10 = ''
    suffixe_01 = ''

with open(dirname+'\\%s\\Data_files\\%s_rectangle_1_Fit_f1s.txt'%(batch_01, prefixe_01+batch_01), 'r') as ff:
    f11_01 = loadtxt(ff)
with open(dirname+'\\%s\\Data_files\\%s_rectangle_1_Fit_f2s.txt'%(batch_01, prefixe_01+batch_01), 'r') as ff:
    f12_01 = loadtxt(ff)
with open(dirname+'\\%s\\Data_files\\%s_rectangle_1_raw_X.txt'%(batch_01, prefixe_01+batch_01), 'r') as ff:
    X_01 = loadtxt(ff)
with open(dirname+'\\%s\\Data_files\\%s_rectangle_1_raw_Y.txt'%(batch_01, prefixe_01+batch_01), 'r') as ff:
    Y_01 = loadtxt(ff)

with open(dirname+'\\%s\\Data_files\\%s_rectangle_1_Fit_f1s.txt'%(batch_10, prefixe_10+batch_10), 'r') as ff:
    f11_10 = loadtxt(ff)
with open(dirname+'\\%s\\Data_files\\%s_rectangle_1_Fit_f2s.txt'%(batch_10, prefixe_10+batch_10), 'r') as ff:
    f12_10 = loadtxt(ff)
with open(dirname+'\\%s\\Data_files\\%s_rectangle_1_raw_X.txt'%(batch_10, prefixe_10+batch_10), 'r') as ff:
    X_10 = loadtxt(ff)
with open(dirname+'\\%s\\Data_files\\%s_rectangle_1_raw_Y.txt'%(batch_10, prefixe_10+batch_10), 'r') as ff:
    Y_10 = loadtxt(ff)

# Correct inversions
correct_inversions(f11_01, f12_01)
correct_inversions(f11_10, f12_10)

# Normalize y by length
# Load X, Y and bitmap
dir_coordinates_01 = dirname + u'\\%s\\Data_files'%batch_01
dir_coordinates_10 = dirname + u'\\%s\\Data_files'%batch_10
if batch_01 == batch_10 and dat >= '20211212':
    if first_direction == '01':
        dir_bitmap_01 = dirname + u'\\%s\\Rectangle %s\\Aller'%(batch_01, rectangle)
        dir_bitmap_10 = dirname + u'\\%s\\Rectangle %s\\Retour'%(batch_10, rectangle)
    elif first_direction == '10':
        dir_bitmap_01 = dirname + u'\\%s\\Rectangle %s\\Retour'%(batch_01, rectangle)
        dir_bitmap_10 = dirname + u'\\%s\\Rectangle %s\\Aller'%(batch_10, rectangle)
else:
    dir_bitmap_01 = dirname + u'\\%s\\Rectangle %s'%(batch_01, rectangle) if dat >= '20211212' else dirname + u'\\%s\\Réctangle %s'%(batch_01, rectangle)
    dir_bitmap_10 = dirname + u'\\%s\\Rectangle %s'%(batch_10, rectangle) if dat >= '20211212' else dirname + u'\\%s\\Réctangle %s'%(batch_10, rectangle)

bitmap_01 = loadtxt(dir_bitmap_01+u'\\Bitmap%s.txt'%suffixe_01)
X_ini_01 = loadtxt(dir_coordinates_01+u'\\%s_rectangle_1_raw_X.txt'%(prefixe_01+batch_01))
Y_ini_01 = loadtxt(dir_coordinates_01+u'\\%s_rectangle_1_raw_Y.txt'%(prefixe_01+batch_01))

bitmap_10 = loadtxt(dir_bitmap_10+u'\\Bitmap%s.txt'%suffixe_10)
X_ini_10 = loadtxt(dir_coordinates_10+u'\\%s_rectangle_1_raw_X.txt'%(prefixe_10+batch_10))
Y_ini_10 = loadtxt(dir_coordinates_10+u'\\%s_rectangle_1_raw_Y.txt'%(prefixe_10+batch_10))

# Normalize Y by length and X by width
centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(bitmap_01, X_01)
X_01 = get_normalized_X(X_01, centers, widths, first_line_found, last_line_found)
length = abs(Y_01[first_line_found, 0] - Y_01[last_line_found, 0])
basis = Y_01[last_line_found, 0]
Y_01 = (basis - Y_01) / length

centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(bitmap_10, X_10)
X_10 = get_normalized_X(X_10, centers, widths, first_line_found, last_line_found)
length = abs(Y_10[first_line_found, 0] - Y_10[last_line_found, 0])
basis = Y_10[last_line_found, 0]
Y_10 = (basis - Y_10) / length

""" Remove aberrant values """
key = 'rectangle_%s'%rectangle

class FakeParams(object):
    """ Used to emulated empty parameters """
    def __init__(self):
        super(FakeParams, self).__init__()
        self.dd_min_freq1 = {}
        self.dd_max_freq1 = {}
        self.dd_min_freq2 = {}
        self.dd_max_freq2 = {}

# For direction 01
parameters_file_path = dirname+u'\\%s\\parameters.py'%batch_01
if os.path.isfile(parameters_file_path):
    spec = importlib.util.spec_from_file_location('parameters', parameters_file_path)
    parameters = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parameters)
else:
    print('No parameters file found for direction 01, emulating empty.')
    parameters = FakeParams()

min_freq1 = parameters.dd_min_freq1[key] if key in parameters.dd_min_freq1 else nanmin(f11_01)-1
max_freq1 = parameters.dd_max_freq1[key] if key in parameters.dd_max_freq1 else nanmax(f11_01)+1

min_freq2 = parameters.dd_min_freq2[key] if key in parameters.dd_min_freq2 else nanmin(f12_01)-1
max_freq2 = parameters.dd_max_freq2[key] if key in parameters.dd_max_freq2 else nanmax(f12_01)+1

to_remove1 = where((f11_01 < min_freq1) | (f11_01 > max_freq1))[0]
to_remove2 = where((f12_01 < min_freq2) | (f12_01 > max_freq2))[0]

f11_01[to_remove1] = nan
f12_01[to_remove1] = nan

# For direction 10
parameters_file_path = dirname+u'\\%s\\parameters.py'%batch_10
if os.path.isfile(parameters_file_path):
    spec = importlib.util.spec_from_file_location('parameters', parameters_file_path)
    parameters = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parameters)
else:
    print('No parameters file found for direction 10, emulating empty.')
    parameters = FakeParams()

min_freq1 = parameters.dd_min_freq1[key] if key in parameters.dd_min_freq1 else nanmin(f11_10)-1
max_freq1 = parameters.dd_max_freq1[key] if key in parameters.dd_max_freq1 else nanmax(f11_10)+1

min_freq2 = parameters.dd_min_freq2[key] if key in parameters.dd_min_freq2 else nanmin(f12_10)-1
max_freq2 = parameters.dd_max_freq2[key] if key in parameters.dd_max_freq2 else nanmax(f12_10)+1

to_remove1 = where((f11_10 < min_freq1) | (f11_10 > max_freq1))[0]
to_remove2 = where((f12_10 < min_freq2) | (f12_10 > max_freq2))[0]

f11_10[to_remove1] = nan
f12_10[to_remove1] = nan

# Reshape freqs
F11_01 = f11_01.reshape(X_01.shape)
F12_01 = f12_01.reshape(X_01.shape)

F11_10 = f11_10.reshape(X_10.shape)
F12_10 = f12_10.reshape(X_10.shape)

# Plot raw data
ax11.plot(Y_01.flatten(), f11_01/1e3, 'o', label=batch_01)
ax12.plot(Y_01.flatten(), f12_01/1e3, 'o', label=batch_01)
ax11.plot(Y_10.flatten(), f11_10/1e3, 'o', label=batch_10)
ax12.plot(Y_10.flatten(), f12_10/1e3, 'o', label=batch_10)

# Plot mean: points kept for sum and difference
ax11.plot(Y_01[:, 0], nanmean(F11_01, axis=1)/1e3, 'sk', mfc='none')
ax12.plot(Y_01[:, 0], nanmean(F12_01, axis=1)/1e3, 'sk', mfc='none')
ax11.plot(Y_10[:, 0], nanmean(F11_10, axis=1)/1e3, 'sk', mfc='none')
ax12.plot(Y_10[:, 0], nanmean(F12_10, axis=1)/1e3, 'sk', mfc='none')

idx_color += 1

""" Harmonize sizes """
# For mode 1
if dat == '20220622' and batch_01 == 'Fil 6_3' and batch_10 == 'Fil 6_4':
    start1_01 = 0
    stop1_01 = F11_01.shape[1]
    start1_10 = 0
    stop1_10 = -1

elif dat == '20220622' and batch_01 == 'Fil 6_1' and batch_10 == 'Fil 6_2':
    start1_01 = 0
    stop1_01 = F11_01.shape[1]
    start1_10 = 0
    stop1_10 = F11_10.shape[1]

elif dat == '20220622' and batch_01 == 'Fil 6_3' and batch_10 == 'Fil 6_2':
    start1_01 = 0
    stop1_01 = F11_01.shape[1]
    start1_10 = 0
    stop1_10 = F11_10.shape[1]

elif dat == '20220622' and batch_01 == 'Fil 5_6' and batch_10 == 'Fil 5_5':
    start1_01 = 1
    stop1_01 = -1
    start1_10 = 0
    stop1_10 = F11_10.shape[1]

elif dat == '20211216' and batch_01 == 'Fil 19_2' and batch_10 == 'Fil 19':
    start1_01 = 4
    stop1_01 = F11_01.shape[1]
    start1_10 = 0
    stop1_10 = -4

else:
    start1_01 = 0
    stop1_01 = F11_01.shape[1]
    start1_10 = 0
    stop1_10 = F11_10.shape[1]


""" Differences and sums """
y_avg = (Y_01[start1_01:stop1_01, 0] + Y_10[start1_10:stop1_10, 0]) / 2

ax_diff1.plot(y_avg, (nanmean(F11_10, axis=1)[start1_10:stop1_10]-nanmean(F11_01, axis=1)[start1_01:stop1_01])/1e3, 'o', label='Peak 1')
ax_diff1.plot(y_avg, (nanmean(F12_10, axis=1)[start1_10:stop1_10]-nanmean(F12_01, axis=1)[start1_01:stop1_01])/1e3, 'o', label='Peak 2')

avg11 = (nanmean(F11_10, axis=1)[start1_10:stop1_10]+nanmean(F11_01, axis=1)[start1_01:stop1_01])/2
avg12 = (nanmean(F12_10, axis=1)[start1_10:stop1_10]+nanmean(F12_01, axis=1)[start1_01:stop1_01])/2

ax_sum1.plot(y_avg, avg11/1e3, 'o', label='Peak 1')
ax_sum1.plot(y_avg, avg12/1e3, 'o', label='Peak 2')

ax11.plot(y_avg, avg11/1e3, 'o', label='Avg')
ax12.plot(y_avg, avg12/1e3, 'o', label='Avg')

""" Fit avg with dOmega_T """
if do_fits:
    def func_to_fit(y, omega0, dTmax):
        if mode == 1:
            return omega0 * (1 + fes.dOmega1_T(y, dTmax))
        elif mode == 2:
            return omega0 * (1 + fes.dOmega2_T(y, dTmax))
    vect_to_fit = vectorize(func_to_fit)

    # 1st Peak
    to_fit11 = 2*pi*avg11[where(logical_not(isnan(avg11)))[0]]
    y11 = Y_01[where(logical_not(isnan(avg11)))[0], 0]
    if len(to_fit11) > 0:
        guess = [2*pi*to_fit11.max(), 100]
        bounds_inf = (to_fit11.min(), 0)
        bounds_sup = (inf, inf)

        fit = curve_fit(vect_to_fit,
                        y11,
                        to_fit11,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O11, dT11 = fit[0]


        yy = linspace(0, 1, 101)
        ax_sum1.plot(yy, vect_to_fit(yy, O11, dT11)/2/pi/1e3, label='Fit peak 1:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O11/2/pi/1e3, dT11))
        ax11.plot(yy, vect_to_fit(yy, O11, dT11)/2/pi/1e3, label='Fit:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O11/2/pi/1e3, dT11))

    # 2nd Peak
    to_fit12 = 2*pi*avg12[where(logical_not(isnan(avg12)))[0]]
    y12 = Y_01[where(logical_not(isnan(avg12)))[0], 0]
    if len(to_fit12) > 0:
        guess = [2*pi*to_fit12.max(), 100]
        bounds_inf = (to_fit12.min(), 0)
        bounds_sup = (inf, inf)

        fit = curve_fit(vect_to_fit,
                        y12,
                        to_fit12,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O12, dT12 = fit[0]


        yy = linspace(0, 1, 101)
        ax_sum1.plot(yy, vect_to_fit(yy, O12, dT12)/2/pi/1e3, label='Fit peak 2:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O12/2/pi/1e3, dT12))
        ax12.plot(yy, vect_to_fit(yy, O12, dT12)/2/pi/1e3, label='Fit:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O12/2/pi/1e3, dT12))


ax11.legend()
ax12.legend()
fig11.tight_layout()
fig12.tight_layout()

ax_diff1.legend()
fig_diff1.tight_layout()
ax_sum1.legend()
fig_sum1.tight_layout()

show()

if savefigs:
    fig11.savefig(figsdir+'\\Mode_11')
    fig12.savefig(figsdir+'\\Mode_12')

    fig_diff1.savefig(figsdir+'\\Difference_1')
    fig_sum1.savefig(figsdir+'\\Sum_1')
    print('Figures saved')
else:
    print('savefigs set to', savefigs)

show()
