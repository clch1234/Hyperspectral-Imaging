# -*- coding: cp1252 -*-
"""
Script to take the difference between to directions of
measurement of a single NW, the remove it from the
frequency evolution and fit the result as the
temperature contribution to the frequency evolution.

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
from scipy.interpolate import interp1d

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

dat = '20220622'
batch_01 = 'Fil 6_3'
batch_10 = 'Fil 6_2'
rectangle = '1'
mode = 1 # nb of mechanical mode
kind = 'quadratic' # Type of interpolation (use linear, quadratic or cubic)

savefigs = False
do_fits = False

dirname = "D:\\Documents\\Boulot\\Grenoble\\Data\\%s"%dat

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

""" Differences and sums """
# Reshape freqs
F11_01 = f11_01.reshape(X_01.shape)
F12_01 = f12_01.reshape(X_01.shape)

F11_10 = f11_10.reshape(X_10.shape)
F12_10 = f12_10.reshape(X_10.shape)

y_avg = (Y_01[:, 0] + Y_10[:, 0]) / 2

diff11 = nanmean(F11_10, axis=1)-nanmean(F11_01, axis=1)
diff12 = nanmean(F12_10, axis=1)-nanmean(F12_01, axis=1)

avg11 = (nanmean(F11_10, axis=1)+nanmean(F11_01, axis=1))/2
avg12 = (nanmean(F12_10, axis=1)+nanmean(F12_01, axis=1))/2

""" Interpolation """
y_01 = Y_01.flatten()
# Remove nan values
yd11 = y_avg[where(logical_not(isnan(diff11)))[0]]
d11 = diff11[where(logical_not(isnan(diff11)))[0]]
# Interpolate
f11 = interp1d(yd11, d11, kind=kind, fill_value='extrapolate')

y_10 = Y_10.flatten()
# Remove nan values
yd12 = y_avg[where(logical_not(isnan(diff12)))[0]]
d12 = diff12[where(logical_not(isnan(diff12)))[0]]
# Interpolate
f12 = interp1d(yd12, d12, kind=kind, fill_value='extrapolate')

""" Fit """
def func_to_fit(y, omega0, dTmax):
    if mode == 1:
        return omega0 * (1 + fes.dOmega1_T(y, dTmax))
    elif mode == 2:
        return omega0 * (1 + fes.dOmega2_T(y, dTmax))
vect_to_fit = vectorize(func_to_fit)


# Mode 11_01
f11_01_without_nan = f11_01[where(logical_not(isnan(f11_01)))[0]]
y11_01_without_nan = y_01[where(logical_not(isnan(f11_01)))[0]]
to_fit11_01 = f11_01_without_nan - f11(y11_01_without_nan)
guess = [2*pi*to_fit11_01.max(), 100]
bounds_inf = (to_fit11_01.min(), 0)
bounds_sup = (inf, inf)
print('Fitting mode 11_01')
fit = curve_fit(vect_to_fit,
                y11_01_without_nan,
                to_fit11_01,
                p0=guess,
                bounds=(bounds_inf, bounds_sup))
O11_01, dT11_01 = fit[0]

# Mode 12_01
f12_01_without_nan = f12_01[where(logical_not(isnan(f12_01)))[0]]
y12_01_without_nan = y_01[where(logical_not(isnan(f12_01)))[0]]
to_fit12_01 = f12_01_without_nan - f12(y12_01_without_nan)
guess = [2*pi*to_fit12_01.max(), 100]
bounds_inf = (to_fit12_01.min(), 0)
bounds_sup = (inf, inf)
print('Fitting mode 12_01')
fit = curve_fit(vect_to_fit,
                y12_01_without_nan,
                to_fit12_01,
                p0=guess,
                bounds=(bounds_inf, bounds_sup))
O12_01, dT12_01 = fit[0]

# Mode 11_10
f11_10_without_nan = f11_10[where(logical_not(isnan(f11_10)))[0]]
y11_10_without_nan = y_10[where(logical_not(isnan(f11_10)))[0]]
to_fit11_10 = f11_10_without_nan + f11(y11_10_without_nan)
guess = [2*pi*to_fit11_10.max(), 100]
bounds_inf = (to_fit11_10.min(), 0)
bounds_sup = (inf, inf)
print('Fitting mode 11_10')
fit = curve_fit(vect_to_fit,
                y11_10_without_nan,
                to_fit11_10,
                p0=guess,
                bounds=(bounds_inf, bounds_sup))
O11_10, dT11_10 = fit[0]

# Mode 12_10
f12_10_without_nan = f12_10[where(logical_not(isnan(f12_10)))[0]]
y12_10_without_nan = y_10[where(logical_not(isnan(f12_10)))[0]]
to_fit12_10 = f12_10_without_nan + f12(y12_10_without_nan)
guess = [2*pi*to_fit12_10.max(), 100]
bounds_inf = (to_fit12_10.min(), 0)
bounds_sup = (inf, inf)
print('Fitting mode 12_10')
fit = curve_fit(vect_to_fit,
                y12_10_without_nan,
                to_fit12_10,
                p0=guess,
                bounds=(bounds_inf, bounds_sup))
O12_10, dT12_10 = fit[0]

""" Plot """
yy = linspace(0, 1, 101)

fig_diff, ax_diff = subplots()
plot(y_avg, diff11/1e3, 'o', label='Mode 11')
plot(y_avg, diff12/1e3, 'o', label='Mode 12')
plot(yy,  f11(yy)/1e3, label='Interp 11 (%s)'%kind)
plot(yy,  f12(yy)/1e3, label='Interp 12 (%s)'%kind)
xlabel('Normalized position $y/L$')
ylabel('Frequency (kHz)')
title('Difference "10" - "01"')
legend()
fig_diff.tight_layout()

fig11, ax11 = subplots()
plot(y_01, f11_01/1e3, 'o', mfc='none', label='Raw data')
plot(y11_01_without_nan, to_fit11_01/1e3, 'o', label='To fit 11')
plot(yy, vect_to_fit(yy, O11_01, dT11_01)/1e3, label='Fit:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O11_01/2/pi/1e3, dT11_01))
plot(y_10, f11_10/1e3, 'o', mfc='none', label='Raw data')
plot(y11_10_without_nan, to_fit11_10/1e3, 'o', label='To fit 11')
plot(yy, vect_to_fit(yy, O11_10, dT11_10)/1e3, label='Fit:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O11_10/2/pi/1e3, dT11_10))
xlabel('Normalized position $y/L$')
ylabel('Frequency (kHz)')
title('Mode 11')
legend()
fig11.tight_layout()

fig12, ax12 = subplots()
plot(y_01, f12_01/1e3, 'o', mfc='none', label='Raw data')
plot(y12_01_without_nan, to_fit12_01/1e3, 'o', label='To fit 12')
plot(yy, vect_to_fit(yy, O12_01, dT12_01)/1e3, label='Fit:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O12_01/2/pi/1e3, dT12_01))
plot(y_10, f12_10/1e3, 'o', mfc='none', label='Raw data')
plot(y12_10_without_nan, to_fit12_10/1e3, 'o', label='To fit 12')
plot(yy, vect_to_fit(yy, O12_10, dT12_10)/1e3, label='Fit:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O12_10/2/pi/1e3, dT12_10))
xlabel('Normalized position $y/L$')
ylabel('Frequency (kHz)')
title('Mode 12')
legend()
fig12.tight_layout()

show()
