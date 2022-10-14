# -*- coding: cp1252 -*-
"""
Script for calculating average and difference of frequency during a round-trip of measurement,
fit average with temperature contribution on frequency shift
then with temperature + back-action contributions.
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
from scipy.integrate import quad

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

dat = '20220622'
batch_01 = 'Fil 6_3'
batch_10 = 'Fil 6_2'
first_direction = '10' # Only useful if batch_10 == batch_01, ignored otherwise
rectangle = '1'
mode = 1 # nb of mechanical mode

savefigs = False
do_fits = True

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

class EmptyParams(object):
    """ Used to emulated empty parameters """
    def __init__(self):
        super(EmptyParams, self).__init__()
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
    parameters = EmptyParams()

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
    parameters = EmptyParams()

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

""" Calculate averages for each position """
y_01 = Y_01[start1_01:stop1_01, 0]
y_10 = Y_10[start1_10:stop1_10, 0]
y_avg = (y_01 + y_10) / 2

f11_10_avg = nanmean(F11_10, axis=1)[start1_10:stop1_10]
f11_01_avg = nanmean(F11_01, axis=1)[start1_01:stop1_01]
f12_10_avg = nanmean(F12_10, axis=1)[start1_10:stop1_10]
f12_01_avg = nanmean(F12_01, axis=1)[start1_01:stop1_01]

""" Left-right difference """
# Left side
x_left_01 = ma.masked_less(X_01, 0)
mask_left_01 = x_left_01.mask
x_left_10 = ma.masked_less(X_10, 0)
mask_left_10 = x_left_10.mask

y_01_left = y_01 #ma.array(Y_01, mask=mask_left)[start1_01:stop1_01, 0]
y_10_left = y_10 #ma.array(Y_10, mask=mask_left)[start1_10:stop1_10, 0]
y_avg_left = (y_01_left + y_10_left) / 2

f11_01_left = ma.array(F11_01, mask=x_left_01.mask)
f11_10_left = ma.array(F11_10, mask=x_left_10.mask)
f12_01_left = ma.array(F12_01, mask=x_left_01.mask)
f12_10_left = ma.array(F12_10, mask=x_left_10.mask)

f11_01_avg_left = nanmean(f11_01_left, axis=1)
f11_10_avg_left = nanmean(f11_10_left, axis=1)
f12_01_avg_left = nanmean(f12_01_left, axis=1)
f12_10_avg_left = nanmean(f12_10_left, axis=1)

# Right side
x_right_01 = ma.masked_greater_equal(X_01, 0)
mask_right_01 = x_right_01.mask
x_right_10 = ma.masked_greater_equal(X_10, 0)
mask_right_10 = x_right_10.mask

y_01_right = y_01 #ma.array(Y_01, mask=mask_right)[start1_01:stop1_01, 0]
y_10_right = y_10 #ma.array(Y_10, mask=mask_right)[start1_10:stop1_10, 0]
y_avg_right = (y_01_right + y_10_right) / 2

f11_01_right = ma.array(F11_01, mask=x_right_01.mask)
f11_10_right = ma.array(F11_10, mask=x_right_10.mask)
f12_01_right = ma.array(F12_01, mask=x_right_01.mask)
f12_10_right = ma.array(F12_10, mask=x_right_10.mask)

f11_01_avg_right = nanmean(f11_01_right, axis=1)
f11_10_avg_right = nanmean(f11_10_right, axis=1)
f12_01_avg_right = nanmean(f12_01_right, axis=1)
f12_10_avg_right = nanmean(f12_10_right, axis=1)


""" Differences and sums """
# All points (both sides)
diff11 = f11_10_avg - f11_01_avg
diff12 = f12_10_avg - f12_01_avg

avg11 = (nanmean(F11_10, axis=1)[start1_10:stop1_10]+nanmean(F11_01, axis=1)[start1_01:stop1_01])/2
avg12 = (nanmean(F12_10, axis=1)[start1_10:stop1_10]+nanmean(F12_01, axis=1)[start1_01:stop1_01])/2

# Left side
diff11_left = f11_10_avg_left-f11_01_avg_left
diff12_left = f12_10_avg_left-f12_01_avg_left

avg11_left = (f11_10_avg_left + f11_01_avg_left)/2
avg12_left = (f12_10_avg_left + f12_01_avg_left)/2

# Right side
diff11_right = f11_10_avg_right-f11_01_avg_right
diff12_right = f12_10_avg_right-f12_01_avg_right

avg11_right = (f11_10_avg_right + f11_01_avg_right)/2
avg12_right = (f12_10_avg_right + f12_01_avg_right)/2

""" Fit avg with dOmega_T """
if do_fits:
    def func_to_fit(y, omega0, dTmax):
        if mode == 1:
            return omega0 * (1 + fes.dOmega1_T(y, dTmax))
        elif mode == 2:
            return omega0 * (1 + fes.dOmega2_T(y, dTmax))
    vect_to_fit = vectorize(func_to_fit)

    # 1st Peak
    print('Fitting mode 11')
    to_fit11 = 2*pi*avg11[where(logical_not(isnan(avg11)))[0]]
    y11 = y_avg[where(logical_not(isnan(avg11)))[0]]
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

    # 1st Peak - left
    print('Fitting mode 11 - left')
    to_fit11_left = 2*pi*avg11_left[where(logical_not(isnan(avg11_left)))[0]]
    y11_left = y_avg_left[where(logical_not(isnan(avg11_left)))[0]]
    if len(to_fit11_left) > 0:
        guess = [2*pi*to_fit11_left.max(), 100]
        bounds_inf = (to_fit11_left.min(), 0)
        bounds_sup = (inf, inf)

        fit = curve_fit(vect_to_fit,
                        y11_left,
                        to_fit11_left,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O11_left, dT11_left = fit[0]

    # 1st Peak - right
    print('Fitting mode 11 - right')
    to_fit11_right = 2*pi*avg11_right[where(logical_not(isnan(avg11_right)))[0]]
    y11_right = y_avg_right[where(logical_not(isnan(avg11_right)))[0]]
    if len(to_fit11_right) > 0:
        guess = [2*pi*to_fit11_right.max(), 100]
        bounds_inf = (to_fit11_right.min(), 0)
        bounds_sup = (inf, inf)

        fit = curve_fit(vect_to_fit,
                        y11_right,
                        to_fit11_right,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O11_right, dT11_right = fit[0]

    # 2nd Peak
    print('Fitting mode 12')
    to_fit12 = 2*pi*avg12[where(logical_not(isnan(avg12)))[0]]
    y12 = y_avg[where(logical_not(isnan(avg12)))[0]]
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

    # 2nd Peak - left
    print('Fitting mode 12 - left')
    to_fit12_left = 2*pi*avg12_left[where(logical_not(isnan(avg12_left)))[0]]
    y12_left = y_avg_left[where(logical_not(isnan(avg12_left)))[0]]
    if len(to_fit12_left) > 0:
        guess = [2*pi*to_fit12_left.max(), 100]
        bounds_inf = (to_fit12_left.min(), 0)
        bounds_sup = (inf, inf)

        fit = curve_fit(vect_to_fit,
                        y12_left,
                        to_fit12_left,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O12_left, dT12_left = fit[0]

    # 2nd Peak - right
    print('Fitting mode 12 - right')
    to_fit12_right = 2*pi*avg12_right[where(logical_not(isnan(avg12_right)))[0]]
    y12_right = y_avg_right[where(logical_not(isnan(avg12_right)))[0]]
    if len(to_fit12_right) > 0:
        guess = [2*pi*to_fit12_right.max(), 100]
        bounds_inf = (to_fit12_right.min(), 0)
        bounds_sup = (inf, inf)

        fit = curve_fit(vect_to_fit,
                        y12_right,
                        to_fit12_right,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O12_right, dT12_right = fit[0]

""" Back-action """
L = 1

# With uniform beta
beta0 = 1
def beta(y, side=1):
    """ side should be 1 or -1 """
    return side * beta0
def kBA(y0, val, side=1):
    """ side should be 1 or -1 """
    if mode == 1:
        phi = fes.phi1
    elif mode == 2:
        phi = fes.phi2
    else:
            raise ValueError('Mode should be 1 or 2')
    alpha = (1+1j)/sqrt(2) * val/L
    int1 = quad(lambda y: beta(y, side)*phi(y)*sinh(alpha*y)/sinh(alpha*y0),
                 0,
                 y0)[0]
    int2 = quad(lambda y: beta(y, side)*phi(y)*cosh(alpha*(y-L))/cosh(alpha*(y0-L)),
                 0,
                 y0)[0]
    return -alpha * sinh(alpha*y0) * cosh(alpha*(y0-L)) / cosh(alpha*L) * phi(y0) / phi(L) * (int1 + int2)

# With Dirac beta
##def kBA(y0, yD=.4*L):
##    if mode == 1:
##        phi = fes.phi1
##    elif mode == 2:
##        phi = fes.phi2
##    else:
##            raise ValueError('Mode should be 1 or 2')
##    int1 = phi(yD)*sinh(alpha*yD)/sinh(alpha*y0)
##    int2 = phi(yD)*cosh(alpha*(yD-L))/cosh(alpha*(y0-L))
##    return -alpha * sinh(alpha*y0) * cosh(alpha*(y0-L)) / cosh(alpha*L) * phi(y0) / phi(L) * (int1 + int2)

if do_fits:
    def func_to_fit_ba(y, omega0, dTmax, A, val):
        if mode == 1:
            return omega0 * (1 + fes.dOmega1_T(y, dTmax)+A*real(kBA(y, val, side)))
        elif mode == 2:
            return omega0 * (1 + fes.dOmega2_T(y, dTmax)+A*real(kBA(y, val, side)))
        else:
            raise ValueError('Mode should be 1 or 2')
    vect_to_fit_ba = vectorize(func_to_fit_ba)

    # 1st Peak
    print('Fitting mode 11 with backaction')
    to_fit11 = 2*pi*avg11[where(logical_not(isnan(avg11)))[0]]
    y11 = y_avg[where(logical_not(isnan(avg11)))[0]]
    if len(to_fit11) > 0:
        side = 1
        guess = [2*pi*to_fit11.max(), 100, 1, 10]
        bounds_inf = (to_fit11.min(), 0, 0, 0)
        bounds_sup = (inf, inf, inf, inf)

        fit = curve_fit(vect_to_fit_ba,
                        y11,
                        to_fit11,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O11_ba, dT11_ba, A11, val11 = fit[0]
        lab11 = '$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$\n$A_{11} = %.2f\\times 10^{-6}$\n$(\\alpha L)_{11} = %.2f$'%(O11_ba/2/pi/1e3, dT11_ba, A11/1e-6, val11)

    # 1st Peak - left
    print('Fitting mode 11 with backaction - left')
    to_fit11_left = 2*pi*avg11_left[where(logical_not(isnan(avg11_left)))[0]]
    y11_left = y_avg_left[where(logical_not(isnan(avg11_left)))[0]]
    if len(to_fit11_left) > 0:
        side = -1
        guess = [2*pi*to_fit11_left.max(), 100, 1, 10]
        bounds_inf = (to_fit11_left.min(), 0, 0, 0)
        bounds_sup = (inf, inf, inf, inf)

        fit = curve_fit(vect_to_fit_ba,
                        y11_left,
                        to_fit11_left,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O11_ba_left, dT11_ba_left, A11_left, val11_left = fit[0]
        lab11_left = '$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$\n$A_{11} = %.2f\\times 10^{-6}$\n$(\\alpha L)_{11} = %.2f$'%(O11_ba_left/2/pi/1e3, dT11_ba_left, A11_left/1e-6, val11_left)

    # 1st Peak - right
    print('Fitting mode 11 with backaction - right')
    to_fit11_right = 2*pi*avg11_right[where(logical_not(isnan(avg11_right)))[0]]
    y11_right = y_avg_right[where(logical_not(isnan(avg11_right)))[0]]
    if len(to_fit11_right) > 0:
        side = 1
        guess = [2*pi*to_fit11_right.max(), 100, 1, 10]
        bounds_inf = (to_fit11_right.min(), 0, 0, 0)
        bounds_sup = (inf, inf, inf, inf)

        fit = curve_fit(vect_to_fit_ba,
                        y11_right,
                        to_fit11_right,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O11_ba_right, dT11_ba_right, A11_right, val11_right = fit[0]
        lab11_right = '$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$\n$A_{11} = %.2f\\times 10^{-6}$\n$(\\alpha L)_{11} = %.2f$'%(O11_ba_right/2/pi/1e3, dT11_ba_right, A11_right/1e-6, val11_right)

    # 2nd Peak
    print('Fitting mode 12 with backaction')
    to_fit12 = 2*pi*avg12[where(logical_not(isnan(avg12)))[0]]
    y12 = y_avg[where(logical_not(isnan(avg12)))[0]]
    if len(to_fit12) > 0:
        side = 1
        guess = [2*pi*to_fit12.max(), 100, 1, 10]
        bounds_inf = (to_fit12.min(), 0, 0, 0)
        bounds_sup = (inf, inf, inf, inf)

        fit = curve_fit(vect_to_fit_ba,
                        y12,
                        to_fit12,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O12_ba, dT12_ba, A12, val12 = fit[0]
        lab12 = '$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$\n$A_{12} = %.2f\\times 10^{-6}$\n$(\\alpha L)_{12} = %.2f$'%(O12_ba/2/pi/1e3, dT12_ba, A12/1e-6, val12)

    # 2nd Peak - left
    print('Fitting mode 12 with backaction - left')
    to_fit12_left = 2*pi*avg12_left[where(logical_not(isnan(avg12_left)))[0]]
    y12_left = y_avg_left[where(logical_not(isnan(avg12_left)))[0]]
    if len(to_fit12_left) > 0:
        side = -1
        guess = [2*pi*to_fit12_left.max(), 100, 1, 10]
        bounds_inf = (to_fit12_left.min(), 0, 0, 0)
        bounds_sup = (inf, inf, inf, inf)

        fit = curve_fit(vect_to_fit_ba,
                        y12_left,
                        to_fit12_left,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O12_ba_left, dT12_ba_left, A12_left, val12_left = fit[0]
        lab12_left = '$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$\n$A_{12} = %.2f\\times 10^{-6}$\n$(\\alpha L)_{12} = %.2f$'%(O12_ba_left/2/pi/1e3, dT12_ba_left, A12_left/1e-6, val12_left)

    # 2nd Peak - right
    print('Fitting mode 12 with backaction - right')
    to_fit12_right = 2*pi*avg12_right[where(logical_not(isnan(avg12_right)))[0]]
    y12_right = y_avg_right[where(logical_not(isnan(avg12_right)))[0]]
    if len(to_fit12_right) > 0:
        side = 1
        guess = [2*pi*to_fit12_right.max(), 100, 1, 10]
        bounds_inf = (to_fit12_right.min(), 0, 0, 0)
        bounds_sup = (inf, inf, inf, inf)

        fit = curve_fit(vect_to_fit_ba,
                        y12_right,
                        to_fit12_right,
                        p0=guess,
                        bounds=(bounds_inf, bounds_sup))
        O12_ba_right, dT12_ba_right, A12_right, val12_right = fit[0]
        lab12_right = '$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$\n$A_{12} = %.2f\\times 10^{-6}$\n$(\\alpha L)_{12} = %.2f$'%(O12_ba_right/2/pi/1e3, dT12_ba_right, A12_right/1e-6, val12_right)

""" Do the plots """
if savefigs:
    figsdir = dirname + '\\Figures\\Sum_diff_%s_%s\\Back-action'%(batch_01, batch_10)
    if not os.path.isdir(figsdir):
        os.makedirs(figsdir)

yy = linspace(0, 1, 101)

""" Figure : Mode 11 """
fig11, ax11 = subplots()
ax11.set_xlabel('y/L')
ax11.set_ylabel('Frequency (kHz)')
ax11.set_title('1st mechanical mode, 1st peak')

# Plot raw data
ax11.plot(Y_01.flatten(), f11_01/1e3, 'o', label=batch_01)
ax11.plot(Y_01[start1_01:stop1_01, :].flatten(), f11_01_left.flatten()/1e3, 'o', label=batch_01+' left')
ax11.plot(Y_01[start1_01:stop1_01, :].flatten(), f11_01_right.flatten()/1e3, 'o', label=batch_01+' right')

ax11.plot(Y_10.flatten(), f11_10/1e3, 'o', label=batch_10)
ax11.plot(Y_10[start1_10:stop1_10, :].flatten(), f11_10_left.flatten()/1e3, 'o', label=batch_10+' left')
ax11.plot(Y_10[start1_10:stop1_10, :].flatten(), f11_10_right.flatten()/1e3, 'o', label=batch_10+' right')

# Plot mean: points kept for sum and difference
ax11.plot(Y_01[:, 0], nanmean(F11_01, axis=1)/1e3, 'sk', mfc='none')
ax11.plot(y_01_left, f11_01_avg_left/1e3, 'dk', mfc='none')
ax11.plot(y_01_right, f11_01_avg_right/1e3, 'pk', mfc='none')

ax11.plot(Y_10[:, 0], nanmean(F11_10, axis=1)/1e3, 'sk', mfc='none')
ax11.plot(y_10_left, f11_10_avg_left/1e3, 'dk', mfc='none')
ax11.plot(y_10_right, f11_10_avg_right/1e3, 'pk', mfc='none')

# Plot averages
ax11.plot(y_avg, avg11/1e3, 'o', label='Avg')
ax11.plot(y_avg_left, avg11_left/1e3, 'o', label='Avg, left')
ax11.plot(y_avg_right, avg11_right/1e3, 'o', label='Avg, right')

# Fits
if do_fits:
    ax11.plot(yy, vect_to_fit(yy, O11, dT11)/2/pi/1e3, label='Fit:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O11/2/pi/1e3, dT11))
    ax11.plot(yy, vect_to_fit_ba(yy, O11_ba, dT11_ba, A11, val11)/2/pi/1e3,
              label='Fit with BA:\n'+lab11)

# Legend
ax11.legend(bbox_to_anchor=(1., 1.))
fig11.tight_layout()

""" Figure : Mode 12 """
fig12, ax12 = subplots()
ax12.set_xlabel('y/L')
ax12.set_ylabel('Frequency (kHz)')
ax12.set_title('1st mechanical mode, 2nd peak')

# Plot raw data
ax12.plot(Y_01.flatten(), f12_01/1e3, 'o', label=batch_01)
ax12.plot(Y_01[start1_01:stop1_01, :].flatten(), f12_01_left.flatten()/1e3, 'o', label=batch_01+' left')
ax12.plot(Y_01[start1_01:stop1_01, :].flatten(), f12_01_right.flatten()/1e3, 'o', label=batch_01+' right')

ax12.plot(Y_10.flatten(), f12_10/1e3, 'o', label=batch_10)
ax12.plot(Y_10[start1_10:stop1_10, :].flatten(), f12_10_left.flatten()/1e3, 'o', label=batch_10+' left')
ax12.plot(Y_10[start1_10:stop1_10, :].flatten(), f12_10_right.flatten()/1e3, 'o', label=batch_10+' right')

# Plot mean: points kept for sum and difference
ax12.plot(Y_01[:, 0], nanmean(F12_01, axis=1)/1e3, 'sk', mfc='none')
ax12.plot(y_01_left, f12_01_avg_left/1e3, 'dk', mfc='none')
ax12.plot(y_01_right, f12_01_avg_right/1e3, 'pk', mfc='none')

ax12.plot(Y_10[:, 0], nanmean(F12_10, axis=1)/1e3, 'sk', mfc='none')
ax12.plot(y_10_left, f12_10_avg_left/1e3, 'dk', mfc='none')
ax12.plot(y_10_right, f12_10_avg_right/1e3, 'pk', mfc='none')

# Plot averages
ax12.plot(y_avg, avg12/1e3, 'o', label='Avg')
ax12.plot(y_avg_left, avg12_left/1e3, 'o', label='Avg, left')
ax12.plot(y_avg_right, avg12_right/1e3, 'o', label='Avg, right')

# Legend
ax12.legend(bbox_to_anchor=(1., 1.))
fig12.tight_layout()

# Fits
if do_fits:
    ax12.plot(yy, vect_to_fit(yy, O12, dT12)/2/pi/1e3, label='Fit:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O12/2/pi/1e3, dT12))
    ax12.plot(yy, vect_to_fit_ba(yy, O12_ba, dT12_ba, A12, val12)/2/pi/1e3,
              label='Fit with BA:\n'+lab12)

""" Figure : Difference Mode 1 """
fig_diff1, ax_diff1 = subplots()
ax_diff1.set_xlabel('y/L')
ax_diff1.set_ylabel('Frequency difference (kHz) Mode 1')
ax_diff1.set_title('Difference "10" - "01"')

ax_diff1.plot(y_avg, diff11/1e3, 'o', label='Peak 1')
ax_diff1.plot(y_avg_left, diff11_left/1e3, 'o', label='Peak 1, left')
ax_diff1.plot(y_avg_right, diff11_right/1e3, 'o', label='Peak 1, right')

ax_diff1.plot(y_avg, diff12/1e3, 'o', label='Peak 2')
ax_diff1.plot(y_avg_left, diff12_left/1e3, 'o', label='Peak 2, left')
ax_diff1.plot(y_avg_right, diff12_right/1e3, 'o', label='Peak 2, right')

# Legend
ax_diff1.legend(bbox_to_anchor=(1., 1.))
fig_diff1.tight_layout()

""" Figure : Average Mode 11 """
fig_sum11, ax_sum11 = subplots()
ax_sum11.set_xlabel('y/L')
ax_sum11.set_ylabel('Frequency difference (kHz)')
ax_sum11.set_title('Average ("10" + "01") / 2 Mode 1')

ax_sum11.plot(y_avg, avg11/1e3, 'o', c='C0', label='Peak 1')
if do_fits:
    ax_sum11.plot(yy, vect_to_fit(yy, O11, dT11)/2/pi/1e3, '--', c='C0', label='Fit peak 1:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O11/2/pi/1e3, dT11))
    ax_sum11.plot(yy, vect_to_fit_ba(yy, O11_ba, dT11_ba, A11, val11)/2/pi/1e3, c='C0', 
                 label='Fit peak 1 with BA:\n'+lab11+'\n')

ax_sum11.plot(y_avg_left, avg11_left/1e3, 'o', c='C1', label='Peak 1, left')
if do_fits:
    ax_sum11.plot(yy, vect_to_fit(yy, O11_left, dT11_left)/2/pi/1e3, '--', c='C1', label='Fit peak 1 left:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O11_left/2/pi/1e3, dT11_left))
    ax_sum11.plot(yy, vect_to_fit_ba(yy, O11_ba_left, dT11_ba_left, A11_left, val11_left)/2/pi/1e3, c='C1', 
                 label='Fit peak 1 left with BA:\n'+lab11_left+'\n')

ax_sum11.plot(y_avg_right, avg11_right/1e3, 'o', c='C2', label='Peak 1, right')
if do_fits:
    ax_sum11.plot(yy, vect_to_fit(yy, O11_right, dT11_right)/2/pi/1e3, '--', c='C2', label='Fit peak 1 right:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O11_right/2/pi/1e3, dT11_right))
    ax_sum11.plot(yy, vect_to_fit_ba(yy, O11_ba_right, dT11_ba_right, A11_right, val11_right)/2/pi/1e3, c='C2', 
                 label='Fit peak 1 right with BA:\n'+lab11_right+'\n')

# Legend
ax_sum11.legend(bbox_to_anchor=(1., 1.))
fig_sum11.tight_layout()

""" Figure : Average Mode 12 """
fig_sum12, ax_sum12 = subplots()
ax_sum12.set_xlabel('y/L')
ax_sum12.set_ylabel('Frequency difference (kHz)')
ax_sum12.set_title('Average ("10" + "01") / 2 Mode 1')

ax_sum12.plot(y_avg, avg12/1e3, 'o', c='C0', label='Peak 2')
if do_fits:
    ax_sum12.plot(yy, vect_to_fit(yy, O12, dT12)/2/pi/1e3, '--', c='C0', label='Fit peak 1:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O12/2/pi/1e3, dT12))
    ax_sum12.plot(yy, vect_to_fit_ba(yy, O12_ba, dT12_ba, A12, val12)/2/pi/1e3, c='C0', 
                 label='Fit peak 1 with BA:\n'+lab12+'\n')

ax_sum12.plot(y_avg_left, avg12_left/1e3, 'o', c='C1', label='Peak 2, left')
if do_fits:
    ax_sum12.plot(yy, vect_to_fit(yy, O12_left, dT12_left)/2/pi/1e3, '--', c='C1', label='Fit peak 1 left:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O12_left/2/pi/1e3, dT12_left))
    ax_sum12.plot(yy, vect_to_fit_ba(yy, O12_ba_left, dT12_ba_left, A12_left, val12_left)/2/pi/1e3, c='C1',
                 label='Fit peak 1 left with BA:\n'+lab12_left+'\n')

ax_sum12.plot(y_avg_right, avg12_right/1e3, 'o', c='C2', label='Peak 2, right')
if do_fits:
    ax_sum12.plot(yy, vect_to_fit(yy, O12_right, dT12_right)/2/pi/1e3, '--', c='C2', label='Fit peak 1 right:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O12_right/2/pi/1e3, dT12_right))
    ax_sum12.plot(yy, vect_to_fit_ba(yy, O12_ba_right, dT12_ba_right, A12_right, val12_right)/2/pi/1e3, c='C2',
                 label='Fit peak 1 right with BA:\n'+lab12_right)

# Legend
ax_sum12.legend(bbox_to_anchor=(1., 1.))
fig_sum12.tight_layout()

""" Figure : splitting """
fig_splitting, ax_splitting = subplots()
ax_splitting.set_xlabel('y/L')
ax_splitting.set_ylabel('Mode Splitting (kHz)')
ax_splitting.set_title('Mode Splitting')

plot(y_avg_left, (f12_01_avg_left-f11_01_avg_left)/1e3, 'o', label='Splitting 01, left')
plot(y_avg_right, (f12_01_avg_right-f11_01_avg_right)/1e3, 'o', label='Splitting 01, right')

plot(y_avg_left, (f12_10_avg_left-f11_10_avg_left)/1e3, 'o', label='Splitting 10, left')
plot(y_avg_right, (f12_10_avg_right-f11_10_avg_right)/1e3, 'o', label='Splitting 10, right')

plot(y_avg, (avg12_left-avg11_left)/1e3, 'o', label='Splitting avg, left')
plot(y_avg, (avg12_right-avg11_right)/1e3, 'o', label='Splitting avg, right')

legend()
fig_splitting.tight_layout()

show()

""" Save figures """
if savefigs:
    fig11.savefig(figsdir+'\\Mode_11')
    fig12.savefig(figsdir+'\\Mode_12')

    fig_diff1.savefig(figsdir+'\\Difference_1')
    fig_sum1.savefig(figsdir+'\\Sum_1')
    print('Figures saved')
else:
    print('savefigs set to', savefigs)

show()
