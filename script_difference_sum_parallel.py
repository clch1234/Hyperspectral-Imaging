# -*- coding: cp1252 -*-
"""
Script for comparing the frequency evolutions for multiple acquisitions on a single NW
Here only for the mode measured in parallel

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

dat = '20220622'
batch_01 = 'Fil 6_1'
batch_10 = 'Fil 6_2'
rectangle = '1'
mode = 2 # nb of mechanical mode

savefigs = True
do_fits = True

dirname = "D:\\Documents\\Boulot\\Grenoble\\Data\\%s"%dat

""" Create figures """
if savefigs:
    figsdir = dirname + '\\Figures\\Sum_diff_%s_%s'%(batch_01, batch_10)
    if not os.path.isdir(figsdir):
        os.makedirs(figsdir)

idx_color = 0

fig21, ax21 = subplots()
ax21.set_xlabel('y/L')
ax21.set_ylabel('Frequency (kHz)')
ax21.set_title('2nd mechanical mode, 1st peak')

fig22, ax22 = subplots()
ax22.set_xlabel('y/L')
ax22.set_ylabel('Frequency (kHz)')
ax22.set_title('2nd mechanical mode, 2nd peak')

fig_diff2, ax_diff2 = subplots()
ax_diff2.set_xlabel('y/L')
ax_diff2.set_ylabel('Frequency difference (kHz)')
ax_diff2.set_title('Difference "10" - "01" ; Mode 2')

fig_sum2, ax_sum2 = subplots()
ax_sum2.set_xlabel('y/L')
ax_sum2.set_ylabel('Frequency difference (kHz)')
ax_sum2.set_title('Average ("10" + "01") / 2 ; Mode 2')

""" Load and plot """
dir_coordinates_01 = dirname + u'\\%s\\Data_files'%batch_01
dir_bitmap_01 = dirname + u'\\%s\\Rectangle %s'%(batch_01, rectangle)
dir_coordinates_10 = dirname + u'\\%s\\Data_files'%batch_10
dir_bitmap_10 = dirname + u'\\%s\\Rectangle %s'%(batch_10, rectangle)

""" Load mode 2 """
dir_mode2_01 = dirname + u'\\%s\\Data_files\\Parallel_acquisition'%batch_01
dir_mode2_10 = dirname + u'\\%s\\Data_files\\Parallel_acquisition'%batch_10

f21s_01 = load(dir_mode2_01+'\\freqs1.npy', allow_pickle=True)
ys21_01 = load(dir_mode2_01+'\\ys1.npy', allow_pickle=True)

f22s_01 = load(dir_mode2_01+'\\freqs2.npy', allow_pickle=True)
ys22_01 = load(dir_mode2_01+'\\ys2.npy', allow_pickle=True)

f21s_10 = load(dir_mode2_10+'\\freqs1.npy', allow_pickle=True)
ys21_10 = load(dir_mode2_10+'\\ys1.npy', allow_pickle=True)

f22s_10 = load(dir_mode2_10+'\\freqs2.npy', allow_pickle=True)
ys22_10 = load(dir_mode2_10+'\\ys2.npy', allow_pickle=True)


# Normalize y by length
# Load X, Y and bitmap
bitmap_01 = loadtxt(dir_bitmap_01+u'\\Bitmap.txt')
X_ini_01 = loadtxt(dir_coordinates_01+u'\\%s_rectangle_1_raw_X.txt'%batch_01)
Y_ini_01 = loadtxt(dir_coordinates_01+u'\\%s_rectangle_1_raw_Y.txt'%batch_01)

bitmap_10 = loadtxt(dir_bitmap_10+u'\\Bitmap.txt')
X_ini_10 = loadtxt(dir_coordinates_10+u'\\%s_rectangle_1_raw_X.txt'%batch_10)
Y_ini_10 = loadtxt(dir_coordinates_10+u'\\%s_rectangle_1_raw_Y.txt'%batch_10)

# Detection of the edges of the NW and left-right separation
centers_01, edges_01, widths_01, first_line_found_01, last_line_found_01 = detect_NW_edges(bitmap_01, X_ini_01)
centers_10, edges_10, widths_10, first_line_found_10, last_line_found_10 = detect_NW_edges(bitmap_10, X_ini_10)

# Deduce basis position and NW length
length_01 = abs(Y_ini_01[first_line_found_01, 0] - Y_ini_01[last_line_found_01, 0])
basis_01 = Y_ini_01[first_line_found_01, 0]

ys21_01 = (ys21_01 - basis_01) / length_01
ys22_01 = (ys22_01 - basis_01) / length_01

length_10 = abs(Y_ini_10[first_line_found_10, 0] - Y_ini_10[last_line_found_10, 0])
basis_10 = Y_ini_10[first_line_found_10, 0]

ys21_10 = (ys21_10 - basis_10) / length_10
ys22_10 = (ys22_10 - basis_10) / length_10

# Plot

ax21.plot(ys21_01, f21s_01/1e3, 'o', label=batch_01)
ax22.plot(ys22_01, f22s_01/1e3, 'o', label=batch_01)

ax21.plot(ys21_10, f21s_10/1e3, 'o', label=batch_10)
ax22.plot(ys22_10, f22s_10/1e3, 'o', label=batch_10)

""" Take averages for each position """
# Mode 2.1
y21_01_new = []
f21_01_new = []
for ii, y in enumerate(ys21_01):
    if y in y21_01_new:
        pass
    else:
        y21_01_new.append(y)
        f = nanmean(f21s_01[where(ys21_01==y)[0]])
        f21_01_new.append(f)
y21_01_new = array(y21_01_new)
f21_01_new = array(f21_01_new)

y21_10_new = []
f21_10_new = []
for ii, y in enumerate(ys21_10):
    if y in y21_10_new:
        pass
    else:
        y21_10_new.append(y)
        f = nanmean(f21s_10[where(ys21_10==y)[0]])
        f21_10_new.append(f)
y21_10_new = array(y21_10_new)
f21_10_new = array(f21_10_new)

# Mode 2.2
y22_01_new = []
f22_01_new = []
for ii, y in enumerate(ys22_01):
    if y in y22_01_new:
        pass
    else:
        y22_01_new.append(y)
        f = nanmean(f22s_01[where(ys22_01==y)[0]])
        f22_01_new.append(f)
y22_01_new = array(y22_01_new)
f22_01_new = array(f22_01_new)

y22_10_new = []
f22_10_new = []
for ii, y in enumerate(ys22_10):
    if y in y22_10_new:
        pass
    else:
        y22_10_new.append(y)
        f = nanmean(f22s_10[where(ys22_10==y)[0]])
        f22_10_new.append(f)
y22_10_new = array(y22_10_new)
f22_10_new = array(f22_10_new)

""" Harmonize sizes """
if batch_01 == 'Fil 6_3' and batch_10 == 'Fil 6_4':
    start21_01 = 0
    stop21_01 = len(y21_01_new)
    start21_10 = 0
    stop21_10 = len(y21_10_new)

    start22_01 = 1
    stop22_01 = len(y22_01_new)
    start22_10 = 1
    stop22_10 = len(y22_10_new)

    # alternative:
    pts_to_ignore_21_01 = None
    pts_to_ignore_21_10 = None
    pts_to_ignore_22_01 = None
    pts_to_ignore_22_10 = None

elif batch_01 == 'Fil 6_1' and batch_10 == 'Fil 6_2':
    start21_01 = -1
    stop21_01 = 0
    start21_10 = -1
    stop21_10 = 0

    start22_01 = 0
    stop22_01 = len(y22_01_new)
    start22_10 = 1
    stop22_10 = len(y22_10_new)

    # alternative:
    pts_to_ignore_21_01 = None
    pts_to_ignore_21_10 = None
    pts_to_ignore_22_01 = [9, 10]
    pts_to_ignore_22_10 = [0, 3, 4]

elif batch_01 == 'Fil 5_11' and batch_10 == 'Fil 5_10':
    start21_01 = 1
    stop21_01 = len(y21_01_new)
    start21_10 = 1
    stop21_10 = len(y21_10_new)

    start22_01 = 2
    stop22_01 = len(y22_01_new)
    start22_10 = 0
    stop22_10 = -1

    # alternative:
    pts_to_ignore_21_01 = None
    pts_to_ignore_21_10 = None
    pts_to_ignore_22_01 = None
    pts_to_ignore_22_10 = None

else:
    start21_01 = 0
    stop21_01 = len(y21_01_new)
    start21_10 = 0
    stop21_10 = len(y21_10_new)

    start22_01 = 0
    stop22_01 = len(y22_01_new)
    start22_10 = 0
    stop22_10 = len(y22_10_new)

    # alternative:
    pts_to_ignore_21_01 = None
    pts_to_ignore_21_10 = None
    pts_to_ignore_22_01 = None
    pts_to_ignore_22_10 = None

# Mode 2.1
if pts_to_ignore_21_01 is None:
    y21_01_new = y21_01_new[start21_01:stop21_01]
    f21_01_new = f21_01_new[start21_01:stop21_01]
else:
    y21_01_new = y21_01_new[[ii for ii in range(len(y21_01_new)) if not ii in pts_to_ignore_21_01]]
    f21_01_new = f21_01_new[[ii for ii in range(len(f21_01_new)) if not ii in pts_to_ignore_21_01]]
if pts_to_ignore_21_10 is None:
    y21_10_new = y21_10_new[start21_10:stop21_10]
    f21_10_new = f21_10_new[start21_10:stop21_10]
else:
    y21_10_new = y21_10_new[[ii for ii in range(len(y21_10_new)) if not ii in pts_to_ignore_21_10]]
    f21_10_new = f21_10_new[[ii for ii in range(len(f21_10_new)) if not ii in pts_to_ignore_21_10]]

# Mode 2.2
if pts_to_ignore_22_01 is None:
    y22_01_new = y22_01_new[start22_01:stop22_01]
    f22_01_new = f22_01_new[start22_01:stop22_01]
else:
    y22_01_new = y22_01_new[[ii for ii in range(len(y22_01_new)) if not ii in pts_to_ignore_22_01]]
    f22_01_new = f22_01_new[[ii for ii in range(len(f22_01_new)) if not ii in pts_to_ignore_22_01]]
if pts_to_ignore_22_10 is None:
    y22_10_new = y22_10_new[start22_10:stop22_10]
    f22_10_new = f22_10_new[start22_10:stop22_10]
else:
    y22_10_new = y22_10_new[[ii for ii in range(len(y22_10_new)) if not ii in pts_to_ignore_22_10]]
    f22_10_new = f22_10_new[[ii for ii in range(len(f22_10_new)) if not ii in pts_to_ignore_22_10]]

""" Differences and sums """
# Invert direction 01 everywhere !
y_avg21 = (y21_01_new[::-1] + y21_10_new)/2
avg21 = (f21_10_new + f21_01_new[::-1])/2
diff21 = f21_10_new - f21_01_new[::-1]

y_avg22 = (y22_01_new[::-1] + y22_10_new)/2
avg22 = (f22_10_new + f22_01_new[::-1])/2
diff22 = f22_10_new - f22_01_new[::-1]

# Plot
ax21.plot(y_avg21, avg21/1e3, 'o')
ax21.plot(y21_01_new, f21_01_new/1e3, 's', mec='k', mfc='none')
ax21.plot(y21_10_new, f21_10_new/1e3, 's', mec='k', mfc='none')
ax_sum2.plot(y_avg21, avg21/1e3, 'o')
ax_diff2.plot(y_avg21, diff21/1e3, 'o')

ax22.plot(y_avg22, avg22/1e3, 'o')
ax22.plot(y22_01_new, f22_01_new/1e3, 's', mec='k', mfc='none')
ax22.plot(y22_10_new, f22_10_new/1e3, 's', mec='k', mfc='none')

ax_sum2.plot(y_avg22, avg22/1e3, 'o')
ax_diff2.plot(y_avg22, diff22/1e3, 'o')
""" Fit avg with dOmega_T """
if do_fits:
    def func_to_fit(y, omega0, dTmax):
        if mode == 1:
            return omega0 * (1 + fes.dOmega1_T(y, dTmax))
        elif mode == 2:
            return omega0 * (1 + fes.dOmega2_T(y, dTmax))
    vect_to_fit = vectorize(func_to_fit)

    # 1st Peak
    to_fit21 = 2*pi*avg21[where(logical_not(isnan(avg21)))[0]]
    if len(to_fit21) > 0:
        guess = [2*pi*to_fit21.max(), 100]
        bounds_inf = (to_fit21.min(), 0)
        bounds_sup = (inf, inf)

        fit21 = curve_fit(vect_to_fit,
                          y_avg21,
                          to_fit21,
                          p0=guess,
                          bounds=(bounds_inf, bounds_sup))
        O21, dT21 = fit21[0]


        yy = linspace(0, 1, 101)
        ax_sum2.plot(yy, vect_to_fit(yy, O21, dT21)/2/pi/1e3, label='Fit peak 1:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O21/2/pi/1e3, dT21))
        ax21.plot(yy, vect_to_fit(yy, O21, dT21)/2/pi/1e3, label='Fit:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O21/2/pi/1e3, dT21))

    # 2nd Peak
    to_fit22 = 2*pi*avg22[where(logical_not(isnan(avg22)))[0]]
    if len(to_fit22) > 0:
        guess = [2*pi*to_fit22.max(), 100]
        bounds_inf = (to_fit22.min(), 0)
        bounds_sup = (inf, inf)

        fit22 = curve_fit(vect_to_fit,
                          y_avg22,
                          to_fit22,
                          p0=guess,
                          bounds=(bounds_inf, bounds_sup))
        O22, dT22 = fit22[0]


        yy = linspace(0, 1, 101)
        ax_sum2.plot(yy, vect_to_fit(yy, O22, dT22)/2/pi/1e3, label='Fit peak 2:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O22/2/pi/1e3, dT22))
        ax22.plot(yy, vect_to_fit(yy, O22, dT22)/2/pi/1e3, label='Fit:\n$\Omega_0/2\pi = %.1f kHz$\n$\\Delta T_{max} = %i$'%(O22/2/pi/1e3, dT22))


ax21.legend()
ax22.legend()
fig21.tight_layout()
fig22.tight_layout()

ax_diff2.legend()
fig_diff2.tight_layout()
ax_sum2.legend()
fig_sum2.tight_layout()

show()

if savefigs:
    fig21.savefig(figsdir+'\\Mode_21')
    fig22.savefig(figsdir+'\\Mode_22')

    fig_diff2.savefig(figsdir+'\\Difference_2')
    fig_sum2.savefig(figsdir+'\\Sum_2')
    print('Figures saved')
else:
    print('savefigs set to', savefigs)

show()
