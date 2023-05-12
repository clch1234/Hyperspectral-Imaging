# -*- coding: cp1252 -*-
"""
Script for the acquisition made in june 2022 where we measured the first two mechanical modes
of an InAs nanowire in parallel, with two spectrum analysers

When using on a new dataset :
    - correspondence and sens need to be set in the parameters file
    - Run script once. It is probably going to end with an error.
    - Use following command :
        plot(times_RSA, signal_RSA)
    - Adjust theshold_NW_edge in parameters file
    - Use following command :
        plot(times_RSA, spectrogram.max(axis=1)-signal_RSA)
    - Adjust threshold_peak in parameters file

    - Restart kernell to reload parameters
    - Run again. If frequencies are mixed between modes 1 & 2 :
      set freq_min1, freq_middle and freq_max2 in the parameters file

    - Restart kernell to reload parameters
    - Run one last time, it should work properly now.

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
import os
#from scipy.stats import linregress
#from params_reconstitution_frequence_et_image_MEB import *
from matplotlib import colors
from scipy.optimize import curve_fit
#from ..dependances_scripts import displacements
import importlib.util
from functions import *
import parameters_parallel_acquisition as ppa

nw = '843'
fil = '6'
spec = '11'
dat = '20220622'
key = 'NW_'+nw+'_fil_'+fil+'_spec_'+spec.zfill(4)

savefigs = False
save_data = False

dirname = "D:\\Documents\\Boulot\\Grenoble\\Data\\%s"%dat
filename_oscillo = 'NW'+nw+'_fil'+fil+'_spec'+spec+'_oscillo.csv'
filename_RSA = 'NW'+nw+'_fil_'+fil+'_spectrogram_'+spec.zfill(4)+'.npy'
filename_times_RSA = 'NW'+nw+'_fil_'+fil+'_times_spectrogram_'+spec.zfill(4)+'.npy'
filename_freqs_RSA = 'NW'+nw+'_fil_'+fil+'_frequencies_spectrogram_'+spec.zfill(4)+'.txt'

""" Get corrrespondence with ViBR files """
corres = ppa.correspondences[key]
savedir = "D:\\Documents\\Boulot\\Grenoble\\Data\\%s\\%s"%(dat, corres)

""" Get dwell time from acquisition parameters """
filepath_params_acquisition = savedir+"\\Rectangle 1\\Parametres spectres.txt"
params_acquisition = {}
with open(filepath_params_acquisition, 'r') as ff:
    for line in ff.readlines():
        if line[-1] == '\n':
            line = line[:-1]
        k, v = line.split('\t')
        try:
            params_acquisition[k] = float(v)
        except ValueError:
            params_acquisition[k] = v
dwell_time = params_acquisition['Dwell Time (s)']

""" Load numerical parameters """
threshold_NW_edges = ppa.dd_threshold_NW_edges[key] if key in ppa.dd_threshold_NW_edges else -60
threshold_peak = ppa.dd_threshold_peak[key] if key in ppa.dd_threshold_peak else 12.5 # Nb of dB above noise level for peak detection
freq_min1 = ppa.dd_freq_min1[key] if key in ppa.dd_freq_min1 else None # Minimum frequency for mode 2.1
freq_middle = ppa.dd_freq_middle[key] if key in ppa.dd_freq_middle else None # Frequency between modes 2.1 and 2.2
freq_max2 = ppa.dd_freq_max2[key] if key in ppa.dd_freq_max2 else None # Maximum frequency for mode 2.2
sens = ppa.dd_sens[key]
if not sens in ('01', '10'):
    raise ValueError('Wrong value for "sens" : %s'%sens)

""" Load data from RSA """
if os.path.isfile(savedir+'\\'+filename_RSA.split('.')[0]+'_without_duplicates'+'.npy'):
    spectrogram = load(savedir+'\\'+filename_RSA.split('.')[0]+'_without_duplicates'+'.npy')
    times_RSA = load(savedir+'\\'+filename_times_RSA.split('.')[0]+'_without_duplicates'+'.npy')
    cleared = True
    print('Data RSA without duplicates loaded')
else:
    spectrogram = load(dirname+'\\'+filename_RSA, allow_pickle=True)
    times_RSA = load(dirname+'\\'+filename_times_RSA, allow_pickle=True)
    cleared = False
    print('Data RSA loaded')

if not cleared:
    print('Removing sprectrums with wrong number of points. Number of spectrums before : %i'%len(spectrogram))
    N = len(spectrogram[-1])
    S = []
    for spectrum in spectrogram:
        if len(spectrum) == N:
            S.append(array(spectrum))
    spectrogram = array(S)
    print('Done. Number of spectrums after : %i'%len(spectrogram))

    # Remove all duplicates :
    idx_to_keep = [0]
    for ii in range(1, len(spectrogram)):
        if False in (spectrogram[ii] == spectrogram[ii-1]):
            idx_to_keep.append(ii)
    print('Removing duplicates from spectrogram. Number of spectrums before : %i'%len(spectrogram))
    spectrogram = spectrogram[idx_to_keep]
    times_RSA = times_RSA[idx_to_keep]
    print('Duplicates removed from spectrogram. Number of spectrums after : %i'%len(spectrogram))

    save(savedir+'\\'+filename_RSA.split('.')[0] + '_without_duplicates'+'.npy', spectrogram)
    save(savedir+'\\'+filename_times_RSA.split('.')[0] + '_without_duplicates'+'.npy', times_RSA)
    print('Spectrogram and times saved without duplicates')

try:
    with open(dirname + '\\' + filename_freqs_RSA, 'r') as ff:
        lines = ff.readlines()

    f_centre = float(lines[0].split(' = ')[1])
    span = float(lines[1].split(' = ')[1])
    ff = linspace(f_centre-span/2, f_centre+span/2, len(spectrogram[0]))
except FileNotFoundError:
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print('/!\ Frequencies file not found, using indexes as frequencies /!\ ')
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    ff = array(range(spectrogram.shape[1]))
    span = spectrogram.shape[1]
    f_centre = spectrogram.shape[1] / 2

if freq_min1 is None:
    freq_min1 = ff.min()
if freq_middle is None:
    freq_middle = ff.mean()
if freq_max2 is None:
    freq_max2 = ff.max()

""" Load files from ViBR """
bitmap = loadtxt(savedir+u'\\Rectangle 1\\Bitmap.txt')
X_ini = loadtxt(savedir+u'\\Data_files\\%s_rectangle_1_raw_X.txt'%corres)
Y_ini = loadtxt(savedir+u'\\Data_files\\%s_rectangle_1_raw_Y.txt'%corres)

""" Get NW edges from bitmap and drift corrections for bitmap """
centers_b, edges_b, widths_b, first_line_found_b, last_line_found_b = detect_NW_edges(bitmap, X_ini)
X_corrected = get_normalized_X(X_ini, centers_b, widths_b, first_line_found_b, last_line_found_b)

""" Get NW edges from 2nd RSA: we use N consecutive points above/below threshold_NW_edges"""
N = 3
signal_RSA = spectrogram[:, :100].mean(axis=1)
edges = []
left_edge = None
right_edge = len(signal_RSA)
for ii, s in enumerate(signal_RSA[:1-N]):
    val_min = min(signal_RSA[ii:ii+N])
    val_max = max(signal_RSA[ii:ii+N])
    if left_edge is None and val_min > threshold_NW_edges:
        left_edge = ii
        right_edge = None
    elif right_edge is None and val_max < threshold_NW_edges:
        right_edge = ii
        edges.append(array([left_edge, right_edge]))
        left_edge = None
edges = array(edges)

""" Get approximate end of line times """
times_edges = array([[times_RSA[edges[ii, 0]], times_RSA[edges[ii, 1]]] for ii in range(len(edges))])
times_end_of_lines = array([(times_edges[ii, 1]+times_edges[ii+1, 0])/2 for ii in range(len(times_edges)-1)] + [(times_edges[-1, 1]+times_RSA.max())/2])

""" Get frequency and time for each peak """
idx_min1 = abs(ff-freq_min1).argmin()
idx_max2 = abs(ff-freq_max2).argmin()

freqs1 = []
times1 = []
f_idxs1 = []
idxs1 = []
freqs2 = []
times2 = []
f_idxs2 = []
idxs2 = []
for ii, spectrum in enumerate(spectrogram):
    tt = times_RSA[ii]
    ss = signal_RSA[ii]

    if type(freq_middle) == type(lambda x:1):
        idx_middle = abs(ff-freq_middle(tt)).argmin()
    else:
        idx_middle = abs(ff-freq_middle).argmin()

    # Search mode mode 2.1
    if spectrum[idx_min1:idx_middle].max() - ss > threshold_peak:
        idx = idx_min1 + spectrum[idx_min1:idx_middle].argmax()
        freqs1.append(ff[idx])
        times1.append(tt)
        f_idxs1.append(idx)
        idxs1.append(ii)
    # Search mode mode 2.2
    if spectrum[idx_middle:idx_max2].max() - ss > threshold_peak:
        idx = idx_middle + spectrum[idx_middle:idx_max2].argmax()
        freqs2.append(ff[idx])
        times2.append(tt)
        f_idxs2.append(idx)
        idxs2.append(ii)
freqs1 = array(freqs1)
times1 = array(times1)
idxs1 = array(idxs1)
f_idxs1 = array(f_idxs1)
freqs2 = array(freqs2)
times2 = array(times2)
idxs2 = array(idxs2)
f_idxs2 = array(f_idxs2)

""" Get Y positions of each frequency """
# Get all positions
y_all = []
idx_line_all = []
for ii in range(len(spectrogram)):
    tt = times_RSA[ii]
    if max(times_end_of_lines) < tt:
        idx_line = times_end_of_lines.argmax()
    else:
        idx_line = min([ii for ii in range(len(times_end_of_lines)) if times_end_of_lines[ii] > tt])

    if sens == '01':
        idx = idx_line + first_line_found_b
        y = Y_ini[idx, 0]
    elif sens == '10':
        idx = last_line_found_b - idx_line
        y = Y_ini[idx, 0]

    y_all.append(y)
    idx_line_all.append(idx_line)
y_all = array(y_all)
idx_line_all = array(idx_line_all)

# For mode 2.1
ys1 = []
idx_line1 = []
for ii in range(len(freqs1)):
    tt = times1[ii]
    if max(times_end_of_lines) < tt:
        idx_line = times_end_of_lines.argmax()
    else:
        idx_line = min([ii for ii in range(len(times_end_of_lines)) if times_end_of_lines[ii] > tt])

    if sens == '01':
        idx = idx_line + first_line_found_b
        y = Y_ini[idx, 0]
    elif sens == '10':
        idx = last_line_found_b - idx_line
        y = Y_ini[idx, 0]
    ys1.append(y)
    idx_line1.append(idx)
ys1 = array(ys1)
idx_line1 = array(idx_line1)

# For mode 2.2
ys2 = []
idx_line2 = []
for ii in range(len(freqs2)):
    tt = times2[ii]
    if max(times_end_of_lines) < tt:
        idx_line = times_end_of_lines.argmax()
    else:
        idx_line = min([ii for ii in range(len(times_end_of_lines)) if times_end_of_lines[ii] > tt])

    if sens == '01':
        idx = idx_line + first_line_found_b
        y = Y_ini[idx, 0]
    elif sens == '10':
        idx = last_line_found_b - idx_line
        y = Y_ini[idx, 0]
    ys2.append(y)
    idx_line2.append(idx)
ys2 = array(ys2)
idx_line2 = array(idx_line2)

""" Images plotting """
fig_spectrogram, ax_spectrogram = subplots()
imshow(spectrogram, interpolation='none', aspect='auto', extent=[ff.min()/1e3, ff.max()/1e3, times_RSA.max(), times_RSA.min()])
colorbar()
plot(freqs1/1e3, times1, 'xr')
plot(freqs2/1e3, times2, '+k')
xlabel('Frequency (kHz)')
ylabel('Time (s)')
title('Spectrogram')
tight_layout()

fig_y_vs_f, ax_y_vs_f = subplots()
plot(ys1, freqs1/1e3, 'o', label='Mode 2.1')
plot(ys2, freqs2/1e3, 'o', label='Mode 2.2')
xlabel('Longitudinal position')
ylabel('Frequency (kHz)')
legend()
title('Frequencies')
tight_layout()

fig_freq1, ax_freq1 = subplots()
plot(ys1, freqs1/1e3, 'o', label='Mode 2.1')
xlabel('Longitudinal position')
ylabel('Frequency (kHz)')
legend()
title('Frequency Mode 2.1')
tight_layout()

fig_freq2, ax_freq2 = subplots()
plot(ys2, freqs2/1e3, 'o', label='Mode 2.2')
xlabel('Longitudinal position')
ylabel('Frequency (kHz)')
legend()
title('Frequency Mode 2.2')
tight_layout()

fig_positions, ax_positions = subplots()
pcolor(X_ini, Y_ini, bitmap, cmap='gray')
plot(centers_b, Y_ini[:, 0], 'r')
if len(idx_line1) > 0:
    plot(centers_b[idx_line1], ys1, 'ob', ms=9, label='Mode 2.1')
if len(idx_line2) > 0:
    plot(centers_b[idx_line2], ys2, 'og', label='Mode 2.2')
legend()
xlabel('X position (V)')
ylabel('Y position (V)')
title('Longitudinal positions at which peaks where found\n(no lateral position yet)')
tight_layout()

show()

""" Save figures """
if savefigs:
    figsdir = savedir+'\\Figures\\Parallel_acquisition'
    if not os.path.isdir(figsdir):
        os.makedirs(figsdir)
        print('Directory %s created.'%figsdir)

    fig_spectrogram.savefig(figsdir+'\\Spectrogram')
    fig_y_vs_f.savefig(figsdir+'\\Frequencies')
    fig_freq1.savefig(figsdir+'\\Frequency1')
    fig_freq2.savefig(figsdir+'\\Frequency2')
    fig_positions.savefig(figsdir+'\\Positions')

    print('Figures saved.')

""" Save data """
if save_data:
    datadir = savedir + '\\Data_files\\Parallel_acquisition'
    if not os.path.isdir(datadir):
        os.makedirs(datadir)
        print('Directory %s created.'%datadir)

    save(datadir+'\\freqs1.npy', freqs1)
    save(datadir+'\\times1.npy', times1)
    save(datadir+'\\idxs1.npy', idxs1)
    save(datadir+'\\f_idxs1.npy', f_idxs1)
    save(datadir+'\\ys1.npy', ys1)

    save(datadir+'\\freqs2.npy', freqs2)
    save(datadir+'\\times2.npy', times2)
    save(datadir+'\\idxs2.npy', idxs2)
    save(datadir+'\\f_idxs2.npy', f_idxs2)
    save(datadir+'\\ys2.npy', ys2)

    save(datadir+'\\y_all.npy', y_all)

    print('Array files saved.')

