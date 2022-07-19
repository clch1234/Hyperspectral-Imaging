# -*- coding: cp1252 -*-
"""
Script for fitting the results the parallel acquisition as lorentzian peaks.

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
import os
from scipy.optimize import least_squares
import importlib.util
from functions import *
import parameters_parallel_acquisition as ppa

nw = '843'
fil = '6'
spec = '17'
dat = '20220622'
key = 'NW_'+nw+'_fil_'+fil+'_spec_'+spec.zfill(4)

savefigs = True
savefiles = True

downsample = 1
df = 10e3
generalWidth = 2*pi*1000

dirname = "D:\\Documents\\Boulot\\Grenoble\\Data\\%s"%dat
filename_oscillo = 'NW'+nw+'_fil'+fil+'_spec'+spec+'_oscillo.csv'
filename_RSA = 'NW'+nw+'_fil_'+fil+'_spectrogram_'+spec.zfill(4)+'.npy'
filename_times_RSA = 'NW'+nw+'_fil_'+fil+'_times_spectrogram_'+spec.zfill(4)+'.npy'
filename_freqs_RSA = 'NW'+nw+'_fil_'+fil+'_frequencies_spectrogram_'+spec.zfill(4)+'.txt'

""" Get corrrespondence with ViBR files """
corres = ppa.correspondences[key]
rectangle = '1'
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
spectrogram = load(savedir+'\\'+filename_RSA.split('.')[0]+'_without_duplicates'+'.npy')
times_RSA = load(savedir+'\\'+filename_times_RSA.split('.')[0]+'_without_duplicates'+'.npy')
print('Data RSA without duplicates loaded')

with open(dirname + '\\' + filename_freqs_RSA, 'r') as ff:
    lines = ff.readlines()

f_centre = float(lines[0].split(' = ')[1])
span = float(lines[1].split(' = ')[1])
ff = linspace(f_centre-span/2, f_centre+span/2, len(spectrogram[0]))
omegas = 2*pi*ff

if freq_min1 is None:
    freq_min1 = ff.min()
if freq_middle is None:
    freq_middle = ff.mean()
if freq_max2 is None:
    freq_max2 = ff.max()

print('Frequencies loaded')

""" Load preliminary results """
datadir = savedir + '\\Data_files\\Parallel_acquisition'

f1s_ini = load(datadir+'\\freqs1.npy', allow_pickle=True)
ys1 = load(datadir+'\\ys1.npy', allow_pickle=True)
times1 = load(datadir+'\\times1.npy', allow_pickle=True)
idxs1 = load(datadir+'\\idxs1.npy', allow_pickle=True)
f_idxs1 = load(datadir+'\\f_idxs1.npy', allow_pickle=True)

f2s_ini = load(datadir+'\\freqs2.npy', allow_pickle=True)
ys2 = load(datadir+'\\ys2.npy', allow_pickle=True)
times2 = load(datadir+'\\times2.npy', allow_pickle=True)
idxs2 = load(datadir+'\\idxs2.npy', allow_pickle=True)
f_idxs2 = load(datadir+'\\f_idxs2.npy', allow_pickle=True)

y_all = load(datadir+'\\y_all.npy', allow_pickle=True)

""" Normalize y by length """
# Load X, Y and bitmap
dir_bitmap = savedir + u'\\Rectangle %s'%rectangle
bitmap = loadtxt(dir_bitmap+u'\\Bitmap.txt')
dir_coordinates = savedir + u'\\Data_files'
X_ini = loadtxt(dir_coordinates+u'\\%s_rectangle_1_raw_X.txt'%corres)
Y_ini = loadtxt(dir_coordinates+u'\\%s_rectangle_1_raw_Y.txt'%corres)

# Detection of the edges of the NW and left-right separation
centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(bitmap, X_ini)

# Deduce basis position and NW length
length = abs(Y_ini[first_line_found, 0] - Y_ini[last_line_found, 0])
basis = Y_ini[first_line_found, 0]

ys1 = (ys1 - basis) / length
ys2 = (ys2 - basis) / length
y_all = (y_all - basis) / length


""" Define fit functions """
def lorentzian(omega, omega0, gamma, a):
    return a * omega0**2 * gamma**2 / ((omega**2 - omega0**2)**2 + omega**2*gamma**2)

def multi_lorentz(x, params):
    off = params[0]
    paramsRest = params[1:]
    assert not (len(paramsRest) % 3)
    return 10*log10(abs(off + sum([lorentzian(x, *paramsRest[ii:ii+3]) for ii in range(0, len(paramsRest), 3)]))**2)

def res_multi_lorentz(params, xData, yData):
    diff = [multi_lorentz(x, params) - y for x, y in zip(xData, yData)]
    return diff

""" Do the fits """
print('Starting fits. Estimated number of spectrums to fit :',
      max(len(f1s_ini), len(f2s_ini)))
success = []
errors = []
no_peaks = []
popts = []
guesses = []
y_fit = []

for ii in range(0, len(spectrogram), downsample):
    spectrum = spectrogram[ii]
    if ii%100 == 0:
        print('Spectrum %i/%i : sucess %i - errors %i - no peaks %i'%(ii+1,
                                                                      spectrogram.shape[0],
                                                                      len(success),
                                                                      len(errors),
                                                                      len(no_peaks)))

    nb_peaks = 0
    if ii in idxs1:
        guess_f1 = f1s_ini[list(idxs1).index(ii)]
        nb_peaks += 1
    if ii in idxs2:
        guess_f2 = f2s_ini[list(idxs2).index(ii)]
        nb_peaks += 1

    if nb_peaks == 0:
        no_peaks.append(ii)
    else:
        maxS = max(abs(spectrum))
        spectrumLoc = abs(spectrum) / maxS
        startValues = [10**(spectrum[:len(spectrum)//4].mean()/20)]
        bounds_inf = [0]
        bounds_sup = [inf]

        minL = argmin(spectrumLoc)
        minS = spectrum[minL]
        omega1 = omegas[minL]
        omega_inf1 = omega1 - 2*pi*df
        omega_sup1 = omega1 + 2*pi*df
        guess_a1 = 10**(max(spectrum)/20) # abs(minY - max(yDataLoc))
        startValues += [omega1, generalWidth, guess_a1]
        bounds_inf += [omega_inf1, 0, 0]
        bounds_sup += [omega_sup1, inf, inf]
        if nb_peaks == 2:
            if abs(startValues[1]/2/pi - guess_f1) < abs(startValues[1]/2/pi - guess_f2):
                idx = abs(ff - guess_f2).argmin()
                omega2 = omegas[idx]
                omega_inf2 = omega2 - 2*pi*df
                omega_sup2 = omega2 + 2*pi*df
                guess_a2 = 10**(spectrum[idx]/20)
            else:
                idx = abs(ff - guess_f1).argmin()
                omega2 = omegas[idx]
                omega_inf2 = omega2 - 2*pi*df
                omega_sup2 = omega2 + 2*pi*df
                guess_a2 = 10**(spectrum[idx]/20)
            startValues += [omega2, generalWidth, guess_a2]
            bounds_inf += [omega_inf2, 0, 0]
            bounds_sup += [omega_sup2, inf, inf]

        try:
            res = least_squares(res_multi_lorentz, startValues, args=(omegas, spectrum), bounds=(bounds_inf, bounds_sup))
            popt = res.x
            fit = array([multi_lorentz(2*pi*f, popt) for f in ff])
            spectrumLoc = [abs(s - multi_lorentz(o, popt)) / maxS for o, s in zip(omegas, spectrum)]

        except ValueError as er:
            print(er)
            for jj in range(len(startValues)):
                print(bounds_inf[jj], '\\', startValues[jj], '\\', bounds_sup[jj])
            errors.append(ii)
            break

        # If 1 peak detected : add nan to fill popt, and discriminate between f1 and f2
        if len(popt) < 7:
            if ii in idxs2:
                # The peak we have is f2
                popt = [popt[0], nan, nan, nan, popt[1], popt[2], popt[3]]
            elif ii in idxs1:
                # The peak we have is f1
                popt = [popt[0], popt[1], popt[2], popt[3], nan, nan, nan]
            else:
                # Arbitrarly keep it as f1
                popt = [popt[0], popt[1], popt[2], popt[3], nan, nan, nan]

        success.append(ii)
        popts.append(popt)
        guesses.append(startValues)
        y_fit.append(y_all[ii])

omegas1, omegas2, gammas1, gammas2, ampls1, ampls2, offsets = [nan]*len(spectrogram), [nan]*len(spectrogram), [nan]*len(spectrogram), [nan]*len(spectrogram), [nan]*len(spectrogram), [nan]*len(spectrogram), [nan]*len(spectrogram)
for nn, popt in enumerate(popts):
    offsets[success[nn]] = popt[0]

    if not isnan(popt[1]) and not isnan(popt[4]):
        if popt[1] < popt[4]:
            omegas1[success[nn]] = popt[1]
            gammas1[success[nn]] = popt[2]
            ampls1[success[nn]] = popt[3]

            omegas2[success[nn]] = popt[4]
            gammas2[success[nn]] = popt[5]
            ampls2[success[nn]] = popt[6]
        else:
            omegas1[success[nn]] = popt[4]
            gammas1[success[nn]] = popt[5]
            ampls1[success[nn]] = popt[6]

            omegas2[success[nn]] = popt[1]
            gammas2[success[nn]] = popt[2]
            ampls2[success[nn]] = popt[3]
    else:
        omegas1[success[nn]] = popt[1]
        gammas1[success[nn]] = popt[2]
        ampls1[success[nn]] = popt[3]

        omegas2[success[nn]] = popt[4]
        gammas2[success[nn]] = popt[5]
        ampls2[success[nn]] = popt[6]

omegas1, omegas2, gammas1, gammas2, ampls1, ampls2, offsets = array(omegas1), array(omegas2), array(gammas1), array(gammas2), array(ampls1), array(ampls2), array(offsets)
f1s = omegas1/2/pi
f2s = omegas2/2/pi
y_fit = array(y_fit)

""" Plot figures """
fig_freqs, ax_freqs = subplots()
plot(y_all, f1s, 'o', label='Mode 1')
plot(y_all, f2s, 'o', label='Mode 2')
legend()
xlabel('Y position')
ylabel('Frequency (Hz)')
tight_layout()

fig_gammas, ax_gammas = subplots()
plot(y_all, gammas1/2/pi, 'o', label='Mode 1')
plot(y_all, gammas2/2/pi, 'o', label='Mode 2')
legend()
xlabel('Y position')
ylabel('$\Gamma / 2pi$ (Hz)')
tight_layout()

fig_ampls, ax_ampls = subplots()
semilogy(y_all, ampls1, 'o', label='Mode 1')
semilogy(y_all, ampls2, 'o', label='Mode 2')
legend()
xlabel('Y position')
ylabel('Amplitude')
tight_layout()

fig_off, ax_of = subplots()
plot(y_all, offsets, 'o')
xlabel('Y position')
ylabel('Offset')
tight_layout()

show()

""" Save figures and files """
if savefigs:
    print('Saving figures...')
    figsdir = savedir+u'\\Figures\\Parallel_acquisition'
    if not os.path.isdir(figsdir):
        os.mkdir(figsdir)

    fig_freqs.savefig(figsdir+'\\Fit_Frequencies')
    fig_gammas.savefig(figsdir+'\\Fit_Gammas')
    fig_ampls.savefig(figsdir+'\\Fit_amplitudes')
    fig_off.savefig(figsdir+'\\Fit_offsets')

    print('Figures saved')


if savefiles:
    print('Saving data files...')
    filesdir = datadir+'\\Fits'
    if not os.path.isdir(filesdir):
        os.mkdir(filesdir)

    files = (f1s, f2s, gammas1, gammas2, ampls1, ampls2, offsets, success)
    names = ('Fit_f1s', 'Fit_f2s', 'Fit_gammas1', 'Fit_gammas2', 'Fit_ampls1', 'Fit_ampls2', 'Fit_offsets', 'Fit_success')

    for ii in range(len(files)):
        name = names[ii]
        savetxt(filesdir+'\\'+name+'.txt',
                files[ii])
    print('Data files saved')
