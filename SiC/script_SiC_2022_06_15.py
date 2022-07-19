# -*- coding: cp1252 -*-
"""
Functions for fitting the spectrums of vibrations of a SiC nanowire with lorentzian shape
We assume that we only have one peak (one direction of vibration visible)

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
from scipy.optimize import least_squares
import os
from donnees_rsa import *
from numerical_values import *

fil = '9'
savefigs = True # Save figures
savefiles = True # Save results of fits
plot_only = False # Use saved results of fits instead of fitting again
if plot_only:
    savefiles = False
    print('Savefiles set to False due to plot_only = True')

##dat = '20211126'
##batch = 'Fil %s manuel'%fil
dat = '20220620'
batch = 'Fil %s Manuel'%fil

df = 100 # Margin left and right of the maximum of the spectrum for peak detection (in Hz)
generalWidth = 5 # Guess for Gamma/2pi (in Hz)

dirname = 'D:\\Documents\\Boulot\\Grenoble\\Data\\%s\\%s'%(dat, batch)
figsdir = dirname+'\\Figures'
if savefigs and not os.path.isdir(figsdir):
    os.mkdir(figsdir)
    print('Figures directory created :', figsdir)
filesdir = dirname + '\\Fits Results'
if savefiles and not os.path.isdir(filesdir):
    os.mkdir(filesdir)
    print('Files directory created :', filesdir)

""" Define functions """
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

""" Load noise background """
filename_bg  = 'Fond_moyenne_10_fois.csv'
if os.path.isfile(filename_bg):
    freqs_bg, spectrum_bg = get_freqs_data(filename_bg, Trace=1, unit='Hz', dirname=dirname)
    omegas_bg = freqs_bg * 2 * pi
    background = True
    print('Noise Background loaded,', filename_bg)
else:
    background = False

""" Load numerical values """
key = 'Fil %s'%fil
length = lengths[key]
diameter = diameters[key]
step = steps[key]
freq_min = freqs_min[key] if key in freqs_min else freqs_bg.min() if background else None
freq_max = freqs_max[key] if key in freqs_max else freqs_bg.max() if background else None
# normalized position of acquisition at the tip
position_0 = (length-positions_0[key])/length if key in positions_0 else 1

""" Load file names """
_, _, filenames = next(os.walk(dirname), (None, None, [])) # Get all files in dirname, excluding sub-directories
if background:
    filenames.remove(filename_bg)
num_position_min = min([float(filename.split('.')[0].split(' ')[1]) for filename in filenames])

if plot_only:
    """ If plot only, load the values """
    popts = loadtxt(filesdir+'\\popts.txt')
    popts_g = loadtxt(filesdir+'\\popts_g.txt')
    popts_d = loadtxt(filesdir+'\\popts_d.txt')

    positions = loadtxt(filesdir+'\\positions.txt')
    positions_g = loadtxt(filesdir+'\\positions_g.txt')
    positions_d = loadtxt(filesdir+'\\positions_d.txt')

else:
    """ Else, do the fits """
    popts = []
    positions = []
    popts_g = []
    positions_g = []
    popts_d = []
    positions_d = []

    for ii, filename in enumerate(sorted(filenames, key=lambda f: float(f.split('.')[0].split(' ')[1]))):
        print()
        print(filename)

        position = float(filename.split('.')[0].split(' ')[1])
        position_norm = position_0 - (position-num_position_min)*step/length

        for trace in (1, 2, ):
            print('Trace', trace)

            try:
                freqs, spectrum = get_freqs_data(filename, Trace=trace, unit='Hz', dirname=dirname)
                omegas = freqs * 2 * pi
                print('Data loaded,', filename)

                positions.append(position_norm)

                if background:
                    spectrum_no_bg = spectrum - spectrum_bg
                else:
                    spectrum_no_bg = spectrum - spectrum[:100].min()
                    if freq_min is None:
                        freq_min = freqs.min()
                    if freq_max is None:
                        freq_max = freqs.max()

                mask_to_fit = where((freqs >= freq_min) & (freqs <= freq_max))[0]
                freqs_to_fit = freqs[mask_to_fit]
                omegas_to_fit = 2*pi*freqs_to_fit
                spectrum_to_fit = spectrum_no_bg[mask_to_fit]

                maxS = max(abs(spectrum_to_fit))
                spectrumLoc = abs(spectrum_to_fit) / maxS
                minL = argmax(spectrumLoc)

                startValues = [10**(spectrum_to_fit[-len(spectrum)//4:].mean()/20)] # Initial guess for fit, with a guess for the offset
                bounds_inf = [0] # Min bounds for the fit, with value for the offset
                bounds_sup = [inf] # Max bounds for the fit, with value for the offset

                omega1 = omegas_to_fit[minL] # Guess for resonant frequency
                omega_inf1 = omega1 - 2*pi*df # Min  for resonant frequency
                omega_sup1 = omega1 + 2*pi*df # Max bound for resonant frequency

                guess_a1 = 10**(max(spectrum_to_fit)/20) # Guess for peak amplitude (in linear units)

                startValues += [omega1, generalWidth, guess_a1]
                bounds_inf += [omega_inf1, 0, 0] # Min bounds for Gamma and amplitude are 0
                bounds_sup += [omega_sup1, inf, inf] # Max bounds for Gamma and amplitude are inf

                try:
                    # Do the fit 
                    print('Fitting', filename)
                    res = least_squares(res_multi_lorentz, startValues, args=(omegas_to_fit, spectrum_to_fit), bounds=(bounds_inf, bounds_sup))
                    popt = res.x
                    popts.append(popt)
                    spectrumLoc = [abs(s - multi_lorentz(o, popt)) / maxS for o, s in zip(omegas_to_fit, spectrum_to_fit)]
                    print('Done.')

                except ValueError as er:
                    # In case something went wrong
                    print(er)
                    for jj in range(len(startValues)):
                        print(bounds_inf[jj], '\\', startValues[jj], '\\', bounds_sup[jj])
                    errors.append(ii)
                    break

                # Left or right
                if 'gauche' in filename:
                    positions_g.append(position_norm)
                    popts_g.append(popt)
                elif 'droite' in filename:
                    positions_d.append(position_norm)
                    popts_d.append(popt)

                """ Plot the fits """
                ff = linspace(freqs.min(), freqs.max(), int(1e6))
                fit = array([multi_lorentz(2*pi*f, popt) for f in ff])
                offset, omega0, gamma, ampl = popt

                fig, ax = subplots()
            ##    ax.plot(freqs, spectrum, label='Data')
            ##    ax.plot(freqs, spectrum_bg+fit, label='$\Omega_0 = %i$\n$\Gamma_0 = %.2f$\n$A = %.3f$'%(omega0, gamma, ampl))
                ax.plot(freqs, spectrum_no_bg, label='Data', c='C%i'%(2*trace))
                ax.plot(ff, fit, label='$f_0 = %i$\n$\Gamma_0 = %.2f$\n$A = %.3f$'%(omega0/2/pi, gamma, ampl), c='C%i'%(2*trace+1))
                xlabel('Frequency (Hz)')
                ylabel('dB')
                title(filename.split('.')[0].title())
                legend()

                if savefigs:
                    fig.savefig(figsdir+'\\%s Trace %s'%(filename.split('.')[0], trace))

                    x1, x2 = ax.get_xlim()
                    ax.set_xlim(omega0/2/pi-max(df, 3*gamma/2/pi), omega0/2/pi+max(df, 3*gamma/2/pi))
                    fig.savefig(figsdir+'\\Zoom %s Trace %s'%(filename.split('.')[0], trace))
                    ax.set_xlim(x1, x2)
                    print('Figure saved,', filename.split('.')[0])

            except ValueError:
                print('%s : Trace %i not found'%(filename, trace))

""" Save files if required """
if savefiles:
    savetxt(filesdir+'\\popts.txt', popts)
    savetxt(filesdir+'\\popts_g.txt', popts_g)
    savetxt(filesdir+'\\popts_d.txt', popts_d)

    savetxt(filesdir+'\\positions.txt', positions)
    savetxt(filesdir+'\\positions_g.txt', positions_g)
    savetxt(filesdir+'\\positions_d.txt', positions_d)

    print('Files saved')

""" Plot figures """
fig_freqs, ax_freqs = subplots()
plot(positions, [popt[1]/2/pi for popt in popts], 'o')
if len(popts_g) > 0:
    plot(positions_g, [popt[1]/2/pi for popt in popts_g], 'o', label='Left')
if len(popts_d) > 0:
    plot(positions_d, [popt[1]/2/pi for popt in popts_d], 'o', label='Right')
xlabel('Normalized position')
ylabel('Frequency (Hz)')
title('Frequency')
legend()

fig_gamma, ax_gamma = subplots()
plot(positions, [popt[2]/2/pi for popt in popts], 'o')
if len(popts_g) > 0:
    plot(positions_g, [popt[2]/2/pi for popt in popts_g], 'o', label='Left')
if len(popts_d) > 0:
    plot(positions_d, [popt[2]/2/pi for popt in popts_d], 'o', label='Right')
xlabel('Normalized position')
ylabel('$\Gamma/2\pi$ (Hz)')
title('Damping rate')
legend()

fig_ampl, ax_ampl = subplots()
semilogy(positions, [popt[3] for popt in popts], 'o')
if len(popts_g) > 0:
    semilogy(positions_g, [popt[3] for popt in popts_g], 'o', label='Left')
if len(popts_d) > 0:
    semilogy(positions_d, [popt[3] for popt in popts_d], 'o', label='Right')
xlabel('Normalized position')
ylabel('Amplitude (arb.)')
title('Peak amplitude')
legend()

""" Save figures if required """
if savefigs:
    fig_freqs.savefig(figsdir+'\\Frequencies')
    fig_gamma.savefig(figsdir+'\\Gamma')
    fig_ampl.savefig(figsdir+'\\Amplitude')
    print('Figures saved')

show()
