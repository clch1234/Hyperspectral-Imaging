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

fil = '4'
savefigs = True

dat = '20211126'
batch = 'Fil %s manuel'%fil

df = 500 # Margin left and right of the maximum of the spectrum for peak detection (in Hz)
generalWidth = 5 # Guess for Gamma/2pi (in Hz)

dirname = 'D:\\Documents\\Boulot\\Grenoble\\Data\\%s\\%s'%(dat, batch)
figsdir = dirname+'\\Figures'
if savefigs and not os.path.isdir(figsdir):
    os.mkdir(figsdir)
    print('Figures directory created :', figsdir)

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

""" Load numerical values """
key = 'Fil %s'%fil
length = lengths[key]
diameter = diameters[key]
step = steps[key]

""" Do the fits """
popts = []
positions = []
_, _, filenames = next(os.walk(dirname), (None, None, [])) # Get all files in dirname, excluding sub-directories
for ii, filename in enumerate(sorted(filenames)):
    print()
    print(filename)

    position = int(filename.split('.')[0].split(' ')[1])
    position_norm = 1 - (position-1)*step/length
    positions.append(position_norm)

    freqs, spectrum = get_freqs_data(filename, Trace=1, unit='Hz', dirname=dirname)
    omegas = freqs * 2 * pi
    print('Data loaded,', filename)

    maxS = max(abs(spectrum))
    spectrumLoc = abs(spectrum) / maxS
    minL = argmin(spectrumLoc)

    startValues = [10**(spectrum[-len(spectrum)//4:].mean()/20)] # Initial guess for fit, with a guess for the offset
    bounds_inf = [0] # Min bounds for the fit, with value for the offset
    bounds_sup = [inf] # Max bounds for the fit, with value for the offset

    omega1 = omegas[minL] # Guess for resonant frequency
    omega_inf1 = omega1 - 2*pi*df # Min  for resonant frequency
    omega_sup1 = omega1 + 2*pi*df # Max bound for resonant frequency

    guess_a1 = 10**(max(spectrum)/20) # Guess for peak amplitude (in linear units)

    startValues += [omega1, generalWidth, guess_a1]
    bounds_inf += [omega_inf1, 0, 0] # Min bounds for Gamma and amplitude are 0
    bounds_sup += [omega_sup1, inf, inf] # Max bounds for Gamma and amplitude are inf

    try:
        # Do the fit 
        print('Fitting', filename)
        res = least_squares(res_multi_lorentz, startValues, args=(omegas, spectrum), bounds=(bounds_inf, bounds_sup))
        popt = res.x
        popts.append(popt)
        fit = array([multi_lorentz(2*pi*f, popt) for f in freqs])
        spectrumLoc = [abs(s - multi_lorentz(o, popt)) / maxS for o, s in zip(omegas, spectrum)]
        print('Done.')

    except ValueError as er:
        # In case something went wrong
        print(er)
        for jj in range(len(startValues)):
            print(bounds_inf[jj], '\\', startValues[jj], '\\', bounds_sup[jj])
        errors.append(ii)
        break

    """ Plot the fit """
    fit = array([multi_lorentz(2*pi*f, popt) for f in freqs])
    offset, omega0, gamma, ampl = popt

    fig, ax = subplots()
    ax.plot(freqs, spectrum, label='Data')
    ax.plot(freqs, fit, label='$\Omega_0 = %i$\n$\Gamma_0 = %.2f$\n$A = %.3f$'%(omega0, gamma, ampl))
    xlabel('Frequency (Hz)')
    ylabel('dB')
    title(filename.split('.')[0].title())
    legend()
    if savefigs:
        fig.savefig(figsdir+'\\'+filename.split('.')[0])

        x1, x2 = ax.get_xlim()
        ax.set_xlim(omega0/2/pi-max(df, 3*gamma/2/pi), omega0/2/pi+max(df, 3*gamma/2/pi))
        fig.savefig(figsdir+'\\Zoom '+filename.split('.')[0])
        ax.set_xlim(x1, x2)
        print('Figure saved,', filename.split('.')[0])

""" Plot figures """
fig_freqs, ax_freqs = subplots()
plot(positions, [popt[1]/2/pi for popt in popts], 'o')
xlabel('Normalized position')
ylabel('Frequency (Hz)')
title('Frequency')

fig_gamma, ax_gamma = subplots()
plot(positions, [popt[2]/2/pi for popt in popts], 'o')
xlabel('Normalized position')
ylabel('$\Gamma/2\pi$ (Hz)')
title('Damping rate')

fig_ampl, ax_ampl = subplots()
plot(positions, [popt[3] for popt in popts], 'o')
xlabel('Normalized position')
ylabel('Amplitude (arb.)')
title('Peak amplitude')

if savefigs:
    fig_freqs.savefig(figsdir+'\\Frequencies')
    fig_gamma.savefig(figsdir+'\\Gamma')
    fig_ampl.savefig(figsdir+'\\Amplitude')
    print('Figures saved')

show()
