# -*- coding: cp1252 -*-
"""
Functions for fitting the spectrums of vibrations of a nanowire with lorentzian shape

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
#from scipy.optimize import leastsq
from scipy.optimize import least_squares
import os
#from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import datetime
import imp
from data_RSA_new import *
from data_selector import DataSelector

import warnings
warnings.filterwarnings("ignore")

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

def fit_spectrums_batch(**kwargs):
    # Data selector
    ds = DataSelector(description="Select a whole acquisition")
    ds.select_directory()

    batch = ds.directory.split('/')[-1]
    dat = ds.directory.split('/')[-2]
    workdir = '/'.join(ds.directory.split('/')[:-3])
    datadir = '/'.join(ds.directory.split('/')[:-2])

    sens = ds.sens_var.get()
    if not sens in ('10', '01'):
        raise ValueError('Sens should be either "10" or "01", %s received.'%sens)

    AR = ds.AR_var.get()
    if AR == 'None':
        AR = None
    if not AR in ('A', 'R', None, 'both'):
        raise ValueError('AR should be either "A", "R" or None, %s received.'%aller_retour)
    if AR == 'both':
        ARs = ('A', 'R')
    else:
        ARs = (AR, )

    savefigs = ds.savefigs_var.get()
    savefiles = ds.savefiles_var.get()

    kwargs['dat'] = dat
    kwargs['batch'] = batch
    kwargs['workdir'] = workdir
    kwargs['datadir'] = datadir
    kwargs['savefigs'] = savefigs
    kwargs['savefiles'] = savefiles

    dirname = workdir+u'\\Data\\%s\\%s'%(dat, batch)

    rectangles = [d.split(' ')[1] for d in os.listdir(dirname) if 'Rectangle' in d or 'Réctangle' in d]

    results = []
    for rectangle in rectangles:
        for ii, aller_retour in enumerate(ARs):
            kwargs['aller_retour'] = aller_retour
            if ii == 0:
                kwargs['sens'] = sens
            elif ii == 1:
                # Flip the direction for 2nd round
                kwargs['sens'] = '01' if sens == '10' else '10'

            tu = fit_spectrums(rectangle, **kwargs)
            results.append(tu)
    return results

def fit_spectrums(rectangle, **kwargs):
    time_start = datetime.datetime.now()

    # Read keayword arguments
    dat = kwargs['dat']# if 'dat' in kwargs else '20210519'
    batch = kwargs['batch']# if 'batch' in kwargs else 'Test 2 (nuit)'
    workdir = kwargs['workdir']# if 'workdir' in kwargs else 'D:/Documents/Boulot/Scripts/'
    datadir = kwargs['datadir']
    generalWidth = kwargs['guess_gamma'] if 'guess_gamma' in kwargs else 2*pi*1000 # Guess Gamma
    downsample = kwargs['downsample'] if 'downsample' in kwargs else 1 # Down sampling factor
    AR = kwargs['aller_retour']# if 'aller_retour' in kwargs else None
    sens = kwargs['sens']# if 'sens' in kwargs else '10'

    savefigs = kwargs['savefigs']# if 'savefigs' in kwargs else True
    savefiles = kwargs['savefiles']# if 'savefiles' in kwargs else True
    old_position_correction = kwargs['old_position_correction'] if 'sition_correction' in kwargs else False
    print('Starting script : rectangle = %s'%rectangle)
    for k, v in kwargs.items():
        print('\t', k, '=', v)

    # Load acquisition parameters
    print('Loading parameters...')
    try:
        """
        Version for python 3.x
        """
        spec = importlib.util.spec_from_file_location('parameters', workdir+u'\\Data\\%s\\%s'%(dat, batch)+'\\parameters.py')
        parameters = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(parameters)
        """
        Version for python 2.7
        """
        #parameters = imp.load_source(workdir+u'\\Data\\%s\\%s'%(dat, batch)+'\\parameters.py', 'parameters')
    except FileNotFoundError: # For python 3.x
    #except IOError: # For python 2.7
        class DummyParams(object):
            def __init__(self):
                self.dd_middle = {}
                self.dd_threshold_peak = {}
                self.dd_df = {}
                self.dd_last_quarter_for_threshold = {}
                self.dd_required_points_above_threshold = {}
        parameters = DummyParams()
        print('\n!!!!! No parameters file found in %s\\Data\\%s\\%s !!!!!\n'%(workdir, dat, batch))

    key = 'rectangle_%s'%rectangle
    middle = parameters.dd_middle[key] if key in parameters.dd_middle else lambda ii:None
    threshold_peak = parameters.dd_threshold_peak[key] if key in parameters.dd_threshold_peak else 3.#-40
    last_quarter_for_threshold = parameters.dd_last_quarter_for_threshold[key] if key in parameters.dd_last_quarter_for_threshold else False

    df = parameters.dd_df[key] if key in parameters.dd_df else 4e3
    print('Parameters loaded.')

    # Prepare directory for saving files
    if savefiles:
        filesdir = datadir+u'\\%s\\%s\\Data_files'%(dat, batch)
        if not os.path.isdir(filesdir):
            os.mkdir(filesdir)
            print('Savefiles detected, directory created')

    # Load files
    print('Loading files...')
    dirname = datadir+u'\\%s\\%s\\Réctangle %s\\'%(dat, batch, rectangle) if dat < '20211213' else datadir+u'\\%s\\%s\\Rectangle %s\\'%(dat, batch, rectangle)

    filename = u'Paramétres spéctres.txt' if dat < '20211212' else u'Parametres spectres.txt'
    params = load_params(dirname, filename)
    numpoints = int(params['Nombre de point par spectre']) if 'Nombre de point par spectre' in params else 4001

    if AR is None:
        spectrogram = loadtxt(dirname+'\\Spectres.txt')
        bitmap = loadtxt(dirname+u'\\Bitmap.txt')

        X_ini = loadtxt(filesdir+'\\'+batch+'_rectangle_'+rectangle+'_raw_X.txt')
        Y_ini = loadtxt(filesdir+'\\'+batch+'_rectangle_'+rectangle+'_raw_Y.txt')
        F1 = loadtxt(filesdir+'\\'+batch+'_rectangle_'+rectangle+'_F1.txt')
        F2 = loadtxt(filesdir+'\\'+batch+'_rectangle_'+rectangle+'_F2.txt')
        F1_ravel, F2_ravel = F1.ravel(), F2.ravel()
    elif AR == 'A':
        if dat < '20211212':
            spectrogram = loadtxt(dirname+'\\Spectres A.txt')
            bitmap = loadtxt(dirname+u'\\Bitmap A.txt')
        else:
            spectrogram = loadtxt(dirname+'\\Aller\\Spectres.txt')
            if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                bitmap = loadtxt(dirname+'\\Aller\\Bitmap.txt')
            elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                bitmap = loadtxt(dirname+'\\Aller\\Bitmap')
            else:
                raise ValueError('No bitmap file found')

        X_ini = loadtxt(filesdir+'\\A_'+batch+'_rectangle_'+rectangle+'_raw_X.txt')
        Y_ini = loadtxt(filesdir+'\\A_'+batch+'_rectangle_'+rectangle+'_raw_Y.txt')
        F1 = loadtxt(filesdir+'\\A_'+batch+'_rectangle_'+rectangle+'_F1.txt')
        F2 = loadtxt(filesdir+'\\A_'+batch+'_rectangle_'+rectangle+'_F2.txt')
        F1_ravel, F2_ravel = F1.ravel(), F2.ravel()
    elif AR == 'R':
        if dat < '20211212':
            spectrogram = loadtxt(dirname+'\\Spectres R.txt')
            bitmap = loadtxt(dirname+u'\\Bitmap R.txt')
        else:
            spectrogram = loadtxt(dirname+'\\Retour\\Spectres.txt')
            if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                bitmap = loadtxt(dirname+'\\Retour\\Bitmap.txt')
            elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                bitmap = loadtxt(dirname+'\\Retour\\Bitmap')
            else:
                raise ValueError('No bitmap file found')

        X_ini = loadtxt(filesdir+'\\R_'+batch+'_rectangle_'+rectangle+'_raw_X.txt')
        Y_ini = loadtxt(filesdir+'\\R_'+batch+'_rectangle_'+rectangle+'_raw_Y.txt')
        F1 = loadtxt(filesdir+'\\R_'+batch+'_rectangle_'+rectangle+'_F1.txt')
        F2 = loadtxt(filesdir+'\\R_'+batch+'_rectangle_'+rectangle+'_F2.txt')
        F1_ravel, F2_ravel = F1.ravel(), F2.ravel()
    else:
        raise ValueError('AR should be one of "A", "R" or None !')

    F = get_corrected_freqs(dirname, params, old_corrections=old_position_correction, AR=AR, sens=sens)
    print('Files loaded.')

    # Load shape and check for interruption in acquisition
    xshape, yshape = [int(float(st)) for st in params[u'Résolution (pixels)'].split('x')] if u'Résolution (pixels)' in params \
                     else [int(float(st)) for st in params[u'Resolution (pixels)'].split('x')]
    if xshape * yshape > spectrogram.shape[0]:
        shape_ini = (xshape, y_shape)
        yshape = int(spectrogram.shape[0] / xshape)
        #X, Y, X_ini, Y_ini, X_corr, Y_corr = X[:yshape], Y[:yshape], X_ini[:yshape], Y_ini[:yshape], X_corr[:yshape], Y_corr[:yshape]
        print('Acquisition interruption detected ! New shape : %ix%i (from %s)'%(xshape, yshape, shape_ini))

    # Prepare directory for fit plots
    plot_fits = kwargs['plot_fits'] if 'plot_fits' in kwargs else False
    if plot_fits:
        fitsdir = dirname+'\\Fits\\'+time_start.strftime("%d-%m-%Y-%H-%M-%S")
        if AR in ('A', 'R'):
            fitsdir += ' - '+AR
        if not os.path.isdir(dirname+'\\Fits'):
            os.mkdir(dirname+'\\Fits')
        os.mkdir(fitsdir)
        print('Plot fits detected, directory created')

    # Do the fits
    print('Starting fits. Estimated number of spectrums to fit :',
          len([ii for ii in range(len(F1_ravel)) if not(isnan(F1_ravel[ii]) and isnan(F2_ravel[ii]))]))
    success = []
    errors = []
    no_peaks = []
    popts = []
    guesses = []
    #for ii, spectrum in enumerate(spectrogram):
    for ii in range(0, len(spectrogram), downsample):
        spectrum = spectrogram[ii]
        if ii%100 == 0:
            print('Spectrum %i/%i : sucess %i - errors %i - no peaks %i'%(ii+1,
                                                                          spectrogram.shape[0],
                                                                          len(success),
                                                                          len(errors),
                                                                          len(no_peaks)))

        nb_peaks = 0
        guess_f1, guess_f2 = F1_ravel[ii], F2_ravel[ii]
        if not isnan(guess_f1):
            nb_peaks += 1
        if not isnan(guess_f2):
            nb_peaks += 1

        if nb_peaks == 0:
            no_peaks.append(ii)
        else:
            ff = F[ii]
            omegas = 2*pi*ff

            maxS = max(abs(spectrum))
            spectrumLoc = abs(spectrum) / maxS
            startValues = [10**(spectrum[-len(spectrum)//4:].mean()/20)] if last_quarter_for_threshold else [10**(spectrum[:len(spectrum)//4].mean()/20)]
            bounds_inf = [0]
            bounds_sup = [inf]
            counter = 0

            """
            diffs = []
            while max(spectrumLoc) - min(spectrumLoc) > .3:
                counter += 1
                if counter > nb_peaks: ### max nb of peaks, emergency break to avoid infinite loop
                    #popts.append([])
                    print('Break : counter =', counter, ' and nb_peaks =', nb_peaks)
                    break
            """
            #while counter < nb_peaks:
            #counter += 1
            #print('Spectrum %i, %i peak(s)'%(ii, counter))
            #print('Spectrum %i, %i peak(s)'%(ii, nb_peaks))
            #if nb_peaks == 1 or len(startValues) < 4:
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
            #popt, ier = leastsq(res_multi_lorentz, startValues, args=(xData, yData))
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
            # Fin boucle while

            # If 1 peak detected : add nan to fill popt, and discriminate between f1 and f2
            if len(popt) < 7:
                if isnan(guess_f1):
                    # The peak we have is f2
                    popt = [popt[0], nan, nan, nan, popt[1], popt[2], popt[3]]
                elif isnan(guess_f2):
                    # The peak we have is f1
                    popt = [popt[0], popt[1], popt[2], popt[3], nan, nan, nan]
                else:
                    # Arbitrarly keep it as f1
                    popt = [popt[0], popt[1], popt[2], popt[3], nan, nan, nan]

            success.append(ii)
            popts.append(popt)
            guesses.append(startValues)

            if plot_fits and ii == success[-1]:
                fig, ax = subplots()
                plot(ff/1e3, spectrum, label='Data')
                plot(ff/1e3, fit, label='$\Omega_1/2\pi = %s$\n$\Gamma_1/2\pi = %s$\n$A_1 = %s$\n$\Omega_2/2\pi = %s$\n$\Gamma_2/2\pi = %s$\n$A_2 = %s$'%(popt[1], popt[2], popt[3], popt[4], popt[5], popt[6]))
                legend()
                xlabel('Frequency (kHz)')
                ylabel('dB')
                title('Spectrum %i'%ii)

                axx = inset_axes(ax, width="30%", height ="30%", loc=2)
                axx.xaxis.set_visible(False)
                axx.yaxis.set_visible(False)
                axx.pcolor(X_ini, Y_ini, bitmap, cmap='gray', shading='nearest')
                axx.invert_yaxis()
                yy, xx = ii//xshape, ii%xshape
                axx.plot(X_ini[yy, xx], Y_ini[yy, xx], 'xr')

                fig.savefig(fitsdir+'\\spectrum_'+str(ii).zfill(5))
                close(fig)


    omegas1, omegas2, gammas1, gammas2, ampls1, ampls2, offsets = [nan]*len(spectrogram), [nan]*len(spectrogram), [nan]*len(spectrogram), [nan]*len(spectrogram), [nan]*len(spectrogram), [nan]*len(spectrogram), [nan]*len(spectrogram)
    for nn, popt in enumerate(popts):
        offsets[success[nn]] = popt[0]

        omegas1[success[nn]] = popt[1]
        gammas1[success[nn]] = popt[2]
        ampls1[success[nn]] = popt[3]

        omegas2[success[nn]] = popt[4]
        gammas2[success[nn]] = popt[5]
        ampls2[success[nn]] = popt[6]
    omegas1, omegas2, gammas1, gammas2, ampls1, ampls2, offsets = array(omegas1), array(omegas2), array(gammas1), array(gammas2), array(ampls1), array(ampls2), array(offsets)
    f1s = omegas1/2/pi
    f2s = omegas2/2/pi

    fig_freqs, ax_freqs = subplots()
    plot(Y_ini.ravel(), f1s, 'o', label='Mode 1')
    plot(Y_ini.ravel(), f2s, 'o', label='Mode 2')
    legend()
    xlabel('Y position')
    ylabel('Frequency (Hz)')
    tight_layout()

    fig_gammas, ax_gammas = subplots()
    plot(Y_ini.ravel(), gammas1/2/pi, 'o', label='Mode 1')
    plot(Y_ini.ravel(), gammas2/2/pi, 'o', label='Mode 2')
    legend()
    xlabel('Y position')
    ylabel('$\Gamma / 2pi$ (Hz)')
    tight_layout()

    fig_ampls, ax_ampls = subplots()
    semilogy(Y_ini.ravel(), ampls1, 'o', label='Mode 1')
    semilogy(Y_ini.ravel(), ampls2, 'o', label='Mode 2')
    legend()
    xlabel('Y position')
    ylabel('Amplitude')
    tight_layout()

    fig_off, ax_of = subplots()
    plot(Y_ini.ravel(), offsets, 'o')
    xlabel('Y position')
    ylabel('Offset')
    tight_layout()

    if savefigs:
        print('Saving figures...')
        figs = [fig_freqs, fig_gammas, fig_ampls, fig_off]
        names = ['Fit_Frequencies', 'Fit_Gammas', 'Fit_amplitudes', 'Fit_offsets']
        figsdir = datadir+u'\\%s\\%s\\Figures'%(dat, batch)
        if not os.path.isdir(figsdir):
            os.mkdir(figsdir)

        for ii in range(len(figs)):
            name = batch+'_rectangle_'+rectangle+'_'+names[ii]
            if AR in ('A', 'R'):
                name = AR + '_' + name
            figs[ii].savefig(figsdir+'\\'+name)
        print('Figures saved')

    if savefiles:
        print('Saving data files...')
        files = (f1s, f2s, gammas1, gammas2, ampls1, ampls2, offsets, success)
        names = ('Fit_f1s', 'Fit_f2s', 'Fit_gammas1', 'Fit_gammas2', 'Fit_ampls1', 'Fit_ampls2', 'Fit_offsets', 'Fit_success')

        for ii in range(len(files)):
            name = batch+'_rectangle_'+rectangle+'_'+names[ii]+'.txt'
            if AR in ('A', 'R'):
                name = AR + '_' + name
            savetxt(filesdir+'\\'+name,
                    files[ii])
        print('Data files saved')

    return f1s, f2s, gammas1, gammas2, ampls1, ampls2, offsets, popts, guesses, success, X_ini, Y_ini
