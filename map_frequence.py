# -*- coding: cp1252 -*-
"""
Functions for extracting and plotting the frequencies of vibration of a nanowire as a function of the SEM beam position

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
#import importlib.util
#import imp
import itertools
import datetime
from data_RSA_new import *
from data_selector import DataSelector

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

def get_value_for_threshold(spectrum, last_quarter_for_threshold):
    value_for_threshold = spectrum[-len(spectrum)//4:].max() if last_quarter_for_threshold else spectrum[:len(spectrum)//4].max()
    return value_for_threshold

def criteria_above_noise_max(spectrum, threshold_peak, required_numpoints, value_for_threshold):
    """
    Detection method for peaks :
        If more than a certain nb of CONSECUTIVE points in the spectrum are above threshold, we consider we have a peak
        Threshold is determined by the max value of the spectrum over a quarter of the data where we KNOW there is no peak,
        plus a certain value sest by threshold_peak (3 dB by default, can be the variance of the noise)

    Returns True if we detect a peak in spectrum, False otherwise
    """
    if threshold_peak == 'var':
        threshold_peak = var(spectrum[-len(spectrum)//4:]) if last_quarter_for_threshold else var(spectrum[:len(spectrum)//4])
    aux = lambda ii: min(spectrum[max(0, ii-required_numpoints//2):min(len(spectrum), ii+required_numpoints//2)])
    return len(where(array([aux(ii) - value_for_threshold for ii in range(len(spectrum))]) > threshold_peak)[0]) > 0

def criteria_above_noise_mean(spectrum, threshold_peak, required_numpoints, value_for_threshold):
    """
    Detection method for peaks :
        If more than a certain nb of CONSECUTIVE points in the spectrum are above threshold, we consider we have a peak
        Threshold is determined by the mean value of the spectrum over a quarter of the data where we KNOW there is no peak,
        plus a certain value sest by threshold_peak (3 dB by default, can be the variance of the noise)

    Returns True if we detect a peak in spectrum, False otherwise
    """
    if threshold_peak == 'var':
        threshold_peak = var(spectrum[-len(spectrum)//4:]) if last_quarter_for_threshold else var(spectrum[:len(spectrum)//4])
    aux = lambda ii: min(spectrum[max(0, ii-required_numpoints//2):min(len(spectrum), ii+required_numpoints//2)])
    return len(where(array([aux(ii) - value_for_threshold for ii in range(len(spectrum))]) > threshold_peak)[0]) > 0

def map_frequence_batch(**kwargs):
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
    aller_retour = ds.AR_var.get()
    if aller_retour == 'None':
        aller_retour = None
    if not aller_retour in ('A', 'R', None):
        raise ValueError('AR should be either "A", "R" or None, %s received.'%aller_retour)

    savefigs = ds.savefigs_var.get()
    savefiles = ds.savefiles_var.get()

    kwargs['dat'] = dat
    kwargs['batch'] = batch
    kwargs['workdir'] = workdir
    kwargs['datadir'] = datadir
    kwargs['sens'] = sens
    kwargs['aller_retour'] = aller_retour
    kwargs['savefigs'] = savefigs
    kwargs['savefiles'] = savefiles

    rectangles = [d.split(' ')[1] for d in os.listdir(ds.directory) if 'Rectangle' in d or 'Réctangle' in d]

    results = []
    for rectangle in rectangles:
        tu = map_frequence(rectangle, **kwargs)
        results.append(tu)
    return results

def map_frequence(rectangle, **kwargs):
    time_start = datetime.datetime.now()

    # Read keayword arguments
    dat = kwargs['dat']# if 'dat' in kwargs else '20210519'
    batch = kwargs['batch']# if 'batch' in kwargs else 'Test 2 (nuit)'
    workdir = kwargs['workdir']# if 'workdir' in kwargs else 'D:/Documents/Boulot/Scripts/'
    datadir = kwargs['datadir']
    sens = kwargs['sens']# if 'sens' in kwargs else '10'
    aller_retour = kwargs['aller_retour']# if 'aller_retour' in kwargs else None # aller_retour = 'A', 'R', or None

    savefigs = kwargs['savefigs']# if 'savefigs' in kwargs else True
    savefiles = kwargs['savefiles']# if 'savefiles' in kwargs else True
    old_position_correction = kwargs['old_position_correction'] if 'sition_correction' in kwargs else False

    df = kwargs['df'] if 'df' in kwargs else 4e3

    threshold_peak = kwargs['threshold_peak'] if 'threshold_peak' in kwargs else 3 # If no value is set in kwargs, we check in the parameters file
    required_numpoints = kwargs['required_points_above_threshold'] if 'required_points_above_threshold' in kwargs else 10 # If no value is set in kwargs, we check in the parameters file
    general_last_quarter_for_threshold = kwargs['last_quarter_for_threshold'] if 'last_quarter_for_threshold' in kwargs else 'auto' # If no value is set in kwargs, we check in the parameters file
    if 'criteria_for_peak_selection' in kwargs:
        if kwargs['criteria_for_peak_selection'] == 'above_noise_max':
            criteria_for_peak_selection = criteria_above_noise_max
        elif kwargs['criteria_for_peak_selection'] == 'above_noise_mean':
            criteria_for_peak_selection = criteria_above_noise_mean
        else:
            raise ValueError('Wrong criteria_for_peak_selection')
    else:
        criteria_for_peak_selection = criteria_above_noise_max
        threshold_peak = 3

    print('Starting script : rectangle = %s'%rectangle)
    for k, v in kwargs.items():
        print('\t', k, '=', v)

    # Load acquisition parameters
    print('Loading parameters...')
    try:
        """
        Version for python 3.x
        Comment if you have python 2.7
        """
        spec = importlib.util.spec_from_file_location('parameters', workdir+u'\\Data\\%s\\%s'%(dat, batch)+'\\parameters.py')
        parameters = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(parameters)
        """
        Version for python 2.7
        Comment if you have python 3.x
        """
    #parameters = imp.load_source(workdir+u'\\Data\\%s\\%s'%(dat, batch)+'\\parameters.py', 'parameters')
    except FileNotFoundError: # Use this if you have python 3.x
    #except IOError: # Use this if you have python 2.7
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
    # middle is used for the old method of peak detection, see later
    middle = parameters.dd_middle[key] if key in parameters.dd_middle and callable(parameters.dd_middle[key]) else lambda ii:None
    if key in parameters.dd_last_quarter_for_threshold:
        general_last_quarter_for_threshold = parameters.dd_last_quarter_for_threshold[key]
        print('general_last_quarter_for_threshold taken for parameter.py instead ok kwargs !')
    if key in parameters.dd_threshold_peak:
        threshold_peak = parameters.dd_threshold_peak[key]
        print('threshold_peak taken for parameter.py instead ok kwargs !')
    if key in parameters.dd_required_points_above_threshold:
        required_numpoints = parameters.dd_required_points_above_threshold[key]
        print('required_numpoints taken for parameter.py instead ok kwargs !')

    if key in parameters.dd_df:
        df = parameters.dd_df[key]
        print('df taken for parameter.py instead ok kwargs !')
    print('Parameters loaded.')

    # Load data files
    print('Loading files...')
    dirname = datadir+u'\\%s\\%s\\Réctangle %s\\'%(dat, batch, rectangle) if dat < '20211212' else datadir+u'\\%s\\%s\\Rectangle %s\\'%(dat, batch, rectangle)
    filename = u'Paramétres spéctres.txt' if dat < '20211212' else u'Parametres spectres.txt'
    params = load_params(dirname, filename)
    numpoints = int(params['Nombre de point par spectre']) if 'Nombre de point par spectre' in params else 4001

    if aller_retour is None:
        spectrogram = loadtxt(dirname+'\\Spectres.txt')
        bitmap = loadtxt(dirname+u'\\Bitmap.txt')
    else:
        if dat < '20211212':
            spectrogram = loadtxt(dirname+'\\Spectres %s.txt'%aller_retour)
            bitmap = loadtxt(dirname+'\\Bitmap %s.txt'%aller_retour)
        else:
            if aller_retour == 'A':
                spectrogram = loadtxt(dirname+'\\Aller\\Spectres.txt')
                if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                    bitmap = loadtxt(dirname+'\\Aller\\Bitmap.txt')
                elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                    bitmap = loadtxt(dirname+'\\Aller\\Bitmap')
                else:
                    raise ValueError('No bitmap file found')
            elif aller_retour == 'R':
                spectrogram = loadtxt(dirname+'\\Retour\\Spectres.txt')
                if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                    bitmap = loadtxt(dirname+'\\Retour\\Bitmap.txt')
                elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                    bitmap = loadtxt(dirname+'\\Retour\\Bitmap')
                else:
                    raise ValueError('No bitmap file found')
            else:
                raise ValueError('aller_retour should be either "A" or "R"')
    bitmap_ravel = bitmap.ravel()
    print('Files loaded.')

    # If we want to plot the spectrums, create the directory to save them
    plot_spectrums = kwargs['plot_spectrums'] if 'plot_spectrums' in kwargs else False
    if plot_spectrums:
        spectrumsdir = dirname+'\\Spectrums\\'+time_start.strftime("%d-%m-%Y-%H-%M-%S")
        if not os.path.isdir(dirname+'\\Spectrums'):
            os.mkdir(dirname+'\\Spectrums')
        os.mkdir(spectrumsdir)
        print('Plot spectrums detected, directory created')

    # Frequencies for each spectrum and indexes of spectrums
    print('Extracting all coordinates and frequencies...')
    F = get_corrected_freqs(dirname, params, old_corrections=old_position_correction, AR=aller_retour, sens=sens)
    I = repeat(range(spectrogram.shape[0]), spectrogram.shape[1]).reshape(spectrogram.shape)

    # Coordinates for each pixel of the bitmap image
    X_corr, Y_corr = get_corrected_coordinates(dirname, params, old_corrections=old_position_correction)
    if aller_retour is None:
        X_ini = loadtxt(dirname+u'\\Tensions X.txt')
        Y_ini = loadtxt(dirname+u'\\Tensions Y.txt')
    else:
        if dat < '20211212':
            X_ini = loadtxt(dirname+'\\Tensions X %s.txt'%aller_retour)
            Y_ini = loadtxt(dirname+'\\Tensions Y %s.txt'%aller_retour)
        else:
            if aller_retour == 'A' :
                if os.path.isfile(dirname+'\\Aller\\Tensions X.txt'):
                    X_ini = loadtxt(dirname+'\\Aller\\Tensions X.txt')
                    Y_ini = loadtxt(dirname+'\\Aller\\Tensions Y.txt')
                elif os.path.isfile(dirname+'\\Aller\\Tensions X'):
                    X_ini = loadtxt(dirname+'\\Aller\\Tensions X')
                    Y_ini = loadtxt(dirname+'\\Aller\\Tensions Y')
                else:
                    raise ValueError('No tensions file found')
            elif aller_retour == 'R':
                if os.path.isfile(dirname+'\\Retour\\Tensions X.txt'):
                    X_ini = loadtxt(dirname+'\\Retour\\Tensions X.txt')
                    Y_ini = loadtxt(dirname+'\\Retour\\Tensions Y.txt')
                elif os.path.isfile(dirname+'\\Retour\\Tensions X'):
                    X_ini = loadtxt(dirname+'\\Retour\\Tensions X')
                    Y_ini = loadtxt(dirname+'\\Retour\\Tensions Y')
                else:
                    raise ValueError('No tensions file found')
    X, Y = X_ini + X_corr, Y_ini + Y_corr

    # Current measured after each correction
    if not old_position_correction:
        corrections = loadtxt(dirname+'\\Corrections.txt', skiprows=1)
        if len(corrections.shape) == 1:
            num_lines_currents = array([corrections[2]])
            currents = array([corrections[5]])
        else:
            num_lines_currents = array(list(itertools.accumulate(corrections[:, 2])))
            currents = corrections[:, 5]

    xshape, yshape = [int(float(st)) for st in params[u'Résolution (pixels)'].split('x')] if dat < '20211212' else [int(float(st)) for st in params[u'Resolution (pixels)'].split('x')]
    if xshape * yshape > spectrogram.shape[0]:
        shape_ini = (xshape, yshape)
        yshape = int(spectrogram.shape[0] / xshape)
        #X, Y, X_ini, Y_ini, X_corr, Y_corr = X[:yshape], Y[:yshape], X_ini[:yshape], Y_ini[:yshape], X_corr[:yshape], Y_corr[:yshape]
        print('Acquisition interruption detected (xshape * yshape > spectrogram.shape[0]) ! New shape : %ix%i (from %s)'%(xshape, yshape, shape_ini))
    if spectrogram.shape[0] > F.shape[0]:
        spectrogram = spectrogram[:F.shape[0]]
        I = I[:F.shape[0]]
        shape_ini = (xshape, yshape)
        yshape = int(spectrogram.shape[0] / xshape)
        print('Acquisition interruption detected (spectrogram.shape[0] > F.shape[0]) ! New shape : %ix%i (from %s)'%(xshape, yshape, shape_ini))

    print('All coordinates and frequencies extracted.')

    # Define spectrum zone if last_quarter_for_threshold is in auto mode
    # Now defined for each spectrum individually
    if general_last_quarter_for_threshold == 'auto':
        print('last_quarter_for_threshold is in auto mode, setting it for each spectrum individually')

    # Get frequencies for each point
    print('Extracting peak frequencies...')
    f1 = []
    f2 = []
    nb_peaks = 0
    nb_false_positives = 0
    false_positives = []
    temp = []
    for ii in range(len(spectrogram)):
        spectrum = spectrogram[ii]
        ff = F[ii]

        if general_last_quarter_for_threshold == 'auto':
            last_quarter_for_threshold = spectrum[:len(spectrum)//4].max() > spectrum[-len(spectrum)//4:].max()
            temp.append(int(last_quarter_for_threshold))
        else:
            last_quarter_for_threshold = general_last_quarter_for_threshold

        # Old method with a middle freq set manually
        if key in parameters.dd_middle and callable(parameters.dd_middle[key]):
            f_middle = middle(ii)
            if f_middle is None:
                f_middle = ff[int(len(ff)/2)]
                idx_middle = int(len(ff)/2)
            else:
                if f_middle > ff.max() or f_middle < ff.min():
                    raise ValueError('f_middle not between ff.min() and ff.max() !')
                idx_middle = abs(ff-f_middle).argmin()
            
            if spectrum[:idx_middle].max() > threshold_peak:
                f1.append(ff[spectrum[:idx_middle].argmax()])
            else:
                f1.append(nan)
            
            if spectrum[idx_middle:].max() > threshold_peak:
                f2.append(ff[idx_middle + spectrum[idx_middle:].argmax()])
            else:
                f2.append(nan)

        # New method
        else:
            value_for_threshold = get_value_for_threshold(spectrum, last_quarter_for_threshold)
            has_peak = criteria_for_peak_selection(spectrum, threshold_peak, required_numpoints, value_for_threshold)
            if not has_peak:
                # No peak detected
                f1.append(nan)
                f2.append(nan)
                f_max = ff[spectrum.argmax()] # Needed if we plot the spectrums even in the middle of the NW

            else:
                # Take max for one of the 2 peaks
                f_max = ff[spectrum.argmax()]

                # Exclude a zone of df on both sides of the max to search for the other peak
                idx_stop1 = max(1, abs(ff - (f_max - df)).argmin())
                idx_start2 = min(len(ff), abs(ff - (f_max + df)).argmin())

                if spectrum[:idx_stop1].max() > spectrum[idx_start2:].max() or f_max+df > ff.max():
                    idx_other_peak = spectrum[:idx_stop1].argmax()

                    # Check in all the search area for a peak
                    has_other_peak = criteria_for_peak_selection(spectrum[:idx_stop1], threshold_peak, required_numpoints, value_for_threshold)
                    if has_other_peak:
                        f1.append(ff[idx_other_peak])
                        f2.append(f_max)
                        only_one_peak = False
                    else:
                        only_one_peak = True

                elif spectrum[:idx_stop1].max() < spectrum[idx_start2:].max() or f_max+df > ff.max():
                    idx_other_peak = idx_start2 + spectrum[idx_start2:].argmax()

                    # Check in all the search area for required_numpoints consecutive above threshold
                    has_other_peak = criteria_for_peak_selection(spectrum[idx_start2:], threshold_peak, required_numpoints, value_for_threshold)
                    if has_other_peak:
                        f2.append(ff[idx_other_peak])
                        f1.append(f_max)
                        only_one_peak = False
                    else:
                        only_one_peak = True

                else:
                   only_one_peak = True

                # If only one peak was detected, determine if it's f1 or f2
                if only_one_peak:
                    mask1 = where(logical_not(isnan(f1)))[0]
                    mask2 = where(logical_not(isnan(f2)))[0]

                    if len(mask1) == 0 or len(mask2) == 0:
                        # Security check, arbitrarly put value in f1
                        f1.append(f_max)
                        f2.append(nan)

##                    # Else : put the value in the array that has the closest last (non-nan) value
##                    # Unless the new value is more than df away from that last value in which case it is probabmy a false positive (not peak)
##                    elif abs(f1[mask1[-1]] - f_max) < abs(f2[mask2[-1]] - f_max):
##                        if 
##                            f1.append(f_max)
##                        else:
##                            f1.append(nan)
##                        f2.append(nan)
##                    elif abs(f1[mask1[-1]] - f_max) >= abs(f2[mask2[-1]] - f_max):
##                        if abs(f_max - f2[mask2[-1]]) < df:
##                            f2.append(f_max)
##                        else:
##                            f2.append(nan)
##                        f1.append(nan)

                    else:
                        m1 = mean([f1[mask1[kk]] for kk in range(max(-10, -len(mask1)), 0, 1)])
                        m2 = mean([f2[mask2[kk]] for kk in range(max(-10, -len(mask2)), 0, 1)])

                        if abs(f_max - m1) < df:
                            # If f_max is close to the last confirmed value of f1 it should always go to f1
                            f1.append(f_max)
                            f2.append(nan)
                        elif abs(f_max - m2) < df:
                            # If f_max is close to the last confirmed value of f2 it should always go to f2
                            f1.append(nan)
                            f2.append(f_max)
                        else:
                            # Reject value as false positive
                            f1.append(nan)
                            f2.append(nan)
                            #print(ii, 'false positive')
                            nb_false_positives += 1
                            false_positives.append((ii, f_max))

        if not isnan(f1[-1]) or not isnan(f2[-1]):
            nb_peaks += 1

        if ii % 100 == 0:
            print(nb_peaks, 'spectrums with peak(s) detected over %i spectrums in total'%(ii+1), '(%i false positives)'%nb_false_positives)

        if plot_spectrums:
            if not isnan(f1[-1]) or not isnan(f2[-1]) or bitmap_ravel[ii] > bitmap.min() + .2*(bitmap.max()-bitmap.min()):
                fig, ax = subplots()
                plot(ff, spectrum)
                axvline(f_max, color='k', linestyle='--')
                axvline(f_max + df, color='c', linestyle='--')
                axvline(f_max - df, color='c', linestyle='--')
                axvline(f1[-1], color='r')
                axvline(f2[-1], color='g')
                axhline(value_for_threshold + threshold_peak, color='y')
                xlabel('Frequency')

                axx = inset_axes(ax, width="30%", height ="30%", loc=2)
                axx.xaxis.set_visible(False)
                axx.yaxis.set_visible(False)
                axx.pcolor(X_ini, Y_ini, bitmap, cmap='gray', shading='nearest')
                axx.invert_yaxis()
                yy, xx = ii//xshape, ii%xshape
                axx.plot(X_ini[yy, xx], Y_ini[yy, xx], 'xr')

                if aller_retour is None:
                    fig.savefig(spectrumsdir+'\\spectrum_'+str(ii).zfill(5))
                else:
                    fig.savefig(spectrumsdir+'\\'+aller_retour+'_spectrum_'+str(ii).zfill(5))
                close(fig)

    print(nb_peaks, 'spectrums with peak(s) detected over %i spectrums in total'%(ii+1))
    f1, f2 = array(f1), array(f2)
    F1 = f1.reshape((yshape, xshape))
    F2 = f2.reshape((yshape, xshape))
    print('Peak frequencies extracted')

    # Plot frequency maps
    print('Plotting...')
    fig_raw1, ax_raw1 = subplots()
    pcolor(X_ini, Y_ini, bitmap, cmap='gray', shading='nearest')
    norm = colors.Normalize(vmin=nanmin(F1), vmax=nanmax(F1))
    sc1 = scatter(X_ini[:yshape], Y_ini[:yshape], c=F1, norm=norm)
    colorbar(sc1)
    xlabel('X position (px)')
    ylabel('Y position (px)')
    title('Frequency map, Mode 1')
    suptitle('Rectangle '+rectangle)
    ax_raw1.invert_yaxis()
    tight_layout()

    fig_raw2, ax_raw2 = subplots()
    pcolor(X_ini, Y_ini, bitmap, cmap='gray', shading='nearest')
    norm = colors.Normalize(vmin=nanmin(F2), vmax=nanmax(F2))
    sc2 = scatter(X_ini[:yshape], Y_ini[:yshape], c=F2, norm=norm)
    colorbar(sc1)
    xlabel('X position (px)')
    ylabel('Y position (px)')
    title('Frequency map, Mode 2')
    suptitle('Rectangle '+rectangle)
    ax_raw2.invert_yaxis()
    tight_layout()

    fig_y_vs_freqs, ax_yvf = subplots()
    yy = Y_ini[:yshape].ravel()
    plot(yy, f1/1e3, 'o', label='Mode 1')
    plot(yy, f2/1e3, 'o', label='Mode 2')
    for ii, f in false_positives:
        plot(yy[ii], f/1e3, 'ok')
    if len(false_positives) > 0:
        plot([], [], 'ok', label='False positives ?') # Dummy plot for legend
    xlabel('Y position (px)')
    ylabel('Frequency (kHz)')
    suptitle('Rectangle '+rectangle)
    legend()
    tight_layout()

    fig_spectrogram, ax_spectrogram = subplots()
    pcolormesh(F/1e3, I, spectrogram, shading='nearest')
    plot(f1/1e3, range(len(f1)), '+r')
    plot(f2/1e3, range(len(f2)), 'xk')
    xlabel('Frequency (kHz)')
    ylabel('Spectrum nb')
    title('Spectrogram')
    suptitle('Rectangle '+rectangle)
    ax_spectrogram.invert_yaxis()
    for ii, f in false_positives:
        plot(f/1e3, ii, 'ow', ms=2, markerfacecolor="None")
    tight_layout()

    """
    fig_lines1, ax_lines1 = subplots()
    for nn in range(F1.shape[0]):
        plot(X_ini[nn], F1[nn], 'o', label='y = %f'%Y_ini[nn, 0])
    legend()
    xlabel('X position (px)')
    ylabel('Frequency (Hz)')
    title('Frequency per line, mode 1')
    suptitle('Rectangle '+rectangle)

    fig_lines2, ax_lines2 = subplots()
    for nn in range(F2.shape[0]):
        plot(X_ini[nn], F2[nn], 'o', label='y = %f'%Y_ini[nn, 0])
    legend()
    xlabel('X position (px)')
    ylabel('Frequency (Hz)')
    title('Frequency per line, mode 2')
    suptitle('Rectangle '+rectangle)
    tight_layout()
    """

    if not old_position_correction:
        fig_cur, ax_cur = subplots()
        plot(num_lines_currents, currents*1e12, 'o')
        xlabel('Line nb. (px)')
        ylabel('Current (pA)')
        title('Currents after each correction')
        suptitle('Rectangle '+rectangle)
        tight_layout()

    fig_temp, ax_temp = subplots()
    plot(temp, 'o')
    title('last_quarter_for_threshold')

    #show()
    print('Plots done.')

    if savefigs:
        print('Saving figures...')
        figs = [fig_raw1, fig_raw2, fig_y_vs_freqs, fig_spectrogram]#, fig_lines1, fig_lines2]
        names = ['Map_raw_mode1', 'Map_raw_mode2', 'Y_vs_Freqs', 'Spectrogram']#, 'Lines_mode1', 'Lines_mode2']
        if not old_position_correction:
            figs.append(fig_cur)
            names.append('Currents')
        figsdir = datadir+u'\\%s\\%s'%(dat, batch)+'\\Figures'
        if not os.path.isdir(figsdir):
            os.mkdir(figsdir)

        for ii in range(len(figs)):
            if aller_retour is None:
                figs[ii].savefig(figsdir+'\\'+batch+'_rectangle_'+rectangle+'_'+names[ii])
            else:
                figs[ii].savefig(figsdir+'\\'+aller_retour+'_'+batch+'_rectangle_'+rectangle+'_'+names[ii])
        print('Figures saved')

    if savefiles:
        print('Saving data files...')
        files = (F1, F2, X_ini, Y_ini)
        names = ('F1', 'F2', 'raw_X', 'raw_Y')
        filesdir = datadir+u'\\%s\\%s'%(dat, batch)+'\\Data_files'
        if not os.path.isdir(filesdir):
            os.mkdir(filesdir)

        for ii in range(len(files)):
            if aller_retour is None:
                savetxt(filesdir+'\\'+batch+'_rectangle_'+rectangle+'_'+names[ii]+'.txt',
                        files[ii])
            else:
                savetxt(filesdir+'\\'+aller_retour+'_'+batch+'_rectangle_'+rectangle+'_'+names[ii]+'.txt',
                        files[ii])
        print('Data files saved')

    print('All good ! \(^_^)/ \n')

    return F, I, spectrogram, F1, F2, f1, f2, false_positives, X_ini, Y_ini, X, Y, bitmap, temp
