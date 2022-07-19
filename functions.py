# -*- coding: cp1252 -*-
"""
Funtions used to treat data in the nanowire vibrations project

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
import os
import importlib
from numpy import *
from matplotlib.pyplot import imread
from data_RSA_new import load_params

#workdir = "C:\\Users\\cchardin\\Documents\\Grenoble new"
#workdir = "D:\\Documents\\Boulot\\Grenoble"

def load_fit_data(savedir, batch, prefixe_AR=''):
    print('Loading data from :')
    print('Savedir =', savedir)
    print('Batch =', batch)
    if not os.path.isdir(savedir):
        raise OSError

    if os.path.isfile(savedir+'\\'+prefixe_AR+batch+'_f1.txt'):
        f1s = loadtxt(savedir+'\\'+prefixe_AR+batch+'_f1.txt')
        gammas1 = loadtxt(savedir+'\\'+prefixe_AR+batch+'_gammas1.txt')
        ampls1 = loadtxt(savedir+'\\'+prefixe_AR+batch+'_ampls1.txt')
        xs1 = loadtxt(savedir+'\\'+prefixe_AR+batch+'_x1_normalized.txt')
        ys1 = loadtxt(savedir+'\\'+prefixe_AR+batch+'_y1_normalized.txt')
    else:
        f1s, gammas1, ampls1, xs1, ys1 = array([]), array([]), array([]), array([]), array([])

    if os.path.isfile(savedir+'\\'+prefixe_AR+batch+'_f2.txt'):
        f2s = loadtxt(savedir+'\\'+prefixe_AR+batch+'_f2.txt')
        gammas2 = loadtxt(savedir+'\\'+prefixe_AR+batch+'_gammas2.txt')
        ampls2 = loadtxt(savedir+'\\'+prefixe_AR+batch+'_ampls2.txt')
        xs2 = loadtxt(savedir+'\\'+prefixe_AR+batch+'_x2_normalized.txt')
        ys2 = loadtxt(savedir+'\\'+prefixe_AR+batch+'_y2_normalized.txt')
    else:
        f2s, gammas2, ampls2, xs2, ys2 = array([]), array([]), array([]), array([]), array([])
    return (f1s, gammas1, ampls1, xs1, ys1,
            f2s, gammas2, ampls2, xs2, ys2)

def correct_inversions(f1s_all, f2s_all, gammas1_all, gammas2_all, ampls1_all, ampls2_all):
    """ Correct inversions in mode 1 and mode 2 """
    nb_corrections = 0
    for ii in range(len(f1s_all)):
        if not isnan(f1s_all[ii]) and not isnan(f2s_all[ii]) and f1s_all[ii] > f2s_all[ii]:
            tempf = f1s_all[ii]
            f1s_all[ii] = f2s_all[ii]
            f2s_all[ii] = tempf

            tempg = gammas1_all[ii]
            gammas1_all[ii] = gammas2_all[ii]
            gammas2_all[ii] = tempg

            tempa = ampls1_all[ii]
            ampls1_all[ii] = ampls2_all[ii]
            ampls2_all[ii] = tempa

            nb_corrections += 1

    ii = 0
    while abs(f2s_all[:50].mean() - f1s_all[0]) < abs(f1s_all[:50].mean() - f1s_all[0]):
        f2s_all.insert(0, f1s_all[0])
        f1s_all = f1s_all[1:]

        gammas2_all.insert(0, gammas1_all[0])
        gammas1_all = gammas1_all[1:]

        ampls2_all.insert(0, ampls1_all[0])
        ampls1_all = ampls1_all[1:]

        ii += 1

    jj= 0
    while abs(f2s_all[:50].mean() - f2s_all[0]) > abs(f1s_all[:50].mean() - f2s_all[0]):
        f1s_all.insert(0, f2s_all[0])
        f2s_all = f2s_all[1:]

        gammas1_all.insert(0, gammas2_all[0])
        gammas2_all = gammas2_all[1:]

        ampls1_all.insert(0, ampls2_all[0])
        ampls2_all = ampls2_all[1:]

        jj += 1

    print('Inversions corrected. (%i corrections made, then %i and %i)'%(nb_corrections, ii, jj))

def detect_NW_edges(bitmap, X_ini):
    """ Detection of the edges of the NW and left-right separation"""
    threshold = bitmap.min() + .5*(bitmap.max() - bitmap.min())
    edges = []
    centers = []
    widths = []
    masks = []
    for ii, xx in enumerate(X_ini):
        mask = where(bitmap[ii] > threshold)[0]
        masks.append(mask)
        if len(mask) > 0:
            left, right = xx[mask.min()], xx[mask.max()]
            center = (left+right)/2
            width = abs(right-left)
            edges.append((left, right))
            centers.append(center)
            widths.append(width)
        else:
            edges.append((nan, nan))
            centers.append(nan)
            widths.append(nan)
    centers = array(centers)
    widths = array(widths)
    edges = array(edges)
    masks = array(masks)

    # Detect 1st and last lines
    first_line_found = None
    last_line_found = None
    for ii, center in enumerate(centers):
        if isnan(center):
            if not first_line_found is None and last_line_found is None:
                last_line_found = ii
        else:
            if first_line_found is None:
                first_line_found = ii
    print('First line detected :', first_line_found)
    print('Last line detected :', last_line_found)
    if last_line_found is None:
        # We did not see the end of the NW, so we conside the last line measured is the end
        last_line_found = len(centers)-1
        print('Last line detected was None, now set to', last_line_found)

    return centers, edges, widths, first_line_found, last_line_found

def get_normalized_X(X_ini, centers, widths, first_line_found, last_line_found):
    X_corrected = zeros(X_ini.shape)
    for ii, xx in enumerate(X_ini):
        if ii < first_line_found:
            # Apply a dummy correction for lines where there is no NW and get NW length
            center_1st = centers[first_line_found]
            width_1st = widths[first_line_found]
            X_corrected[ii] = (X_ini[ii] - center_1st)/width_1st
        elif ii >= last_line_found:
            # Apply a dummy correction for lines where there is no NW and get NW length
            center_last = centers[last_line_found-1]
            width_last = widths[last_line_found-1]
            X_corrected[ii] = (X_ini[ii] - center_last)/width_last
        else:
            X_corrected[ii] = (X_ini[ii] - centers[ii])/widths[ii]
    return X_corrected

def popts_to_fit_data(datadir,
                      dat,
                      batch,
                      rectangle,
                      sens,
                      remove_aberrant=True,
                      AR=None,
                      normalize_from_global_image=False,
                      only_two_peaks=False,
                      select_from_gradient=None, # Either the minimal value for gradient or None
                      ):

    dirname = datadir+u'\\%s\\%s\\Réctangle %s\\'%(dat, batch, rectangle) if dat < '20211213' else datadir+u'\\%s\\%s\\Rectangle %s\\'%(dat, batch, rectangle)
    filesdir = datadir+u'\\%s\\%s\\Data_files'%(dat, batch)
    savedir = filesdir+u'\\\\Fits_%s_%s'%(batch, rectangle)
    if normalize_from_global_image:
        savedir += '_Norm_from_global_image'
    if not os.path.isdir(savedir):
        os.mkdir(savedir)

    """ Load array data """
    prefixe_AR = AR+'_' if AR in ('A', 'R') else ''
    sufixe_AR = ' '+AR if AR in ('A', 'R') else ''
    f1s_all = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_f1s.txt'%(batch, rectangle))
    gammas1_all = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_gammas1.txt'%(batch, rectangle))
    ampls1_all = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_ampls1.txt'%(batch, rectangle))
    f2s_all = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_f2s.txt'%(batch, rectangle))
    gammas2_all = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_gammas2.txt'%(batch, rectangle))
    ampls2_all = loadtxt(filesdir+'\\'+prefixe_AR+'%s_rectangle_%s_Fit_ampls2.txt'%(batch, rectangle))

    X_ini = loadtxt(filesdir+'\\'+prefixe_AR+batch+'_rectangle_'+rectangle+'_raw_X.txt')
    Y_ini = loadtxt(filesdir+'\\'+prefixe_AR+batch+'_rectangle_'+rectangle+'_raw_Y.txt')
    if AR is None or dat < '20211212':
        bitmap = loadtxt(dirname+u'\\Bitmap'+sufixe_AR+'.txt')
    elif AR == 'A':
        if dat < '20211212':
            bitmap = loadtxt(dirname+u'\\Bitmap A.txt')
        else:
            if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                bitmap = loadtxt(dirname+'\\Aller\\Bitmap.txt')
            elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                bitmap = loadtxt(dirname+'\\Aller\\Bitmap')
            else:
                raise ValueError('No bitmap file found')
    elif AR == 'R':
        if dat < '20211212':
            bitmap = loadtxt(dirname+u'\\Bitmap R.txt')
        else:
            if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                bitmap = loadtxt(dirname+'\\Retour\\Bitmap.txt')
            elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                bitmap = loadtxt(dirname+'\\Retour\\Bitmap')
            else:
                raise ValueError('No bitmap file found')
        

    """ Correct inversions in mode 1 and mode 2 """
    correct_inversions(f1s_all, f2s_all, gammas1_all, gammas2_all, ampls1_all, ampls2_all)

    """ Detection of the edges of the NW and left-right separation"""
    centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(bitmap, X_ini)

    # Normalize X by width
    X_corrected = get_normalized_X(X_ini, centers, widths, first_line_found, last_line_found)
    
    # Normalize Y by length
    if normalize_from_global_image:
        print('Normalizing from global image')
        Y_corrected = get_normalized_Y(datadir, dat, batch, rectangle, sens)#, color_box='r')
    else:
        print('Normalizing from acquisition image')
        length = abs(Y_ini[first_line_found, 0] - Y_ini[last_line_found, 0])
        basis = Y_ini[last_line_found, 0]
        Y_corrected = (basis - Y_ini) / length

    ys = Y_corrected.flatten()
    xs = X_corrected.flatten()

    # Remove nan values and select from gradient if required
    if not select_from_gradient is None:
        f1s, f2s, ys1, ys2, xs1, xs2 = [], [], [], [], [], []
        gammas1, gammas2, ampls1, ampls2 = [], [], [], []

        f1s_all = f1s_all.reshape(bitmap.shape)
        f2s_all = f2s_all.reshape(bitmap.shape)
        gammas1_all = gammas1_all.reshape(bitmap.shape)
        gammas2_all = gammas2_all.reshape(bitmap.shape)
        ampls1_all = ampls1_all.reshape(bitmap.shape)
        ampls2_all = ampls2_all.reshape(bitmap.shape)

        grad = gradient(bitmap, axis=1)
        rows, cols = where(abs(grad) > select_from_gradient)

        for nn in range(len(rows)):
            ii, jj = rows[nn], cols[nn]
            if not isnan(f1s_all[ii, jj]):
                f1s.append(f1s_all[ii, jj])
                xs1.append(X_corrected[ii, jj])
                ys1.append(Y_corrected[ii, jj])
                gammas1.append(gammas1_all[ii, jj])
                ampls1.append(ampls1_all[ii, jj])

            if not isnan(f2s_all[ii, jj]):
                f2s.append(f2s_all[ii, jj])
                xs2.append(X_corrected[ii, jj])
                ys2.append(Y_corrected[ii, jj])
                gammas2.append(gammas2_all[ii, jj])
                ampls2.append(ampls2_all[ii, jj])

        f1s, f2s = array(f1s), array(f2s)
        xs1, xs2 = array(xs1), array(xs2)
        ys1, ys2 = array(ys1), array(ys2)
        gammas1, gammas2 = array(gammas1), array(gammas2)
        ampls1, ampls2 = array(ampls1), array(ampls2)

    else:
        mask_nan_1 = where(logical_not(isnan(f1s_all)))[0]
        mask_nan_2 = where(logical_not(isnan(f2s_all)))[0]
        if only_two_peaks:
            mask = intersect1d(mask_nan_1, mask_nan_2, assume_unique=True)
            mask_nan_1 = mask
            mask_nan_2 = mask
            print('Selecting only two peaks')
            print('mask two peaks :', mask)

        f1s = f1s_all[mask_nan_1]
        gammas1 = gammas1_all[mask_nan_1]
        ampls1 = ampls1_all[mask_nan_1]
        ys1 = ys[mask_nan_1]
        xs1 = xs[mask_nan_1]

        f2s = f2s_all[mask_nan_2]
        gammas2 = gammas2_all[mask_nan_2]
        ampls2 = ampls2_all[mask_nan_2]
        ys2 = ys[mask_nan_2]
        xs2 = xs[mask_nan_2]

    """ Remove aberrant values """
    if remove_aberrant:
        print('Removing aberrant values')
        try:
            key = 'rectangle_%s'%rectangle if AR is None else 'rectangle_%s_%s'%(rectangle, AR)
            # Loading parameters file
            spec = importlib.util.spec_from_file_location('parameters', datadir+u'\\%s\\%s'%(dat, batch)+'\\parameters.py')
            parameters = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(parameters)

            min_freq1 = parameters.dd_min_freq1[key] if key in parameters.dd_min_freq1 else f1s.min()-1
            max_freq2 = parameters.dd_max_freq2[key] if key in parameters.dd_max_freq2 else f2s.max()+1
            try:
                max_freq1 = parameters.dd_max_freq1[key] if key in parameters.dd_max_freq1 else f1s.max()+1
                min_freq2 = parameters.dd_min_freq2[key] if key in parameters.dd_min_freq2 else f2s.min()-1
            except AttributeError:
                middle_freq = parameters.dd_middle[key] if key in parameters.dd_middle else (f1s.max()+f2s.min())/2
                max_freq1 = middle_freq
                min_freq2 = middle_freq
            if not min_freq1 <= max_freq1 and max_freq1 <= min_freq2 and min_freq2 <= max_freq2:
                raise ValueError('Wrong values if min/middle/max freqs, check parameters file.')

            to_remove1 = where((f1s < min_freq1) | (f1s > max_freq1))[0]
            to_remove2 = where((f2s < min_freq2) | (f2s > max_freq2))[0]
            if only_two_peaks:
                to_remove = union1d(to_remove1, to_remove2)
                to_remove1 = to_remove
                to_remove2 = to_remove

            print(len(f1s), 'points in f1s initially')
            f1s = delete(f1s, to_remove1)
            gammas1 = delete(gammas1, to_remove1)
            ampls1 = delete(ampls1, to_remove1)
            ys1 = delete(ys1, to_remove1)
            xs1 = delete(xs1, to_remove1)
            print(len(to_remove1), 'points removed from f1s.')

            print(len(f2s), 'points in f2s initially')
            f2s = delete(f2s, to_remove2)
            gammas2 = delete(gammas2, to_remove2)
            ampls2 = delete(ampls2, to_remove2)
            ys2 = delete(ys2, to_remove2)
            xs2 = delete(xs2, to_remove2)
            print(len(to_remove2), 'points removed from f2s.')
        except FileNotFoundError: # Use this if you have python 3.x
            print('No parameter file, did not remove anything.')
    return (f1s, gammas1, ampls1, xs1, ys1,
            f2s, gammas2, ampls2, xs2, ys2)


def save_fit_data(datadir,
                  dat,
                  batch,
                  rectangle,
                  sens,
                  remove_aberrant=True,
                  AR=None,
                  normalize_from_global_image=False,
                  select_from_gradient=None, # Either the minimal value for gradient or None
                  only_two_peaks=False,
                  ):
    f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = popts_to_fit_data(datadir,
                                                                                       dat,
                                                                                       batch,
                                                                                       rectangle,
                                                                                       sens,
                                                                                       remove_aberrant=remove_aberrant,
                                                                                       AR=AR,
                                                                                       normalize_from_global_image=normalize_from_global_image,
                                                                                       only_two_peaks=only_two_peaks,
                                                                                       select_from_gradient=select_from_gradient)

    dirname = datadir+u'\\%s\\%s\\Réctangle %s\\'%(dat, batch, rectangle) if dat < '20211213' else datadir+u'\\%s\\%s\\Rectangle %s\\'%(dat, batch, rectangle)
    filesdir = datadir+u'\\%s\\%s\\Data_files'%(dat, batch)
    savedir = filesdir+u'\\Fits_%s_%s'%(batch, rectangle)
    if not select_from_gradient is None:
        savedir += '_Select_from_gradient'
    if normalize_from_global_image:
        savedir += '_Norm_from_global_image'
    if only_two_peaks:
        savedir += '_Only_two_peaks'

    prefixe_AR = AR+'_' if AR in ('A', 'R') else ''

    if not os.path.isdir(savedir):
        os.mkdir(savedir)

    """ Save files """
    savetxt(savedir+'\\'+prefixe_AR+batch+'_f1.txt', f1s)
    savetxt(savedir+'\\'+prefixe_AR+batch+'_gammas1.txt', gammas1)
    savetxt(savedir+'\\'+prefixe_AR+batch+'_ampls1.txt', ampls1)
    savetxt(savedir+'\\'+prefixe_AR+batch+'_x1_normalized.txt', xs1)
    savetxt(savedir+'\\'+prefixe_AR+batch+'_y1_normalized.txt', ys1)

    savetxt(savedir+'\\'+prefixe_AR+batch+'_f2.txt', f2s)
    savetxt(savedir+'\\'+prefixe_AR+batch+'_gammas2.txt', gammas2)
    savetxt(savedir+'\\'+prefixe_AR+batch+'_ampls2.txt', ampls2)
    savetxt(savedir+'\\'+prefixe_AR+batch+'_x2_normalized.txt', xs2)
    savetxt(savedir+'\\'+prefixe_AR+batch+'_y2_normalized.txt', ys2)

    print('Files saved')

    if len(ys1) > 0:
        print('Min position for mode 1 : y =', round(ys1.min()*100), '%')
    else:
        print('Min position for mode 1 : No points !')
    if len(ys2) > 0:
        print('Min position for mode 2 : y =', round(ys2.min()*100), '%')
    else:
        print('Min position for mode 2 : No points !')

def detect_box(image, color='r'):
    """ Get the rectangle on which the measurement was made """
    if color == 'r':
        idx_color = 0
    elif color == 'g':
        idx_color = 1
    elif color == 'b':
        idx_color = 2

    ar = image[:, :, idx_color]
    rows, cols = where(ar==1)

    return (rows.min(), rows.max()), (cols.min(), cols.max())

def rgb2gray(rgb):
    return dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])

def get_closeup(datadir, dat, batch, rectangle, color_box='r'):
    dirname = datadir+u'\\%s\\%s\\Réctangle %s\\'%(dat, batch, rectangle) if dat < '20211213' else datadir+u'\\%s\\%s\\Rectangle %s\\'%(dat, batch, rectangle)
    image_path = dirname+'\\Image global.png'

    """ Load 'Image global.png' for the given measurement """
    im = imread(image_path)

    (i1, i2), (j1, j2) = detect_box(im, color=color_box)
    closeup = rgb2gray(im[i1:i2, j1+1:j2, :])

    return closeup, i1, i2, j1, j2

def get_normalized_Y(datadir, dat, batch, rectangle, sens, color_box='r'):
    """ Gives the normalized Y corrdinates of the measurement, using the global image for normalization """
    closeup, i1, i2, j1, j2 = get_closeup(datadir, dat, batch, rectangle, color_box=color_box)

    """ Detect the NW """
    Y, X = mgrid[i1:i2, j1+1:j2]
    centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(closeup, X)

    """ Get the actual pixels positions with the right measurement resolution """
    dirname = datadir+u'\\%s\\%s\\Réctangle %s\\'%(dat, batch, rectangle) if dat < '20211213' else datadir+u'\\%s\\%s\\Rectangle %s\\'%(dat, batch, rectangle)
    params = load_params(dirname)
    xshape, yshape = [int(float(st)) for st in params[u'Résolution (pixels)'].split('x')] if u'Résolution (pixels)' in params \
                     else [int(float(st)) for st in params[u'Resolution (pixels)'].split('x')]

    x_positions = linspace(j1, j2, xshape)
    y_positions = linspace(i1, i2, yshape)
    X_meas, Y_meas = meshgrid(x_positions, y_positions)

    if sens == '01':
        basis = Y[first_line_found, 0]
        end = Y[last_line_found, 0]
    elif sens == '10':
        basis = Y[last_line_found, 0]
        end = Y[first_line_found, 0]
    length_px = end - basis
    Y_meas_norm = (Y_meas - basis)/length_px

    return Y_meas_norm

def get_ncols_nrows(params):
    st = params['Résolution (pixels)'] if u'Résolution (pixels)' in params \
         else [int(float(st)) for st in params[u'Resolution (pixels)'].split('x')]
    tu = st.split('x')
    ncols = int(float(tu[0]))
    nrows = int(float(tu[1]))
    return ncols, nrows

def get_dimensions(params):
    st = params['Dimensions (m)']
    tu = st.split('x')
    dimx = float(tu[0])
    dimy = float(tu[1])
    return dimx, dimy

def get_real_length_diameter_from_rectangle(bitmap, X_ini, Y_ini, params):
    centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(bitmap, X_ini)

    dimensions_image = get_dimensions(params)
    ncols, nrows = get_ncols_nrows(params)

    pixel_size_x = dimensions_image[0]/ncols
    pixel_size_y = dimensions_image[1]/nrows

    # Extract diameter value for each line
    diameters = []
    for ii, (left, right) in enumerate(edges):
        xx = X_ini[ii]

        idx_left = abs(xx-left).argmin()
        idx_right = abs(xx-right).argmin()

        diameter = pixel_size_x * (idx_right - idx_left)
        diameters.append(diameter)

    # Extract length value
    length = pixel_size_y * abs(last_line_found - first_line_found)

    return diameters, length

def get_real_length_diameter_from_closeup(closeup, params):
    X, Y = meshgrid(range(closeup.shape[1]), range(closeup.shape[0]))
    centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(closeup, X)

    dimensions_image = get_dimensions(params)
    ncols, nrows = get_ncols_nrows(params)

    pixel_size_x = dimensions_image[0]/ncols
    pixel_size_y = dimensions_image[1]/nrows

    # Extract diameter value for each line
    diameters = []
    for ii, (left, right) in enumerate(edges):
        diameter = pixel_size_x * abs(right - left)
        diameters.append(diameter)

    # Extract length value
    length = pixel_size_y * abs(last_line_found - first_line_found)

    return diameters, length
