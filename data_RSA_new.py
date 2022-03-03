# -*- coding: cp1252 -*-
"""
Functions for extracted useful data from the raw file created by the acquisition program ViBr

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
import io

def load_params(dirname, filename=u'Parametres spectres.txt'):
    if dirname[-1:] in ('\\', '/'):
        complete_path = dirname + filename
    else:
        complete_path = dirname + '\\' + filename
    with open(complete_path, 'r') as ff:
        lines = ff.readlines()
##    txt = io.open(dirname+'\\'+filename, 'rb')#, encoding='utf-8')
##    lines = txt.readlines()

    dd = {}
    for line in lines:
        if line[-1] == '\n':
            line = line[:-1]
        if len(line) > 0:
            ll = line.split('\t')
            dd[ll[0]] = ll[1]

    return dd

def get_corrected_freqs(dirname,
                        params,
                        filename_corrections='Corrections.txt',
                        old_corrections=False,
                        AR=None,
                        sens='10'):
    assert AR in ('A', 'R', None)
    assert sens in ('10', '01')

    skiprows = 0 if old_corrections else 1
    corrections = loadtxt(dirname+'\\'+filename_corrections, skiprows=skiprows)

    if 'SPAN (Hz)' in params:
        span = float(params['SPAN (Hz)'])
    else:
        span = .2 * corrections[4, 0] if old_corrections else .2 * corrections[0, 4]

    nb_points = int(params['Nombre de point par spectre']) if 'Nombre de point par spectre' in params else 4001

    xshape, yshape = [int(float(st)) for st in params['Résolution (pixels)'].split('x')] if 'Résolution (pixels)' in params \
                     else [int(float(st)) for st in params['Resolution (pixels)'].split('x')]

    F = []
    my_range = corrections.shape[1] if old_corrections else  \
               corrections.shape[0] if len(corrections.shape) > 1 else 1
    for ii in range(my_range):
        if old_corrections:
            corr_x, corr_y, nb_px_x, nb_px_y, center_freq = corrections[:, ii]
        else:
            if my_range == 1:
                corr_x, corr_y, nb_px_x, nb_px_y, center_freq, current = corrections
            else:
                corr_x, corr_y, nb_px_x, nb_px_y, center_freq, current = corrections[ii]
        for jj in range(int(nb_px_x * nb_px_y)):
            ff = linspace(center_freq - span/2, center_freq + span/2, nb_points)
            F.append(ff)
    F = array(F)

    if AR == 'A':
        idx = F.shape[0] // 2
        F = F[:idx]
    elif AR == 'R':
        idx = F.shape[0] // 2
        F = F[idx:]

    if sens == '01':
        F = flip(F, axis=0)

    return F

def get_corrected_coordinates(dirname,
                              params,
                              filename_corrections='Corrections.txt',
                              old_corrections=False):
    skiprows = 0 if old_corrections else 1
    corrections = loadtxt(dirname+'\\'+filename_corrections, skiprows=skiprows)
    xshape, yshape = [int(float(st)) for st in params['Résolution (pixels)'].split('x')] if 'Résolution (pixels)' in params \
                     else [int(float(st)) for st in params['Resolution (pixels)'].split('x')]

    X = zeros((yshape, xshape))
    Y = zeros((yshape, xshape))
    idx = 0
    my_range = corrections.shape[1] if old_corrections else \
               corrections.shape[0] if len(corrections.shape) > 1 else 1
    for ii in range(my_range):
        if old_corrections:
            corr_x, corr_y, nb_px_x, nb_px_y, center_freq = corrections[:, ii]
        else:
            if my_range == 1:
                corr_x, corr_y, nb_px_x, nb_px_y, center_freq, current = corrections
            else:
                corr_x, corr_y, nb_px_x, nb_px_y, center_freq, current = corrections[ii]
        X[idx:idx+int(nb_px_x)] += corr_x
        Y[idx:idx+int(nb_px_x)] += corr_y
        idx += int(nb_px_x)

    return X, Y
