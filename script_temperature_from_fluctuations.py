# -*- coding: cp1252 -*-
"""
Script for plotting the temperature defined as the integral of the lorentzian peak

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from data_RSA_new import *
from functions import *
from data_selector import DataSelector
from fit_spectres import multi_lorentz
import functions_engraving_study as fes
from scipy.optimize import curve_fit

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

# Arguments that still need to be set
remove_aberrant = True # Remove aberrant frequency values
normalize_from_global_image = False # Use global image taken before acquisition for position normalization
old_position_correction = False # Normally never needs to be changed

# Data selector
ds = DataSelector(description="Select a single rectangle")
ds.savefiles_var.set(False)
ds.savefiles_check['state'] = 'disabled'
ds.select_directory()

rectangle = ds.directory.split('/')[-1].split(' ')[1]
batch = ds.directory.split('/')[-2]
dat = ds.directory.split('/')[-3]
datadir = '/'.join(ds.directory.split('/')[:-3])
dirname = datadir+u'\\%s\\%s\\Réctangle %s\\'%(dat, batch, rectangle) if dat < '20211212' else datadir+u'\\%s\\%s\\Rectangle %s\\'%(dat, batch, rectangle)

sens = ds.sens_var.get()
if not sens in ('10', '01'):
    raise ValueError('Sens should be either "10" or "01", %s received.'%sens)
AR = ds.AR_var.get()
if AR == 'None':
    AR = None
if not AR in ('A', 'R', None):
    raise ValueError('AR should be either "A", "R" or None, %s received.'%AR)

savefigs = ds.savefigs_var.get()

print('Starting plot fits script.')
print('Batch :', batch, '\nRectangle :', rectangle)

""" Load data """
filesdir = datadir+u'\\%s\\%s\\Data_files'%(dat, batch)
figsdir = datadir+u'\\%s\\%s\\Figures\\Normalized_coordinates'%(dat, batch)
savedir = filesdir+u'\\Fits_%s_%s'%(batch, rectangle)

prefixe_AR = AR+'_' if AR in ('A', 'R') else ''
sufixe_AR = ' '+AR if AR in ('A', 'R') else ''
if os.path.isfile(savedir+'\\'+prefixe_AR+batch+'_gammas1.txt'):
    f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = load_fit_data(savedir, batch, prefixe_AR)
else:
    save_fit_data(datadir, dat, batch, rectangle, sens, remove_aberrant, AR, normalize_from_global_image)
    f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = load_fit_data(savedir, batch, prefixe_AR)

if AR is None:
    spectrogram = loadtxt(dirname+'\\Spectres.txt')
    bitmap = loadtxt(dirname+u'\\Bitmap.txt')
    X_ini = loadtxt(dirname+u'\\Tensions X.txt')
    Y_ini = loadtxt(dirname+u'\\Tensions Y.txt')
else:
    if dat < '20211212':
        spectrogram = loadtxt(dirname+'\\Spectres %s.txt'%aller_retour)
        bitmap = loadtxt(dirname+'\\Bitmap %s.txt'%aller_retour)
        X_ini = loadtxt(dirname+'\\Tensions X %s.txt'%aller_retour)
        Y_ini = loadtxt(dirname+'\\Tensions Y %s.txt'%aller_retour)
    else:
        if AR == 'A':
            spectrogram = loadtxt(dirname+'\\Aller\\Spectres.txt')
            if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                bitmap = loadtxt(dirname+'\\Aller\\Bitmap.txt')
                X_ini = loadtxt(dirname+'\\Aller\\Tensions X.txt')
                Y_ini = loadtxt(dirname+'\\Aller\\Tensions Y.txt')
            elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                bitmap = loadtxt(dirname+'\\Aller\\Bitmap')
                X_ini = loadtxt(dirname+'\\Aller\\Tensions X')
                Y_ini = loadtxt(dirname+'\\Aller\\Tensions Y')
            else:
                raise ValueError('No bitmap file found')
        elif AR == 'R':
            spectrogram = loadtxt(dirname+'\\Retour\\Spectres.txt')
            if os.path.isfile(dirname+'\\Aller\\Bitmap.txt'):
                bitmap = loadtxt(dirname+'\\Retour\\Bitmap.txt')
                X_ini = loadtxt(dirname+'\\Retour\\Tensions X.txt')
                Y_ini = loadtxt(dirname+'\\Retour\\Tensions Y.txt')
            elif os.path.isfile(dirname+'\\Aller\\Bitmap'):
                bitmap = loadtxt(dirname+'\\Retour\\Bitmap')
                X_ini = loadtxt(dirname+'\\Retour\\Tensions X')
                Y_ini = loadtxt(dirname+'\\Retour\\Tensions Y')
            else:
                raise ValueError('No bitmap file found')
        else:
            raise ValueError('aller_retour should be either "A" or "R"')
bitmap_ravel = bitmap.ravel()
filename = u'Paramétres spéctres.txt' if dat < '20211212' else u'Parametres spectres.txt'
params = load_params(dirname, filename)
F = get_corrected_freqs(dirname, params, old_corrections=old_position_correction, AR=AR, sens=sens)
I = repeat(range(spectrogram.shape[0]), spectrogram.shape[1]).reshape(spectrogram.shape)
xshape, yshape = [int(float(st)) for st in params[u'Résolution (pixels)'].split('x')] if dat < '20211212' else [int(float(st)) for st in params[u'Resolution (pixels)'].split('x')]

centers, edges, widths, first_line_found, last_line_found = detect_NW_edges(bitmap, X_ini)
X_corrected = get_normalized_X(X_ini, centers, widths, first_line_found, last_line_found)

if normalize_from_global_image:
    print('Normalizing from global image')
    Y_corrected = get_normalized_Y(datadir, dat, batch, rectangle, sens)#, color_box='r')
else:
    print('Normalizing from acquisition image')
    length = abs(Y_ini[first_line_found, 0] - Y_ini[last_line_found, 0])
    basis = Y_ini[last_line_found, 0]
    Y_corrected = (basis - Y_ini) / length

print('Files loaded.')

""" Calculate temperature """
t1 = ampls1 * gammas1**2
t2 = ampls2 * gammas2*2

""" Extract values """
xs = X_corrected.ravel()
ys = Y_corrected.ravel()

mask_1l = where(xs1<0)[0]
mask_1r = where(xs1>0)[0]
mask_2l = where(xs2<0)[0]
mask_2r = where(xs2>0)[0]

# Mode 1 left
x1l = xs1[mask_1l]
y1l = ys1[mask_1l]
t1l = t1[mask_1l]
a1l = ampls1[mask_1l]
g1l = gammas1[mask_1l]

yy1l = []
tt1l = []
aa1l = []
gg1l = []
for y in y1l:
    if not y in yy1l:
        idx = (a1l[where(y1l == y)[0]]).argmax()
        t = t1l[where(y1l == y)[0]][idx]
        a = a1l[where(y1l == y)[0]][idx]
        g = g1l[where(y1l == y)[0]][idx]
        yy1l.append(y)
        tt1l.append(t)
        aa1l.append(a)
        gg1l.append(g)
yy1l = array(yy1l)
tt1l = array(tt1l)
aa1l = array(aa1l)
gg1l = array(gg1l)

# Mode 1 right
x1r = xs1[mask_1r]
y1r = ys1[mask_1r]
t1r = t1[mask_1r]
a1r = ampls1[mask_1r]
g1r = gammas1[mask_1r]

yy1r = []
tt1r = []
aa1r = []
gg1r = []
for y in y1r:
    if not y in yy1r:
        idx = (a1r[where(y1r == y)[0]]).argmax()
        t = t1r[where(y1r == y)[0]][idx]
        a = a1r[where(y1r == y)[0]][idx]
        g = g1r[where(y1r == y)[0]][idx]
        yy1r.append(y)
        tt1r.append(t)
        aa1r.append(a)
        gg1r.append(g)
yy1r = array(yy1r)
tt1r = array(tt1r)
aa1r = array(aa1r)
gg1r = array(gg1r)

# Mode 2 left
x2l = xs2[mask_2l]
y2l = ys2[mask_2l]
t2l = t2[mask_2l]
a2l = ampls2[mask_2l]
g2l = gammas2[mask_2l]

yy2l = []
tt2l = []
aa2l = []
gg2l = []
for y in y2l:
    if not y in yy2l:
        idx = (a2l[where(y2l == y)[0]]).argmax()
        t = t2l[where(y2l == y)[0]][idx]
        a = a2l[where(y2l == y)[0]][idx]
        g = g2l[where(y2l == y)[0]][idx]
        yy2l.append(y)
        tt2l.append(t)
        aa2l.append(a)
        gg2l.append(g)
yy2l = array(yy2l)
tt2l = array(tt2l)
aa2l = array(aa2l)
gg2l = array(gg2l)

# Mode 2 right
x2r = xs2[mask_2r]
y2r = ys2[mask_2r]
t2r = t2[mask_2r]
a2r = ampls2[mask_2r]
g2r = gammas2[mask_2r]

yy2r = []
tt2r = []
aa2r = []
gg2r = []
for y in y2r:
    if not y in yy2r:
        idx = (a2r[where(y2r == y)[0]]).argmax()
        t = t2r[where(y2r == y)[0]][idx]
        a = a2r[where(y2r == y)[0]][idx]
        g = g2r[where(y2r == y)[0]][idx]
        yy2r.append(y)
        tt2r.append(t)
        aa2r.append(a)
        gg2r.append(g)
yy2r = array(yy2r)
tt2r = array(tt2r)
aa2r = array(aa2r)
gg2r = array(gg2r)

""" Fit """
def func_to_fit(y, c):
    return c*fes.phi1(y)
vect_to_fit = vectorize(func_to_fit)

# Mode 1 left
guess = (tt1l.max(), )
bounds_inf = (0, )
bounds_sup = (inf, )
fit = curve_fit(vect_to_fit,
                yy1l,
                tt1l,
                p0=guess,
                bounds=(bounds_inf, bounds_sup))
c1l = fit[0]

# Mode 1 right
guess = (tt1r.max(), )
bounds_inf = (0, )
bounds_sup = (inf, )
fit = curve_fit(vect_to_fit,
                yy1r,
                tt1r,
                p0=guess,
                bounds=(bounds_inf, bounds_sup))
c1r = fit[0]

# Mode 2 left
guess = (tt2l.max(), )
bounds_inf = (0, )
bounds_sup = (inf, )
fit = curve_fit(vect_to_fit,
                yy2l,
                tt2l,
                p0=guess,
                bounds=(bounds_inf, bounds_sup))
c2l = fit[0]

# Mode 2 right
guess = (tt2r.max(), )
bounds_inf = (0, )
bounds_sup = (inf, )
fit = curve_fit(vect_to_fit,
                yy2r,
                tt2r,
                p0=guess,
                bounds=(bounds_inf, bounds_sup))
c2r = fit[0]

""" Plot figures """
def mass_eff(y):
    return fes.phi1(1)**2/fes.phi1(y)**2

fig_map1, ax_map1 = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Mode 1')
suptitle('Rectangle '+rectangle)
if len(ys1) > 0:
    norm = colors.Normalize(vmin=nanmin(t1), vmax=nanmax(t1))
    sc = scatter(xs1, ys1, c=t1, norm=norm)
    colorbar(sc, label='Temperature (arb.)')
xlim(max(-5, nanmin(xs[xs != -inf])), min(5, nanmax(xs[xs != inf])))
ylim(nanmin(ys[ys != -inf]), nanmax(ys[ys != inf]))

fig_map2, ax_map2 = subplots()
pcolor(X_corrected, Y_corrected, bitmap, cmap='gray', shading='nearest')
xlabel('Normalized X position')
ylabel('Normalized Y position')
title('Mode 2')
suptitle('Rectangle '+rectangle)
if len(ys2) > 0:
    norm = colors.Normalize(vmin=nanmin(t2), vmax=nanmax(t2))
    sc = scatter(xs2, ys2, c=t2, norm=norm)
    colorbar(sc, label='Temperature (arb.)')
xlim(max(-5, nanmin(xs[xs != -inf])), min(5, nanmax(xs[xs != inf])))
ylim(nanmin(ys[ys != -inf]), nanmax(ys[ys != inf]))

fig1, ax1 = subplots()
plot(ys1[mask_1l], t1[mask_1l], 'o', label='Mode 1 left')
plot(ys1[mask_1r], t1[mask_1r], 'o', label='Mode 1 right')
xlabel('Normalized y position')
ylabel('Temperature ($A\Gamma^2$) (arb.)')
legend()
title('Temperature mode 1, all points')
tight_layout()

fig2, ax2 = subplots()
plot(ys2[mask_2l], t2[mask_2l], 'o', label='Mode 2 left')
plot(ys2[mask_2r], t2[mask_2r], 'o', label='Mode 2 right')
xlabel('Normalized y position')
ylabel('Temperature ($A\Gamma^2$) (arb.)')
legend()
title('Temperature mode 2, all points')
tight_layout()

fig_fit1, ax_fit1 = subplots()
plot(yy1l, tt1l, 'o', label='Mode 1 left')
plot(yy1l, func_to_fit(yy1l, c1l), label='Fit 1 left')
plot(yy1r, tt1r, 'o', label='Mode 1 right')
plot(yy1r, func_to_fit(yy1r, c1r), label='Fit 1 right')
xlabel('Normalized y position')
ylabel('Temperature ($A\Gamma^2$) (arb.)')
legend()
title('Temperature mode 1, points of max amplitude')
tight_layout()

fig_fit2, ax_fit2 = subplots()
plot(yy2l, tt2l, 'o', label='Mode 2 left')
plot(yy2l, func_to_fit(yy2l, c2l), label='Fit 2 left')
plot(yy2r, tt2r, 'o', label='Mode 2 right')
plot(yy2r, func_to_fit(yy2r, c2r), label='Fit 2 right')
xlabel('Normalized y position')
ylabel('Temperature ($A\Gamma^2$) (arb.)')
legend()
title('Temperature mode 2, points of max amplitude')
tight_layout()

fig_g1, ax_g1 = subplots()
plot(yy1l, gg1l/2/pi, 'o', label='Mode 1 left')
plot(yy1r, gg1r/2/pi, 'o', label='Mode 1 right')
xlabel('Normalized y position')
ylabel('$\Gamma/2\pi$')
legend()
title('Damping rate mode 1, points of max amplitude')
tight_layout()

fig_g2, ax_g2 = subplots()
plot(yy2l, gg2l/2/pi, 'o', label='Mode 2 left')
plot(yy2r, gg2r/2/pi, 'o', label='Mode 2 right')
xlabel('Normalized y position')
ylabel('$\Gamma/2\pi$')
legend()
title('Damping rate mode 2, points of max amplitude')
tight_layout()

fig_a1, ax_a1 = subplots()
semilogy(yy1l, aa1l/2/pi, 'o', label='Mode 1 left')
semilogy(yy1r, aa1r/2/pi, 'o', label='Mode 1 right')
xlabel('Normalized y position')
ylabel('Amplitude (arb.)')
legend()
title('Amplitude mode 1, points of max amplitude')
tight_layout()

fig_a2, ax_a2 = subplots()
semilogy(yy2l, aa2l/2/pi, 'o', label='Mode 2 left')
semilogy(yy2r, aa2r/2/pi, 'o', label='Mode 2 right')
xlabel('Normalized y position')
ylabel('Amplitude (arb.)')
legend()
title('Amplitude mode 2, points of max amplitude')
tight_layout()

show()
