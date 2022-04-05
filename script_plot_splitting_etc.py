"""
Script for plotting mode splitting, difference of gamma and amplitude ratio.

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
from numpy import *
from matplotlib.pyplot import *
import os
from functions import *
from data_selector import DataSelector
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

# To disable matplotlib's UserWarning
import warnings
warnings.filterwarnings("ignore")

# Arguments that still need to be set
remove_aberrant = True # Remove aberrant frequency values
normalize_from_global_image = False # Use global image taken before acquisition for position normalization

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
savefiles = ds.savefiles_var.get()

print('Starting plot data and maps script.')
print('Batch :', batch, '\nRectangle :', rectangle)

""" Get the data """
f1s, gammas1, ampls1, xs1, ys1, f2s, gammas2, ampls2, xs2, ys2 = popts_to_fit_data(datadir,
                                                                                   dat,
                                                                                   batch,
                                                                                   rectangle,
                                                                                   sens,
                                                                                   remove_aberrant=True,
                                                                                   AR=AR,
                                                                                   normalize_from_global_image=False,
                                                                                   only_two_peaks=True,
                                                                                   select_from_gradient=None,
                                                                                   )
df = f2s - f1s
dg = gammas2 - gammas1
ra = ampls2/ampls1

maskl = where(xs1 < 0)
maskr = where(xs1 > 0)

""" Plot figures """
fig_split, ax_split = subplots()
plot(ys1[maskl], df[maskl], 'o', label='Left')
plot(ys1[maskr], df[maskr], 'o', label='Right')
xlabel('y/L')
ylabel('Mode splitting (Hz)')
suptitle('Mode splitting')
title('%s - Rectangle %s'%(batch, rectangle))
legend()

fig_gamma, ax_gamma = subplots()
plot(ys1[maskl], dg[maskl], 'o', label='Left')
plot(ys1[maskr], dg[maskr], 'o', label='Right')
axhline(0, color='k', linestyle='--')
xlabel('y/L')
ylabel('$\Gamma_2/2\pi - \Gamma_1/2\pi$ (Hz)')
suptitle('Difference of damping')
title('%s - Rectangle %s'%(batch, rectangle))
legend()

fig_ampl, ax_ampl = subplots()
semilogy(ys1[maskl], ra[maskl], 'o', label='Left')
semilogy(ys1[maskr], ra[maskr], 'o', label='Right')
xlabel('y/L')
ylabel('$A_2 / A_1$')
suptitle('Amplitude ratio')
title('%s - Rectangle %s'%(batch, rectangle))
legend()

fig_pts, ax_pts = subplots()
plot(xs1[maskl], ys1[maskl], 'o', label='Left')
plot(xs1[maskr], ys1[maskr], 'o', label='Right')
rects = [Rectangle((-.5, 0), 1, 1)]
pc = PatchCollection(rects, facecolor='gray', alpha=.5, edgecolor='k')
ax_pts.add_collection(pc)
legend()
xlabel('x/D')
ylabel('y/D')
suptitle('Position of points')
title('%s - Rectangle %s'%(batch, rectangle))

show()

figsdir = datadir+u'\\%s\\%s\\Figures\\Splitting_etc'%(dat, batch)
if not os.path.isdir(figsdir):
    os.mkdir(figsdir)
if savefigs:
    fig_split.savefig(figsdir+'\\splitting')
    fig_gamma.savefig(figsdir+'\\gamma_difference')
    fig_ampl.savefig(figsdir+'\\amplitude_ratio')
    fig_pts.savefig(figsdir+'\\points_position')
