"""
Détection des pics: si la méthode par défaut est utilisée (above_noise_mean) :
    - On sélectionne une zone où on sait qu'il n'y a pas de pics
        par défaut, c'est fait en comparant le 1er quart de chaque spectre avec le dernier quart
        et en gardant celui dont le maximum est le plus bas (ça suppose qu'un de ces deux quarts n'a pas de pics)

    - On prend le max de ce bout de spectre

    - On recherche dans le reste du spectre N points consécutifs au-dessus de cette valeur +  threshold_peak
        - N est définie par required_numpoints
            - defaut=10
            - peut être modifié :
                soit dans script_map_frequence_batch en ajoutant 'required_numpoints'=nouvelle_valeur dans kwargs
                soit dans le fichier parameters.py de chaque acquisition
        - threshold_peak peut être modifié avec les deux mêmes méthodes (défaut = 3 dB)

Ce scripte aide à vérifier si les pics devraient être détectés, mais CE N'EST QU'UN APERCU !
Ici on ne regarde que le max des spectres (ie required_numpoints = 1)


Sinon on peut aussi tracer le spectrogram et voir si il y a des pics manifestement ratés
"""
from numpy import *
from matplotlib.pyplot import *

# Load spectrogram : dirname is the directory
dirname = 'D:\\Documents\\Boulot\\Grenoble\\Data\\20220622\\Fil 6_3\\Rectangle 1'
spectrogram = loadtxt(dirname+'\\Spectres.txt')
print('Spectrogram loaded')

# Valeur de threshold_peak qu'on veut tester (défaut = 3)
threshold_peak = 3

threshold_first_quarter = spectrogram[:, :spectrogram.shape[1]//4].max(axis=1) + threshold_peak
peaks_first = where(spectrogram.max(axis=1) > threshold_first_quarter)[0]
print('Nb peaks detected with 1st quarter :', len(peaks_first))

threshold_last_quarter = spectrogram[:, 3*spectrogram.shape[1]//4:].max(axis=1) + threshold_peak
peaks_last = where(spectrogram.max(axis=1) > threshold_last_quarter)[0]
print('Nb peaks detected with last quarter :', len(peaks_last))

idx = array(range(len(spectrogram)))

fig, (ax1, ax2) = subplots(ncols=2)
ax1.plot(spectrogram.max(axis=1), label='Max of spectra')
ax1.plot(threshold_first_quarter, label='Threshold 1st quarter')
ax1.plot(idx[peaks_first], spectrogram.max(axis=1)[peaks_first], 'o', label='Peaks detected')
ax1.legend()

ax2.plot(spectrogram.max(axis=1), label='Max of spectra')
ax2.plot(threshold_last_quarter, label='Threshold last quarter')
ax2.plot(idx[peaks_last], spectrogram.max(axis=1)[peaks_last], 'o', label='Peaks detected')
ax2.legend()

fig_sp, ax_sp = subplots()
imshow(spectrogram)
title('Spectrogram')

show()
