# The goal of this script is to construct a linear model of the background spectrum of the NANODROP microspectrometer
# when it is measuring the solutes in 500 uL of solvent in the 2mL glass vials. To this end, we will use the
# 2023-03-25_14-33-50__plate0000003__- folder, which contains 27 spectra of 500 uL of solvent in the 2mL glass vials.
# Then, we will perform principal component analysis on the background spectra, and save the two first principal
# components for later use in a linear model of the background spectra.
import importlib

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from sklearn.decomposition import PCA as pca
import os
# from process_wellplate_spectra import SpectraProcessor, create_folder_unless_it_exists
process_wellplate_spectra = importlib.import_module("uv-vis-absorption-spectroscopy.process_wellplate_spectra")

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
# sp.nanodrop_lower_cutoff_of_wavelengths = 220
# load background
run_name = 'simple-reactions/2023-09-06-run01'
calibration_folder = data_folder + f'{run_name}/' + 'microspectrometer_data/calibration/'
# bkg_spectrum = np.load(calibration_folder + 'background/bkg_spectrum.npy')

# Load all background spectra
spectra = sp.load_all_spectra(data_folder + 'nanodrop-spectrophotometer-measurements/reference_for_E1/MeCN_54times/2023-09-09_20-46-58_MeCN_54times.csv')
bkg_spectrum = np.mean([spectrum for spectrum in spectra], axis=0)
np.save(calibration_folder + 'background/bkg_spectrum.npy', bkg_spectrum)
cut_from = 0
do_save = True

process_wellplate_spectra.create_folder_unless_it_exists(data_folder + f'{run_name}/microspectrometer_data/background_model/')

# Principal component analysis of the background spectra
fig = plt.figure(1, figsize=(10, 5))
nc = 2
pca1 = pca(n_components=nc)
nfeat1 = [spectrum[cut_from:, 1] - bkg_spectrum[cut_from:, 1] for spectrum in spectra]
nfeat1 = [x - np.mean(x) for x in nfeat1]
X1 = pca1.fit(nfeat1)
expl_var_1 = X1.explained_variance_ratio_
sv = X1.components_.T

# plt.plot(np.cumsum(pca1.explained_variance_ratio_))
# plt.show()

for spectrum in spectra:
    plt.plot(spectrum[cut_from:, 0], spectrum[cut_from:, 1] - bkg_spectrum[cut_from:, 1], alpha=0.2, color='grey')
colors = ['C0', 'C1', 'C2']
for i in range(nc):
    plt.plot(spectra[0][cut_from:, 0], sv[:, i], label=f'PC {i+1}', color=colors[i], alpha=0.3)

# Smooth the components and plot them
for i in range(nc):
    smoothed_component = savgol_filter(sv[:, i], polyorder=1, window_length=41)
    plt.plot(spectra[0][cut_from:, 0], smoothed_component, linestyle='--', label=f'smoothed PC {i+1}', color=colors[i])
    for_saving = np.zeros_like(bkg_spectrum[:, 1])
    for_saving[cut_from:] = smoothed_component
    if not do_save:
        continue
    np.save(
       data_folder + f'{run_name}/microspectrometer_data/background_model/component_{i}.npy',
        for_saving)
plt.title('PCA of background spectra')
plt.xlabel('Wavelength, nm')
plt.ylabel('Absorbance units,\nspectra are shifted to be around zero')
plt.legend()
plt.tight_layout()

if do_save:
    fig.savefig(data_folder + f'{run_name}/microspectrometer_data/background_model/pca.png', dpi=300)
    fig.savefig(data_folder + f'{run_name}/microspectrometer_data/background_model/pca.eps')

plt.show()

if do_save:
    np.save(data_folder + f'simple-reactions/2023-07-05-run01/microspectrometer_data/background_model/component_{i}.npy',
        for_saving)
