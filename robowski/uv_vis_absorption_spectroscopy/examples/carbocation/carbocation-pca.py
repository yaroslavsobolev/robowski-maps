# The goal of this script is to construct a linear model of the background spectrum of the CRAIC microspectrometer
# when it is measuring the solutes in 500 uL of DMF in the 2mL glass vials. To this end, we will use the
# 2023-03-25_14-33-50__plate0000003__- folder, which contains 27 spectra of 500 uL of DMF in the 2mL glass vials.
# Then, we will perform principal component analysis on the background spectra, and save the two first principal
# components for later use in a linear model of the background spectra.

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from sklearn.decomposition import PCA as pca
import os
import importlib
process_wellplate_spectra = importlib.import_module("uv_vis_absorption_spectroscopy.process_wellplate_spectra")

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'
sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')

litdata = np.loadtxt('misc_scripts/literature_spectra/carbocation_Yuichi_Nishimae_2004.txt', delimiter='\t')
litdata_e12 = np.loadtxt('misc_scripts/literature_spectra/ammer-sailer-riedle-2012-E12.txt', delimiter='\t')
litdata_e3 = np.loadtxt('misc_scripts/literature_spectra/ammer-sailer-riedle-2012-E3.txt', delimiter='\t')

# load background
plate_folder = craic_folder + '2023-03-25_14-33-50__plate0000003__-' + '/'
calibration_folder = data_folder + 'multicomp-reactions/2023-01-18-run01/' + 'microspectrometer_data/calibration/'
calibrant_shortnames = ['IIO029A', 'ald001']
calibrant_upper_bounds = [np.inf, 1e-10]

calibrants = []
for calibrant_shortname in calibrant_shortnames:
    dict_here = dict()
    dict_here['coeff_to_concentration_interpolator'], dict_here['reference_interpolator'], dict_here[
        'bkg_spectrum'] = \
        sp.load_calibration_for_one_calibrant(calibrant_shortname, calibration_folder)
    calibrants.append(dict_here.copy())
bkg_spectrum = calibrants[0]['bkg_spectrum']

wavelengths = bkg_spectrum[:, 0]

def wavelength_id_by_value(wavelengths, value):
    return np.argmin(np.abs(wavelengths - value))
# interm_wav_id = wavelength_id_by_value(wavelengths, 510)
# interm_wav_id2 = wavelength_id_by_value(wavelengths, 645)

interm_wav_id = wavelength_id_by_value(wavelengths, 575)
interm_wav_id2 = wavelength_id_by_value(wavelengths, 645)

# print(interm_wav_id)
# print(interm_wav_id2)
# print(wavelength_id_by_value(wavelengths, 416))
# print(wavelength_id_by_value(wavelengths, 438))

# Load all spectra
plates = ['2023-04-08_11-09-22__plate0000006__2023-04-07-run01/',
          '2023-04-08_11-43-19__plate0000009__2023-04-07-run01/',
          '2023-04-08_12-56-26__plate0000011__2023-04-07-run01/'
          ]
# plates = ['2023-04-12_20-35-01__plate0000020__simple_reactions_2023-04-12-run01-diluted/']

spectra = sp.load_all_spectra(craic_folder + plates[0])
for i in range(1, len(plates)):
    spectra.extend(sp.load_all_spectra(craic_folder + plates[i]))

# empty wells mask
empty_vial_mask = np.array([spectrum[-1, 1] < 0.3 for spectrum in spectra])
spectra = [spectrum for spectrum in np.array(spectra)[empty_vial_mask, :, :]]

# carbocat_mask = np.array([spectrum[interm_wav_id, 1] - spectrum[interm_wav_id2, 1] > 0.16 for spectrum in spectra])
# carbocat_mask = np.array([spectrum[interm_wav_id, 1] - spectrum[interm_wav_id2, 1] > 0.16 for spectrum in spectra])
carbocat_mask = np.array([spectrum[interm_wav_id, 1] - spectrum[interm_wav_id2, 1] > 0.037 for spectrum in spectra])
cut_from = 160
spectra_masked = np.array(spectra)[carbocat_mask, :, :]
spectra_masked = [spectrum for spectrum in spectra_masked]
spectra_masked = spectra_masked[:-1]

# plt.axvline(x=wavelengths[interm_wav_id], color='black')

spectra_unmasked = np.array(spectra)[~carbocat_mask, :, :]
spectra_unmasked = [spectrum for spectrum in spectra_unmasked]

f1 = plt.figure(figsize=(4, 4))

# for spectrum in spectra_unmasked:
#     ys = spectrum[cut_from:, 1] - bkg_spectrum[cut_from:, 1]
#     ys = ys - ys[700-350 + cut_from + 60]
#     plt.plot(spectrum[cut_from:, 0], ys, alpha=0.1, color='grey')
# for i, spectrum in enumerate(spectra_masked):
#     ys = spectrum[cut_from:, 1] - bkg_spectrum[cut_from:, 1]
#     ys = ys - ys[700-350 + cut_from + 60]
#     if i == 0:
#         plt.plot(spectrum[cut_from:, 0], ys, alpha=0.5, color='C1', linewidth=2, label='Spectra with anomaly')
#     else:
#         plt.plot(spectrum[cut_from:, 0], ys, alpha=0.5, color='C1', linewidth=2,)


plt.plot(litdata[:, 0], litdata[:, 1]/0.575*0.2, '--', color='C0', linewidth=2, alpha=0.5, label='Nishimae et al. (2004)')
plt.plot(litdata_e3[:, 0], litdata_e3[:, 1]/0.45*0.2, '--', color='C2', linewidth=2, alpha=1, label='Ammer et al. (2012), E3+')
plt.plot(litdata_e12[:, 0], litdata_e12[:, 1]/0.15*0.2, color='C1', linewidth=2, alpha=0.5, label='Ammer et al. (2012), E12+')

carboca_highest_fom = np.load(data_folder + 'Yaroslav/mystery_prod/pink_spectrum.npy')
spectrum_with_lowest_diff = np.load(data_folder + 'Yaroslav/mystery_prod/nonpink_spectrum.npy')
plt.plot(carboca_highest_fom[:, 0], carboca_highest_fom[:, 1] - spectrum_with_lowest_diff[:, 1], color='deeppink', label='Unexpected compound')

nms = np.load('misc_scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_nms.npy')
spectrum = np.load('misc_scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_mean.npy')
perc1 = np.load('misc_scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_perc1.npy')
perc2 = np.load('misc_scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_perc2.npy')

plt.plot(nms, spectrum, label=f'Theory for compound   ', linestyle='-',
         color='black', zorder=10)
plt.fill_between(nms, perc1, perc2, alpha=0.4,
                 label=f'Theory, 68% confidence interval', color='grey', zorder=-10)

plt.xlabel('Wavelength, nm')
plt.ylabel('Absorbance')

plt.xlim(420, 700)
plt.ylim(-0.01, 0.3)

# plt.legend()

plt.tight_layout()
print(f'cwd: {os.getcwd()}')
plt.gcf().savefig('misc_scripts/figures/pink-spectrum-comparisons.png', dpi=300)

plt.show()
#
# spectra = spectra_unmasked
# # Principal component analysis of the background spectra
# nc = 3
# pca1 = pca(n_components=nc)
# nfeat1 = [spectrum[cut_from:, 1] - bkg_spectrum[cut_from:, 1] for spectrum in spectra]
# nfeat1 = [x - np.mean(x) for x in nfeat1]
# X1 = pca1.fit(nfeat1)
# expl_var_1 = X1.explained_variance_ratio_
# sv = X1.components_.T
# for spectrum in spectra:
#     plt.plot(spectrum[cut_from:, 0], spectrum[cut_from:, 1] - bkg_spectrum[cut_from:, 1], alpha=0.2, color='grey')
# # colors = ['C0', 'C1', 'C2']
# for i in range(nc):
#     plt.plot(spectra[0][cut_from:, 0], sv[:, i], label=f'PC {i+1}', alpha=1)
#
# # # Smooth the components and plot them
# # for i in range(nc):
# #     smoothed_component = savgol_filter(sv[:, i], polyorder=1, window_length=71)
# #     plt.plot(spectra[0][cut_from:, 0], smoothed_component, linestyle='--', label=f'smoothed PC {i+1}', color=colors[i])
# #     for_saving = np.zeros_like(bkg_spectrum[:, 1])
# #     for_saving[cut_from:] = smoothed_component
#
# plt.title('PCA components')
# plt.xlabel('Wavelength, nm')
# plt.ylabel('Absorbance units,\nspectra are shifted to be around zero')
# plt.show()

# # Save the components
# np.save(data_folder + f'multicomp-reactions/2023-03-20-run01/microspectrometer_data/background_model/component_{i}.npy',
#         for_saving)
