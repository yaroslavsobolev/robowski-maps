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
process_wellplate_spectra = importlib.import_module("uv-vis-absorption-spectroscopy.process_wellplate_spectra")

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'
sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')

litdata = np.loadtxt('misc-scripts/literature_spectra/carbocation_Yuichi_Nishimae_2004.txt', delimiter='\t')
litdata_e12 = np.loadtxt('misc-scripts/literature_spectra/ammer-sailer-riedle-2012-E12.txt', delimiter='\t')
litdata_e3 = np.loadtxt('misc-scripts/literature_spectra/ammer-sailer-riedle-2012-E3.txt', delimiter='\t')

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
foms = [spectrum[interm_wav_id, 1] - spectrum[interm_wav_id2, 1] for spectrum in spectra]
# find the spectrum with highest fom
highest_fom_id = np.argmax(foms)
carboca_highest_fom = spectra[highest_fom_id]
# subtract the mean of the last 100 points
carboca_highest_fom[:, 1] -= np.mean(carboca_highest_fom[-100:, 1])

carbocat_mask = np.array([spectrum[interm_wav_id, 1] - spectrum[interm_wav_id2, 1] > 0.037 for spectrum in spectra])
spectra_unmasked = np.array(spectra)[~carbocat_mask, :, :]
spectra_unmasked = [spectrum for spectrum in spectra_unmasked]

# f1 = plt.figure(dpi=300, figsize=(4, 4))

spectral_diffs = []
for spectrum in spectra_unmasked:
    # limits for comparison
    limit_wav_1 = wavelength_id_by_value(wavelengths, 430)
    limit_wav_2 = wavelength_id_by_value(wavelengths, 640)
    limit_wav_3 = wavelength_id_by_value(wavelengths, 410)
    limit_wav_4 = wavelength_id_by_value(wavelengths, 700)
    # compare the spectrum to carbocat_highest_fom where wav_id is below limit_wav_1 or above limit_wav_2
    # wave_id_mask = np.logical_or(np.arange(len(wavelengths)) < limit_wav_1, np.arange(len(wavelengths)) > limit_wav_2)
    # wave_id_mask = np.logical_and(wave_id_mask, np.arange(len(wavelengths)) > limit_wav_3)
    # wave_id_mask = np.logical_and(wave_id_mask, np.arange(len(wavelengths)) < limit_wav_4)
    spectrum[:, 1] = spectrum[:, 1] - np.mean(spectrum[-100:, 1])
    wave_id_mask = np.logical_and(np.arange(len(wavelengths)) < limit_wav_1, np.arange(len(wavelengths)) >= limit_wav_3)
    spectral_diff = np.sqrt(np.mean((spectrum[wave_id_mask, 1] - carboca_highest_fom[wave_id_mask, 1])**2))
    spectral_diffs.append(spectral_diff)

# plot the spectra_unmasked that have highest spectral_diff
index_of_spectrum_with_lowest_spectral_diff = np.argmin(spectral_diffs)
spectrum_with_lowest_diff = spectra_unmasked[index_of_spectrum_with_lowest_spectral_diff]
spectrum_with_lowest_diff[:, 1] *= 0.416 / 0.450 * 0.90

difference_spectrum = carboca_highest_fom[:, 1] - spectrum_with_lowest_diff[:, 1]
wavelengths = spectrum_with_lowest_diff[:, 0]
for_saving = np.stack((wavelengths, difference_spectrum)).T

# make a figure with two subplots, xharex=true
factor = 0.85
fig, ax = plt.subplots(1, sharex=True, figsize=(4.4*factor, 4*factor))
axs = [ax]
simpleaxis(axs[0])
# axs[0].fill_between(x=wavelengths, y1=carboca_highest_fom[:, 1], y2=spectrum_with_lowest_diff[:, 1], color='deeppink', alpha=0.4)
# axs[0].fill_between(x=wavelengths, y1=0, y2=spectrum_with_lowest_diff[:, 1], color='grey', alpha=0.4)

axs[0].plot(spectrum_with_lowest_diff[:, 0], spectrum_with_lowest_diff[:, 1], linestyle='--', label='Spectrum possible for a mix of substrate and products', color='black')
# plot carboca_highest_fom
axs[0].plot(carboca_highest_fom[:, 0], carboca_highest_fom[:, 1], color='deeppink', label='Experimental spectrum of a crude with unknown (pink) compound')

# add second subplot below the first, with sharex = true

# axs[1].plot(wavelengths, difference_spectrum, color='deeppink', label='Difference spectrum')
# axs[1].fill_between(x=wavelengths, y1=0, y2=difference_spectrum, color='black', alpha=0.4)
# save (wavelengths, difference_spectrum).T to numpy file
# np.save(data_folder + 'Yaroslav/mystery_prod/extracted_pink_spectrum.npy', for_saving)
np.save(data_folder + 'Yaroslav/mystery_prod/nonpink_spectrum.npy', spectrum_with_lowest_diff)
np.save(data_folder + 'Yaroslav/mystery_prod/pink_spectrum.npy', carboca_highest_fom)


nms = np.load('misc-scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_nms.npy')
spectrum = np.load('misc-scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_mean.npy')
perc1 = np.load('misc-scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_perc1.npy')
perc2 = np.load('misc-scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_perc2.npy')

# resample spectrum_with_lowest_diff[:, 0], spectrum_with_lowest_diff[:, 1] to nms
spectrum_with_lowest_diff_resampled = np.interp(nms, spectrum_with_lowest_diff[:, 0], spectrum_with_lowest_diff[:, 1])

plt.plot(nms, spectrum + spectrum_with_lowest_diff_resampled, label=f'Theoretical', linestyle='-',
         color='C0', zorder=-10)
plt.fill_between(nms, perc1 + spectrum_with_lowest_diff_resampled, perc2+spectrum_with_lowest_diff_resampled, alpha=0.2,
                 label=f'Theoretical 1$\sigma$-confidence interval', color='C0', zorder=-10)

# axs[0].legend()
# axs[1].legend()
axs[0].set_xlabel('Wavelength, nm')
axs[0].set_ylabel('Absorbance')
for ax in axs:
    ax.set_ylabel('Absorbance')
    ax.set_xlim(400, 700)
    ax.set_ylim(-0.01, 0.44)
plt.tight_layout()
fig.savefig(data_folder + 'Yaroslav/mystery_prod/pink_spectrum_with_dft.png', dpi=300)
plt.show()