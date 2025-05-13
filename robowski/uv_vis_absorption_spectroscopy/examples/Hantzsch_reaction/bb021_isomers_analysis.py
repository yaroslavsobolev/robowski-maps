import importlib
import os

import numpy as np
from matplotlib import pyplot as plt

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

process_wellplate_spectra = importlib.import_module('uv_vis_absorption_spectroscopy.process_wellplate_spectra')

sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220
well_id = 7
# substances_for_fitting = ['bb021', 'bb021_f2']
substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'dm36', 'dm37', 'dm40_12', 'dm40_10', 'ethyl_acetoacetate', 'EAB', 'bb017', 'bb021', 'dm70', 'dm053', 'dm088_4', 'bb021_f2']
# cut_from = 40
cut_from = 1
plate_folder = data_folder + 'BPRF/2024-01-17-run01/calibrations/2024_04_0_UV-Vis_bb021_F2.csv'
spectrum1 = sp.load_single_nanodrop_spectrum(plate_folder=plate_folder, well_id=well_id)[:, 1]
spectrum2 = spectrum1

plt.plot(spectrum1 * 1)
# plt.plot(spectrum2 * 200)
print('len of spectrum1', len(spectrum1))
plt.show()

# set logging level to debug

# spectrum1 = spectrum2 * 1.1

# concentrations = sp.spectrum_to_concentration(target_spectrum_input=spectrum2,
#                                                    calibration_folder=data_folder + 'BPRF/2023-11-08-run01/' + 'microspectrometer_data/calibration/',
#                                                    calibrant_shortnames=substances_for_fitting,
#                                                    background_model_folder=data_folder + 'simple-reactions/2023-09-06-run01/microspectrometer_data/background_model/',
#                                                    upper_bounds=[np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf],
#                                                    do_plot=True, cut_from=cut_from,
#                                                    ignore_abs_threshold=True, ignore_pca_bkg=True)

concentrations = sp.multispectrum_to_concentration(target_spectrum_inputs=[spectrum1, spectrum2],
                                                   dilution_factors=[20, 20],
                                                   calibration_folder=data_folder + 'BPRF/2024-01-17-run01/' + 'microspectrometer_data/calibration/',
                                                   calibrant_shortnames=substances_for_fitting,
                                                   background_model_folder=data_folder + 'BPRF/cross_conamination_and_backgound_test/ethanol_background_model/',
                                                   upper_bounds=[np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
                                                                 np.inf, np.inf, np.inf],
                                                   do_plot=True, cut_from=cut_from, cut_to=250,
                                                   ignore_abs_threshold=False, ignore_pca_bkg=False,
                                                   plot_calibrant_references=True,
                                                   upper_limit_of_absorbance=0.95)