import importlib
import os

import numpy as np
from matplotlib import pyplot as plt

SpectraProcessor = importlib.import_module('uv-vis-absorption-spectroscopy.process_wellplate_spectra').SpectraProcessor
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

sp = SpectraProcessor(folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220

# well_id = 44
well_id = 33
# Condition 154
# plate_folder = data_folder + 'BPRF/2024-01-08-run01/nanodrop_spectra/2024-01-10_12-51-07_UV-Vis_plate_71.csv'
# plate_folder = data_folder + 'BPRF/2024-01-08-run02/nanodrop_spectra/2024-01-10_17-10-28_UV-Vis_plate_61.csv'
# plate_folder = data_folder + 'BPRF/2024-01-17-run01/nanodrop_spectra/2024-01-19_12-22-47_UV-Vis_plate_66.csv'
plate_folder = data_folder + 'BPRF/2024-01-17-run01/nanodrop_spectra/2024-01-19_12-22-47_UV-Vis_plate_66.csv'
spectrum1 = sp.load_msp_by_id(
    plate_folder=plate_folder,
    well_id=well_id)[:, 1]

# plate_folder = data_folder + 'BPRF/2024-01-08-run01/nanodrop_spectra/2024-01-10_13-48-13_UV-Vis_plate_73.csv'
# plate_folder = data_folder + 'BPRF/2024-01-08-run02/nanodrop_spectra/2024-01-10_17-55-20_UV-Vis_plate_66.csv'
# plate_folder = data_folder + 'BPRF/2024-01-17-run01/nanodrop_spectra/2024-01-19_13-00-17_UV-Vis_plate_67.csv'
plate_folder = data_folder + 'BPRF/2024-01-17-run01/nanodrop_spectra/2024-01-19_12-22-47_UV-Vis_plate_66.csv'
spectrum2 = sp.load_msp_by_id(
    plate_folder=plate_folder,
    well_id=well_id)[:, 1]


# dedicated volume check
plate_folder = data_folder + 'BPRF/volume_check_20240119/nanodrop_spectra/2024-01-19_20-46-33_UV-Vis_volume_check.csv'
bkg = np.mean(np.array([sp.load_msp_by_id(
    plate_folder=plate_folder,
    well_id=i)[:, 1] for i in [15]]), axis=0)


spectrum1 = np.mean(np.array([sp.load_msp_by_id(
    plate_folder=plate_folder,
    well_id=i)[:, 1] for i in [0, 1, 2, 3]]), axis=0)
spectrum1 -= bkg
spectrum1 = spectrum1 - np.mean(spectrum1[300:])
spectrum2 = np.mean(np.array([sp.load_msp_by_id(
    plate_folder=plate_folder,
    well_id=i)[:, 1] for i in [5,6,7,8,9]]), axis=0)
spectrum2 -= bkg
spectrum2 = spectrum2 - np.mean(spectrum2[300:])

plt.plot(spectrum1)
plt.plot(spectrum2)

#
# plate_folder = data_folder + 'BPRF/volume_check_20240119/nanodrop_spectra/2024-01-19_20-46-33_UV-Vis_volume_check.csv'
# spectrum1 = np.mean(np.array([sp.load_msp_by_id(
#     plate_folder=plate_folder,
#     well_id=i)[:, 1] for i in [5,6,7,8,9]]), axis=0)
# spectrum1 -= bkg
# spectrum1 = spectrum1 - np.mean(spectrum1[300:])
# spectrum2 = np.mean(np.array([sp.load_msp_by_id(
#     plate_folder=plate_folder,
#     well_id=i)[:, 1] for i in [10, 12, 13, 14]]), axis=0)
# spectrum2 -= bkg
# spectrum2 = spectrum2 - np.mean(spectrum2[300:])

# plt.plot(spectrum1)
# plt.plot(spectrum2)
#
# plt.plot(220 + np.arange(spectrum1.shape[0]), spectrum1, label='Measured at dilution 20x')
# plt.plot(220 + np.arange(spectrum1.shape[0]), spectrum2*200/20, label='Measured at dilution 200x, then multiplied by 10')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorbance')
plt.legend()

plt.show()