import os
import numpy as np
import pandas as pd
import importlib
from matplotlib import pyplot as plt
import logging
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import pickle

from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# set level to info
logging.basicConfig(level=logging.INFO)
st = importlib.import_module('uv_vis_absorption_spectroscopy.spectraltools')
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

def average_all_cycles(cary_file, column_name_without_cycle_suffix):
    # read first line of cary_file and split by comma delimiter
    with open(cary_file, 'r') as file:
        first_line = file.readline()
    column_names = first_line.split(',')
    column_names = [x for x in column_names if x != '']
    # find all the columns that start with column_name_without_cycle_suffix and end with '_1', '_2' and so on,
    # until the first number that is not present
    spectra = []
    for cycle_number in range(1, 20):
        if f'{column_name_without_cycle_suffix}_{cycle_number}' not in column_names:
            # print(f'Cycle {cycle_number} not found in {cary_file}')
            break
        wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=f'{column_name_without_cycle_suffix}_{cycle_number}')
        if cycle_number == 1:
            wav_true = np.copy(wav)
        # print the smallest wav
        # print(f'{cycle_number}: minwav: {min(wav)}')
        if min(wav) > 201:
            print(f'min wav: {min(wav)}')
            break
        spectra.append(spec)
    print(f'{column_name_without_cycle_suffix}, largest_cycle_number = ', cycle_number)
    # average the spectra
    if len(spectra) == 1:
        average_spectrum = spectra[0]
    else:
        average_spectrum = np.mean(np.array(spectra), axis=0)
    return wav_true, average_spectrum


# Acetic acid is added in microliters (from 0 to 60 uL) of the following stock:
# 85.5 uL of acetic acid
# 4.914 mL of Ethanol
# The concentration of acetic acid in stock is 0.2989 mol/L
# The maximum concentration of acetic acid in cuvette, achieved with 60 uL of acetic acid/ethanol stock is 0.005978 mol/L
# Given that the smallest dilution factor is 20x, 20 times the 0.005978 mol/L is 0.11956 mol/L

cary_file = data_folder + 'Yaroslav/Hantzsch_acetic_influence/HRP01_acetic_acid_and_pure_ethanol.csv'
cary_column_name = f'pure_ethanol_60ul_acetic_acid_rep1_1'
wav_acetic, spec_acetic = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
plt.plot(wav_acetic, spec_acetic)
plt.show()
spec_acetic_interp = interpolate.interp1d(wav_acetic, spec_acetic, kind='linear', fill_value='extrapolate')

cary_file = data_folder + 'Yaroslav/Hantzsch_acetic_influence/Run2_bb017_with_acetic_acid_and_pure_ethanol_1.csv'
# cary_column_name = f'aminal_bb017_4.9ug_per_mL_and_0_acetic_acid_rep1_1'
# wav, _ = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
acetic_acid_volumes = np.array([0, 1, 2] + list(range(5, 65, 5)))
# remove the 25 value
acetic_acid_volumes = acetic_acid_volumes[acetic_acid_volumes != 25]
acetic_acid_concentrations = acetic_acid_volumes / 60 * 0.005978

cary_column_name = f'fresh_bb017_20ugPERmL_and_0uL_acetic_rep1'
wav, zero_concentration_spectrum = average_all_cycles(cary_file, cary_column_name)
# make zero_concentration_spectrum_interpolator
zero_concentration_spectrum_interpolator = interpolate.interp1d(wav, zero_concentration_spectrum, kind='linear', fill_value='extrapolate')
# fit the zero-concentration spectrum to the reference loaded from file, and get the fitting coefficient
calibration_folder=data_folder + 'BPRF/2024-01-17-run01/' + 'microspectrometer_data/calibration/'
calibrant_shortname = 'bb017'
wavelengths_of_reference = np.load(calibration_folder + f'references/{calibrant_shortname}/bkg_spectrum.npy')[:, 0]
ref_spectrum = np.load(calibration_folder + f'references/{calibrant_shortname}/ref_spectrum.npy')
# resample the ref_spectrum into the wav
# ref_spectrum_resampled = np.interp(wav, wavelengths_of_reference, ref_spectrum)
# make a fit of the zero-concentration spectrum multiplied by fitting scale factor to the ref_spectrum with curve_fit
# plt.plot(wavelengths_of_reference, ref_spectrum, '--', label='reference')
# plt.plot(wav, ref_spectrum_resampled, '--', label='reference resampled')



# Resample and save the calibration for acetic acid
acetic_spec_for_saving = interpolate.interp1d(wav_acetic, spec_acetic, kind='linear', fill_value='extrapolate')(wavelengths_of_reference)
# plt.plot(wavelengths_of_reference, acetic_spec_for_saving, label='acetic acid', alpha=0.5)
# plt.plot(wav_acetic, spec_acetic, label='acetic acid original', alpha=0.5)
# plt.legend()
# plt.show()
# csv_filename_for_acetic = data_folder + 'BPRF/2024-01-17-run01/calibrations/acetic_acid_2.csv'
csv_filename_for_acetic = data_folder + 'BPRF/2024-01-17-run01/calibrations/acetic_acid_2.csv'
np.savetxt(csv_filename_for_acetic, np.column_stack((wavelengths_of_reference, acetic_spec_for_saving)), delimiter=',', header=',1')

def func(x, a):
    return a * zero_concentration_spectrum_interpolator(x)

popt, _ = curve_fit(func, wavelengths_of_reference, ref_spectrum)
scaling_factor = popt[0]
print(f'scaling factor = {scaling_factor}')
plt.plot(wav, func(wav, *popt), label='fit')
plt.plot(wavelengths_of_reference, ref_spectrum, label='ref')
plt.legend()
# plt.plot(wav, zero_concentration_spectrum * scaling_factor, 'o', label='zero concentration')
plt.show()

wavelengths_grid, acetic_grid = np.meshgrid(wavelengths_of_reference, acetic_acid_volumes, indexing='ij')
wavelengths_grid = wavelengths_grid.astype(float)
acetic_grid = acetic_grid.astype(float)

spectral_2d_grid = np.zeros_like(acetic_grid, dtype=float)

for acetic_acid_volume_added in acetic_acid_volumes:
    cary_column_name = f'fresh_bb017_20ugPERmL_and_{acetic_acid_volume_added}uL_acetic_rep1'
    wav, spec = average_all_cycles(cary_file, cary_column_name)
    spec_acetic_res = spec_acetic_interp(wav)
    spec = spec - spec_acetic_res / 60 * acetic_acid_volume_added
    spec = spec * scaling_factor
    spec = interpolate.interp1d(wav, spec, kind='linear', fill_value='extrapolate')(wavelengths_of_reference)
    # for wavelengths above 320, replace with values of ref_spectrum
    mask = wavelengths_of_reference > 320
    spec[mask] = ref_spectrum[mask]
    # spec = spec - np.mean(spec[-100:])
    # spec[spec < 0.04] = savgol_filter(spec, 9, 3)[spec < 0.04]
    plt.plot(wavelengths_of_reference, spec, label=f'{acetic_acid_volume_added/ 60 * 0.005978*1000:.2f} mM acetic acid', color=plt.cm.viridis(acetic_acid_volume_added / 60))
    acetic_index = np.where(acetic_acid_volumes == acetic_acid_volume_added)[0][0]
    spectral_2d_grid[:, acetic_index] = spec

#the data structure is
array_of_wavelengths = wavelengths_of_reference
data_for_pickling = [array_of_wavelengths, acetic_acid_concentrations, spectral_2d_grid]
# pickle to data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/bb017/acetic_acid_influence.pkl'
# with open(data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/bb017/acetic_acid_influence.pkl', 'wb') as f:
#     pickle.dump(data_for_pickling, f)

# make 2d interpolator
interp = RegularGridInterpolator((array_of_wavelengths, acetic_acid_concentrations), spectral_2d_grid,
                                 bounds_error=False, fill_value=None)

# evaluate the spectrum at one concentration of acetic acid -- just for fun
c = 3 / 60 * 0.005978
spec = interp((wavelengths_of_reference, np.ones_like(wavelengths_of_reference)*c))
plt.plot(wavelengths_of_reference, spec, label=f'{c*1000:.2f} mM acetic acid', color='red')

plt.xlabel('Wavelength, nm')
plt.ylabel('Absorbance')
plt.legend()
plt.show()