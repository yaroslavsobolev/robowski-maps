import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import importlib
from scipy import interpolate
from scipy.optimize import curve_fit
process_wellplate_spectra = importlib.import_module("uv-vis-absorption-spectroscopy.process_wellplate_spectra")

dilution_factor = 1

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
run_shortname = '2023-08-21-run01'
experiment_name = f'simple-reactions/{run_shortname}/'
calibrant_shortname = 'HBr'
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'
plate_folder = data_folder + 'nanodrop-spectrophotometer-measurements/reference_for_simple_SN1/2023-09-07_22-46-02_E1_ref_and_etoh_hbr_aac.csv'
one_calibrant_df = pd.read_csv(data_folder + 'nanodrop-spectrophotometer-measurements/reference_for_simple_SN1/'
                                             '2023_09_07_UV-Vis_E1_ref_and_etoh_hbr_aac.txt')
one_calibrant_df['nanodrop_col_name'] = one_calibrant_df['nanodrop_col_name'].astype(str)
one_calibrant_df = one_calibrant_df[one_calibrant_df['substance'] == 'HBR_48%']

ref_concentration=0.0035759999999999998/dilution_factor
do_plot=True
lower_limit_of_absorbance = 0.007

process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + 'references')
process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + 'background')
process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}')

# load background from different file
bkg_spectrum = np.load(calibration_folder + f'background//bkg_spectrum.npy')

nanodrop_df = sp.load_nanodrop_csv_for_one_plate(plate_folder)

# This piece of code gets rid of Yankai's excessive creativity in naming columns
for i in range(28):
    # average columns with name like "1", "1.1", "1.2" and write to column "1"
    nanodrop_df[f'{i}'] = nanodrop_df[[f'{i}', f'{i}.1', f'{i}.2']].mean(axis=1)
    # delete columns "1.1" and "1.2"
    nanodrop_df.drop(columns=[f'{i}.1', f'{i}.2'], inplace=True)

wavelengths = nanodrop_df['wavelength'].to_numpy()

def load_spectrum_by_df_row(row):
    wavelengths = nanodrop_df['wavelength'].to_numpy()
    absorbances = nanodrop_df[row['nanodrop_col_name']].to_numpy()
    spectrum = np.array([wavelengths, absorbances]).T
    spectrum[:, 1] -= bkg_spectrum[:, 1]
    spectrum[:, 1] -= np.mean(spectrum[-100:, 1])
    return spectrum

concentration_column_name = 'concentration'

# make sure that only one well for this calibrant has concentration equal to ref_concentration
assert one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == ref_concentration].shape[0] == 1
ref_spectrum = load_spectrum_by_df_row(
    one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == ref_concentration].iloc[0])[:, 1]
ref_spectrum -= np.mean(ref_spectrum[-100:])
if do_plot:
    plt.plot(wavelengths, ref_spectrum)
    plt.title('Ref spectrum')
    plt.show()

wavelength_indices = np.arange(ref_spectrum.shape[0])
reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')

concentrations = sorted(one_calibrant_df[concentration_column_name].to_list())
concentrations = [0, concentrations[2], concentrations[3]]

process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}/concentration_fits')
# cut_from = 115
cut_from=50
coeffs = []
coeff_errs = []
for concentration in concentrations:
    if concentration == 0:
        coeffs.append(0)
        coeff_errs.append(0)
        continue

    df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]
    target_spectrum = load_spectrum_by_df_row(df_row_here)[:, 1]
    mask = wavelength_indices > cut_from
    # mask = np.logical_and(mask, target_spectrum > np.min(target_spectrum) + lower_limit_of_absorbance)

    def func(xs, a, b):
        return a * reference_interpolator(xs) + b

    p0 = (concentration / ref_concentration, 0)
    bounds = ([-1e-10, -np.inf], [np.inf, np.inf])
    popt, pcov = curve_fit(func, wavelength_indices[mask], target_spectrum[mask],
                           p0=p0, bounds=bounds)
    # sigma=noise_std*np.ones_like(target_spectrum),
    # absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    slope = popt[0]
    slope_error = perr[0]
    coeffs.append(slope)
    coeff_errs.append(slope_error)

    fig1 = plt.figure(1)
    plt.plot(target_spectrum, label='data', color='C0', alpha=0.5)
    mask_illustration = np.ones_like(target_spectrum) * np.max(target_spectrum)
    mask_illustration[mask] = 0
    plt.fill_between(x=wavelength_indices, y1=0, y2=mask_illustration, color='yellow', alpha=0.3,
                     label='ignored (masked) data')
    plt.plot(func(wavelength_indices, *popt), color='r', label='fit', alpha=0.5)
    plt.plot(func(wavelength_indices, popt[0], 0), color='C1', label='reference', alpha=0.5)
    plt.ylim(-0.03,
             np.max((func(wavelength_indices, *popt)[mask])) * 2)
    plt.title(
        f"conc {df_row_here[concentration_column_name]}, well {df_row_here['nanodrop_col_name']}")
    plt.legend()
    fig1.savefig(
        calibration_folder + f"references/{calibrant_shortname}/concentration_fits/{df_row_here[concentration_column_name]}_fit.png")
    if do_plot:
        plt.show()
    else:
        plt.clf()

fig3 = plt.figure(3)
plt.loglog(coeffs, concentrations, 'o-')
plt.xlabel('Fit coefficients')
plt.ylabel('Concentrations, mol/liter')
fig3.savefig(calibration_folder + f"references/{calibrant_shortname}/concentration-vs-coeff.png", dpi=300)
if do_plot:
    plt.show()
    plt.plot(coeffs, concentrations, 'o-')
    plt.xlabel('Fit coefficients')
    plt.ylabel('Concentrations, mol/liter')
    plt.show()
else:
    plt.clf()

coeff_to_concentration_interpolator = interpolate.interp1d(coeffs, concentrations,
                                                           fill_value='extrapolate')

np.save(calibration_folder + f'references/{calibrant_shortname}/bkg_spectrum.npy', bkg_spectrum)
np.save(calibration_folder + f'background//bkg_spectrum.npy', bkg_spectrum)
np.save(calibration_folder + f'references/{calibrant_shortname}/ref_spectrum.npy', ref_spectrum)
np.save(calibration_folder + f'references/{calibrant_shortname}/interpolator_coeffs.npy', np.array(coeffs))
np.save(calibration_folder + f'references/{calibrant_shortname}/interpolator_concentrations.npy',
        concentrations)