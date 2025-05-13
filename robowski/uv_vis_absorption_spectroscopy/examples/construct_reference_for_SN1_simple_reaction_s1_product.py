from robowski.settings import *
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from scipy import interpolate
from scipy.optimize import curve_fit
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra


data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'
calibration_plate_folder = '2023-05-23_01-51-15__plate0000020__simple-reactions-2023-05-22-run01_calibration/'

sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
run_shortname = '2023-07-05-run01'
experiment_name = f'simple-reactions/{run_shortname}/'
calibrant_shortname = 'SN1Br01s1'
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'

do_plot=False
lower_limit_of_absorbance = 0.005
ref_concentration = 0.002585158

process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + 'references')
process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + 'background')
process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}')

# bkg_spectra = np.array(
#     sp.load_all_spectra(plate_folder=craic_folder + '2023-06-15_15-49-38__plate0000037__pure-dmf-bkg-test/'))

bkg_spectra = np.array(
    sp.load_all_spectra(plate_folder=craic_folder + '2023-05-23_01-51-15__plate0000020__simple-reactions-2023-05-22-run01_calibration/'))
bkg_spectrum = np.mean(bkg_spectra[:3,:], axis=0)
# bkg_spectrum[:, 1] = 0
if do_plot:
    plt.plot(bkg_spectrum[:, 1])
    plt.title('Bkg spectrum')
    plt.show()

# one_calibrant_df = df_excel.loc[df_excel['c#HBr'] == 0]

def load_spectrum_by_id(id):
    spectrum = sp.load_msp_by_id(
        plate_folder=craic_folder + calibration_plate_folder,
        well_id=id)
    spectrum[:, 1] -= bkg_spectrum[:, 1]
    return spectrum

ref_spectrum = load_spectrum_by_id(34)[:, 1]
ref_spectrum -= 0.47 * bkg_spectrum[:, 1]
ref_spectrum -= np.mean(ref_spectrum[-100:])
if do_plot:
    plt.plot(ref_spectrum)
    plt.title('Ref spectrum')
    plt.show()

wavelength_indices = np.arange(ref_spectrum.shape[0])
reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')


# advanced background model
background_model_folder=data_folder + 'multicomp-reactions/2023-03-20-run01/microspectrometer_data/background_model/'
background_interpolators = [interpolate.interp1d(wavelength_indices,
                                                 np.load(background_model_folder + f'component_{i}.npy'),
                                                 fill_value='extrapolate')
                            for i in range(2)]

# thresh_w_indices = [0, 25, 127, 2000]
# thresh_as = [0.67, 0.75, 1.6, 1.6]
# threshold_interpolator = interpolate.interp1d(thresh_w_indices, thresh_as, fill_value='extrapolate')
# concentrations = sorted(one_calibrant_df[concentration_column_name].to_list())
concentrations = np.loadtxt(f'{data_folder}simple-reactions/2023-05-22-run01/outVandC/SN1Br01s1.txt')
unique_concentrations = sorted(list(set(concentrations)))

process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}/concentration_fits')
cut_from = 79
coeffs = []
coeff_errs = []
for concentration in unique_concentrations:
    if concentration == 0:
        coeffs.append(0)
        coeff_errs.append(0)
        continue

    well_ids = np.where(concentrations == concentration)[0]
    slopes_for_this_concentration = []
    for well_id in well_ids:
        target_spectrum = load_spectrum_by_id(well_id)[:, 1]
        mask = wavelength_indices > cut_from
        mask = np.logical_and(mask, target_spectrum > target_spectrum[-1] + lower_limit_of_absorbance)

        def func(xs, a, b, c, d, e):
            return a * reference_interpolator(xs) + b \
                       + c*xs + d * background_interpolators[0](xs) + e * background_interpolators[1](xs)

        p0 = (concentration / ref_concentration, 0, 0, 0, 0)
        bounds = ([-1e-10, -np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf])
        popt, pcov = curve_fit(func, wavelength_indices[mask], target_spectrum[mask],
                               p0=p0, bounds=bounds)
        perr = np.sqrt(np.diag(pcov))
        slopes_for_this_concentration.append(popt[0])

        fig1 = plt.figure(1)
        plt.plot(target_spectrum, label='data', color='C0', alpha=0.5)
        mask_illustration = np.ones_like(target_spectrum) * np.max(target_spectrum)
        mask_illustration[mask] = 0
        plt.fill_between(x=wavelength_indices, y1=0, y2=mask_illustration, color='yellow', alpha=0.3,
                         label='ignored (masked) data')
        plt.plot(func(wavelength_indices, *popt), color='r', label='fit', alpha=0.5)
        plt.plot(func(wavelength_indices, popt[0], 0, 0, 0, 0), color='C1', label='reference', alpha=0.5)
        plt.ylim(-0.1,
                 np.max((func(wavelength_indices, *popt)[mask])) * 2)
        plt.title(
            f"conc {concentration}, well {well_id}")
        plt.legend()
        fig1.savefig(
            calibration_folder + f"references/{calibrant_shortname}/concentration_fits/{well_id}_fit.png")
        if do_plot:
            plt.show()
        else:
            plt.clf()

    slope = np.mean(np.array(slopes_for_this_concentration))
    slope_error = perr[0]
    coeffs.append(slope)
    coeff_errs.append(slope_error)

fig3 = plt.figure(3)
plt.loglog(coeffs, unique_concentrations, 'o-')
plt.xlabel('Fit coefficients')
plt.ylabel('Concentrations, mol/liter')
fig3.savefig(calibration_folder + f"references/{calibrant_shortname}/concentration-vs-coeff.png", dpi=300)
if do_plot:
    plt.show()
else:
    plt.clf()

coeff_to_concentration_interpolator = interpolate.interp1d(coeffs, unique_concentrations,
                                                           fill_value='extrapolate')

np.save(calibration_folder + f'references/{calibrant_shortname}/bkg_spectrum.npy', bkg_spectrum)
np.save(calibration_folder + f'background//bkg_spectrum.npy', bkg_spectrum)
np.save(calibration_folder + f'references/{calibrant_shortname}/ref_spectrum.npy', ref_spectrum)
np.save(calibration_folder + f'references/{calibrant_shortname}/interpolator_coeffs.npy', np.array(coeffs))
np.save(calibration_folder + f'references/{calibrant_shortname}/interpolator_concentrations.npy',
        unique_concentrations)