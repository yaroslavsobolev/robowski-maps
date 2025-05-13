import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import importlib
from scipy import interpolate
from scipy.optimize import curve_fit
process_wellplate_spectra = importlib.import_module("uv_vis_absorption_spectroscopy.process_wellplate_spectra")

dilution_factor = 200

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'

sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
run_shortname = '2023-07-05-run01'
experiment_name = f'simple-reactions/{run_shortname}/'
calibrant_shortname = 'SN1OH01'
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'

ref_concentration=0.01146145002002/dilution_factor
do_plot=True
do_reference_refinements=True
lower_limit_of_absorbance = 0.02

process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + 'references')
process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + 'background')
process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}')

# calibration_sequence_df = pd.read_csv(calibration_folder + 'calibration_sequence_dataframe.csv')
# one_calibrant_df = calibration_sequence_df.loc[calibration_sequence_df['shortname'] == calibrant_shortname]

# bkg_row = one_calibrant_df.loc[one_calibrant_df['concentration'] == 0].iloc[0]

bkg_spectra = np.array(
    sp.load_all_spectra(plate_folder=craic_folder + '2023-06-15_15-49-38__plate0000037__pure-dmf-bkg-test/'))
bkg_spectrum = np.mean(bkg_spectra, axis=0)
bkg_spectrum[:, 1] = 0
# if do_plot:
#     plt.plot(bkg_spectrum[:, 1])
#     plt.title('Bkg spectrum')
#     plt.show()

# make a dataframe with data for one calibrant. This means that we find conditions where this calibrant has
# non-zero concentraion, while all other substances have zero concentration.
# load excel with conditions (volumes) used by the pipetter
# df_excel = pd.read_excel(data_folder + experiment_name + f'{run_shortname}.xlsx', sheet_name=0)

df_excel = pd.read_csv(data_folder + experiment_name + f'/results/run_structure.csv', header=1)

# load excel table with stock solution compositions
df_stock = pd.read_excel(data_folder + experiment_name + f'outVandC/stock_solutions.xlsx', sheet_name=0)
# iterate over rows in df_excel
for substance in df_stock.columns[1:]:
    df_excel[f'c#{substance}'] = 0

for index, row in df_excel.iterrows():
    # iterate over columns in df_excel
    substances = df_stock.columns[1:]
    net_volume = 0
    for stock_solution_name in df_excel.columns:
        # if this is a column with a stock solution
        if stock_solution_name in df_stock['stock_solution'].to_list():
            # get the stock solution concentration
            stock_concentration = df_stock.loc[df_stock['stock_solution'] == stock_solution_name]
            # get the volume of this stock solution used by the pipetter
            volume = row[stock_solution_name]
            for substance in df_stock.columns[1:]:
                this_stock_solution_concentration = stock_concentration[substance].iloc[0]
                df_excel.loc[index, f'c#{substance}'] = df_excel.loc[index, f'c#{substance}'] + volume * this_stock_solution_concentration
            net_volume += volume
    for substance in df_stock.columns[1:]:
        df_excel.loc[index, f'c#{substance}'] = df_excel.loc[index, f'c#{substance}'] / net_volume

# divide all the concentrations by dilution factor, because here we are loading spectra from the diluted plate
df_excel.loc[:, df_excel.columns.str.startswith('c#')] = df_excel.loc[:, df_excel.columns.str.startswith('c#')] / dilution_factor

one_calibrant_df = df_excel.loc[df_excel['c#HBr'] == 0]

def load_spectrum_by_df_row(row):
    spectrum = sp.load_msp_by_id(
        plate_folder=craic_folder + row['craic_folder'] + '/',
        well_id=row['vial_id'])
    # spectrum[:, 1] -= bkg_spectrum[:, 1]
    return spectrum

concentration_column_name = f'c#{calibrant_shortname}'

# make sure that only one well for this calibrant has concentration equal to ref_concentration
assert one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == ref_concentration].shape[0] == 1
ref_spectrum = load_spectrum_by_df_row(
    one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == ref_concentration].iloc[0])[:, 1]
ref_spectrum -= np.mean(ref_spectrum[-100:])
if do_plot:
    plt.plot(ref_spectrum)
    plt.title('Ref spectrum')
    plt.show()

wavelength_indices = np.arange(ref_spectrum.shape[0])
reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')

thresh_w_indices = [0, 25, 127, 2000]
thresh_as = [0.67, 0.75, 1.6, 1.6]
threshold_interpolator = interpolate.interp1d(thresh_w_indices, thresh_as, fill_value='extrapolate')
concentrations = sorted(one_calibrant_df[concentration_column_name].to_list())
process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}/concentration_fits')
# cut_from = 115
cut_from=79
coeffs = []
coeff_errs = []
for concentration in concentrations:
    if concentration == 0:
        coeffs.append(0)
        coeff_errs.append(0)
        continue

    df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]
    target_spectrum = load_spectrum_by_df_row(df_row_here)[:, 1]
    mask = np.logical_and(target_spectrum < threshold_interpolator(wavelength_indices),
                          wavelength_indices > cut_from)
    mask = np.logical_and(mask, target_spectrum > np.min(target_spectrum) + lower_limit_of_absorbance)

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
    plt.ylim(-0.3,
             np.max((func(wavelength_indices, *popt)[mask])) * 2)
    plt.title(
        f"conc {df_row_here[concentration_column_name]}, well {df_row_here['vial_id']}, diluted plate {df_row_here['diluted_plate_id']:04d}")
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