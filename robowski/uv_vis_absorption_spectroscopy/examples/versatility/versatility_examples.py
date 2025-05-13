import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import importlib
from scipy import interpolate
from scipy.optimize import curve_fit
process_wellplate_spectra = importlib.import_module("uv_vis_absorption_spectroscopy.process_wellplate_spectra")
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

def construct_calibrant(
                        cut_from,
                        lower_limit_of_absorbance,
                        concentration_column_name,
                        do_plot,
                        experiment_name,
                        calibration_source_filename,
                        calibrant_shortnames,
                        ref_concentrations,
                        max_concentrations,
                        custom_bkg_spectrum_npy_file=None,
                        no_right_edge_subtraction=False,
                        ):

    plate_folder = f'{data_folder}{experiment_name}{calibration_source_filename}.csv'
    all_calibrants_df = pd.read_csv(f'{data_folder}{experiment_name}{calibration_source_filename}.txt')

    sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                         '2022-12-01/interpolator-dataset/')
    sp.nanodrop_lower_cutoff_of_wavelengths = 220

    calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'

    process_wellplate_spectra.create_folder_unless_it_exists(data_folder + experiment_name + 'microspectrometer_data')
    process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder)
    process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + 'references')
    process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + 'background')

    all_calibrants_df['nanodrop_col_name'] = all_calibrants_df['nanodrop_col_name'].astype(str)

    # load background from different file
    # bkg_spectrum = np.load(calibration_folder + f'background/bkg_spectrum.npy')
    nanodrop_df = sp.load_nanodrop_csv_for_one_plate(plate_folder)
    wavelengths = nanodrop_df['wavelength'].to_numpy()

    # load background from the row where concentration is 0
    if custom_bkg_spectrum_npy_file is not None:
        bkg_spectrum = np.load(custom_bkg_spectrum_npy_file)
    else:
        col_name_with_background = all_calibrants_df.loc[all_calibrants_df[concentration_column_name] == 0].iloc[0]['nanodrop_col_name']
        absorbances = nanodrop_df[col_name_with_background].to_numpy()
        bkg_spectrum = np.array([wavelengths, absorbances]).T

    if do_plot:
        # plot bkg_spectrum
        plt.plot(bkg_spectrum[:, 0], bkg_spectrum[:, 1])
        plt.title('Background spectrum')
        plt.show()

    def load_spectrum_by_df_row(row):
        wavelengths = nanodrop_df['wavelength'].to_numpy()
        absorbances = nanodrop_df[row['nanodrop_col_name']].to_numpy()
        spectrum = np.array([wavelengths, absorbances]).T
        spectrum[:, 1] -= bkg_spectrum[:, 1]
        if not no_right_edge_subtraction:
            spectrum[:, 1] -= np.mean(spectrum[-100:, 1])
        return spectrum

    np.save(calibration_folder + f'background//bkg_spectrum.npy', bkg_spectrum)

    def reference_for_one_calibrant(calibrant_shortname, ref_concentration, min_concentration=0, max_concentration=1000,
                                    do_plot=True):
        process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}')
        one_calibrant_df = all_calibrants_df[all_calibrants_df['substance'] == calibrant_shortname]

        # make sure that only one well for this calibrant has concentration equal to ref_concentration
        assert one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == ref_concentration].shape[0] == 1
        ref_spectrum = load_spectrum_by_df_row(
            one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == ref_concentration].iloc[0])[:, 1]
        if not no_right_edge_subtraction:
            ref_spectrum -= np.mean(ref_spectrum[-100:])
        if do_plot:
            plt.plot(wavelengths, ref_spectrum)
            plt.title('Ref spectrum')
            plt.show()

        wavelength_indices = np.arange(ref_spectrum.shape[0])
        reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')

        concentrations = sorted([0] + one_calibrant_df[concentration_column_name].to_list())
        concentrations = [x for x in concentrations if min_concentration <= x <= max_concentration]

        process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}/concentration_fits')

        coeffs = []
        coeff_errs = []
        spectra = []
        for concentration in concentrations:
            if concentration == 0:
                coeffs.append(0)
                coeff_errs.append(0)
                continue

            df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]
            target_spectrum = load_spectrum_by_df_row(df_row_here)[:, 1]
            spectra.append(np.copy(target_spectrum))
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

        fig3, axarr = plt.subplots(1, 2, figsize=(10, 5))
        ax1, ax2 = axarr
        ax2.plot(coeffs, concentrations, 'o-')
        ax2.set_xlabel('Best-fit scale coefficient')
        ax2.set_ylabel('Concentration, mol/liter')
        fig3.savefig(calibration_folder + f"references/{calibrant_shortname}/concentration-vs-coeff.png", dpi=300)

        if do_plot:
            plt.show()
            plt.loglog(coeffs, concentrations, 'o-')
            plt.xlabel('Best-fit scale coefficient')
            plt.ylabel('Concentration, mol/liter')
            plt.show()
        else:
            plt.clf()

        coeff_to_concentration_interpolator = interpolate.interp1d(coeffs, concentrations,
                                                                   fill_value='extrapolate')
        # Saving
        np.save(calibration_folder + f'references/{calibrant_shortname}/bkg_spectrum.npy', bkg_spectrum)
        np.save(calibration_folder + f'background//bkg_spectrum.npy', bkg_spectrum)
        np.save(calibration_folder + f'references/{calibrant_shortname}/ref_spectrum.npy', ref_spectrum)
        np.save(calibration_folder + f'references/{calibrant_shortname}/interpolator_coeffs.npy', np.array(coeffs))
        np.save(calibration_folder + f'references/{calibrant_shortname}/interpolator_concentrations.npy',
                concentrations)

    for i, calibrant_shortname in enumerate(calibrant_shortnames):
        reference_for_one_calibrant(calibrant_shortname, ref_concentration=ref_concentrations[i],
                                    max_concentration=max_concentrations[i], do_plot=do_plot)


def process_plate(sp, dilution_factor,
                  plate_folder,
                  well_ids,
                  calibrant_shortnames,
                  calibration_folder,
                  experiment_name,
                  cut_from,
                  cut_to,
                  do_plot=False,
                  use_line=True,
                  ignore_bkg_pca=True,
                  get_errors_from_fit=False,
                  std_calib = 0.0105,
                  upper_bounds='auto'):
    if upper_bounds == 'auto':
        upper_bounds = [np.inf] * len(calibrant_shortnames)
    plate_name = plate_folder.split('/')[-1]
    df = pd.DataFrame(columns=['well_id'] + calibrant_shortnames, dtype=object)
    for well_id in well_ids:
        spectrum = sp.load_msp_by_id(
            plate_folder=plate_folder,
            well_id=well_id)[:, 1]
        process_wellplate_spectra.create_folder_unless_it_exists(data_folder + experiment_name + 'results')
        process_wellplate_spectra.create_folder_unless_it_exists(data_folder + experiment_name + 'results/uv-vis-fits')
        concentrations_here = sp.spectrum_to_concentration(target_spectrum_input=spectrum,
                                                             calibration_folder=calibration_folder,
                                                             calibrant_shortnames=calibrant_shortnames,
                                                             fig_filename=data_folder + experiment_name + f'results/uv-vis-fits/{plate_name}-well{well_id:02d}',
                                                             do_plot=do_plot,
                                                             background_model_folder=data_folder + 'multicomp-reactions/2023-03-20-run01/microspectrometer_data/background_model/',
                                                             upper_bounds=upper_bounds,
                                                             cut_from=cut_from,
                                                             cut_to=cut_to,
                                                             ignore_abs_threshold=True,
                                                             ignore_pca_bkg=ignore_bkg_pca,
                                                             use_line=use_line,
                                                             return_errors=get_errors_from_fit)
        if get_errors_from_fit:
            concentrations_here, concentration_errors = concentrations_here
            concentrations_here = np.array(concentrations_here) * dilution_factor
            concentration_errors = np.array(concentration_errors) * dilution_factor
        else:
            concentrations_here = np.array(concentrations_here) * dilution_factor
        print(f'Well {well_id:02d} concentrations: {[calibrant_shortnames[i] + ": " + str(concentrations_here[i]) for i in range(len(calibrant_shortnames))]}')
        # add to dataframe
        df.loc[len(df)] = [well_id] + list(concentrations_here)
        # print the mean ± std

    for column in df.columns[1:]:
        factor = 1000

        if not get_errors_from_fit:
            std_here = df[column].std()
        else:
            std_here = concentration_errors[calibrant_shortnames.index(column)]
        mean_here = df[column].mean()
        # additional relative error due to uncertainty of calibration solutions
        std_here = mean_here * ( (std_here / mean_here) + std_calib )

        print(f'{column}: ({factor * df[column].mean():.2f} ± {factor * std_here:.2f}) mM, n = {len(df)}')
    print(df)

    return df