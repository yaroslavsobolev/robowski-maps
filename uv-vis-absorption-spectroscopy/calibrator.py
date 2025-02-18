import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import importlib
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
process_wellplate_spectra = importlib.import_module("uv-vis-absorption-spectroscopy.process_wellplate_spectra")
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

nanodrop_errorbar_folder = data_folder + 'nanodrop-spectrophotometer-measurements/nanodrop_errorbar_folder_2024-03-16/raw_residuals/'

st = importlib.import_module('uv-vis-absorption-spectroscopy.spectraltools')

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
                        min_concentrations=None,
                        custom_bkg_spectrum_npy_file=None,
                        no_right_edge_subtraction=False,
                        upper_limit_of_absorbance=1000,
                        do_reference_stitching=False,
                        artefact_generating_upper_limit_of_absorbance=1.5,
                        cut_to=None,
                        bkg_multiplier=1,
                        do_smoothing_at_low_absorbance=0.005,
                        savgol_window=31,
                        forced_reference_from_agilent_cary_file=None,
                        cary_column_name = None,
                        nanodrop_wavelength_shift=0,
                        do_record_residuals=False,
                        do_not_save_data=False,
                        skip_concentrations=tuple([]),
                        dont_save_residuals_below_cut_to=False
):
    if min_concentrations is None:
        min_concentrations = np.zeros(len(calibrant_shortnames))

    plate_folder = f'{data_folder}{experiment_name}{calibration_source_filename}.csv'
    all_calibrants_df = pd.read_csv(f'{data_folder}{experiment_name}{calibration_source_filename}.txt')

    sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                                         '2022-12-01/interpolator-dataset/')
    sp.nanodrop_lower_cutoff_of_wavelengths = 220 - nanodrop_wavelength_shift
    sp.nanodrop_upper_cutoff_of_wavelengths = 600 - nanodrop_wavelength_shift

    calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'

    process_wellplate_spectra.create_folder_unless_it_exists(data_folder + experiment_name + 'microspectrometer_data')
    process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder)
    process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + 'references')
    process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + 'background')

    all_calibrants_df['nanodrop_col_name'] = all_calibrants_df['nanodrop_col_name'].astype(str)

    # load background from different file
    # bkg_spectrum = np.load(calibration_folder + f'background/bkg_spectrum.npy')
    nanodrop_df = sp.load_nanodrop_csv_for_one_plate(plate_folder)
    nanodrop_df['wavelength'] = nanodrop_df['wavelength'] + nanodrop_wavelength_shift
    wavelengths = nanodrop_df['wavelength'].to_numpy()

    # load background from the row where concentration is 0
    if custom_bkg_spectrum_npy_file is not None:
        bkg_spectrum = np.load(custom_bkg_spectrum_npy_file)
    else:
        col_name_with_background = all_calibrants_df.loc[all_calibrants_df[concentration_column_name] == 0].iloc[0]['nanodrop_col_name']
        absorbances = nanodrop_df[col_name_with_background].to_numpy()
        bkg_spectrum = np.array([wavelengths, absorbances]).T
    bkg_spectrum[:, 1] *= bkg_multiplier

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

    if not do_not_save_data:
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

        concentrations = one_calibrant_df[concentration_column_name].to_list()
        concentrations = [x for x in concentrations if (min_concentration <= x <= max_concentration) and (x not in skip_concentrations)]
        concentrations = sorted([0] + concentrations)

        process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}/concentration_fits')

        if forced_reference_from_agilent_cary_file is not None:
            wavelengths_cary, ys = st.read_cary_agilent_csv_spectrum(forced_reference_from_agilent_cary_file, column_name=cary_column_name)
            ref_spectrum = ys
            plt.semilogy(wavelengths_cary, ref_spectrum, label='forced reference')
            plt.title(f'Ref spectrum, forced reference')
            reference_interpolator = interpolate.interp1d(wavelengths_cary, ref_spectrum, fill_value='extrapolate')
            ref_spectrum = reference_interpolator(wavelengths)
            reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')
            plt.semilogy(wavelengths, ref_spectrum, label='forced reference, resampled')
            plt.show()

        if do_reference_stitching:
            linefit_parameters = []
            for concentration in concentrations:
                if concentration < ref_concentration:
                    continue

                df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]
                target_spectrum = load_spectrum_by_df_row(df_row_here)[:, 1]
                mask = wavelength_indices > cut_from
                if cut_to is not None:
                    mask = np.logical_and(mask, wavelength_indices < cut_to)
                mask = np.logical_and(mask, target_spectrum < upper_limit_of_absorbance)

                # find the largest index where target_spectrum is above the value 1.5
                if len(np.where(target_spectrum > artefact_generating_upper_limit_of_absorbance)[0]) == 0:
                    largest_index_above_2 = -1
                else:
                    largest_index_above_2 = np.max(np.where(target_spectrum > artefact_generating_upper_limit_of_absorbance)[0])
                # mask all the indices smaller than largest_index_above_1.5
                mask2 = wavelength_indices > largest_index_above_2
                mask = np.logical_and(mask, mask2)

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
                linefit_parameters.append([popt[0], popt[1]])

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
                if do_plot:
                    plt.show()
                else:
                    plt.clf()

                ref_spectrum[mask] = (target_spectrum[mask] - popt[1])/popt[0]
                ref_spectrum = ref_spectrum - np.min(ref_spectrum)
                reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum,
                                                              fill_value='extrapolate')
                # plot new ref spectrum in semilog scale
                plt.semilogy(wavelengths, ref_spectrum)
                plt.title(f'Ref spectrum, concentration: {concentration}')
                plt.show()

        if do_smoothing_at_low_absorbance is not None:
            savgol_smoothed_signal = savgol_filter(ref_spectrum, window_length=savgol_window, polyorder=4)
            # make exponential weight between the smoothed data and the original
            exponential_decay_constant = do_smoothing_at_low_absorbance * np.max(ref_spectrum)
            exponential_weight = np.exp(-1*ref_spectrum/exponential_decay_constant)
            plt.semilogy(wavelengths, ref_spectrum, label='original spectrum')
            plt.semilogy(wavelengths, savgol_smoothed_signal, label='smoothed spectrum')
            ref_spectrum = savgol_smoothed_signal * exponential_weight + ref_spectrum * (1-exponential_weight)
            ref_spectrum = ref_spectrum - np.min(ref_spectrum)
            reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum,
                                                          fill_value='extrapolate')
            # plot new ref spectrum in semilog scale
            plt.semilogy(wavelengths, ref_spectrum, label='hybrid spectrum')
            plt.title(f'Ref spectrum, savgol_smoothed')
            plt.legend()
            plt.show()

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
            if cut_to is not None:
                mask = np.logical_and(mask, wavelength_indices < cut_to)
            mask = np.logical_and(mask, target_spectrum < upper_limit_of_absorbance)

            # find the largest index where target_spectrum is above the value 1.5
            if len(np.where(target_spectrum > artefact_generating_upper_limit_of_absorbance)[0]) == 0:
                largest_index_above_2 = -1
            else:
                largest_index_above_2 = np.max(np.where(target_spectrum > artefact_generating_upper_limit_of_absorbance)[0])
            # mask all the indices smaller than largest_index_above_1.5
            mask2 = wavelength_indices > largest_index_above_2
            mask = np.logical_and(mask, mask2)

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

            if do_record_residuals:
                residuals = target_spectrum - func(wavelength_indices, *popt)
                # stack wavelengths, target_spectrum, resoduals into a single 3xN array for saving
                if dont_save_residuals_below_cut_to:
                    residuals_for_saving = np.vstack((wavelengths[mask], target_spectrum[mask], residuals[mask])).T
                else:
                    residuals_for_saving = np.vstack((wavelengths[mask2], target_spectrum[mask2], residuals[mask2])).T
                filename_from_calibration_source = calibration_source_filename.split('/')[-1]
                np.save(f'{nanodrop_errorbar_folder}residuals_{filename_from_calibration_source}__colname{df_row_here["nanodrop_col_name"]}.npy', residuals_for_saving)

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

        xs = coeffs
        ys = concentrations
        popt, pcov = curve_fit(lambda x, a: a * x, xs, ys, p0=(1))
        new_xs = np.array([0, 1*max(xs)])
        new_ys = np.array([0, popt[0]*max(xs)])
        ax2.plot(new_xs, new_ys, label='linear fit', color='C1')

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

        if not do_not_save_data:
            # Saving
            np.save(calibration_folder + f'references/{calibrant_shortname}/bkg_spectrum.npy', bkg_spectrum)
            np.save(calibration_folder + f'background//bkg_spectrum.npy', bkg_spectrum)
            np.save(calibration_folder + f'references/{calibrant_shortname}/ref_spectrum.npy', ref_spectrum)
            np.save(calibration_folder + f'references/{calibrant_shortname}/interpolator_coeffs.npy', np.array(coeffs))
            np.save(calibration_folder + f'references/{calibrant_shortname}/interpolator_concentrations.npy',
                    concentrations)

    for i, calibrant_shortname in enumerate(calibrant_shortnames):
        reference_for_one_calibrant(calibrant_shortname, ref_concentration=ref_concentrations[i],
                                    min_concentration=min_concentrations[i],
                                    max_concentration=max_concentrations[i],
                                    do_plot=do_plot)


def construct_degenerate_calibrant():
    pass


def take_median_of_nanodrop_spectra(plate_folder, nanodrop_lower_cutoff_of_wavelengths = 220,
        nanodrop_upper_cutoff_of_wavelengths = 600):

    # repeated measurements
    nanodrop_df = pd.read_csv(plate_folder)

    # rename first column to "wavelength" and make it float type
    nanodrop_df = nanodrop_df.rename(columns={nanodrop_df.columns[0]: "wavelength"})

    # remove rows where wavelength is lower than nanodrop_lower_cutoff_of_wavelengths
    nanodrop_df = nanodrop_df[nanodrop_df["wavelength"] >= nanodrop_lower_cutoff_of_wavelengths]

    # remove rows where wavelength is higher than nanodrop_upper_cutoff_of_wavelengths
    nanodrop_df = nanodrop_df[nanodrop_df["wavelength"] <= nanodrop_upper_cutoff_of_wavelengths]

    nanodrop_df["wavelength"] = nanodrop_df["wavelength"].astype(float)

    # make a new df with same wavelength column
    nanodrop_df_medianned = nanodrop_df[["wavelength"]]

    sample_ids = []
    for column in nanodrop_df.columns[1:]:
        sample_id = column.split('-')[0]
        if sample_id not in sample_ids:
            sample_ids.append(sample_id)
    for sample_id in sample_ids:
        # find all columns which start with the sample_id before underscore, such as 1_0, 1_2, 1_3 if sample id is 1
        columns_to_average = [column for column in nanodrop_df.columns if column.startswith(f'{sample_id}-')]
        # take median of these columns
        stacked_cols = [nanodrop_df[column].to_numpy() for column in columns_to_average]
        ys = np.median(np.array(stacked_cols), axis=0)
        nanodrop_df_medianned[f'{sample_id}'] = ys
        # for col in stacked_cols:
        #     plt.plot(nanodrop_df["wavelength"], col, label='original')
        # plt.plot(nanodrop_df["wavelength"], ys, label='medianned')
        # plt.legend()
        # plt.show()

    # # plot all columns against column 'wavelength' and use labels as columns
    # for column in nanodrop_df.columns[1:]:
    #     plt.plot(nanodrop_df["wavelength"], nanodrop_df[column], label=column)
    #
    # # plot medianned df
    # for column in nanodrop_df_medianned.columns[1:]:
    #     plt.plot(nanodrop_df_medianned["wavelength"], nanodrop_df_medianned[column], label=column, linestyle='--')
    #
    # plt.legend()
    # plt.show()

    # save medianned df to csv
    nanodrop_df_medianned.to_csv(plate_folder[:-4] + '_medianned.csv', index=False)


if __name__ == '__main__':
    pass