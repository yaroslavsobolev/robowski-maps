import logging
from robowski.settings import *
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra
from dataclasses import dataclass, field
from typing import Optional, Tuple, List, Dict


nanodrop_errorbar_folder = data_folder + 'nanodrop-spectrophotometer-measurements/nanodrop_errorbar_folder_2024-03-16/raw_residuals/'
import robowski.uv_vis_absorption_spectroscopy.spectraltools as st


@dataclass
class CalibrantConfig:
    """Configuration for processing a single calibrant."""
    shortname: str
    ref_concentration: float
    min_concentration: float = 0
    max_concentration: float = 1000
    do_plot: bool = True


def construct_calibrant(
                        cut_from,
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
                        dont_save_residuals_below_cut_to=False,
                        lower_limit_of_absorbance=0.007
):
    if min_concentrations is None:
        min_concentrations = np.zeros(len(calibrant_shortnames))

    plate_folder = f'{data_folder}{experiment_name}{calibration_source_filename}.csv'
    all_calibrants_df = pd.read_csv(f'{data_folder}{experiment_name}{calibration_source_filename}.txt')

    sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
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
        _plot_diagnostic_spectrum(bkg_spectrum[:, 0], bkg_spectrum[:, 1], title="Background spectrum")

    if not do_not_save_data:
        np.save(calibration_folder + f'background//bkg_spectrum.npy', bkg_spectrum)

    for i, calibrant_shortname in enumerate(calibrant_shortnames):
        calibrant_config = CalibrantConfig(
            shortname=calibrant_shortname,
            ref_concentration=ref_concentrations[i],
            min_concentration=min_concentrations[i],
            max_concentration=max_concentrations[i],
            do_plot=do_plot
        )

        _process_single_calibrant(
            calibrant_config=calibrant_config,
            all_calibrants_df=all_calibrants_df,
            nanodrop_df=nanodrop_df,
            bkg_spectrum=bkg_spectrum,
            calibration_folder=calibration_folder,
            concentration_column_name=concentration_column_name,
            cut_from=cut_from,
            no_right_edge_subtraction=no_right_edge_subtraction,
            upper_limit_of_absorbance=upper_limit_of_absorbance,
            artefact_generating_upper_limit_of_absorbance=artefact_generating_upper_limit_of_absorbance,
            do_reference_stitching=do_reference_stitching,
            do_smoothing_at_low_absorbance=do_smoothing_at_low_absorbance,
            savgol_window=savgol_window,
            forced_reference_from_agilent_cary_file=forced_reference_from_agilent_cary_file,
            cary_column_name=cary_column_name,
            do_record_residuals=do_record_residuals,
            do_not_save_data=do_not_save_data,
            skip_concentrations=skip_concentrations,
            dont_save_residuals_below_cut_to=dont_save_residuals_below_cut_to,
            lower_limit_of_absorbance=lower_limit_of_absorbance,
            cut_to=cut_to
        )


def _process_single_calibrant(
        calibrant_config: CalibrantConfig,
        all_calibrants_df: pd.DataFrame,
        nanodrop_df: pd.DataFrame,
        bkg_spectrum: np.ndarray,
        calibration_folder: str,
        concentration_column_name: str,
        cut_from: int,
        no_right_edge_subtraction: bool,
        upper_limit_of_absorbance: float,
        artefact_generating_upper_limit_of_absorbance: float,
        do_reference_stitching: bool,
        do_smoothing_at_low_absorbance: Optional[float],
        savgol_window: int,
        forced_reference_from_agilent_cary_file: Optional[str],
        cary_column_name: Optional[str],
        do_record_residuals: bool,
        do_not_save_data: bool,
        skip_concentrations: Tuple,
        dont_save_residuals_below_cut_to: bool,
        lower_limit_of_absorbance: float,
        cut_to: Optional[int]
) -> None:
    """
    Process a single calibrant to create reference spectrum and calibration curve.

    This function validates calibrant data, creates reference spectra, optionally
    applies stitching and smoothing, and establishes concentration-coefficient
    relationships for spectral unmixing.
    """

    calibrant_shortname = calibrant_config.shortname
    ref_concentration = calibrant_config.ref_concentration
    min_concentration = calibrant_config.min_concentration
    max_concentration = calibrant_config.max_concentration
    do_plot = calibrant_config.do_plot

    process_wellplate_spectra.create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}')

    one_calibrant_df, concentrations = _validate_and_filter_calibrant_data(
        all_calibrants_df, calibrant_shortname, concentration_column_name,
        ref_concentration, min_concentration, max_concentration, skip_concentrations)

    ref_spectrum = _load_and_process_spectrum_by_metadata_row(
        one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == ref_concentration].iloc[0],
        nanodrop_df, bkg_spectrum, no_right_edge_subtraction)[:, 1]
    if not no_right_edge_subtraction:
        ref_spectrum -= np.mean(ref_spectrum[-100:])
    if do_plot:
        _plot_diagnostic_spectrum(nanodrop_df['wavelength'].to_numpy(), ref_spectrum,
                                  title=f"Reference spectrum for {calibrant_shortname}")

    wavelength_indices = np.arange(ref_spectrum.shape[0])
    reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')

    process_wellplate_spectra.create_folder_unless_it_exists(
        calibration_folder + f'references/{calibrant_shortname}/concentration_fits')

    if forced_reference_from_agilent_cary_file is not None:
        wavelengths_cary, ys = st.read_cary_agilent_csv_spectrum(forced_reference_from_agilent_cary_file,
                                                                 column_name=cary_column_name)
        ref_spectrum = ys
        # copy it for later plotting
        ref_spectrum_before_resampling = np.copy(ref_spectrum)
        reference_interpolator = interpolate.interp1d(wavelengths_cary, ref_spectrum, fill_value='extrapolate')
        ref_spectrum = reference_interpolator(nanodrop_df['wavelength'].to_numpy())
        reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')
        if do_plot:
            _plot_about_resampling(nanodrop_df['wavelength'].to_numpy(), ref_spectrum, ref_spectrum_before_resampling,
                                   wavelengths_cary)

    if do_reference_stitching:
        ref_spectrum, reference_interpolator = _apply_reference_stitching(
            ref_spectrum=ref_spectrum,
            reference_interpolator=reference_interpolator,
            concentrations=concentrations,
            ref_concentration=ref_concentration,
            one_calibrant_df=one_calibrant_df,
            nanodrop_df=nanodrop_df,
            bkg_spectrum=bkg_spectrum,
            concentration_column_name=concentration_column_name,
            cut_from=cut_from,
            cut_to=cut_to,
            upper_limit_of_absorbance=upper_limit_of_absorbance,
            artefact_generating_upper_limit_of_absorbance=artefact_generating_upper_limit_of_absorbance,
            no_right_edge_subtraction=no_right_edge_subtraction,
            calibrant_shortname=calibrant_shortname,
            do_plot=do_plot
        )

    if do_smoothing_at_low_absorbance is not None:
        original_spectrum = np.copy(ref_spectrum)
        savgol_smoothed_signal = savgol_filter(ref_spectrum, window_length=savgol_window, polyorder=4)
        # make exponential weight between the smoothed data and the original
        exponential_decay_constant = do_smoothing_at_low_absorbance * np.max(ref_spectrum)
        exponential_weight = np.exp(-1 * ref_spectrum / exponential_decay_constant)

        ref_spectrum = savgol_smoothed_signal * exponential_weight + ref_spectrum * (1 - exponential_weight)
        ref_spectrum = ref_spectrum - np.min(ref_spectrum)
        reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum,
                                                      fill_value='extrapolate')
        if do_plot:
            _plot_smoohting_comparison(nanodrop_df['wavelength'].to_numpy(), original_spectrum, ref_spectrum,
                                       savgol_smoothed_signal)

    coeffs, coeff_errs, spectra = _calculate_concentration_coefficients(
        concentrations=concentrations,
        one_calibrant_df=one_calibrant_df,
        nanodrop_df=nanodrop_df,
        bkg_spectrum=bkg_spectrum,
        reference_interpolator=reference_interpolator,
        concentration_column_name=concentration_column_name,
        cut_from=cut_from,
        cut_to=cut_to,
        upper_limit_of_absorbance=upper_limit_of_absorbance,
        artefact_generating_upper_limit_of_absorbance=artefact_generating_upper_limit_of_absorbance,
        no_right_edge_subtraction=no_right_edge_subtraction,
        do_record_residuals=do_record_residuals,
        dont_save_residuals_below_cut_to=dont_save_residuals_below_cut_to,
        calibrant_shortname=calibrant_shortname,
        calibration_folder=calibration_folder,
        do_plot=do_plot
    )

    xs = coeffs
    ys = concentrations
    popt, pcov = curve_fit(lambda x, a: a * x, xs, ys, p0=(1))
    new_xs = np.array([0, 1 * max(xs)])
    new_ys = np.array([0, popt[0] * max(xs)])

    _plot_calibration_curve(coeffs, concentrations, do_plot, new_xs, new_ys,
                            savefigpath=calibration_folder + f"references/{calibrant_shortname}/concentration-vs-coeff.png")

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


def _apply_reference_stitching(
        ref_spectrum: np.ndarray,
        reference_interpolator: callable,
        concentrations: List[float],
        ref_concentration: float,
        one_calibrant_df: pd.DataFrame,
        nanodrop_df: pd.DataFrame,
        bkg_spectrum: np.ndarray,
        concentration_column_name: str,
        cut_from: int,
        cut_to: Optional[int],
        upper_limit_of_absorbance: float,
        artefact_generating_upper_limit_of_absorbance: float,
        no_right_edge_subtraction: bool,
        calibrant_shortname: str,
        do_plot: bool
) -> Tuple[np.ndarray, callable]:
    """
    Apply reference spectrum stitching using higher concentration samples.

    Stitching improves the reference spectrum by using data from higher
    concentration samples where the signal is stronger and more reliable.

    Returns:
        Tuple of (modified_ref_spectrum, updated_reference_interpolator)
    """
    wavelength_indices = np.arange(ref_spectrum.shape[0])
    linefit_parameters = []

    for concentration in concentrations:
        if concentration < ref_concentration:
            continue

        df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]
        target_spectrum = _load_and_process_spectrum_by_metadata_row(
            df_row_here, nanodrop_df, bkg_spectrum, no_right_edge_subtraction)[:, 1]

        # Create mask for fitting
        mask = wavelength_indices > cut_from
        if cut_to is not None:
            mask = np.logical_and(mask, wavelength_indices < cut_to)
        mask = np.logical_and(mask, target_spectrum < upper_limit_of_absorbance)

        # Handle artifacts at very high absorbance
        if len(np.where(target_spectrum > artefact_generating_upper_limit_of_absorbance)[0]) == 0:
            largest_index_above_2 = -1
        else:
            largest_index_above_2 = np.max(np.where(target_spectrum > artefact_generating_upper_limit_of_absorbance)[0])
        mask2 = wavelength_indices > largest_index_above_2
        mask = np.logical_and(mask, mask2)

        # Fit target spectrum to reference
        def func(xs, a, b):
            return a * reference_interpolator(xs) + b

        p0 = (concentration / ref_concentration, 0)
        bounds = ([-1e-10, -np.inf], [np.inf, np.inf])
        popt, pcov = curve_fit(func, wavelength_indices[mask], target_spectrum[mask],
                               p0=p0, bounds=bounds)
        perr = np.sqrt(np.diag(pcov))
        linefit_parameters.append([popt[0], popt[1]])

        # Plot the fit
        _plot_concentration_fit(df_row_here, do_plot, func, mask, popt, target_spectrum,
                                wavelength_indices, concentration_column_name)

        # Update reference spectrum with stitched data
        ref_spectrum[mask] = (target_spectrum[mask] - popt[1]) / popt[0]
        ref_spectrum = ref_spectrum - np.min(ref_spectrum)
        reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum,
                                                      fill_value='extrapolate')

        # Plot the updated reference
        if do_plot:
            _plot_diagnostic_spectrum(
                nanodrop_df['wavelength'].to_numpy(),
                ref_spectrum,
                title=f"Reference spectrum for calibrant {calibrant_shortname}, concentration: {concentration}, after stitching",
                semilog=True
            )

    return ref_spectrum, reference_interpolator


def _calculate_concentration_coefficients(
        concentrations: List[float],
        one_calibrant_df: pd.DataFrame,
        nanodrop_df: pd.DataFrame,
        bkg_spectrum: np.ndarray,
        reference_interpolator: callable,
        concentration_column_name: str,
        cut_from: int,
        cut_to: Optional[int],
        upper_limit_of_absorbance: float,
        artefact_generating_upper_limit_of_absorbance: float,
        no_right_edge_subtraction: bool,
        do_record_residuals: bool,
        dont_save_residuals_below_cut_to: bool,
        calibrant_shortname: str,
        calibration_folder: str,
        do_plot: bool
) -> Tuple[List[float], List[float], List[np.ndarray]]:
    """
    Calculate coefficients relating concentration to spectral scaling factors.

    Returns:
        Tuple of (coefficients, coefficient_errors, spectra)
    """
    coeffs = []
    coeff_errs = []
    spectra = []
    wavelength_indices = np.arange(bkg_spectrum.shape[0])

    for concentration in concentrations:
        if concentration == 0:
            coeffs.append(0)
            coeff_errs.append(0)
            continue

        df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]
        target_spectrum = _load_and_process_spectrum_by_metadata_row(df_row_here, nanodrop_df, bkg_spectrum,
                                                                     no_right_edge_subtraction)[:, 1]
        spectra.append(np.copy(target_spectrum))

        # Create mask for fitting
        mask = wavelength_indices > cut_from
        if cut_to is not None:
            mask = np.logical_and(mask, wavelength_indices < cut_to)
        mask = np.logical_and(mask, target_spectrum < upper_limit_of_absorbance)

        # Handle artifacts at high absorbance
        if len(np.where(target_spectrum > artefact_generating_upper_limit_of_absorbance)[0]) == 0:
            largest_index_above_2 = -1
        else:
            largest_index_above_2 = np.max(np.where(target_spectrum > artefact_generating_upper_limit_of_absorbance)[0])
        mask2 = wavelength_indices > largest_index_above_2
        mask = np.logical_and(mask, mask2)

        # Fit spectrum to reference
        def func(xs, a, b):
            return a * reference_interpolator(xs) + b

        p0 = (concentration / 0.006, 0)  # ref_concentration was typically 0.006
        bounds = ([-1e-10, -np.inf], [np.inf, np.inf])
        popt, pcov = curve_fit(func, wavelength_indices[mask], target_spectrum[mask],
                               p0=p0, bounds=bounds)
        perr = np.sqrt(np.diag(pcov))

        coeffs.append(popt[0])
        coeff_errs.append(perr[0])

        # Handle residuals recording
        if do_record_residuals:
            residuals = target_spectrum - func(wavelength_indices, *popt)
            if dont_save_residuals_below_cut_to:
                residuals_for_saving = np.vstack(
                    (nanodrop_df['wavelength'].to_numpy()[mask], target_spectrum[mask], residuals[mask])).T
            else:
                residuals_for_saving = np.vstack(
                    (nanodrop_df['wavelength'].to_numpy()[mask2], target_spectrum[mask2], residuals[mask2])).T
            filename_from_calibration_source = 'dummy_filename'
            np.save(
                f'{nanodrop_errorbar_folder}residuals_{filename_from_calibration_source}__colname{df_row_here["nanodrop_col_name"]}.npy',
                residuals_for_saving)

        # Plot concentration fit
        _plot_concentration_fit(df_row_here, do_plot, func, mask, popt, target_spectrum, wavelength_indices,
                                concentration_column_name,
                                savefigpath=calibration_folder + f"references/{calibrant_shortname}/concentration_fits/{df_row_here[concentration_column_name]}_fit.png")

    return coeffs, coeff_errs, spectra


def _plot_diagnostic_spectrum(wavelengths, spectrum, title="Background spectrum", semilog=False):
    """Plot background spectrum."""
    if semilog:
        plt.semilogy(wavelengths, spectrum)
    else:
        plt.plot(wavelengths, spectrum)
    plt.title(title)
    plt.show()


def _plot_about_resampling(wavelengths, ref_spectrum, ref_spectrum_before_resampling, wavelengths_cary):
    plt.semilogy(wavelengths_cary, ref_spectrum_before_resampling, label='forced reference')
    plt.title(f'Ref spectrum, forced reference')
    plt.semilogy(wavelengths, ref_spectrum, label='forced reference, resampled')
    plt.show()


def _plot_concentration_fit(df_row_here, do_plot, func, mask, popt, target_spectrum, wavelength_indices, concentration_column_name,
                            savefigpath=None):
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
    if savefigpath is not None:
        fig1.savefig(savefigpath)
        logging.debug(f"Saved figure to {os.path.abspath(savefigpath)}")
    if do_plot:
        plt.show()
    else:
        plt.clf()


def _plot_smoohting_comparison(wavelengths, original_spectrum, ref_spectrum, savgol_smoothed_signal):
    # plot new ref spectrum in semilog scale
    plt.semilogy(wavelengths, original_spectrum, label='original spectrum')
    plt.semilogy(wavelengths, savgol_smoothed_signal, label='smoothed spectrum')
    plt.semilogy(wavelengths, ref_spectrum, label='hybrid spectrum')
    plt.title(f'Ref spectrum, savgol_smoothed')
    plt.legend()
    plt.show()


def _plot_calibration_curve(coeffs, concentrations, do_plot, new_xs, new_ys,
                            savefigpath):
    fig3, axarr = plt.subplots(1, 2, figsize=(10, 5))
    ax1, ax2 = axarr
    ax2.plot(coeffs, concentrations, 'o-')
    ax2.plot(new_xs, new_ys, label='linear fit', color='C1')
    ax2.set_xlabel('Best-fit scale coefficient')
    ax2.set_ylabel('Concentration, mol/liter')
    fig3.savefig(savefigpath, dpi=300)
    if do_plot:
        plt.show()
        plt.loglog(coeffs, concentrations, 'o-')
        plt.xlabel('Best-fit scale coefficient')
        plt.ylabel('Concentration, mol/liter')
        plt.show()
    else:
        plt.clf()


def _load_and_process_spectrum_by_metadata_row(row, nanodrop_df, bkg_spectrum, no_right_edge_subtraction=False):
    """
    Load and process a spectrum for a specific metadata row.

    This function loads the raw spectrum data, subtracts the background,
    and optionally applies right edge correction.

    Parameters
    ----------
    row : pandas.Series
        Metadata row containing 'nanodrop_col_name' column
    nanodrop_df : pandas.DataFrame
        DataFrame with spectral data
    bkg_spectrum : numpy.ndarray
        Background spectrum to subtract
    no_right_edge_subtraction : bool, optional
        Whether to skip right edge baseline correction

    Returns
    -------
    numpy.ndarray
        Processed spectrum as 2D array [wavelengths, absorbances]
    """
    wavelengths = nanodrop_df['wavelength'].to_numpy()
    absorbances = nanodrop_df[row['nanodrop_col_name']].to_numpy()
    spectrum = np.array([wavelengths, absorbances]).T
    spectrum[:, 1] -= bkg_spectrum[:, 1]
    if not no_right_edge_subtraction:
        spectrum[:, 1] -= np.mean(spectrum[-100:, 1])
    return spectrum


def _validate_and_filter_calibrant_data(all_calibrants_df, calibrant_shortname,
                                        concentration_column_name, ref_concentration,
                                        min_concentration, max_concentration,
                                        skip_concentrations):
    """
    Validate and filter calibrant data for a specific calibrant.

    This function extracts data for a specific calibrant, validates that the reference
    concentration exists exactly once, and filters concentrations to the specified range.

    Parameters
    ----------
    all_calibrants_df : pandas.DataFrame
        DataFrame containing all calibrant metadata
    calibrant_shortname : str
        Name of the calibrant to process
    concentration_column_name : str
        Name of the concentration column
    ref_concentration : float
        Reference concentration that must exist exactly once
    min_concentration : float
        Minimum concentration to include
    max_concentration : float
        Maximum concentration to include
    skip_concentrations : tuple
        Concentrations to exclude from processing

    Returns
    -------
    tuple
        (one_calibrant_df, filtered_concentrations)
        - one_calibrant_df: DataFrame with data for this calibrant only
        - filtered_concentrations: Sorted list of concentrations to process

    Raises
    ------
    AssertionError
        If reference concentration doesn't exist exactly once
    """
    # Extract data for this calibrant
    one_calibrant_df = all_calibrants_df[all_calibrants_df['substance'] == calibrant_shortname]

    # Validate that reference concentration exists exactly once
    ref_conc_count = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == ref_concentration].shape[0]
    assert ref_conc_count == 1, f"Reference concentration {ref_concentration} for {calibrant_shortname} must exist exactly once, found {ref_conc_count} times"

    # Filter and sort concentrations
    concentrations = one_calibrant_df[concentration_column_name].to_list()
    concentrations = [x for x in concentrations if
                      (min_concentration <= x <= max_concentration) and (x not in skip_concentrations)]
    concentrations = sorted([0] + concentrations)

    return one_calibrant_df, concentrations


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