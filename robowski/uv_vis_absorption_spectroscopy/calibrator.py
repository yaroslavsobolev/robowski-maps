"""
Spectral calibration utilities for UV-Vis absorption spectroscopy.

This module provides functionality for creating calibration references used in spectral
unmixing of UV-Vis absorption spectra from multi-component chemical mixtures. The main
workflow involves processing known standard samples to establish the relationship between
component concentrations and their spectra.

Key Functions:
    construct_calibrant: Main entry point for creating calibration data from standard samples
    _process_single_calibrant: Core processing workflow for individual calibrants
    _calculate_concentration_coefficients: Establishes concentration-to-coefficient relationships
    _apply_reference_stitching: Improves reference spectra using high-concentration calibration samples
    _apply_reference_smoothing: Applies noise reduction to reference spectra

The calibration process typically involves:
1. Loading spectral data from NanoDrop spectrophotometer measurements
2. Background subtraction and spectral preprocessing
3. Creating reference spectra from known standard concentrations
4. Optional spectrum enhancement (stitching, smoothing, external reference integration)
5. Fitting concentration-coefficient relationships for each calibrant
6. Saving calibration data for use in spectral unmixing workflows

The resulting calibration files are used by the spectral unmixing algorithms in
`process_wellplate_spectra.py` to determine component concentrations in unknown mixtures.

Example:
    calibrator.construct_calibrant(
        cut_from=5,
        concentration_column_name='concentration',
        do_plot=False,
        experiment_name='my_experiment/',
        calibration_source_filename='standards_data',
        calibrant_shortnames=['compound_A', 'compound_B'],
        ref_concentrations=[0.001, 0.002],
        max_concentrations=[0.01, 0.02]
    )

Note:
    This module requires external spectral data files and follows the data organization
    structure defined in the robowski package documentation.

    This module is intentionally written in functional style, instead of Object-Oriented style.
    **Reasons:**

    1. **Clear dependencies**: Each function signature explicitly shows what data it needs
    2. **No hidden coupling**: Easy to see exactly what each function depends on
    3. **Simple testing**: Pass individual values rather than constructing objects
    4. **Functional purity**: Functions remain pure with explicit inputs/outputs
    5. **Easy modification**: Can change individual parameters without object restructuring

"""

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
from robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra import create_folder_unless_it_exists
from typing import Optional, Tuple, List
import robowski.uv_vis_absorption_spectroscopy.spectraltools as st

nanodrop_errorbar_folder = (data_folder + 'nanodrop-spectrophotometer-measurements/' +
                            'nanodrop_errorbar_folder_2024-03-16/raw_residuals/')


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
        artefactogenic_upper_limit_of_absorbance=1.5,
        cut_to=None,
        bkg_multiplier=1,
        do_smoothing_at_low_absorbance=0.005,
        savgol_window=31,
        forced_reference_from_agilent_cary_file=None,
        cary_column_name=None,
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

    sp = process_wellplate_spectra.SpectraProcessor(
        folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                        '2022-12-01/interpolator-dataset/')
    sp.nanodrop_lower_cutoff_of_wavelengths = 220 - nanodrop_wavelength_shift
    sp.nanodrop_upper_cutoff_of_wavelengths = 600 - nanodrop_wavelength_shift

    calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'

    create_folder_unless_it_exists(data_folder + experiment_name + 'microspectrometer_data')
    create_folder_unless_it_exists(calibration_folder)
    create_folder_unless_it_exists(calibration_folder + 'references')
    create_folder_unless_it_exists(calibration_folder + 'background')

    all_calibrants_df['nanodrop_col_name'] = all_calibrants_df['nanodrop_col_name'].astype(str)
    nanodrop_df = sp.load_nanodrop_csv_for_one_plate(plate_folder)
    nanodrop_df['wavelength'] = nanodrop_df['wavelength'] + nanodrop_wavelength_shift
    wavelengths = nanodrop_df['wavelength'].to_numpy()

    # load background from the row where concentration is 0
    if custom_bkg_spectrum_npy_file is not None:
        bkg_spectrum = np.load(custom_bkg_spectrum_npy_file)
    else:
        col_name_with_background = all_calibrants_df.loc[all_calibrants_df[concentration_column_name] == 0].iloc[0][
            'nanodrop_col_name']
        absorbances = nanodrop_df[col_name_with_background].to_numpy()
        bkg_spectrum = np.array([wavelengths, absorbances]).T
    bkg_spectrum[:, 1] *= bkg_multiplier

    if do_plot:
        _plot_diagnostic_spectrum(bkg_spectrum[:, 0], bkg_spectrum[:, 1], title="Background spectrum")

    if not do_not_save_data:
        np.save(calibration_folder + f'background//bkg_spectrum.npy', bkg_spectrum)

    for i, calibrant_shortname in enumerate(calibrant_shortnames):
        _process_single_calibrant(
            calibrant_shortname=calibrant_shortname,
            ref_concentration=ref_concentrations[i],
            min_concentration=min_concentrations[i],
            max_concentration=max_concentrations[i],
            do_plot=do_plot,
            all_calibrants_df=all_calibrants_df,
            nanodrop_df=nanodrop_df,
            bkg_spectrum=bkg_spectrum,
            calibration_folder=calibration_folder,
            concentration_column_name=concentration_column_name,
            cut_from=cut_from,
            no_right_edge_subtraction=no_right_edge_subtraction,
            upper_limit_of_absorbance=upper_limit_of_absorbance,
            artefactogenic_upper_limit_of_absorbance=artefactogenic_upper_limit_of_absorbance,
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
        calibrant_shortname,
        ref_concentration,
        min_concentration: float,
        max_concentration: float,
        do_plot: bool,
        all_calibrants_df: pd.DataFrame,
        nanodrop_df: pd.DataFrame,
        bkg_spectrum: np.ndarray,
        calibration_folder: str,
        concentration_column_name: str,
        cut_from: int,
        no_right_edge_subtraction: bool,
        upper_limit_of_absorbance: float,
        artefactogenic_upper_limit_of_absorbance: float,
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
    create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}')

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

    create_folder_unless_it_exists(
        calibration_folder + f'references/{calibrant_shortname}/concentration_fits')

    if forced_reference_from_agilent_cary_file is not None:
        ref_spectrum, reference_interpolator = _load_reference_from_cary_file(
            forced_reference_from_agilent_cary_file=forced_reference_from_agilent_cary_file,
            cary_column_name=cary_column_name,
            wavelengths=nanodrop_df['wavelength'].to_numpy(),
            do_plot=do_plot
        )

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
            artefactogenic_upper_limit_of_absorbance=artefactogenic_upper_limit_of_absorbance,
            no_right_edge_subtraction=no_right_edge_subtraction,
            calibrant_shortname=calibrant_shortname,
            do_plot=do_plot
        )

    if do_smoothing_at_low_absorbance is not None:
        ref_spectrum, reference_interpolator = _apply_reference_smoothing(
            ref_spectrum=ref_spectrum,
            do_smoothing_at_low_absorbance=do_smoothing_at_low_absorbance,
            savgol_window=savgol_window,
            wavelengths=nanodrop_df['wavelength'].to_numpy(),
            do_plot=do_plot
        )

    coeffs = _calculate_concentration_coefficients(
        concentrations=concentrations,
        one_calibrant_df=one_calibrant_df,
        nanodrop_df=nanodrop_df,
        bkg_spectrum=bkg_spectrum,
        reference_interpolator=reference_interpolator,
        concentration_column_name=concentration_column_name,
        cut_from=cut_from,
        cut_to=cut_to,
        upper_limit_of_absorbance=upper_limit_of_absorbance,
        artefactogenic_upper_limit_of_absorbance=artefactogenic_upper_limit_of_absorbance,
        no_right_edge_subtraction=no_right_edge_subtraction,
        do_record_residuals=do_record_residuals,
        dont_save_residuals_below_cut_to=dont_save_residuals_below_cut_to,
        calibrant_shortname=calibrant_shortname,
        calibration_folder=calibration_folder,
        do_plot=do_plot
    )

    _plot_calibration_curve(coeffs, concentrations, do_plot,
                            savefigpath=calibration_folder + f"references/{calibrant_shortname}/concentration-vs-coeff.png")

    if not do_not_save_data:
        _save_calibration_data(
            calibration_folder=calibration_folder,
            calibrant_shortname=calibrant_shortname,
            bkg_spectrum=bkg_spectrum,
            ref_spectrum=ref_spectrum,
            coeffs=coeffs,
            concentrations=concentrations
        )


def _load_reference_from_cary_file(
        forced_reference_from_agilent_cary_file: str,
        cary_column_name: str,
        wavelengths: np.ndarray,
        do_plot: bool
) -> Tuple[np.ndarray, callable]:
    """
    Load and resample reference spectrum from Agilent Cary CSV file.

    This replaces the NanoDrop-derived reference spectrum with one from
    a high-quality Agilent Cary spectrophotometer, resampled to match
    the NanoDrop wavelength grid.

    Returns:
        Tuple of (resampled_ref_spectrum, updated_reference_interpolator)
    """
    wavelengths_cary, ys = st.read_cary_agilent_csv_spectrum(
        forced_reference_from_agilent_cary_file, column_name=cary_column_name)
    ref_spectrum = ys

    # Store for plotting comparison
    ref_spectrum_before_resampling = np.copy(ref_spectrum)

    # Resample Cary spectrum to NanoDrop wavelength grid
    reference_interpolator = interpolate.interp1d(wavelengths_cary, ref_spectrum, fill_value='extrapolate')
    ref_spectrum = reference_interpolator(wavelengths)

    # Create new interpolator for the resampled spectrum
    wavelength_indices = np.arange(ref_spectrum.shape[0])
    reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')

    if do_plot:
        _plot_about_resampling(wavelengths, ref_spectrum, ref_spectrum_before_resampling, wavelengths_cary)

    return ref_spectrum, reference_interpolator


def _apply_reference_smoothing(
        ref_spectrum: np.ndarray,
        do_smoothing_at_low_absorbance: float,
        savgol_window: int,
        wavelengths: np.ndarray,
        do_plot: bool
) -> Tuple[np.ndarray, callable]:
    """
    Apply Savitzky-Golay smoothing to reference spectrum at low absorbance values.

    Uses exponential weighting to blend smoothed and original spectra,
    applying more smoothing at lower absorbance values where noise is more problematic.

    Returns:
        Tuple of (smoothed_ref_spectrum, updated_reference_interpolator)
    """
    original_spectrum = np.copy(ref_spectrum)
    savgol_smoothed_spectrum = savgol_filter(ref_spectrum, window_length=savgol_window, polyorder=4)

    # Create exponential weight between smoothed and original data
    exponential_decay_constant = do_smoothing_at_low_absorbance * np.max(ref_spectrum)
    exponential_weight = np.exp(-1 * ref_spectrum / exponential_decay_constant)

    ref_spectrum = savgol_smoothed_spectrum * exponential_weight + ref_spectrum * (1 - exponential_weight)
    ref_spectrum = ref_spectrum - np.min(ref_spectrum)

    wavelength_indices = np.arange(ref_spectrum.shape[0])
    reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum,
                                                  fill_value='extrapolate')

    if do_plot:
        _plot_smoohting_comparison(wavelengths, original_spectrum, ref_spectrum, savgol_smoothed_spectrum)

    return ref_spectrum, reference_interpolator


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
        artefactogenic_upper_limit_of_absorbance: float,
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

    for concentration in concentrations:
        if concentration < ref_concentration:
            continue

        df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]
        target_spectrum = _load_and_process_spectrum_by_metadata_row(
            df_row_here, nanodrop_df, bkg_spectrum, no_right_edge_subtraction)[:, 1]

        mask = _create_spectrum_mask(
            wavelength_indices, target_spectrum, cut_from, cut_to,
            upper_limit_of_absorbance, artefactogenic_upper_limit_of_absorbance
        )

        popt, pcov = _fit_spectrum_to_reference(
            wavelength_indices, target_spectrum, reference_interpolator,
            mask, initial_scale_guess=concentration / ref_concentration  # or appropriate ratio
        )

        # Plot the fit
        _plot_concentration_fit(df_row_here, do_plot, lambda xs, a, b: a * reference_interpolator(xs) + b,
                                mask, popt, target_spectrum,
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
        artefactogenic_upper_limit_of_absorbance: float,
        no_right_edge_subtraction: bool,
        do_record_residuals: bool,
        dont_save_residuals_below_cut_to: bool,
        calibrant_shortname: str,
        calibration_folder: str,
        do_plot: bool
) -> List[float]:
    """
    Calculate coefficients relating concentration to spectral scaling factors.

    Returns:
        Tuple of (coefficients, coefficient_errors, spectra)
    """
    coeffs = []
    wavelength_indices = np.arange(bkg_spectrum.shape[0])

    for concentration in concentrations:
        if concentration == 0:
            coeffs.append(0)
            continue

        df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]
        target_spectrum = _load_and_process_spectrum_by_metadata_row(df_row_here, nanodrop_df, bkg_spectrum,
                                                                     no_right_edge_subtraction)[:, 1]

        mask = _create_spectrum_mask(
            wavelength_indices, target_spectrum, cut_from, cut_to,
            upper_limit_of_absorbance, artefactogenic_upper_limit_of_absorbance
        )

        # Get artifact mask separately for residuals recording
        artifact_mask = _create_artifact_mask(
            wavelength_indices, target_spectrum, artefactogenic_upper_limit_of_absorbance
        )

        popt, pcov = _fit_spectrum_to_reference(
            wavelength_indices, target_spectrum, reference_interpolator,
            mask, initial_scale_guess=concentration / 0.006
        )

        coeffs.append(popt[0])

        _record_residuals_if_needed(
            target_spectrum=target_spectrum,
            wavelength_indices=wavelength_indices,
            wavelengths=nanodrop_df['wavelength'].to_numpy(),
            reference_interpolator=reference_interpolator,
            popt=popt,
            mask=mask,
            artifact_mask=artifact_mask,
            do_record_residuals=do_record_residuals,
            dont_save_residuals_below_cut_to=dont_save_residuals_below_cut_to,
            df_row=df_row_here
        )

        # Plot concentration fit
        _plot_concentration_fit(df_row_here, do_plot, lambda xs, a, b: a * reference_interpolator(xs) + b,
                                mask, popt, target_spectrum, wavelength_indices,
                                concentration_column_name,
                                savefigpath=calibration_folder + f"references/{calibrant_shortname}/concentration_fits/{df_row_here[concentration_column_name]}_fit.png")

    return coeffs


def _record_residuals_if_needed(
        target_spectrum: np.ndarray,
        wavelength_indices: np.ndarray,
        wavelengths: np.ndarray,
        reference_interpolator: callable,
        popt: np.ndarray,
        mask: np.ndarray,
        artifact_mask: np.ndarray,
        do_record_residuals: bool,
        dont_save_residuals_below_cut_to: bool,
        df_row: pd.Series
) -> None:
    """
    Record residuals from spectrum fitting if requested.

    Saves residuals data for uncertainty analysis, using either the full mask
    or just the artifact mask depending on configuration.
    """
    if not do_record_residuals:
        return

    # Recreate the fitting function
    def func(xs, a, b):
        return a * reference_interpolator(xs) + b

    # Calculate residuals
    residuals = target_spectrum - func(wavelength_indices, *popt)

    # Choose which mask to use for saving
    if dont_save_residuals_below_cut_to:
        save_mask = mask
    else:
        save_mask = artifact_mask

    # Stack data for saving
    residuals_for_saving = np.vstack((
        wavelengths[save_mask],
        target_spectrum[save_mask],
        residuals[save_mask]
    )).T

    # Save residuals
    filename_from_calibration_source = 'dummy_filename'
    np.save(
        f'{nanodrop_errorbar_folder}residuals_{filename_from_calibration_source}__colname{df_row["nanodrop_col_name"]}.npy',
        residuals_for_saving)


def _save_calibration_data(
        calibration_folder: str,
        calibrant_shortname: str,
        bkg_spectrum: np.ndarray,
        ref_spectrum: np.ndarray,
        coeffs: List[float],
        concentrations: List[float]
) -> None:
    """
    Save all calibration data files to disk for later use in spectral unmixing.

    Saves background spectrum, reference spectrum, coefficients, and concentrations
    as numpy arrays in the calibration folder structure.
    """
    np.save(calibration_folder + f'references/{calibrant_shortname}/bkg_spectrum.npy', bkg_spectrum)
    np.save(calibration_folder + f'background//bkg_spectrum.npy', bkg_spectrum)
    np.save(calibration_folder + f'references/{calibrant_shortname}/ref_spectrum.npy', ref_spectrum)
    np.save(calibration_folder + f'references/{calibrant_shortname}/interpolator_coeffs.npy', np.array(coeffs))
    np.save(calibration_folder + f'references/{calibrant_shortname}/interpolator_concentrations.npy', concentrations)


def _create_artifact_mask(
        wavelength_indices: np.ndarray,
        target_spectrum: np.ndarray,
        artefactogenic_upper_limit_of_absorbance: float
) -> np.ndarray:
    """
    Create a mask that excludes wavelength regions contaminated by measurement artifacts.

    Finds the highest wavelength index where absorbance exceeds the artifact threshold
    and excludes all wavelengths up to that point.
    """
    artifact_indices = np.where(target_spectrum > artefactogenic_upper_limit_of_absorbance)[0]
    if len(artifact_indices) == 0:
        largest_artifact_index = -1
    else:
        largest_artifact_index = np.max(artifact_indices)

    artifact_mask = wavelength_indices > largest_artifact_index
    return artifact_mask


def _create_spectrum_mask(
        wavelength_indices: np.ndarray,
        target_spectrum: np.ndarray,
        cut_from: int,
        cut_to: Optional[int],
        upper_limit_of_absorbance: float,
        artefactogenic_upper_limit_of_absorbance: float
) -> np.ndarray:
    """
    Create a boolean mask for spectrum fitting, excluding problematic regions.
    """
    # Basic wavelength and absorbance filtering
    mask = wavelength_indices > cut_from
    if cut_to is not None:
        mask = np.logical_and(mask, wavelength_indices < cut_to)
    mask = np.logical_and(mask, target_spectrum < upper_limit_of_absorbance)

    # Exclude artifact-contaminated regions
    artifact_mask = _create_artifact_mask(
        wavelength_indices, target_spectrum, artefactogenic_upper_limit_of_absorbance
    )
    mask = np.logical_and(mask, artifact_mask)

    return mask


def _fit_spectrum_to_reference(
        wavelength_indices: np.ndarray,
        target_spectrum: np.ndarray,
        reference_interpolator: callable,
        mask: np.ndarray,
        initial_scale_guess: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit target spectrum to reference spectrum using linear scaling.

    Returns:
        Tuple of (fit_parameters, parameter_covariance)
    """

    def func(xs, a, b):
        return a * reference_interpolator(xs) + b

    p0 = (initial_scale_guess, 0)
    bounds = ([-1e-10, -np.inf], [np.inf, np.inf])
    return curve_fit(func, wavelength_indices[mask], target_spectrum[mask],
                           p0=p0, bounds=bounds)


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


def _plot_concentration_fit(df_row_here, do_plot, func, mask, popt, target_spectrum, wavelength_indices,
                            concentration_column_name,
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


def _plot_calibration_curve(coeffs, concentrations, do_plot,
                            savefigpath):
    xs = coeffs
    ys = concentrations
    popt, pcov = curve_fit(lambda x, a: a * x, xs, ys, p0=(1))
    new_xs = np.array([0, 1 * max(xs)])
    new_ys = np.array([0, popt[0] * max(xs)])
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


def take_median_of_nanodrop_spectra(plate_folder, nanodrop_lower_cutoff_of_wavelengths=220,
                                    nanodrop_upper_cutoff_of_wavelengths=600):
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
