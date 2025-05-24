"""
Spectral calibration utilities for UV-Vis absorption spectroscopy.

This module provides functionality for creating calibration references used in spectral
unmixing of UV-Vis absorption spectra from multi-component chemical mixtures. The main
workflow involves processing known standard samples to establish the relationship between
component concentrations and their spectra.

Key Functions:
    perform_calibration: Main entry point for creating calibration data from standard samples
    _process_single_calibrant: Core processing workflow for individual calibrants
    _calculate_concentration_coefficients: Establishes concentration-to-coefficient relationships
    _apply_reference_stitching: Improves reference spectra using high-concentration calibration samples
    _apply_reference_smoothing: Applies noise reduction to reference spectra

The workflow of this module consists of several steps:
1. Loading spectral data from NanoDrop spectrophotometer measurements.
2. Background subtraction and spectral preprocessing.
3. Creating reference spectra from spectra of known standard concentrations, or loading it from an external NumPy file.
4. Optional enhancement of reference spectrum (stitching, smoothing)
5. For each concentration, finding the coefficients by which the reference spectrum must be multiplied in order
        to match the absorbance spectrum of the sample at that concentration. These "coefficients"/"coeffs" as they
        are called throughout the code, are reflecting the magnitude of the spectrum of a specific substance in the
        mix, are used to determine the concentration of the substance in the mixture.
6. Saving calibration data for use in spectral unmixing workflows

The resulting calibration files are used by the spectral unmixing algorithms in
`process_wellplate_spectra.py` to determine component concentrations in unknown mixtures.

Example:
    calibrator.perform_calibration(
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
import os
from typing import Optional, Tuple, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra
import robowski.uv_vis_absorption_spectroscopy.spectraltools as st
from robowski.settings import *
from robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra import create_folder_unless_it_exists

nanodrop_errorbar_folder = (data_folder + 'nanodrop-spectrophotometer-measurements/' +
                            'nanodrop_errorbar_folder_2024-03-16/raw_residuals/')


def perform_calibration(
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
        lower_limit_of_absorbance=0.007,
        npoints_for_right_edge_subtraction=-100
):
    """
    Create calibration data for spectral unmixing from known standard samples.
    --------------------

    This is the main entry point for the calibration workflow. It processes
    multiple calibrants from an input dataset, creating reference spectra and
    conversion from the spectral magnitudes (coefficients) to the concentrations. The results are
    saved to disk for later use in spectral unmixing.

    The calibration process is based on the assumption that absorption spectra of a pure substance
    scale linearly with concentration (Beer-Lambert law). The "coefficient" represents the scaling
    factor by which the reference spectrum must be multiplied to match the spectrum
    measured at a given concentration (known in the case of calibration samples, unknown in the case of
    unknown mixtures). This coefficient is a direct reflection of the magnitude of the substance's
    spectrum at that concentration relative to the reference concentration.
    This coefficient-to-concentration relationship is later used in spectral unmixing to determine
    unknown concentrations from mixed spectra.

    The workflow for each calibrant:
    1. Load data for the calibrant from the input files
    2. Creating reference spectra from spectra of known standard concentrations, or loading it from an external NumPy file.
    3. Optionally enhance the reference using stitching or smoothing
    4. For each concentration, finding the coefficients by which the reference spectrum must be multiplied in order
        to match the absorbance spectrum of the sample at that concentration. These "coefficients"/"coeffs" as they
        are called throughout the code, are reflecting the magnitude of the spectrum of a specific substance in the
        mix, are used to determine the concentration of the substance in the mixture.
    5. Evaluate linear fit of coefficient-to-concentration relationships across all concentrations
    6. Save all calibration data to files for later use, following the specified folder structure

    Structure of input files
    --------------------
    Two input files are required, both with the same base name (calibration_source_filename) but different extensions:

    **Spectral Data File (.csv):**
    Contains the UV-Vis absorption spectra for all samples. The structure is:

    - First column (unnamed): Wavelength values in nanometers (e.g., 190, 191, 192, ...). See note below about the wavelength range.
    - Subsequent columns (0, 1, 2, 3, ...): Absorbance values for each sample (in absorbance units)
    - Each row represents one wavelength point
    - Each numbered column represents one sample measurement

    Example CSV structure (can be found, e.g., in `BPRF\2024-01-17-run01\calibrations` subdirectory of data folder)::

    >     ,0,1,2,3,4,5,6,7,8,9,10,11,...
    >     190,0,0.421,0.705,0.817,0.656,0.868,0.855,0.484,0.248,0.293,0.221,0,...
    >     191,0,1.274,1.556,1.704,1.576,1.744,1.677,1.247,0.898,0.891,0.858,0,...
    >     192,0,2.406,2.57,2.637,2.572,2.618,2.565,2.132,1.784,1.759,1.747,0,...

    **Metadata File (.txt):**
    Contains sample information linking column numbers to chemical identities and concentrations.
    Required columns:
    - nanodrop_col_name: Column number in the CSV file (0, 1, 2, ...)
    - substance: Chemical substance name (must match calibrant_shortnames)
    - concentration: Concentration value in mol/L
    - solvent: Solvent used (for reference)
    Optional columns:
    - dilution_factor: Dilution factor applied. Its value is not currently used in the code.

    For samples containing the pure solvent, the value of `concentration` column should be equal to 0

    Example TXT structure (can be found, e.g., in `BPRF\2024-01-17-run01\calibrations` subdirectory of data folder)::

    >    nanodrop_col_name,substance,dilution_factor,solvent,concentration
    >    0,solvent_Ethanol,1,Ethanol,0
    >    1,methoxybenzaldehyde,1,Ethanol,0.0000025
    >    2,methoxybenzaldehyde,1,Ethanol,0.000005
    >    3,methoxybenzaldehyde,1,Ethanol,0.00001
    >    ...

    The linkage is: CSV column "1" contains the spectrum for the sample described in TXT row 2
    (methoxybenzaldehyde at 0.0000025 M concentration).

    **Integration with Spectral Unmixing:**

    The saved calibration data is later loaded by functions in `process_wellplate_spectra.py`,
    specifically `load_calibration_for_one_calibrant()` and `load_concentration_to_coeff_for_one_calibrant()`,
    to enable determination of component concentrations in unknown mixture spectra.

    Parameters
    ----------
    cut_from : int
        Wavelength index to start analysis from. Earlier wavelengths are ignored. It is counting points from the left
        edge of the spectrum. For example, in a typical nanodrop spectrum, the first point is at 220 nm, so if you set
        cut_from=5, it will start processing from 225 nm onwards.
    concentration_column_name : str
        Name of the column containing concentration values in the metadata file. Normally, its value is 'concentration'.
    do_plot : bool
        Whether to generate and display diagnostic plots during processing. Useful for quality control
        and method development but should be False for automated processing.
    experiment_name : str
        Path to experiment folder relative to data_folder.
    calibration_source_filename : str
        Base filename (without extension) for calibration data files.
        Should exist as both .csv (spectra) and .txt (metadata) files.
    calibrant_shortnames : list of str
        Short names identifying each calibrant in the metadata.
    ref_concentrations : list of float
        List of reference concentration (mol/L) - one for each calibrant in the calibrant_shortnames list.
        In the metadata, must be one and only one row with this concentration value and the respective calibrant shortname.
        This concentration should be chosen to provide good signal-to-noise ratio while avoiding instrumental artifacts
        and the breakdown of the Beer-Lambert law at high absorbance values (if the reference spectrum stays below
        the absorbance 1, it is usually fine).
    max_concentrations : list of float
        Maximum concentration to include for each calibrant (mol/L).
    min_concentrations : list of float, optional
        Minimum concentration to include for each calibrant (mol/L).
        Defaults to zeros if not provided.
    custom_bkg_spectrum_npy_file : str, optional
        Path to custom background spectrum file (in the NumPy binary format). If None, uses zero-concentration
        sample from the provided dataset.
    no_right_edge_subtraction : bool, optional
        Whether to skip right-edge baseline correction. Default is False, which means that the mean value of the
        last 100 points of the spectrum is subtracted from the whole spectrum. The number `100` can be changed by
        setting the `npoints_for_right_edge_subtraction` parameter.
    upper_limit_of_absorbance : float, optional
        Upper limit for absorbance values in fitting. Points above this value are ignored (masked).
        Default is 1000, which in practice means "no limit".
    do_reference_stitching : bool, optional
        Whether to improve reference spectra using higher-concentration samples. Instead of using the reference
        spectrum at the reference concentration, it stitches together the reference spectrum from the spectra of
        higher concentrations, if available. This is useful for improving the quality of the reference spectrum.
        In a way, this synthetically increases the dynamic range of the reference spectrum in terms of absorbance values.
        Important for multispectrum unmixing (when the parent sample is diluted by several different factors and
        the spectra of these diluted samples are used simultaneously for unmixing). See the
        documentation of the `_apply_reference_stitching()` function for details. Default is False.
    artefactogenic_upper_limit_of_absorbance : float, optional
        Threshold of absorbance above which the artefacts are generated by NanoDrop spectrometer. This is not related
        to the breakdown of the Beer-Lambert law, but rather to the artefacts generated by the spectrophotometer
        when almost no light reaches the detector. The algorithm finds the highest
        wavelength where absorbance exceeds this value and masks all lower (bluer) wavelengths, removing
        artifact-contaminated regions caused by insufficient light transmission.
        Default 1.5.
    cut_to : int, optional
        Wavelelengths with index above this integer value are ignored. If None, then all wavelengths with indices
        above cut_from are used. Useful for excluding
        the long sections of almost zero absorbance at the red end of the spectrum, which are not contributing
        anything but noise to the calibration.
    bkg_multiplier : float, optional
        Scaling factor for background spectrum. Default is 1.
    do_smoothing_at_low_absorbance : float, optional
        Absorbance threshold below which to apply Savitzky-Golay smoothing (see method _smooth_the_reference_spectrum()).
        If None, no smoothing is applied. Default 0.005.
    savgol_window : int, optional
        Window length (in points) for Savitzky-Golay filter. Must be an odd integer. Default 31.
    forced_reference_from_agilent_cary_file : str, optional
        Path to CSV file containing a spectrum (measured and saved to CSV by Agilent Cary spectrophotometer)
        to use as reference spectrum instead of the NanoDrop-derived reference spectrum.
    cary_column_name : str, optional
        Column name in the Agilent Cary file to use as the reference spectrum. Only relevant if
        forced_reference_from_agilent_cary_file is provided.
    nanodrop_wavelength_shift : float, optional
        Wavelength offset to apply to all the NanoDrop data (nm). Default is 0.
    do_record_residuals : bool, optional
        Whether to save fitting residuals for uncertainty analysis. Residuals are saved
        to enable later mapping of spectrophotometer errors as a function of wavelength and absorbance
        (see `robowski/uv_vis_absorption_spectroscopy/absorbance_errorbar_model.py`
         and the respective section in the Supplementary Information document of the accompanying research article).
        Default is False.
    do_not_save_data : bool, optional
        Whether to skip saving calibration results to files. Useful for testing. Default is False.
    skip_concentrations : tuple, optional
        Concentrations to exclude from processing. Default is an empty tuple.
    dont_save_residuals_below_cut_to : bool, optional
        Whether to exclude low-wavelength residuals from saved residuals. Only relevant if do_record_residuals is True.
        Default False.
    lower_limit_of_absorbance : float, optional
        Deprecated parameter retained for backward compatibility. Previously used to exclude
        very low absorbance values from fitting.
        Default is 0.007.
    npoints_for_right_edge_subtraction : int, optional
        Negative of the bumber of points from the red end of the spectrum to be used if the 'no_right-edge-subtraction' option is set to
        False. This is the number of points to be used for calculating the mean value of the spectrum at the right edge
        of the spectrum. The mean value is then subtracted from the whole spectrum. This is useful for removing the
        vertical offset of the spectrum.
        Default is (-100), which means "last 100 points".

    Notes
    -----
    This function creates the following directory structure in the calibration folder:

    - background/bkg_spectrum.npy
    - references/{calibrant}/ref_spectrum.npy
    - references/{calibrant}/interpolator_coeffs.npy
    - references/{calibrant}/interpolator_concentrations.npy
    - references/{calibrant}/concentration_fits/{conc}_fit.png (if do_plot=True)

    The saved calibration data is used by the spectral unmixing functions in
    process_wellplate_spectra.py to determine component concentrations in unknown
    mixtures.

    Note that in the specific case of
     NanoDrop spectrophotometer, the output spectrum may contain wavelengths from 190 nm to 220 nm, but sometimes it does not:
     it varies from sample to sample. The guaranteed wavelength range is from 220 nm upwards. For consistent workflow,
     the wavelengths below 220 nm are removed upon loading the data from file.

    Examples
    --------
    >>> perform_calibration(
    ...     cut_from=5,
    ...     concentration_column_name='concentration',
    ...     do_plot=False,
    ...     experiment_name='my_experiment/',
    ...     calibration_source_filename='standards_data',
    ...     calibrant_shortnames=['methoxybenzaldehyde', 'ethyl_acetoacetate'],
    ...     ref_concentrations=[0.001, 0.002],
    ...     max_concentrations=[0.01, 0.02]
    ... )
    """

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
        _calibrate_a_single_calibrant(calibrant_shortname=calibrant_shortname, ref_concentration=ref_concentrations[i],
                                      min_concentration=min_concentrations[i], max_concentration=max_concentrations[i],
                                      skip_concentrations=skip_concentrations, cut_from=cut_from, cut_to=cut_to,
                                      do_reference_stitching=do_reference_stitching,
                                      all_calibrants_df=all_calibrants_df, nanodrop_df=nanodrop_df,
                                      calibration_folder=calibration_folder, cary_column_name=cary_column_name,
                                      concentration_column_name=concentration_column_name, bkg_spectrum=bkg_spectrum,
                                      forced_reference_from_agilent_cary_file=forced_reference_from_agilent_cary_file,
                                      no_right_edge_subtraction=no_right_edge_subtraction,
                                      upper_limit_of_absorbance=upper_limit_of_absorbance,
                                      artefactogenic_upper_limit_of_absorbance=artefactogenic_upper_limit_of_absorbance,
                                      do_smoothing_at_low_absorbance=do_smoothing_at_low_absorbance,
                                      savgol_window=savgol_window, do_plot=do_plot, do_not_save_data=do_not_save_data,
                                      do_record_residuals=do_record_residuals,
                                      dont_save_residuals_below_cut_to=dont_save_residuals_below_cut_to,
                                      lower_limit_of_absorbance=lower_limit_of_absorbance,
                                      npoints_for_right_edge_subtraction=npoints_for_right_edge_subtraction)


def _calibrate_a_single_calibrant(calibrant_shortname, ref_concentration, min_concentration: float,
                                  max_concentration: float, skip_concentrations: Tuple, cut_from: int,
                                  cut_to: Optional[int], do_reference_stitching: bool, all_calibrants_df: pd.DataFrame,
                                  nanodrop_df: pd.DataFrame, calibration_folder: str, cary_column_name: Optional[str],
                                  concentration_column_name: str, bkg_spectrum: np.ndarray,
                                  forced_reference_from_agilent_cary_file: Optional[str],
                                  no_right_edge_subtraction: bool, upper_limit_of_absorbance: float,
                                  artefactogenic_upper_limit_of_absorbance: float,
                                  do_smoothing_at_low_absorbance: Optional[float], savgol_window: int, do_plot: bool,
                                  do_not_save_data: bool, do_record_residuals: bool,
                                  dont_save_residuals_below_cut_to: bool, lower_limit_of_absorbance: float,
                                  npoints_for_right_edge_subtraction=-100) -> None:
    """
    Process complete calibration workflow for a single chemical component.

    This function is called by `perform_calibration()` for each calibrant and handles the entire
    calibration pipeline for one substance. It isolates the calibrant's data from the full dataset,
    creates and optionally enhances the reference spectrum, establishes the quantitative relationship
    between the spectral coefficients (i.e. the magnitude of the spectrum of a specific substance at a given
    concentration relative to the magnitude at the reference concentration) and the concentrations, and saves
    all calibration results to files for later use.

    The calibration workflow for each substance:

    1. **Data isolation**: Extract all concentration data for this specific calibrant from the metadata
    2. **Reference spectrum creation**: Load the spectrum at the reference concentration and apply background subtraction
    3. **Optional spectrum enhancement**:
       - Replace with high-quality Agilent Cary spectrum if provided
       - Apply reference stitching using higher-concentration samples
       - Apply Savitzky-Golay smoothing to reduce noise at low absorbance
    4. **Coefficient calculation**: For each concentration, determine the scaling factor (coefficient) needed to
       match the scaled reference spectrum to the spectrum of same substance measured at this concentration.
    5. **Calibration curve fitting**: Establish the linear relationship between coefficients and concentrations
    6. **Data persistence**: Save all calibration results following the standardized folder structure

    Scientific Background
    ---------------------
    The calibration process is based on the assumption that absorption spectra of a pure substance
    scale linearly with concentration (Beer-Lambert law). The "coefficient" represents the scaling
    factor by which the reference spectrum must be multiplied to match the spectrum
    measured at a given concentration (known in the case of calibration samples, unknown in the case of
    unknown mixtures). This coefficient is a direct reflection of the magnitude of the substance's
    spectrum at that concentration relative to the reference concentration.
    This coefficient-to-concentration relationship is later used in spectral unmixing to determine
    unknown concentrations from mixed spectra.

    Reference stitching (when enabled) improves the reference spectrum by incorporating data from
    higher concentration samples where the signal-to-noise ratio is better. This creates a synthetic
    reference spectrum with improved dynamic range, particularly important for multispectrum unmixing
    where the same "parent" sample is diluted by different factors and the spectra of these diluted samples
    are used simultaneously for unmixing and understanding the concentration of the substance in the "parent" mixture.

    Parameters
    ----------
    calibrant_shortname : str
        Chemical substance identifier that must exactly match the `substance` column values in `all_calibrants_df`.
    ref_concentration : float
        Reference concentration (mol/L) that serves as the basis for the reference spectrum. Must exist
        exactly once in the dataset for this calibrant. This concentration should be chosen to provide good
        signal-to-noise ratio while avoiding instrumental artifacts and the breakdown of the Beer-Lambert law at high
        absorbance values (if the reference spectrum stays below the absorbance 1, it is usually fine).
    min_concentration : float
        Minimum concentration (mol/L) to include in calibration curve fitting. Concentrations below this
        threshold are excluded from the coefficient-concentration relationship.
    max_concentration : float
        Maximum concentration (mol/L) to include in calibration curve fitting. Concentrations above this
        threshold are excluded to avoid non-linear deviation from Beer-Lambert law or instrumental artifacts.
    skip_concentrations : tuple of float
        Specific concentration values to exclude from processing, typically problematic measurements
        identified during quality control.
    cut_from : int
        Wavelength index for analysis start. Wavelengths with smaller indices are masked during fitting.
        This parameter removes noisy or artifact-prone regions at the blue end of the spectrum.
        Foe example, if cut_from=5, the analysis starts from the 6th wavelength point (e.g., from 225 nm if the
        first point of the spectrum is at 220 nm, which is the case for NanoDrop spectrophotometer).
    cut_to : int or None
        Wavelength index for analysis end. If None, analysis continues to the red end of the spectrum.
        Used to exclude regions with poor signal or instrumental artifacts. Useful for excluding
        the long sections of almost zero absorbance at the red end of the spectrum, which are not contributing
        anything but noise to the calibration.
    do_reference_stitching : bool
        Whether to enhance the reference spectrum using higher-concentration samples. When True,
        the algorithm iteratively improves the reference spectrum by incorporating data from samples
        with concentrations higher than ref_concentration, where signal quality is typically better.This is useful for
        improving the quality of the reference spectrum.
        In a way, this synthetically increases the dynamic range of the reference spectrum in terms of absorbance values.
        Important for multispectrum unmixing (when the parent sample is diluted by several different factors and
        the spectra of these diluted samples are used simultaneously for unmixing). See the
        documentation of the `_apply_reference_stitching()` function for details.
    all_calibrants_df : pandas.DataFrame
        Complete metadata dataframe containing sample information for all calibrants. Must contain
        columns: 'substance', concentration_column_name, 'nanodrop_col_name'. See the documentation of
        `perform_calibration()` for details of the expected structure.
    nanodrop_df : pandas.DataFrame
        Spectral data from NanoDrop measurements with wavelengths in first column and sample
        absorbances in numbered columns corresponding to nanodrop_col_name values. See the documentation of
        `perform_calibration()` for details of the expected structure.
    calibration_folder : str
        Root directory path where calibration results will be saved. The function creates
        subdirectories following the pattern: references/{calibrant_shortname}/
        See the documentation of `perform_calibration()` for details of the directory structure.
    cary_column_name : str or None
        Column identifier in the Agilent Cary CSV file when using forced_reference_from_agilent_cary_file.
        The Cary file typically contains multiple columns with sample identifiers.
    concentration_column_name : str
        Name of the concentration column in all_calibrants_df, typically 'concentration'.
    bkg_spectrum : numpy.ndarray
        Background spectrum to subtract from all measurements, shape (n_wavelengths, 2) with
        columns [wavelengths, absorbances]. Usually derived from solvent-only measurements.
    forced_reference_from_agilent_cary_file : str or None
        Path to CSV file containing high-quality reference spectrum from Agilent Cary spectrophotometer.
        When provided, this replaces the NanoDrop-derived reference spectrum.
    no_right_edge_subtraction : bool
        Whether to skip baseline correction using the red end of the spectrum. When False (default),
        the mean of the last npoints_for_right_edge_subtraction points is subtracted from the entire spectrum.
    upper_limit_of_absorbance : float
        Maximum absorbance value for inclusion in fitting. Points above this threshold are masked.
        Prevents fitting in regions where Beer-Lambert law breaks down or instrumental artifacts occur.
    artefactogenic_upper_limit_of_absorbance : float
        Absorbance threshold for detecting spectrophotometer artifacts. The algorithm finds the highest
        wavelength where absorbance exceeds this value and masks all lower wavelengths, removing
        artifact-contaminated regions caused by insufficient light transmission.
    do_smoothing_at_low_absorbance : float or None
        Absorbance threshold below which Savitzky-Golay smoothing is applied to the reference spectrum.
        Smoothing is weighted exponentially, with more smoothing at lower absorbance values where
        noise dominates. If None, no smoothing is performed.
    savgol_window : int
        Window length for Savitzky-Golay filter (must be odd). Larger values provide more smoothing
        but may remove genuine spectral features. Default 31 is suitable for typical NanoDrop spectra.
    do_plot : bool
        Whether to generate and display diagnostic plots during processing. Useful for quality control
        and method development but should be False for automated processing.
    do_not_save_data : bool
        Whether to skip saving calibration results to disk. Used primarily for testing and method
        development where persistent storage is not needed.
    do_record_residuals : bool
        Whether to save fitting residuals for instrumental uncertainty analysis. Residuals are saved
        to enable later mapping of spectrophotometer errors as a function of wavelength and absorbance
        (see `robowski/uv_vis_absorption_spectroscopy/absorbance_errorbar_model.py`).
    dont_save_residuals_below_cut_to : bool
        Whether to exclude low-wavelength residuals from saved data. Only relevant if do_record_residuals is True.
        When True, only residuals above the cut_to wavelength are saved, focusing analysis on the useful spectral range.
    lower_limit_of_absorbance : float
        Deprecated parameter retained for backward compatibility. Previously used to exclude
        very low absorbance values from fitting.
    npoints_for_right_edge_subtraction : int
        Number of points from the red end of the spectrum (negative value) used for baseline correction.
        Default -100 uses the last 100 wavelength points to calculate the baseline offset.

    Raises
    ------
    AssertionError
        If the reference concentration doesn't exist exactly once in the dataset for this calibrant.
        This validation ensures the calibration has a well-defined reference point.

    Notes
    -----
    **Output Files Created:**

    For each calibrant, the following files are saved in calibration_folder/references/{calibrant_shortname}/:

    - `bkg_spectrum.npy`: Background spectrum used for this calibrant
    - `ref_spectrum.npy`: Final reference spectrum (after enhancement if applied)
    - `interpolator_coeffs.npy`: Array of scaling coefficients for each concentration
    - `interpolator_concentrations.npy`: Array of corresponding concentrations
    - `concentration_fits/{concentration}_fit.png`: Diagnostic plots for each fitted concentration (if do_plot=True)

    **Relationship to Parent Function:**

    This function is called once for each calibrant listed in calibrant_shortnames by the parent
    `perform_calibration()` function. The parent function handles dataset loading and iteration,
    while this function focuses on the detailed processing of individual calibrants.

    **Integration with Spectral Unmixing:**

    The saved calibration data is later loaded by functions in `process_wellplate_spectra.py`,
    specifically `load_calibration_for_one_calibrant()` and `load_concentration_to_coeff_for_one_calibrant()`,
    to enable determination of component concentrations in unknown mixture spectra.
    """
    create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}')

    one_calibrant_df, concentrations = _isolate_data_for_a_single_calibrant_from_dataframe(all_calibrants_df,
                                                                                           calibrant_shortname,
                                                                                           concentration_column_name,
                                                                                           ref_concentration,
                                                                                           min_concentration,
                                                                                           max_concentration,
                                                                                           skip_concentrations)

    ref_spectrum = _load_and_process_spectrum_by_row_in_metadata_dataframe(
        one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == ref_concentration].iloc[0], nanodrop_df,
        bkg_spectrum, no_right_edge_subtraction)[:, 1]
    if not no_right_edge_subtraction:
        ref_spectrum -= np.mean(ref_spectrum[npoints_for_right_edge_subtraction:])
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
        ref_spectrum, reference_interpolator = _stitch_reference_spectrum_from_spectra_at_many_concentrations(
            ref_spectrum=ref_spectrum, reference_interpolator=reference_interpolator, concentrations=concentrations,
            ref_concentration=ref_concentration, one_calibrant_df=one_calibrant_df, nanodrop_df=nanodrop_df,
            bkg_spectrum=bkg_spectrum, concentration_column_name=concentration_column_name, cut_from=cut_from,
            cut_to=cut_to, upper_limit_of_absorbance=upper_limit_of_absorbance,
            artefactogenic_upper_limit_of_absorbance=artefactogenic_upper_limit_of_absorbance,
            no_right_edge_subtraction=no_right_edge_subtraction, calibrant_shortname=calibrant_shortname,
            do_plot=do_plot)

    if do_smoothing_at_low_absorbance is not None:
        ref_spectrum, reference_interpolator = _smooth_the_reference_spectrum(ref_spectrum=ref_spectrum,
                                                                              do_smoothing_at_low_absorbance=do_smoothing_at_low_absorbance,
                                                                              savgol_window=savgol_window,
                                                                              wavelengths=nanodrop_df[
                                                                                  'wavelength'].to_numpy(),
                                                                              do_plot=do_plot)

    coeffs = _calculate_coefficients_vs_concentration(concentrations=concentrations, one_calibrant_df=one_calibrant_df,
                                                      nanodrop_df=nanodrop_df, bkg_spectrum=bkg_spectrum,
                                                      reference_interpolator=reference_interpolator,
                                                      concentration_column_name=concentration_column_name,
                                                      cut_from=cut_from, cut_to=cut_to,
                                                      upper_limit_of_absorbance=upper_limit_of_absorbance,
                                                      artefactogenic_upper_limit_of_absorbance=artefactogenic_upper_limit_of_absorbance,
                                                      no_right_edge_subtraction=no_right_edge_subtraction,
                                                      do_record_residuals=do_record_residuals,
                                                      dont_save_residuals_below_cut_to=dont_save_residuals_below_cut_to,
                                                      calibrant_shortname=calibrant_shortname,
                                                      calibration_folder=calibration_folder, do_plot=do_plot)

    _plot_calibration_curve(coeffs, concentrations, do_plot,
                            savefigpath=calibration_folder + f"references/{calibrant_shortname}/concentration-vs-coeff.png")

    if not do_not_save_data:
        _save_calibration_results(calibration_folder=calibration_folder, calibrant_shortname=calibrant_shortname,
                                  bkg_spectrum=bkg_spectrum, ref_spectrum=ref_spectrum, coeffs=coeffs,
                                  concentrations=concentrations)


def _load_reference_from_cary_file(
        forced_reference_from_agilent_cary_file: str,
        cary_column_name: str,
        wavelengths: np.ndarray,
        do_plot: bool
) -> Tuple[np.ndarray, callable]:
    """
    Replace NanoDrop-derived reference spectrum with high-quality Agilent Cary spectrophotometer data.

    This function loads a reference spectrum from an Agilent Cary 5000 spectrophotometer CSV file
    and resamples it to match the NanoDrop wavelength grid. The Cary instrument provides superior
    spectral quality with better wavelength accuracy, higher resolution, and lower noise compared
    to NanoDrop measurements, making it ideal for creating high-quality reference spectra.

    The workflow:

    1. **Load Cary spectrum**: Read the specified column from the Agilent Cary CSV file
    2. **Wavelength alignment**: Interpolate the Cary spectrum onto the NanoDrop wavelength grid
    3. **Interpolator creation**: Generate a new reference interpolator function for the resampled spectrum
    4. **Optional visualization**: Display comparison plots if requested

    Scientific Rationale
    --------------------
    The Agilent Cary 5000 is a research-grade double-beam UV-Vis spectrophotometer with significantly
    better performance characteristics than the NanoDrop 2000c:

    - **Wavelength accuracy**: ±0.08 nm (Cary) vs ±1 nm (NanoDrop)
    - **Spectral bandwidth**: 0.1-5 nm variable (Cary) vs ~1.5 nm fixed (NanoDrop)
    - **Noise level**: <0.0002 Abs (Cary) vs ~0.002-0.005 Abs (NanoDrop)
    - **Stray light**: <0.00003% (Cary) vs <0.02% (NanoDrop)

    Using Cary-derived reference spectra improves calibration accuracy, especially for:
    - Substances with narrow spectral features that may be broadened by NanoDrop
    - Low-concentration calibrations where noise reduction is critical
    - Multispectrum unmixing requiring high-precision reference spectra

    Parameters
    ----------
    forced_reference_from_agilent_cary_file : str
        Path to CSV file containing Agilent Cary spectrophotometer data. The file should follow
        the standard Agilent export format with wavelengths in one column and absorbances in
        named sample columns. The file typically contains a header section followed by the data.
    cary_column_name : str
        Exact column name in the Cary CSV file containing the desired reference spectrum.
        Cary files often contain multiple sample measurements with descriptive column names corresponding
        to the sample names chosen by the user during the measurement,
        like 'sample_1_baseline_corrected' or 'standard_0.001M_rep1'.
    wavelengths : numpy.ndarray
        Target wavelength array from NanoDrop measurements (in nm). The Cary spectrum will be
        interpolated onto this wavelength grid to ensure compatibility with the NanoDrop calibration spectra.
        Typically spans 220-600 nm for NanoDrop data, with steps of 1 nm.
    do_plot : bool
        Whether to display diagnostic plots showing the original Cary spectrum and the resampled
        version. Useful for quality control and verifying that interpolation preserves spectral features.

    Returns
    -------
    ref_spectrum : numpy.ndarray
        Reference spectrum resampled to the NanoDrop wavelength grid. Shape (n_wavelengths,)
        with absorbance values corresponding to the input wavelengths array.
    reference_interpolator : callable
        Interpolation function that takes wavelength indices (not wavelengths) and returns
        absorbance values. Compatible with existing calibration workflow functions.
        Function signature: f(wavelength_indices) -> absorbances

    Notes
    -----
    **File Format Requirements:**

    The Agilent Cary CSV file must follow the standard export format from Cary WinUV software:
    - Header rows containing instrument parameters (automatically skipped)
    - Data section with wavelength column and sample columns
    - Wavelength values in nm, absorbance values in standard units

    **Wavelength Range Considerations:**

    The Cary spectrum should cover at least the full range of the NanoDrop wavelengths (typically
    220-600 nm). If the Cary range is narrower, extrapolation will be used, which may introduce
    artifacts. If the Cary range is broader, only the overlapping region will be used effectively.

    **Interpolation Method:**

    Linear interpolation is used for resampling, which is appropriate for the smooth spectral
    features typically encountered in UV-Vis absorption spectra. The linear interpolation
    avoids any artefacts that could arise from more complex methods like spline fitting, which may
    introduce oscillations in the spectrum. Normally, it preserves
    peak positions and relative intensities while adapting to the NanoDrop wavelength grid.

    **Integration with Calibration Workflow:**

    This function is called by `_calibrate_a_single_calibrant()` when
    `forced_reference_from_agilent_cary_file` is provided. It replaces the NanoDrop-derived
    reference spectrum but maintains full compatibility with subsequent processing steps including
    reference stitching and smoothing.

    Examples
    --------
    Loading a Cary reference spectrum::

        >>> wavelengths = np.arange(220, 601, 1)  # NanoDrop wavelength grid from 220 nm to 600 nm
        >>> ref_spectrum, ref_interpolator = _load_reference_from_cary_file(
        ...     forced_reference_from_agilent_cary_file='path/to/cary_data.csv',
        ...     cary_column_name='methoxybenzaldehyde_0.001M',
        ...     wavelengths=wavelengths,
        ...     do_plot=True
        ... )
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
        _plot_about_resampling(wavelengths_cary, ref_spectrum_before_resampling, wavelengths, ref_spectrum)

    return ref_spectrum, reference_interpolator


def _smooth_the_reference_spectrum(
        ref_spectrum: np.ndarray,
        do_smoothing_at_low_absorbance: float,
        savgol_window: int,
        wavelengths: np.ndarray,
        do_plot: bool
) -> Tuple[np.ndarray, callable]:
    """
    Apply adaptive Savitzky-Golay smoothing to reference spectrum based on absorbance magnitude.

    This function reduces noise in reference spectra by applying Savitzky-Golay filtering with
    exponentially-weighted blending between smoothed and original data. More smoothing is applied
    at lower absorbance values where noise typically dominates the signal, while preserving
    genuine spectral features at higher absorbances where signal-to-noise ratio is better.

    The smoothing workflow:

    1. **Generate smoothed spectrum**: Apply Savitzky-Golay filter to the entire reference spectrum
    2. **Calculate exponential weights**: Determine blending weights based on local absorbance values
    3. **Weighted combination**: Blend original and smoothed spectra according to exponential weights
    4. **Baseline correction**: Subtract the minimum value to ensure non-negative absorbances
    5. **Interpolator update**: Create new interpolation function for the smoothed spectrum

    Scientific Rationale
    --------------------
    At low absorbance values (typically <0.01 Abs), instrumental noise becomes comparable to or
    larger than the genuine absorption signal, particularly for NanoDrop measurements. Traditional
    uniform smoothing would blur genuine spectral features at higher absorbances where signal
    quality is good. This adaptive approach addresses both issues:

    - **Noise reduction**: Aggressive smoothing at low absorbances reduces instrumental noise
    - **Feature preservation**: Minimal smoothing at high absorbances preserves peak shapes and fine structure
    - **Smooth transition**: Exponential weighting provides gradual transition between smoothing levels

    The exponential weighting function is: weight = exp(-absorbance / threshold)
    where threshold = do_smoothing_at_low_absorbance × max(spectrum)

    This approach is particularly beneficial for:

    - Reference spectra derived from low-concentration standards
    - Multispectrum unmixing where reference quality directly impacts accuracy
    - Substances with weak absorption features that are noise-limited

    Parameters
    ----------
    ref_spectrum : numpy.ndarray
        Input reference spectrum to be smoothed. Shape (n_wavelengths,) containing absorbance
        values at each wavelength point. Should be background-subtracted and baseline-corrected
        before smoothing.
    do_smoothing_at_low_absorbance : float
        Absorbance threshold parameter controlling the transition between smoothed and original data.
        Lower values result in more aggressive smoothing. Typical values range from 0.001 to 0.01.
        The actual threshold is calculated as this value times the maximum absorbance in the spectrum.
    savgol_window : int
        Window length for Savitzky-Golay filter in number of wavelength points. Must be an odd
        integer. Larger values provide more smoothing but may remove genuine spectral features.
        Typical values: 11-51 for NanoDrop spectra (1 nm spacing), with 31 being a good default.
    wavelengths : numpy.ndarray
        Wavelength array in nm corresponding to the reference spectrum. Used for diagnostic
        plotting and should match the original NanoDrop wavelength grid (typically 220-600 nm).
    do_plot : bool
        Whether to display diagnostic plots comparing original, smoothed, and final hybrid spectra.
        Essential for method validation and quality control during calibration development.

    Returns
    -------
    ref_spectrum : numpy.ndarray
        Adaptively smoothed reference spectrum with shape (n_wavelengths,). Baseline-corrected
        to have minimum value of zero. Preserves original magnitude scaling.
    reference_interpolator : callable
        Updated interpolation function for the smoothed spectrum. Takes wavelength indices
        (not wavelength values) and returns absorbance values. Compatible with subsequent
        calibration workflow functions.

    Notes
    -----
    **Savitzky-Golay Filter Parameters:**

    The function uses a 4th-order polynomial for Savitzky-Golay filtering (polyorder=4), which
    provides good balance between noise reduction and feature preservation for typical UV-Vis
    absorption spectra. The polynomial order is fixed to avoid overfitting artifacts that can
    occur with higher orders.

    **Exponential Weighting Calculation:**

    The exponential decay constant is calculated as:
    decay_constant = do_smoothing_at_low_absorbance × max(ref_spectrum)

    This adaptive scaling ensures consistent behavior across spectra with different maximum
    absorbances while maintaining the user-specified threshold interpretation.

    **Baseline Correction:**

    The final spectrum is baseline-corrected by subtracting its minimum value to ensure all
    absorbances are non-negative. This prevents artifacts in subsequent concentration calculations
    where negative absorbances would be unphysical.

    **Integration with Calibration Workflow:**

    This function is called by `_calibrate_a_single_calibrant()` when `do_smoothing_at_low_absorbance`
    is not None. It operates after optional reference stitching but before coefficient calculation,
    ensuring the smoothed reference is used for all subsequent fitting operations.

    **Quality Control Considerations:**

    Visual inspection of the diagnostic plots is recommended to verify that:

    - Noise is adequately reduced at low absorbances
    - Genuine spectral features are preserved at high absorbances
    - The transition between regions appears smooth and natural
    - No smoothing artifacts (oscillations, peak shifts) are introduced

    Examples
    --------
    Applying adaptive smoothing to a noisy reference spectrum::

        >>> wavelengths = np.arange(220, 601, 1)  # Standard NanoDrop grid
        >>> ref_spectrum = np.array([...])  # Noisy reference spectrum
        >>> smoothed_spectrum, interpolator = _smooth_the_reference_spectrum(
        ...     ref_spectrum=ref_spectrum,
        ...     do_smoothing_at_low_absorbance=0.005,  # 0.5% of max absorbance threshold
        ...     savgol_window=31,  # 31-point window for 1 nm spacing
        ...     wavelengths=wavelengths,
        ...     do_plot=True  # Show before/after comparison
        ... )
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


def _stitch_reference_spectrum_from_spectra_at_many_concentrations(
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
        target_spectrum = _load_and_process_spectrum_by_row_in_metadata_dataframe(df_row_here, nanodrop_df,
                                                                                  bkg_spectrum,
                                                                                  no_right_edge_subtraction)[:, 1]

        mask = _create_spectrum_mask(
            wavelength_indices, target_spectrum, cut_from, cut_to,
            upper_limit_of_absorbance, artefactogenic_upper_limit_of_absorbance
        )

        popt, pcov = _fit_reference_model_to_target_spectrum(wavelength_indices, target_spectrum,
                                                             reference_interpolator, mask,
                                                             initial_scale_guess=concentration / ref_concentration)

        # Plot the fit
        _plot_concentration_fit(wavelength_indices, target_spectrum,
                                lambda xs, a, b: a * reference_interpolator(xs) + b, popt, mask, df_row_here,
                                concentration_column_name, do_plot)

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


def _calculate_coefficients_vs_concentration(
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
        target_spectrum = _load_and_process_spectrum_by_row_in_metadata_dataframe(df_row_here, nanodrop_df,
                                                                                  bkg_spectrum,
                                                                                  no_right_edge_subtraction)[:, 1]

        mask = _create_spectrum_mask(
            wavelength_indices, target_spectrum, cut_from, cut_to,
            upper_limit_of_absorbance, artefactogenic_upper_limit_of_absorbance
        )

        # Get artifact mask separately for residuals recording
        artifact_mask = _create_artifact_mask(
            wavelength_indices, target_spectrum, artefactogenic_upper_limit_of_absorbance
        )

        popt, pcov = _fit_reference_model_to_target_spectrum(wavelength_indices, target_spectrum,
                                                             reference_interpolator, mask,
                                                             initial_scale_guess=concentration / 0.006)

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
        _plot_concentration_fit(wavelength_indices, target_spectrum,
                                lambda xs, a, b: a * reference_interpolator(xs) + b, popt, mask, df_row_here,
                                concentration_column_name, do_plot,
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
    def model_function(xs, a, b):
        return a * reference_interpolator(xs) + b

    # Calculate residuals
    residuals = target_spectrum - model_function(wavelength_indices, *popt)

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


def _save_calibration_results(
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


def _fit_reference_model_to_target_spectrum(
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

    def model_function(xs, a, b):
        return a * reference_interpolator(xs) + b

    p0 = (initial_scale_guess, 0)
    bounds = ([-1e-10, -np.inf], [np.inf, np.inf])
    return curve_fit(model_function, wavelength_indices[mask], target_spectrum[mask],
                     p0=p0, bounds=bounds)


def _plot_diagnostic_spectrum(wavelengths, spectrum, title="Background spectrum", semilog=False):
    """Plot background spectrum."""
    if semilog:
        plt.semilogy(wavelengths, spectrum)
    else:
        plt.plot(wavelengths, spectrum)
    plt.title(title)
    plt.show()


def _plot_about_resampling(wavelengths_before_resampling, spectrum_before_resampling, resampled_wavelengths,
                           resampled_spectrum):
    plt.semilogy(wavelengths_before_resampling, spectrum_before_resampling, label='Original spectrum from Cary')
    plt.title(f'Resampling of forced reference spectrum from Agilent Cary spectrophotometer')
    plt.semilogy(resampled_wavelengths, resampled_spectrum, label='Resampled spectrum')
    plt.show()


def _plot_concentration_fit(wavelength_indices, target_spectrum, model_function, popt, mask, df_row_here,
                            concentration_column_name, do_plot, savefigpath=None):
    fig1 = plt.figure(1)
    plt.plot(target_spectrum, label='data', color='C0', alpha=0.5)
    mask_illustration = np.ones_like(target_spectrum) * np.max(target_spectrum)
    mask_illustration[mask] = 0
    plt.fill_between(x=wavelength_indices, y1=0, y2=mask_illustration, color='yellow', alpha=0.3,
                     label='ignored (masked) data')
    plt.plot(model_function(wavelength_indices, *popt), color='r', label='fit', alpha=0.5)
    plt.plot(model_function(wavelength_indices, popt[0], 0), color='C1', label='reference', alpha=0.5)
    plt.ylim(-0.03,
             np.max((model_function(wavelength_indices, *popt)[mask])) * 2)
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


def _load_and_process_spectrum_by_row_in_metadata_dataframe(row, nanodrop_df, bkg_spectrum,
                                                            no_right_edge_subtraction=False):
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


def _isolate_data_for_a_single_calibrant_from_dataframe(all_calibrants_df, calibrant_shortname,
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
    assert ref_conc_count == 1, (f"Reference concentration {ref_concentration} for {calibrant_shortname} must "
                                 f"exist exactly once, found {ref_conc_count} times")

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

        ######### PLOTTING; UNCOMMENT IF DEBUGGING #########
        # for col in stacked_cols:
        #     plt.plot(nanodrop_df["wavelength"], col, label='original')
        # plt.plot(nanodrop_df["wavelength"], ys, label='medianned')
        # plt.legend()
        # plt.show()

    ######### PLOTTING; UNCOMMENT IF DEBUGGING #########
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
