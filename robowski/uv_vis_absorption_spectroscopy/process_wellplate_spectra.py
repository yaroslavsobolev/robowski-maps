"""
Utilities for spectral processing and unmixing of UV-Vis absorption spectra from well plates.

This module provides functionality for:
1. Loading spectral data from different sources (CRAIC microspectrometer, NanoDrop spectrophotometer)
2. Processing and correcting raw spectral data
3. Calibrating using reference spectra and known concentrations
4. Unmixing multi-component spectra to determine individual component concentrations
5. Calculating stoichiometric relationships between products and substrates

The main workflow involves:
1. Creating a SpectraProcessor object
2. Loading spectral data from files
3. Creating calibration references using known standards
4. Processing unknown samples to determine component concentrations
"""

from robowski.settings import *
import logging
import pickle

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.signal import savgol_filter
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.ticker as mticker
import statsmodels.api as sm
import time

# matplotlib.use('Agg')
plt.ioff()

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'

def create_folder_unless_it_exists(path):
    """
    Create a folder at the given path if it doesn't already exist.

    Parameters
    ----------
    path : str
        Directory path to create
    """
    if not os.path.exists(path):
        os.makedirs(path)


def load_msp_file(experimental_data_filename, cut=False):
    """
    Load a CRAIC microspectrometer data file (.msp format) containing spectral data.

    Parameters
    ----------
    experimental_data_filename : str
        Path to the .msp file
    cut : tuple or False, optional
        If specified as (min_id, max_id), cuts the spectrum to this wavelength index range.
        If False, no cutting is performed.

    Returns
    -------
    numpy.ndarray
        2D array with shape (n_wavelengths, 2), where the first column contains
        wavelengths (nm) and the second column contains absorbance values
    """
    input_spectrum = np.loadtxt(experimental_data_filename, skiprows=10,
                                delimiter='\t')
    input_spectrum = np.transpose(input_spectrum)
    # input_spectrum_cut = np.flipud(input_spectrum)
    if cut:
        min_id = cut[0]
        max_id = cut[1]
        input_spectrum_cut = input_spectrum[min_id:max_id, :]
    return input_spectrum


def get_spectra_file_list(target_folder, prefix='spectrum_'):
    """
    Get a list of spectrum files (.MSP) extension in the target folder that match the given prefix.
    This excludes the files that have 'rep2' anywhere in their filename, because these files are
    for second repetition of the measurement in CRAIC.

    Parameters
    ----------
    target_folder : str
        Path to folder containing spectrum files
    prefix : str, optional
        Prefix of spectrum files to match, defaults to 'spectrum_'

    Returns
    -------
    list
        List of filenames (excluding any '-3D.msp' files which contain 2D maps, and file containing 'rep2' substring.
    """
    os.chdir(target_folder)
    file_list = glob.glob(f"{prefix}*.msp")
    try:
        file_list.remove(f'{prefix}-3D.msp')  # this file contains the 2D map of absorbance at single fixed wavelength
    except ValueError:
        pass
    return [filename for filename in file_list if 'rep2' not in filename]


def construct_interpolators_for_absorbance_correction(
        nd_names=(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0),
        nd_names_used=(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.5, 3.0, 3.5, 4.0),
        microspec_folder=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/2022-12-01/2-inch-nd-calibrations',
        folder_for_saving_interpolator_datasets=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/2022-12-01/interpolator-dataset/'):
    """
    Construct interpolators for correcting CRAIC microspectrometer absorbance measurements.

    This function generates a wavelength-dependent correction for CRAIC microspectrometer
    absorbance measurements by comparing measurements of identical neutral density (ND)
    filters on both the CRAIC microspectrometer and the reference Agilent Cary 5000
    spectrophotometer. The Agilent instrument is considered the ground truth due to its
    higher accuracy and precision.

    The function:
    1. Loads absorbance spectra for ND filters from both instruments
    2. Applies Savitzky-Golay filtering to smooth the spectra
    3. Aligns the wavelength domains using interpolation
    4. For each wavelength, sorts the filter measurements by absorbance value
    5. Saves the aligned and sorted datasets as numpy arrays for later use

    These saved datasets can be used to create an interpolation function at each wavelength
    that maps from CRAIC absorbance values to the corresponding Agilent values, effectively
    calibrating the CRAIC instrument against the Agilent reference.

    Parameters
    ----------
    nd_names : tuple
        Complete list of neutral density filter optical densities used in Agilent measurements.
        These values correspond to column indices in the Agilent data file.
    nd_names_used : tuple
        Subset of nd_names that are actually used in the correction. This allows excluding
        certain filters that might have measurement issues.
    microspec_folder : str
        Path to the folder containing CRAIC microspectrometer measurements of ND filters.
        Each file should be named according to the filter density (e.g., "0p5.msp" for 0.5 OD).
    folder_for_saving_interpolator_datasets : str
        Path where the processed numpy arrays will be saved for later use in the correction.

    Returns
    -------
    None
        The function doesn't return values directly but saves two numpy arrays:
        - 'craic_data.npy': Array of sorted CRAIC absorbance values for each wavelength
        - 'agilent_data.npy': Array of corresponding Agilent absorbance values

    Notes
    -----
    The saved arrays can be loaded by the `load_dataset_for_absorbance_correction` function
    and used by the `apply_correction` function to correct any CRAIC measurement.

    The filenames in the microspec_folder should follow the pattern "{ND_value}.msp" where
    decimal points are replaced with 'p' (e.g., "0p5.msp" for 0.5 OD).

    The first row in both the craic_data and agilent_data arrays is filled with zeros to
    ensure that a zero absorbance in CRAIC corresponds to a zero absorbance in Agilent.
    """
    microspec_absorbances = dict()

    example_data = load_msp_file(experimental_data_filename=microspec_folder +
                                                            f'/{nd_names_used[0]:.1f}'.replace('.', 'p') + '.msp')
    craic_data = [np.zeros_like(example_data[:, 1])]
    for nd_name in nd_names_used:
        data = load_msp_file(experimental_data_filename=microspec_folder + f'/{nd_name:.1f}'.replace('.', 'p') + '.msp')
        data[:, 1] = savgol_filter(data[:, 1], window_length=31, polyorder=2)
        microspec_absorbances[nd_name] = data[:, 1]
        wavelengths = data[:, 0]
        craic_data.append(data[:, 1])
    craic_data = np.stack(craic_data)

    spectrophotometer_file = repo_data_path + 'uv_vis_absorption_spectroscopy/spectrophotometer-references/' \
                             '2-inch-nd-filters/nd-filters.csv'
    df = pd.read_csv(spectrophotometer_file, skiprows=[0], nrows=451)
    wavelengths_agilent = df.loc[:, 'Wavelength (nm)']
    agilent_data = [np.zeros_like(example_data[:, 1])]
    for nd_name in nd_names_used:
        absorbances_agilent = df.loc[:, f'Abs.{nd_names.index(nd_name) + 2}']
        absorbances_agilent = savgol_filter(absorbances_agilent, window_length=7, polyorder=2)
        agilent_interpolator = interpolate.interp1d(wavelengths_agilent, absorbances_agilent)
        agilent_data.append(agilent_interpolator(wavelengths))
    agilent_data = np.stack(agilent_data)

    for wavelength_id, wavelength in enumerate(wavelengths):
        agilent_here = np.copy(agilent_data[:, wavelength_id])
        craic_here = np.copy(craic_data[:, wavelength_id])
        sorting_ids = craic_here.argsort()
        agilent_here = agilent_here[sorting_ids]
        craic_here = craic_here[sorting_ids]
        agilent_data[:, wavelength_id] = np.copy(agilent_here)
        craic_data[:, wavelength_id] = np.copy(craic_here)

    np.save(folder_for_saving_interpolator_datasets + 'craic_data.npy', craic_data)
    np.save(folder_for_saving_interpolator_datasets + 'agilent_data.npy', agilent_data)


def load_dataset_for_absorbance_correction(target_folder=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                         '2022-12-01/interpolator-dataset/'):
    craic_data = np.load(target_folder + 'craic_data.npy')
    agilent_data = np.load(target_folder + 'agilent_data.npy')
    return [craic_data, agilent_data]


def apply_correction(input_craic_spectrum, absorbance_correction_dataset):
    craic_data, agilent_data = absorbance_correction_dataset
    result = np.zeros_like(input_craic_spectrum)
    for wavelength_id in range(input_craic_spectrum.shape[0]):
        f = interpolate.interp1d(craic_data[:, wavelength_id], agilent_data[:, wavelength_id],
                                 fill_value='extrapolate')
        result[wavelength_id] = f(input_craic_spectrum[wavelength_id])
    return result


def well_id_to_file_id(well_id):
    """
    Convert well ID to file ID based on the scanning pattern of the microspectrometer.

    The microspectrometer scans from the bottom right corner. Each scan line goes up.
    The next scan line is to the left of the previous one. Well ID is counted from the
    top left corner in left-to-right lines, with each next line below the previous one.

    Parameters
    ----------
    well_id : int
        Well ID (0-based, counted from top left in left-to-right, top-to-bottom pattern)

    Returns
    -------
    int
        File ID (1-based, corresponding to the scan sequence numbering)
    """
    # well_id to i,j. Index i is left-to-right. Index j is top-to-bottom. Both start from zero.
    j = well_id // 9
    i = well_id % 9

    # i,j to file_id
    file_id = (8 - i) * 6 + (5 - j) + 1  # plus one because file id is counted from one, not from zero
    return file_id


def load_raw_msp_by_id(plate_folder, well_id, prefix='spectrum_', suffix=''):
    """
    Load a raw microspectrometer data file for a specific well ID.

    Parameters
    ----------
    plate_folder : str
        Path to the folder containing spectrum files
    well_id : int
        Well ID to load
    prefix : str, optional
        Prefix of spectrum files, defaults to 'spectrum_'
    suffix : str, optional
        Suffix to add to the filename, defaults to ''

    Returns
    -------
    numpy.ndarray
        2D array with shape (n_wavelengths, 2), where the first column contains
        wavelengths (nm) and the second column contains absorbance values
    """
    data = load_msp_file(plate_folder + prefix + f'-{well_id_to_file_id(well_id)}{suffix}.msp')
    return data


def diluted_vials_only(list_of_vials_on_plate):
    """
    Returns the indices of diluted vials on the plate.

    Every second row of the plate was not filled with reaction mixtures and is
    later used to hold a diluted reaction mixture.

    Parameters
    ----------
    list_of_vials_on_plate : array-like
        List or array of data for all vials on the plate

    Returns
    -------
    array-like
        Subset of the input containing only the vials that are diluted
    """
    return list_of_vials_on_plate[[i + j for i in [9, 27, 45] for j in range(9)]]


class SpectraProcessor:
    """
    Process and analyze spectral data from microspectrometers and spectrophotometers.

    This class provides methods for loading, processing, and analyzing spectral data
    from various sources, including CRAIC microspectrometers and NanoDrop spectrophotometers.
    It provides functionality for spectral correction, calibration, and concentration determination
    through spectral unmixing.

    The primary workflow is:
    1. Initialize the processor with correction dataset
    2. Load spectral data from files
    3. Apply corrections and calibrations
    4. Calculate concentrations through spectral unmixing

    Attributes
    ----------
    absorbance_correction_dataset : list
        Dataset for correcting absorbance measurements
    spectrum_data_type : str
        Type of spectral data being processed ('craic' by default)
    nanodrop_lower_cutoff_of_wavelengths : float
        Lower wavelength cutoff for NanoDrop data
    nanodrop_upper_cutoff_of_wavelengths : float
        Upper wavelength cutoff for NanoDrop data
    thresh_w_indices : list
        Wavelength indices for threshold interpolation
    thresh_as : list
        Absorbance thresholds for spectral masking
    lower_limit_of_absorbance : float
        Lower limit for absorbance values
    use_instrumental_sigmas : bool
        Whether to use instrument-specific uncertainty values
    sigma_interpolator : callable
        Function for interpolating measurement uncertainties
    uncertainty_of_stoichiometric_overspending_ratio : float
        Uncertainty factor for stoichiometric calculations
    filepath_of_csv_stoichiometry_table : str
        Path to CSV file with stoichiometric data
    substrates : tuple
        Names of substrate chemicals
    df_stoich : pandas.DataFrame
        DataFrame containing stoichiometric information
    """

    def __init__(self, folder_with_correction_dataset, spectrum_data_type='craic',
                 sigma_interpolator_filename=f'{data_folder}nanodrop-spectrophotometer-measurements/'
                                             f'nanodrop_errorbar_folder_2024-03-16/bivariate_spline_interpolator.pkl',
                 filepath_of_csv_stoichiometry_table='BPRF/misc/Hnamesstechiometry3.csv',
                 substrates = ('methoxybenzaldehyde', 'ethyl_acetoacetate', 'ammonium_acetate')):
        """
        Initialize the SpectraProcessor with correction dataset and parameters.

        Parameters
        ----------
        folder_with_correction_dataset : str
            Path to the folder containing absorbance correction data (craic_data.npy and agilent_data.npy)
        spectrum_data_type : str, optional
            Type of spectral data ('craic' by default)
        sigma_interpolator_filename : str, optional
            Path to the pickle file containing the uncertainty interpolator
        filepath_of_csv_stoichiometry_table : str, optional
            Path to the CSV file with stoichiometric information
        substrates : tuple, optional
            Names of substrate chemicals
        """
        self.absorbance_correction_dataset = load_dataset_for_absorbance_correction(
            target_folder=folder_with_correction_dataset)
        self.spectrum_data_type = spectrum_data_type

        self.nanodrop_lower_cutoff_of_wavelengths = 250
        self.nanodrop_upper_cutoff_of_wavelengths = 600

        self.thresh_w_indices = [0, 2000]
        # self.thresh_as = [100.0, 100.0]
        self.thresh_as = [1.0, 1.0]
        self.lower_limit_of_absorbance = 0
        self.use_instrumental_sigmas = False
        with open(sigma_interpolator_filename, 'rb') as f:
            self.sigma_interpolator = pickle.load(f)

        # uncertainty associated with stoichiometric overspending ratio
        self.uncertainty_of_stoichiometric_overspending_ratio = 0.1

        self.filepath_of_csv_stoichiometry_table = filepath_of_csv_stoichiometry_table
        self.substrates = substrates
        self.load_df_stoich()

    def load_df_stoich(self):
        """
        Load the stoichiometry dataframe from the CSV file.

        This method reads the stoichiometric coefficients table from the specified CSV file. Adds
        rows for each substrate. with a 1 in the column corresponding to the substrate
        and 0s in all other substrate columns.
        """
        self.df_stoich = pd.read_csv(data_folder + self.filepath_of_csv_stoichiometry_table)
        # substrates = ['methoxybenzaldehyde', 'ethyl_acetoacetate', 'ammonium_acetate']

        for i, s in enumerate(self.substrates):
            # add row with string s in 'Names', string f'SUB{i}' in Short_names column, number 1 in the column called s and zeros in other columns
            dict_to_add = {'Names': s, 'Short_names': f'SUB{i}', s: 1}
            dict_to_add.update({x: 0 for x in self.substrates if x != s})
            self.df_stoich = self.df_stoich.append(dict_to_add, ignore_index=True)

    def load_nanodrop_csv_for_one_plate(self, plate_folder,
                                        ):
        """
        Load NanoDrop spectrophotometer data from a CSV file into a DataFrame.

        This method loads the CSV file with NanoDrop measurements. The first column is wavelength,
        and the remaining columns are absorbances for each sample (numerated as well/vial).
        It handles various formats including UUIDs attached to column names.

        Parameters
        ----------
        plate_folder : str
            Path to the CSV file with NanoDrop data

        Returns
        -------
        pandas.DataFrame
            DataFrame with NanoDrop measurements, filtered to the wavelength range specified
            by nanodrop_lower_cutoff_of_wavelengths and nanodrop_upper_cutoff_of_wavelengths
        """
        nanodrop_df = pd.read_csv(plate_folder)

        # rename first column to "wavelength" and make it float type
        nanodrop_df = nanodrop_df.rename(columns={nanodrop_df.columns[0]: "wavelength"})

        # print the lowest value of 'wavelength' column
        # print(f"Lowest wavelength in nanodrop file: {nanodrop_df['wavelength'].min()}")

        # remove rows where wavelength is lower than nanodrop_lower_cutoff_of_wavelengths
        nanodrop_df = nanodrop_df[nanodrop_df["wavelength"] >= self.nanodrop_lower_cutoff_of_wavelengths]

        # remove rows where wavelength is higher than nanodrop_upper_cutoff_of_wavelengths
        nanodrop_df = nanodrop_df[nanodrop_df["wavelength"] <= self.nanodrop_upper_cutoff_of_wavelengths]

        nanodrop_df["wavelength"] = nanodrop_df["wavelength"].astype(float)

        # Remove underscore from the column names and everything after it.
        # This is because Yankai has added the UUID of each comdition into the column names -- a good idea, because
        # it allows to cross-validate the relation between spectra and the list of conditions.
        nanodrop_df.columns = nanodrop_df.columns.str.split('_').str[0]

        return nanodrop_df


    def load_single_nanodrop_spectrum(self, plate_folder, well_id):
        """
        Loads the Nanodrop spectrum for a single well.

        The code can also process the nanodrop's CSV files whose nanodrop column names are, for instance, like so:
        `0_4BhYCtsRm6MY7wqCRnBh43,1_iEd9wJZKzJFfF63WizGXsJ,2_bicD5pn9i6yKwsuEwLX59r,3_cHXhtGgQtdSgyBZTx6942M, ...`
        i.e. there is a UUID added after the underscore. In principle, this commit allows for any string to be can
        be added after the underscore, it will be loaded successfully by this code.

        The code retains reverse compatibility to old nanodrop's CSV files that don't contain underscores or UUIDs.

        Parameters
        ----------
        plate_folder: str
            Path to the folder with the Nanodrop measurements.
        well_id: int
            Number of the well.

        Returns
        -------
        nanodrop_spectrum: np.array
            Array with the Nanodrop spectrum.
        """
        nanodrop_df = self.load_nanodrop_csv_for_one_plate(plate_folder=plate_folder)
        wavelengths = nanodrop_df["wavelength"].to_numpy()
        # get the column whose name is equal to well_id
        absorbances = nanodrop_df[str(well_id)].to_numpy()
        nanodrop_spectrum = np.array([wavelengths, absorbances]).T
        return nanodrop_spectrum

    def load_craic_spectrum_by_id(self, plate_folder, well_id, prefix='spectrum_', do_show=False, ignore_second_repetition=False):
        """
        Load and process a CRAIC microspectrometer spectrum by well ID.

        This method loads the spectrum file, applies absorbance correction, and optionally
        applies photo-bleaching correction if a second repetition file exists.

        Parameters
        ----------
        plate_folder : str
            Path to the folder containing spectrum files
        well_id : int
            Well ID to load
        prefix : str, optional
            Prefix of spectrum files, defaults to 'spectrum_'
        do_show : bool, optional
            Whether to display plots of the spectrum, defaults to False
        ignore_second_repetition : bool, optional
            Whether to ignore second repetition files, defaults to False

        Returns
        -------
        numpy.ndarray
            2D array with wavelengths and corrected absorbance values
        """
        spectrum = load_raw_msp_by_id(plate_folder=plate_folder, well_id=well_id, prefix=prefix)

        # if the file of the same name but suffix '_rep2' exists, then load it and apply the 'zero-dose extrapolation'
        # to correct for photobleaching
        if (not ignore_second_repetition) and \
                (os.path.isfile(plate_folder + prefix + f'-{well_id_to_file_id(well_id)}_rep2.msp')):
            try:
                spectrum_rep2 = load_raw_msp_by_id(plate_folder=plate_folder, well_id=well_id, prefix=prefix,
                                                   suffix='_rep2')
                spectrum[:, 1] = spectrum[:, 1] - 0.5*(spectrum_rep2[:, 1] - spectrum[:, 1])
            except FileNotFoundError:
                pass
        if do_show:
            plt.plot(spectrum[:, 0], spectrum[:, 1])
        spectrum[:, 1] = apply_correction(spectrum[:, 1],
                                          absorbance_correction_dataset=self.absorbance_correction_dataset)
        if do_show:
            plt.plot(spectrum[:, 0], spectrum[:, 1], label='corr')
            plt.legend()
            plt.show()
        return spectrum

    def load_msp_by_id(self, plate_folder, well_id, prefix='spectrum_', do_show=False, ignore_second_repetition=False):
        """
        Load a spectrum by well ID, automatically determining the data source type.

        This is a versatile method that can load spectra from both NanoDrop and CRAIC sources.
        It determines the source type from the plate_folder path.

        Parameters
        ----------
        plate_folder : str
            Path to the folder or file containing spectral data
        well_id : int
            Well ID to load
        prefix : str, optional
            Prefix of spectrum files (for CRAIC data), defaults to 'spectrum_'
        do_show : bool, optional
            Whether to display plots of the spectrum, defaults to False
        ignore_second_repetition : bool, optional
            Whether to ignore second repetition files (for CRAIC data), defaults to False

        Returns
        -------
        numpy.ndarray
            2D array with wavelengths and absorbance values
        """

        # if plate folder contains string "nanodrop", then treat the plate_folder as path to nanodrop CSV file
        if 'nanodrop' in plate_folder:
            spectrum = self.load_single_nanodrop_spectrum(plate_folder=plate_folder, well_id=well_id)
        else: # plate_folder is a folder with CRAIC spectra. Use load_craic_spectrum_by_id and pass all args
            spectrum = self.load_craic_spectrum_by_id(plate_folder=plate_folder, well_id=well_id, prefix=prefix,
                                                      do_show=do_show, ignore_second_repetition=ignore_second_repetition)
        return spectrum


    def load_all_spectra(self, plate_folder, prefix='spectrum_-'):
        """
        Load all spectra from a plate.

        Parameters
        ----------
        plate_folder : str
            Path to the folder or file containing spectral data
        prefix : str, optional
            Prefix of spectrum files (for CRAIC data), defaults to 'spectrum_-'

        Returns
        -------
        list
            List of spectra (numpy.ndarray) for all wells
        """
        if 'nanodrop' in plate_folder:
            # load the nanodrop csv file and count the columns
            nanodrop_df = self.load_nanodrop_csv_for_one_plate(plate_folder=plate_folder)
            well_id = 0
            resulting_array = []
            while str(well_id) in nanodrop_df.columns:
                resulting_array.append(self.load_msp_by_id(plate_folder, well_id))
                well_id += 1
            # make a warning if the length of resulting array is higher than 54
            if len(resulting_array) > 54:
                logging.warning(f'Warning: the number of wells is {len(resulting_array)}, '
                             f'which is higher than 54. Check the Nanodrop file.')
            return resulting_array
        else:
            filelist = get_spectra_file_list(plate_folder)
            return [self.load_msp_by_id(plate_folder, well_id) for well_id in range(len(filelist))]


    def show_all_spectra(self, plate_folder, prefix='spectrum_-', specific_well_ids=None):
        """
        Plot all spectra from a plate.

        Parameters
        ----------
        plate_folder : str
            Path to the folder or file containing spectral data
        prefix : str, optional
            Prefix of spectrum files (for CRAIC data), defaults to 'spectrum_-'
        specific_well_ids : list or None, optional
            If specified, only plot spectra for these well IDs
        """
        for well_id, spectrum in enumerate(self.load_all_spectra(plate_folder, prefix=prefix)):
            if specific_well_ids is None or well_id in specific_well_ids:
                plt.plot(spectrum[:, 0], spectrum[:, 1], alpha=0.3, label=f'{well_id}')
            print(f'{well_id}: max {np.max(spectrum[:, 1])}, min {np.min(spectrum[:, 1])}')
        plt.ylabel('Absorbance')
        plt.xlabel('Wavelength, nm')

    def show_all_spectra_for_one_calibrant(self, calibrant_shortname, calibration_sequence_df,
                                           experiment_name, subtract_red_end=True):
        """
        Plot all spectra for a specific calibrant from a calibration sequence.

        Parameters
        ----------
        calibrant_shortname : str
            Short name of the calibrant
        calibration_sequence_df : pandas.DataFrame
            DataFrame with calibration sequence information
        subtract_red_end : bool, optional
            Whether to subtract the median of the last 10 points from the spectrum, defaults to True

        Returns
        -------
        numpy.ndarray
            The last spectrum loaded
        """
        one_calibrant_df = calibration_sequence_df.loc[calibration_sequence_df['shortname'] == calibrant_shortname]
        for index, row in one_calibrant_df.iterrows():
            spectrum = self.load_msp_by_id(
                plate_folder=data_folder + experiment_name + f"microspectrometer_data/calibration/plate-{row['plate_id']:04d}/",
                well_id=row['well_id'])
            if subtract_red_end:
                spectrum[:, 1] -= np.median(spectrum[-10:, 1])
            plt.plot(spectrum[:, 0], spectrum[:, 1], alpha=0.5,
                     label=f"{row['concentration']:.5f} M, well_id:{row['well_id']}")
        plt.legend()
        plt.show()
        return spectrum


    def construct_reference_for_calibrant(self, calibrant_shortname,
                                          calibration_folder, ref_concentration,
                                          do_plot=True, lower_limit_of_absorbance=0.05, do_reference_refinements=True):
        """
        Construct a reference spectrum for a calibrant and save it for later use.

        This method creates a reference spectrum, optionally refines it, and establishes
        the relationship between concentration and absorption coefficient. The results are
        saved to disk for later use in concentration calculations.

        Parameters
        ----------
        calibrant_shortname : str
            Short name of the calibrant
        calibration_folder : str
            Path to the folder containing calibration data
        ref_concentration : float
            Reference concentration for the calibrant
        do_plot : bool, optional
            Whether to display plots during processing, defaults to True
        lower_limit_of_absorbance : float, optional
            Lower limit for absorbance values, defaults to 0.05
        do_reference_refinements : bool, optional
            Whether to perform reference spectrum refinements, defaults to True

        Returns
        -------
        tuple
            (coeff_to_concentration_interpolator, reference_interpolator, bkg_spectrum)
        """
        create_folder_unless_it_exists(calibration_folder + 'references')
        create_folder_unless_it_exists(calibration_folder + 'background')
        create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}')
        calibration_sequence_df = pd.read_csv(calibration_folder + 'calibration_sequence_dataframe.csv')
        one_calibrant_df = calibration_sequence_df.loc[calibration_sequence_df['shortname'] == calibrant_shortname]

        # make sure that only one well for this calibrant has zero concentration. Otherwise it's weird.
        assert one_calibrant_df.loc[one_calibrant_df['concentration'] == 0].shape[0] == 1
        bkg_row = one_calibrant_df.loc[one_calibrant_df['concentration'] == 0].iloc[0]
        bkg_spectrum = self.load_msp_by_id(
            plate_folder=calibration_folder + f"plate-{bkg_row['plate_id']:04d}/",
            well_id=bkg_row['well_id'])

        def load_spectrum_by_df_row(row):
            spectrum = self.load_msp_by_id(
                plate_folder=calibration_folder + f"plate-{row['plate_id']:04d}/",
                well_id=row['well_id'])
            spectrum[:, 1] -= bkg_spectrum[:, 1]
            return spectrum

        # make sure that only one well for this calibrant has concentration equal to ref_concentration
        assert one_calibrant_df.loc[one_calibrant_df['concentration'] == ref_concentration].shape[0] == 1
        ref_spectrum = load_spectrum_by_df_row(
            one_calibrant_df.loc[one_calibrant_df['concentration'] == ref_concentration].iloc[0])[:, 1]
        ref_spectrum -= np.mean(ref_spectrum[-100:])
        if do_plot:
            plt.plot(ref_spectrum)
            plt.title('Ref spectrum')
            plt.show()

        all_spectra = [self.load_msp_by_id(
            plate_folder=calibration_folder + f"plate-{row['plate_id']:04d}/",
            well_id=row['well_id']) for index, row in one_calibrant_df.iterrows()]

        wavelength_indices = np.arange(ref_spectrum.shape[0])
        reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')

        thresh_w_indices = [0, 25, 127, 2000]
        thresh_as = [0.67, 0.75, 1.6, 1.6]
        threshold_interpolator = interpolate.interp1d(thresh_w_indices, thresh_as, fill_value='extrapolate')

        concentrations = sorted(one_calibrant_df['concentration'].to_list())

        def refine_reference(cut_from, row, do_plot=True, use_first_n_points_after_masking=100):
            create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}/refinements')
            target_spectrum = load_spectrum_by_df_row(row)[:, 1]
            mask_containing_entire_tail = np.logical_and(target_spectrum < threshold_interpolator(wavelength_indices),
                                                         wavelength_indices > cut_from)
            first_index_where_data_is_not_ignored = np.argmax(mask_containing_entire_tail)
            mask = np.logical_and(mask_containing_entire_tail,
                                  wavelength_indices < first_index_where_data_is_not_ignored + use_first_n_points_after_masking)

            def func(xs, a, b):
                return a * reference_interpolator(xs) + b

            p0 = (0.5, 0)
            bounds = ([0, -np.inf], [np.inf, np.inf])
            popt, pcov = curve_fit(func, wavelength_indices[mask], target_spectrum[mask],
                                   p0=p0, bounds=bounds)
            # sigma=noise_std*np.ones_like(target_spectrum),
            # absolute_sigma=True)
            perr = np.sqrt(np.diag(pcov))
            slope = popt[0]
            slope_error = perr[0]

            new_ref_spectrum = np.copy(ref_spectrum)
            new_ref_spectrum[mask_containing_entire_tail] = (target_spectrum[mask_containing_entire_tail] - popt[
                1]) / slope
            new_ref_spectrum -= new_ref_spectrum[-1]

            ### PLOTTING
            fig1 = plt.figure(1)
            plt.plot(target_spectrum, label='data', color='C0', alpha=0.5)

            mask_illustration = np.ones_like(target_spectrum) * 4
            mask_illustration[mask] = 0
            plt.fill_between(x=wavelength_indices, y1=0, y2=mask_illustration, color='yellow', alpha=0.3,
                             label='Data is ignored')
            plt.plot(func(wavelength_indices, *popt), color='r', label='fit', alpha=0.5)
            plt.plot(func(wavelength_indices, popt[0], 0), color='C1', label='reference', alpha=0.5)
            plt.ylim(np.min((func(wavelength_indices, *popt)[mask_containing_entire_tail])),
                     np.max((func(wavelength_indices, *popt)[mask_containing_entire_tail])) * 2)
            plt.title(
                f"conc {row['concentration']}, well {row['well_id']}, plate {row['plate_id']:04d}")
            plt.legend()
            fig1.savefig(
                calibration_folder + f"references/{calibrant_shortname}/refinements/{row['concentration']}_refinement_fit.png",
                dpi=300)
            if do_plot:
                plt.show()
            else:
                plt.clf()


            fig2 = plt.figure(2)
            plt.title(
                f"Refined reference, well {row['concentration']}, well {row['well_id']}, plate {row['plate_id']:04d} was used")
            plt.plot(new_ref_spectrum - np.min(new_ref_spectrum), color='black', alpha=0.5, label='new reference')
            plt.plot(ref_spectrum - np.min(new_ref_spectrum), color='C0', alpha=0.5, label='old reference')
            plt.yscale('log')
            plt.legend()
            fig2.savefig(
                calibration_folder + f"references/{calibrant_shortname}/refinements/{row['concentration']}_refined_result.png",
                dpi=300)
            if do_plot:
                plt.show()
            else:
                plt.clf()

            return new_ref_spectrum, interpolate.interp1d(wavelength_indices, new_ref_spectrum,
                                                          fill_value='extrapolate')

        if do_reference_refinements:
            for concentration in concentrations:
                if concentration <= ref_concentration:
                    # Right tail of absorption band is better only in spectra having higher concentrations than the reference
                    continue
                df_row_here = one_calibrant_df.loc[one_calibrant_df['concentration'] == concentration].iloc[0]
                ref_spectrum, reference_interpolator = refine_reference(cut_from=250, row=df_row_here, do_plot=False)

        create_folder_unless_it_exists(calibration_folder + f'references/{calibrant_shortname}/concentration_fits')
        # cut_from = 115
        cut_from=200
        coeffs = []
        coeff_errs = []
        for concentration in concentrations:
            if concentration == 0:
                coeffs.append(0)
                coeff_errs.append(0)
                continue

            df_row_here = one_calibrant_df.loc[one_calibrant_df['concentration'] == concentration].iloc[0]
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
                f"conc {df_row_here['concentration']}, well {df_row_here['well_id']}, plate {df_row_here['plate_id']:04d}")
            plt.legend()
            fig1.savefig(
                calibration_folder + f"references/{calibrant_shortname}/concentration_fits/{df_row_here['concentration']}_fit.png")
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

        return coeff_to_concentration_interpolator, reference_interpolator, bkg_spectrum

    def load_calibration_for_one_calibrant(self, calibrant_shortname, calibration_folder, use_line_fit=False,
                                           do_savgol_filtering=False, ignore_acetic_dependence=True):
        """
        Load calibration data for a calibrant from saved files.

        Parameters
        ----------
        calibrant_shortname : str
            Short name of the calibrant
        calibration_folder : str
            Path to the folder containing calibration data
        use_line_fit : bool, optional
            Whether to use a linear fit for the calibration curve, defaults to False
        do_savgol_filtering : bool, optional
            Whether to apply Savitzky-Golay filtering to the reference spectrum, defaults to False
        ignore_acetic_dependence : bool, optional
            Whether to ignore acetic acid influence on the spectrum, defaults to True

        Returns
        -------
        tuple
            (coeff_to_concentration_interpolator, reference_interpolator, bkg_spectrum)
        """
        bkg_spectrum = np.load(calibration_folder + f'references/{calibrant_shortname}/bkg_spectrum.npy')
        # assert len(bkg_spectrum) == 381

        coeffs = np.load(calibration_folder + f'references/{calibrant_shortname}/interpolator_coeffs.npy')
        concentrations = np.load(
            calibration_folder + f'references/{calibrant_shortname}/interpolator_concentrations.npy')
        if not use_line_fit:
            coeff_to_concentration_interpolator = interpolate.interp1d(coeffs, concentrations,
                                                                       fill_value='extrapolate')
        else:
            # make a line fit with zero intercept
            xs = coeffs
            ys = concentrations
            popt, pcov = curve_fit(lambda x, a: a * x, xs, ys, p0=(1))
            new_xs = np.array([0, 1])
            new_ys = np.array([0, popt[0]])
            coeff_to_concentration_interpolator = interpolate.interp1d(new_xs, new_ys,
                                                                          fill_value='extrapolate')

        # if there is a file called 'acetic_acid_influence.pkl' then open it
        if os.path.isfile(calibration_folder + f'references/{calibrant_shortname}/acetic_acid_influence.pkl') and (not ignore_acetic_dependence):
            with open(calibration_folder + f'references/{calibrant_shortname}/acetic_acid_influence.pkl', 'rb') as f:
                array_of_wavelengths, acetic_acid_concentrations, spectral_2d_grid = pickle.load(f)
            reference_interpolator = RegularGridInterpolator((array_of_wavelengths, acetic_acid_concentrations), spectral_2d_grid,
                                             bounds_error=False, fill_value=None)
        else:
            ref_spectrum = np.load(calibration_folder + f'references/{calibrant_shortname}/ref_spectrum.npy')
            # assert len(ref_spectrum) == 381
            wavelength_indices = np.arange(ref_spectrum.shape[0])
            if do_savgol_filtering:
                plt.plot(ref_spectrum)
                ref_spectrum = savgol_filter(ref_spectrum, 7, 3)
                plt.plot(ref_spectrum)
                plt.show()
            reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')
        return coeff_to_concentration_interpolator, reference_interpolator, bkg_spectrum

    def load_concentration_to_coeff_for_one_calibrant(self, calibrant_shortname, calibration_folder, use_line_fit=False):
        """
        Load calibration data and create a concentration-to-coefficient interpolator.

        This is the inverse of the coefficient-to-concentration interpolator.

        Parameters
        ----------
        calibrant_shortname : str
            Short name of the calibrant
        calibration_folder : str
            Path to the folder containing calibration data
        use_line_fit : bool, optional
            Whether to use a linear fit for the calibration curve, defaults to False

        Returns
        -------
        tuple
            (concentration_to_coeff_interpolator, reference_interpolator, bkg_spectrum)
        """
        bkg_spectrum = np.load(calibration_folder + f'references/{calibrant_shortname}/bkg_spectrum.npy')

        coeffs = np.load(calibration_folder + f'references/{calibrant_shortname}/interpolator_coeffs.npy')
        concentrations = np.load(
            calibration_folder + f'references/{calibrant_shortname}/interpolator_concentrations.npy')
        # assert that there is zero concentration
        assert 0 in concentrations, f"Zero concentration is not in the list of concentrations for {calibrant_shortname}"
        assert 0 in coeffs # assert that there is zero coefficient

        if not use_line_fit:
            concentration_to_coeff_interpolator = interpolate.interp1d(concentrations, coeffs,
                                                                       fill_value='extrapolate')
        else:
            # make a line fit with zero intercept
            xs = coeffs
            ys = concentrations
            popt, pcov = curve_fit(lambda x, a: a * x, xs, ys, p0=(1))
            new_xs = np.array([0, 1])
            new_ys = np.array([0, popt[0]])
            concentration_to_coeff_interpolator = interpolate.interp1d(new_ys, new_xs,
                                                                          fill_value='extrapolate')

        ref_spectrum = np.load(calibration_folder + f'references/{calibrant_shortname}/ref_spectrum.npy')
        wavelength_indices = np.arange(ref_spectrum.shape[0])
        reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')
        return concentration_to_coeff_interpolator, reference_interpolator, bkg_spectrum


    def spectrum_to_concentration(self, target_spectrum_input, calibration_folder, calibrant_shortnames,
                                  background_model_folder,
                                  lower_limit_of_absorbance=-0.2, fig_filename='temp', do_plot=False, #lower_limit_of_absorbance=0.02
                                  upper_bounds=[np.inf, np.inf], use_line=False, cut_from = 200, ignore_abs_threshold=False,
                                  cut_to = False, ignore_pca_bkg=False, return_errors=False): #upper_bounds=[np.inf, np.inf]
        """
        Calculate the concentrations of multiple components from a single spectrum.

        This method performs spectral unmixing to determine the concentrations of multiple
        components in a mixture based on their reference spectra.

        Parameters
        ----------
        target_spectrum_input : numpy.ndarray
            1D array of absorbance values
        calibration_folder : str
            Path to the folder containing calibration data
        calibrant_shortnames : list
            List of short names of the calibrants to include in the unmixing
        background_model_folder : str
            Path to the folder containing background model data
        lower_limit_of_absorbance : float, optional
            Lower limit for absorbance values, defaults to -0.2
        fig_filename : str, optional
            Base filename for saving figures, defaults to 'temp'
        do_plot : bool, optional
            Whether to display plots during processing, defaults to False
        upper_bounds : list, optional
            Upper bounds for calibrant concentrations, defaults to [np.inf, np.inf]
        use_line : bool, optional
            Whether to use a linear term in the model, defaults to False
        cut_from : int, optional
            Wavelength index to start the analysis from, defaults to 200
        ignore_abs_threshold : bool, optional
            Whether to ignore absorbance thresholds, defaults to False
        cut_to : int or False, optional
            Wavelength index to end the analysis at, defaults to False
        ignore_pca_bkg : bool, optional
            Whether to ignore PCA background components, defaults to False
        return_errors : bool, optional
            Whether to return concentration errors, defaults to False

        Returns
        -------
        list or tuple
            List of calculated concentrations, or tuple of (concentrations, errors) if return_errors=True
        """
        calibrants = []
        for calibrant_shortname in calibrant_shortnames:
            dict_here = dict()
            dict_here['coeff_to_concentration_interpolator'], dict_here['reference_interpolator'], dict_here[
                'bkg_spectrum'] = \
                self.load_calibration_for_one_calibrant(calibrant_shortname, calibration_folder)
            calibrants.append(dict_here.copy())

        bkg_spectrum = calibrants[0]['bkg_spectrum']
        wavelengths = bkg_spectrum[:, 0]
        target_spectrum = target_spectrum_input - bkg_spectrum[:, 1]
        wavelength_indices = np.arange(calibrants[0]['bkg_spectrum'].shape[0])

        thresh_w_indices = [0, 25, 127, 2000]
        thresh_as = [0.67, 0.75, 1.6, 1.6]
        threshold_interpolator = interpolate.interp1d(thresh_w_indices, thresh_as, fill_value='extrapolate')

        if not ignore_abs_threshold:
            mask = np.logical_and(target_spectrum < threshold_interpolator(wavelength_indices),
                                  wavelength_indices > cut_from)
        else:
            mask = wavelength_indices > cut_from

        if cut_to:
            mask = np.logical_and(mask, wavelength_indices <= cut_to)

        mask = np.logical_and(mask,
                              target_spectrum > np.min(target_spectrum) + lower_limit_of_absorbance)

        if not ignore_pca_bkg:
            xxx = np.load(background_model_folder + f'component_0.npy')
            yyy = np.load(background_model_folder + f'component_1.npy')
            background_interpolators = [interpolate.interp1d(wavelength_indices,
                                                          np.load(background_model_folder + f'component_{i}.npy'),
                                                          fill_value='extrapolate')
                                     for i in range(2)]
        else:
            background_interpolators = [interpolate.interp1d(wavelength_indices,
                                                             np.ones_like(wavelength_indices),
                                                             fill_value='extrapolate')
                                        for i in range(2)]

        if len(wavelength_indices[mask]) == 0:
            print('There is no data that is within mask. Returning zeros.')
            return [0 for i in range(4)]

        ## old implementation
        # if len(calibrant_shortnames) == 2:
        #     def func(xs, a, b, c, d, e, f):
        #         return a * calibrants[0]['reference_interpolator'](xs) + b * calibrants[1]['reference_interpolator'](xs) + c \
        #                + d*xs + e * background_interpolators[0](xs) + f * background_interpolators[1](xs)
        # elif len(calibrant_shortnames) == 3:
        #     def func(xs, a1, a2, a3, c, d, e, f):
        #         return a1 * calibrants[0]['reference_interpolator'](xs) + \
        #                a2 * calibrants[1]['reference_interpolator'](xs) + \
        #                a3 * calibrants[2]['reference_interpolator'](xs)\
        #                + c + d * xs + e * background_interpolators[0](xs) + f * background_interpolators[1](xs)
        # else:
        #     raise NotImplementedError

        ## New implementation
        def func(*args):
            xs = args[0]
            c,d,e,f = args[-4:]
            calibrant_coefficients = args[1:-4]
            return sum([calibrant_coefficients[i] * calibrants[i]['reference_interpolator'](xs) for i in range(len(calibrant_coefficients))]) \
                      + c + d * xs + e * background_interpolators[0](xs) + f * background_interpolators[1](xs)

        # p0 = tuple(0.5 if upper_bounds[0] is np.inf else upper_bounds[0],
        #       0.5 if upper_bounds[1] is np.inf else upper_bounds[1],
        #       0,
        #       0,
        #       0,
        #       0)
        p0 = tuple([0.5 if upper_bound is np.inf else upper_bound for upper_bound in upper_bounds] + [0] * 4)
        if use_line:
            linebounds = [-np.inf, np.inf]
        else:
            linebounds = [-1e-15, 1e-15]

        if ignore_pca_bkg:
            bkg_comp_limit = 1e-12
        else:
            bkg_comp_limit = np.inf
        bounds = ([-1e-20] * len(calibrant_shortnames) + [-np.inf, linebounds[0], -1*bkg_comp_limit, -1*bkg_comp_limit],
                  upper_bounds + [np.inf, linebounds[1], bkg_comp_limit, bkg_comp_limit])
        popt, pcov = curve_fit(func, wavelength_indices[mask], target_spectrum[mask],
                               p0=p0, bounds=bounds)
        perr = np.sqrt(np.diag(pcov))  # errors of the fitted coefficients

        concentrations_here = [calibrants[calibrant_index]['coeff_to_concentration_interpolator'](fitted_coeff)
                               for calibrant_index, fitted_coeff in enumerate(popt[:-4])]


        # #plot covariance matrix
        # plt.figure(figsize=(5, 10))
        # pcov_to_plot = pcov[:len(calibrant_shortnames), :len(calibrant_shortnames)]
        # plt.imshow(pcov_to_plot, vmin=-1*max(np.abs(pcov_to_plot).flatten()), vmax=max(np.abs(pcov_to_plot).flatten()),
        #            cmap='RdBu_r')
        # # make tick labels from calibrant_shortnames
        # plt.yticks(range(len(calibrant_shortnames)), calibrant_shortnames)
        # plt.xticks(range(len(calibrant_shortnames)), calibrant_shortnames, rotation=90)
        # plt.colorbar(orientation='vertical', fraction=0.046)
        # plt.tight_layout()
        # plt.show()


        fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
        ax = ax1
        ax.plot(wavelengths, target_spectrum_input, label='Raw data', color='grey', alpha=0.2)
        ax.plot(wavelengths, target_spectrum, label='Data minus bkg.', color='black', alpha=0.5)
        mask_illustration = np.ones_like(target_spectrum) * np.max(target_spectrum)
        mask_illustration[mask] = 0
        ax.fill_between(x=wavelengths, y1=0, y2=mask_illustration, color='yellow', alpha=0.3,
                         label='Masked data')
        ax.plot(wavelengths, func(wavelength_indices, *popt), color='r', label='Fit', alpha=0.5)
        for calibrant_index in range(len(calibrant_shortnames)):
            cpopt = [x if i == calibrant_index else 0 for i, x in enumerate(popt)]
            ax.plot(wavelengths, func(wavelength_indices, *cpopt), label=calibrant_shortnames[calibrant_index], alpha=0.5)
        # make a list where only the third from the end item is the same as in popt, while the other ones are zero
        if use_line:
            cpopt = [x if i == len(popt) - 3 else 0 for i, x in enumerate(popt)]
            ax.plot(wavelengths, func(wavelength_indices, *cpopt), label='Line', alpha=0.5)
        if not ignore_pca_bkg:
            cpopt = [x if i == len(popt) - 2 else 0 for i, x in enumerate(popt)]
            ax.plot(wavelengths, func(wavelength_indices, *cpopt), label='Bkg. PC1', alpha=0.5)
            cpopt = [x if i == len(popt) - 1 else 0 for i, x in enumerate(popt)]
            ax.plot(wavelengths, func(wavelength_indices, *cpopt), label='Bkg. PC2', alpha=0.5)
        # plt.ylim(-0.3,
        #          np.max((func(wavelength_indices, *popt)[mask])) * 3)
        title_str = f'Concentrations:\n'
        for i in range(len(concentrations_here)):
            title_str += f'{np.array(concentrations_here)[i]:.6f} M ({calibrant_shortnames[i]})\n '
        fig1.suptitle(title_str[:-2])
        ax.set_ylabel('Absorbance')
        ax.legend()
        # Residuals subplot
        ax = ax2
        ax.plot(wavelengths[mask], target_spectrum[mask] - func(wavelength_indices[mask], *popt), color='black', alpha=0.5,
                label='residuals')
        ax.legend()
        ax.set_xlabel('Wavelength, nm')
        ax.set_ylabel('Absorbance')
        fig1.savefig(f"{fig_filename}.png")

        if do_plot:
            plt.show()
        else:
            plt.close(fig1)
            plt.close('all')
            plt.clf()

        if return_errors:
            # convert coefficient errors into concentration errors
            upper_confidence_limit = [calibrants[calibrant_index]['coeff_to_concentration_interpolator'](fitted_coeff + perr[calibrant_index])
                               for calibrant_index, fitted_coeff in enumerate(popt[:-4])]
            concentration_errors = [upper_confidence_limit[i] - concentrations_here[i] for i in range(len(concentrations_here))]
            return concentrations_here, concentration_errors

        return concentrations_here


    def concentrations_for_one_plate(self, experiment_folder, plate_folder,
                                      calibration_folder, calibrant_shortnames, calibrant_upper_bounds,
                                     background_model_folder, do_plot=False, return_all_substances=False,
                                     cut_from = 200, cut_to=False, ignore_abs_threshold=False, ignore_pca_bkg=False):
        """
        Calculate concentrations for all wells in a plate.

        Parameters
        ----------
        experiment_folder : str
            Path to the experiment folder
        plate_folder : str
            Path to the folder or file containing spectral data
        calibration_folder : str
            Path to the folder containing calibration data
        calibrant_shortnames : list
            List of short names of the calibrants to include in the unmixing
        calibrant_upper_bounds : list
            Upper bounds for calibrant concentrations
        background_model_folder : str
            Path to the folder containing background model data
        do_plot : bool, optional
            Whether to display plots during processing, defaults to False
        return_all_substances : bool, optional
            Whether to return concentrations for all substances, defaults to False
        cut_from : int, optional
            Wavelength index to start the analysis from, defaults to 200
        cut_to : int or False, optional
            Wavelength index to end the analysis at, defaults to False
        ignore_abs_threshold : bool, optional
            Whether to ignore absorbance thresholds, defaults to False
        ignore_pca_bkg : bool, optional
            Whether to ignore PCA background components, defaults to False

        Returns
        -------
        numpy.ndarray
            Array of calculated concentrations for all wells
        """
        plate_name = plate_folder.split('/')[-1]
        create_folder_unless_it_exists(experiment_folder + 'results')
        create_folder_unless_it_exists(experiment_folder + f'results/uv-vis-fits')
        # input_compositions = pd.read_csv(path_to_input_compositions_csv)
        concentrations = []
        # for index, row in input_compositions.iterrows():
        #     plate_id = index // 54
        #     well_id = index % 54
        #     print(f'{plate_id}-{well_id}')
        if 'nanodrop' in plate_folder:
            # load the nanodrop csv file and count the columns
            nanodrop_df = self.load_nanodrop_csv_for_one_plate(plate_folder=plate_folder)
            well_id = 0
            range_of_wells = []
            while str(well_id) in nanodrop_df.columns:
                range_of_wells.append(well_id)
                well_id += 1
            # make a warning if the length of resulting array is higher than 54
            if len(range_of_wells) > 54:
                logging.warning(f'Warning: the number of wells is {len(range_of_wells)}, '
                             f'which is higher than 54. Check the Nanodrop file.')
        else:
            range_of_wells = range(54)

        for well_id in range_of_wells:
            spectrum = self.load_msp_by_id(
                plate_folder=plate_folder,
                well_id=well_id)[:, 1]
            concentrations_here = self.spectrum_to_concentration(target_spectrum_input=spectrum,
                                                                 calibration_folder=calibration_folder,
                                                                 calibrant_shortnames=calibrant_shortnames,
                                                                 fig_filename=experiment_folder + f'results/uv-vis-fits/{plate_name}-well{well_id:02d}.png',
                                                                 do_plot=do_plot,
                                                                 background_model_folder=background_model_folder,
                                                                 upper_bounds=calibrant_upper_bounds, cut_from=cut_from,
                                                                 cut_to=cut_to,
                                                                 ignore_abs_threshold=ignore_abs_threshold,
                                                                 ignore_pca_bkg=ignore_pca_bkg)
            if return_all_substances:
                concentrations.append(concentrations_here)
            else:
                concentrations.append(concentrations_here[0])
        # input_compositions[calibrant_shortnames[0]] = concentrations
        # input_compositions.to_csv(
        #     data_folder + experiment_name + f'results/timepoint{timepoint_id:03d}-reaction_results.csv', index=False)
        return np.array(concentrations)


    def multispectrum_concentrations_for_one_plate(self, experiment_folder, plate_folder_1,
                                                   plate_folder_2, dilution_factors,
                                      calibration_folder, calibrant_shortnames, calibrant_upper_bounds,
                                     background_model_folder, do_plot=False, return_all_substances=False,
                                     cut_from = 200, cut_to=False, ignore_abs_threshold=False, ignore_pca_bkg=False,
                                                   return_report=False, upper_limit_of_absorbance=1000,
                                                   list_of_starting_concentration_dicts=None,
                                                   obey_stoichiometric_inequalities=True):
        """
        Calculate concentrations from multiple spectra of the same wells with different dilutions.

        This method is used when a plate has been measured twice with different dilution factors
        to capture both high and low concentration components accurately.

        Parameters
        ----------
        experiment_folder : str
            Path to the experiment folder
        plate_folder_1 : str
            Path to the first plate data
        plate_folder_2 : str
            Path to the second plate data
        dilution_factors : list
            List of dilution factors for each plate
        calibration_folder : str
            Path to the folder containing calibration data
        calibrant_shortnames : list
            List of short names of the calibrants to include in the unmixing
        calibrant_upper_bounds : list
            Upper bounds for calibrant concentrations
        background_model_folder : str
            Path to the folder containing background model data
        do_plot : bool, optional
            Whether to display plots during processing, defaults to False
        return_all_substances : bool, optional
            Whether to return concentrations for all substances, defaults to False
        cut_from : int, optional
            Wavelength index to start the analysis from, defaults to 200
        cut_to : int or False, optional
            Wavelength index to end the analysis at, defaults to False
        ignore_abs_threshold : bool, optional
            Whether to ignore absorbance thresholds, defaults to False
        ignore_pca_bkg : bool, optional
            Whether to ignore PCA background components, defaults to False
        return_report : bool, optional
            Whether to return a detailed report of the unmixing, defaults to False
        upper_limit_of_absorbance : float, optional
            Upper limit for absorbance values, defaults to 1000
        list_of_starting_concentration_dicts : list or None, optional
            List of dictionaries with starting concentrations for each well
        obey_stoichiometric_inequalities : bool, optional
            Whether to enforce stoichiometric constraints, defaults to True

        Returns
        -------
        numpy.ndarray or tuple
            Array of calculated concentrations, or tuple of (concentrations, reports) if return_report=True
        """
        if list_of_starting_concentration_dicts is None:
            list_of_starting_concentration_dicts = [None] * len(calibrant_shortnames)
        plate_name = plate_folder_1.split('/')[-1]
        create_folder_unless_it_exists(experiment_folder + 'results')
        create_folder_unless_it_exists(experiment_folder + f'results/uv-vis-fits')
        # input_compositions = pd.read_csv(path_to_input_compositions_csv)
        concentrations = []
        reports = []
        # for index, row in input_compositions.iterrows():
        #     plate_id = index // 54
        #     well_id = index % 54
        #     print(f'{plate_id}-{well_id}')
        if 'nanodrop' in plate_folder_1:
            # load the nanodrop csv file and count the columns
            nanodrop_df = self.load_nanodrop_csv_for_one_plate(plate_folder=plate_folder_1)
            well_id = 0
            range_of_wells = []
            while str(well_id) in nanodrop_df.columns:
                range_of_wells.append(well_id)
                well_id += 1
            # make a warning if the length of resulting array is higher than 54
            if len(range_of_wells) > 54:
                logging.warning(f'Warning: the number of wells is {len(range_of_wells)}, '
                             f'which is higher than 54. Check the Nanodrop file.')
        else:
            range_of_wells = range(54)

        for well_id in range_of_wells:
            spectrum_1 = self.load_msp_by_id(
                plate_folder=plate_folder_1,
                well_id=well_id)[:, 1]
            spectrum_2 = self.load_msp_by_id(
                plate_folder=plate_folder_2,
                well_id=well_id)[:, 1]
            print(f'>>>>>>>>>>>>>>>>>> Well {well_id}')
            concentrations_here = self.multispectrum_to_concentration(
                                            target_spectrum_inputs=[spectrum_1, spectrum_2],
                                            dilution_factors=dilution_factors,
                                            calibration_folder=calibration_folder,
                                            calibrant_shortnames=calibrant_shortnames,
                                            fig_filename=experiment_folder + f'results/uv-vis-fits/{plate_name}-well{well_id:02d}.png',
                                            do_plot=do_plot,
                                            background_model_folder=background_model_folder,
                                            upper_bounds=calibrant_upper_bounds, cut_from=cut_from,
                                            cut_to=cut_to,
                                            ignore_abs_threshold=ignore_abs_threshold,
                                            ignore_pca_bkg=ignore_pca_bkg,
                                            upper_limit_of_absorbance=upper_limit_of_absorbance,
                                            return_report=return_report,
                                            starting_concentration_dict=list_of_starting_concentration_dicts[well_id],
                                            obey_stoichiometric_inequalities=obey_stoichiometric_inequalities)
            if return_report:
                concentrations_here, report = concentrations_here
                reports.append(report)

            if return_all_substances:
                concentrations.append(concentrations_here)
            else:
                concentrations.append(concentrations_here[0])
        # input_compositions[calibrant_shortnames[0]] = concentrations
        # input_compositions.to_csv(
        #     data_folder + experiment_name + f'results/timepoint{timepoint_id:03d}-reaction_results.csv', index=False)
        if not return_report:
            return np.array(concentrations)
        else:
            return np.array(concentrations), reports

    def concentrations_for_all_plates(self, timepoint_id, experiment_folder,
                                      calibration_folder,
                                      calibrant_shortnames,
                                      path_to_input_compositions_csv,
                                      calibrant_upper_bounds, do_plot=False):
        """
        Calculate concentrations for all wells in all plates for a specific timepoint.

        Parameters
        ----------
        timepoint_id : int
            ID of the timepoint to process
        experiment_folder : str
            Path to the experiment folder
        calibration_folder : str
            Path to the folder containing calibration data
        calibrant_shortnames : list
            List of short names of the calibrants to include in the unmixing
        path_to_input_compositions_csv : str
            Path to the CSV file with input compositions
        calibrant_upper_bounds : list
            Upper bounds for calibrant concentrations
        do_plot : bool, optional
            Whether to display plots during processing, defaults to False

        Returns
        -------
        pandas.DataFrame
            DataFrame with calculated concentrations for all plates and wells
        """
        create_folder_unless_it_exists(experiment_folder + 'results')
        create_folder_unless_it_exists(experiment_folder + f'results/uv-vis-fits')
        input_compositions = pd.read_csv(path_to_input_compositions_csv)
        concentrations = []
        for index, row in input_compositions.iterrows():
            plate_id = index // 54
            well_id = index % 54
            print(f'{plate_id}-{well_id}')
            spectrum = sp.load_msp_by_id(
                plate_folder=experiment_folder + f"microspectrometer_data/timepoint_{timepoint_id:03d}/plate-{plate_id:02d}/",
                well_id=well_id)[:, 1]
            concentrations_here = self.spectrum_to_concentration(target_spectrum_input=spectrum,
                                                                 calibration_folder=calibration_folder,
                                                                 calibrant_shortnames=calibrant_shortnames,
                                                                 fig_filename=experiment_folder + f'results/uv-vis-fits/plate{plate_id:04d}-well{well_id:02d}.png',
                                                                 do_plot=do_plot,
                                                                 upper_bounds=calibrant_upper_bounds)
            concentrations.append(concentrations_here[0])
        input_compositions[calibrant_shortnames[0]] = concentrations
        input_compositions.to_csv(
            experiment_folder + f'results/timepoint{timepoint_id:03d}-reaction_results.csv', index=False)
        return input_compositions


    def get_absorbance_at_single_wavelength_for_one_plate(self, plate_folder, wavelength=None, ref_wavelengths=None,
                                                          wavelength_id = 100, ref_wavelength_id=[500]):
        """
        Get the absorbance at a specific wavelength for all wells in a plate.

        This method calculates the difference between absorbance at the target wavelength
        and the mean absorbance at reference wavelengths.

        Parameters
        ----------
        plate_folder : str
            Path to the folder or file containing spectral data
        wavelength : float or None, optional
            Target wavelength (nm), defaults to None
        ref_wavelengths : list or None, optional
            List of reference wavelengths (nm), defaults to None
        wavelength_id : int, optional
            Target wavelength index, used if wavelength is None, defaults to 100
        ref_wavelength_id : list, optional
            List of reference wavelength indices, used if ref_wavelengths is None, defaults to [500]

        Returns
        -------
        numpy.ndarray
            Array of absorbance values for all wells
        """
        if not (wavelength is None):
            wavelengths = self.load_msp_by_id(plate_folder=plate_folder, well_id=0)[:, 0]
            wavelength_id = np.absolute(wavelengths - wavelength).argmin()
            ref_wavelength_id = [np.absolute(wavelengths - ref_wavelength).argmin() for ref_wavelength in ref_wavelengths]

        wavelengths = self.load_msp_by_id(plate_folder=plate_folder, well_id=0)[:, 0]
        print(f'Wavelength for wavelength_id is: {wavelengths[wavelength_id]}')
        for ref_wav in ref_wavelength_id:
            print(f'Reference wavelength for ref_wavelength_id is: {wavelengths[ref_wav]}')

        concentrations = []
        for well_id in range(54):
            spectrum = self.load_msp_by_id(plate_folder=plate_folder, well_id=well_id)[:, 1]
            concentrations.append(spectrum[wavelength_id] - np.mean(np.array([spectrum[ref_wav] for ref_wav in ref_wavelength_id])))
        return np.array(concentrations)


    def mask_the_spectrum(self, wavelength_indices, target_spectrum,
                          cut_from, ignore_abs_threshold=False, cut_to=False):
        """
        Apply masking to a spectrum based on wavelength range and absorbance thresholds.

        Parameters
        ----------
        wavelength_indices : numpy.ndarray
            Array of wavelength indices
        target_spectrum : numpy.ndarray
            Array of absorbance values
        cut_from : int
            Wavelength index to start the analysis from
        ignore_abs_threshold : bool, optional
            Whether to ignore absorbance thresholds, defaults to False
        cut_to : int or False, optional
            Wavelength index to end the analysis at, defaults to False

        Returns
        -------
        tuple
            (masked_wavelength_indices, masked_target_spectrum)
        """
        threshold_interpolator = interpolate.interp1d(self.thresh_w_indices, self.thresh_as,
                                                      fill_value='extrapolate')

        if not ignore_abs_threshold:
            mask = np.logical_and(target_spectrum < threshold_interpolator(wavelength_indices),
                                  wavelength_indices > cut_from)
        else:
            mask = wavelength_indices > cut_from

        if cut_to:
            mask = np.logical_and(mask, wavelength_indices <= cut_to)

        mask = np.logical_and(mask,
                              target_spectrum > np.min(target_spectrum) + self.lower_limit_of_absorbance)

        return wavelength_indices[mask], target_spectrum[mask]


    def mask_multispectrum(self, wavelength_indices, target_spectrum,
                           cut_from, ignore_abs_threshold=False, cut_to=False,
                           upper_limit_of_absorbance=0.95,
                           artefact_generating_upper_limit_of_absorbance=1.5):
        """
        Apply masking to a spectrum for multispectrum analysis.

        This masking is more sophisticated than mask_the_spectrum, handling artifacts
        that can occur at very high absorbance values.

        Parameters
        ----------
        wavelength_indices : numpy.ndarray
            Array of wavelength indices
        target_spectrum : numpy.ndarray
            Array of absorbance values
        cut_from : int
            Wavelength index to start the analysis from
        ignore_abs_threshold : bool, optional
            Whether to ignore absorbance thresholds, defaults to False
        cut_to : int or False, optional
            Wavelength index to end the analysis at, defaults to False
        upper_limit_of_absorbance : float, optional
            Upper limit for absorbance values, defaults to 0.95
        artefact_generating_upper_limit_of_absorbance : float, optional
            Threshold for identifying artefact-generating regions, defaults to 1.5

        Returns
        -------
        tuple
            (masked_wavelength_indices, masked_target_spectrum)
        """
        mask = wavelength_indices > cut_from
        if cut_to is not None:
            mask = np.logical_and(mask, wavelength_indices < cut_to)
        mask = np.logical_and(mask, target_spectrum < upper_limit_of_absorbance)

        # find the largest index where target_spectrum is above the artefact_generating_upper_limit_of_absorbance
        if len(np.where(target_spectrum > artefact_generating_upper_limit_of_absorbance)[0]) == 0:
            largest_index_above_2 = -1
        else:
            largest_index_above_2 = np.max(np.where(target_spectrum > artefact_generating_upper_limit_of_absorbance)[0])
        # mask all the indices smaller than largest_index_above_1.5
        mask2 = wavelength_indices > largest_index_above_2
        mask = np.logical_and(mask, mask2)

        return wavelength_indices[mask], target_spectrum[mask]


    def multispectrum_to_concentration(self, target_spectrum_inputs, calibration_folder, calibrant_shortnames,
                                       background_model_folder, dilution_factors,
                                       upper_limit_of_absorbance=1000, fig_filename='temp', do_plot=False,  #lower_limit_of_absorbance=0.02
                                       upper_bounds=[np.inf, np.inf], use_line=False, cut_from = 200, ignore_abs_threshold=False,
                                       cut_to = False, ignore_pca_bkg=False, return_errors=False,
                                       use_linear_calibration=True, plot_calibrant_references=False,
                                       return_report=False, sigma_of_absorbance=0.01,
                                       starting_concentration_dict=None, second_dilution_factor_bound_range=0.1,
                                       maximum_wavelength_offset=1.5,
                                       obey_stoichiometric_inequalities=True):
        """
        Calculate concentrations from multiple spectra of the same sample with different dilutions.

        This is a key method for determining component concentrations in complex mixtures by
        using multiple spectra with different dilution factors to capture both high and low
        concentration components accurately. The method performs a sophisticated spectral unmixing
        with constraints and can account for wavelength offsets between spectra.

        Parameters
        ----------
        target_spectrum_inputs : list
            List of input spectra (absorbance arrays)
        calibration_folder : str
            Path to the folder containing calibration data
        calibrant_shortnames : list
            List of short names of the calibrants to include in the unmixing
        background_model_folder : str
            Path to the folder containing background model data
        dilution_factors : list
            List of dilution factors for each spectrum
        upper_limit_of_absorbance : float, optional
            Upper limit for absorbance values, defaults to 1000
        fig_filename : str, optional
            Base filename for saving figures, defaults to 'temp'
        do_plot : bool, optional
            Whether to display plots during processing, defaults to False
        upper_bounds : list, optional
            Upper bounds for calibrant concentrations, defaults to [np.inf, np.inf]
        use_line : bool, optional
            Whether to use a linear term in the model, defaults to False
        cut_from : int, optional
            Wavelength index to start the analysis from, defaults to 200
        ignore_abs_threshold : bool, optional
            Whether to ignore absorbance thresholds, defaults to False
        cut_to : int or False, optional
            Wavelength index to end the analysis at, defaults to False
        ignore_pca_bkg : bool, optional
            Whether to ignore PCA background components, defaults to False
        return_errors : bool, optional
            Whether to return concentration errors, defaults to False
        use_linear_calibration : bool, optional
            Whether to use linear calibration, defaults to True
        plot_calibrant_references : bool, optional
            Whether to plot calibrant reference spectra, defaults to False
        return_report : bool, optional
            Whether to return a detailed report of the unmixing, defaults to False
        sigma_of_absorbance : float, optional
            Standard deviation of absorbance measurements, defaults to 0.01
        starting_concentration_dict : dict or None, optional
            Dictionary with starting concentrations
        second_dilution_factor_bound_range : float, optional
            Allowed range for second dilution factor, defaults to 0.1
        maximum_wavelength_offset : float, optional
            Maximum allowed wavelength offset between spectra, defaults to 1.5
        obey_stoichiometric_inequalities : bool, optional
            Whether to enforce stoichiometric constraints, defaults to True

        Returns
        -------
        numpy.ndarray or tuple
            Array of calculated concentrations, or tuple with additional information
            depending on return_errors and return_report values
        """
        t0 = time.time()
        calibrants = []
        for calibrant_shortname in calibrant_shortnames:
            dict_here = dict()
            dict_here['coeff_to_concentration_interpolator'], dict_here['reference_interpolator'], dict_here[
                'bkg_spectrum'] = \
                self.load_calibration_for_one_calibrant(calibrant_shortname, calibration_folder,
                                                        use_line_fit=use_linear_calibration,
                                                        do_savgol_filtering=False)
            dict_here['concentration_to_coeff_interpolator'], _, _ = \
                self.load_concentration_to_coeff_for_one_calibrant(calibrant_shortname, calibration_folder,
                                                                     use_line_fit=use_linear_calibration)
            calibrants.append(dict_here.copy())

        if plot_calibrant_references:
            for i, calibrant in enumerate(calibrants):
                if i<10:
                    linestyle='-'
                else:
                    linestyle='--'
                xs = 220+np.linspace(0, 400, 400)
                # xs = 220 + np.linspace(150, 400, 400)
                plt.plot(xs, calibrant['reference_interpolator'](xs-220),
                         label=calibrant_shortnames[i], linestyle=linestyle)
            plt.legend()
            plt.xlabel('Wavelength, nm')
            plt.ylabel('Absorbance')
            plt.xlim(220, 500)
            plt.ylim(0, 1.4)
            plt.show()
        print(f'N calibrants: {len(calibrants)}')

        if ignore_pca_bkg:
            bkg_spectrum = np.mean(np.array([calibrant['bkg_spectrum'] for calibrant in calibrants]), axis=0)
        else:
            bkg_spectrum = np.load(background_model_folder + 'bkg_spectrum.npy')
        for target_spectrum_input in target_spectrum_inputs:
            assert len(bkg_spectrum) == len(target_spectrum_input), \
                'Length of background spectrum is not the same as the length of the target spectrum.' \
                'This may be because the wavelengths are not aligned.'

        wavelengths = bkg_spectrum[:, 0]
        target_spectra = [target_spectrum_input - bkg_spectrum[:, 1] for target_spectrum_input in target_spectrum_inputs]
        wavelength_indices = np.arange(calibrants[0]['bkg_spectrum'].shape[0])

        target_spectra_wavelength_indices_masked = []
        target_spectra_amplitudes_masked = []
        for i, target_spectrum in enumerate(target_spectra):
            target_spectrum_wavelengths_masked, target_spectrum_amplitudes_masked = self.mask_multispectrum(
                wavelength_indices, target_spectrum, cut_from, upper_limit_of_absorbance=upper_limit_of_absorbance, cut_to=cut_to)
            target_spectra_wavelength_indices_masked.append(target_spectrum_wavelengths_masked)
            target_spectra_amplitudes_masked.append(target_spectrum_amplitudes_masked)

            if len(target_spectrum_wavelengths_masked) == 0:
                print(f'There is no data that is within mask for spectrum #{i}. Returning zeros.')
                return [0 for i in range(len(calibrant_shortnames))]

        comboX = np.concatenate(target_spectra_wavelength_indices_masked)
        comboY = np.concatenate(target_spectra_amplitudes_masked)

        # uncertainties
        if not self.use_instrumental_sigmas:
            combo_sigmas = np.ones_like(comboY) * sigma_of_absorbance
        else:
            combo_sigmas = []
            for i, wavelength_index in enumerate(comboX):
                wavelength_here = wavelengths[wavelength_index]
                absorbance_here = comboY[i]
                sigma_here_here = self.uncertainty_of_measured_absorbance(wavelength_here, absorbance_here)
                combo_sigmas.append(sigma_here_here)
            combo_sigmas = np.array(combo_sigmas)

        indices_for_splitting = np.cumsum([len(target_spectrum_wavelengths_masked) for target_spectrum_wavelengths_masked in target_spectra_wavelength_indices_masked])[:-1]
        number_of_calibrants = len(calibrant_shortnames)
        number_of_spectra = len(target_spectrum_inputs)

        if not ignore_pca_bkg:
            background_interpolators = [interpolate.interp1d(wavelength_indices,
                                                          np.load(background_model_folder + f'component_{i}.npy'),
                                                          fill_value='extrapolate')
                                     for i in range(1)]
        else:
            background_interpolators = [interpolate.interp1d(wavelength_indices,
                                                             np.ones_like(wavelength_indices),
                                                             fill_value='extrapolate')
                                        for i in range(1)]

        def func_prelim(*args):
            xs = args[0]
            separate_spectra = np.split(xs, indices_for_splitting)
            calibrants_concentrations = args[1:number_of_calibrants + 1]
            dilutions_factors_here = [dilution_factors[0]] + list(args[number_of_calibrants + 1: number_of_calibrants + 1 + number_of_spectra - 1])
            assert len(dilutions_factors_here) == number_of_spectra
            offsets = args[number_of_calibrants + 1 + number_of_spectra - 1:
                           number_of_calibrants + 1 + number_of_spectra - 1 + number_of_spectra]
            bkg_pca_weights = args[number_of_calibrants + 1 + number_of_spectra - 1 + number_of_spectra:
                                   number_of_calibrants + 1 + number_of_spectra - 1 + number_of_spectra + number_of_spectra]
            wavelength_offsets = args[number_of_calibrants + 1 + number_of_spectra - 1 + number_of_spectra + number_of_spectra:
                                      number_of_calibrants + 1 + number_of_spectra - 1 + number_of_spectra + number_of_spectra + number_of_spectra]
            separate_predicted_spectra = []
            for spectrum_index, wavelengths in enumerate(separate_spectra):
                wavelengths = wavelengths + wavelength_offsets[spectrum_index]
                dilution_factor_for_this_spectrum = dilutions_factors_here[spectrum_index]
                calibrants_concentrations_for_this_spectrum = [x / dilution_factor_for_this_spectrum for x in calibrants_concentrations]

                calibrants_coeffs_for_this_spectrum = [np.asscalar(calibrants[i]['concentration_to_coeff_interpolator'](calibrants_concentrations_for_this_spectrum[i]))
                                                       for i in range(number_of_calibrants)]

                predicted_spectrum = np.zeros_like(wavelengths)

                for i in range(number_of_calibrants):
                    predicted_spectrum += calibrants_coeffs_for_this_spectrum[i] * calibrants[i]['reference_interpolator'](wavelengths)

                predicted_spectrum += offsets[spectrum_index] + background_interpolators[0](wavelengths) * bkg_pca_weights[spectrum_index]
                separate_predicted_spectra.append(predicted_spectrum)
            comboY = np.concatenate(separate_predicted_spectra)

            # if stoich_callable is not None:
            #     stoich_penalization_of_cost = stoich_callable(calibrants_concentrations)
            #     comboY[-2] += stoich_penalization_of_cost * 1e-2
            #     comboY[-1] -= stoich_penalization_of_cost * 1e-2

            return comboY

        def func(*args):
            calibrants_concentrations = args[1:number_of_calibrants + 1]
            comboY = func_prelim(*args)
            if starting_concentration_dict is not None:
                stoich_penalization_of_cost = self.stoich_cost(calibrants_concentrations, calibrant_shortnames, starting_concentration_dict)
                comboY[-2] += stoich_penalization_of_cost * combo_sigmas[-2] / self.uncertainty_of_stoichiometric_overspending_ratio
                comboY[-1] -= stoich_penalization_of_cost * combo_sigmas[-1] / self.uncertainty_of_stoichiometric_overspending_ratio

            return comboY

        p0 = []
        lower_bounds = []
        upper_bounds = []
        x_scale = []

        # it stoichiometry is obeyed, limit the upper bound of calibrant concentration to two times the possible
        # concentration if all the starting material is used to make this calibrant only
        if obey_stoichiometric_inequalities:
            calibrant_concentration_upper_bounds = []
            for i, substance_for_fitting in enumerate(calibrant_shortnames):
                limits_here = []
                for s in self.substrates:
                    if self.df_stoich.loc[self.df_stoich['Names'] == substance_for_fitting, s].values[0] == 0:
                        continue
                    limit_by_this_substrate = starting_concentration_dict[s] / self.df_stoich.loc[self.df_stoich['Names'] == substance_for_fitting, s].values[0]
                    limits_here.append(limit_by_this_substrate)
                calibrant_concentration_upper_bounds.append(min(limits_here))
            print(f'Limits: {calibrant_concentration_upper_bounds}')


        for i in range(number_of_calibrants):
            # if calibrant_shortnames[i] == 'bb017':
            #     p0.append(0.0068)
            #     lower_bounds.append(0.0068)
            #     upper_bounds.append(1)
            # else:
            p0.append(1e-7)
            x_scale.append(1e-3)
            # if calibrant_shortnames[i] == 'dm37':
            #     lower_bounds.append(-1*np.inf)
            # else:
            lower_bounds.append(0)
            if obey_stoichiometric_inequalities:
                upper_bounds.append(max([2e-7, 2 * calibrant_concentration_upper_bounds[i]]))
            else:
                upper_bounds.append(np.inf)

        for i in range(number_of_spectra - 1):
            p0.append(dilution_factors[i+1])
            x_scale.append(dilution_factors[i+1])
            lower_bounds.append(dilution_factors[i+1]*(1-second_dilution_factor_bound_range))
            upper_bounds.append(dilution_factors[i+1]*(1+second_dilution_factor_bound_range))

        for i in range(number_of_spectra): # these are offsets
            p0.append(0)
            x_scale.append(0.01)
            lower_bounds.append(-np.inf)
            upper_bounds.append(np.inf)

        if ignore_pca_bkg:
            max_bkg_pca_weight = 1e-12
        else:
            max_bkg_pca_weight = 0.5

        for i in range(number_of_spectra): # these are weights for PCA componetns of the background
            p0.append(0)
            x_scale.append(max_bkg_pca_weight)
            lower_bounds.append(-1*max_bkg_pca_weight)
            upper_bounds.append(max_bkg_pca_weight)

        for i in range(number_of_spectra): # these are wavelength offsets (shifts) of individual spectra
            p0.append(0)
            x_scale.append(1)
            lower_bounds.append(-1*maximum_wavelength_offset)
            upper_bounds.append(maximum_wavelength_offset)

        bounds = (lower_bounds, upper_bounds)
        popt, pcov = curve_fit(func_prelim, comboX, comboY, method='trf',
                               p0=p0, bounds=bounds, sigma=combo_sigmas, absolute_sigma=True,
                               maxfev=100000, ftol=1e-15, xtol=1e-15, gtol=1e-15, verbose=1, x_scale=x_scale)
        if obey_stoichiometric_inequalities:
            print('prelim done')
            try:
                popt, pcov = curve_fit(func, comboX, comboY,
                                       p0=popt, bounds=bounds, sigma=combo_sigmas, absolute_sigma=True,
                                       maxfev=100, ftol=1e-15, xtol=1e-15, gtol=1e-15, verbose=1,
                                       x_scale=x_scale)
            # if max_nfev is reached
            except RuntimeError:
                print(f'RuntimeError, hopefully max_nfev')
                print('Trying again with 1e-12 tolerances')
                try:
                    popt, pcov = curve_fit(func, comboX, comboY, method='trf',
                                           p0=popt, bounds=bounds, sigma=combo_sigmas, absolute_sigma=True,
                                           maxfev=100, ftol=1e-12, xtol=1e-12, gtol=1e-12, verbose=1,
                                           x_scale=x_scale)
                except RuntimeError:
                    print(f'RuntimeError, hopefully max_nfev')
                    print('Trying again with 1e-10 tolerances')
                    try:
                        popt, pcov = curve_fit(func, comboX, comboY, method='trf',
                                               p0=popt, bounds=bounds, sigma=combo_sigmas, absolute_sigma=True,
                                               maxfev=100, ftol=1e-10, xtol=1e-10, gtol=1e-10, verbose=1,
                                               x_scale=x_scale)
                    except RuntimeError:
                        print(f'RuntimeError, hopefully max_nfev')
                        print('Trying again with 1e-6 tolerances')
                        popt, pcov = curve_fit(func, comboX, comboY, method='trf',
                                               p0=popt, bounds=bounds, sigma=combo_sigmas, absolute_sigma=True,
                                               maxfev=100, ftol=1e-6, xtol=1e-6, gtol=1e-6, verbose=1,
                                               x_scale=x_scale)

        # print(f'infodict nfev: {infodict["nfev"]}')
        perr = np.sqrt(np.diag(pcov))  # errors of the fitted coefficients

        concentrations_here = popt[0:number_of_calibrants]

        if obey_stoichiometric_inequalities:
            required_subs = self.product_concentrations_to_required_substrates(concentrations_here, calibrant_shortnames)
            os_string = ''
            for s in self.substrates:
                overspending_ratio = required_subs[s] / starting_concentration_dict[s]
                string_here = f'{s} osr: {overspending_ratio-1:.1%}\n'
                os_string += string_here
            os_string = os_string[:-1]
            print(os_string)
        else:
            os_string = ''

        fitted_dilution_factors = popt[number_of_calibrants: number_of_calibrants + number_of_spectra - 1]
        fitted_offsets = popt[number_of_calibrants + number_of_spectra - 1: number_of_calibrants + number_of_spectra - 1 + number_of_spectra]
        fitted_wavelength_offsets = popt[number_of_calibrants + number_of_spectra - 1 + number_of_spectra + number_of_spectra:
                                      number_of_calibrants + number_of_spectra - 1 + number_of_spectra + number_of_spectra + number_of_spectra]
        print(f'Fitted wavelength offsets: {fitted_wavelength_offsets}')
        logging.debug('Fitted concentrations:', concentrations_here)
        logging.debug('Fitted dilution factors:', fitted_dilution_factors)
        logging.debug('Fitted offsets:', fitted_offsets)
        logging.debug('Popt: ', popt)
        logging.debug('p0: ', p0)

        if return_errors:
            # convert coefficient errors into concentration errors
            upper_confidence_limit = [calibrants[calibrant_index]['coeff_to_concentration_interpolator'](fitted_coeff + perr[calibrant_index])
                               for calibrant_index, fitted_coeff in enumerate(popt[:-4])]
            concentration_errors = [upper_confidence_limit[i] - concentrations_here[i] for i in range(len(concentrations_here))]
            return concentrations_here, concentration_errors

        # plot the fit vs the data

        # make number of subplots equal to number of spectra
        predicted_comboY = func_prelim(comboX, *popt)

        fit_report = dict()
        for calibrant_id, calibrant_shortname in enumerate(calibrant_shortnames):
            fit_report[f'pcerr#{calibrant_shortname}'] = perr[calibrant_id]

        fit_report['rmse'] = np.sqrt(np.mean((predicted_comboY - comboY) ** 2))
        fit_report['maxresidual'] = np.max(np.abs(predicted_comboY - comboY))

        separate_predicted_spectra = np.split(predicted_comboY, indices_for_splitting)
        separate_sigmas = np.split(combo_sigmas, indices_for_splitting)
        fig1, axs = plt.subplots(len(target_spectrum_inputs), 1, figsize=(10, 10), sharex=True)
        for spectrum_index in range(number_of_spectra):
            axs[spectrum_index].plot(220+target_spectra_wavelength_indices_masked[spectrum_index],
                                     target_spectra_amplitudes_masked[spectrum_index], color='red',
                                     label='Data', alpha=0.7)
            axs[spectrum_index].fill_between(x=220 + target_spectra_wavelength_indices_masked[spectrum_index],
                                     y1=target_spectra_amplitudes_masked[spectrum_index] - separate_sigmas[spectrum_index],
                                     y2=target_spectra_amplitudes_masked[spectrum_index] + separate_sigmas[spectrum_index],
                                     color='gold', alpha=0.15)
            axs[spectrum_index].plot(220+target_spectra_wavelength_indices_masked[spectrum_index],
                                     separate_predicted_spectra[spectrum_index], color='black', label='Fit', alpha=0.7)
            residuals_here = separate_predicted_spectra[spectrum_index] - target_spectra_amplitudes_masked[spectrum_index]
            lag = len(residuals_here) // 5
            lb_df = sm.stats.acorr_ljungbox(residuals_here, lags=[lag])
            # print(f'spectrum index {spectrum_index}, LB_pvalue: {lb_df.loc[lag, "lb_pvalue"]}, lag: {lag}')
            # print(lb_df)
            if len(lb_df) == 1:
                # take values from first row of dataframe lb_pvalue
                fit_report[f'LB_pvalue_dil_{spectrum_index}'] = lb_df.loc[lag, 'lb_pvalue']
                fit_report[f'LB_stat_dil_{spectrum_index}'] = lb_df.loc[lag, 'lb_stat']
            # axs[spectrum_index].legend()

        for calibrant_index in range(len(calibrant_shortnames)):
            if calibrant_index <= 9:
                linestyle_here = '-'  # solid line
            else:
                linestyle_here = '--' # dashed line
            cpopt = popt.copy()
            for i in range(number_of_calibrants):
                if i != calibrant_index:
                    cpopt[i] = 0
            predicted_comboY = func(comboX, *cpopt)
            separate_predicted_spectra = np.split(predicted_comboY, indices_for_splitting)
            for spectrum_index in range(number_of_spectra):
                axs[spectrum_index].plot(220+target_spectra_wavelength_indices_masked[spectrum_index], separate_predicted_spectra[spectrum_index],
                                         label=calibrant_shortnames[calibrant_index], linestyle=linestyle_here)

        weights_of_background_components = popt[number_of_calibrants + number_of_spectra - 1 + number_of_spectra:
                                                number_of_calibrants + number_of_spectra - 1 + number_of_spectra + number_of_spectra]
        print(f'weights_of_background_components: {weights_of_background_components}')
        for spectrum_index in range(number_of_spectra):
            # plot background
            axs[spectrum_index].plot(220+target_spectra_wavelength_indices_masked[spectrum_index],
                                     background_interpolators[0](target_spectra_wavelength_indices_masked[spectrum_index]) * weights_of_background_components[spectrum_index],
                                     label='Bkg. PC1', linestyle=':')
        # for spectrum_index in range(number_of_spectra):
        #     axs[spectrum_index].legend()
        axs[0].legend()

        plt.xlabel('Wavelength, nm')
        plt.ylabel('Absorbance')

        axs[0].set_title(os_string)

        fit_report['fitted_dilution_factor_2'] = fitted_dilution_factors[0]

        # plot covariance matrix
        # # plt.legend()
        # plt.show()
        #
        # plt.figure(figsize=(5, 10))
        # pcov_to_plot = pcov[:len(calibrant_shortnames), :len(calibrant_shortnames)]
        # plt.imshow(pcov_to_plot, vmin=-1*max(np.abs(pcov_to_plot).flatten()), vmax=max(np.abs(pcov_to_plot).flatten()),
        #            cmap='RdBu_r')
        # # make tick labels from calibrant_shortnames
        # plt.yticks(range(len(calibrant_shortnames)), calibrant_shortnames)
        # plt.xticks(range(len(calibrant_shortnames)), calibrant_shortnames, rotation=90)
        # plt.colorbar(orientation='vertical', fraction=0.046)
        # plt.tight_layout()
        # plt.show()

        fig1.savefig(f"{fig_filename}.png")

        print(f'Time spent for this unmixing: {time.time() - t0} seconds.')
        if do_plot:
            plt.show()
        else:
            plt.close(fig1)
            plt.close('all')
            plt.clf()

        if return_report:
            return concentrations_here, fit_report
        else:
            return concentrations_here

    def uncertainty_of_measured_absorbance(self, wavelength, absorbance, lower_threshold_of_sigma=0.005):
        """
        Calculate the uncertainty of a measured absorbance value.

        Parameters
        ----------
        wavelength : float
            Wavelength (nm)
        absorbance : float
            Absorbance value
        lower_threshold_of_sigma : float, optional
            Minimum uncertainty threshold, defaults to 0.005

        Returns
        -------
        float
            Estimated uncertainty of the absorbance measurement
        """
        variance = self.sigma_interpolator(wavelength, absorbance)[0][0]
        if variance < lower_threshold_of_sigma**2:
            variance = lower_threshold_of_sigma**2
        # if wavelength <= 270:
        #     variance = variance * (1 + 200*np.exp(-1*(wavelength - 220)/20))
        return np.sqrt(variance)

    def product_concentrations_to_required_substrates(self, concentrations, calibrant_shortnames):
        """
        Calculate the required substrate concentrations for given product concentrations,
        based on stoichiometric coefficients.

        Parameters
        ----------
        concentrations : list
            List of product concentrations
        calibrant_shortnames : list
            List of short names of the products

        Returns
        -------
        dict
            Dictionary mapping substrate names to required concentrations
        """
        dict_of_substrate_concentrations = {s: 0 for s in self.substrates}
        for i, substance_for_fitting in enumerate(calibrant_shortnames):
            for s in self.substrates:
                # find the cell in self.df_stoich where Name is equal to substance_for_fitting and column is s
                dict_of_substrate_concentrations[s] += concentrations[i] * self.df_stoich.loc[
                    self.df_stoich['Names'] == substance_for_fitting, s].values[0]
        return dict_of_substrate_concentrations

    def stoich_cost(self, concentrations, calibrant_shortnames, starting_concentrations_dict):
        """
        Calculate the stoichiometric cost (penalization) for a set of concentrations.

        This is used to enforce stoichiometric constraints in the spectral unmixing algorithm.

        Parameters
        ----------
        concentrations : list
            List of product concentrations
        calibrant_shortnames : list
            List of short names of the products
        starting_concentrations_dict : dict
            Dictionary mapping substrate names to starting concentrations

        Returns
        -------
        float
            Stoichiometric cost (penalization value)
        """
        def smooth_step(x):
            if x <= 0:
                return 0
            else:
                return x
            # return np.exp(x)

        required_subs = self.product_concentrations_to_required_substrates(concentrations, calibrant_shortnames)
        final_penalization = 0
        for s in self.substrates:
            if starting_concentrations_dict[s] > 0:
                overspending_ratio = required_subs[s] / starting_concentrations_dict[s]
                final_penalization += smooth_step((overspending_ratio - 1))
            # final_penalization += smooth_step((required_subs[s] - starting_concentrations_dict[s])/0.001)

        # find concentration of acetic acid
        # acetic_mixmatch = np.abs(required_subs['ammonium_acetate'] - concentrations[calibrant_shortnames.index('acetic_acid')])
        # final_penalization += acetic_mixmatch/0.001

        return final_penalization


def plot_differential_absorbances_for_plate(craic_exp_name,
                                            wavelength,
                                            ref_wavelengths,
                                            ):
    """
    Plots the difference between absorbance at the target wavelength and the mean absorbance at reference wavelengths
    from ref_wavelength list.

    Parameters
    ----------
    craic_exp_name: str
        Name of the folder with CRAIC microspectrometer measurements.
    wavelength
        Target wavelength at which the absorbance is calculated.
    ref_wavelengths
        List of reference wavelengths. Mean absorbance at these wavelengths is subtracted from the absorbance at the
        target wavelength.

    Returns
    -------
    diff: np.array
        Array of differential absorbances.
    """
    sp = SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                         '2022-12-01/interpolator-dataset/')
    craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'
    sp.show_all_spectra(craic_folder + craic_exp_name + '/')
    plt.show()
    diff = sp.get_absorbance_at_single_wavelength_for_one_plate(craic_folder + craic_exp_name + '/',
                                                                wavelength=wavelength,
                                                                ref_wavelengths=ref_wavelengths)
    diluted_indices = [i + j for i in [9, 27, 45] for j in range(9)]
    undiluted_indices = [i + j for i in [0, 18, 36] for j in range(9)]
    diff = diff[diluted_indices]
    print(f'rel.std {np.std(diff) / np.mean(diff)}')
    plt.plot(diff)
    plt.xlabel('Vial ID')
    plt.ylabel(f'Absorbance at {wavelength} nm minus absorbance at wavelengths {ref_wavelengths} nm')
    plt.title(f'{craic_exp_name}.\nRel. STD: {np.std(diff) / np.mean(diff)}')
    plt.show()
    return diff




if __name__ == '__main__':
    pass


