from robowski.settings import *
import importlib
import pickle
from distutils import dir_util
import pandas as pd
from pytest import fixture, mark
import os
import numpy as np
import matplotlib.pyplot as plt

plt.ioff()  # Turn off interactive plotting for testing
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as pwc
import shutil
import json
import pytest

# Global flag to control saving of expected outputs
# Set this to True only when you intentionally want to generate/update expected outputs
# Can be controlled via environment variable
# SAVE_EXPECTED_OUTPUTS = os.environ.get('SAVE_EXPECTED_OUTPUTS', '').lower() == 'true'
# For now, it is controlled by this hardcoded value:
SAVE_EXPECTED_OUTPUTS = False


@fixture
def datadir(tmpdir, request):
    '''
    Fixture responsible for searching a folder with the same name of test module and, if available, moving all
    contents to a temporary directory so tests can use them freely.
    '''
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, str(tmpdir))

    return tmpdir


@fixture
def setup_processor(datadir):
    """
    Set up a SpectraProcessor instance for testing.
    """
    with datadir.as_cwd():
        sp = pwc.SpectraProcessor(
            folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                            '2022-12-01/interpolator-dataset/')
        sp.nanodrop_lower_cutoff_of_wavelengths = 250
        sigma_interpolator_filename = (f'nanodrop-spectrophotometer-measurements/'
                                       f'nanodrop_errorbar_folder_2024-03-16/bivariate_spline_interpolator.pkl')
        with open(sigma_interpolator_filename, 'rb') as f:
            sp.sigma_interpolator = pickle.load(f)
        sp.use_instrumental_sigmas = True

        # Create output directories if they don't exist
        os.makedirs('expected_outputs', exist_ok=True)

        return sp


def save_and_verify_expected_output(data, filename, output_dir='expected_outputs'):
    """
    Helper function to save expected outputs for regression testing AND verify they can be loaded back

    Only saves if SAVE_EXPECTED_OUTPUTS is True, but always verifies against existing files if they exist
    """
    filepath = os.path.join(output_dir, filename)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save data if we're in generation mode
    if SAVE_EXPECTED_OUTPUTS:
        if isinstance(data, np.ndarray):
            np.save(filepath, data)
        elif isinstance(data, pd.DataFrame):
            data.to_pickle(filepath)
        else:
            with open(filepath, 'wb') as f:
                pickle.dump(data, f)

        print(f"Saved expected output to {filepath}")

    # Always verify against existing file if it exists
    if os.path.exists(filepath):
        if isinstance(data, np.ndarray):
            expected_data = np.load(filepath, allow_pickle=True)
            assert np.allclose(data, expected_data, equal_nan=True), f"Data mismatch with {filepath}"
        elif isinstance(data, pd.DataFrame):
            expected_data = pd.read_pickle(filepath)
            # reload expected data, but use index saved in the file as the index of the dataframe
            pd.testing.assert_frame_equal(data, expected_data)
        else:
            with open(filepath, 'rb') as f:
                expected_data = pickle.load(f)

            if isinstance(data, dict):
                assert data.keys() == expected_data.keys(), f"Dict keys mismatch with {filepath}"
                for key in data:
                    if isinstance(data[key], np.ndarray):
                        assert np.allclose(data[key], expected_data[key],
                                           equal_nan=True), f"Dict value mismatch for key {key}"
                    else:
                        assert data[key] == expected_data[key], f"Dict value mismatch for key {key}"
            else:
                assert data == expected_data, f"Data mismatch with {filepath}"

        return expected_data
    else:
        if not SAVE_EXPECTED_OUTPUTS:
            pytest.fail(f"Expected output file {filepath} does not exist and SAVE_EXPECTED_OUTPUTS=False")
        return data


def test_concentrations_for_one_plate_end2end(datadir):
    """
    Test uses a fixture that copies all structure from `tests/test_process_wellplate_spectra` directory into a temporary
    directory, which is later treated as the data_folder (that is normally in the Dropbox, but not for tests). Test
    loads the concentrations for one plate from respective locations in the
    `simple-reactions/2023-09-06-run01/` in the temporary folder and then checks the results against an expected
    concentrations loaded from `expected_outputs/concentrations_for_simple_reactions_2023-09-06-run01_plate50.npy`.

    Parameters
    ----------
    datadir: pytest fixture
        Temporary directory with the same structure as `tests/process_wellplate_spectra` directory.
    """
    data_folder = ''
    with datadir.as_cwd():
        sp = pwc.SpectraProcessor(
            folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                            '2022-12-01/interpolator-dataset/')
        sp.nanodrop_lower_cutoff_of_wavelengths = 250
        sp.spectrum_data_type = 'nanodrop'
        run_name = 'simple-reactions/2023-09-06-run01/'
        concentrations_here = sp.concentrations_for_one_plate(experiment_folder=data_folder + run_name,
                                                              plate_folder=data_folder + 'simple-reactions/2023-09-06-run01/nanodrop_spectra/2023-09-06_20-29-24_plate_50.csv',
                                                              calibration_folder=data_folder + 'simple-reactions/2023-09-06-run01/' + 'microspectrometer_data/calibration/',
                                                              calibrant_shortnames=['E1DB02', 'E1OH02'],
                                                              background_model_folder=data_folder + 'simple-reactions/2023-09-06-run01/microspectrometer_data/background_model/',
                                                              calibrant_upper_bounds=[np.inf, np.inf],
                                                              do_plot=False, return_all_substances=True,
                                                              cut_from=50, cut_to=False,
                                                              ignore_abs_threshold=True, ignore_pca_bkg=True)
        expected_concentrations_here = np.load(
            'expected_outputs/concentrations_for_simple_reactions_2023-09-06-run01_plate50.npy')
        # use np.isclose() to assert equality
        if not np.isclose(concentrations_here, expected_concentrations_here).all():
            # print the differences
            print("Differences:")
            print(concentrations_here - expected_concentrations_here)

        assert np.isclose(concentrations_here, expected_concentrations_here).all()


def test_multispectrum_to_concentration_end2end(datadir):
    """
    Test uses a fixture that copies all structure from `tests/test_process_wellplate_spectra` directory into a temporary
    directory, which is later treated as the data_folder (that is normally in the Dropbox, but not for tests). Test
    loads the concentrations for one plate from respective locations in the
    `simple-reactions/2023-09-06-run01/` in the temporary folder and then checks the results against an expected
    concentrations loaded from `expected_outputs/concentrations_for_simple_reactions_2023-09-06-run01_plate50.npy`.

    Parameters
    ----------
    datadir: pytest fixture
        Temporary directory with the same structure as `tests/process_wellplate_spectra` directory.
    """
    data_folder = ''
    with datadir.as_cwd():
        sp = pwc.SpectraProcessor(
            folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                            '2022-12-01/interpolator-dataset/')
        sigma_interpolator_filename = (f'{data_folder}nanodrop-spectrophotometer-measurements/'
                                       f'nanodrop_errorbar_folder_2024-03-16/bivariate_spline_interpolator.pkl')
        with open(sigma_interpolator_filename, 'rb') as f:
            sp.sigma_interpolator = pickle.load(f)
        sp.nanodrop_lower_cutoff_of_wavelengths = 220
        sp.use_instrumental_sigmas = True

        # well_id = 44
        well_id = 9
        substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'dm36', 'dm37', 'dm40_12',
                                  'dm40_10', 'ethyl_acetoacetate', 'EAB', 'bb017', 'bb021', 'dm70', 'dm053', 'dm088_4',
                                  'bb021_f2']
        # cut_from = 40
        cut_from = 0
        # Condition 154
        # plate_folder = data_folder + 'BPRF/2024-01-08-run01/nanodrop_spectra/2024-01-10_12-51-07_UV-Vis_plate_71.csv'
        # plate_folder = data_folder + 'BPRF/2024-01-08-run02/nanodrop_spectra/2024-01-10_17-10-28_UV-Vis_plate_61.csv'
        # plate_folder = data_folder + 'BPRF/2024-01-17-run01/nanodrop_spectra/2024-01-19_12-22-47_UV-Vis_plate_66.csv'
        plate_folder = data_folder + 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-22-22_UV-Vis_plate_92.csv'
        # plate_folder = data_folder + 'BPRF/2024-02-16-run01/nanodrop_spectra/2024-02-18_17-48-07_UV-Vis_plate74.csv'
        spectrum1 = sp.load_msp_by_id(
            plate_folder=plate_folder,
            well_id=well_id)[:, 1]

        # plate_folder = data_folder + 'BPRF/2024-01-08-run01/nanodrop_spectra/2024-01-10_13-48-13_UV-Vis_plate_73.csv'
        # plate_folder = data_folder + 'BPRF/2024-01-08-run02/nanodrop_spectra/2024-01-10_17-55-20_UV-Vis_plate_66.csv'
        # plate_folder = data_folder + 'BPRF/2024-01-17-run01/nanodrop_spectra/2024-01-19_13-00-17_UV-Vis_plate_67.csv'
        plate_folder = data_folder + 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-44-22_UV-Vis_plate_93.csv'
        # plate_folder = data_folder + 'BPRF/2024-02-16-run01/nanodrop_spectra/2024-02-18_18-02-42_UV-Vis_plate76.csv'
        spectrum2 = sp.load_msp_by_id(
            plate_folder=plate_folder,
            well_id=well_id)[:, 1]

        concentrations_here = sp.multispectrum_to_concentration(target_spectrum_inputs=[spectrum1, spectrum2],
                                                                dilution_factors=[20, 200],
                                                                calibration_folder=data_folder + 'BPRF/2024-01-17-run01/' + 'microspectrometer_data/calibration/',
                                                                calibrant_shortnames=substances_for_fitting,
                                                                background_model_folder=data_folder + 'BPRF/cross_conamination_and_backgound_test/ethanol_background_model/',
                                                                upper_bounds=[np.inf, np.inf, np.inf, np.inf, np.inf,
                                                                              np.inf,
                                                                              np.inf, np.inf, np.inf, np.inf],
                                                                do_plot=False, cut_from=cut_from, cut_to=250,
                                                                ignore_abs_threshold=False, ignore_pca_bkg=False,
                                                                plot_calibrant_references=False,
                                                                upper_limit_of_absorbance=0.95,
                                                                obey_stoichiometric_inequalities=False)

        expected_concentrations_here = np.load(
            'expected_outputs/multispectrum_to_concentration_BPRF_2024-01-17-run01_plates_92_and_93.npy')
        # use np.isclose() to assert equality
        if not np.isclose(concentrations_here, expected_concentrations_here).all():
            # print the differences
            print("Differences:")
            print(concentrations_here - expected_concentrations_here)

        assert np.isclose(concentrations_here, expected_concentrations_here).all()


# =============== NEW TESTS BELOW THIS LINE ===============


def test_load_msp_file(datadir, setup_processor):
    """Test the load_msp_file function"""
    with datadir.as_cwd():
        # Test loading an MSP file
        experimental_data_filename = 'uv-vis-absorption-spectroscopy/microspectrometer-calibration/2022-12-01/2-inch-nd-calibrations/0p1.msp'

        # Since we're testing with numpy files instead of .msp files, let's modify to load numpy directly
        input_spectrum = pwc.load_msp_file(experimental_data_filename)

    # Save and verify against expected output
    expected_spectrum = save_and_verify_expected_output(
        input_spectrum,
        'load_msp_file_HRP01_ref_spectrum.npy'
    )

    # Additional verification
    assert np.array_equal(input_spectrum, expected_spectrum)
    assert isinstance(input_spectrum, np.ndarray)
    assert len(input_spectrum.shape) >= 1


def test_load_nanodrop_csv_for_one_plate(datadir, setup_processor):
    """Test loading nanodrop CSV data"""
    with datadir.as_cwd():
        sp = setup_processor

        # Test loading a nanodrop CSV file
        plate_folder = 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-22-22_UV-Vis_plate_92.csv'
        nanodrop_df = sp.load_nanodrop_csv_for_one_plate(plate_folder=plate_folder)

        # Save and verify against expected output
        expected_df = save_and_verify_expected_output(
            nanodrop_df,
            'df_nanodrop_df_plate_92.pickle'
        )

        # Additional verification
        assert isinstance(nanodrop_df, pd.DataFrame)
        assert 'wavelength' in nanodrop_df.columns
        assert nanodrop_df['wavelength'].min() >= sp.nanodrop_lower_cutoff_of_wavelengths
        assert nanodrop_df['wavelength'].max() <= sp.nanodrop_upper_cutoff_of_wavelengths
        pd.testing.assert_frame_equal(nanodrop_df, expected_df)


def test_load_single_nanodrop_spectrum(datadir, setup_processor):
    """Test loading a single nanodrop spectrum"""
    with datadir.as_cwd():
        sp = setup_processor

        # Test loading a single spectrum
        plate_folder = 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-22-22_UV-Vis_plate_92.csv'
        well_id = 9
        spectrum = sp.load_single_nanodrop_spectrum(plate_folder=plate_folder, well_id=well_id)

        # Save and verify against expected output
        expected_spectrum = save_and_verify_expected_output(
            spectrum,
            f'nanodrop_spectrum_plate_92_well_{well_id}.npy'
        )

        # Additional verification
        assert np.array_equal(spectrum, expected_spectrum)
        assert isinstance(spectrum, np.ndarray)
        assert spectrum.shape[1] == 2  # [wavelengths, absorbances]


def test_load_msp_by_id(datadir, setup_processor):
    """Test loading a spectrum by well ID"""
    with datadir.as_cwd():
        sp = setup_processor

        # Test with CRAIC data
        plate_folder = 'craic_microspectrometer_measurements/absorbance/2023-06-29_07-32-24__plate0000036__multicomp_reactions_2023-06-28-run01/'
        well_id = 9
        spectrum = sp.load_msp_by_id(plate_folder=plate_folder, well_id=well_id)

        # Save and verify against expected output
        expected_spectrum = save_and_verify_expected_output(
            spectrum,
            f'msp_spectrum_2023-06-29_07-32-24__plate0000036__multicomp_reactions_2023-06-28-run01_well_{well_id}.npy'
        )

        # Additional verification
        assert np.array_equal(spectrum, expected_spectrum)
        assert isinstance(spectrum, np.ndarray)
        assert spectrum.shape[1] == 2  # [wavelengths, absorbances]

        # Test with nanodrop data
        plate_folder = 'simple-reactions/2023-09-06-run01/nanodrop_spectra/2023-09-06_21-50-04_plate_51.csv'
        well_id = 9
        spectrum = sp.load_msp_by_id(plate_folder=plate_folder, well_id=well_id)

        # Save and verify against expected output
        expected_spectrum = save_and_verify_expected_output(
            spectrum,
            f'nanodrop_spectrum_plate_nanodrop_spectra_2023-09-06_21-50-04_plate_51_{well_id}.npy'
        )

        # Additional verification
        assert np.array_equal(spectrum, expected_spectrum)
        assert isinstance(spectrum, np.ndarray)
        assert spectrum.shape[1] == 2  # [wavelengths, absorbances]


def test_mask_multispectrum(datadir, setup_processor):
    """Test the spectrum masking functionality"""
    with datadir.as_cwd():
        sp = setup_processor

        # Load a spectrum
        plate_folder = 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-22-22_UV-Vis_plate_92.csv'
        well_id = 9
        spectrum = sp.load_msp_by_id(plate_folder=plate_folder, well_id=well_id)

        # Create wavelength indices and get the target spectrum
        wavelength_indices = np.arange(spectrum.shape[0])
        target_spectrum = spectrum[:, 1]

        # Test masking with different parameters
        cut_from = 50
        upper_limit = 0.95

        masked_wavelengths, masked_spectrum = sp.mask_multispectrum(
            wavelength_indices,
            target_spectrum,
            cut_from=cut_from,
            upper_limit_of_absorbance=upper_limit
        )

        # Save and verify against expected output
        result_data = {
            'masked_wavelengths': masked_wavelengths,
            'masked_spectrum': masked_spectrum,
            'original_wavelengths': wavelength_indices,
            'original_spectrum': target_spectrum
        }

        expected_data = save_and_verify_expected_output(
            result_data,
            f'masked_spectrum_plate_92_well_{well_id}.pkl'
        )

        # Additional verification
        assert np.array_equal(masked_wavelengths, expected_data['masked_wavelengths'])
        assert np.array_equal(masked_spectrum, expected_data['masked_spectrum'])
        assert len(masked_wavelengths) == len(masked_spectrum)
        assert len(masked_wavelengths) <= len(wavelength_indices)
        assert np.all(masked_wavelengths >= cut_from)
        assert np.all(masked_spectrum <= upper_limit)


@mark.parametrize("well_id", [9, 10, 11])
def test_multispectrum_to_concentration_multiple_wells(datadir, setup_processor, well_id):
    """Test multispectrum_to_concentration with multiple well IDs"""
    with datadir.as_cwd():
        sp = setup_processor
        sp.nanodrop_lower_cutoff_of_wavelengths = 220

        substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'dm36', 'dm37', 'dm40_12',
                                  'dm40_10', 'ethyl_acetoacetate', 'EAB', 'bb017', 'bb021', 'dm70', 'dm053',
                                  'dm088_4',
                                  'bb021_f2']
        cut_from = 0

        # Load spectra
        plate_folder1 = 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-22-22_UV-Vis_plate_92.csv'
        plate_folder2 = 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-44-22_UV-Vis_plate_93.csv'

        spectrum1 = sp.load_msp_by_id(plate_folder=plate_folder1, well_id=well_id)[:, 1]
        spectrum2 = sp.load_msp_by_id(plate_folder=plate_folder2, well_id=well_id)[:, 1]

        # Calculate concentrations
        concentrations = sp.multispectrum_to_concentration(
            target_spectrum_inputs=[spectrum1, spectrum2],
            dilution_factors=[20, 200],
            calibration_folder='BPRF/2024-01-17-run01/microspectrometer_data/calibration/',
            calibrant_shortnames=substances_for_fitting,
            background_model_folder='BPRF/cross_conamination_and_backgound_test/ethanol_background_model/',
            upper_bounds=[np.inf] * 10,
            do_plot=False,
            cut_from=cut_from,
            cut_to=250,
            ignore_abs_threshold=False,
            ignore_pca_bkg=False,
            plot_calibrant_references=False,
            upper_limit_of_absorbance=0.95,
            obey_stoichiometric_inequalities=False
        )

        # Save and verify against expected output
        expected_concentrations = save_and_verify_expected_output(
            concentrations,
            f'multispectrum_concentration_well_{well_id}.npy'
        )

        # Additional verification
        assert np.allclose(concentrations, expected_concentrations)
        assert isinstance(concentrations, np.ndarray)
        assert len(concentrations) == len(substances_for_fitting)


@mark.parametrize("return_report", [False, True])
def test_multispectrum_to_concentration_with_report(datadir, setup_processor, return_report):
    """Test multispectrum_to_concentration with report option"""
    with datadir.as_cwd():
        sp = setup_processor
        sp.nanodrop_lower_cutoff_of_wavelengths = 220

        substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'dm36']
        cut_from = 0
        well_id = 9

        # Load spectra
        plate_folder1 = 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-22-22_UV-Vis_plate_92.csv'
        plate_folder2 = 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-44-22_UV-Vis_plate_93.csv'

        spectrum1 = sp.load_msp_by_id(plate_folder=plate_folder1, well_id=well_id)[:, 1]
        spectrum2 = sp.load_msp_by_id(plate_folder=plate_folder2, well_id=well_id)[:, 1]

        # Calculate concentrations with or without report
        result = sp.multispectrum_to_concentration(
            target_spectrum_inputs=[spectrum1, spectrum2],
            dilution_factors=[20, 200],
            calibration_folder='BPRF/2024-01-17-run01/microspectrometer_data/calibration/',
            calibrant_shortnames=substances_for_fitting,
            background_model_folder='BPRF/cross_conamination_and_backgound_test/ethanol_background_model/',
            upper_bounds=[np.inf] * 5,
            do_plot=False,
            cut_from=cut_from,
            cut_to=250,
            ignore_abs_threshold=False,
            ignore_pca_bkg=False,
            plot_calibrant_references=False,
            upper_limit_of_absorbance=0.95,
            obey_stoichiometric_inequalities=False,
            return_report=return_report
        )

        # Save and verify expected output based on return type
        if return_report:
            concentrations, report = result
            expected_concentrations = save_and_verify_expected_output(
                concentrations,
                f'multispectrum_concentration_with_report_well_{well_id}.npy'
            )
            expected_report = save_and_verify_expected_output(
                report,
                f'multispectrum_report_well_{well_id}.pkl'
            )

            # Verify report structure and content
            assert report == expected_report
            assert isinstance(report, dict)
            assert 'rmse' in report
            assert 'maxresidual' in report
        else:
            concentrations = result
            expected_concentrations = save_and_verify_expected_output(
                concentrations,
                f'multispectrum_concentration_without_report_well_{well_id}.npy'
            )

        # Verify concentrations
        assert np.allclose(concentrations, expected_concentrations)
        assert isinstance(concentrations, np.ndarray)
        assert len(concentrations) == len(substances_for_fitting)


def test_uncertainty_of_measured_absorbance(datadir, setup_processor):
    """Test uncertainty calculation for absorbance measurements"""
    with datadir.as_cwd():
        sp = setup_processor

        # Test various wavelength and absorbance combinations
        test_points = [
            (220, 0.1),
            (300, 0.5),
            (400, 1.0),
            (500, 1.5),
            (600, 2.0)
        ]

        results = {}
        for wavelength, absorbance in test_points:
            uncertainty = sp.uncertainty_of_measured_absorbance(wavelength, absorbance)
            results[(wavelength, absorbance)] = uncertainty

        # Save and verify expected output
        expected_results = save_and_verify_expected_output(
            results,
            'uncertainty_measurements.pkl'
        )

        # Verify results
        assert results == expected_results
        for point, uncertainty in results.items():
            assert uncertainty >= 0  # Uncertainty should be non-negative
            assert isinstance(uncertainty, float)


def test_get_absorbance_at_single_wavelength_for_one_plate(datadir, setup_processor):
    """Test getting absorbance at a single wavelength"""
    with datadir.as_cwd():
        sp = setup_processor

        # Test with nanodrop data
        plate_folder = 'simple-reactions/2023-09-06-run01/nanodrop_spectra/2023-09-06_21-50-04_plate_51.csv'

        # Test with specific wavelength
        wavelength = 300
        ref_wavelengths = [500, 550]
        absorbances = sp.get_absorbance_at_single_wavelength_for_one_plate(
            plate_folder=plate_folder,
            wavelength=wavelength,
            ref_wavelengths=ref_wavelengths
        )

        # Save and verify expected output
        expected_absorbances = save_and_verify_expected_output(
            absorbances,
            f'absorbance_at_{wavelength}nm_plate_51.npy'
        )

        # Test with wavelength ID
        wavelength_id = 100
        ref_wavelength_id = [300]
        absorbances_by_id = sp.get_absorbance_at_single_wavelength_for_one_plate(
            plate_folder=plate_folder,
            wavelength_id=wavelength_id,
            ref_wavelength_id=ref_wavelength_id
        )

        # Save and verify expected output
        expected_absorbances_by_id = save_and_verify_expected_output(
            absorbances_by_id,
            f'absorbance_at_id_{wavelength_id}_plate_51.npy'
        )

        # Verify outputs
        assert np.array_equal(absorbances, expected_absorbances)
        assert np.array_equal(absorbances_by_id, expected_absorbances_by_id)
        assert isinstance(absorbances, np.ndarray)
        assert isinstance(absorbances_by_id, np.ndarray)


def test_product_concentrations_to_required_substrates(datadir, setup_processor):
    """Test the stoichiometry calculation function"""
    with datadir.as_cwd():
        sp = setup_processor

        # Test with some example concentrations and calibrants
        concentrations = [0.01, 0.02, 0.005]
        calibrant_shortnames = ['HRP01', 'dm35_8', 'dm35_9']

        # Calculate required substrates
        required_substrates = sp.product_concentrations_to_required_substrates(
            concentrations, calibrant_shortnames
        )

        # Save and verify expected output
        expected_substrates = save_and_verify_expected_output(
            required_substrates,
            'required_substrates_example.pkl'
        )

        # Verify results
        assert required_substrates == expected_substrates
        assert isinstance(required_substrates, dict)
        assert all(substrate in required_substrates for substrate in sp.substrates)


def generate_validation_report(datadir):
    """Generate a validation report of all saved expected outputs"""
    with datadir.as_cwd():
        output_dir = 'expected_outputs'
        if not os.path.exists(output_dir):
            return "No expected outputs directory found"

        files = os.listdir(output_dir)

        report = {
            "timestamp": pd.Timestamp.now().isoformat(),
            "total_files": len(files),
            "file_details": []
        }

        for file in files:
            file_path = os.path.join(output_dir, file)
            file_info = {
                "filename": file,
                "size_bytes": os.path.getsize(file_path),
                "last_modified": pd.Timestamp.fromtimestamp(os.path.getmtime(file_path)).isoformat()
            }

            # Try to get summary info based on file type
            if file.endswith('.npy'):
                try:
                    data = np.load(file_path, allow_pickle=True)
                    file_info["shape"] = data.shape if hasattr(data, 'shape') else "Not an array"
                    file_info["dtype"] = str(data.dtype) if hasattr(data, 'dtype') else "N/A"
                except Exception as e:
                    file_info["error"] = str(e)

            elif file.endswith('.csv'):
                try:
                    df = pd.read_csv(file_path)
                    file_info["rows"] = len(df)
                    file_info["columns"] = list(df.columns)
                except Exception as e:
                    file_info["error"] = str(e)

            elif file.endswith('.pkl'):
                try:
                    with open(file_path, 'rb') as f:
                        data = pickle.load(f)
                    file_info["type"] = type(data).__name__
                    if isinstance(data, dict):
                        file_info["keys"] = list(data.keys())
                    elif isinstance(data, list):
                        file_info["length"] = len(data)
                except Exception as e:
                    file_info["error"] = str(e)

            report["file_details"].append(file_info)

        # Save the report
        report_path = os.path.join(output_dir, 'validation_report.json')
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)

        print(f"Validation report saved to {report_path}")
        return report


def test_generate_validation_report(datadir):
    """Test generating a validation report for all expected outputs"""
    with datadir.as_cwd():
        report = generate_validation_report(datadir)
        assert report is not None
