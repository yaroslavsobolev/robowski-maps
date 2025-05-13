from robowski.settings import *
import importlib
import pickle
from distutils import dir_util
import pandas as pd
from pytest import fixture
import os
import numpy as np
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as pwc


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
        expected_concentrations_here = np.load('expected_outputs/concentrations_for_simple_reactions_2023-09-06-run01_plate50.npy')
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
                                                           upper_bounds=[np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
                                                                         np.inf, np.inf, np.inf, np.inf],
                                                           do_plot=False, cut_from=cut_from, cut_to=250,
                                                           ignore_abs_threshold=False, ignore_pca_bkg=False,
                                                           plot_calibrant_references=False,
                                                           upper_limit_of_absorbance=0.95,
                                                           obey_stoichiometric_inequalities=False)

        expected_concentrations_here = np.load('expected_outputs/multispectrum_to_concentration_BPRF_2024-01-17-run01_plates_92_and_93.npy')
        # use np.isclose() to assert equality
        if not np.isclose(concentrations_here, expected_concentrations_here).all():
            # print the differences
            print("Differences:")
            print(concentrations_here - expected_concentrations_here)

        assert np.isclose(concentrations_here, expected_concentrations_here).all()