import importlib
import pickle
from distutils import dir_util
import pandas as pd
from pytest import fixture
import os
import numpy as np
pwc = importlib.import_module("uv-vis-absorption-spectroscopy.process_wellplate_spectra")


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
            folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
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
        assert np.isclose(concentrations_here, expected_concentrations_here).all()