from robowski.settings import *
from distutils import dir_util
import pandas as pd
from pytest import fixture, mark
from pandas.testing import assert_frame_equal
import os
import matplotlib.pyplot as plt
import robowski.uv_vis_absorption_spectroscopy.examples.versatility.versatility_examples as versatility_examples

plt.ioff()  # Turn off interactive plotting for testing
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra
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

def test_unmixing_for_claisen(datadir):
    with datadir.as_cwd():
        data_folder = ''
        experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Claisen_WaiShing/'

        versatility_examples.construct_calibrant(
            cut_from=5,
            lower_limit_of_absorbance=0.007,
            concentration_column_name='concentration',
            do_plot=False,
            calibration_source_filename='2023-10-11_16-26-58_UV-Vis_methoxychalcone',
            calibrant_shortnames=['methoxychalcone'],
            ref_concentrations=[0.0002],
            max_concentrations=[0.0006],
            experiment_name=experiment_name,
        )

        versatility_examples.construct_calibrant(
            cut_from=5,
            lower_limit_of_absorbance=0.007,
            concentration_column_name='concentration',
            do_plot=False,
            calibration_source_filename='2023-10-11_19-50-36_UV-Vis_anisaldehyde',
            calibrant_shortnames=['anisaldehyde', 'acetophenone'],
            ref_concentrations=[0.0002, 0.0003],
            max_concentrations=[0.001, 0.001],
            experiment_name=experiment_name,
        )

        sp = process_wellplate_spectra.SpectraProcessor(
            folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                            '2022-12-01/interpolator-dataset/')
        sp.nanodrop_lower_cutoff_of_wavelengths = 220
        calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'
        df = versatility_examples.process_plate(sp, dilution_factor=500,
                           plate_folder=f'{data_folder}{experiment_name}2023-10-12_14-51-59_UV-Vis_crude.csv',
                           well_ids=range(5),
                           cut_from=5,
                           cut_to=False,
                           calibrant_shortnames=['methoxychalcone', 'anisaldehyde', 'acetophenone'],
                           calibration_folder=calibration_folder,
                           experiment_name=experiment_name,
                           do_plot=False)

        if SAVE_EXPECTED_OUTPUTS:
            # save df as 'Claisen_WaiShing_results.pickle'
            df.to_pickle(f'Claisen_WaiShing_results.pickle')

        # assert that the df and the one loaded from pickle are equal
        df_loaded = pd.read_pickle(f'Claisen_WaiShing_results.pickle')
        assert_frame_equal(df, df_loaded, check_exact=False, check_like=True)
