from robowski.settings import *
import importlib
import pickle
from distutils import dir_util
import pandas as pd
from pytest import fixture
import os
import robowski.misc_scripts.organize_run_results as organize_run_results


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


def test_load_df_from_run_info(datadir):
    """
    Test uses a fixture that copies all structure from `tests/test_organize_run_results` directory into a temporary
    directory, which is later treated as the data_folder (that is normally in the Dropbox, but not for tests). Test
    loads the run_info from respective locations in the `multicomp-reactions/2023-06-20-run01/` in the temporary
    folder and then checks the results against an expected dataframe loaded from
    `expected_outputs/run_info_no_version.pkl`. Then it does the same for `multicomp-reactions/2023-07-04-run01/` and
    checks against `expected_outputs/run_info_V1.00.pkl` to check for run_info version 1.00

    Parameters
    ----------
    datadir: pytest fixture
        Temporary directory with the same structure as `tests/test_organize_run_results` directory.
    """
    data_folder = ''
    with datadir.as_cwd():
        experiment_name = 'multicomp-reactions/2023-06-20-run01/'
        pd.testing.assert_frame_equal(
            organize_run_results.load_df_from_run_info(experiment_name + 'pipetter_io/run_info.csv'),
            pd.read_pickle('expected_outputs/run_info_no_version.pkl'))

        experiment_name = 'multicomp-reactions/2023-07-04-run01/'
        pd.testing.assert_frame_equal(
            organize_run_results.load_df_from_run_info(experiment_name + 'pipetter_io/run_info.csv'),
            pd.read_pickle('expected_outputs/run_info_v1.00.pkl'))


def test_load_df_from_dilution_info(datadir):
    """
    Test uses a fixture that copies all structure from `tests/test_organize_run_results` directory into a temporary
    directory, which is later treated as the data_folder (that is normally in the Dropbox, but not for tests). Test
    loads the dataframe with information about dilution from respective locations in the
    `multicomp-reactions/2023-06-20-run01/` in the temporary folder and then checks the results against an expected
    dataframe loaded from `expected_outputs/dilution_info.pkl`.

    Parameters
    ----------
    datadir: pytest fixture
        Temporary directory with the same structure as `tests/test_organize_run_results` directory.
    """
    data_folder = ''
    with datadir.as_cwd():
        pd.testing.assert_frame_equal(
            organize_run_results.load_df_from_dilution_info('multicomp-reactions/2023-06-20-run01/'),
            pd.read_pickle('expected_outputs/dilution_info.pkl'))


def test_join_data_from_runs(datadir):
    """
    Test uses a fixture that copies all structure from `tests/test_organize_run_results` directory into a temporary
    directory, which is later treated as the data_folder (that is normally in the Dropbox, but not for tests).
    Test joins data from runs `multicomp-reactions/2023-06-20-run01/`, `multicomp-reactions/2023-06-21-run01/` and
    `multicomp-reactions/2023-06-21-run02/` and checks the results against an expected dataframe loaded from
    `expected_outputs/joined_data_from_runs.pkl`.

    Parameters
    ----------
    datadir: pytest fixture
        Temporary directory with the same structure as `tests/test_organize_run_results` directory.
    """
    data_folder = ''
    with datadir.as_cwd():
        frame1 = organize_run_results.join_data_from_runs(['multicomp-reactions/2023-06-20-run01/',
                                                      'multicomp-reactions/2023-06-21-run01/',
                                                      'multicomp-reactions/2023-06-21-run02/'
                                                     ])
        # frame1.to_pickle('C:/Users/ICanSeeYourPixels/Desktop/temp/joined_data_from_runs.pkl')
        frame2 = pd.read_pickle('expected_outputs/joined_data_from_runs.pkl')
        pd.testing.assert_frame_equal(frame1, frame2)


def test_get_excel_file_version(datadir):
    with datadir.as_cwd():
        assert organize_run_results.get_excel_file_version(
            'simple-reactions/2023-08-21-run01/2023-08-21-run01.xlsx') == 1
        assert organize_run_results.get_excel_file_version(
            'simple-reactions/2023-07-16-run01/2023-07-16-run01.xlsx') == 0


def test_locate_condition_by_operation_datetime_and_plate_id(datadir):
    with datadir.as_cwd():
        df_conditions = pd.read_excel('simple-reactions/2099-99-99-run01/2099-99-99-run01.xlsx',
                                      sheet_name='reactions_with_run_info')
        d = pickle.load(open('expected_outputs/locate_condition_by_operation_datetime_and_plate_id.pkl', 'rb'))
    for timestamp in d.keys():
        res = organize_run_results.locate_condition_by_operation_datetime_and_plate_id(
            timestamp=timestamp,
            plate_id=51,
            dataframe_with_conditions=df_conditions,
            column_name_for_plate_id='plate_barcodes_for_dilution',
            column_name_for_timestamp='timestamp_dilution',
            column_name_for_container_id='container_id'
        )
        assert res == tuple(d[timestamp])

    res = organize_run_results.locate_condition_by_operation_datetime_and_plate_id(
        timestamp=timestamp,
        plate_id=152,
        dataframe_with_conditions=df_conditions,
        column_name_for_plate_id='plate_barcodes_for_dilution',
        column_name_for_timestamp='timestamp_dilution',
        column_name_for_container_id='container_id'
    )
    assert len(res) == 0