from datetime import datetime
import logging
import os
import pickle
import re
import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from xlrd.biffh import XLRDError

logging.basicConfig(level=logging.INFO)

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'


def remove_one_outlier_and_average_rest(x, verbose=False):
    # x is a numpy array
    sorted_list = np.sort(x)
    diff = np.diff(sorted_list)
    if np.abs(diff[0]) > np.abs(diff[-1]):
        if verbose:
            print('Removing the first value.')
        result = np.mean(sorted_list[1:])
    else:
        if verbose:
            print('Removing the last value.')
        result = np.mean(sorted_list[:-1])
    return result


def round_to_nearest(df_new, df_reference, column_names):
    """
    Round the values in the new dataframe to the nearest values in the reference dataframe.
    For example, if column in reference dataframe has values [0.1, 0.2, 0.3] and the new dataframe has values
    [5, 8, 11, 0.099999999999999999998, 11, 0.2000000000000000001, 0.35], the new dataframe values will be rounded to
    [5, 8, 11, 0.1, 0.2, 0.3. 11. 0.2, 0.35]
    and returned.

    Parameters
    ----------
    df_new: pd.DataFrame
        The dataframe with the new values that may be slightly different from the reference values due to rounding
        errors.

    df_reference: pd.DataFrame
        The dataframe with the reference values.

    column_names: list of str
        Only these columns will be changed in the new dataframe. Columns of same names are compared between
        the dataframes.

    Returns
    -------
    df_new: pd.DataFrame
        The dataframe with the new values rounded to the nearest values in the reference dataframe.
    """
    for column_name in column_names:
        new_values = df_new[column_name].to_numpy()
        unique_values_from_reference_df = df_reference[column_name].unique()
        for i, new_value in enumerate(new_values):
            for unique_value_from_reference in unique_values_from_reference_df:
                if np.isclose(new_value, unique_value_from_reference):
                    new_values[i] = unique_value_from_reference
        df_new[column_name] = new_values
    return df_new


def join_data_from_runs(run_names, round_on_columns=('ic001', 'am001', 'ald001', 'ptsa'), subfolder_containing_csv='results'):
    """
    Loads the `results/product_concentration.csv` from multiple runs and joins them into one dataframe.

    Parameters
    ----------
    run_names: list
        The names of the runs whose results will be loaded and joined. These names are the names of the subfolders
        of the data folder.

    round_on_columns
        The columns that will be rounded to the nearest values in the reference dataframe. This is done to avoid
        rounding errors when joining the dataframes.

    Returns
    -------
    df_result: pd.DataFrame
        The dataframe with the joined data. The sequence of the rows is the same as the sequence of the runs.
    """
    df_result = pd.read_csv(data_folder + run_names[0] + f'{subfolder_containing_csv}/product_concentration.csv')
    df_result['run_name'] = run_names[0]
    logging.info(f'Load the run, number of rows: {df_result.shape[0]}')
    if 'Unnamed: 0' in df_result.columns:
        df_result.drop('Unnamed: 0', inplace=True, axis=1)
    for run_name in run_names[1:]:
        df_temporary = pd.read_csv(data_folder + run_name + f'{subfolder_containing_csv}/product_concentration.csv')
        logging.info(f'Loaded the run, number of rows: {df_temporary.shape[0]}')
        df_temporary['run_name'] = run_name
        if 'Unnamed: 0' in df_temporary.columns:
            df_temporary.drop('Unnamed: 0', inplace=True, axis=1)
        if round_on_columns is not None:
            df_temporary = round_to_nearest(df_temporary, df_result, round_on_columns)
        df_result = df_result.append(df_temporary, ignore_index=True)

    logging.info(f'Number of rows in the joined dataframe: {df_result.shape[0]}')
    return df_result


def load_df_from_run_info(path_to_run_info_file):
    """
    Load the `pipetter_io/run_info.csv` file into a dataframe.
    Parameters
    ----------
    path_to_run_info_file: str
        The name of the `pipetter_io/run_info.csv` file.

    Returns
    -------
    df_pipetter: pd.DataFrame
        The dataframe with the run info.

    """
    with open(path_to_run_info_file, 'r') as f:
        first_line = f.readline()
    # take care of the version of run_info:
    if first_line.startswith('#version'):
        if first_line.startswith('#version: 1.00'):
            df_pipetter = pd.read_csv(path_to_run_info_file, delimiter=',', header=0, index_col= False,
                                      names=['plate_code', 'experiment_name', 'start_time_unix',
                                             'start_time_string', 'finish_time_unix', 'finish_time_string', 'note'])
        else:
            raise ValueError(
                f'Unknown version of run_info.csv: {first_line[:-1]}. Good luck coding your own damn loader.')
    else:
        df_pipetter = pd.read_csv(path_to_run_info_file, delimiter=', ', header=None,
                                  names=['plate_code', 'experiment_name', 'start_time_unix',
                                         'start_time_string', 'finish_time_unix', 'finish_time_string', 'note'])
    return df_pipetter


def load_df_from_dilution_info(experiment_name):
    # check if 'dilution' subfolder is present in the experiment folder
    path_to_dilution_info_file = data_folder + experiment_name + 'dilution/dilution_info.csv'
    if not os.path.exists(path_to_dilution_info_file):
        raise ValueError(
            f'The dilution info file {path_to_dilution_info_file} does not exist. '
            f'Make it or code your own damn loader.')
    else:
        with open(path_to_dilution_info_file, 'r') as f:
            first_line = f.readline()
        if first_line.startswith('#version'):
            if first_line.startswith('#version: 1.00'):
                df_dilution = pd.read_csv(path_to_dilution_info_file, header=1)
                # make sure that there are only two columns: 'from' and 'to'
                assert len(df_dilution.columns) == 2, 'The dilution_info.csv file has more than two columns.'
                assert df_dilution.columns[0] == 'from', 'The first column of dilution_info.csv is not "from".'
                assert df_dilution.columns[1] == 'to', 'The second column of dilution_info.csv is not "to".'
            else:
                raise ValueError(
                    f'Unknown version of dilution_info.csv: {first_line[:-1]}. Good luck coding your own damn loader.')
        else:
            raise ValueError(f'No version version in first line of dilution_info.csv: {first_line[:-1]}. '
                             f'Good luck coding your own damn loader.')
    return df_dilution


def check_run_data_consistency(list_of_runs):
    """
    Check if the data in the completed runs is consistent.

    Parameters
    ----------
    list_of_runs: list of str
        The names of the runs whose consistency will be checked.

    Returns
    -------
    bool
        True if the data is consistent, False otherwise.
    """
    global craic_folder
    global data_folder

    logging.info('Checking consistency of the CRAIC database.')
    df_craic = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
    # iterate over rows of df_craic
    for index, row in df_craic.iterrows():
        if row['timestamp'] < 1686989327:
            continue
        # check if the existing folder corresponds to the experiment_name
        assert row['exp_name'] in row['folder'], \
            'The exp_name is not in folder name.'
        # check that the folder actually exists on the disk
        assert os.path.exists(craic_folder + row['folder']), \
            f'Folder {row["folder"]} is in database but does not exist on drive: {craic_folder + row["folder"]}'

    for experiment_name in list_of_runs:
        logging.info(f'Checking consistency of {experiment_name}')
        run_name = experiment_name.split('/')[1]
        run_type = experiment_name.split('/')[0]

        # make sure that there exists file 'outC.csv' in the experiment folder 'outVandC'
        # assert os.path.exists(data_folder + experiment_name + 'outVandC/outC.csv'), 'File outC.csv does not exist.'
        # assert os.path.exists(data_folder + experiment_name + 'outVandC/outV.csv'), 'File outV.csv does not exist.'


        # open dilution_info.csv as dataframe
        df_dilution = load_df_from_dilution_info(experiment_name)

        # open run_info.csv as dataframe
        df_pipetter = load_df_from_run_info(data_folder + experiment_name + 'pipetter_io/run_info.csv')

        # open the Excel file with the volumes, use first sheet
        df_structure = pd.read_excel(data_folder + experiment_name + f'{run_name}.xlsx', sheet_name=0)

        df_craic = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
        exp_names_craic = [f'{run_type.replace("-", "_")}_{run_name}']
        exp_names_craic_dil = [f'{run_type.replace("-", "_")}_{run_name}_dil']
        if df_craic['exp_name'].isin(exp_names_craic_dil).any():
            df_craic = df_craic.loc[df_craic['exp_name'].isin(exp_names_craic_dil)].copy().reset_index()
        else:
            df_craic = df_craic.loc[df_craic['exp_name'].isin(exp_names_craic)].copy().reset_index()
        assert len(
            df_craic) > 0, f"No plates found in CRAIC database for {experiment_name}. Good luck learning to type."
        # logging.info(f'Found {len(df_craic)} plates in CRAIC database for {experiment_name}: '
        #              f'plate IDs: {df_craic["plate_id"].values}')

        assert len(df_dilution) == len(df_pipetter) == len(df_structure) / 54 == len(df_craic), \
            'numbers of rows in the dataframes are inconsistent'

        assert len(df_structure) % 54 == 0, 'number of conditions is not divisible by 54'

        # iterate over the plates of df_pipetter
        for row_id, row in enumerate(df_pipetter.itertuples()):
            # verify that the logic of the data is correct
            reaction_plate = row.plate_code
            assert df_dilution.iloc[row_id]['from'] == reaction_plate, 'Time sequence of plates is wrong.'
            plate_for_dilution = df_dilution.iloc[row_id]['to']
            assert df_craic.iloc[row_id]['plate_id'] == plate_for_dilution, 'Time sequence of plates is wrong.'

    logging.info('Consistency is OK.')
    return True

def organize_run_structure(experiment_name, version='1.00'):
    if version == '1.00' or version == '':
        organize_run_structure_v1_00(experiment_name)
    elif version == '2.00':
        organize_run_structure_v2_00(experiment_name)
    elif version == '3.00':
        return organize_run_structure_v3_00(experiment_name)
    else:
        raise ValueError(f'Version {version} of organize_run_structure is not implemented.')


def organize_run_structure_v1_00(experiment_name):
    """
    Automatically combine `run_info.csv`, `dilution_info.csv` and CRAIC plate database into
    a unified "run structure" table saved into `results/run_structure.csv`.

    The output table indicates for each condition the vial_id, the reaction plate id, the id of plate it was
    diluted into, and the CRAIC folder name that contains the spectra for this plate.

    Parameters
    ----------
    experiment_name: str
        The name of the experiment, e.g. `multicomp-reactions/2023-06-26-run01/`
    """
    version = '1.00'
    global craic_folder
    global data_folder

    run_name = experiment_name.split('/')[1]
    run_type = experiment_name.split('/')[0]

    # if there is not 'results' folder in the run folder, create it
    target_folder = data_folder + experiment_name + 'results'
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # open dilution_info.csv as dataframe
    df_dilution = load_df_from_dilution_info(experiment_name)

    # open run_info.csv as dataframe
    df_pipetter = load_df_from_run_info(data_folder + experiment_name + 'pipetter_io/run_info.csv')

    # open the Excel file with the volumes, use first sheet
    df_structure = pd.read_excel(data_folder + experiment_name + f'{run_name}.xlsx', sheet_name=0)
    df_structure['vial_id'] = 0
    df_structure['reaction_plate_id'] = 0
    df_structure['diluted_plate_id'] = 0
    df_structure['craic_folder'] = ''
    df_structure['is_outlier'] = 0

    exp_names_craic = [f'multicomp_reactions_{run_name}']
    df_craic = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
    exp_names_craic = [f'{run_type.replace("-", "_")}_{run_name}']
    exp_names_craic_dil = [f'{run_type.replace("-", "_")}_{run_name}_dil']
    if df_craic['exp_name'].isin(exp_names_craic_dil).any():
        df_craic = df_craic.loc[df_craic['exp_name'].isin(exp_names_craic_dil)].copy().reset_index()
    else:
        df_craic = df_craic.loc[df_craic['exp_name'].isin(exp_names_craic)].copy().reset_index()
    assert len(
        df_craic) > 0, f"No plates found in CRAIC database for {experiment_name}. Good luck learning to type."
    # logging.info(f'Found {len(df_craic)} plates in CRAIC database for {experiment_name}: '
    #              f'plate IDs: {df_craic["plate_id"].values}')

    # verify than numbers of rows in all the dataframes are the same
    assert len(df_dilution) == len(df_pipetter) == len(df_structure) / 54 == len(df_craic)

    # verify that number of rows is divisible by 54
    assert len(df_structure) % 54 == 0

    # iterate over the plates of df_pipetter
    for row_id, row in enumerate(df_pipetter.itertuples()):
        # verify that the logic of the data is correct
        reaction_plate = row.plate_code
        assert df_dilution.iloc[row_id]['from'] == reaction_plate
        plate_for_dilution = df_dilution.iloc[row_id]['to']
        assert df_craic.iloc[row_id]['plate_id'] == plate_for_dilution
        craic_folder_here = df_craic.iloc[row_id]['folder']

        # populate the structure dataframe with the plate codes and the craic folder
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'reaction_plate_id'] = reaction_plate
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'diluted_plate_id'] = plate_for_dilution
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'craic_folder'] = craic_folder_here
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'vial_id'] = tuple(range(54))

    # save the dataframe as csv
    structure_csv_filepath = data_folder + experiment_name + 'results/run_structure.csv'
    with open(structure_csv_filepath, 'w') as f:
        f.write(f'version: {version}\n')
    df_structure.to_csv(structure_csv_filepath, index=False, mode='a')


def organize_run_structure_v2_00(experiment_name):
    """
     Automatically combine `run_info.csv`, `dilution_info.csv` and CRAIC plate database into
     a unified "run structure" table saved into `results/run_structure.csv`.

     The output table indicates for each condition the vial_id, the reaction plate id, the id of plate it was
     diluted into, and the CRAIC folder name that contains the spectra for this plate.

     Parameters
     ----------
     experiment_name: str
         The name of the experiment, e.g. `multicomp-reactions/2023-06-26-run01/`
     """
    version = '2.00'
    global craic_folder
    global data_folder

    run_name = experiment_name.split('/')[1]
    run_type = experiment_name.split('/')[0]

    # if there is not 'results' folder in the run folder, create it
    target_folder = data_folder + experiment_name + 'results'
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # open dilution_info.csv as dataframe
    df_dilution = load_df_from_dilution_info(experiment_name)

    # open run_info.csv as dataframe
    df_pipetter = load_df_from_run_info(data_folder + experiment_name + 'pipetter_io/run_info.csv')

    # open the Excel file with the volumes, use first sheet
    df_structure = pd.read_excel(data_folder + experiment_name + f'{run_name}.xlsx', sheet_name=0)
    df_structure['vial_id'] = 0
    df_structure['reaction_plate_id'] = 0
    df_structure['diluted_plate_id'] = 0
    df_structure['craic_folder'] = ''
    df_structure['is_outlier'] = 0

    exp_names_craic = [f'multicomp_reactions_{run_name}']
    df_craic = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
    exp_names_craic = [f'{run_type.replace("-", "_")}_{run_name}']
    exp_names_craic_dil = [f'{run_type.replace("-", "_")}_{run_name}_dil']
    if df_craic['exp_name'].isin(exp_names_craic_dil).any():
        df_craic = df_craic.loc[df_craic['exp_name'].isin(exp_names_craic_dil)].copy().reset_index()
    else:
        df_craic = df_craic.loc[df_craic['exp_name'].isin(exp_names_craic)].copy().reset_index()
    assert len(
        df_craic) > 0, f"No plates found in CRAIC database for {experiment_name}. Good luck learning to type."
    # logging.info(f'Found {len(df_craic)} plates in CRAIC database for {experiment_name}: '
    #              f'plate IDs: {df_craic["plate_id"].values}')

    # verify than numbers of rows in all the dataframes are the same
    assert len(df_dilution) == len(df_pipetter) == len(df_structure) / 54 == len(df_craic)

    # verify that number of rows is divisible by 54
    assert len(df_structure) % 54 == 0

    # iterate over the plates of df_pipetter
    for row_id, row in enumerate(df_pipetter.itertuples()):
        # verify that the logic of the data is correct
        reaction_plate = row.plate_code
        assert df_dilution.iloc[row_id]['from'] == reaction_plate
        plate_for_dilution = df_dilution.iloc[row_id]['to']
        assert df_craic.iloc[row_id]['plate_id'] == plate_for_dilution
        craic_folder_here = df_craic.iloc[row_id]['folder']

        # populate the structure dataframe with the plate codes and the craic folder
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'reaction_plate_id'] = reaction_plate
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'diluted_plate_id'] = plate_for_dilution
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'craic_folder'] = craic_folder_here
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'vial_id'] = tuple(range(54))

    # if there is a file 'outVandC/stock_solutions.xlsx', load it and use to determine concentrations

    # load excel table with stock solution compositions
    df_stock = pd.read_excel(data_folder + experiment_name + f'outVandC/stock_solutions.xlsx', sheet_name=0)
    # iterate over rows in df_excel
    for substance in df_stock.columns[1:]:
        df_structure[f'c#{substance}'] = 0

    for index, row in df_structure.iterrows():
        # iterate over columns in df_structure
        substances = df_stock.columns[1:]
        net_volume = 0
        for stock_solution_name in df_structure.columns:
            # if this is a column with a stock solution
            if stock_solution_name in df_stock['stock_solution'].to_list():
                # get the stock solution concentration
                stock_concentration = df_stock.loc[df_stock['stock_solution'] == stock_solution_name]
                # get the volume of this stock solution used by the pipetter
                volume = row[stock_solution_name]
                for substance in df_stock.columns[1:]:
                    this_stock_solution_concentration = stock_concentration[substance].iloc[0]
                    df_structure.loc[index, f'c#{substance}'] = df_structure.loc[
                                                                index, f'c#{substance}'] + volume * this_stock_solution_concentration
                net_volume += volume
        for substance in df_stock.columns[1:]:
            df_structure.loc[index, f'c#{substance}'] = df_structure.loc[index, f'c#{substance}'] / net_volume

    # save the dataframe as csv
    structure_csv_filepath = data_folder + experiment_name + 'results/run_structure.csv'
    with open(structure_csv_filepath, 'w') as f:
        f.write(f'version: {version}\n')
    df_structure.to_csv(structure_csv_filepath, index=False, mode='a')


def organize_run_structure_v2_01(experiment_name):
    """
     Automatically combine `run_info.csv`, `dilution_info.csv` and CRAIC plate database into
     a unified "run structure" table saved into `results/run_structure.csv`.

     The output table indicates for each condition the vial_id, the reaction plate id, the id of plate it was
     diluted into, and the CRAIC folder name that contains the spectra for this plate.

     This version also adds the nondiluted spectra folders into the final table

     Parameters
     ----------
     experiment_name: str
         The name of the experiment, e.g. `multicomp-reactions/2023-06-26-run01/`
     """
    version = '2.00'
    global craic_folder
    global data_folder

    run_name = experiment_name.split('/')[1]
    run_type = experiment_name.split('/')[0]

    # if there is not 'results' folder in the run folder, create it
    target_folder = data_folder + experiment_name + 'results'
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # open dilution_info.csv as dataframe
    df_dilution = load_df_from_dilution_info(experiment_name)

    # open run_info.csv as dataframe
    df_pipetter = load_df_from_run_info(data_folder + experiment_name + 'pipetter_io/run_info.csv')

    # open the Excel file with the volumes, use first sheet
    df_structure = pd.read_excel(data_folder + experiment_name + f'{run_name}.xlsx', sheet_name=0)
    df_structure['vial_id'] = 0
    df_structure['reaction_plate_id'] = 0
    df_structure['diluted_plate_id'] = 0
    df_structure['craic_folder'] = ''
    df_structure['is_outlier'] = 0

    exp_names_craic = [f'multicomp_reactions_{run_name}']
    df_craic = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
    exp_names_craic = [f'{run_type.replace("-", "_")}_{run_name}']
    exp_names_craic_dil = [f'{run_type.replace("-", "_")}_{run_name}_dil']
    if df_craic['exp_name'].isin(exp_names_craic_dil).any():
        df_craic = df_craic.loc[df_craic['exp_name'].isin(exp_names_craic_dil)].copy().reset_index()
    else:
        df_craic = df_craic.loc[df_craic['exp_name'].isin(exp_names_craic)].copy().reset_index()
    assert len(
        df_craic) > 0, f"No plates found in CRAIC database for {experiment_name}. Good luck learning to type."

    df_craic2 = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
    df_craic_undil = df_craic2.loc[df_craic2['exp_name'].isin(exp_names_craic)].copy().reset_index()

    # logging.info(f'Found {len(df_craic)} plates in CRAIC database for {experiment_name}: '
    #              f'plate IDs: {df_craic["plate_id"].values}')

    # verify than numbers of rows in all the dataframes are the same
    assert len(df_dilution) == len(df_pipetter) == len(df_structure) / 54 == len(df_craic)

    # verify that number of rows is divisible by 54
    assert len(df_structure) % 54 == 0

    # iterate over the plates of df_pipetter
    for row_id, row in enumerate(df_pipetter.itertuples()):
        # verify that the logic of the data is correct
        reaction_plate = row.plate_code
        assert df_dilution.iloc[row_id]['from'] == reaction_plate
        plate_for_dilution = df_dilution.iloc[row_id]['to']
        assert df_craic.iloc[row_id]['plate_id'] == plate_for_dilution
        craic_folder_here = df_craic.iloc[row_id]['folder']
        assert df_craic_undil.iloc[row_id]['plate_id'] == reaction_plate
        undil_folder_here = df_craic_undil.iloc[row_id]['folder']

        # populate the structure dataframe with the plate codes and the craic folder
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'reaction_plate_id'] = reaction_plate
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'diluted_plate_id'] = plate_for_dilution
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'craic_folder'] = craic_folder_here
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'craic_folder_undil'] = undil_folder_here
        df_structure.at[row_id * 54:(row_id + 1) * 54 - 1, 'vial_id'] = tuple(range(54))

    # if there is a file 'outVandC/stock_solutions.xlsx', load it and use to determine concentrations

    # load excel table with stock solution compositions
    df_stock = pd.read_excel(data_folder + experiment_name + f'outVandC/stock_solutions.xlsx', sheet_name=0)
    # iterate over rows in df_excel
    for substance in df_stock.columns[1:]:
        df_structure[f'c#{substance}'] = 0

    for index, row in df_structure.iterrows():
        # iterate over columns in df_structure
        substances = df_stock.columns[1:]
        net_volume = 0
        for stock_solution_name in df_structure.columns:
            # if this is a column with a stock solution
            if stock_solution_name in df_stock['stock_solution'].to_list():
                # get the stock solution concentration
                stock_concentration = df_stock.loc[df_stock['stock_solution'] == stock_solution_name]
                # get the volume of this stock solution used by the pipetter
                volume = row[stock_solution_name]
                for substance in df_stock.columns[1:]:
                    this_stock_solution_concentration = stock_concentration[substance].iloc[0]
                    df_structure.loc[index, f'c#{substance}'] = df_structure.loc[
                                                                index, f'c#{substance}'] + volume * this_stock_solution_concentration
                net_volume += volume
        for substance in df_stock.columns[1:]:
            df_structure.loc[index, f'c#{substance}'] = df_structure.loc[index, f'c#{substance}'] / net_volume

    # save the dataframe as csv
    structure_csv_filepath = data_folder + experiment_name + 'results/run_structure.csv'
    with open(structure_csv_filepath, 'w') as f:
        f.write(f'version: {version}\n')
    df_structure.to_csv(structure_csv_filepath, index=False, mode='a')


def organize_run_structure_v3_00(experiment_name):
    """
    This function is for nanodrop spectrometer version. It takes the run_info.csv and dilution_info.csv files and
    creates a run_structure.csv file that contains all the information about the run.

     Parameters
     ----------
     experiment_name: str
         The name of the experiment, e.g. `multicomp-reactions/2023-06-26-run01/`
     """
    compatible_versions_of_excel_file = [1.00]
    output_version = '3.00'
    global data_folder

    run_name = experiment_name.split('/')[1]
    run_type = experiment_name.split('/')[0]

    # if there is not 'results' folder in the run folder, create it
    target_folder = data_folder + experiment_name + 'results'
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # open the Excel file with the volumes, use first sheet
    path_to_excel_file = data_folder + experiment_name + f'{run_name}.xlsx'
    assert get_excel_file_version(path_to_excel_file) in compatible_versions_of_excel_file, \
        f"Excel file version is not compatible with this function. " \
        f"Compatible versions are: {compatible_versions_of_excel_file}"
    df_structure = pd.read_excel(path_to_excel_file, sheet_name=0)
    # drop columns that contain 'Unnamed'
    df_structure = df_structure.loc[:, ~df_structure.columns.str.contains('^Unnamed')]
    df_structure['is_outlier'] = 0

    # if there is a file 'outVandC/stock_solutions.xlsx', load it and use to determine concentrations

    # load excel table with stock solution compositions
    df_stock = pd.read_excel(path_to_excel_file, sheet_name='stock_solutions_concentration', skiprows=[0])
    # iterate over rows in df_excel
    for substance in df_stock.columns[1:]:
        df_structure[f'c#{substance}'] = 0

    names_of_stock_solutions_used = [stock_solution_name.replace('vol#','') for stock_solution_name
                                     in df_structure.columns if 'vol#' in stock_solution_name]

    for stock_solution_name in names_of_stock_solutions_used:
        if stock_solution_name not in df_stock['stock_solution'].to_list():
            raise ValueError(f"Stock solution {stock_solution_name} is used in conditions list but is not in the stock solutions table.")
    for index, row in df_structure.iterrows():
        total_volume = 0
        for stock_solution_name in names_of_stock_solutions_used:
            # if this is a column with a stock solution
            if stock_solution_name in df_stock['stock_solution'].to_list():
                concentrations_of_substances_in_this_stock_solution = df_stock.loc[df_stock['stock_solution'] == stock_solution_name]
                pipetted_volume_of_this_stock_solution = row[f'vol#{stock_solution_name}']
                for substance in df_stock.columns[1:]:
                    concentration_of_this_substance = concentrations_of_substances_in_this_stock_solution[substance].iloc[0]
                    df_structure.loc[index, f'c#{substance}'] += pipetted_volume_of_this_stock_solution * concentration_of_this_substance
                total_volume += pipetted_volume_of_this_stock_solution
        for substance in df_stock.columns[1:]:
            df_structure.loc[index, f'c#{substance}'] = df_structure.loc[index, f'c#{substance}'] / total_volume

    # assign the nanodrop spectra to rows of the df_structure based on plate id and well id
    nanodrop_spectra_folder = data_folder + experiment_name + 'nanodrop_spectra/'
    # Use regular expressions to get the plate id and unix timestamp for each file by parsing the filename. For example,
    # filename '2023-08-23_23-52-13_plate_51.csv' means for 2023-08-23 date, 23-52-13 hour-minute-second, and plate 51
    ## Old regexp:
    # regexp_pattern = re.compile(r'(?P<timestamp_string>\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})_plate_(?P<plate_id>\d+).csv')

    # Changing regexp because Yankai is very creative
    pattern = re.compile(r'(?P<timestamp_string>\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})(?:.*?)_plate_(?P<plate_id>\d+).csv',
                         re.DOTALL)
    for filename in os.listdir(nanodrop_spectra_folder):
        if not filename.endswith('.csv'):
            continue
        match = pattern.match(filename)
        if not match:
            logging.warning(f'Could not parse nanodrop spectrum filename {filename} for timestamp and id')
            continue
        plate_id = int(match.group('plate_id'))
        unix_timestamp = datetime.strptime(match.group('timestamp_string'), '%Y-%m-%d_%H-%M-%S').timestamp()
        logging.info(
            f'Found nanodrop spectrum file {filename} with plate id {plate_id} '
            f'and unix timestamp {unix_timestamp}')

        # locate the conditions in the df_structure for each file using the locate_condition_by_operation... function
        # and write the filepath into the df_structure column 'nanodrop_filepath"
        indices = locate_condition_by_operation_datetime_and_plate_id(unix_timestamp, plate_id, df_structure,
                                                                      column_name_for_plate_id='plate_barcodes_for_dilution',
                                                                      column_name_for_timestamp='timestamp_dilution',
                                                                      column_name_for_container_id='container_id')
        if len(indices) == 0:
            logging.warning(f'Could not find conditions in run structure for nanodrop spectrum file {filename}')
            continue
        elif len(indices) > 0:
            logging.info(f'Found {len(indices)} conditions in run structure for nanodrop spectrum file {filename}')
            df_structure.loc[indices, 'nanodrop_filepath'] = experiment_name + 'nanodrop_spectra/' + filename
            df_structure.loc[indices, 'timestamp_nanodrop'] = unix_timestamp

    # save the dataframe as csv
    structure_csv_filepath = data_folder + experiment_name + 'results/run_structure.csv'
    with open(structure_csv_filepath, 'w') as f:
        f.write(f'version: {output_version}\n')
    df_structure.to_csv(structure_csv_filepath, index=False, mode='a')
    return df_structure


def load_run_structure(experiment_name):
    """
    Load the run structure from the csv file into a dataframe.

    Returns
    -------
    df: pandas.DataFrame
        The run structure dataframe.
    """
    global data_folder
    path_to_file = data_folder + experiment_name + 'results/run_structure.csv'
    # check if the first line of text file contains 'version'
    with open(path_to_file, 'r') as f:
        first_line = f.readline()
    if 'version' in first_line:
        version = first_line.split(' ')[1].strip()      # get the version number
        df = pd.read_csv(path_to_file, skiprows=1)
    else:
        df = pd.read_csv(path_to_file)
    # if there are columns containing string "Unnamed", remove them
    for column in df.columns:
        if 'Unnamed' in column:
            df.drop(column, axis=1, inplace=True)
    return df


def outV_to_outC_by_lookup(experiment_name, lookup_run):
    df_lookup_C = pd.read_csv(data_folder + lookup_run + 'outVandC/outC.csv')
    df_lookup_C.drop('Unnamed: 0', inplace=True, axis=1)
    df_lookup_V = pd.read_csv(data_folder + lookup_run + 'outVandC/outV.csv')
    df_lookup_V.drop('Unnamed: 0', inplace=True, axis=1)

    # load Excel file from experiment_name run into dataframe
    run_name = experiment_name.split('/')[1]
    df_excel = pd.read_excel(data_folder + experiment_name + f'{run_name}.xlsx', sheet_name=0)

    df_outC = pd.DataFrame(dtype=object).reindex_like(df_lookup_C)[0:0]
    # iterate over df_excel, look for row of df_lookup_C with index equal to 'reactions' and populate df_outC
    for row_id, row in df_excel.iterrows():
        # Zero-filled rows in Excel are used as padding for matching to 54 conditions on each plate).
        # In this case, fill row of concentrations in df_outC with zeros.
        if (row.to_numpy()[1:] == 0).all():
            df_outC.loc[row_id] = df_lookup_C.loc[0] * 0
        else:
            # find the row with same id in df_lookup_C and write to df_outC
            df_outC.loc[row_id] = df_lookup_C.loc[row['reactions']]
            # Lake sure that the volumes in this row are the same as in the row in lookup_V, unless the row is all zeros
            same_row_in_lookup_V = df_lookup_V.loc[row['reactions']].to_numpy()
            assert np.isclose(same_row_in_lookup_V, row.to_numpy()[1:]).all()

    # save outC to /outVandC/outC.csv
    # if there is not 'outVansC' folder in the run folder, create it
    target_folder = data_folder + experiment_name + 'outVandC'
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)
    df_outC.to_csv(data_folder + experiment_name + 'outVandC/outC.csv')


def merge_repeated_outliers(original_run, outlier_runs,
                            output_csv_filename='product_concentration_after_substituting_outliers.csv'):
    # find rows in outliers run that have same conditions as in the original run
    df_results_original = pd.read_csv(data_folder + original_run + 'results/product_concentration.csv')
    if 'Unnamed: 0' in df_results_original.columns:
        df_results_original.drop('Unnamed: 0', inplace=True, axis=1)
    df_results_outliers = join_data_from_runs(outlier_runs, round_on_columns=tuple())
    if 'Unnamed: 0' in df_results_outliers.columns:
        df_results_outliers.drop('Unnamed: 0', inplace=True, axis=1)
    df_results_concatenated = pd.concat([df_results_original, df_results_outliers], ignore_index=True)
    # remove known outliers
    df_results_concatenated = df_results_concatenated[df_results_concatenated['is_outlier'] == 0]
    df_results_concatenated.reset_index(inplace=True, drop=True)
    substrate_names = ['ic001', 'ald001', 'am001', 'ptsa']
    substrate_vectors_concatenated = df_results_concatenated[substrate_names].to_numpy()
    distances = pairwise_distances(substrate_vectors_concatenated,
                                   substrate_vectors_concatenated, metric='euclidean')
    distance_threshold = 1e-7
    mindistance_by_first_index = (distances + np.eye(distances.shape[0])).min(axis=1)
    yield_lists = []
    indices_that_were_processed = []
    for first_index in range(substrate_vectors_concatenated.shape[0]):
        if mindistance_by_first_index[first_index] < distance_threshold:
            indices_that_match = np.where(distances[first_index, :] < distance_threshold)[0]
            all_concentrations_match = True
            for another_index in indices_that_match[1:]:
                for substance in substrate_names:
                    if not np.isclose(df_results_concatenated.iloc[indices_that_match[0]][substance],
                                      df_results_concatenated.iloc[another_index][substance]):
                        all_concentrations_match = False
            if not all_concentrations_match:
                print(f'Not all concentrations match in the following rows:{indices_that_match}')
                print(df_results_concatenated.iloc[indices_that_match])
                continue

            if any([index_here in indices_that_were_processed for index_here in indices_that_match]):
                print('Some of the indices were already processed. Skipping.')
                continue

            indices_that_were_processed.extend(indices_that_match.tolist())

            yields_here = df_results_concatenated.iloc[indices_that_match]['yield'].to_numpy()
            yield_lists.append(yields_here)
            mean_yield_without_outlier = remove_one_outlier_and_average_rest(yields_here)
            df_results_concatenated.at[indices_that_match[0], 'yield'] = mean_yield_without_outlier
            df_results_concatenated.at[indices_that_match[1:], 'is_outlier'] = 1
    df_output = df_results_concatenated[df_results_concatenated['is_outlier'] == 0]
    df_output.to_csv(data_folder + original_run + 'results/' + output_csv_filename, index=False)
    return df_output, yield_lists


def get_excel_file_version(excel_file_path):
    """
    Get the version of the data structure used in excel file of experimental run.
    Version value is stored in the first row of first column of the sheet named 'version'
    If the file does not have a version, return 0.

    Parameters
    ----------
    excel_file_path: str
        Path to the excel file

    Returns
    -------
    version: float

    """
    try:
        version_df = pd.read_excel(excel_file_path, sheet_name='version')
    except XLRDError:
        version = 0
    else:
        # get name of the first column
        version = float(version_df.columns[0])

    return version


def locate_condition_by_operation_datetime_and_plate_id(timestamp, plate_id, dataframe_with_conditions,
                                                        column_name_for_plate_id='plate_barcodes_for_dilution',
                                                        column_name_for_timestamp='timestamp_dilution',
                                                        column_name_for_container_id='container_id') -> tuple:
    """
    Given the dataframe of conditions, the time of operation and the barcode of the plate, find the condition that
    correspond to each container on which the operation was performed. For example, if the operation is measurement
    of spectra of containers in this plate, then it will find rows of conditions that correspond to each measured
    container (vial, well) of the plate.
    The logic of this method is: for each container of a given plate (by plate barcode),
    the method finds the latest condition whose timestamp precedes the input timestamp.

    Returns empty list if all the timestamps in the dataframe are later than the input timestamp.
    Returns empty list if the plate barcode is not found in the dataframe.

    Parameters
    ----------
    timestamp: int
        Timestamp of the operation (UNIX time)
    plate_id: int or str
        Barcode of the plate. This is normally an integer.
    dataframe_with_conditions: pd.DataFrame
        Dataframe with conditions. This dataframe should have columns with names specified in the next three arguments.
    column_name_for_plate_id: str
        Name of the column in dataframe_with_conditions that contains plate barcodes.
    column_name_for_timestamp: str
        Name of the column in dataframe_with_conditions that contains UNIXTIME timestamps of operations that are
        logically assumed to precede the current operation.
    column_name_for_container_id: str
        Name of the column in dataframe_with_conditions that contains container ids (vial or well ids).

    Returns
    -------
    indices_of_found_conditions: list
        List of dataframe indices of conditions that correspond to each container of the plate that was operated on.
    """
    if plate_id in dataframe_with_conditions[column_name_for_plate_id].unique():
        # get all rows for this particular plate
        df = dataframe_with_conditions[(dataframe_with_conditions[column_name_for_plate_id] == plate_id)].copy()
    else:
        logging.warning(f'Plate barcode {plate_id} is not found in the dataframe. Dataframe has plates '
                        f'with barcodes {dataframe_with_conditions[column_name_for_plate_id].unique()}')
        return []

    # make sure that column_name_for_timestamp column is of integer type
    df[column_name_for_timestamp] = df[column_name_for_timestamp].astype(int)

    # iterate over containers in this plate and look at all the sorted timestamps for each container.
    # If the input timestamp is between two timestamps, return the condition as the sought-after condition.
    indices_of_found_conditions = []
    for container_id in df[column_name_for_container_id].unique():
        # get all rows for this particular container
        df_container = df[df[column_name_for_container_id] == container_id].copy()
        # sort these rows df by column_name_for_timestamp column
        df_container.sort_values(by=column_name_for_timestamp, inplace=True)
        # find the position just before the insertion point for the input timestamp in the sorted array of timestamps
        i = df_container[column_name_for_timestamp].searchsorted(timestamp) - 1
        # if the input timestamp is earlier than the earliest timestamp in the dataframe, skip this container
        if i < 0:
            logging.warning(f'For container {container_id} of plate {plate_id}, '
                            f'the input timestamp {timestamp} is earlier than the earliest timestamp '
                            f'{df_container[column_name_for_timestamp].min()} in the dataframe. '
                            f'This container will be unassigned (skipped).')
            continue
        indices_of_found_conditions.append(df_container.index.values[i])

    return tuple(indices_of_found_conditions)


def find_measurement_files_and_compiile_a_table(target_folder, extension, regexp_pattern,
                                                datetime_string_format='%Y-%m-%d_%H-%M-%S'):
    # make an empty dateframe with columns for filename, plate_id, timestamp
    df_files = pd.DataFrame(columns=['filename', 'plate_id', 'timestamp'], dtype=object)
    for filename in os.listdir(target_folder):
        if not filename.endswith(extension):
            continue
        match = regexp_pattern.match(filename)
        if not match:
            logging.warning(f'Could not parse measurement filename {filename} for timestamp and id')
            continue
        plate_id = int(match.group('plate_id'))
        unix_timestamp = datetime.strptime(match.group('timestamp_string'), datetime_string_format).timestamp()
        logging.info(f'Found measurement file {filename} with plate id {plate_id} and unix timestamp {unix_timestamp}')
        df_files = df_files.append({'filename': filename, 'plate_id': plate_id, 'timestamp': unix_timestamp}, ignore_index=True)
    return df_files


def find_nanodrop_files_and_compiile_a_table(target_folder):
    regexp_pattern = re.compile(
        r'(?P<timestamp_string>\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})(?:.*?)_plate(_?)(?P<plate_id>\d+).csv',
        re.DOTALL)
    return find_measurement_files_and_compiile_a_table(target_folder, extension='.csv', regexp_pattern=regexp_pattern)


def locate_measurement_file_by_condition_timestamp_and_plate_id(df_files, condition_timestamp, plate_id):
    df_files_single_plate = df_files[df_files['plate_id'] == plate_id].copy()
    df_files_single_plate.sort_values(by='timestamp', inplace=True)
    i = df_files_single_plate['timestamp'].searchsorted(condition_timestamp)
    if i > len(df_files_single_plate) - 1:
        # logging.error(f'No measurement file found for plate {plate_id} and condition timestamp {condition_timestamp}, which is later than the latest measurement file.')
        raise ValueError(f'No measurement file found for plate {plate_id} and condition timestamp {condition_timestamp}, which is later than the latest measurement file ({df_files["timestamp"].max()})')
    return df_files_single_plate.index.values[i]


def organize_run_structure_v4_00(experiment_name):
    """
    This function is for nanodrop spectrometer version. It takes the run_info.csv and dilution_info.csv files and
    creates a run_structure.csv file that contains all the information about the run.

     Parameters
     ----------
     experiment_name: str
         The name of the experiment, e.g. `multicomp-reactions/2023-06-26-run01/`
     """
    compatible_versions_of_excel_file = [1.02]
    output_version = '4.00'
    global data_folder

    run_name = experiment_name.split('/')[1]
    run_type = experiment_name.split('/')[0]

    # if there is not 'results' folder in the run folder, create it
    target_folder = data_folder + experiment_name + 'results'
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # open the Excel file with the volumes, use first sheet
    path_to_excel_file = data_folder + experiment_name + f'{run_name}.xlsx'
    assert get_excel_file_version(path_to_excel_file) in compatible_versions_of_excel_file, \
        f"Excel file version is not compatible with this function. " \
        f"Compatible versions are: {compatible_versions_of_excel_file}"
    df_structure = pd.read_excel(path_to_excel_file, sheet_name=0)
    # drop columns that contain 'Unnamed'
    df_structure = df_structure.loc[:, ~df_structure.columns.str.contains('^Unnamed')]
    df_structure['is_outlier'] = 0

    # if there is a file 'outVandC/stock_solutions.xlsx', load it and use to determine concentrations

    # load excel table with stock solution compositions
    df_stock = pd.read_excel(path_to_excel_file, sheet_name='stock_solutions_concentration', skiprows=[0])
    # iterate over rows in df_excel
    for substance in df_stock.columns[1:]:
        df_structure[f'c#{substance}'] = 0

    names_of_stock_solutions_used = [stock_solution_name.replace('vol#','') for stock_solution_name
                                     in df_structure.columns if 'vol#' in stock_solution_name]

    for stock_solution_name in names_of_stock_solutions_used:
        if stock_solution_name not in df_stock['stock_solution'].to_list():
            raise ValueError(f"Stock solution {stock_solution_name} is used in conditions list but is not in the stock solutions table.")
    for index, row in df_structure.iterrows():
        total_volume = 0
        for stock_solution_name in names_of_stock_solutions_used:
            # if this is a column with a stock solution
            if stock_solution_name in df_stock['stock_solution'].to_list():
                concentrations_of_substances_in_this_stock_solution = df_stock.loc[df_stock['stock_solution'] == stock_solution_name]
                pipetted_volume_of_this_stock_solution = row[f'vol#{stock_solution_name}']
                for substance in df_stock.columns[1:]:
                    concentration_of_this_substance = concentrations_of_substances_in_this_stock_solution[substance].iloc[0]
                    df_structure.loc[index, f'c#{substance}'] += pipetted_volume_of_this_stock_solution * concentration_of_this_substance
                total_volume += pipetted_volume_of_this_stock_solution
        for substance in df_stock.columns[1:]:
            df_structure.loc[index, f'c#{substance}'] = df_structure.loc[index, f'c#{substance}'] / total_volume

    # # assign the nanodrop spectra to rows of the df_structure based on plate id and well id
    nanodrop_spectra_folder = data_folder + experiment_name + 'nanodrop_spectra/'
    df_files = find_nanodrop_files_and_compiile_a_table(nanodrop_spectra_folder)
    dilution_related_colnames = [['plate_barcodes_for_dilution', 'timestamp_dilution'],
                                 ['plate_barcodes_for_dilution_2', 'timestamp_dilution_2']]
    for dilution_id in range(len(dilution_related_colnames)):
        this_dilution_plate_colname = dilution_related_colnames[dilution_id][0]
        this_dilution_timestamp_colname = dilution_related_colnames[dilution_id][1]
        suffix_for_dilution_id = f'_{dilution_id+1}' if dilution_id > 0 else ''
        for index, row in df_structure.iterrows():
            # try:
            file_index = locate_measurement_file_by_condition_timestamp_and_plate_id(df_files,
                                                                                     row[this_dilution_timestamp_colname],
                                                                                     row[this_dilution_plate_colname])
            df_structure.loc[index, f'nanodrop_filepath{suffix_for_dilution_id}'] = df_files.loc[file_index, 'filename']
            # except ValueError:
            #     df_structure.loc[index, f'nanodrop_filepath{suffix_for_dilution_id}'] = np.nan

    # nanodrop_spectra_folder = data_folder + experiment_name + 'nanodrop_spectra/'
    # # Use regular expressions to get the plate id and unix timestamp for each file by parsing the filename. For example,
    # # filename '2023-08-23_23-52-13_plate_51.csv' means for 2023-08-23 date, 23-52-13 hour-minute-second, and plate 51
    # ## Old regexp:
    # # regexp_pattern = re.compile(r'(?P<timestamp_string>\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})_plate_(?P<plate_id>\d+).csv')
    #
    # # Changing regexp because Yankai is very creative
    # pattern = re.compile(r'(?P<timestamp_string>\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})(?:.*?)_plate_(?P<plate_id>\d+).csv',
    #                      re.DOTALL)
    # for filename in os.listdir(nanodrop_spectra_folder):
    #     if not filename.endswith('.csv'):
    #         continue
    #     match = pattern.match(filename)
    #     if not match:
    #         logging.warning(f'Could not parse nanodrop spectrum filename {filename} for timestamp and id')
    #         continue
    #     plate_id = int(match.group('plate_id'))
    #     unix_timestamp = datetime.strptime(match.group('timestamp_string'), '%Y-%m-%d_%H-%M-%S').timestamp()
    #     logging.info(
    #         f'Found nanodrop spectrum file {filename} with plate id {plate_id} '
    #         f'and unix timestamp {unix_timestamp}')
    #
    #     # locate the conditions in the df_structure for each file using the locate_condition_by_operation... function
    #     # and write the filepath into the df_structure column 'nanodrop_filepath"
    #     indices = locate_condition_by_operation_datetime_and_plate_id(unix_timestamp, plate_id, df_structure,
    #                                                                   column_name_for_plate_id='plate_barcodes_for_dilution',
    #                                                                   column_name_for_timestamp='timestamp_dilution',
    #                                                                   column_name_for_container_id='container_id')
    #     if len(indices) == 0:
    #         logging.warning(f'Could not find conditions in run structure for nanodrop spectrum file {filename}')
    #         continue
    #     elif len(indices) > 0:
    #         logging.info(f'Found {len(indices)} conditions in run structure for nanodrop spectrum file {filename}')
    #         df_structure.loc[indices, 'nanodrop_filepath'] = experiment_name + 'nanodrop_spectra/' + filename
    #         df_structure.loc[indices, 'timestamp_nanodrop'] = unix_timestamp

    # save the dataframe as csv
    structure_csv_filepath = data_folder + experiment_name + 'results/run_structure.csv'
    with open(structure_csv_filepath, 'w') as f:
        f.write(f'version: {output_version}\n')
    df_structure.to_csv(structure_csv_filepath, index=False, mode='a')
    return df_structure

if __name__ == '__main__':

    # df_files = find_nanodrop_files_and_compiile_a_table(target_folder=data_folder + 'BPRF/2024-01-17-run01/nanodrop_spectra')
    # print(locate_measurement_file_by_condition_timestamp_and_plate_id(df_files, condition_timestamp=1705636817.00000, plate_id=71))


    # target_folder = 'BPRF/2024-03-20-run01/'
    # target_folder = 'C:/Users/ICanSeeYourPixels/Desktop/temp/2024-03-20-run01/'
    # target_folder = 'BPRF/2024-03-06-run02/'
    # target_folder = 'BPRF/2024-03-04-run02/'
    # target_folder = 'BPRF/2024-03-20-run01/'
    # target_folder = 'BPRF/2024-04-17-run01/'
    # target_folder = 'BPRF/2024-04-17-run02/'
    target_folder = 'BPRF/2024-09-30-run02/'
    organize_run_structure_v4_00(experiment_name=target_folder)

    # list_of_runs = tuple(['2024-01-29-run01',
    #                       '2024-01-29-run02',
    #                       '2024-01-30-run01'])
    # for run in list_of_runs:
    #     organize_run_structure_v4_00(experiment_name=f'BPRF/{run}/')

    # list_of_runs = tuple(['2024-02-16-run01',
    #                       '2024-02-17-run01',
    #                       '2024-02-17-run02'])
    # for run in list_of_runs:
    #     organize_run_structure_v4_00(experiment_name=f'BPRF/{run}/')

    # list_of_runs = tuple(['2024-01-16-run01',
    #                       '2024-01-16-run02',
    #                       '2024-01-17-run01'])
    # for run in list_of_runs:
    #     organize_run_structure_v4_00(experiment_name=f'BPRF/{run}/')

    # df_conditions = pd.read_excel('tests/test_organize_run_results/simple-reactions/2099-99-99-run01/2099-99-99-run01.xlsx',
    #                               sheet_name='reactions_with_run_info')
    #
    # timestamp = 1692793066 + 10
    # res = locate_condition_by_operation_datetime_and_plate_id(
    #     timestamp=timestamp,
    #     plate_id=51,
    #     dataframe_with_conditions=df_conditions,
    #     column_name_for_plate_id='plate_barcodes_for_dilution',
    #     column_name_for_timestamp='timestamp_dilution',
    #     column_name_for_container_id='container_id'
    # )
    # print(res)

    # list_of_runs = tuple([
    #                       '2023-11-28-run01',
    #                       '2023-11-29-run01',
    #                       '2023-11-29-run02',
    #                       '2023-12-02-run01',
    #                       '2023-12-04-run01',
    #                       '2023-12-04-run02'])
    # for run_shortname in list_of_runs:
    #     logging.info(f'Organizing run {run_shortname}')
    #     organize_run_structure(f'simple-reactions/{run_shortname}/', version='3.00')

    # list_of_runs = tuple([
    #                       '2023-11-08-run01',
    #                       '2023-11-13-run01',
    #                       '2023-11-14-run01',
    #                       '2023-11-21-run01'
    # ])
    # for run_shortname in list_of_runs:
    #     logging.info(f'Organizing run {run_shortname}')
    #     organize_run_structure(f'BPRF/{run_shortname}/', version='3.00')

    # experiment_name = 'simple-reactions/2023-08-21-run01/'
    # df_structure = organize_run_structure(experiment_name, version='3.00')

    # list_of_runs = tuple([
    #                       '2023-09-06-run01'])
    # for run_shortname in list_of_runs:
    #     organize_run_structure(f'simple-reactions/{run_shortname}/', version='3.00')

    # list_of_runs = tuple([
    #                       '2023-09-14-run01',
    #                       '2023-09-15-run01',
    #                       '2023-09-18-run01',
    #                       '2023-09-19-run01',
    #                       '2023-09-20-run01'
    # ])
    # for run_shortname in list_of_runs:
    #     logging.info(f'Organizing run {run_shortname}')
    #     organize_run_structure(f'simple-reactions/{run_shortname}/', version='3.00')

    # list_of_runs = tuple([
    #                       # '2023-08-21-run01',
    #                       '2023-08-22-run01',
    #                       '2023-08-22-run02',
    #                       '2023-08-28-run01',
    #                       '2023-08-29-run01',
    #                       '2023-08-29-run02'
    # ])
    # for run_shortname in list_of_runs:
    #     organize_run_structure(f'simple-reactions/{run_shortname}/', version='3.00')


    # list_of_runs = tuple([
    #                       '2023-07-05-run01',
    #                       '2023-07-06-run01',
    #                       '2023-07-07-run01',
    #                       '2023-07-10-run01',
    #                       '2023-07-10-run02',
    #                       '2023-07-11-run01',
    #                       '2023-07-11-run02'])
    # for run_shortname in list_of_runs:
    #     organize_run_structure_v2_01(f'simple-reactions/{run_shortname}/')

    # for run_shortname in list_of_runs:
    #     organize_run_structure_v2_00(f'simple-reactions/{run_shortname}/')

    # list_of_runs = tuple([
    #                       '2023-07-05-run01',
    #                       '2023-07-06-run01',
    #                       '2023-07-07-run01',
    #                       '2023-07-10-run01',
    #                       '2023-07-10-run02',
    #                       '2023-07-11-run01',
    #                       '2023-07-11-run02',
    #                       # '2023-07-13-run01' # this run should not be included.
    #                       ])
    # check_run_data_consistency([f'simple-reactions/{run_name}/' for run_name in list_of_runs])

    # ############################# E1
    # list_of_runs = tuple([
    #                       '2023-09-06-run01',
    #                       '2023-09-07-run01'
    # ])
    # for run_shortname in list_of_runs:
    #     organize_run_structure(f'simple-reactions/{run_shortname}/', version='3.00')
