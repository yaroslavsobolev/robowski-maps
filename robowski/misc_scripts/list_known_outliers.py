from robowski.settings import *
import os

import pandas as pd

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

def join_structures_from_runs(run_names):
    """
    Loads the `results/run_structure.csv` from multiple runs and joins them into one dataframe.

    Parameters
    ----------
    run_names: list
        The names of the runs whose results will be loaded and joined. These names are the names of the subfolders
        of the data folder.

    Returns
    -------
    df_result: pd.DataFrame
        The dataframe with the joined data. The sequence of the rows is the same as the sequence of the runs.
    """
    df_result = pd.read_csv(data_folder + run_names[0] + f'results/run_structure.csv')
    if 'Unnamed: 0' in df_result.columns:
        df_result.drop('Unnamed: 0', inplace=True, axis=1)
    for run_name in run_names[1:]:
        df_temporary = pd.read_csv(data_folder + run_name + f'results/run_structure.csv')
        if 'Unnamed: 0' in df_temporary.columns:
            df_temporary.drop('Unnamed: 0', inplace=True, axis=1)
        df_result = df_result.append(df_temporary, ignore_index=True)
    return df_result

destination_run = 'multicomp-reactions/2023-06-19-run01/'

list_of_runs = tuple(['2023-06-20-run01',
                    '2023-06-21-run01',
                    '2023-06-21-run02',
                    '2023-06-22-run01',
                    '2023-06-22-run02',
                    '2023-06-22-run03',
                    '2023-06-23-run01',
                    '2023-06-23-run02',
                    '2023-06-26-run01',
                    '2023-06-26-run02',
                    '2023-06-27-run01',
                    '2023-06-27-run02',
                    '2023-06-27-run03',
                    '2023-06-28-run01',
                    '2023-06-28-run02',
                    '2023-06-28-run03'])

df_results = join_structures_from_runs([f'multicomp-reactions/{run}/' for run in list_of_runs])

# remove padding rows
substances = ['ic001','am001','ald001','ptsa']
padding_rows_count = (df_results[substances] == 0).all(axis=1).sum()
print(f"There are {padding_rows_count} padding rows (with zero concentrations of substrates).")
df_results = df_results[(df_results[substances] != 0).any(axis=1)]

# save df_results to /results/run_structure.csv
df_results.to_csv(data_folder + destination_run + 'results/run_structure.csv')

print(f"There are {df_results[df_results['is_outlier'] == 1].shape[0]} outliers.")
df_outliers = df_results[df_results['is_outlier'] == 1]

target_folder = data_folder + destination_run + 'results/outliers'
if not os.path.exists(target_folder):
    os.makedirs(target_folder)
df_outliers.to_csv(data_folder + destination_run + 'results/outliers/known_outliers.csv', index=False)

# columns_for_outV_excel = ['reactions', 'DMF', 'ald001', 'ptsa', 'ptsa_dil_x5', 'am001', 'ic001']
# df_outliers[columns_for_outV_excel].to_excel(data_folder + destination_run + 'results/outliers/outV.xlsx',
#                                              sheetname='reactions', index=False)
