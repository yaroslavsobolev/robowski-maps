import pandas as pd
import numpy as np
import os
import importlib
organize_run_results = importlib.import_module("misc_scripts.organize_run_results")

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

run_name = 'multicomp-reactions/2023-06-19-run01/'

df_data = pd.read_csv(data_folder + run_name + 'results/run_structure.csv')
if 'Unnamed: 0' in df_data.columns:
    df_data.drop('Unnamed: 0', inplace=True, axis=1)
df_results = pd.read_csv(data_folder + run_name + 'results/product_concentration.csv')

assert df_data.shape[0] == df_results.shape[0]

def row_index_by_indices_in_unique_value_lists(df_data, param_values_by_index,
                                         column_names=('ic001', 'am001', 'ald001')):

    target = {column_name: sorted(df_data[column_name].unique())[param_values_by_index[i]]
              for i, column_name in enumerate(column_names)}

    indices = df_data[(df_data[column_names[0]] == target[column_names[0]]) &
                      (df_data[column_names[1]] == target[column_names[1]]) &
                      (df_data[column_names[2]] == target[column_names[2]]) &
                      (df_data[column_names[3]] == target[column_names[3]])
                      ].index

    return indices

def select_rows_by_indices_in_unique_value_lists(df_data, param_values_by_index,
                                         column_names=('ic001', 'am001', 'ald001')):

    return df_data.loc[row_index_by_indices_in_unique_value_lists(df_data=df_data,
                                                                  param_values_by_index=param_values_by_index,
                                                                  column_names=column_names),
                       :]

substances = ['ic001','am001','ald001','ptsa']
product = 'IIO029A'

for substance in substances:
    sorted_unique = sorted(df_results[substance].unique())
    print(f'{len(sorted_unique)} unique concentrations for {substance}')
    print(sorted_unique)

for substance in substances:
    sorted_unique = sorted(df_data[substance].unique())
    print(f'{len(sorted_unique)} unique volumes for {substance}')
    print(sorted_unique)

# df_manual_outliers = pd.DataFrame(dtype=object).reindex_like(df_data)[0:0]

def make_outliers_at_given_locations(locations, output_filename, substances = ('ic001','am001','ald001','ptsa')):

    df_manual_outliers = pd.concat(
            [df_data.loc[row_index_by_indices_in_unique_value_lists(df_results,
                                                                    location,
                                                                    column_names=substances)]
             for location in locations],
        ignore_index=True,
        sort=False)

    # remove known outliers
    df_manual_outliers = df_manual_outliers[df_manual_outliers['is_outlier'] == 0]

    #save manual outliers to csv
    target_folder = data_folder + run_name + 'results/outliers'
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)
    df_manual_outliers.to_csv(data_folder + run_name + f'results/outliers/{output_filename}.csv', index=False)

if __name__ == '__main__':
    # locations = [
    #     [-1, 1, 6, 10],
    #     [-1, 1, 5, 10],
    #     [-1, 0, 4, 10],
    #     [-1, 0, 5, 10],
    #     [-1, 2, 4, 10],
    #     [2, 0, 5, 10],
    #     [-1, 1, 2, 10],
    #     [4, 0, 6, 4],
    #     [4, 0, 4, 4],
    #     [6, 0, 4, 4],
    #     [3, 0, 5, 6],
    #     [6, 0, 4, 7],
    #     [6, 0, 5, 7],
    #     [5, 0, 5, 7],
    #     [-1, 4, 5, 7],
    #     [-1, 1, 5, 8],
    #     [-1, 0, 5, 8],
    #     [3, 0, 2, 8],
    #     [5, 0, 6, 8],
    #     [5, 0, 5, 8],
    #     [5, 0, 4, 8],
    #     [5, 0, 3, 8],
    #     [4, 0, 2, 8],
    #     [-1, 0, 0, 8],
    #     [-2, 0, 2, 9],
    #     [-1, 0, 6, 9],
    #     [-1, 1, 6, 11],
    #     [-1, 1, 1, 11],
    #     [-1, 2, 5, 11],
    #     [-1, 0, 6, 12],
    #     [-1, 1, 5, 12],
    #     [-1, 2, 5, 13],
    #     [5, 4, 0, 14],
    #     [2, 6, 6, 16],
    #     [-1, 4, 5, 16],
    #     [5, 5, 0, 16],
    #     [-1, 5, 6, 17],
    #     [-1, 6, 5, 17],
    #     [-1, 5, 3, 17],
    #     [-1, 5, 6, 18],
    #     [-1, 5, 3, 18],
    #     [4, 6, 0, 18],
    #     [2, 6, 6, 20],
    #     [0, 6, 6, 20],
    #     [2, 6, 6, 20],
    #     [-1, 6, 2, 20]
    # ]
    # make_outliers_at_given_locations(locations=locations, output_filename = 'manual_outliers')

    locations = [
    [-1, 0, 6, 4],
    [-1, 0, 5, 4],
    [-1, 0, 4, 4],
    [ 5, 0, 6, 4],
    [ 5, 0, 5, 4],
    [ 4, 0, 6, 4],
    [-1, 0, 5, 9],
    [-1, 0, 4, 9],
    [-1, 0, 3, 9],
    [-1, 0, 2, 9],
    [ 5, 0, 6, 9],
    [ 5, 0, 5, 9],
    [ 5, 0, 4, 9],
    [ 5, 0, 3, 9],
    [ 4, 0, 6, 9],
    [ 6, 1 ,4, 11],
    [-1, 0, 4, 10],
    [-1, 5, 0, 18],
    ]
    make_outliers_at_given_locations(locations=locations, output_filename='manual_outliers_2')