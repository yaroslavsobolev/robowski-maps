import pandas as pd
import numpy as np
import os
import importlib
visualize_results = importlib.import_module("visualize-results.visualize_results")

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

run_name = 'multicomp-reactions/2023-06-19-run01/'

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

df_results = visualize_results.join_data_from_runs([f'multicomp-reactions/{run}/' for run in list_of_runs])

print(f"There are {df_results[df_results['is_outlier'] == 1].shape[0]} outliers.")
df_results = df_results[df_results['is_outlier'] == 0]


# df_results = pd.read_csv(data_folder + experiment_name + f'results/interpolated_product_concentration.csv')
# replace negative values in yield column with zeros
df_results['yield'] = df_results['yield'].apply(lambda x: 0 if x < 0 else x)

substances = ['ic001','am001','ald001','ptsa']

padding_rows_count = (df_results[substances] == 0).all(axis=1).sum()
print(f"There are {padding_rows_count} padding rows (with zero concentrations of substrates).")
df_results = df_results[(df_results[substances] != 0).any(axis=1)]


substances = ['ic001','am001','ald001','ptsa']
product = 'IIO029A'

substrate_cs = []
for substance in substances:
    substrate_cs.append(df_results[substance].to_numpy())

xs0, ys0, zs0, cats = substrate_cs

print('Max concentrations of substrates: ')
for x in [xs0, ys0, zs0]:
    print(max(x))

minimal_concentration_of_substrates = np.min(np.array([xs0, ys0, zs0]))

unique_cats = sorted(list(set(list(cats))))
print(f'Unique cats: {unique_cats}')


# make interpolation between two conditions
for substance in substances:
    print(substance)
    print(df_results[substance].unique())

# linear interpolation between two conditions
def linear_interp_between_df_rows(number_of_points_between, df_original, targets_ids,
                                  only_use_these=None):
    interp_conditions = df_original.iloc[:0,:].copy()
    for column in interp_conditions.columns:
        if only_use_these is not None:
            interp_conditions[column] = np.linspace(df_original[column].iloc[targets_ids[0]],
                                                    df_original[column].iloc[targets_ids[1]],
                                                    number_of_points_between)[1:-1][only_use_these]
        else:
            interp_conditions[column] = np.linspace(df_original[column].iloc[targets_ids[0]],
                                                   df_original[column].iloc[targets_ids[1]],
                                                   number_of_points_between)[1:-1]
    interp_conditions = interp_conditions.fillna(0)
    return interp_conditions

# Target concentrations, point 1
new_concentrations = df_results.iloc[:0, :].copy()
new_volumes = df_volumes.iloc[:0, :].copy()


def append_interpolation(targets, number_of_points_between, only_use_these=None):
    global new_concentrations, new_volumes
    targets_ids = [
        df_results[(df_results['ic001'] == target['ic001']) & (df_results['am001'] == target['am001'])
                   & (df_results['ald001'] == target['ald001']) & (df_results['ptsa'] == target['ptsa'])].index[0]
        for target in targets
    ]
    interpolated_points_concentration = linear_interp_between_df_rows(number_of_points_between=number_of_points_between,
                                                                      df_original=df_results.copy(),
                                                                      targets_ids=targets_ids,
                                                                      only_use_these=only_use_these)
    new_concentrations = new_concentrations.append(interpolated_points_concentration, ignore_index=True)
    interpolated_points_volumes = linear_interp_between_df_rows(number_of_points_between=number_of_points_between,
                                                                df_original=df_volumes.copy(),
                                                                targets_ids=targets_ids,
                                                                only_use_these=only_use_these)
    new_volumes = new_volumes.append(interpolated_points_volumes, ignore_index=True)

def append_interpolation_by_indices(param_values_by_index_start, param_values_by_index_stop,
                                    number_of_points_between, cats_range_indices,
                                    column_names=('ic001', 'am001', 'ald001'),
                                    catalyst_name='ptsa', only_use_these=None):
    global new_concentrations, new_volumes
    for cats_id in cats_range_indices:
        targets = []
        for param_value_indices in [param_values_by_index_start, param_values_by_index_stop]:
            target = {column_name: sorted(df_results[column_name].unique())[param_value_indices[i]]
                      for i, column_name in enumerate(column_names)}
            target[catalyst_name] = unique_cats[cats_id]
            targets.append(target)
        append_interpolation(targets=targets, number_of_points_between=number_of_points_between,
                             only_use_these=only_use_these)