from robowski.settings import *
from traits.api import HasTraits, Range, Instance, \
        on_trait_change
from traitsui.api import View, Item, Group
from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor, \
                MlabSceneModel
import pandas as pd
from scipy.interpolate import Rbf
import numpy as np
import os

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

# timepoint_id = 1
# experiment_name = 'multicomp-reactions/2023-01-18-run01/'
# df_results = pd.read_csv(data_folder + experiment_name + f'results/timepoint{timepoint_id:03d}-reaction_results.csv')

experiment_name = 'multicomp-reactions/2023-03-20-run01/'
df_results = pd.read_csv(data_folder + experiment_name + f'results/product_concentration.csv')
df_results.drop('Unnamed: 0', inplace=True, axis=1)

df_volumes = pd.read_csv(data_folder + experiment_name + f'outVandC/outVRF038202303201421.csv')
df_volumes.drop('Unnamed: 0', inplace=True, axis=1)

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

# # Sampling for 2023-03-29-run01
# # # Interpolating between two maxima at large catalysts
# append_interpolation_by_indices(param_values_by_index_start=(2, 0, 2),
#                                 param_values_by_index_stop=(2, 2, 0),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=16, stop=30))
#
#
# # Interpolating near the corner (max IC, max ALD), side 1
# append_interpolation_by_indices(param_values_by_index_start=(1, 0, 2),
#                                 param_values_by_index_stop=(2, 0, 2),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=11, stop=18))
#
# # Interpolating near the corner (max IC, max ALD), side 2
# append_interpolation_by_indices(param_values_by_index_start=(2, 0, 1),
#                                 param_values_by_index_stop=(2, 0, 2),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=11, stop=18))
#
#
# # Interpolating near the corner (max IC, max ALD), side 3
# append_interpolation_by_indices(param_values_by_index_start=(2, 0, 2),
#                                 param_values_by_index_stop=(2, 1, 2),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=11, stop=18))
#
# # Interpolating near the corner (max IC, min ALD),span 1
# append_interpolation_by_indices(param_values_by_index_start=(2, 0, 0),
#                                 param_values_by_index_stop=(2, 1, 0),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=11, stop=18))
#
# # Interpolating near the corner (max IC, min ALD),span 1
# append_interpolation_by_indices(param_values_by_index_start=(2, 1, 0),
#                                 param_values_by_index_stop=(2, 2, 0),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=11, stop=18))
#
# # Interpolating between large maxima, diagonal
# append_interpolation_by_indices(param_values_by_index_start=(2, 0, 2),
#                                 param_values_by_index_stop=(2, 1, 0),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=16, stop=30))
#
# output_folder = data_folder + 'multicomp-reactions/2023-03-29-run01/outVandC/'

# Adding points to the maximum IC plane
# append_interpolation_by_indices(param_values_by_index_start=(2, 0, 0),
#                                 param_values_by_index_stop= (2, 1, 1),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=16, stop=30))
# append_interpolation_by_indices(param_values_by_index_start=(2, 2, 2),
#                                 param_values_by_index_stop= (2, 1, 1),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=16, stop=30))
# append_interpolation_by_indices(param_values_by_index_start=(2, 1, 2),
#                                 param_values_by_index_stop= (2, 1, 1),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=16, stop=30))
# append_interpolation_by_indices(param_values_by_index_start=(2, 1, 0),
#                                 param_values_by_index_stop= (2, 1, 1),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=16, stop=30))
# append_interpolation_by_indices(param_values_by_index_start=(2, 2, 1),
#                                 param_values_by_index_stop= (2, 1, 1),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=16, stop=30))
# append_interpolation_by_indices(param_values_by_index_start=(2, 0, 1),
#                                 param_values_by_index_stop= (2, 1, 1),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=16, stop=30))
# append_interpolation_by_indices(param_values_by_index_start=(2, 1, 0),
#                                 param_values_by_index_stop= (1, 1, 0),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=16, stop=30))
# append_interpolation_by_indices(param_values_by_index_start=(2, 1, 1),
#                                 param_values_by_index_stop= (1, 1, 0),
#                                 number_of_points_between=3, cats_range_indices=np.arange(start=20, stop=30))



# # Interpolating between two maxima at lower half of catalyst range, only using two points of diagonal
append_interpolation_by_indices(param_values_by_index_start=(2, 0, 2),
                                param_values_by_index_stop=(2, 2, 0),
                                number_of_points_between=10, cats_range_indices=np.arange(start=0, stop=16),
                                only_use_these=[0, 1])

# One point on this diagonal near center at high concentrations of PTSA
append_interpolation_by_indices(param_values_by_index_start=(2, 0, 2),
                                param_values_by_index_stop=(2, 2, 0),
                                number_of_points_between=10, cats_range_indices=np.arange(start=27, stop=30),
                                only_use_these=[3])

# one point at max-IC plane at low catalysts
append_interpolation_by_indices(param_values_by_index_start=(2, 0, 1),
                                param_values_by_index_stop= (2, 1, 1),
                                number_of_points_between=3, cats_range_indices=np.arange(start=0, stop=16))


output_folder = data_folder + 'multicomp-reactions/2023-04-11-run01/outVandC/'

assert np.isclose(new_volumes.sum(axis=1).to_numpy(), 200).all()
new_volumes['DMF'] = 0
new_volumes['DMF'] = 200 - new_volumes.sum(axis=1)
assert (new_volumes.sum(axis=1) == 200).all()

new_volumes.to_csv(output_folder + 'outV.csv')
new_concentrations.drop(['yield', product], inplace=False, axis=1).to_csv(output_folder + 'outC.csv')