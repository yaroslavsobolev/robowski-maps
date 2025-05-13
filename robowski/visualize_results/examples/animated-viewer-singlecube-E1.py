from robowski.settings import *

import os
import numpy as np
import robowski.misc_scripts.organize_run_results as organize_run_results
import robowski.visualize_results.animated_viewer_static as avs
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

experiment_name = 'simple-reactions/2023-09-06-run01/'
list_of_runs = tuple([
    '2023-09-06-run01',
    '2023-09-07-run01'
])
# column_to_plot = 'HBr_relative_change'
column_to_plot = 'conversion'

substances = ['c#E1OH02', 'c#HBr', 'temperature']
substance_titles = ['Alcohol', 'HBr', 'Temperature']
substrates = ['c#E1OH02', 'c#HBr']

df_results = organize_run_results.join_data_from_runs([f'simple-reactions/{x}/' for x in list_of_runs],
                                 round_on_columns=None)

# set yields to nan at the minimum of c#SN1OH01 column
# df_results.loc[df_results[substances[0]].round(4) == df_results[substances[0]].round(4).min(), column_to_plot] = np.nan

# # set yields to nan where the temperature is 6 and yield is above 0.9
# df_results.loc[(df_results['temperature'] == 6) & (df_results[column_to_plot] > 0.9), column_to_plot] = np.nan

# set yields to nan where the HBr is zero
# df_results.loc[df_results['c#HBr'] == 0, column_to_plot] = np.nan

# divide yields by 5
# df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: x/5)

# remove all rows where 'yield' columns is nan
# df_results.dropna(subset=[column_to_plot], inplace=True)
# remove all rows where 'yield' columns is inf
# df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]

# df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: 0.5 if x <= 0.5 else 1)

# df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: x if x > 1e-10 else 0)
# df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: x if x <= 1 else 1)

## drop rows where 'is_outlier' columns is equal to 1
df_results.drop(df_results[df_results['is_outlier'] == 1].index, inplace=True)

# convert from mol/L to mM
for substrate in substrates:
    df_results[substrate] = df_results[substrate].apply(lambda x: x*1000 if x>1e-10 else 0)
    # df_results[substrate] = df_results[substrate].round(4)

xs = df_results[substrates[0]].to_numpy()
ys = df_results[substrates[1]].to_numpy()
zs = df_results['temperature'].to_numpy()
yields = df_results[column_to_plot].to_numpy()

print(f'Min concentrations of substrates: {[np.min(x) for x in [xs, ys, zs]]}')
print(f'Max concentrations of substrates: {[np.max(x) for x in [xs, ys, zs]]}')
print(f'Yields - min: {min(yields)}, max: {max(yields)}')

avs.plot_3d_dataset_as_cube(xs, ys, zs, yields,
                            substance_titles=('Alcohol,\nmM', 'HBr,\nmM', 'Temperature,\n°C'),
                            colorbar_title=column_to_plot,
                            npoints=50, sparse_npoints=7, rbf_epsilon=1,
                            rbf_smooth=0.05,
                            interpolator_choice='rbf',
                            data_for_spheres='raw',
                            rbf_function='multiquadric',
                            axes_ticks_format='%.0f',
                            axes_font_factor=1.3,
                            contours=[0.2, 0.4, 0.7, 0.85],
                            contour_opacity=0.1)

# df_results.loc[df_results['yield#SN1Br01s1'] < 0, 'yield#SN1Br01s1'] = 0
# df_results.loc[df_results['yield#SN1Br01s1'] > 1, 'yield#SN1Br01s1'] = 1
# df_results.loc[df_results['c#SN1OH01'].round(4) == df_results['c#SN1OH01'].round(4).min(), 'yield#SN1Br01s1'] = 0
#
# avs.plot_3d_dataset_as_cube(xs, ys, zs, df_results['yield#SN1Br01s1'].to_numpy(),
#                             substance_titles=('Alcohol,\nmM', 'HBr,\nmM', 'temperature,\n°C'),
#                             colorbar_title='yield of SN1Br01s1',
#                             npoints=50, sparse_npoints=7, rbf_epsilon=1,
#                             rbf_smooth=0.05,
#                             interpolator_choice='rbf',
#                             data_for_spheres='interpolated',
#                             rbf_function='multiquadric',
#                             axes_ticks_format='%.0f',
#                             axes_font_factor=1.3,
#                             contours=[0.15, 0.2])