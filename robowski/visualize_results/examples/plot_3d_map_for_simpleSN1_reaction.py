from robowski.settings import *
import logging
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

import robowski.misc_scripts.organize_run_results as organize_run_results
import robowski.visualize_results.animated_viewer_static as avs
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

# Third batch from december
experiment_name = 'simple-reactions/2023-12-11-run01/'
list_of_runs = tuple([
    '2023-12-11-run01',
    '2023-12-11-run02',
    '2023-12-12-run01',
    '2023-12-12-run02',
    '2023-12-16-run01',
    '2023-12-16-run02'])


column_to_plot = 'yield'

substances = ['c#SN1OH03', 'c#HBr', 'temperature']
substance_titles = ['Alcohol', 'HBr', 'Temperature']
substrates = ['c#SN1OH03', 'c#HBr']

df_results = organize_run_results.join_data_from_runs([f'simple-reactions/{x}/' for x in list_of_runs],
                                 round_on_columns=substances)

for i, row in df_results.iterrows():
    df_results.loc[i, 'c#H2O'] = row['c#HBr'] / 4.5 * 21.88934517
    df_results.loc[i, 'c#acetic_acid'] = row['c#HBr'] / 4.5 * 8.583416744
    df_results.loc[i, 'product_sum'] = row['pc#SN1OH03'] + row['pc#SN1Br03']

df_results.dropna(subset=[column_to_plot], inplace=True)
df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]

# round concentrations "c#SN1OH03" and "c#HBr" to 6 decimal places
df_results['c#SN1OH03'] = df_results['c#SN1OH03'].round(6)
df_results['c#HBr'] = df_results['c#HBr'].round(6)

# iterate over rows and compute column yield2
for i, row in df_results.iterrows():
    limiting_substrate = min(row['pc#SN1OH03'] + row['pc#SN1Br03'], row['c#HBr'] / row['c#SN1OH03'] * (row['pc#SN1OH03'] + row['pc#SN1Br03']))
    if limiting_substrate == 0:
        df_results.loc[i, 'yield2'] = 0
    else:
        df_results.loc[i, 'yield2'] = row['pc#SN1Br03'] / limiting_substrate


for substrate in substrates:
    df_results[substrate] = df_results[substrate].apply(lambda x: x*1000 if x>1e-10 else 0)
xs = df_results[substrates[0]].to_numpy()
ys = df_results[substrates[1]].to_numpy()
zs = df_results['temperature'].to_numpy()

column_to_plot = 'yield'
column_to_plot = 'yield2'

yields = df_results[column_to_plot].to_numpy()
yields[yields > 1] = 1

print(f'Min concentrations of substrates: {[np.min(x) for x in [xs, ys, zs]]}')
print(f'Max concentrations of substrates: {[np.max(x) for x in [xs, ys, zs]]}')
print(f'Yields - min: {min(yields)}, max: {max(yields)}')

avs.plot_3d_dataset_as_cube(xs, ys, zs, yields,
                            substance_titles=('Alcohol,\nmM', 'HBr,\nmM', 'Temperature,\n°C'),
                            colorbar_title=column_to_plot,
                            npoints=50, sparse_npoints=5, rbf_epsilon=0.2,
                            rbf_smooth=0.05,
                            interpolator_choice='rbf',
                            data_for_spheres='raw',#'interpolated',
                            rbf_function='multiquadric',
                            axes_ticks_format='%.0f',
                            axes_font_factor=1.5,
                            contours=[0.1, 0.70, 0.93, 0.98], contour_opacity=0.5) # [0.2, 0.4, 0.55, 0.7, 0.85]

