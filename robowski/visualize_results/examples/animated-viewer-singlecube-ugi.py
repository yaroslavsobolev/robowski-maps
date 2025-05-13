from robowski.settings import *

import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

# what_to_plot = 'model'
what_to_plot = 'data'

import robowski.misc_scripts.organize_run_results as organize_run_results
import robowski.visualize_results.animated_viewer_static as avs
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

experiment_name = 'multicomp-reactions/2023-06-19-run01/'
if what_to_plot == 'model':
    df_results = pd.read_csv(data_folder + experiment_name + f'results/interpolated_product_concentration_model.csv')
elif what_to_plot == 'data':
    df_results = pd.read_csv(data_folder + experiment_name + f'results/interpolated_product_concentration.csv')


# raw data
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

df_results = organize_run_results.join_data_from_runs([f'multicomp-reactions/{run}/' for run in list_of_runs])

sorted_unique_values_of_ptsa_column = df_results['ptsa'].unique()
sorted_unique_values_of_ptsa_column.sort()
print(sorted_unique_values_of_ptsa_column)

ptsa_targets = [0.05, 0.12, 0.17, 0.24, 0.298]
ptsa_target = ptsa_targets[1]
# first index where ptsa is greater than
ith_ptsa = next(i for i, x in enumerate(sorted_unique_values_of_ptsa_column) if x > ptsa_target)
# ith_ptsa = len(sorted_unique_values_of_ptsa_column) - 1

print(f'i={ith_ptsa}, PTSA: {sorted_unique_values_of_ptsa_column[ith_ptsa]}')
df_results = df_results[df_results['ptsa'] == sorted_unique_values_of_ptsa_column[ith_ptsa]]
# drop nans in the yield column
df_results.dropna(subset=['yield'], inplace=True)

# substrates = ['c#SN1OH03', 'c#HBr']
substances = ['am001','ic001','ald001']
substance_titles = ['Amine', 'Isocyanide', 'Aldehyde']

column_to_plot = 'yield'
# df_results.dropna(subset=[column_to_plot], inplace=True)
# df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]

for substance in substances:
    df_results[substance] = df_results[substance].round(6).astype(np.float64)


# convert from mol/L to mM
for substrate in substances:
    df_results[substrate] = df_results[substrate].apply(lambda x: x*1000 if x>1e-10 else 0)

xs = df_results[substances[0]].to_numpy()
ys = df_results[substances[1]].to_numpy()
zs = df_results[substances[2]].to_numpy()
yields = df_results[column_to_plot].to_numpy() # + df_results['pc#dm40_10'].to_numpy()

print(f'Min concentrations of substrates: {[np.min(x) for x in [xs, ys, zs]]}')
print(f'Max concentrations of substrates: {[np.max(x) for x in [xs, ys, zs]]}')
print(f'Yields - min: {min(yields)}, max: {max(yields)}')

# MAYAVI glitches when plotting degenerate points, so we add a small amount of noise to the data
xs += np.random.normal(0, 1e-12, len(xs))
ys += np.random.normal(0, 1e-12, len(ys))
zs += np.random.normal(0, 1e-12, len(zs))

# avs.plot_3d_dataset_as_cube(xs, ys, zs, yields,
#                             substance_titles=substance_titles,
#                             colorbar_title=column_to_plot,
#                             npoints=50, sparse_npoints=7, rbf_epsilon=1,
#                             rbf_smooth=0.05,
#                             interpolator_choice='linear',
#                             data_for_spheres='raw',
#                             rbf_function='multiquadric',
#                             axes_ticks_format='%.0f',
#                             axes_font_factor=1.3,
#                             contours=[2])

def reorient_callable(scene):
    # scene.scene.camera.position = [3.027169769132524, 2.9982106643412405, 3.0085160168502623]
    # scene.scene.camera.focal_point = [0.5206201337277889, 0.4916610289365053, 0.5019663814455271]
    # scene.scene.camera.view_angle = 30.0
    # scene.scene.camera.view_up = [0.0, 0.0, 1.0]
    # scene.scene.camera.clipping_range = [2.450147039354766, 6.732730810611246]
    # scene.scene.camera.compute_view_plane_normal()
    #
    # scene.scene.camera.position = [3.5535451925675186, 3.524586087776235, 3.5348914402852567]
    # scene.scene.camera.focal_point = [0.5206201337277889, 0.4916610289365053, 0.5019663814455271]
    # scene.scene.camera.view_angle = 30.0
    # scene.scene.camera.view_up = [0.0, 0.0, 1.0]
    # scene.scene.camera.clipping_range = [3.3527389268273087, 7.658115422514912]
    # scene.scene.camera.compute_view_plane_normal()
    #
    # scene.scene.camera.position = [2.5571363102886355, 4.601804632974648, 3.062283677095859]
    # scene.scene.camera.focal_point = [0.5206201337277889, 0.4916610289365053, 0.5019663814455271]
    # scene.scene.camera.view_angle = 30.0
    # scene.scene.camera.view_up = [-0.21520257220714775, -0.4373014770797476, 0.8731868477360951]
    # scene.scene.camera.clipping_range = [3.42782919845877, 7.5635922413913095]
    # scene.scene.camera.compute_view_plane_normal()

    scene.scene.camera.position = [3.5022834678305252, 6.509322279608751, 4.250526934007179]
    scene.scene.camera.focal_point = [0.5206201337277889, 0.4916610289365053, 0.5019663814455271]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [-0.21520257220714775, -0.4373014770797476, 0.8731868477360951]
    scene.scene.camera.clipping_range = [5.8418018839435515, 10.04878256955108]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()

# figure_filename = frepo_data_path + 'misc_scripts/figures/cubes/Ugi-smoothed_{column_to_plot}_ptsa{ith_ptsa}_ptsaval{ptsa_target:.3f}_{what_to_plot}.png'

figure_filename = frepo_data_path + 'misc_scripts/figures/cubes/Ugi-raw_{column_to_plot}_ptsa{ith_ptsa}_ptsaval{ptsa_target:.3f}_{what_to_plot}.png'

avs.plot_3d_dataset_as_cube(xs, ys, zs, yields,
                            substance_titles=('', '', ''),
                            colorbar_title=column_to_plot,
                            npoints=50, sparse_npoints=5, rbf_epsilon=1,
                            rbf_smooth=0.01,
                            interpolator_choice='rbf',
                            data_for_spheres='raw',
                            rbf_function='multiquadric',
                            axes_ticks_format='%.0f',
                            axes_font_factor=1.3,
                            contours=[0.005, 0.05, 0.01, 0.015], # for HRP01
                            # contours=[0.1, 0.9, 2],
                            contour_opacity=0.3,
                            forced_kmax=0.2,
                            dont_mlabshow=False,
                            reorient_callable=reorient_callable,
                            savefig=figure_filename,
                            transparent=True)

