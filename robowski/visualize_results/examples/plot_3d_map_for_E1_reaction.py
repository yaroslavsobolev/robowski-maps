from robowski.settings import *
import robowski.misc_scripts.organize_run_results as organize_run_results
import robowski.visualize_results.animated_viewer_static as avs
import os
import logging

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

logging.basicConfig(level=logging.INFO)


experiment_name = 'simple-reactions/2023-09-14-run01/'
list_of_runs = tuple([
    '2023-09-14-run01',
    '2023-09-15-run01',
    '2023-09-18-run01',
    '2023-09-19-run01',
    '2023-09-20-run01'
])

column_to_plot = 'conversion'

substances = ['c#E1OH02', 'c#HBr', 'temperature']
substance_titles = ['Alcohol', 'HBr', 'Temperature']
substrates = ['c#E1OH02', 'c#HBr']

df_results = organize_run_results.join_data_from_runs([f'simple-reactions/{x}/' for x in list_of_runs],
                                 round_on_columns=None)

## drop rows where 'is_outlier' columns is equal to 1
df_results.drop(df_results[df_results['is_outlier'] == 1].index, inplace=True)

# convert from mol/L to mM
for substrate in substrates:
    df_results[substrate] = df_results[substrate].apply(lambda x: x*1000 if x>1e-10 else 0)

xs = df_results[substrates[0]].to_numpy()
ys = df_results[substrates[1]].to_numpy()
zs = df_results['temperature'].to_numpy()
yields = df_results[column_to_plot].to_numpy()

avs.plot_3d_dataset_as_cube(xs, ys, zs, yields,
                            substance_titles=('Alcohol,\nmM', 'HBr,\nmM', 'Temperature,\n°C'),
                            colorbar_title=column_to_plot,
                            npoints=50, sparse_npoints=5, rbf_epsilon=1,
                            rbf_smooth=0.05,
                            interpolator_choice='rbf',
                            data_for_spheres='raw',
                            rbf_function='multiquadric',
                            axes_ticks_format='%.0f',
                            axes_font_factor=1.3,
                            contours=[0.2, 0.4, 0.7, 0.85],
                            contour_opacity=0.4)