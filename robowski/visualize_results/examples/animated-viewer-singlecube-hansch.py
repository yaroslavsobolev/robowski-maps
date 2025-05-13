from robowski.settings import *

import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

import robowski.misc_scripts.organize_run_results as organize_run_results
import robowski.visualize_results.animated_viewer_static as avs
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

# list_of_runs = tuple([
#     '2023-11-08-run01',
#     '2023-11-13-run01',
#     '2023-11-14-run01',
#     '2023-11-21-run01'])
column_to_plot = 'yield'

# 80 degrees
#
# list_of_runs = tuple(['2024-01-29-run01',
#                       '2024-01-29-run02',
#                       '2024-01-30-run01'
#                       ])

list_of_runs = tuple(['2024-02-16-run01',
                      '2024-02-17-run01',
                      '2024-02-17-run02'])

# Room temp
# list_of_runs = tuple(['2024-01-16-run01',
#                       '2024-01-16-run02',
#                       '2024-01-17-run01'])

substances = ['c#ethyl_acetoacetate',  'c#methoxybenzaldehyde', 'c#ammonium_acetate']
substance_titles = ['Acetoacetate', 'Methoxy', 'Ammonium acetate']
# substrates = ['c#SN1OH03', 'c#HBr']

df_results = organize_run_results.join_data_from_runs([f'BPRF/{x}/' for x in list_of_runs],
                                 round_on_columns=None)
df_results = df_results[df_results['c#ammonium_acetate'] > 0.05]
# df_results = df_results[df_results['c#ammonium_acetate'] > 0.01]

print(len(df_results))
# column_to_plot = 'yield#HRP01'
# df_results.dropna(subset=[column_to_plot], inplace=True)
# df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]

for substance in substances:
    df_results[substance] = df_results[substance].round(6).astype(np.float64)


# find top 2 percentile of the 'rmse' column
rmse_thresh = 1000 #df_results['rmse'].quantile(0.95)
# df_results = df_results[df_results['rmse'] < rmse_thresh]

# target_substance = 'HRP01'
# column_to_plot = 'yield#HRP01'

# df_results = df_results[df_results[column_to_plot] <= 1]

# column_to_plot = 'yield#dm37'
# df_results = df_results[df_results[column_to_plot] <= 0.285]

# target_substance = 'bb017'
# column_to_plot = 'yield#bb017'

# column_to_plot = 'yield#dm70'
# df_results = df_results[df_results[column_to_plot] <= 0.1]

# column_to_plot = 'yield#dm035_8_dm35_9'

# column_to_plot = 'pc#dm40_12'
# df_results = df_results[df_results[column_to_plot] <= 0.01]
# column_to_plot = 'rmse'
# column_to_plot = 'fitted_dilution_factor_2'

# column_to_plot = 'pc#bb021'
# df_results = df_results[df_results[column_to_plot] <= 0.020]


# target_substance = 'dm40'
column_to_plot = 'yield#HRP01'
column_to_plot = 'yield#bb017'

substrates = ['ethyl_acetoacetate',  'methoxybenzaldehyde', 'ammonium_acetate']
# product_name = 'bb021'
# for index, row in df_results.iterrows():
#     product_concentration = df_results.loc[index, f'pc#{product_name}']
#     coefficients_dict = {'methoxybenzaldehyde': 2, 'ethyl_acetoacetate': 1, 'ammonium_acetate': 0}
#     candidate_yields = [
#         product_concentration / (df_results.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name])
#         for substrate_name in substrates if substrate_name != 'ammonium_acetate']
#     df_results.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)

# product_name = 'dm40'
# for index, row in df_results.iterrows():
#     product_concentration = df_results.loc[index, f'pc#dm40_10'] + df_results.loc[index, f'pc#dm40_12']
#     coefficients_dict = {'methoxybenzaldehyde': 1, 'ethyl_acetoacetate': 1, 'ammonium_acetate': 0}
#     candidate_yields = [
#         product_concentration / (df_results.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name])
#         for substrate_name in substrates if substrate_name != 'ammonium_acetate']
#     df_results.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)

# column_to_plot = 'yield#dm40'

# df_results.dropna(subset=[column_to_plot], inplace=True)
# df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]
# limit the df to only where 'c#ammonium_acetate" is greater than 0.01
# df_results['fitted_dilution_factor_2'] = df_results['fitted_dilution_factor_2'].apply(lambda x: x/200)
# if yield is above 1, replace with 1

df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: 1 if x>1 else x)

# df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: 1 if x>1 else x)
# negative yields are omitted
# df_results = df_results[df_results[column_to_plot] >= 0]
# df_results = df_results[df_results[column_to_plot] <= 1.3]
# df_results = df_results[df_results[column_to_plot] <= 1.55]
# df_results = df_results[df_results[column_to_plot] <= 2.15]
# df_results = df_results[df_results[column_to_plot] <= 0.35]
# df_results = df_results[df_results[column_to_plot] <= 0.04]

# Smoothing across ammonium acetate concentrations


def smooth_across_concentrations(x, y, mask, do_plot=True):

    def custom_fit_func(x, a, b, c, d, e):
        # return a * (1 - np.exp(-1*x*b)) + c * x + d * x**2 + e * x**3 # for bb017
        # return a * x**4 + b * x**5 + c * x + d * x ** 2 + e * x ** 3 # for hrp01
        return a * np.exp(-1*x*b) + c * x + d * x**2 + e * x**3 # for bb017

    bounds = ([0, 0, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf])
    def produce_fit(x, y):
        popt, pcov = curve_fit(custom_fit_func, x, y, p0=[0.5, 1, 0, 0, 0], maxfev=100000,
                               loss='soft_l1', f_scale=1, bounds=bounds)
        best_f = lambda x: custom_fit_func(x, *popt)
        return best_f

    # def produce_fit(x, y):
        ### Polynomial fit version
        # z = np.polyfit(x, y, 4)
        # f = np.poly1d(z)
        # return f

    x = np.array(x)
    y = np.array(y)
    if do_plot:
        plt.scatter(x, y, color='grey')
        plt.scatter(x[mask], y[mask], color='C0')
    f = produce_fit(x[mask], y[mask])
    if do_plot:
        plt.plot(np.sort(x), f(np.sort(x)), '--', color='C1')

    return f(x)

def smooth_col(column_to_smooth_over, show_smoothing_plots = False):
    # iterate over unique values of "c#SN1OH03" and "temperature" and smooth the data across the 'c#HBr' values
    for ethyl_acetoacetate_concentration in df_results['c#ethyl_acetoacetate'].unique():
        for methoxy_concentration in df_results['c#methoxybenzaldehyde'].unique():
            print(ethyl_acetoacetate_concentration, methoxy_concentration)
            indices = (df_results['c#methoxybenzaldehyde'] == methoxy_concentration) & \
                        (df_results['c#ethyl_acetoacetate'] == ethyl_acetoacetate_concentration)
            # mask where rows at these indices have rmse values above the threshold
            mask = (df_results.loc[indices, 'rmse'] < rmse_thresh)
            df_results.loc[indices, column_to_smooth_over] = \
                smooth_across_concentrations(df_results.loc[indices, 'c#ammonium_acetate'],
                                                 df_results.loc[indices, column_to_smooth_over],
                                                 mask=mask,
                                                 do_plot=show_smoothing_plots)
            if show_smoothing_plots:
                plt.show()

# smooth_col(f'pc#HRP01', show_smoothing_plots=True)
# smooth_col('pc#dm40_10')
# smooth_col('pc#dm40_12')

substrates = ['ethyl_acetoacetate',  'methoxybenzaldehyde', 'ammonium_acetate']
for index, row in df_results.iterrows():
    # HRP01
    product_name = 'HRP01'
    product_concentration = df_results.loc[index, f'pc#{product_name}']
    coefficients_dict = {'methoxybenzaldehyde': 1, 'ethyl_acetoacetate': 2, 'ammonium_acetate': 1}
    candidate_yields = [product_concentration / (df_results.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name]) for substrate_name in substrates]
    df_results.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)

    # bb017
    product_name = 'bb017'
    product_concentration = df_results.loc[index, f'pc#{product_name}']
    coefficients_dict = {'methoxybenzaldehyde': 2, 'ethyl_acetoacetate': 1, 'ammonium_acetate': 2}
    candidate_yields = [product_concentration / (df_results.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name]) for substrate_name in substrates]
    df_results.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)

    # product_name = 'dm40'
    # # for index, row in df_results.iterrows():
    # product_concentration = df_results.loc[index, f'pc#dm40_10'] + df_results.loc[index, f'pc#dm40_12']
    # coefficients_dict = {'methoxybenzaldehyde': 1, 'ethyl_acetoacetate': 1, 'ammonium_acetate': 0}
    # candidate_yields = [
    #     product_concentration / (df_results.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name])
    #     for substrate_name in substrates if substrate_name != 'ammonium_acetate']
    # df_results.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)

# plotting

# df_results = df_results[df_results[column_to_plot] <= 0.009]

print(len(df_results))

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

figure_filename = repo_data_path + f'misc_scripts/figures/cubes/Hansch-80deg-2024-02-16-run01_{column_to_plot}.png'

avs.plot_3d_dataset_as_cube(xs, ys, zs, yields,
                            substance_titles=('', '', ''),
                            colorbar_title=column_to_plot,
                            npoints=50, sparse_npoints=5, rbf_epsilon=1,
                            rbf_smooth=0.01,
                            interpolator_choice='rbf',
                            data_for_spheres='interpolated',
                            rbf_function='multiquadric',
                            axes_ticks_format='%.0f',
                            axes_font_factor=1.3,
                            contours=[0.1, 0.2, 0.4, 0.6, 2], # for HRP01
                            # contours=[0.1, 0.9, 2],
                            contour_opacity=0.3,
                            forced_kmax=1,
                            dont_mlabshow=False,
                            reorient_callable=reorient_callable,
                            savefig=figure_filename,
                            transparent=True)

