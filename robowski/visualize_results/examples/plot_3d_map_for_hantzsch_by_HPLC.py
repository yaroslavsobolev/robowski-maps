import importlib
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.integrate import odeint
import itertools
from scipy.optimize import curve_fit

organize_run_results = importlib.import_module("misc-scripts.organize_run_results")
avs = importlib.import_module("visualize-results.animated_viewer_static")


data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

# HPLC_80 deg
file_path = r'D:\Dropbox\robochem\data\BPRF\HPLC_results\Hantzsch_80\average_results_HPLCb.csv'
# # HPLC_26 deg
# file_path = r'C:\Users\UNIST\Dropbox\robochem\data\BPRF\HPLC_results\Hantzsch_26\average_results_HPLC_26.csv'
# Load the CSV file into a DataFrame
resultst_df = pd.read_csv(file_path)
# column_to_plot = 'yield'
column_to_plot = 'S'
# 'HE+THP-19k', 'HE+19k-THP', '19k+THP-HE'
# ['HE_c_HPLC','THP_c_HPLC', '19k_c_HPLC', '19j_c_HPLC', '19p_c_HPLC']
# ['HE_ae_HPLC','THP_ae_HPLC', '19k_ae_HPLC', 'total_ae_HPLC']
# substance_titles = ['19c','19a','19b']
substance_titles = ['','','']
xs = resultst_df['c#methoxybenzaldehyde'].to_numpy()*1000
ys = resultst_df['c#ethyl_acetoacetate'].to_numpy()*1000
zs = resultst_df['c#ammonium_acetate'].to_numpy()*1000

# yields
yields = resultst_df[column_to_plot].to_numpy()
contours = np.linspace(min(yields), max(yields), 15)[1:-1].tolist()
# ['HE_y_HPLC','THP_y_HPLC', '19k_y_HPLC', '19j_y_HPLC', '19p_y_HPLC' ]

print(f'Min concentrations of substrates: {[np.min(x) for x in [xs, ys, zs]]}')
print(f'Max concentrations of substrates: {[np.max(x) for x in [xs, ys, zs]]}')
print(f'Yields - min: {min(yields)}, max: {max(yields)}')

avs.plot_3d_dataset_as_cube(xs, ys, zs, yields,
                                substance_titles=substance_titles,
                                colorbar_title=column_to_plot,
                                npoints=50, sparse_npoints=3, rbf_epsilon=0.01,
                                rbf_smooth=0.01,
                                interpolator_choice='rbf',
                                data_for_spheres='interpolated',
                                contour_opacity=0.8,
                                forced_kmax=None,
                                # axes_ticks_format='',
                                axes_ticks_format='%.0f',
                                axes_font_factor=1.0,
                                # contours=[0.2,0.3,0.4])
                                # contours = [0.1,0.2, 0.3, 0.4, 0.5, 0.6])
                                contours=contours)

# # concentrations
# concentrations = resultst_df[column_to_plot].to_numpy()*1000
# print(f'Min concentrations of substrates: {[np.min(x) for x in [xs, ys, zs]]}')
# print(f'Max concentrations of substrates: {[np.max(x) for x in [xs, ys, zs]]}')
# print(f'Product_c - min: {min(concentrations)}, max: {max(concentrations)}')
#
# avs.plot_3d_dataset_as_cube(xs, ys, zs, concentrations,
#                                 substance_titles=substance_titles,
#                                 colorbar_title=column_to_plot,
#                                 npoints=50, sparse_npoints=3, rbf_epsilon=0.01,
#                                 rbf_smooth=0.01,
#                                 interpolator_choice='rbf',
#                                 data_for_spheres='interpolated',
#                                 contour_opacity=0.3,
#                                 forced_kmax=31,
#                                 axes_ticks_format='',
#                                 # axes_ticks_format='%.0f',
#                                 axes_font_factor=1.3,
#                                 # contours=[0.2,0.3,0.4])
#                                 # contours = [0.2, 0.3, 0.4, 0.5, 0.6])
#                                 contours=[5, 10, 15, 20, 25, 30])
# ae
# ae = resultst_df[column_to_plot].to_numpy()
# print(f'Min concentrations of substrates: {[np.min(x) for x in [xs, ys, zs]]}')
# print(f'Max concentrations of substrates: {[np.max(x) for x in [xs, ys, zs]]}')
# print(f'Product_c - min: {min(ae)}, max: {max(ae)}')
#
# avs.plot_3d_dataset_as_cube(xs, ys, zs, ae,
#                                 substance_titles=substance_titles,
#                                 colorbar_title=column_to_plot,
#                                 npoints=50, sparse_npoints=3, rbf_epsilon=0.01,
#                                 rbf_smooth=0.01,
#                                 interpolator_choice='rbf',
#                                 data_for_spheres='interpolated',
#                                 contour_opacity=0.3,
#                                 forced_kmax=0.06,
#                                 axes_ticks_format='',
#                                 # axes_ticks_format='%.0f',
#                                 axes_font_factor=1.3,
#                                 # contours=[0.2,0.3,0.4])
#                                 # contours = [0.2, 0.3, 0.4, 0.5, 0.6])
#                                 contours=[0.02, 0.04])