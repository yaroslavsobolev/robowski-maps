from robowski.settings import *
import numpy as np
import pandas as pd
import robowski.visualize_results.animated_viewer_static as avs

# Load the CSV file into a DataFrame
resultst_df = pd.read_csv(data_folder + r'\BPRF\HPLC_results\Hantzsch_80\average_results_HPLCb.csv')

column_to_plot = 'S' # Entropy

substance_titles = ['','','']

# convert to mM
xs = resultst_df['c#methoxybenzaldehyde'].to_numpy()*1000
ys = resultst_df['c#ethyl_acetoacetate'].to_numpy()*1000
zs = resultst_df['c#ammonium_acetate'].to_numpy()*1000

# values_to_plot
values_to_plot = resultst_df[column_to_plot].to_numpy()
contours = np.linspace(min(values_to_plot), max(values_to_plot), 15)[1:-1].tolist()

avs.plot_3d_dataset_as_cube(xs, ys, zs, values_to_plot,
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