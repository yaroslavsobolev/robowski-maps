import importlib
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
# import plotly.graph_objects as go
from scipy.interpolate import griddata
from mayavi import mlab
import pandas as pd
from scipy.ndimage import gaussian_filter

# what_to_plot = 'model'
what_to_plot = 'data'
do_smooth = False

organize_run_results = importlib.import_module("misc-scripts.organize_run_results")
avs = importlib.import_module("visualize-results.animated_viewer_static")
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

df_results = pd.read_csv(data_folder + 'BPRF/2024-10-01-run01/results/product_concentrations.csv')

substances = ['c#ethyl_acetoacetate',  'c#methoxybenzaldehyde', 'c#ammonium_acetate']
substance_titles = ['Acetoacetate', 'Methoxy', 'Ammonium acetate']
# substrates = ['c#SN1OH03', 'c#HBr']

if what_to_plot == 'model':
    target_folder = 'visualize-results/examples/kinetics_models/'
    df_results = pd.read_hdf(f'{target_folder}hntz_df_results_model_p4.hdf', key='df')

substrates = ['ethyl_acetoacetate',  'methoxybenzaldehyde', 'ammonium_acetate']

# # round values in substance columns to 6 significant digits
# for substance in substances:
#     df_results[substance] = df_results[substance].round(6)

# column_to_plot = 'yield#bb017'
# df_results.dropna(subset=[column_to_plot], inplace=True)
# df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]

# find the 95% perfentile of rmse, LB_stat_dil_0, LB_stat_dil_1
# percentile_to_target = 0.9
# rmse_95 = df_results['rmse'].quantile(percentile_to_target)
# LB_stat_dil_0_95 = df_results['LB_stat_dil_0'].quantile(percentile_to_target)
# LB_stat_dil_1_95 = df_results['LB_stat_dil_1'].quantile(percentile_to_target)

# use only df_results with values below these percentiles
# df_results = df_results[df_results['rmse'] < rmse_95]
# df_results = df_results[df_results['LB_stat_dil_0'] < LB_stat_dil_0_95]
# df_results = df_results[df_results['LB_stat_dil_1'] < LB_stat_dil_1_95]

# df_results = df_results[df_results['yield#HRP01'] < 0.3]

xs = df_results['c#ethyl_acetoacetate']
ys = df_results['c#methoxybenzaldehyde']
plt.scatter(xs, ys)
plt.show()

unique_ethyl_acetoacetate = df_results['c#ethyl_acetoacetate'].unique()
unique_methoxybenzaldehyde = df_results['c#methoxybenzaldehyde'].unique()

def get_xys_for_one_col(column_to_plot, df_results):
    xs = []
    ys = []
    zs = []
    zserr = []
    for ethyl_acetoacetate in unique_ethyl_acetoacetate:
        for methoxybenzaldehyde in unique_methoxybenzaldehyde:
            df = df_results[(df_results['c#ethyl_acetoacetate'] == ethyl_acetoacetate) &
                            (df_results['c#methoxybenzaldehyde'] == methoxybenzaldehyde) &
                            (df_results['c#ammonium_acetate'] > 0.3)]
            if len(df) > 0:
                xs.append(ethyl_acetoacetate)
                ys.append(methoxybenzaldehyde)
                # average values of column to plot
                zs.append(df[column_to_plot].median())
                zserr.append(df[column_to_plot.replace('yield', 'yielderr')].median())
            else:
                print(f'No data for {ethyl_acetoacetate} and {methoxybenzaldehyde}')


    # # use griddata to interpolate to uniform xs ys grid
    # plotsteps = 18
    # grid_x, grid_y = np.mgrid[min(xs):max(xs):plotsteps*1j, min(ys):max(ys):plotsteps*1j]
    # points = np.array([xs, ys]).T
    # grid_z = griddata(points, zs, (grid_x, grid_y), method='linear')
    #
    # # smooth the grid_z with gaussian kernel 5x5
    # if (what_to_plot == 'data') and do_smooth:
    #     grid_z = gaussian_filter(grid_z, sigma=1)
    #
    # gridzerr = griddata(points, zserr, (grid_x, grid_y), method='linear')
    # xs_for_plot = np.linspace(min(xs), max(xs), plotsteps)
    # ys_for_plot = np.linspace(min(ys), max(ys), plotsteps)
    # return grid_x, grid_y, grid_z, gridzerr

    return np.array(xs), np.array(ys), np.array(zs), np.array(zserr)

##### 3D PLOT
data = []

######## MAIN PRODUCTS
data.append(get_xys_for_one_col(column_to_plot = 'yield#HRP01',
                                df_results = df_results[df_results['yield#HRP01'] <= 1]))

# data.append(get_xys_for_one_col(column_to_plot = 'yield#HRP01',
#                                 df_results = df_results[df_results['yield#HRP01'] <= 1]))


# all the rows where 'yield#bb017' is above 1 replace with 1
# df_results['yield#bb017'] = df_results['yield#bb017'].apply(lambda x: 1 if x>1 else x)
# df_results['yield#bb017'] = df_results['yield#bb017']/1.06
data.append(get_xys_for_one_col(column_to_plot = 'yield#bb017',
                                df_results = df_results[df_results['yield#bb017'] <= 1.55]))

# data.append(get_xys_for_one_col(column_to_plot = 'yield#dm053',
#                                 df_results = df_results))


# surfcolors = [(31/255, 119/255, 180/255), (1, 127/255, 14/255)]
surfcolors = [tuple(np.array((33, 64, 154))/255), tuple(np.array((243, 185, 26))/255)]
yieldmax = 100

# ########## intermediates:
# data.append(get_xys_for_one_col(column_to_plot = 'yield#dm40',
#                                 df_results = df_results[df_results['yield#dm40'] <= 1]))
# data.append(get_xys_for_one_col(column_to_plot = 'yield#bb021',
#                                 df_results = df_results[df_results['yield#bb021'] <= 1]))
# surfcolors = [(44/255, 160/255, 44/255), (214/255, 39/255, 40/255)]
# yieldmax = 15

mlab.figure(size=(1224, 1024), bgcolor=(1, 1, 1), fgcolor=(0.2, 0.2, 0.2))

for dataset_id, data_point in enumerate(data):
    grid_x, grid_y, grid_z, gridzerr = data_point
    print(f'min z: {np.min(grid_z)}, max z: {np.max(grid_z)}')
    grid_x *= 1000
    grid_y *= 1000
    grid_z *= 100
    gridzerr *= 100

    # if (what_to_plot == 'model') and (dataset_id == 1):
    #     # this zero-radius point is just to make the z scale up to 60.3 -- same for model and data
    #     grid_x = np.append(grid_x, 10)
    #     grid_y = np.append(grid_y, 100)
    #     grid_z = np.append(grid_z, 60.3)
    #     gridzerr = np.append(gridzerr, 0.1)


    xs0 = (grid_x - np.min(grid_x)) / (np.max(grid_x) - np.min(grid_x))
    ys0 = (grid_y - np.min(grid_y)) / (np.max(grid_y) - np.min(grid_y))
    # zs0 = (grid_z - np.min(grid_z)) / (np.max(grid_z) - np.min(grid_z))
    zs0 = grid_z/yieldmax
    rel_error = gridzerr/grid_z
    max_xs0 = np.max(xs0)
    max_ys0 = np.max(ys0)
    max_zs0 = np.max(zs0)

    print(f'zs0 range: {np.min(zs0)} - {np.max(zs0)}')

    # This makes array of colors based on error: the low error color is more vivid, the high error color is more white
    rel_error = rel_error.flatten()
    print(f'min rel_error: {np.min(rel_error)}, max rel_error: {np.max(rel_error)}')
    error_based_colors = []
    vivid_color_here = surfcolors[dataset_id]
    max_rel_error = np.percentile(rel_error, 60)
    min_rel_error = np.percentile(rel_error, 30)
    for error_here in rel_error:
        # make a linear transition from white to vivid_color_here based on error_here
        error_here = min(error_here, max_rel_error)
        error_here = max(error_here, min_rel_error)
        error_here = error_here - min_rel_error
        error_here = error_here/max_rel_error
        error_based_colors.append(tuple(
            np.array([1, 1, 1]) * error_here + np.array(vivid_color_here) * (1 - error_here)
        ))

    ## Plotting a surface. Discontinued.
    # use mayavi mlab surf
    # mlab.figure(bgcolor=(1, 1, 1))
    # plot = mlab.surf(xs0, ys0, zs0, color=surfcolors[dataset_id])
    # plot = mlab.surf(xs0, ys0, zs0, color=surfcolors[dataset_id], extent=[np.min(xs0), np.max(xs0),
    #                                                                          np.min(ys0), np.max(ys0),
    #                                                                          0, 1],
    #                  opacity=1)


    plot = mlab.points3d(xs0.flatten(), ys0.flatten(), zs0.flatten(), color=vivid_color_here, scale_factor=0.08)




# set(plot,'FaceColor',[1 0 0],'FaceAlpha',0.5)
# mlab.colorbar()
# mlab.xlabel('c#ethyl_acetoacetate')
# mlab.ylabel('c#methoxybenzaldehyde')
# mlab.show()

# plot.actor.actor.property.ambient = 0
# # for i in range(3):
# #     start = np.array([np.min(xs0), np.min(ys0), np.min(zs0)])
# #     end = np.array([np.min(xs0), np.min(ys0), np.min(zs0)])
# #     end[i] = list([max_xs0, max_ys0, max_zs0])[i]
# #     arr = Arrow_From_A_to_B(start[0], start[1], start[2], end[0], end[1], end[2])
# # arr_temp = Arrow_From_A_to_B(np.max(xs0), np.min(ys0), np.min(zs0),
# #                                   np.max(xs0), np.min(ys0), np.max(zs0))
# # plot = mlab.surf(xs0, ys0, (zs0-np.min(zs0))/(np.max(zs0)-np.min(zs0)), color=surfcolors[dataset_id], opacity=0)
#
# plot = mlab.surf(xs0, ys0, (zs0-np.min(zs0))/(np.max(zs0)-np.min(zs0)), color=surfcolors[dataset_id])
#
sparse_npoints=4
ax1 = mlab.axes(color=(0.5, 0.5, 0.5), nb_labels=sparse_npoints, ranges=[np.min(grid_x), np.max(grid_x),
                                                                         np.min(grid_y), np.max(grid_y),
                                                                         0, np.max(grid_z)])

substance_titles = ['Acetoacetate', 'Methoxy', 'Yield, %']
substance_titles = ['']*3
mlab.xlabel(f'{substance_titles[0]}')
mlab.ylabel(f'{substance_titles[1]}')
mlab.zlabel(f'{substance_titles[2]}')

axes_ticks_format='%.0f'
axes_font_factor=1.3
# substance_titles = ['Acetoacetate', 'Methoxy', 'Ammonium acetate']
# mlab.outline(plot)
# cb = mlab.colorbar(object=plot, title=colorbar_title, orientation='horizontal', nb_labels=5)
# cb.scalar_bar.unconstrained_font_size = True
# cb.label_text_property.font_size = 19

# ax1.axes.font_factor = axes_font_factor
# ax1.axes.label_format = axes_ticks_format
# ax1.axes.corner_offset = 0.05

scene = mlab.get_engine().scenes[0]

# scene.scene.camera.position = [2.5924475454278446, 4.014092072475744, 3.196806857037632]
# scene.scene.camera.focal_point = [0.5, 0.5, 0.5]
# scene.scene.camera.view_angle = 30.0
# scene.scene.camera.view_up = [-0.3157787130206195, -0.4521410923093301, 0.834177581242967]
# scene.scene.camera.clipping_range = [3.1635435457316055, 7.0953401952622475]
# scene.scene.camera.compute_view_plane_normal()

# scene.scene.camera.position = [3.5022834678305252, 6.509322279608751, 4.250526934007179]
# scene.scene.camera.focal_point = [0.5206201337277889, 0.4916610289365053, 0.5019663814455271]
# scene.scene.camera.view_angle = 30.0
# scene.scene.camera.view_up = [-0.21520257220714775, -0.4373014770797476, 0.8731868477360951]
# scene.scene.camera.clipping_range = [5.8418018839435515, 10.04878256955108]
# scene.scene.camera.compute_view_plane_normal()

scene.scene.parallel_projection = True
scene.scene.camera.position = [3.534194063285721, 6.570709849223479, 4.1265976964057804]
scene.scene.camera.focal_point = [0.5525307291829846, 0.5530485985512339, 0.37803714384412934]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.21520257220714775, -0.4373014770797476, 0.8731868477360951]
scene.scene.camera.clipping_range = [5.8418018839435515, 10.04878256955108]
scene.scene.camera.compute_view_plane_normal()
if do_smooth:
    smooth_suffix = '-smooth'
else:
    smooth_suffix = ''
scene.scene.save(f'misc-scripts/figures/hantzsch-points3d-NMR-sparse-{what_to_plot}{smooth_suffix}.png')

# scene.scene.light_manager.lights[1].intensity = 0
# scene.scene.light_manager.lights[2].intensity = 0
# scene.scene.light_manager.lights[3].intensity = 0
# camera_light = scene.scene.light_manager.lights[0]
# camera_light.elevation = 47
# camera_light.azimuth = -5.0
scene.scene.render()

# mlab.savefig('misc-scripts/figures/hansch-mainprod.obj')

mlab.show()


# data = []
# data.append(get_xys_for_one_col(column_to_plot = 'yield#HRP01',
#                                 df_results = df_results[df_results['yield#HRP01'] <= 1]))
# # df_results['yield#bb017'] = df_results['yield#bb017'].apply(lambda x: 1 if x>1 else x)
# data.append(get_xys_for_one_col(column_to_plot = 'yield#bb017',
#                                 df_results = df_results[df_results['yield#bb017'] <= 10.55]))
#
# # data.append(get_xys_for_one_col(column_to_plot = 'yield#HRP02',
# #                                 df_results = df_results[df_results['yield#bb017'] <= 10.55]))
#
# # data.append(get_xys_for_one_col(column_to_plot = 'yield#HRI03',
# #                                 df_results = df_results[df_results['yield#bb017'] <= 10.55]))
# # data.append(get_xys_for_one_col(column_to_plot = 'yield#dm40',
# #                                 df_results = df_results[df_results['yield#dm40'] <= 1]))
# # data.append(get_xys_for_one_col(column_to_plot = 'yield#bb021',
# #                                 df_results = df_results[df_results['yield#bb021'] <= 1]))
#
# # make figure with double x axes
# fig = plt.figure(figsize=(4,2.5), dpi=300)
# ax1 = fig.add_subplot(111)
# ax2 = ax1.twiny()
#
# def tick_function(X):
#     V = 110 - X
#     return ["%.0f" % z for z in V]
#
# labels = ['HRP01', 'HRP04']
# colors = [f'C{i}' for i in range(4)]
# colors[1] = 'gold'
# for id, d in enumerate(data):
#     # iterate over unique x and find z where y=x
#     xs, ys, zs, zserr = d
#     xxs = []
#     zzs = []
#     # for x in np.unique(xs):
#
#     # find indices of ys where xs==ys==x
#     indices = np.where(np.isclose(xs+ys, 110/1000))
#     # find z where xs==ys==x
#     z = zs[indices]
#     x = xs[indices]
#     # ax1.plot(110 - x*1000, z*100, label=labels[id], color=colors[id])
#     # plot with zserr errorbars
#     ax1.errorbar(110 - x*1000, z*100, yerr=zserr[indices]*100, fmt='o-', label=labels[id], color=colors[id],
#                  alpha=0.5, capsize=3, capthick=1, elinewidth=1, markersize=3)
#
# ax1.set_xlim(10, 100)
# ax1Ticks = ax1.get_xticks()
# ax2Ticks = ax1Ticks
# ax2.set_xticks(ax2Ticks)
# ax2.set_xbound(ax1.get_xbound())
# ax2.set_xticklabels(tick_function(ax2Ticks))
#
# ax2.set_xlabel("EAA, mM")
# ax1.set_xlabel('Aldehyde, mM')
#
# ax1.set_ylabel('Yield, %')
# ax1.legend()
#
# plt.show()