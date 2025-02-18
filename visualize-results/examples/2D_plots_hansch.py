import importlib
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import plotly.graph_objects as go
from scipy.interpolate import griddata

organize_run_results = importlib.import_module("misc-scripts.organize_run_results")
avs = importlib.import_module("visualize-results.animated_viewer_static")
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

# list_of_runs = tuple(['2024-01-29-run01',
#                       '2024-01-29-run02',
#                       '2024-01-30-run01'
#                       ])

list_of_runs = tuple(['2024-02-16-run01',
                      '2024-02-17-run01',
                      '2024-02-17-run02'])


substances = ['c#ethyl_acetoacetate',  'c#methoxybenzaldehyde', 'c#ammonium_acetate']
substance_titles = ['Acetoacetate', 'Methoxy', 'Ammonium acetate']
# substrates = ['c#SN1OH03', 'c#HBr']

df_results = organize_run_results.join_data_from_runs([f'BPRF/{x}/' for x in list_of_runs],
                                 round_on_columns=None)

# round values in substance columns to 6 significant digits
for substance in substances:
    df_results[substance] = df_results[substance].round(6)

# column_to_plot = 'yield#HRP01'

# column_to_plot = 'yield#dm37'
# df_results = df_results[df_results[column_to_plot] <= 0.285]

column_to_plot = 'yield#bb017'

# column_to_plot = 'yield#dm70'
# column_to_plot = 'yield#dm035_8_dm35_9'
# column_to_plot = 'rmse'
# column_to_plot = 'fitted_dilution_factor_2'

df_results.dropna(subset=[column_to_plot], inplace=True)
df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]
# limit the df to only where 'c#ammonium_acetate" is greater than 0.01
df_results = df_results[df_results['c#ammonium_acetate'] > 0.01]
# df_results['fitted_dilution_factor_2'] = df_results['fitted_dilution_factor_2'].apply(lambda x: x/200)
# if yield is above 1, replace with 1
# df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: 1 if x>1 else x)
# negative yields are omitted
# df_results = df_results[df_results[column_to_plot] >= 0]
# df_results = df_results[df_results[column_to_plot] <= 1.3]

# df_results = df_results[df_results[column_to_plot] <= 2.15]
# df_results = df_results[df_results[column_to_plot] <= 0.35]
# df_results = df_results[df_results[column_to_plot] <= 0.04]

unique_ethyl_acetoacetate = df_results['c#ethyl_acetoacetate'].unique()
unique_methoxybenzaldehyde = df_results['c#methoxybenzaldehyde'].unique()

def get_xys_for_one_col(column_to_plot, df_results):
    xs = []
    ys = []
    zs = []
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
            else:
                print(f'No data for {ethyl_acetoacetate} and {methoxybenzaldehyde}')

    # use griddata to interpolate to uniform xs ys grid
    plotsteps = 30
    grid_x, grid_y = np.mgrid[min(xs):max(xs):plotsteps*1j, min(ys):max(ys):plotsteps*1j]
    points = np.array([xs, ys]).T
    grid_z = griddata(points, zs, (grid_x, grid_y), method='linear')
    xs_for_plot = np.linspace(min(xs), max(xs), plotsteps)
    ys_for_plot = np.linspace(min(ys), max(ys), plotsteps)
    return xs_for_plot, ys_for_plot, grid_z


data = []
data.append(get_xys_for_one_col(column_to_plot = 'yield#HRP01',
                                df_results = df_results[df_results['yield#HRP01'] <= 1]))
data.append(get_xys_for_one_col(column_to_plot = 'yield#dm37',
                                df_results = df_results[df_results['yield#dm37'] <= 0.285]))
data.append(get_xys_for_one_col(column_to_plot = 'yield#bb017',
                                df_results = df_results[df_results['yield#bb017'] <= 1.55]))

# # use imshow, and plot points as well
# fig, ax = plt.subplots()
# ax.imshow(grid_z, extent=(min(xs), max(xs), min(ys), max(ys)), origin='lower', aspect='auto')
# ax.scatter(xs, ys, c=zs, s=100, edgecolors='k', cmap='viridis', alpha=0.7)
# ax.set_xlabel('c#ethyl_acetoacetate')
# ax.set_ylabel('c#methoxybenzaldehyde')
# plt.colorbar(ax.scatter(xs, ys, c=zs, s=100, edgecolors='k', cmap='viridis', alpha=0.7))
# plt.show()

# make 2d surface with plotly
grid_z = data[0][2]
cmap = plt.get_cmap("viridis")
colorscale = [[0, 'rgb' + str(cmap(0)[0:3])],
              [0.5, 'rgb' + str(cmap(0.5)[0:3])],
              [1, 'rgb' + str(cmap(0.9)[0:3])]]
colors1 = np.zeros(shape=grid_z.shape)
colors2 = 0.5*np.ones(shape=grid_z.shape)
colors3 = np.ones(shape=grid_z.shape)

opacity = 0.3
fig = go.Figure(data=[
    go.Surface(z=data[0][2], x=data[0][0]*1000, y=data[0][1]*1000, surfacecolor=colors1, colorscale=colorscale, cmin=0, cmax=1, opacity=opacity),
    go.Surface(z=data[1][2], x=data[1][0]*1000, y=data[1][1]*1000, surfacecolor=colors2, colorscale=colorscale, cmin=0, cmax=1, opacity=opacity),
    go.Surface(z=data[2][2], x=data[2][0]*1000, y=data[2][1]*1000, surfacecolor=colors3, colorscale=colorscale, cmin=0, cmax=1, opacity=opacity),
])
fig.update_layout(scene = dict(
                    xaxis_title='methoxybenzaldehyde',
                    yaxis_title='ethyl_acetoacetate',
                    zaxis_title='Yield'),
                    margin=dict(l=0, r=0, b=0, t=0))
# save to html
fig.write_html("misc-scripts/figures/2d_plot_surface-hansch-smooth.html")
fig.show()