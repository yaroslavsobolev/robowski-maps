from robowski.settings import *
import os
import numpy as np
import robowski.misc_scripts.organize_run_results as organize_run_results
import robowski.visualize_results.animated_viewer_static as avs
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, brentq
import matplotlib.ticker as mtick
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
#
# # convert from mol/L to mM
# for substrate in substrates:
#     df_results[substrate] = df_results[substrate].apply(lambda x: x*1000 if x>1e-10 else 0)
#     # df_results[substrate] = df_results[substrate].round(4)
#
# xs = df_results[substrates[0]].to_numpy()
# ys = df_results[substrates[1]].to_numpy()
# zs = df_results['temperature'].to_numpy()
# yields = df_results[column_to_plot].to_numpy()
# #
# # logging.info(f'Min concentrations of substrates: {[np.min(x) for x in [xs, ys, zs]]}')
# #
# # [ic(np.min(x)) for x in [xs, ys, zs]]
# #
# # print(f'Max concentrations of substrates: {[np.max(x) for x in [xs, ys, zs]]}')
# # print(f'Yields - min: {min(yields)}, max: {max(yields)}')
# #
# avs.plot_3d_dataset_as_cube(xs, ys, zs, yields,
#                             substance_titles=('Alcohol,\nmM', 'HBr,\nmM', 'Temperature,\n°C'),
#                             colorbar_title=column_to_plot,
#                             npoints=50, sparse_npoints=5, rbf_epsilon=1,
#                             rbf_smooth=0.05,
#                             interpolator_choice='rbf',
#                             data_for_spheres='raw',
#                             rbf_function='multiquadric',
#                             axes_ticks_format='%.0f',
#                             axes_font_factor=1.3,
#                             contours=[0.2, 0.4, 0.7, 0.85],
#                             contour_opacity=0.4)


def model_of_yield_for_one_condition(index_in_df, K_1, k_forward, k_backward):
    # get row of df_results with index_in_df
    row = df_results.loc[index_in_df]
    npoints = 200
    reaction_time_in_hours = 4
    dt = reaction_time_in_hours / npoints
    # get concentrations of substrates
    c_alcohol = row['c#E1OH02']
    c_hbr = row['c#HBr']
    c_product = 0
    for time in np.linspace(0, reaction_time_in_hours, npoints):
        # Concentration of C-OH2+ is assumed to be in dynamic equilibrium with the concentration of alcohol
        def left_hand_side(x):
            return x**2 - K_1*(c_alcohol - x)*(c_hbr - x)
        # find roots using brentq method
        c_oh2_plus = brentq(left_hand_side, 0, min(c_alcohol, c_hbr))

        # find the reaction rates at this timepoint
        rate_of_forward_reaction = k_forward * c_oh2_plus
        rate_of_backward_reaction = k_backward * c_product
        overall_rate = rate_of_forward_reaction - rate_of_backward_reaction

        # update concentrations of alcohol and the product
        c_alcohol -= overall_rate * dt
        c_product += overall_rate * dt

    reaction_yield = c_product / row['c#E1OH02']
    return reaction_yield

def model_of_yield_for_many_conditions(indices, K_1, k_forward, k_backward):
    print('Evaluating yield model for many conditions...')
    return [model_of_yield_for_one_condition(index_in_df=i,
                                             K_1=K_1,
                                             k_forward=k_forward,
                                             k_backward=k_backward)
            for i in indices]

predicted_yield = model_of_yield_for_one_condition(index_in_df=0, K_1 = 100, k_forward=1, k_backward=1)
print(f'predicted yield {predicted_yield:.2f}')

# filtering
# df_results.drop(df_results[df_results['c#HBr'] < 0.065].index, inplace=True)
# df_results.drop(df_results[df_results['yield'] > 0.96].index, inplace=True)
# df_results.drop(df_results[df_results['is_outlier'] == 1].index, inplace=True)
df_results.dropna(subset=[column_to_plot], inplace=True)
df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]
# df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: x if x > 1e-10 else 0)
# df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: x if x <= 100 else 100)

# round concentrations "c#SN1OH03" and "c#HBr" to 6 decimal places
df_results['c#E1OH02'] = df_results['c#E1OH02'].round(6)
df_results['c#HBr'] = df_results['c#HBr'].round(6)


def fit_kinetic_model(indices_here, do_plot=False):

    def produce_fit(x, y):
        popt, pcov = curve_fit(model_of_yield_for_many_conditions, x, y, p0=[16, 0.7, 0.3],
                               max_nfev=1000, bounds=([0]*3, [np.inf]*3),
                               loss='soft_l1', f_scale=0.05)
        best_f = lambda x: model_of_yield_for_many_conditions(x, *popt)
        print(f'popt = {popt}')
        return best_f, popt

    # def produce_fit(x, y):
        ### Polynomial fit version
        # z = np.polyfit(x, y, 4)
        # f = np.poly1d(z)
        # return f

    # make unique values of the alcohol concentrations for indices_here
    unique_alcohol_concentrations = df_results.loc[indices_here, 'c#E1OH02'].unique()
    # sort it
    unique_alcohol_concentrations = np.sort(unique_alcohol_concentrations)
    # make a list of colors based on the id of the alcohol concentration
    colors = [f'C{i}' for i in range(len(unique_alcohol_concentrations))]
    colors_to_plot = df_results.loc[indices_here, 'c#E1OH02'].apply(lambda x: colors[np.where(unique_alcohol_concentrations == x)[0][0]])
    xs_to_plot = df_results.loc[indices_here, 'c#HBr']
    measured_yields = df_results.loc[indices_here, 'yield']
    # plt.scatter(xs_to_plot, measured_yields, color='yellow', marker='o')
    f, keq_fit = produce_fit(indices_here, measured_yields)
    if do_plot:
        plt.scatter(xs_to_plot, measured_yields, color=colors_to_plot, alpha=0.5, s=10)
        for c_alc in unique_alcohol_concentrations:
            color_here = colors[np.where(unique_alcohol_concentrations == c_alc)[0][0]]
            # find df_indices among indices_here where alcolhol concentration is c_alc
            indices_where_mask_is_true = df_results.loc[indices_here, 'c#E1OH02'] == c_alc
            # sort indices by HBr concentration
            xs_here = df_results.loc[indices_here[indices_where_mask_is_true], 'c#HBr']
            ys_here = f(indices_here[indices_where_mask_is_true])
            # sort xs and ys by increasing xs
            xs_here, ys_here = zip(*sorted(zip(xs_here, ys_here)))
            plt.plot(xs_here, ys_here, color=color_here, label=f'{c_alc:.3f} M')
        # plt.scatter(xs_to_plot, ys_to_plot, color=colors_to_plot, marker='x')
        simpleaxis(plt.gca())
        plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
        plt.ylabel('Yield with respect to alcohol')
        plt.xlabel('Starting concentration of HBr, M')
        plt.tight_layout()
        plt.legend(title="Starting alcohol\nconcentration")

    return keq_fit

do_plot = True
keq_fits = []
temperatures = df_results['temperature'].unique()
temperatures = np.sort(temperatures)
for temperature in temperatures:
    fig1 = plt.figure(figsize=(4, 3.9), dpi=300)
    print(f'Fitting model at temperature = {temperature}')
    mask = (df_results['temperature'] == temperature)
    indices_where_mask_is_true = df_results[mask].index.to_numpy()
    keq_fit = fit_kinetic_model(indices_where_mask_is_true, do_plot=True)
    keq_fits.append(keq_fit)
    if do_plot:
        plt.title(f'Temperature {temperature} °C')
        simpleaxis(plt.gca())
        plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
        plt.ylabel('Yield with respect to alcohol')
        plt.xlabel('Starting concentration of HBr, M')
        plt.tight_layout()
        plt.gcf().savefig(f'{data_folder}simple-reactions/2023-11-28-run01/results/kinetics/figures/temperature_{temperature}C.png', dpi=300)
        plt.xlim(-0.001, 0.011)
        plt.show()

# save keq_fits as numpy array
# np.save(f'{data_folder}simple-reactions/2023-09-14-run01/results/kinetics/keq_fits.npy', np.array(keq_fits))
