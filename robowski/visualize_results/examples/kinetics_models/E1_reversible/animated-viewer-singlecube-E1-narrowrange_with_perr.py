import importlib
import logging
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit, brentq
import matplotlib.ticker as mtick
from _decimal import *
from scipy.integrate import solve_ivp
from multiprocessing import Pool
from tqdm import tqdm


def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

if 'ROBOCHEM_DATA_PATH' in os.environ:
    abe = importlib.import_module("visualize_results.examples.kinetics_models.acid_base_equilibrium")

    logging.basicConfig(level=logging.INFO)

    organize_run_results = importlib.import_module("misc_scripts.organize_run_results")
    avs = importlib.import_module("visualize_results.animated_viewer_static")
    data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

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
    df_results.drop(df_results[df_results['is_outlier'] == 1].index, inplace=True)
    df_results.dropna(subset=[column_to_plot], inplace=True)
    df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]
    # df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: x if x > 1e-10 else 0)
    # df_results[column_to_plot] = df_results[column_to_plot].apply(lambda x: x if x <= 100 else 100)

    # round concentrations "c#SN1OH03" and "c#HBr" to 6 decimal places
    df_results['c#E1OH02'] = df_results['c#E1OH02'].round(6)
    df_results['c#HBr'] = df_results['c#HBr'].round(6)

    df_results.to_csv('e1results.csv', index=False)
    njobs = 4
else:
    abe = importlib.import_module("visualize_results.examples.kinetics_models.acid_base_equilibrium")
    df_results = pd.read_csv('e1results.csv')
    njobs = 77

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


def model_of_yield_for_one_condition(index_in_df, pKa_HBr, pKa_COH2plus, K_tilde, k_forward_over_k_backward, k_backward):

    # 10**(-1*(pKa_carbocat - pKa_COH2plus)) * k_forward_over_k_backward = K_tilde
    # therefore, we can find pKa_carbocat
    pKa_carbocat = pKa_COH2plus - np.log10(K_tilde / k_forward_over_k_backward)

    # get row of df_results with index_in_df
    row = df_results.loc[index_in_df]
    # npoints = 1000
    reaction_time_in_hours = 4
    # dt = reaction_time_in_hours / npoints
    # get concentrations of substrates
    c_alcohol_0 = row['c#E1OH02']
    c_hbr = row['c#HBr']
    c_product_0 = 0
    fixed_added_water = 0.231
    moles_of_water_per_moles_of_HBr = 4.7759

    def right_hand_side(t, y):
        c_alcohol, c_product = y
        c_h2o = fixed_added_water + c_hbr * moles_of_water_per_moles_of_HBr + c_product
        # Concentration of C-OH2+ is assumed to be in dynamic equilibrium with the concentration of alcohol
        if c_alcohol < 1e-12:
            c_alcohol = 1e-12

        # OLD VERSION
        # def left_hand_side(x):
        #     return x**2 - K_1*(c_alcohol - x)*(c_hbr - x)
        # # find roots using brentq method
        # c_oh2_plus = brentq(left_hand_side, 0, min(c_alcohol, c_hbr))
        # c_br_minus = c_oh2_plus
        # c_remaining_HBr = c_hbr - c_br_minus

        # NEW VERSION
        # Figuring our the pH and the concentrations of different ions
        substances = ({'pKa': pKa_HBr, 'conc': c_hbr, 'charge': 0}, # HBr + H2O -> H3O+ + Br-
                      {'pKa': pKa_COH2plus, 'conc': c_alcohol, 'charge': 1}, # COH2+ + H2O -> H3O+ + COH
                        {'pKa': pKa_carbocat, 'conc': c_product, 'charge': 1} # Carbocation+ + H2O -> H3O+ + Product
                      )
        Kw_here = Decimal(10) ** Decimal(-14) / (abe.PURE_WATER_MOLARITY) ** 2
        try:
            ph_here = abe.solve_for_zero_charge(water_concentration=c_h2o,
                                                kW=Kw_here, substances=substances)
        except Exception as e:
            print(f'Failed to solve for zero charge at index {index_in_df}.')
            raise e
        subs, remaining_water, oh_minus = abe.concentrations_by_ph(water_concentration=c_h2o, Kw=Kw_here,
                                                               substances=substances,
                                                               ph=ph_here, return_Decimals=False)
        c_oh2_plus = subs[1]
        c_product_after_prot_exchange = subs[2]
        # find the reaction rates at this timepoint
        k_forward = k_forward_over_k_backward * k_backward
        rate_of_forward_reaction = k_forward * c_oh2_plus['conc_prot']
        rate_of_backward_reaction = k_backward * c_product_after_prot_exchange['conc_prot'] * remaining_water
        overall_rate = rate_of_forward_reaction - rate_of_backward_reaction
        return np.array((-1*overall_rate, overall_rate))

    y0 = np.array([c_alcohol_0, c_product_0])

    # use solve_ivp to solve the system of ODEs
    sol = solve_ivp(right_hand_side, (0, reaction_time_in_hours), y0, t_eval=[reaction_time_in_hours],
                    method='RK45', first_step=0.000001, max_step=1, atol=1e-8, rtol=1e-7)
    # sol = solve_ivp(right_hand_side, (0, reaction_time_in_hours), y0, t_eval=[reaction_time_in_hours],
    #                 method='RK45', first_step=0.01, max_step=2, atol=1e-4, rtol=1e-3)
    c_product = sol.y[1][-1]

    reaction_yield = c_product / row['c#E1OH02']
    return reaction_yield


def work_wrapper(packed_args):
    i, pKa_HBr, pKa_COH2plus, K_tilde, k_forward_over_k_backward, k_backward = packed_args
    return model_of_yield_for_one_condition(i, pKa_HBr, pKa_COH2plus, K_tilde, k_forward_over_k_backward, k_backward)


def model_of_yield_for_many_conditions(indices, pKa_HBr, pKa_COH2plus, K_tilde, k_forward_over_k_backward, k_backward):
    print(f'Evaluating yield model for many conditions... {pKa_HBr}, {pKa_COH2plus}, {K_tilde}, {k_forward_over_k_backward}, {k_backward}')
    # # nonparallel version
    # return [model_of_yield_for_one_condition(index_in_df=i, pKa_HBr=pKa_HBr, pKa_COH2plus=pKa_COH2plus,
    #                                          K_tilde=K_tilde, k_forward_over_k_backward=k_forward_over_k_backward,
    #                                          k_backward=k_backward)
    #         for i in indices]

    # parallel version
    list_of_parameter_lists = [(i, pKa_HBr, pKa_COH2plus, K_tilde, k_forward_over_k_backward, k_backward) for i in indices]
    with Pool(njobs) as pool:
        model = list(tqdm(pool.imap(work_wrapper, list_of_parameter_lists), total=len(indices)))
    return model

# predicted_yield = model_of_yield_for_one_condition(index_in_df=0, K_1 = 100, k_forward_over_k_backward=1, k_backward=1)
# print(f'predicted yield {predicted_yield:.2f}')

# filtering
# df_results.drop(df_results[df_results['c#HBr'] < 0.065].index, inplace=True)
# df_results.drop(df_results[df_results['yield'] > 0.96].index, inplace=True)
# df_results.drop(df_results[df_results['is_outlier'] == 1].index, inplace=True)


# # iterate over unique values of "c#SN1OH03" and "temperature" and smooth the data across the 'c#HBr' values
# for alcohol_concentration in df_results['c#SN1OH03'].unique():
#     for temperature in df_results['temperature'].unique():
#         print(f'c#SN1OH03 = {alcohol_concentration}, temperature = {temperature}')
#         indices = (df_results['c#SN1OH03'] == alcohol_concentration) & \
#                                                    (df_results['temperature'] == temperature)
#         df_results.loc[indices, column_to_plot] = \
#             smooth_across_HBr_concentrations(df_results.loc[indices, 'c#HBr'],
#                                              df_results.loc[indices, column_to_plot]
#                                              )

def fit_kinetic_model(indices_here, do_plot=False,
                      p0=tuple([-9.00000001, 0.48564436, 4.05143725, 0.06418353, 7.76078559])):

    def produce_fit(x, y):
        # parameters are pKa_HBr, pKa_COH2plus, K_tilde, k_forward_over_k_backward, k_backward
        lower_bounds = [-10, -9, 1e-9, 1e-9, 1e-9]
        upper_bounds = [2, 6, 1e6, 1e6, 1e9]
        popt, pcov = curve_fit(model_of_yield_for_many_conditions, x, y, p0=p0,
                               max_nfev=1000, bounds=(lower_bounds, upper_bounds), verbose=2,
                               x_scale=[0.1, 0.01, 1, 0.01, 1], diff_step=[0.0001, 0.001, 0.0001, 0.0001, 0.0001],
                               loss='soft_l1', f_scale=0.05, sigma=0.04*y + 0.0001, absolute_sigma=True)
        best_f = lambda x: model_of_yield_for_many_conditions(x, *popt)
        print(f'popt = {popt}')
        perr = np.sqrt(np.diag(pcov))
        # print the rmse of the residuals
        residuals = y - best_f(x)
        rmse = np.sqrt(np.mean(residuals**2))
        print(f'RMSE = {rmse}')
        return best_f, popt, perr

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
    f, keq_fit, keq_err = produce_fit(indices_here, measured_yields)
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

    return keq_fit, keq_err


def plot_kinetic_model(indices_here, popt):
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
    colors_to_plot = df_results.loc[indices_here, 'c#E1OH02'].apply(
        lambda x: colors[np.where(unique_alcohol_concentrations == x)[0][0]])
    xs_to_plot = df_results.loc[indices_here, 'c#HBr']
    measured_yields = df_results.loc[indices_here, 'yield']
    # plt.scatter(xs_to_plot, measured_yields, color='yellow', marker='o')
    f = best_f = lambda x: model_of_yield_for_many_conditions(x, *popt)
    keq_fit = popt
    do_plot = True
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

if __name__ == '__main__':

    do_plot = True
    keq_fits = []
    temperatures = df_results['temperature'].unique()
    temperatures = np.sort(temperatures)
    keq_fit = tuple([-9.000000007982367, 0.7800766169293273, 5.247668539433714, 0.19537521789277223, 2.2963657721096467])
    for temperature in temperatures:
        fig1 = plt.figure(figsize=(4, 3.9), dpi=300)
        print(f'Fitting model at temperature = {temperature}')
        mask = (df_results['temperature'] == temperature)
        indices_where_mask_is_true = df_results[mask].index.to_numpy()

        keq_fit, keq_err = fit_kinetic_model(indices_where_mask_is_true, do_plot=True, p0=keq_fit)
        keq_fits.append((keq_fit, keq_err))

        # plot_kinetic_model(indices_where_mask_is_true, keq_fit)

        if do_plot:
            plt.title(f'Temperature {temperature} °C')
            simpleaxis(plt.gca())
            plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
            plt.ylabel('Yield with respect to alcohol')
            plt.xlabel('Starting concentration of HBr, M')
            plt.tight_layout()
            if 'ROBOCHEM_DATA_PATH' in os.environ:
                plt.gcf().savefig(f'{data_folder}simple-reactions/2023-11-28-run01/results/kinetics/figures/temperature_{temperature}C.png', dpi=300)
                plt.xlim(-0.001, 0.011)
                plt.show()
            else:
                plt.gcf().savefig(
                    f'temperature_{temperature}C.png',
                    dpi=300)
                # close all figs
                plt.close('all')

    print(np.array(keq_fits))

    # save keq_fits as numpy array
    if 'ROBOCHEM_DATA_PATH' in os.environ:
        np.save(f'{data_folder}simple-reactions/2023-09-14-run01/results/kinetics/keq_fits_with_errs.npy', np.array(keq_fits))
    else:
        np.save('keq_fits_with_errs.npy', np.array(keq_fits))
