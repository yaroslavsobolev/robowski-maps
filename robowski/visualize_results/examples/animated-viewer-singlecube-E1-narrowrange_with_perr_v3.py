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

R_gas = 8.31446261815324  # J/(K*mol)

if 'ROBOCHEM_DATA_PATH' in os.environ:
    activity_folder = 'misc_scripts/activity_data/e1_hanna/'
else:
    activity_folder = 'e1_hanna/'
df_water_activities = pd.read_csv(f'{activity_folder}activity_coeffs_vs_T_for_water.csv')
df_substrate_activity = pd.read_csv(f'{activity_folder}activity_coeffs_vs_T_for_alcohol_neutral.csv')
df_product_activity = pd.read_csv(f'{activity_folder}activity_coeffs_vs_T_for_product_neutral.csv')
acetonitrile_molarity = 19.145982 # mol / L

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

    # df_results.to_csv('e1results.csv', index=False)
    njobs = 4
else:
    abe = importlib.import_module("visualize_results.examples.kinetics_models.acid_base_equilibrium")
    df_results = pd.read_csv('e1results.csv')
    njobs = 77

# PURE_WATER_MOLARITY is actually replaced here with pure acetone molarity
abe.PURE_WATER_MOLARITY = Decimal(19.1473812) # mol/L
abe.LOG10_PURE_WATER_MOLARITY = abe.PURE_WATER_MOLARITY.log10()
abe.LOG10_PURE_WATER_MOLARITY_FLOAT = np.log10(19.1473812)

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


def model_of_yield_for_one_condition(index_in_df, pKa_COH2plus_deltaS, pKa_COH2plus_deltaH, K_tilde_deltaS, K_tilde_deltaH, k_forward_over_k_backward_deltaS, k_forward_over_k_backward_deltaH, k_backward_kappa, k_backward_deltaS, k_backward_deltaH, HBr_B):
    # NOTE: in the research article, the signs in the definitions of enthalpies except K_tilde_deltaH are reversed
    # in comparison to this code.

    # get row of df_results with index_in_df
    row = df_results.loc[index_in_df]
    # pKa_HBr = 6.6


    # npoints = 1000
    reaction_time_in_hours = 4
    # dt = reaction_time_in_hours / npoints
    # get concentrations of substrates
    c_alcohol_0 = row['c#E1OH02']
    c_hbr = row['c#HBr']
    c_product_0 = 0
    temperature_here = row['temperature']
    fixed_added_water = 0.231
    moles_of_water_per_moles_of_HBr = 4.7759

    # temperature dependence for equilibrium constant of the entire reaction
    x = 1000/(273.15 + temperature_here)
    a = K_tilde_deltaS / R_gas
    b = -1 * K_tilde_deltaH / R_gas
    # a = 28.39501394273001, b = -7.8753586965650655
    corrected_K_tilde = np.exp(a + b * x)
    # but, corrected_keq = K_tilde / acetonitrile_molarity * water_activity_here * product_activity_here / substrate_activity_here
    # therefore, we can find K_tilde
    water_activity_here = df_water_activities.loc[df_water_activities['T_C'] == temperature_here, 'gamma'].values[0]
    substrate_activity_here = df_substrate_activity.loc[df_substrate_activity['T_C'] == temperature_here, 'gamma'].values[0]
    product_activity_here = df_product_activity.loc[df_product_activity['T_C'] == temperature_here, 'gamma'].values[0]
    K_tilde = corrected_K_tilde * acetonitrile_molarity / water_activity_here / product_activity_here * substrate_activity_here

    # temperature dependence for hydration equilibrium
    a = k_forward_over_k_backward_deltaS
    b = k_forward_over_k_backward_deltaH
    # a = 22.97029428800321, b = -6.769846371083024
    k_forward_over_k_backward = np.exp(a + b * x)

    # pKa of COH2+ is temperature dependent
    a = pKa_COH2plus_deltaS
    b = pKa_COH2plus_deltaH
    # a = 44.47454362481631, b = -10.490899871702174
    pKa_COH2plus = a + b * x

    # 10**(-1*(pKa_carbocat - pKa_COH2plus)) * k_forward_over_k_backward = K_tilde
    # therefore, we can find pKa_carbocat
    pKa_carbocat = pKa_COH2plus - np.log10(K_tilde / k_forward_over_k_backward)

    # pKa of HBr.
    # It is known that pKa_HBr is 6.6 at 25C. So if pKa_HBr = a + b * 1000/T, then we can find 'a' for a given 'b'
    b = HBr_B
    a = 6.6 - b * 1000/(273.15 + 25)
    pKa_HBr = a + b * x

    # Backward equilibrium constant
    # kappa_eff = 5.180674455597947e-06, a = -4.617352877847709, b = 3.084571995450058
    kappa_eff = k_backward_kappa
    a = k_backward_deltaS
    b = k_backward_deltaH
    k_backward = kappa_eff * (1000/x) * np.exp(a + b * x)

    def right_hand_side(t, y):
        c_alcohol, c_product = y
        c_h2o = fixed_added_water + c_hbr * moles_of_water_per_moles_of_HBr + c_product
        # Concentration of C-OH2+ is assumed to be in dynamic equilibrium with the concentration of alcohol
        if c_alcohol < 1e-12:
            c_alcohol = 1e-12
        if c_product < 0:
            c_product = 0

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

        # Autoprotolysis constant taken from the following article:
        # Barbosa, J. & Sanz-Nebot, V. Autoprotolysis constants and standardization of the glass electrode in
        # acetonitrile-water mixtures. Effect of solvent composition. Anal. Chim. Acta 244, 183–191 (1991).

        # PURE_WATER_MOLARITY is actually replaced here with pure acetone molarity
        Kw_here = Decimal(10) ** Decimal(-33.58) / (abe.PURE_WATER_MOLARITY) ** 2
        try:
            ph_here = abe.solve_for_zero_charge(water_concentration=abe.PURE_WATER_MOLARITY,
                                                kW=Kw_here, substances=substances)
        except Exception as e:
            print(f'Failed to solve for zero charge at index {index_in_df}.')
            raise e
        subs, remaining_water, oh_minus = abe.concentrations_by_ph(water_concentration=abe.PURE_WATER_MOLARITY, Kw=Kw_here,
                                                               substances=substances,
                                                               ph=ph_here, return_Decimals=False)
        c_oh2_plus = subs[1]
        c_product_after_prot_exchange = subs[2]
        # find the reaction rates at this timepoint
        k_forward = k_forward_over_k_backward * k_backward
        rate_of_forward_reaction = k_forward * c_oh2_plus['conc_prot']
        rate_of_backward_reaction = k_backward * c_product_after_prot_exchange['conc_prot'] * c_h2o
        overall_rate = rate_of_forward_reaction - rate_of_backward_reaction
        return np.array((-1*overall_rate, overall_rate))

    y0 = np.array([c_alcohol_0, c_product_0])

    # use solve_ivp to solve the system of ODEs
    sol = solve_ivp(right_hand_side, (0, reaction_time_in_hours), y0, t_eval=[reaction_time_in_hours],
                    method='RK45', first_step=0.00000001, max_step=1, atol=1e-8, rtol=1e-7)
    # sol = solve_ivp(right_hand_side, (0, reaction_time_in_hours), y0, t_eval=[reaction_time_in_hours],
    #                 method='RK45', first_step=0.01, max_step=2, atol=1e-4, rtol=1e-3)
    c_product = sol.y[1][-1]

    reaction_yield = c_product / row['c#E1OH02']
    return reaction_yield


def work_wrapper(packed_args):
    i, pKa_COH2plus_deltaS, pKa_COH2plus_deltaH, K_tilde_deltaS, K_tilde_deltaH, k_forward_over_k_backward_deltaS, k_forward_over_k_backward_deltaH, k_backward_kappa, k_backward_deltaS, k_backward_deltaH, HBr_B = packed_args
    return model_of_yield_for_one_condition(i, pKa_COH2plus_deltaS, pKa_COH2plus_deltaH, K_tilde_deltaS, K_tilde_deltaH, k_forward_over_k_backward_deltaS, k_forward_over_k_backward_deltaH, k_backward_kappa, k_backward_deltaS, k_backward_deltaH, HBr_B)


def model_of_yield_for_many_conditions(indices, pKa_COH2plus_deltaS, pKa_COH2plus_deltaH, K_tilde_deltaS, K_tilde_deltaH, k_forward_over_k_backward_deltaS, k_forward_over_k_backward_deltaH, k_backward_kappa, k_backward_deltaS, k_backward_deltaH, HBr_B):
    print(f'Evaluating yield model for many conditions, params are {pKa_COH2plus_deltaS, pKa_COH2plus_deltaH, K_tilde_deltaS, K_tilde_deltaH, k_forward_over_k_backward_deltaS, k_forward_over_k_backward_deltaH, k_backward_kappa, k_backward_deltaS, k_backward_deltaH, HBr_B}')
    # # nonparallel version
    # return [model_of_yield_for_one_condition(index_in_df=i, pKa_HBr=pKa_HBr, pKa_COH2plus=pKa_COH2plus,
    #                                          K_tilde=K_tilde, k_forward_over_k_backward=k_forward_over_k_backward,
    #                                          k_backward=k_backward)
    #         for i in indices]

    # parallel version
    list_of_parameter_lists = [(i, pKa_COH2plus_deltaS, pKa_COH2plus_deltaH, K_tilde_deltaS, K_tilde_deltaH, k_forward_over_k_backward_deltaS, k_forward_over_k_backward_deltaH, k_backward_kappa, k_backward_deltaS, k_backward_deltaH, HBr_B) for i in indices]
    with Pool(njobs) as pool:
        model = list(tqdm(pool.imap(work_wrapper, list_of_parameter_lists, chunksize=6), total=len(indices)))
    return model


def fit_kinetic_model(indices_here, do_plot=False,
                      p0=tuple([1]* 10)):

    def produce_fit(x, y):

        # # bounds for optimization:
        # lower_bounds = [-100, -100, # pKa_COH2plus_deltaS, pKa_COH2plus_deltaH
        #        -100*R_gas, -100*R_gas, #K_tilde_deltaS, K_tilde_deltaH
        #        -100, -100, # k_forward_over_k_backward_deltaS, k_forward_over_k_backward_deltaH
        #        -1e-3, -100, -100,  # k_backward_kappa, k_backward_deltaS, k_backward_deltaH
        #        -200 # HBr_B
        #         ]
        # upper_bounds = [-1*lower_bound for lower_bound in lower_bounds]

        # bounds for bootstrapping:
        p_default = [ 4.08008635e+01, -4.04065529e+00,  2.68463787e+02,  8.08157922e+01, 1.90497493e+01, -5.60271559e+00,  5.38416819e-06, -3.46293065e+00,  2.71691036e+00,  1.63131631e+02]
        lower_bounds_factor = 0.8
        upper_bounds_factor = 1.2
        # if p_default is negative, lower bound is p_default * upper_bounds_factor, upper bound is p_default * lower_bounds_factor
        # if p_default is positive, lower bound is p_default * lower_bounds_factor, upper bound is p_default * upper_bounds_factor
        lower_bounds = [p_default[i] * upper_bounds_factor if p_default[i] < 0 else p_default[i] * lower_bounds_factor for i in range(len(p_default))]
        upper_bounds = [p_default[i] * lower_bounds_factor if p_default[i] < 0 else p_default[i] * upper_bounds_factor for i in range(len(p_default))]
        # lower_bounds[0] = p_default[0] * 0.95
        # upper_bounds[0] = p_default[0] * 1.05
        # lower_bounds[1] = p_default[1] * 1.05
        # upper_bounds[1] = p_default[1] * 0.95
        # lower_bounds[-3] = p_default[-3] * 1.1
        # upper_bounds[-3] = p_default[-3] * 0.9
        # lower_bounds[-2] = p_default[-2] * 0.9
        # upper_bounds[-2] = p_default[-2] * 1.1
        popt, pcov = curve_fit(model_of_yield_for_many_conditions, x, y, p0=p0,
                               max_nfev=1000, bounds=(lower_bounds, upper_bounds), verbose=2,
                               x_scale=[x for x in [45, 10, 250, 80, 20, 5, 1e-6, 5, 3, 180]],
                               diff_step=[x/10 for x in ([0.02] + [0.03] + [0.03] + [0.03] + [0.01] + [0.01] + [0.01] + [0.2] + [0.2] + [0.1])],
                               loss='soft_l1', f_scale=0.1, sigma=0.04*np.ones_like(y), absolute_sigma=True, jac='3-point', # sigma=0.04*y + 0.0001
                               xtol=None, ftol=1e-6)
        best_f = lambda x: model_of_yield_for_many_conditions(x, *popt)
        print(f'popt = {popt}')
        perr = np.sqrt(np.diag(pcov))
        # print the rmse of the residuals
        residuals = y - best_f(x)
        rmse = np.sqrt(np.mean(residuals**2))
        print(f'RMSE = {rmse}')
        return best_f, popt, perr

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
    temperatures = df_results['temperature'].unique()
    temperatures = np.sort(temperatures)

    keq_initials = np.load('E1_3temps/keq_fits_with_errs_tdep_homosc_2024-12-16a_temps16-21-26__from_medians_exact__size249_rep6.npy')

    # only indices where temperature is below 30
    mask = df_results['temperature'] < 30
    indices_where_mask_is_true = df_results[mask].index.to_numpy()

    p0 = keq_initials[0]

    ## make a bootstrapped sampling (with replacement) of the indices_where_mask_is_true
    N_bootstrapped_samples = 1000
    bootstrapped_keq_fits = []
    for i in range(N_bootstrapped_samples):
        print(f'Bootstrapped sample {i+1} of {N_bootstrapped_samples}')
        indices_where_mask_is_true_bootstrapped = np.random.choice(indices_where_mask_is_true, size=len(indices_where_mask_is_true),
                                                                   replace=True)
        # multiply elements of p0 by random number uniformly distributed between 0.85 and 1.15
        p0_here = np.copy(p0)
        p0_here[2] = p0[2] * np.random.uniform(0.85, 1.15)
        p0_here[3] = p0[3] * np.random.uniform(0.85, 1.15)
        print(f'p0_here = {p0_here}')
        keq_fit, keq_err = fit_kinetic_model(indices_where_mask_is_true_bootstrapped, do_plot=False, p0=p0_here)
        bootstrapped_keq_fits.append(keq_fit)
        # pickle bootstrapped_keq_fits
        np.save(f'bootstrapped_keq_fits_size{i+1}_2024-12-24a_temps16-21-26.npy', np.array(bootstrapped_keq_fits))

    print(f'keq_fit = {keq_fit}')
    print(f'keq_err = {keq_err}')

    # # save keq_fits as numpy array
    # output_fname = 'keq_fits_with_errs_tdep_homosc_2024-12-16a_temps16-21-26__from_medians_exact__size249_rep6.npy'
    # if 'ROBOCHEM_DATA_PATH' in os.environ:
    #     np.save(f'{data_folder}simple-reactions/2023-09-14-run01/results/kinetics/{output_fname}', np.array((keq_fit, keq_err)))
    # else:
    #     np.save(output_fname, np.array((keq_fit, keq_err)))

    for temp_index, temperature in enumerate(temperatures):
        fig1 = plt.figure(figsize=(4, 3.9), dpi=300)
        print(f'Plotting model at temperature = {temperature}')
        mask = (df_results['temperature'] == temperature)
        indices_where_mask_is_true = df_results[mask].index.to_numpy()
        plt.title(f'Temperature {temperature} °C')
        plot_kinetic_model(indices_where_mask_is_true, keq_fit)
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

    print(np.array(keq_fit))


