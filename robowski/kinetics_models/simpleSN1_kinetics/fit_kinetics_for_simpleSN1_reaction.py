from robowski.settings import *

import logging
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit, brentq
from scipy import interpolate

import robowski.misc_scripts.organize_run_results as organize_run_results
import robowski.visualize_results.animated_viewer_static as avs
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

dioxane_density = 1.034  # g/mL
dioxane_molar_mass = 88.11  # g/mol
water_molar_mass = 18.01528  # g/mol
hbr_molar_mass = 80.9119  # g/mol
acetic_acid_molar_mass = 60.05  # g/mol


def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


def water_activity_function(do_plot=False):
    data = np.loadtxt(repo_data_path + 'misc_scripts/activity_data/water-dioxane/'
                      'kogan-fridman-kafarov-1966/data_activity_of_water_in_dioxane.txt',
                      skiprows=1, delimiter=',')
    mole_fraction = data[:, 0]
    activity = data[:, 1]

    def custom_fit_func_3(x, a, b, c, d, e):
        return a + b * np.exp(-1 * c * x) + d * x + e * x ** 2

    def produce_fit(x, y):
        popt, pcov = curve_fit(custom_fit_func_3, x, y, p0=[1, -1, 1, 1, 1], maxfev=100000)
        best_f = lambda x: custom_fit_func_3(x, *popt)
        return best_f

    best_f = produce_fit(mole_fraction, activity)
    xs = np.linspace(0, 1, 100)
    ys_fit = best_f(xs)
    if do_plot:
        plt.scatter(mole_fraction, activity)
        plt.xlabel('Water mole fraction in dioxane-water mixture')
        plt.ylabel('Water activity')
        plt.plot(xs, ys_fit)
        plt.show()
    return best_f


# Third batch from december
experiment_name = 'simple-reactions/2023-12-11-run01/'
list_of_runs = tuple([
    '2023-12-11-run01',
    '2023-12-11-run02',
    '2023-12-12-run01',
    '2023-12-12-run02',
    '2023-12-16-run01',
    '2023-12-16-run02'])


# column_to_plot = 'HBr_relative_change'
column_to_plot = 'yield'

substances = ['c#SN1OH03', 'c#HBr', 'temperature']
substance_titles = ['Alcohol', 'HBr', 'Temperature']
substrates = ['c#SN1OH03', 'c#HBr']

df_results = organize_run_results.join_data_from_runs([f'simple-reactions/{x}/' for x in list_of_runs],
                                 round_on_columns=substances)

for i, row in df_results.iterrows():
    df_results.loc[i, 'c#H2O'] = row['c#HBr'] / 4.5 * 21.88934517
    df_results.loc[i, 'c#acetic_acid'] = row['c#HBr'] / 4.5 * 8.583416744
    df_results.loc[i, 'product_sum'] = row['pc#SN1OH03'] + row['pc#SN1Br03']

# # save dataframe into the summary of reaction yield data
# df_results.to_csv(f'summary_of_reaction_yield_data/SN1_simple/raw_yields.csv', index=False)

water_activity = water_activity_function()

def model_of_yield_for_one_condition(index_in_df, target_equilibrium_constant):
    # get row of df_results with index_in_df
    row = df_results.loc[index_in_df]
    # initial concentrations of substrates, HBR, water
    c_SN1OH03 = row['c#SN1OH03']
    c_HBr = row['c#HBr']
    c_H2O = row['c#H2O']
    logging.debug(f'c_SN1OH03 {c_SN1OH03:.2e}, c_HBr {c_HBr:.2e}, c_H2O {c_H2O:.2e}')
    # if one of concentrations is zero, return zero yield
    if c_SN1OH03 == 0 or c_HBr == 0:
        return 0

    dioxane_volume = row['vol#Dioxane'] + row['vol#SN1OH03']  # in microliters
    dioxane_moles = dioxane_volume / 1000 * dioxane_density / dioxane_molar_mass
    dioxane_molar_concentration = dioxane_moles / (500e-6)
    solution_density_g_per_L = row['c#H2O'] * water_molar_mass + \
                               dioxane_molar_concentration * dioxane_molar_mass + \
                               row['c#HBr'] * hbr_molar_mass + \
                               row['c#acetic_acid'] * acetic_acid_molar_mass
    dioxane_mass_fraction = dioxane_molar_concentration * dioxane_molar_mass / (solution_density_g_per_L)

    def eq_const_function(c_prod):
        # resulting concentrations of substrates, HBR, water
        c_SN1OH03_result = c_SN1OH03 - c_prod
        c_HBr_result = c_HBr - c_prod
        c_H2O_result = c_H2O + c_prod

        # activities of substrates, HBR, water
        water_molar_fraction = c_H2O_result / (c_H2O_result + dioxane_molar_concentration)
        water_activity_here = water_activity(water_molar_fraction)
        maximum_possible_water_molar_concetration = 55.56
        effective_water_molar_concentration = water_activity_here * c_H2O_result

        c_HBr_result_molality = c_HBr_result / solution_density_g_per_L * 1000
        HBr_activity_coefficient = 1
        return (c_HBr_result * HBr_activity_coefficient * c_SN1OH03_result) / (c_prod * effective_water_molar_concentration)

    # solve the equation eq_const_function(c_prod) = target_equilibrium_constant
    c_prod = brentq(lambda x: eq_const_function(x) - target_equilibrium_constant, a=0, b=min(c_SN1OH03, c_HBr))
    reaction_yield = c_prod / c_SN1OH03
    return reaction_yield

def model_of_yield_for_many_conditions(indices, target_equilibrium_constant):
    return [model_of_yield_for_one_condition(index_in_df=i,
                                             target_equilibrium_constant=target_equilibrium_constant)
            for i in indices]

predicted_yield = model_of_yield_for_one_condition(index_in_df=0, target_equilibrium_constant=1e-5)
print(f'predicted yield {predicted_yield:.2f}')

def smooth_across_HBr_concentrations(x, y, fraction_of_outliers=0.15, do_plot=True):

    def custom_fit_func(x, a, b, c, d):
        return a + b * np.exp(-1*c*x) + d*x

    def produce_fit(x, y):
        popt, pcov = curve_fit(custom_fit_func, x, y, p0=[1, -1, 4, 0], maxfev=100000)
        best_f = lambda x: custom_fit_func(x, *popt)
        return best_f

    x = np.array(x)
    y = np.array(y)
    if do_plot:
        plt.scatter(x, y)
    f = produce_fit(x, y)
    if do_plot:
        plt.plot(np.sort(x), f(np.sort(x)), '--', color='C1')
    # remove the 10% of the points furthest from the fit
    diff = np.abs(f(x) - y)
    indices_to_keep = np.argsort(diff)[:-int(len(diff) * fraction_of_outliers)]
    # if point with lowest x is not in the indices to keep, add it
    if np.argmin(x) not in indices_to_keep:
        indices_to_keep = np.append(indices_to_keep, np.argmin(x))
    x2 = x[indices_to_keep]
    y2 = y[indices_to_keep]
    if do_plot:
        plt.scatter(x2, y2, color='C2', marker='x')
    # fit polynomial again
    f = produce_fit(x2, y2)
    if do_plot:
        plt.plot(np.sort(x), f(np.sort(x)), color='C1')
        # plt.ylim(0, 1)
        plt.xlabel('Concentration of HBr')
        plt.ylabel('Allegedly, the equilibrium constant')
        plt.show()
    return f(x)

df_results.dropna(subset=[column_to_plot], inplace=True)
df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]

# round concentrations "c#SN1OH03" and "c#HBr" to 6 decimal places
df_results['c#SN1OH03'] = df_results['c#SN1OH03'].round(6)
df_results['c#HBr'] = df_results['c#HBr'].round(6)


def fit_kinetic_model(indices_here, do_plot=False):

    def produce_fit(x, y):
        popt, pcov = curve_fit(model_of_yield_for_many_conditions, x, y, p0=[1e-2],
                               max_nfev=100000, bounds=([0], [np.inf]),
                               loss='soft_l1', f_scale=0.05)
        best_f = lambda x: model_of_yield_for_many_conditions(x, *popt)
        return best_f, popt[0]

    # make unique values of the alcohol concentrations for indices_here
    unique_alcohol_concentrations = df_results.loc[indices_here, 'c#SN1OH03'].unique()
    # sort it
    unique_alcohol_concentrations = np.sort(unique_alcohol_concentrations)
    # make a list of colors based on the id of the alcohol concentration
    colors = [f'C{i}' for i in range(len(unique_alcohol_concentrations))]
    colors_to_plot = df_results.loc[indices_here, 'c#SN1OH03'].apply(lambda x: colors[np.where(unique_alcohol_concentrations == x)[0][0]])
    xs_to_plot = df_results.loc[indices_here, 'c#HBr']
    measured_yields = df_results.loc[indices_here, 'yield']
    # plt.scatter(xs_to_plot, measured_yields, color='yellow', marker='o')
    f, keq_fit = produce_fit(indices_here, measured_yields)
    if do_plot:
        plt.scatter(xs_to_plot, measured_yields, s=10, color=colors_to_plot, alpha=0.5)
        for c_alc in unique_alcohol_concentrations:
            color_here = colors[np.where(unique_alcohol_concentrations == c_alc)[0][0]]
            # find df_indices among indices_here where alcolhol concentration is c_alc
            indices_where_mask_is_true = df_results.loc[indices_here, 'c#SN1OH03'] == c_alc
            # sort indices by HBr concentration
            xs_here = df_results.loc[indices_here[indices_where_mask_is_true], 'c#HBr']
            ys_here = f(indices_here[indices_where_mask_is_true])
            # sort xs and ys by increasing xs
            xs_here, ys_here = zip(*sorted(zip(xs_here, ys_here)))
            plt.plot(xs_here, ys_here, color=color_here, label=f'{c_alc:.3f} M')
        # plt.scatter(xs_to_plot, ys_to_plot, color=colors_to_plot, marker='x')
        plt.ylabel('Yield')
        plt.xlabel('Initial concentration of HBr')
        plt.legend(title="Starting alcohol\nconcentration")

    return keq_fit

do_plot = True
keq_fits = []
temperatures = df_results['temperature'].unique()
temperatures = np.sort(temperatures)
print(f'Temperature values: {temperatures}')
for temperature in temperatures:
    fig1 = plt.figure(figsize=(4, 3.9))
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
        plt.show()

xs = 1000/(273.15 + temperatures)
ys = -1*np.log(keq_fits)

np.savetxt(f'{data_folder}{experiment_name}results/kinetics/1000_over_t.txt', xs)
np.savetxt(f'{data_folder}{experiment_name}results/kinetics/logK.txt', ys)