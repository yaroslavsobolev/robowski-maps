import importlib
import os
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy import interpolate
from scipy.optimize import curve_fit
organize_run_results = importlib.import_module("misc_scripts.organize_run_results")
avs = importlib.import_module("visualize_results.animated_viewer_static")
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'


# experiment_name = 'simple-reactions/2023-08-21-run01/'
# list_of_runs = tuple([
#     '2023-08-21-run01',
#     '2023-08-22-run01',
#     '2023-08-22-run02',
#     '2023-08-28-run01',
#     '2023-08-29-run01',
#     '2023-08-29-run02'])
# column_to_plot = 'HBr_relative_change'
# column_to_plot = 'conversion'
# substances = ['c#SN1OH03', 'c#HBr', 'temperature']
# df_results = organize_run_results.join_data_from_runs([f'simple-reactions/{x}/' for x in list_of_runs],
#                                  round_on_columns=substances)
# for i, row in df_results.iterrows():
#     df_results.loc[i, 'c#H2O'] = row['c#HBr'] / 4.5 * 21.88934517
#     df_results.loc[i, 'c#acetic_acid'] = row['c#HBr'] / 4.5 * 8.583416744
#     df_results.loc[i, 'product_sum'] = row['pc#SN1OH03'] + row['pc#SN1Br03']
#
#
# for i, row in df_results.iterrows():
#     dioxane_molar_mass = 88.11 # g/mol
#     water_molar_mass = 18.01528 # g/mol
#     hbr_molar_mass = 80.9119 # g/mol
#     acetic_acid_molar_mass = 60.05 # g/mol
#     # concentration of product
#     C_RBr = row['pc#SN1Br03']/(row['pc#SN1OH03'] + row['pc#SN1Br03']) * row['c#SN1OH03']
#     # concentration of water at the end of reaction
#     C_H2O = row['c#H2O'] + C_RBr
#
#     # in this section we compute the activity of water in the reaction mixture
#     dioxane_density = 1.034 # g/mL
#     dioxane_volume = row['vol#Dioxane'] + row['vol#SN1OH03'] # in microliters
#     dioxane_moles = dioxane_volume / 1000 * dioxane_density / dioxane_molar_mass
#     dioxane_molar_concentration = dioxane_moles / (500e-6)
#     solution_density_g_per_L = row['c#H2O'] * water_molar_mass + \
#                                dioxane_molar_concentration * dioxane_molar_mass + \
#                                row['c#HBr'] * hbr_molar_mass + \
#                                row['c#acetic_acid'] * acetic_acid_molar_mass
#     dioxane_mass_fraction = dioxane_molar_concentration * dioxane_molar_mass / (solution_density_g_per_L)
#     df_results['dioxane_mass_fraction'] = dioxane_mass_fraction
#
#     C_HBr_remaining = row['c#HBr'] - C_RBr
#     C_HBr_remaining_molality = C_HBr_remaining / solution_density_g_per_L * 1000
#     df_results.loc[i, 'HBr_remaining_molality'] = C_HBr_remaining_molality
#     plt.scatter(dioxane_mass_fraction, C_HBr_remaining_molality, color='black')
#
# plt.xlabel('dioxane mass fraction')
# plt.ylabel('HBr molality')
# plt.show()

def make_interpolator_of_hbr_activity(do_plot = True):
    if do_plot:
        fig2 = plt.figure(2)
    df = pd.read_csv('misc_scripts/activity_data/hbr-water-dioxane/mussini-et-al-1971-electroanalchem.csv')

    interp_dict = dict()
    for i, row in df.iterrows():
        molality = row['molality']
        xs = [int(x) for x in df.columns[1:]]
        ys = np.array(row[1:])
        if do_plot:
            plt.plot(xs, ys, 'o', label=f'{molality:.3f} mol/kg')

        def custom_fit_func(x, a, b, c, d, e):
            return a + b/(1 + np.exp(-1 * c * (x - d))) + e*x

        def produce_fit(x, y):
            lower_bounds = [-np.inf] * 5
            upper_bounds = [np.inf] * 5
            upper_bounds[4] = 0
            upper_bounds[0] = 1
            lower_bounds[0] = 0
            # upper_bounds[1] = 2
            popt, pcov = curve_fit(custom_fit_func, x, y, p0=[0, 1, 1/20, 70, 0], bounds=(lower_bounds, upper_bounds),
                                   maxfev=100000)
            best_f = lambda x: custom_fit_func(x, *popt)
            return best_f

        best_f = produce_fit(xs, ys)
        interp_dict[molality] = best_f
        xs2 = np.linspace(0, 100, 100)
        ys_fit = best_f(xs2)
        if do_plot:
            # plt.scatter(xs, ys)
            plt.plot(xs2, ys_fit, color='grey', linestyle='--')
            plt.ylim(0, 1)
            # plt.show()
    plt.xlabel('Dioxane mass fraction, %')
    plt.ylabel('HBr molal activity coefficient')
    plt.legend()
    plt.show()

    dio_mass_fracts = np.linspace(0, 100, 100)
    hbr_molalities = np.linspace(0, 1.2, 500)

    xx, yy = np.meshgrid(dio_mass_fracts, hbr_molalities)
    zz = np.zeros_like(xx)

    for i, dio_mass_frac in enumerate(dio_mass_fracts):
        ms = []
        acts = []
        for molality in interp_dict.keys():
            best_f = interp_dict[molality]
            ms.append(molality)
            acts.append(best_f(dio_mass_frac))
        if do_plot:
            plt.plot(ms, acts, 'o', label=f'{dio_mass_frac:.2f}% dioxane')

        # def custom_fit_func_2(x, a, b, c):
        #     return a + b*np.exp(-1 * c * x)

        def custom_fit_func_2(x, a, b, c, d):
            return a*np.exp(-1*b*x) + c * np.exp(-1 * d * x)

        def produce_fit_2(x, y):
            lower_bounds = [-np.inf] * 3
            upper_bounds = [np.inf] * 3

            # upper_bounds[1] = 2
            # p0 = [0, 1, 1/0.02]
            # p0 = [0, 1/0.04, 0.02]
            p0 = (0.5, 1/0.04, 0.5, 1/0.5)
            popt, pcov = curve_fit(custom_fit_func_2, x, y, p0=p0, bounds=(lower_bounds, upper_bounds),
                                   maxfev=100000)
            best_f = lambda x: custom_fit_func_2(x, *popt)
            return best_f

        best_f = produce_fit_2(ms, acts)
        xs2 = np.linspace(0, 1, 500)
        ys_fit = best_f(xs2)
        if do_plot:
            # plt.scatter(xs, ys)
            plt.plot(xs2, ys_fit, color='grey', linestyle='--')
            plt.ylim(0, 1)
            # plt.show()

        # find the indices of zz where the dio_mass_frac is xx and the molality is yy
        # then set zz to the activity
        for j, hbr_molality in enumerate(hbr_molalities):
            assert xx[j, i] == dio_mass_frac
            assert yy[j, i] == hbr_molality
            zz[j, i] = best_f(hbr_molality)

    if do_plot:
        plt.xlabel('HBr molality')
        plt.ylabel('HBr molal activity coefficient')
        plt.show()


    # make an interpolator from xx, yy and zz using interp2d
    interp2d_here = interpolate.interp2d(dio_mass_fracts, hbr_molalities, zz, kind='cubic')

    if do_plot:
        # plor a pcolor plot of xx, yy, zz
        plt.figure(3)
        plt.pcolor(xx, yy, zz)
        plt.colorbar()
        plt.xlabel('Dioxane mass fraction')
        plt.ylabel('HBr molality')
        plt.show()

        # do a second plot, now with the interpolator
        hbr_molality = hbr_molalities[50]
        zs = [interp2d_here(x, hbr_molality) for x in dio_mass_fracts]
        plt.plot(dio_mass_fracts, zs, color='black')
            # plot original data
        plt.plot(dio_mass_fracts, zz[50, :], 'o', color='black')
        plt.xlabel('Dioxane mass fraction')
        plt.ylabel('HBr molal activity coefficient')
        plt.show()

    return interp2d_here

xxx = make_interpolator_of_hbr_activity(do_plot=True)


