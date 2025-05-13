import importlib
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.signal import savgol_filter

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

def get_derivs_3d(xs, ys, zs, ks):
    unique_xs = np.sort(np.unique(xs))
    unique_ys = np.sort(np.unique(ys))
    unique_zs = np.sort(np.unique(zs))

    numunique_xs = len(unique_xs)
    numunique_ys = len(unique_ys)
    numunique_zs = len(unique_zs)

    # # interpolate to a cube with unique values of xs, ys, zs
    # xnew, ynew, znew = np.meshgrid(unique_xs, unique_ys, unique_zs, indexing='ij')

    # mgrid with the same number of points as the original data, but the points are evenly spaced between max unique and min unique
    xnew, ynew, znew = np.mgrid[unique_xs[0]:unique_xs[-1]:numunique_xs * 1j,
                       unique_ys[0]:unique_ys[-1]:numunique_ys * 1j,
                       unique_zs[0]:unique_zs[-1]:numunique_zs * 1j]
    interp_here = LinearNDInterpolator((xs, ys, zs), ks)
    wnew = interp_here((xnew, ynew, znew))

    # for every dimension, get the partial derivative with respect to that dimension at all points of the grid.
    # Use 3-point scheme.

    # differentiating with respect to x
    deriv_cube_x = np.zeros_like(wnew)
    # dx is the spacing between x values
    dx = xnew[1, 0, 0] - xnew[0, 0, 0]
    # for every combination of y and z, get the array of wnew values
    for i in range(numunique_ys):
        for j in range(numunique_zs):
            wnew_slice = wnew[:, i, j]
            # use savgol filter with derivative order 1 and window size 3
            wnew_slice_deriv = savgol_filter(wnew_slice, 3, 2, deriv=1, delta=dx, mode='interp')
            deriv_cube_x[:, i, j] = wnew_slice_deriv

    # differentiating with respect to y
    deriv_cube_y = np.zeros_like(wnew)
    dy = ynew[0, 1, 0] - ynew[0, 0, 0]
    for i in range(numunique_xs):
        for j in range(numunique_zs):
            wnew_slice = wnew[i, :, j]
            wnew_slice_deriv = savgol_filter(wnew_slice, 3, 2, deriv=1, delta=dy, mode='interp')
            deriv_cube_y[i, :, j] = wnew_slice_deriv

    # differentiating with respect to z
    deriv_cube_z = np.zeros_like(wnew)
    dz = znew[0, 0, 1] - znew[0, 0, 0]
    for i in range(numunique_xs):
        for j in range(numunique_ys):
            wnew_slice = wnew[i, j, :]
            wnew_slice_deriv = savgol_filter(wnew_slice, 3, 2, deriv=1, delta=dz, mode='interp')
            deriv_cube_z[i, j, :] = wnew_slice_deriv

    # concatenate flattened deriv cubes
    derivs = np.concatenate([deriv_cube_x.flatten(), deriv_cube_y.flatten(), deriv_cube_z.flatten()])

    return derivs


def get_derivs_2d(xs, ys, ks):
    unique_xs = np.sort(np.unique(xs))
    unique_ys = np.sort(np.unique(ys))

    numunique_xs = len(unique_xs)
    numunique_ys = len(unique_ys)

    # # interpolate to a cube with unique values of xs, ys, zs
    # xnew, ynew, znew = np.meshgrid(unique_xs, unique_ys, unique_zs, indexing='ij')

    # mgrid with the same number of points as the original data, but the points are evenly spaced between max unique and min unique
    xnew, ynew = np.mgrid[unique_xs[0]:unique_xs[-1]:numunique_xs * 1j,
                       unique_ys[0]:unique_ys[-1]:numunique_ys * 1j]
    interp_here = LinearNDInterpolator((xs, ys), ks)
    wnew = interp_here((xnew, ynew))

    # for every dimension, get the partial derivative with respect to that dimension at all points of the grid.
    # Use 3-point scheme.

    # differentiating with respect to x
    deriv_cube_x = np.zeros_like(wnew)
    # dx is the spacing between x values
    dx = xnew[1, 0] - xnew[0, 0]
    # for every combination of y, get the array of wnew values
    for i in range(numunique_ys):
            wnew_slice = wnew[:, i]
            # use savgol filter with derivative order 1 and window size 3
            wnew_slice_deriv = savgol_filter(wnew_slice, 3, 2, deriv=1, delta=dx, mode='interp')
            deriv_cube_x[:, i] = wnew_slice_deriv

    # differentiating with respect to y
    deriv_cube_y = np.zeros_like(wnew)
    dy = ynew[0, 1] - ynew[0, 0]
    for i in range(numunique_xs):
            wnew_slice = wnew[i, :]
            wnew_slice_deriv = savgol_filter(wnew_slice, 3, 2, deriv=1, delta=dy, mode='interp')
            deriv_cube_y[i, :] = wnew_slice_deriv

    # concatenate flattened deriv cubes
    derivs = np.concatenate([deriv_cube_x.flatten(), deriv_cube_y.flatten()])

    print(f'dx = {dx}, dy = {dy}')

    return derivs

def make_hist_fot_ugi():
    # what_to_plot = 'model'
    what_to_plot = 'data'

    data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

    experiment_name = 'multicomp-reactions/2023-06-19-run01/'
    if what_to_plot == 'model':
        df_results_all = pd.read_csv(data_folder + experiment_name + f'results/interpolated_product_concentration_model.csv')
    elif what_to_plot == 'data':
        df_results_all = pd.read_csv(data_folder + experiment_name + f'results/interpolated_product_concentration.csv')

    sorted_unique_values_of_ptsa_column = df_results_all['ptsa'].unique()
    sorted_unique_values_of_ptsa_column.sort()
    print(sorted_unique_values_of_ptsa_column)

    # substrates = ['c#SN1OH03', 'c#HBr']
    substances = ['am001','ic001','ald001']
    substance_titles = ['Amine', 'Isocyanide', 'Aldehyde']

    # ptsa_targets = [0.05, 0.12, 0.17, 0.24, 0.298]
    # ptsa_target = ptsa_targets[4]

    def get_derivs_for_one_ptsa(df_results_all, ptsa_target):
        # first index where ptsa is greater than
        ith_ptsa = next(i for i, x in enumerate(sorted_unique_values_of_ptsa_column) if x > ptsa_target)
        # ith_ptsa = len(sorted_unique_values_of_ptsa_column) - 1

        print(f'i={ith_ptsa}, PTSA: {sorted_unique_values_of_ptsa_column[ith_ptsa]}')
        df_results = df_results_all[df_results_all['ptsa'] == sorted_unique_values_of_ptsa_column[ith_ptsa]]
        # drop nans in the yield column
        df_results.dropna(subset=['yield'], inplace=True)
        column_to_plot = 'yield'

        for substance in substances:
            df_results[substance] = df_results[substance].round(6).astype(np.float64)

        xs = df_results[substances[0]].to_numpy()
        ys = df_results[substances[1]].to_numpy()
        zs = df_results[substances[2]].to_numpy()
        yields = df_results[column_to_plot].to_numpy() # + df_results['pc#dm40_10'].to_numpy()
        # for every data point, find the smallest of the xs, ys, za. The length of the array of limiting reactants is the same as lengths of xs, ys, zs
        limiting_reactant_c = np.zeros_like(xs)
        for i in range(len(xs)):
            limiting_reactant_c[i] = min(xs[i], ys[i], zs[i])
        product_concentrations = yields * limiting_reactant_c

        derivs = get_derivs_3d(xs, ys, zs, product_concentrations)
        return derivs

    all_derivs_list = []
    derivs = get_derivs_for_one_ptsa(df_results_all, ptsa_target=0)
    all_derivs_list.append(derivs)
    for ptsa_target in sorted_unique_values_of_ptsa_column[:-1]:
        derivs = get_derivs_for_one_ptsa(df_results_all, ptsa_target)
        all_derivs_list.append(derivs)

    all_derivs_list = np.concatenate(all_derivs_list)

    print('Max and median:')
    print(f'{np.max(np.abs(all_derivs_list)):.3f}\t{np.median(np.abs(all_derivs_list)):.3f}')

    output_filename = 'experimental_hist_Ugi'
    f1 = plt.figure(figsize=(5, 1.8))
    # set right margin to zero
    plt.subplots_adjust(right=0.99)
    # set bottom margin to larger value
    plt.subplots_adjust(bottom=0.25)
    plt.hist(np.abs(all_derivs_list), bins=100, density=True, color='grey')
    plt.axvline(np.max(np.abs(all_derivs_list)), color='red')
    plt.xlabel('Partial derivative')
    plt.ylabel('Probability density')
    # plt.tight_layout()
    f1.savefig(f'misc_scripts/smoothness_theory/figures/{output_filename}.png', dpi=300)
    plt.show()

def make_hist_for_e1():
    data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
    df_results = pd.read_csv('visualize_results/examples/final_datasets/e1results.csv')

    temps = np.array([16, 21, 26, 31, 36])
    all_derivs_list = []
    for temperature in temps:
        df_one_temp = df_results[df_results['temperature'] == temperature]
        xs = df_one_temp['c#E1OH02'].to_numpy()
        ys = df_one_temp['c#HBr'].to_numpy()
        ks = df_one_temp['pc#E1DB02'].to_numpy()
        derivs = get_derivs_2d(xs, ys, ks)
        all_derivs_list.append(derivs)

    all_derivs_list = np.concatenate(all_derivs_list)

    print('Max and median:')
    print(f'{np.max(np.abs(all_derivs_list)):.3f}\t{np.median(np.abs(all_derivs_list)):.3f}')

    output_filename = 'experimental_hist_E1'
    f1 = plt.figure(figsize=(5, 1.8))
    # set right margin to zero
    plt.subplots_adjust(right=0.99)
    # set bottom margin to larger value
    plt.subplots_adjust(bottom=0.25)
    plt.hist(np.abs(all_derivs_list), bins=100, density=True, color='grey')
    plt.axvline(np.max(np.abs(all_derivs_list)), color='red')

    # # plot noise.
    # absolute_noise = 0.001 # 1 mM at the highest alcohol concentration
    # dx = 0.003375
    # dy = 0.0003333333333333333
    # # multiplication by sqrt(2) is because the derivative has subtraction in the numerator,
    # # and the variance of the difference of two independent random variables is the sum of their variances
    # plt.axvline(absolute_noise/dx * np.sqrt(2), color='blue')

    plt.xlabel('Partial derivative')
    plt.ylabel('Probability density')
    # plt.tight_layout()
    f1.savefig(f'misc_scripts/smoothness_theory/figures/{output_filename}.png', dpi=300)
    plt.show()


def make_hist_for_sn1():
    organize_run_results = importlib.import_module("misc_scripts.organize_run_results")
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

    df_results.dropna(subset=[column_to_plot], inplace=True)
    df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]

    # round concentrations "c#SN1OH03" and "c#HBr" to 6 decimal places
    df_results['c#SN1OH03'] = df_results['c#SN1OH03'].round(6)
    df_results['c#HBr'] = df_results['c#HBr'].round(6)

    column_to_plot = 'yield'
    # convert from mol/L to mM

    # iterate over rows and compute column yield2
    for i, row in df_results.iterrows():
        limiting_substrate = min(row['pc#SN1OH03'] + row['pc#SN1Br03'], row['c#HBr'] / row['c#SN1OH03'] * (row['pc#SN1OH03'] + row['pc#SN1Br03']))
        # limiting_substrate = min(row['pc#SN1OH03'] + row['pc#SN1Br03'], row['c#HBr'])
        # limiting_substrate = min(row['c#SN1OH03'], row['c#HBr'])
        if limiting_substrate == 0:
            df_results.loc[i, 'yield2'] = 0
        else:
            df_results.loc[i, 'yield2'] = row['pc#SN1Br03'] / limiting_substrate

    column_to_plot = 'yield2'

    column_to_plot = 'pc#SN1Br03'
    temps = df_results['temperature'].unique()
    all_derivs_list = []
    for temperature in temps:
        df_one_temp = df_results[df_results['temperature'] == temperature]
        xs = df_one_temp['c#SN1OH03'].to_numpy()
        ys = df_one_temp['c#HBr'].to_numpy()
        ks = df_one_temp[column_to_plot].to_numpy()
        derivs = get_derivs_2d(xs, ys, ks)
        all_derivs_list.append(derivs)

    all_derivs_list = np.concatenate(all_derivs_list)

    print('Max and median:')
    print(f'{np.max(np.abs(all_derivs_list)):.3f}\t{np.median(np.abs(all_derivs_list)):.3f}')

    output_filename = 'experimental_hist_SN1'
    f1 = plt.figure(figsize=(5, 1.8))
    # set right margin to zero
    plt.subplots_adjust(right=0.99)
    # set bottom margin to larger value
    plt.subplots_adjust(bottom=0.25)
    plt.hist(np.abs(all_derivs_list), bins=100, density=True, color='grey')
    plt.axvline(np.max(np.abs(all_derivs_list)), color='red')

    plt.xlabel('Partial derivative')
    plt.ylabel('Probability density')
    # plt.tight_layout()
    f1.savefig(f'misc_scripts/smoothness_theory/figures/{output_filename}.png', dpi=300)
    plt.show()


def make_hist_for_sn1_carbocat():
    organize_run_results = importlib.import_module("misc_scripts.organize_run_results")
    # timepoint_id = 1
    experiment_name = 'simple-reactions/2023-07-05-run01/'
    # df_results = pd.read_csv(data_folder + experiment_name + f'results/timepoint{timepoint_id:03d}-reaction_results.csv')

    list_of_runs = tuple([
        '2023-07-05-run01',
        '2023-07-06-run01',
        '2023-07-07-run01',
        '2023-07-10-run01',
        '2023-07-10-run02',
        '2023-07-11-run01',
        '2023-07-11-run02'])

    df_results = organize_run_results.join_data_from_runs([f'simple-reactions/{x}/' for x in list_of_runs],
                                     round_on_columns=('c#SN1OH01', 'c#HBr', 'Temperature'))

    # set yields to nan at the minimum of c#SN1OH01 column
    df_results.loc[df_results['c#SN1OH01'].round(4) == df_results['c#SN1OH01'].round(4).min(), 'yield'] = np.nan

    # set yields to nan where the temperature is 6 and yield is above 0.9
    df_results.loc[(df_results['Temperature'] == 6) & (df_results['yield'] > 0.9), 'yield'] = np.nan

    # remove all rows where 'yield' columns is nan
    df_results.dropna(subset=['yield'], inplace=True)

    df_results['yield'] = df_results['yield'].apply(lambda x: x if x>1e-10 else 0)
    # df_results.drop(indices_of_outliers, inplace=True)

    substances = ['c#SN1OH01', 'c#HBr', 'Temperature']
    substance_titles = ['Alcohol', 'HBr', 'Temperature']
    substrates = ['c#SN1OH01', 'c#HBr']

    for substance in substances:
        df_results[substance] = df_results[substance].round(6).astype(np.float64)

    temps = df_results['Temperature'].unique()
    all_derivs_list = []
    for yieldname in ['yield#SN1Br01', 'yield#SN1Br01s1']:
        for temperature in temps:
            df_one_temp = df_results[df_results['Temperature'] == temperature]
            xs = df_one_temp['c#SN1OH01'].to_numpy()
            ys = df_one_temp['c#HBr'].to_numpy()
            ks = xs * df_one_temp[yieldname].to_numpy()

            derivs = get_derivs_2d(xs, ys, ks)
            all_derivs_list.append(derivs)

    all_derivs_list = np.concatenate(all_derivs_list)

    print('Max and median:')
    print(f'{np.max(np.abs(all_derivs_list)):.3f}\t{np.median(np.abs(all_derivs_list)):.3f}')

    output_filename = 'experimental_hist_SN1carbocat'
    f1 = plt.figure(figsize=(5, 1.8))
    # set right margin to zero
    plt.subplots_adjust(right=0.99)
    # set bottom margin to larger value
    plt.subplots_adjust(bottom=0.25)
    plt.hist(np.abs(all_derivs_list), bins=100, density=True, color='grey')
    plt.axvline(np.max(np.abs(all_derivs_list)), color='red')

    plt.xlabel('Partial derivative')
    plt.ylabel('Probability density')
    # plt.tight_layout()
    f1.savefig(f'misc_scripts/smoothness_theory/figures/{output_filename}.png', dpi=300)
    plt.show()

if __name__ == '__main__':
    pass
    # make_hist_fot_ugi()
    # make_hist_for_e1()
    # make_hist_for_sn1()
    # make_hist_for_sn1_carbocat()