import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import importlib
from tqdm import tqdm
import seaborn as sns
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib import pyplot as plt, ticker as mticker

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

date_dict = {
    'dm35_8': '2023-12-20',
    'dm35_9': '2023-12-20',
    'dm36': '2023-12-20',
    'dm37': '2023-12-20',
    'EAB': '2023-12-26',
    'dm40_12': '2024-01-16',
    'dm40_10': '2024-01-16',
    'dm053': '2024-01-16',
    'bb017': '2024-01-22',
    'bb021': '2024-02-02',
    'dm70': '2024-02-02',
    'dm088_4': '2024-03-14',
    'bb021_f2': '2024-04-03'
    }

calibrant_sets = {
    0: {'calibrant_set': ['methoxybenzaldehyde', 'HRP01', 'ethyl_acetoacetate'], 'date': '2023-11-08'},
    1: {'calibrant_set': ['dm35_8', 'dm35_9', 'dm36', 'dm37'], 'date': '2023-12-20'},
    2: {'calibrant_set': ['EAB'], 'date': '2023-12-26'},
    3: {'calibrant_set': ['dm40_12', 'dm40_10', 'dm053'], 'date': '2024-01-16'},
    4: {'calibrant_set': ['bb017'], 'date': '2024-01-22'},
    5: {'calibrant_set': ['bb021', 'dm70'], 'date': '2024-02-02'},
    6: {'calibrant_set': ['dm088_4'], 'date': '2024-03-14'},
    7: {'calibrant_set': ['bb021_f2'], 'date': '2024-04-03'}
}

process_wellplate_spectra = importlib.import_module("uv-vis-absorption-spectroscopy.process_wellplate_spectra")
organize_run_results = importlib.import_module("misc-scripts.organize_run_results")
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
plt.ioff()

import logging
logging.basicConfig(level=logging.INFO)


# substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'dm36', 'dm37', 'dm40_12', 'dm40_10', 'ethyl_acetoacetate', 'EAB', 'bb017', 'bb021', 'dm70', 'dm053', 'dm088_4', 'bb021_f2']

# substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'dm36', 'dm37', 'dm40_12', 'dm40_10',
#                           'ethyl_acetoacetate', 'bb017', 'bb021', 'dm70', 'dm053', 'dm088_4', 'bb021_f2']
# try removing: 35_8, 35_9, dm37, dm_70

# substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'ethyl_acetoacetate', 'EAB', 'bb017', 'dm088_4', 'dm053', 'dm70'] # should inlclude dm37 which is oxidized dm70

# ############### best benchmarks so far. This is what was used before 2024-10-01 for the 2024-04-17-run1-2 data
# substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'ethyl_acetoacetate', 'EAB', 'bb017', 'dm088_4', 'dm053', 'dm70']
#
# ############### what was used for the 2024-10-01 reprocessing of the 2024-04-17-run1-2 data
# substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'ethyl_acetoacetate', 'EAB', 'bb017',
#                           'dm088_4', 'dm053', 'dm70', 'dm40_12', 'dm40_10']

def process_run_by_shortname(run_shortname, substances_for_fitting, folder_for_prod_conc, cut_from=1):
    substrates = ['methoxybenzaldehyde', 'ethyl_acetoacetate', 'ammonium_acetate']
    product_name = 'HRP01'
    run_name = f'BPRF/{run_shortname}/'
    print(f'Processing run {run_name}...')

    sp = process_wellplate_spectra.SpectraProcessor(
        folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
        '2022-12-01/interpolator-dataset/')
    sp.nanodrop_lower_cutoff_of_wavelengths = 220
    sp.use_instrumental_sigmas = True

    df_structure = organize_run_results.load_run_structure(run_name)

    # unique filepaths in the column 'plate_barcodes_for_dilution'
    nanodrop_filepaths = df_structure['nanodrop_filepath'].unique()

    for filepath_1 in nanodrop_filepaths:
        filepath_for_first_dilution = data_folder + run_name + '/nanodrop_spectra/' + filepath_1
        filepath_for_second_dilution =  data_folder + run_name + '/nanodrop_spectra/' + \
                                        df_structure[df_structure['nanodrop_filepath'] == filepath_1].iloc[0]['nanodrop_filepath_2']
        logging.info(f'Processing nanodrop file (1st dilution): {filepath_1} and (2nd dilution): {filepath_for_second_dilution}')

        # number of vials in this plate
        number_of_vials = df_structure[df_structure['nanodrop_filepath'] == filepath_1].shape[0]
        list_of_starting_concentration_dicts = []
        for vial_id in range(number_of_vials):
            index_of_this_vial = df_structure.index[(df_structure['container_id'] == vial_id)
                                                    & (df_structure['nanodrop_filepath'] == filepath_1)][0]
            print({s: df_structure.loc[index_of_this_vial, f'c#{s}'] for s in substrates})
            list_of_starting_concentration_dicts.append({s: float(df_structure.loc[index_of_this_vial, f'c#{s}']) for s in substrates})

        concentrations_here, reports = sp.multispectrum_concentrations_for_one_plate(experiment_folder=data_folder + run_name,
                                                                                     plate_folder_1=filepath_for_first_dilution,
                                                                                     plate_folder_2=filepath_for_second_dilution,
                                                                                     dilution_factors=[20, 200],
                                                                                     calibration_folder=data_folder + 'BPRF/2024-01-17-run01/' + 'microspectrometer_data/calibration/',
                                                                                     calibrant_shortnames=substances_for_fitting,
                                                                                     background_model_folder=data_folder + 'BPRF/cross_conamination_and_backgound_test/ethanol_background_model/',
                                                                                     calibrant_upper_bounds=[np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf],
                                                                                     do_plot=False, cut_from=cut_from, cut_to=350,
                                                                                     ignore_abs_threshold=False, ignore_pca_bkg=False,
                                                                                     return_all_substances=True,
                                                                                     upper_limit_of_absorbance=0.95,
                                                                                     return_report=True,
                                                                                     list_of_starting_concentration_dicts=list_of_starting_concentration_dicts,
                                                                                     obey_stoichiometric_inequalities=True)

        for vial_id, product_concentrations in enumerate(concentrations_here):
            # index of this vial in the concentrations_df dataframe
            try:
                index_of_this_vial = df_structure.index[(df_structure['container_id'] == vial_id)
                                                        & (df_structure['nanodrop_filepath'] == filepath_1)][0]
            except IndexError:
                logging.warning(f'Container_id {vial_id} not present in conditions dataframe but was measured in nanodrop file {filepath_1}, '
                                f'second dilution at {df_structure[df_structure["nanodrop_filepath"] == filepath_1].iloc[0]["nanodrop_filepath_2"]}')
                continue

            for i, substance_name in enumerate(substances_for_fitting):
                df_structure.loc[index_of_this_vial, f'pc#{substance_name}'] = product_concentrations[i]

            for key in reports[vial_id]:
                df_structure.loc[index_of_this_vial, key] = reports[vial_id][key]

            required_subs = process_wellplate_spectra.product_concentrations_to_required_substrates(product_concentrations, substances_for_fitting)
            for substrate in substrates:
                df_structure.loc[index_of_this_vial, f'req#{substrate}'] = required_subs[substrate]
                df_structure.loc[index_of_this_vial, f'os#{substrate}'] = (required_subs[substrate] - float(df_structure.loc[index_of_this_vial, f'c#{substrate}'])) / float(df_structure.loc[index_of_this_vial, f'c#{substrate}']) * 100

    logging.debug(f'global id min {df_structure["global_index"].min()}')
    concentrations_df = df_structure[df_structure[f'pc#{product_name}'].notna()]
    logging.debug(f'cdf global id min {concentrations_df["global_index"].min()}')

    # Calculate yields
    for index, row in concentrations_df.iterrows():
        # HRP01
        product_name = 'HRP01'
        product_concentration = concentrations_df.loc[index, f'pc#{product_name}']
        product_error = concentrations_df.loc[index, f'pcerr#{product_name}']
        coefficients_dict = {'methoxybenzaldehyde': 1, 'ethyl_acetoacetate': 2, 'ammonium_acetate': 1}
        candidate_yields = [product_concentration / (concentrations_df.loc[index, f'c#{substrate_name}'] / coefficients_dict[substrate_name]) for substrate_name in substrates]
        candidate_errs = [product_error / (concentrations_df.loc[index, f'c#{substrate_name}'] / coefficients_dict[substrate_name]) for substrate_name in substrates]
        concentrations_df.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)
        concentrations_df.loc[index, f'yielderr#{product_name}'] = candidate_errs[np.argmax(candidate_yields)]

        # bb017
        product_name = 'bb017'
        if product_name in substances_for_fitting:
            product_concentration = concentrations_df.loc[index, f'pc#{product_name}']
            product_error = concentrations_df.loc[index, f'pcerr#{product_name}']
            coefficients_dict = {'methoxybenzaldehyde': 2, 'ethyl_acetoacetate': 1, 'ammonium_acetate': 2}
            candidate_yields = [product_concentration / (concentrations_df.loc[index, f'c#{substrate_name}'] / coefficients_dict[substrate_name]) for substrate_name in substrates]
            candidate_errs = [product_error / (concentrations_df.loc[index, f'c#{substrate_name}'] / coefficients_dict[substrate_name]) for substrate_name in substrates]
            concentrations_df.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)
            concentrations_df.loc[index, f'yielderr#{product_name}'] = candidate_errs[np.argmax(candidate_yields)]

    target_folder = data_folder + run_name + 'results/historical'
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    target_folder = data_folder + run_name + folder_for_prod_conc
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)
    concentrations_df.to_csv(data_folder + run_name + folder_for_prod_conc + '/product_concentration.csv', index=False)


def plot_timeline(list_of_runs, list_of_calibrant_set_ids, param_to_plot='rmse'):
    # make a dataframe with columns 'param', 'set_index'
    overall_df = pd.DataFrame(columns=['param', 'set_index'], dtype=object)
    # fig, axarr = plt.subplots(1, len(list_of_calibrant_set_ids), figsize=(15, 5), sharey=True, dpi=150)
    for set_index in list_of_calibrant_set_ids:
        df_results = organize_run_results.join_data_from_runs([f'BPRF/{x}/' for x in list_of_runs],
                                                          round_on_columns=None,
                                                          subfolder_containing_csv=f'results/historical/calibrantset_{set_index}')
        param_values = df_results[param_to_plot].values
        # add all the param values to the overall_df with column 'set_index' containing the set_index
        overall_df = pd.concat([overall_df, pd.DataFrame({'param': param_values, 'set_index': [set_index] * len(param_values)}, dtype=object)])
        # make a histogram of the param values

        # # linear scale
        # axarr[list_of_calibrant_set_ids.index(set_index)].hist(param_values, orientation='horizontal')

        # log scale
        # axarr[list_of_calibrant_set_ids.index(set_index)].hist(param_values, orientation='horizontal', bins=np.logspace(np.log10(np.min(param_values)),np.log10(np.max(param_values)), 30))
        # axarr[list_of_calibrant_set_ids.index(set_index)].set_yscale('log')

    fig, ax = plt.subplots(figsize=(5, 4))
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.set_style('whitegrid')
    ax = sns.boxplot(ax=ax, x='set_index', y='param', data=overall_df, color='grey', notch=True, whis=[5, 95])
    ax = sns.stripplot(ax=ax, x="set_index", y="param", data=overall_df, alpha=0.35, color='C0', jitter=0.2, size=3)
    # set y axis to log scale
    ax.set_yscale('log')
    ax.set_ylim(5e-3, 1)
    simpleaxis(ax)
    # disable gridlines
    ax.grid(False)
    # plt.xlabel('Components discovered')
    plt.xlabel('')
    # plt.ylabel('Max. residual, absorbance units')
    # plt.ylabel('')

    # Noise level is displayed as the median (over spectra) of the maximum (over wavelengths) absolute deviation of
    # absorbance between the true and measured spectrum, assuming that the standard error of absorbance at each
    # wavelength is 0.01 absorbance units.
    sigma = 0.01
    N_wavelengths = 349
    max_residual_samples = []
    for M in range(500):
        # sample from normal distribution with sigma = 0.03 and zero mean
        x = np.random.normal(0, sigma, N_wavelengths)
        # calculate the max residual
        max_residual_samples.append(np.max(np.abs(x)))
    max_residual_samples = np.array(max_residual_samples)
    median_maxresidual = np.median(max_residual_samples)
    plt.axhline(y=median_maxresidual, color='C2', linestyle='--', zorder=15, linewidth=3)

    # set x axis labels to the dates of the calibrant sets
    ax.set_xticklabels([f"{sum([len(calibrant_sets[i]['calibrant_set']) for i in range(k+1)])}" for k in list_of_calibrant_set_ids])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_locator(mticker.LogLocator(numticks=999))
    ax.yaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))

    plt.tight_layout()
    fig.savefig('misc-scripts/figures/hntz_discovery_timeline_withstoich1.png', dpi=300)

    plt.show()


if __name__ == '__main__':
    list_of_runs = tuple(['2024-02-16-run01',
                          '2024-02-17-run01',
                          '2024-02-17-run02'])

    ### Uncomment to recalibrate the analysis
    # for set_index in [0, 1, 2, 3, 4, 5, 6, 7]:
    #     print(f'Processing calibrant set {set_index}...')
    #     combined_set = []
    #     for i in range(set_index+1):
    #         combined_set += calibrant_sets[i]['calibrant_set']
    #     print(f'Combined set: {combined_set}')
    #     folder_for_prod_conc = f'results/historical/calibrantset_{set_index}'
    #
    #     for i, run_shortname in enumerate(list_of_runs):
    #         process_run_by_shortname(run_shortname, combined_set, folder_for_prod_conc)

    plot_timeline(list_of_runs, [0, 1, 2, 3, 4, 5, 6, 7], param_to_plot='maxresidual')
