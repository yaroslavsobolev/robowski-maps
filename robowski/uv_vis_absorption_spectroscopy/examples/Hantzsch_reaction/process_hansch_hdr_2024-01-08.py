from robowski.settings import *
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd


import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra
import robowski.misc_scripts.organize_run_results as organize_run_results
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
plt.ioff()



import logging
logging.basicConfig(level=logging.INFO)

def yield_by_name(product_name, product_concentration, substrate_concentrations_dictionary={'methoxybenzaldehyde': 1, 'ethyl_acetoacetate': 1, 'ammonium_acetate': 1}):
    """
    Calculate the yield of the reaction.
    :param product_name: str
    :param product_concentration: float
    :param substrate_concentrations_dictionary: dict
    :return: float
    """
    possible_yelds = []
    substrate_molar_factors={'methoxybenzaldehyde': 1, 'ethyl_acetoacetate': 2, 'ammonium_acetate': 1}
    if product_name == 'HRP01':
        for substrate_name in substrate_concentrations_dictionary:
            possible_yelds.append(product_concentration / (substrate_concentrations_dictionary[substrate_name] / substrate_molar_factors[substrate_name]))
    else:
        return product_concentration / substrate_concentrations_dictionary[product_name]


def process_run_by_shortname(run_shortname, cut_from=1, dilution_factor=200):
    substrates = ['methoxybenzaldehyde', 'ethyl_acetoacetate', 'ammonium_acetate']
    product_name = 'HRP01'
    run_name = f'BPRF/{run_shortname}/'
    print(f'Processing run {run_name}...')


    # substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'dm36', 'dm37', 'dm40_12', 'dm40_10', 'ethyl_acetoacetate', 'EAB', 'bb017', 'bb021', 'dm70', 'dm053', 'dm088_4', 'bb021_f2']

    # substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'dm36', 'dm37', 'dm40_12', 'dm40_10',
    #                           'ethyl_acetoacetate', 'bb017', 'bb021', 'dm70', 'dm053', 'dm088_4', 'bb021_f2']
    # try removing: 35_8, 35_9, dm37, dm_70

    # substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'ethyl_acetoacetate', 'EAB', 'bb017', 'dm088_4', 'dm053', 'dm70'] # should inlclude dm37 which is oxidized dm70

    # ############### best benchmarks so far. This is what was used before 2024-10-01 for the 2024-04-17-run1-2 data
    # substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'ethyl_acetoacetate', 'EAB', 'bb017', 'dm088_4', 'dm053', 'dm70']

    ############### what was used for the 2024-10-01 reprocessing of the 2024-04-17-run1-2 data
    substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'ethyl_acetoacetate', 'EAB', 'bb017',
                              'dm088_4', 'dm053', 'dm70', 'dm40_12', 'dm40_10']


    sp = process_wellplate_spectra.SpectraProcessor(
        folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
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
                                                                                     list_of_starting_concentration_dicts=list_of_starting_concentration_dicts)

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

        # # dm037
        # product_name = 'dm37'
        # product_concentration = concentrations_df.loc[index, f'pc#{product_name}']
        # coefficients_dict = {'methoxybenzaldehyde': 2, 'ethyl_acetoacetate': 2, 'ammonium_acetate': 1}
        # candidate_yields = [product_concentration / (concentrations_df.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name]) for substrate_name in substrates]
        # concentrations_df.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)

        # bb017
        product_name = 'bb017'
        product_concentration = concentrations_df.loc[index, f'pc#{product_name}']
        product_error = concentrations_df.loc[index, f'pcerr#{product_name}']
        coefficients_dict = {'methoxybenzaldehyde': 2, 'ethyl_acetoacetate': 1, 'ammonium_acetate': 2}
        candidate_yields = [product_concentration / (concentrations_df.loc[index, f'c#{substrate_name}'] / coefficients_dict[substrate_name]) for substrate_name in substrates]
        candidate_errs = [product_error / (concentrations_df.loc[index, f'c#{substrate_name}'] / coefficients_dict[substrate_name]) for substrate_name in substrates]
        concentrations_df.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)
        concentrations_df.loc[index, f'yielderr#{product_name}'] = candidate_errs[np.argmax(candidate_yields)]

        # product_name = 'bb021'
        # product_concentration = concentrations_df.loc[index, f'pc#{product_name}']
        # product_error = concentrations_df.loc[index, f'pcerr#{product_name}']
        # coefficients_dict = {'methoxybenzaldehyde': 2, 'ethyl_acetoacetate': 1, 'ammonium_acetate': 0}
        # candidate_yields = [product_concentration / (concentrations_df.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name]) for substrate_name in substrates
        #                     if coefficients_dict[substrate_name] > 0]
        # candidate_errs = [product_error / (concentrations_df.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name]) for substrate_name in substrates]
        # concentrations_df.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)
        # concentrations_df.loc[index, f'yielderr#{product_name}'] = candidate_errs[np.argmax(candidate_yields)]

        # product_name = 'HRP02'
        # product_concentration = concentrations_df.loc[index, f'pc#dm35_9']
        # product_error = concentrations_df.loc[index, f'pcerr#dm35_9']
        # coefficients_dict = {'methoxybenzaldehyde': 2, 'ethyl_acetoacetate': 2, 'ammonium_acetate': 1}
        # candidate_yields = [product_concentration / (concentrations_df.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name]) for substrate_name in substrates]
        # candidate_errs = [
        #     product_error / (concentrations_df.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name])
        #     for substrate_name in substrates]
        # concentrations_df.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)
        # concentrations_df.loc[index, f'yielderr#{product_name}'] = candidate_errs[np.argmax(candidate_yields)]

        # product_name = 'HRI03'
        # product_concentration = concentrations_df.loc[index, f'pc#EAB']
        # product_error = concentrations_df.loc[index, f'pcerr#EAB']
        # coefficients_dict = {'methoxybenzaldehyde': 0, 'ethyl_acetoacetate': 1, 'ammonium_acetate': 1}
        # candidate_yields = [product_concentration / (concentrations_df.loc[index, f'c#{substrate_name}'] / coefficients_dict[substrate_name]) for substrate_name in substrates if substrate_name != 'methoxybenzaldehyde']
        # candidate_errs = [
        #     product_error / (concentrations_df.loc[index, f'c#{substrate_name}'] / coefficients_dict[substrate_name])
        #     for substrate_name in substrates if substrate_name != 'methoxybenzaldehyde']
        # concentrations_df.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)
        # concentrations_df.loc[index, f'yielderr#{product_name}'] = candidate_errs[np.argmax(candidate_yields)]

        # product_name = 'dm70'
        # product_concentration = concentrations_df.loc[index, f'pc#{product_name}']
        # coefficients_dict = {'methoxybenzaldehyde': 2, 'ethyl_acetoacetate': 2, 'ammonium_acetate': 1}
        # candidate_yields = [product_concentration / (concentrations_df.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name]) for substrate_name in substrates]
        # concentrations_df.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)

        # product_name = 'dm40'
        # product_concentration = concentrations_df.loc[index, f'pc#dm40_10'] + concentrations_df.loc[index, f'pc#dm40_12']
        # product_error = concentrations_df.loc[index, f'pcerr#dm40_10'] + concentrations_df.loc[index, f'pcerr#dm40_12']
        # coefficients_dict = {'methoxybenzaldehyde': 1, 'ethyl_acetoacetate': 1, 'ammonium_acetate': 0}
        # candidate_yields = [
        #                     product_concentration / (concentrations_df.loc[index, f'c#{substrate_name}'] * coefficients_dict[substrate_name])
        #                     for substrate_name in substrates if substrate_name != 'ammonium_acetate']
        # concentrations_df.loc[index, f'yield#{product_name}'] = np.max(candidate_yields)
        # concentrations_df.loc[index, f'yielderr#{product_name}'] = candidate_errs[np.argmax(candidate_yields)]

        # substrate_concentrations = [concentrations_df.loc[index, f'c#{substrate_name}'] for substrate_name in substrates]
        # # multiply the concentration of acetoacetate by 2 because two moles are needed to produce one mole of product
        # substrate_concentrations[0] = substrate_concentrations[0] * 2
        # substrate_concentration_min = min(substrate_concentrations)
        # product_concentration = concentrations_df.loc[index, f'pc#{product_name}']
        # # if the product concentration is negative, replace it with zero
        # product_concentration = max(product_concentration, 0)
        # concentrations_df.loc[index, 'yield'] = product_concentration / substrate_concentration_min

    concentrations_df.to_csv(data_folder + run_name + 'results/' + 'product_concentration.csv', index=False)


def plot_all_spectra_by_shortname(run_shortname):
    run_name = f'BPRF/{run_shortname}/'

    sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                         '2022-12-01/interpolator-dataset/')
    sp.nanodrop_lower_cutoff_of_wavelengths = 220

    df_structure = organize_run_results.load_run_structure(run_name)
    for nanodrop_filepath in df_structure.nanodrop_filepath.unique():
        sp.show_all_spectra(plate_folder=data_folder + nanodrop_filepath)


if __name__ == '__main__':
    # list_of_runs = tuple([
    #                       '2024-01-17-run01'
    #                       # '2023-11-13-run01',
    #                       # '2023-11-14-run01',
    #                       # '2023-11-21-run01'
    # ])

    # list_of_runs = tuple(['2024-01-29-run01',
    #                       '2024-01-29-run02',
    #                       '2024-01-30-run01'
    #                       ])
    #
    # for i, run_shortname in enumerate(list_of_runs):
    #     process_run_by_shortname(run_shortname)
    #
    # list_of_runs = tuple(['2024-02-16-run01',
    #                       '2024-02-17-run01',
    #                       '2024-02-17-run02'])
    #
    # for i, run_shortname in enumerate(list_of_runs):
    #     process_run_by_shortname(run_shortname)
    #
    # list_of_runs = tuple(['2024-01-16-run01',
    #                       '2024-01-16-run02',
    #                       '2024-01-17-run01'])

    # list_of_runs = tuple(['2024-03-04-run01',
    #                       '2024-03-04-run02'])
    # list_of_runs = tuple(['2024-03-06-run01'])
    # list_of_runs = tuple(['2024-03-06-run02'])
    # list_of_runs = tuple(['2024-03-12-run01'])
    # list_of_runs = tuple(['2024-03-20-run01'])

    # list_of_runs = tuple(['2024-04-17-run01',
    #                       '2024-04-17-run02'])
    # list_of_runs = tuple(['2024-04-17-run01'])
    # list_of_runs = tuple(['2024-04-17-run02'])
    list_of_runs = tuple(['2024-09-30-run02'])
    for i, run_shortname in enumerate(list_of_runs):
        process_run_by_shortname(run_shortname)
        # plot_all_spectra_by_shortname(run_shortname)

    # plt.show()