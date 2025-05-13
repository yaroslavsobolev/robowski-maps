from robowski.settings import *
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra
import robowski.misc_scripts.organize_run_results as organize_run_results

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'

def process_run_by_shortname(run_shortname):
    dilution_factor = 200
    product_name = 'IIO029A'
    run_name = f'multicomp-reactions/{run_shortname}/'

    sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                         '2022-12-01/interpolator-dataset/')

    df_structure = pd.read_csv(data_folder + run_name + 'results/run_structure.csv')

    exp_names_craic = tuple(set(df_structure['craic_folder'].to_list()))
    print(f'{len(exp_names_craic)} craic folders in this run.')

    df = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
    df = df.loc[df['folder'].isin(exp_names_craic)].copy().reset_index()
    concentrations_df = pd.read_csv(data_folder + run_name + 'outVandC/' + 'outC.csv')

    # # make sure that number of rows in concentrations dataframe is number of rows in df times 27 experiments per plate
    # assert len(concentrations_df.index) == len(df.index) * 27

    # add a column for the product concentration, fill it with zeros, then fill with measured values
    concentrations_df[product_name] = np.NaN
    # Fill this column with nans
    # concentrations_df[product_name] = concentrations_df[product_name].fillna(0)

    for index, row in df.iterrows():
        print('Processing CRAIC folder ', row['folder'])
        concentrations_here = sp.concentrations_for_one_plate(experiment_folder=data_folder + run_name,
                                                              plate_folder=craic_folder + row['folder'] + '/',
                                                              calibration_folder=data_folder + 'multicomp-reactions/2023-01-18-run01/' + 'microspectrometer_data/calibration/',
                                                              calibrant_shortnames=[product_name, 'ald001'],
                                                              background_model_folder=data_folder + 'multicomp-reactions/2023-03-20-run01/microspectrometer_data/background_model/',
                                                              calibrant_upper_bounds=[np.inf, 1e-10],
                                                              do_plot=False)
        # this accounts for the dilution
        concentrations_here = concentrations_here * dilution_factor
        for vial_id, product_concentration in enumerate(concentrations_here):
            # indes of this vial in the concentrations_df dataframe
            index_of_this_vial = df_structure.index[(df_structure['vial_id'] == vial_id)
                                                    & (df_structure['craic_folder'] == row['folder'])][0]
            concentrations_df.at[index_of_this_vial, product_name] = product_concentration
        # if true_sequence_of_plates[(row['plate_id'], row['folder'])] == 'skip':
        #     continue
        # id_of_range_of_conditions_for_this_plate = true_sequence_of_plates[(row['plate_id'], row['folder'])]
        # concentrations_df.at[id_of_range_of_conditions_for_this_plate * 27:(id_of_range_of_conditions_for_this_plate + 1) * 27 - 1, product_name] = diluted_vials

    concentrations_df = concentrations_df[concentrations_df[product_name].notna()]

    substrates = ['ald001', 'am001', 'ic001']
    concentrations_df['yield'] = concentrations_df[product_name] * 0
    for index, row in concentrations_df.iterrows():
        substrate_concentrations_min = min([concentrations_df.at[index, substrate] for substrate in substrates])
        yield_here = concentrations_df.at[index, product_name] / substrate_concentrations_min
        concentrations_df.at[index, 'yield'] = yield_here

    for column in ['reaction_plate_id', 'diluted_plate_id', 'vial_id', 'is_outlier']:
        concentrations_df[column] = df_structure[column]

    concentrations_df.to_csv(data_folder + run_name + 'results/' + 'product_concentration.csv', index=False)

if __name__ == '__main__':
    run_shortname = '2023-07-04-run01'
    process_run_by_shortname(run_shortname)

    # list_of_runs = tuple(['2023-06-20-run01',
    #                       '2023-06-21-run01',
    #                       '2023-06-21-run02',
    #                       '2023-06-22-run01',
    #                       '2023-06-22-run02',
    #                       '2023-06-22-run03',
    #                       '2023-06-23-run01',
    #                       '2023-06-23-run02',
    #                       '2023-06-26-run01',
    #                       '2023-06-26-run02',
    #                       '2023-06-27-run01',
    #                       '2023-06-27-run02',
    #                       '2023-06-27-run03',
    #                       '2023-06-28-run01',
    #                       '2023-06-28-run02',
    #                       '2023-06-28-run03'])
    # #
    # # for i, run_name in enumerate(list_of_runs):
    # #     if i < 13:
    # #         print('skipping')
    # #         continue
    # #     print('Processing run ' + run_name)
    # #     process_run_by_shortname(run_name)
    #
    # df_results = organize_run_results.join_data_from_runs([f'multicomp-reactions/{run}/' for run in list_of_runs])
    # substances = ['ic001', 'am001', 'ald001', 'ptsa']
    # substance_titles = ['Isocyanide', 'Amine', 'Aldehyde', 'p-TSA']
    # # substance_titles = ['', '', '', '']
    #
    # padding_rows_count = (df_results[substances] == 0).all(axis=1).sum()
    # print(f"There are {padding_rows_count} padding rows (with zero concentrations of substrates).")
    # df_results = df_results[(df_results[substances] != 0).any(axis=1)]
    #
    # destination_run = 'multicomp-reactions/2023-06-19-run01/'
    # df_results.to_csv(data_folder + destination_run + 'results/' + 'product_concentration.csv', index=False)
