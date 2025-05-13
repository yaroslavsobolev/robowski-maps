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

def process_run_by_shortname(run_shortname):
    dilution_factor = 100
    substrate_name = 'SN1OH03'
    product_name = 'SN1Br03'
    byproduct_name = 'SN1OH03'
    run_name = f'simple-reactions/{run_shortname}/'

    # Find the data about temperature in the readme.md file of the run. Its first line looks like
    # "temperature: 26C", where 26 is the temperature in degrees Celcius
    with open(data_folder + run_name + 'readme.md', 'r') as f:
        first_line = f.readline()
        temperature = first_line.split(' ')[-1].strip('\n').strip('C')

    sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                         '2022-12-01/interpolator-dataset/')

    # df_structure = pd.read_csv(data_folder + run_name + 'results/run_structure.csv')
    df_structure = organize_run_results.load_run_structure(run_name)

    # concentrations_df = pd.read_csv(data_folder + run_name + 'outVandC/' + 'outC.csv')

    # # make sure that number of rows in concentrations dataframe is number of rows in df times 27 experiments per plate
    # assert len(concentrations_df.index) == len(df.index) * 27

    # add a column for the product concentration, fill it with zeros, then fill with measured values
    # concentrations_df[product_name] = np.NaN
    # Fill this column with nans
    # concentrations_df[product_name] = concentrations_df[product_name].fillna(0)

    for nanodrop_filepath in df_structure.nanodrop_filepath.unique():
        print('Processing nanodrop filepath', nanodrop_filepath)
        concentrations_here = sp.concentrations_for_one_plate(experiment_folder=data_folder + run_name,
                                                              plate_folder=data_folder + nanodrop_filepath,
                                                              calibration_folder=data_folder + 'simple-reactions/2023-08-21-run01/' + 'microspectrometer_data/calibration/',
                                                              calibrant_shortnames=[product_name, byproduct_name, 'HBr'],
                                                              background_model_folder=data_folder + 'simple-reactions/2023-08-21-run01/microspectrometer_data/background_model/',
                                                              calibrant_upper_bounds=[np.inf, np.inf, 0.01],
                                                              do_plot=False, return_all_substances=True, cut_from=42, cut_to=180,
                                                              ignore_abs_threshold=True, ignore_pca_bkg=True)
        concentrations_here = concentrations_here * dilution_factor
        for vial_id, product_concentrations in enumerate(concentrations_here):
            # index of this vial in the concentrations_df dataframe
            try:
                index_of_this_vial = df_structure.index[(df_structure['container_id'] == vial_id)
                                                        & (df_structure['nanodrop_filepath'] == nanodrop_filepath)][0]
            except IndexError:
                logging.warning(f'Container_id {vial_id} not present in conditions dataframe but was measured in nanodrop file {nanodrop_filepath}')
                continue


            df_structure.loc[index_of_this_vial, f'pc#{product_name}'] = product_concentrations[0]
            df_structure.loc[index_of_this_vial, f'pc#{byproduct_name}'] = product_concentrations[1]
            df_structure.loc[index_of_this_vial, f'pc#HBr'] = product_concentrations[2]

    concentrations_df = df_structure[df_structure[f'pc#{product_name}'].notna()]

    # Calculate yields
    for index, row in concentrations_df.iterrows():
        substrate_concentration_before_reaction = concentrations_df.loc[index, f'c#{substrate_name}']
        product_concentration = concentrations_df.loc[index, f'pc#{product_name}']
        product_concentration = max(product_concentration, 0)
        substrate_concentration_after_reaction = concentrations_df.loc[index, f'pc#{byproduct_name}']

        if substrate_concentration_before_reaction == 0:
            concentrations_df.loc[index, 'yield'] = 0
            concentrations_df.loc[index, 'conversion'] = 0
        else:
            concentrations_df.loc[index, 'yield'] = product_concentration / (product_concentration + substrate_concentration_after_reaction)
            conversion = (substrate_concentration_before_reaction - substrate_concentration_after_reaction) / \
                         substrate_concentration_before_reaction
            # concentrations_df.loc[index, 'yield'] = np.clip(conversion, 0, 1)
            # if conversion is below 0 or above 1, replace with 0 or 1, respectively, and write to the dataframe
            concentrations_df.loc[index, 'conversion'] = np.clip(conversion, 0, 1)
            if concentrations_df.loc[index, f'c#HBr'] == 0:
                concentrations_df.loc[index, 'HBr_relative_change'] = 1
            else:
                concentrations_df.loc[index, 'HBr_relative_change'] = concentrations_df.loc[index, f'pc#HBr'] / concentrations_df.loc[index, f'c#HBr']

    concentrations_df.to_csv(data_folder + run_name + 'results/' + 'product_concentration.csv', index=False)

if __name__ == '__main__':
    # run_shortname = '2023-08-21-run01'
    # process_run_by_shortname(run_shortname)

    list_of_runs = tuple([
                          '2023-12-11-run01',
                          '2023-12-11-run02',
                          '2023-12-12-run01',
                          '2023-12-12-run02',
                          '2023-12-16-run01',
                          '2023-12-16-run02'])
    for i, run_shortname in enumerate(list_of_runs):
        process_run_by_shortname(run_shortname)