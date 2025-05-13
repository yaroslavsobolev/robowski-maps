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
    product_name = 'SN1OH01'
    byproduct_name = 'SN1Br01s1'
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

    df_structure['Temperature'] = temperature

    exp_names_craic = tuple(set(df_structure['craic_folder'].to_list()))
    print(f'{len(exp_names_craic)} craic folders in this run.')

    df = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
    df = df.loc[df['folder'].isin(exp_names_craic)].copy().reset_index()
    # concentrations_df = pd.read_csv(data_folder + run_name + 'outVandC/' + 'outC.csv')

    # # make sure that number of rows in concentrations dataframe is number of rows in df times 27 experiments per plate
    # assert len(concentrations_df.index) == len(df.index) * 27

    # add a column for the product concentration, fill it with zeros, then fill with measured values
    # concentrations_df[product_name] = np.NaN
    # Fill this column with nans
    # concentrations_df[product_name] = concentrations_df[product_name].fillna(0)

    for index, row in df.iterrows():
        print('Processing CRAIC folder ', row['folder'])
        concentrations_here = sp.concentrations_for_one_plate(experiment_folder=data_folder + run_name,
                                                              plate_folder=craic_folder + row['folder'] + '/',
                                                              calibration_folder=data_folder + 'simple-reactions/2023-07-05-run01/' + 'microspectrometer_data/calibration/',
                                                              calibrant_shortnames=[product_name, byproduct_name],
                                                              background_model_folder=data_folder + 'simple-reactions/2023-07-05-run01/microspectrometer_data/background_model/',
                                                              calibrant_upper_bounds=[np.inf, np.inf],
                                                              do_plot=False, return_all_substances=True, cut_from=79)
        concentrations_here = concentrations_here * dilution_factor

        for vial_id, product_concentrations in enumerate(concentrations_here):
            # index of this vial in the concentrations_df dataframe
            index_of_this_vial = df_structure.index[(df_structure['vial_id'] == vial_id)
                                                    & (df_structure['craic_folder'] == row['folder'])][0]
            df_structure.loc[index_of_this_vial, f'pc#{product_name}'] = product_concentrations[0]
            df_structure.loc[index_of_this_vial, f'pc#{byproduct_name}'] = product_concentrations[1]


    ##################### CARBOCATION ABUNDANCE CALCULATION #####################
    exp_names_craic = tuple(set(df_structure['craic_folder_undil'].to_list()))
    print(f'{len(exp_names_craic)} craic folders in this run.')

    df = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
    df = df.loc[df['folder'].isin(exp_names_craic)].copy().reset_index()

    for index, row in df.iterrows():
        print('Processing CRAIC folder ', row['folder'])
        concentrations_here = sp.get_absorbance_at_single_wavelength_for_one_plate(
            plate_folder=craic_folder + row['folder'] + '/',
            wavelength_id=363,
            ref_wavelength_id=[756])

        for vial_id, product_concentrations in enumerate(concentrations_here):
            if product_concentrations > 0.06:
                wavelengths = sp.load_msp_by_id(plate_folder=craic_folder + row['folder'] + '/', well_id=0)[:, 0]
                spectrum = sp.load_msp_by_id(plate_folder=craic_folder + row['folder'] + '/', well_id=vial_id)[:, 1]
                plt.plot(wavelengths, spectrum)
                print(f'product_concentrations = {product_concentrations} for vial_id = {vial_id} in folder {row["folder"]}')
                plt.show()
                # load the spectrum from the plate and show it

            # index of this vial in the concentrations_df dataframe
            index_of_this_vial = df_structure.index[(df_structure['vial_id'] == vial_id)
                                                    & (df_structure['craic_folder_undil'] == row['folder'])][0]
            df_structure.loc[index_of_this_vial, f'pc#carbocat'] = product_concentrations

    concentrations_df = df_structure[df_structure[f'pc#{product_name}'].notna()]

    concentrations_df['yield'] = concentrations_df[f'c#{product_name}'] * 0
    for index, row in concentrations_df.iterrows():
        substrate_concentration_before_reaction = concentrations_df.at[index, 'c#SN1OH01']
        substrate_concentration_after_reaction = concentrations_df.at[index, f'pc#{product_name}']
        byproduct_concentration = concentrations_df.at[index, f'pc#{byproduct_name}']

        # If any of these concentrations are negative, use zero instead
        substrate_concentration_after_reaction = max(0, substrate_concentration_after_reaction)
        byproduct_concentration = max(0, byproduct_concentration)

        sum_concentration_of_all_products = substrate_concentration_before_reaction - substrate_concentration_after_reaction
        main_product_concentration = sum_concentration_of_all_products - byproduct_concentration
        if substrate_concentration_before_reaction == 0:
            concentrations_df.loc[index, 'yield'] = 0
            concentrations_df.loc[index, 'yield#SN1Br01'] = 0
            concentrations_df.loc[index, 'yield#SN1Br01s1'] = 0
        else:
            concentrations_df.loc[index, 'yield'] = sum_concentration_of_all_products / substrate_concentration_before_reaction
            concentrations_df.loc[index, 'yield#SN1Br01'] = main_product_concentration / substrate_concentration_before_reaction
            concentrations_df.loc[index, 'yield#SN1Br01s1'] = byproduct_concentration / substrate_concentration_before_reaction

    for column in ['Temperature', 'reaction_plate_id', 'diluted_plate_id', 'vial_id', 'is_outlier']:
        concentrations_df[column] = df_structure[column]

    concentrations_df.to_csv(data_folder + run_name + 'results/' + 'product_concentration.csv', index=False)

if __name__ == '__main__':
    # run_shortname = '2023-07-05-run01'
    # process_run_by_shortname(run_shortname)

    list_of_runs = tuple([
                          '2023-07-05-run01',
                          '2023-07-06-run01',
                          '2023-07-07-run01',
                          '2023-07-10-run01',
                          '2023-07-10-run02',
                          '2023-07-11-run01',
                          '2023-07-11-run02'])
    for i, run_shortname in enumerate(list_of_runs):
        process_run_by_shortname(run_shortname)