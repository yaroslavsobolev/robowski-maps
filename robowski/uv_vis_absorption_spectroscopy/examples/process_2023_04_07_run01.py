from robowski.settings import *
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'

temps_by_run_name: dict = { 'simple-reactions/2023-04-07-run01/': 26,
                            'simple-reactions/2023-04-11-run01/': 26,
                            'simple-reactions/2023-04-12-run01/': 36,
                            'simple-reactions/2023-04-14-run01/': 16,
                            'simple-reactions/2023-04-14-run02/': 21,
                            'simple-reactions/2023-04-15-run01/': 31,
                            'simple-reactions/2023-04-15-run02/': 11}
product_name = 'prod1'
run_name = 'simple-reactions/2023-04-07-run01/'
temperature = temps_by_run_name[run_name]

sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
exp_names_craic = ['2023-04-07-run01-diluted']
true_sequence_of_plates = {( 21, '2023-04-08_16-06-36__plate0000021__2023-04-07-run01-diluted'): 0,
                           ( 16, '2023-04-08_16-54-40__plate0000016__2023-04-07-run01-diluted'): 1,
                           ( 19, '2023-04-08_18-17-57__plate0000019__2023-04-07-run01-diluted'): 2}

df = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
df = df.loc[df['exp_name'].isin(exp_names_craic)].copy().reset_index()
concentrations_df = pd.read_csv(data_folder + run_name + 'outVandC/' + 'outC.csv')

# # make sure that number of rows in concentrations dataframe is number of rows in df times 27 experiments per plate
# assert len(concentrations_df.index) == len(df.index) * 27

# add a column for the product concentration, fill it with zeros, then fill with measured values
concentrations_df[product_name] = np.NaN
# Fill this column with nans
# concentrations_df[product_name] = concentrations_df[product_name].fillna(0)

for index, row in df.iterrows():
    concentrations_here = sp.get_absorbance_at_single_wavelength_for_one_plate(
                                                          plate_folder=craic_folder + row['folder'] + '/',
    wavelength_id=98,
    ref_wavelength_id=[198])
    # carbocat_concentrations = sp.get_absorbance_at_single_wavelength_for_one_plate(plate_folder=craic_folder + row['folder'] + '/',
    # wavelength_id=363,
    # ref_wavelength_id=756)
    if true_sequence_of_plates[(row['plate_id'], row['folder'])] == 'skip':
        continue
    id_of_range_of_conditions_for_this_plate = true_sequence_of_plates[(row['plate_id'], row['folder'])]
    last_index = (id_of_range_of_conditions_for_this_plate + 1) * 54 - 1
    if last_index > len(concentrations_df.index):
        number_of_vials_outside_of_range = (last_index - len(concentrations_df.index) + 1)
        last_index = len(concentrations_df.index)
    else:
        number_of_vials_outside_of_range = 0
    concentrations_df.at[id_of_range_of_conditions_for_this_plate * 54:last_index, product_name] = concentrations_here[
                                                                                                   :concentrations_here.shape[
                                                                                                        0] - number_of_vials_outside_of_range]
    # concentrations_df.at[id_of_range_of_conditions_for_this_plate * 54:last_index, 'carbocat'] = carbocat_concentrations[
    #                                                                                                :concentrations_here.shape[
    #                                                                                                     0] - number_of_vials_outside_of_range]


concentrations_df = concentrations_df[concentrations_df[product_name].notna()]

substrates = ['SN1OH01', 'HBr']
concentrations_df['yield'] = concentrations_df[product_name] * 0
for index, row in concentrations_df.iterrows():
    substrate_concentrations = [concentrations_df.at[index, substrate] for substrate in substrates]
    substrate_concentrations_min = min(substrate_concentrations)
    print(f'Substrate that has minimal concentration is: {substrates[substrate_concentrations.index(substrate_concentrations_min)]}')

    yield_here = concentrations_df.at[index, product_name] / substrate_concentrations_min
    concentrations_df.at[index, 'yield'] = yield_here

concentrations_df.to_csv(data_folder + run_name + 'results/' + 'product_concentration.csv', index=False)

# Carbocation concentration

product_name = 'carbocat'
exp_names_craic = ['2023-04-07-run01']
true_sequence_of_plates = {( 6,  '2023-04-08_11-09-22__plate0000006__2023-04-07-run01'): 0,
                           ( 9,  '2023-04-08_11-43-19__plate0000009__2023-04-07-run01'): 1,
                           ( 11, '2023-04-08_12-56-26__plate0000011__2023-04-07-run01'): 2}

df = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
df = df.loc[df['exp_name'].isin(exp_names_craic)].copy().reset_index()

# add a column for the product concentration, fill it with zeros, then fill with measured values
concentrations_df[product_name] = np.NaN

for index, row in df.iterrows():
    concentrations_here = sp.get_absorbance_at_single_wavelength_for_one_plate(
                                                          plate_folder=craic_folder + row['folder'] + '/',
    wavelength_id=363,
    ref_wavelength_id=[756])
    if true_sequence_of_plates[(row['plate_id'], row['folder'])] == 'skip':
        continue
    id_of_range_of_conditions_for_this_plate = true_sequence_of_plates[(row['plate_id'], row['folder'])]
    last_index = (id_of_range_of_conditions_for_this_plate + 1) * 54 - 1
    if last_index > len(concentrations_df.index):
        number_of_vials_outside_of_range = (last_index - len(concentrations_df.index) + 1)
        last_index = len(concentrations_df.index)
    else:
        number_of_vials_outside_of_range = 0
    concentrations_df.at[id_of_range_of_conditions_for_this_plate * 54:last_index, product_name] = concentrations_here[
                                                                                                   :concentrations_here.shape[
                                                                                                        0] - number_of_vials_outside_of_range]
concentrations_df['Temperature'] = temperature
concentrations_df.to_csv(data_folder + run_name + 'results/' + 'product_concentration.csv', index=False)
