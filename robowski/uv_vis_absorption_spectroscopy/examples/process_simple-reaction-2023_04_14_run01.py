from robowski.settings import *
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'

product_name = 'prod1'
run_name = 'simple-reactions/2023-04-14-run01/'
temperature = 16
exp_names_craic = ['simple-reactions-2023-04-14-run01-diluted']
true_sequence_of_plates = {( 27, '2023-04-14_15-51-30__plate0000027__simple-reactions-2023-04-14-run01-diluted'): 0,
                           ( 19, '2023-04-14_16-59-58__plate0000019__simple-reactions-2023-04-14-run01-diluted'): 1,
                           ( 34, '2023-04-14_18-31-00__plate0000034__simple-reactions-2023-04-14-run01-diluted'): 2,
                           ( 34, '2023-04-14_18-19-58__plate0000034__simple-reactions-2023-04-14-run01-diluted'): 'skip',
                           ( 34, '2023-04-14_17-30-53__plate0000034__simple-reactions-2023-04-14-run01-diluted'): 'skip'}


sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')

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
exp_names_craic = ['simple-reactions-2023-04-14-run01']
true_sequence_of_plates = {( 21, '2023-04-14_15-06-39__plate0000021__simple-reactions-2023-04-14-run01'): 0,
                           ( 18, '2023-04-14_16-46-13__plate0000018__simple-reactions-2023-04-14-run01'): 2,
                           ( 22, '2023-04-14_16-10-39__plate0000022__simple-reactions-2023-04-14-run01'): 1}

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

product_name = 'unknown-product'
for index, row in df.iterrows():
    concentrations_here = sp.get_absorbance_at_single_wavelength_for_one_plate(
                                                          plate_folder=craic_folder + row['folder'] + '/',
    wavelength_id=149,
    ref_wavelength_id=[199])
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
