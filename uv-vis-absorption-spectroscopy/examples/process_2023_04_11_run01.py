import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import importlib
process_wellplate_spectra = importlib.import_module("uv-vis-absorption-spectroscopy.process_wellplate_spectra")

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'

product_name = 'prod1'
run_name = 'simple-reactions/2023-04-11-run01/'
temperature = 26

sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
exp_names_craic = ['2023-04-11-run01-diluted']
true_sequence_of_plates = {( 9, '2023-04-11_19-51-37__plate0000009__2023-04-11-run01-diluted'): 0,
                           ( 20, '2023-04-11_20-47-50__plate0000020__2023-04-11-run01-diluted'): 1,
                           ( 18, '2023-04-11_21-15-15__plate0000018__2023-04-11-run01-diluted'): 2}

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
exp_names_craic = ['2023-04-11-run01']
true_sequence_of_plates = {( 17,  '2023-04-11_19-39-38__plate0000017__2023-04-11-run01'): 0,
                           ( 27,  '2023-04-11_20-22-49__plate0000027__2023-04-11-run01'): 1,
                           ( 34, '2023-04-11_21-05-31__plate0000034__2023-04-11-run01'): 2}

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
