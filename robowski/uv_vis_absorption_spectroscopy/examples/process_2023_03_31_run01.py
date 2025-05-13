from robowski.settings import *
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'

dilution_factor = 200
product_name = 'IIO029A'
run_name = 'multicomp-reactions/2023-03-31-run01/'

sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
exp_names_craic = ['2023-03-31-run01']
true_sequence_of_plates = {(17, '2023-04-04_16-13-11__plate0000017__2023-03-31-run01'): 'skip',
                           ( 3, '2023-04-05_12-17-13__plate0000003__2023-03-31-run01'): 1,
                           ( 6, '2023-04-05_13-10-20__plate0000006__2023-03-31-run01'): 2,
                           ( 9, '2023-04-05_13-24-55__plate0000009__2023-03-31-run01'): 3}

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
    concentrations_here = sp.concentrations_for_one_plate(experiment_folder=data_folder + run_name,
                                                          plate_folder=craic_folder + row['folder'] + '/',
                                                          calibration_folder=data_folder + 'multicomp-reactions/2023-01-18-run01/' + 'microspectrometer_data/calibration/',
                                                          calibrant_shortnames=[product_name, 'ald001'],
                                                          background_model_folder=data_folder + 'multicomp-reactions/2023-03-20-run01/microspectrometer_data/background_model/',
                                                          calibrant_upper_bounds=[np.inf, 1e-10],
                                                          do_plot=False)
    diluted_vials = process_wellplate_spectra.diluted_vials_only(concentrations_here) * dilution_factor
    if true_sequence_of_plates[(row['plate_id'], row['folder'])] == 'skip':
        continue
    id_of_range_of_conditions_for_this_plate = true_sequence_of_plates[(row['plate_id'], row['folder'])]
    concentrations_df.at[id_of_range_of_conditions_for_this_plate * 27:(id_of_range_of_conditions_for_this_plate + 1) * 27 - 1, product_name] = diluted_vials

concentrations_df = concentrations_df[concentrations_df[product_name].notna()]

substrates = ['ald001', 'am001', 'ic001']
concentrations_df['yield'] = concentrations_df[product_name] * 0
for index, row in concentrations_df.iterrows():
    substrate_concentrations_min = min([concentrations_df.at[index, substrate] for substrate in substrates])
    yield_here = concentrations_df.at[index, product_name] / substrate_concentrations_min
    concentrations_df.at[index, 'yield'] = yield_here

concentrations_df.to_csv(data_folder + run_name + 'results/' + 'product_concentration.csv', index=False)
