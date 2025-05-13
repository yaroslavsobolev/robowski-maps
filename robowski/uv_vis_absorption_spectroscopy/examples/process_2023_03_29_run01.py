import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import importlib
process_wellplate_spectra = importlib.import_module("uv_vis_absorption_spectroscopy.process_wellplate_spectra")

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'

dilution_factor = 200
product_name = 'IIO029A'
run_name = 'multicomp-reactions/2023-03-29-run01/'

sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
exp_names_craic = ['multicomponent_0329_zoom_in', '2023-03-29-run01']
true_sequence_of_plates = {(17, '2023-03-30_15-59-11__plate0000017__multicomponent_0329_zoom_in'): 1,
                           (18, '2023-03-30_16-24-16__plate0000018__multicomponent_0329_zoom_in'): 2,
                           (11, '2023-03-30_17-01-03__plate0000011__multicomponent_0329_zoom_in'): 0,
                           (11, '2023-04-04_14-07-38__plate0000011__2023-03-29-run01') :           'skip',
                           (16, '2023-04-04_14-36-47__plate0000016__2023-03-29-run01') :           'skip',
                           (21, '2023-04-05_15-02-19__plate0000021__2023-03-29-run01') :           3,
                           (22, '2023-04-05_15-26-12__plate0000022__2023-03-29-run01') :           4,
                           (27, '2023-04-05_15-57-05__plate0000027__2023-03-29-run01') :           5}

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
    if true_sequence_of_plates[(row['plate_id'], row['folder'])] == 'skip':
        continue
    concentrations_here = sp.concentrations_for_one_plate(experiment_folder=data_folder + run_name,
                                                          plate_folder=craic_folder + row['folder'] + '/',
                                                          calibration_folder=data_folder + 'multicomp-reactions/2023-01-18-run01/' + 'microspectrometer_data/calibration/',
                                                          calibrant_shortnames=[product_name, 'ald001'],
                                                          background_model_folder=data_folder + 'multicomp-reactions/2023-03-20-run01/microspectrometer_data/background_model/',
                                                          calibrant_upper_bounds=[np.inf, 1e-10],
                                                          do_plot=False)
    diluted_vials = process_wellplate_spectra.diluted_vials_only(concentrations_here) * dilution_factor
    id_of_range_of_conditions_for_this_plate = true_sequence_of_plates[(row['plate_id'], row['folder'])]
    last_index = (id_of_range_of_conditions_for_this_plate + 1) * 27 - 1
    if last_index > len(concentrations_df.index):
        number_of_vials_outside_of_range = (last_index - len(concentrations_df.index) + 1)
        last_index = len(concentrations_df.index)
    else:
        number_of_vials_outside_of_range = 0
    concentrations_df.at[id_of_range_of_conditions_for_this_plate * 27:last_index, product_name] = diluted_vials[:diluted_vials.shape[0] - number_of_vials_outside_of_range]

concentrations_df = concentrations_df[concentrations_df[product_name].notna()]

substrates = ['ald001', 'am001', 'ic001']
concentrations_df['yield'] = concentrations_df[product_name] * 0
for index, row in concentrations_df.iterrows():
    substrate_concentrations_min = min([concentrations_df.at[index, substrate] for substrate in substrates])
    yield_here = concentrations_df.at[index, product_name] / substrate_concentrations_min
    concentrations_df.at[index, 'yield'] = yield_here

concentrations_df.to_csv(data_folder + run_name + 'results/' + 'product_concentration.csv', index=False)
