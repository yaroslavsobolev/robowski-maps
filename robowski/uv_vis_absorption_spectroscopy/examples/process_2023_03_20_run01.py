from robowski.settings import *
from ..process_wellplate_spectra import *

dilution_factor = 200
experiment_name = 'multicomp-reactions/2023-03-20-run01/'

sp = SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')

df = pd.read_csv(craic_folder + 'database_about_these_folders.csv')
df = df.loc[df['exp_name'] == '-'].copy().reset_index()

concentrations_df = pd.read_csv(data_folder + experiment_name + 'outVandC/' + 'outCRF038202303201421.csv')

# make sure that number of rows in concentrations dataframe is number of rows in df times 27 experiments per plate
assert len(concentrations_df.index) == len(df.index) * 27

# add a column for the product concentration, fill it with zeros, then fill with measured values
concentrations_df['IIO029A'] = concentrations_df['ptsa'] * 0
for index, row in df.iterrows():
    concentrations_here = sp.concentrations_for_one_plate(experiment_folder=data_folder + experiment_name,
                                                          plate_folder=craic_folder + row['folder'] + '/',
                                                          calibration_folder=data_folder + 'multicomp-reactions/2023-01-18-run01/' + 'microspectrometer_data/calibration/',
                                                          calibrant_shortnames=['IIO029A', 'ald001'],
                                                          background_model_folder=data_folder + 'multicomp-reactions/2023-03-20-run01/microspectrometer_data/background_model/',
                                                          calibrant_upper_bounds=[np.inf, 1e-10],
                                                          do_plot=False)
    diluted_vials = diluted_vials_only(concentrations_here) * dilution_factor
    concentrations_df.at[index * 27:(index + 1) * 27 - 1, 'IIO029A'] = diluted_vials

concentrations_df.to_csv(data_folder + experiment_name + 'results/' + 'product_concentration.csv', index=False)

substrates = ['ald001', 'am001', 'ic001']
concentrations_df['yield'] = concentrations_df['IIO029A'] * 0
for index, row in concentrations_df.iterrows():
    substrate_concentrations_min = min([concentrations_df.at[index, substrate] for substrate in substrates])
    yield_here = concentrations_df.at[index, 'IIO029A'] / substrate_concentrations_min
    concentrations_df.at[index, 'yield'] = yield_here

concentrations_df.to_csv(data_folder + experiment_name + 'results/' + 'product_concentration.csv', index=False)
