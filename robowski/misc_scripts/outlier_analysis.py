import importlib
import os

organize_run_results = importlib.import_module("misc_scripts.organize_run_results")
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

list_of_runs = tuple(['2024-01-29-run01',
                      '2024-01-29-run02',
                      '2024-01-30-run01'
                      ])

substances = ['c#ethyl_acetoacetate',  'c#methoxybenzaldehyde', 'c#ammonium_acetate']
substance_titles = ['Acetoacetate', 'Methoxy', 'Ammonium acetate']
# substrates = ['c#SN1OH03', 'c#HBr']

df_results = organize_run_results.join_data_from_runs([f'BPRF/{x}/' for x in list_of_runs],
                                 round_on_columns=[])
target_folder_for_images = data_folder + 'BPRF/2024-01-29-run01/results/outlier_spectra'

df_outliers = df_results[(df_results['fitted_dilution_factor_2'] < 190) | (df_results['fitted_dilution_factor_2'] > 210)]
for index, row in df_outliers.iterrows():
    well_id = row['container_id']
    plate_name = row['nanodrop_filepath']
    experiment_folder = data_folder + row['run_name']
    spectrum_picture_filepath = experiment_folder + f'results/uv-vis-fits/{plate_name}-well{well_id:02d}.png.png'
    # copy this file into the target_folder_for_images
    import shutil
    shutil.copy(spectrum_picture_filepath, target_folder_for_images + f'/dilfacfit{row["fitted_dilution_factor_2"]:.0f}_{plate_name}-well{well_id:02d}.png')

df_outliers.to_csv(data_folder + 'BPRF/2024-01-29-run01/results/outlier_dataframe.csv')