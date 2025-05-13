from robowski.settings import *
from versatility_examples import *

experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Suzuki_JC/'
calibration_source_filename = '2023-10-11_15-28-13_sub_and_prod_calib_curve_good'
calibrant_shortnames = ['CH3COPhBr', 'CH3CObiPhOMe', 'MeOPhB(OH)2']
ref_concentrations = [0.0002, 0.0002, 0.0004]
max_concentrations = [0.0002, 0.0003, 0.0004]

construct_calibrant(
    cut_from=5,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=False,
    calibration_source_filename=calibration_source_filename,
    calibrant_shortnames=calibrant_shortnames,
    ref_concentrations=ref_concentrations,
    max_concentrations=max_concentrations,
    experiment_name=experiment_name,
)

sp = process_wellplate_spectra.SpectraProcessor(
    folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                   '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'
process_plate(sp, dilution_factor=78.65,
              plate_folder=f'{data_folder}{experiment_name}2023-10-11_16-15-52_crude_after_reaction_78.65_df.csv',
              well_ids=range(7),
              cut_from=5,
              cut_to=False,
              calibrant_shortnames=calibrant_shortnames,
              calibration_folder=calibration_folder,
              experiment_name=experiment_name,
              std_calib=0.0265)
