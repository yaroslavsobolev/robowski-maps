from versatility_examples import *

experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Click_DM/'

cut_from = 3
# cut_from = 52

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=False,
    calibration_source_filename='2023-10-18_20-31-42__substrates_product_crude',
    calibrant_shortnames=['1-ethynyl-4-nitrobenzene', 'BnN3', 'triazole'],
    ref_concentrations=[0.0005, 0.0005, 0.0005],
    max_concentrations=[0.1, 0.1, 0.1],
    experiment_name=experiment_name
)

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='2023-10-27_16-48-05_UV-Vis_crude4',
    calibrant_shortnames=['CuI'],
    ref_concentrations=[0.00025],
    max_concentrations=[0.1],
    experiment_name=experiment_name
)


sp = process_wellplate_spectra.SpectraProcessor(
    folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                   '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'

# df = process_plate(sp, dilution_factor=2000,
#               plate_folder=f'{data_folder}{experiment_name}2023-10-23_14-05-35_crude1.csv',
#               well_ids=[1],
#               cut_from=cut_from,
#               cut_to=240,
#               calibrant_shortnames=['1-ethynyl-4-nitrobenzene', 'BnN3', 'triazole'],
#               calibration_folder=calibration_folder,
#               experiment_name=experiment_name,
#               do_plot=True,
#               use_line=False,
#               get_errors_from_fit=True,
#               upper_bounds=[np.inf, 1e-10, np.inf])

df = process_plate(sp, dilution_factor=500,
              plate_folder=f'{data_folder}{experiment_name}2023-10-26_15-38-47_UV-Vis_crude2.csv',
              well_ids=[2],
              cut_from=cut_from,
              cut_to=300,
              calibrant_shortnames=['1-ethynyl-4-nitrobenzene', 'triazole', 'CuI'],
              calibration_folder=calibration_folder,
              experiment_name=experiment_name,
              do_plot=True,
              use_line=False,
              get_errors_from_fit=True,
              upper_bounds=[np.inf, np.inf, np.inf])


# df = process_plate(sp, dilution_factor=1000,
#               plate_folder=f'{data_folder}{experiment_name}2023-10-26_15-38-47_UV-Vis_crude2.csv',
#               well_ids=[3],
#               cut_from=cut_from,
#               cut_to=240,
#               calibrant_shortnames=['1-ethynyl-4-nitrobenzene', 'BnN3', 'triazole'],
#               calibration_folder=calibration_folder,
#               experiment_name=experiment_name,
#               do_plot=True,
#               use_line=False,
#               get_errors_from_fit=True,
#               upper_bounds=[np.inf, 1e-10, np.inf])



# df = process_plate(sp, dilution_factor=500,
#               plate_folder=f'{data_folder}{experiment_name}2023-10-27_13-14-40_UV-Vis_crude3.csv',
#               well_ids=[1],
#               cut_from=cut_from,
#               cut_to=240,
#               calibrant_shortnames=['1-ethynyl-4-nitrobenzene', 'BnN3', 'triazole'],
#               calibration_folder=calibration_folder,
#               experiment_name=experiment_name,
#               do_plot=True,
#               use_line=False,
#               get_errors_from_fit=True,
#               upper_bounds=[np.inf, 1e-10, np.inf])


# df = process_plate(sp, dilution_factor=1000,
#               plate_folder=f'{data_folder}{experiment_name}2023-10-27_13-14-40_UV-Vis_crude3.csv',
#               well_ids=[2],
#               cut_from=cut_from,
#               cut_to=240,
#               calibrant_shortnames=['1-ethynyl-4-nitrobenzene', 'BnN3', 'triazole'],
#               calibration_folder=calibration_folder,
#               experiment_name=experiment_name,
#               do_plot=True,
#               use_line=False,
#               get_errors_from_fit=True,
#               upper_bounds=[np.inf, 1e-10, np.inf])

# df = process_plate(sp, dilution_factor=1000,
#               plate_folder=f'{data_folder}{experiment_name}2023-10-27_15-01-43_UV-Vis_mixture.csv',
#               well_ids=[1],
#               cut_from=cut_from,
#               cut_to=240,
#               calibrant_shortnames=['1-ethynyl-4-nitrobenzene', 'BnN3', 'triazole'],
#               calibration_folder=calibration_folder,
#               experiment_name=experiment_name,
#               do_plot=True,
#               use_line=False,
#               get_errors_from_fit=True,
#               upper_bounds=[np.inf, 1e-10, np.inf])