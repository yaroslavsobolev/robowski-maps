from robowski.settings import *
from versatility_examples import *

experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Glaser_WaiShing/'

cut_from = 70

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='2023-10-16_12-50-44_UV-Vis_phenylacetylene_MeCN',
    calibrant_shortnames=['phenylacetylene'],
    ref_concentrations=[0.0005],
    max_concentrations=[0.1],
    experiment_name=experiment_name,
    no_right_edge_subtraction=True
)


construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='2023-10-16_15-42-27_UV-Vis_diphenylbutadiyne_MeCN',
    calibrant_shortnames=['diphenylbutadiyne'],
    ref_concentrations=[0.0003],
    max_concentrations=[0.00045],
    experiment_name=experiment_name,
    no_right_edge_subtraction=True,
    custom_bkg_spectrum_npy_file=data_folder + 'nanodrop-spectrophotometer-measurements/versatility_test/Glaser_WaiShing/microspectrometer_data/calibration/references/phenylacetylene/bkg_spectrum.npy'
)


construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='2023-10-19_17-26-29_UV-Vis_iodine_MeCN',
    calibrant_shortnames=['iodine'],
    ref_concentrations=[0.000610],
    max_concentrations=[0.610],
    experiment_name=experiment_name,
    no_right_edge_subtraction=False,
    custom_bkg_spectrum_npy_file=data_folder + 'nanodrop-spectrophotometer-measurements/versatility_test/Glaser_WaiShing/microspectrometer_data/calibration/references/phenylacetylene/bkg_spectrum.npy'
)


sp = process_wellplate_spectra.SpectraProcessor(
    folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                   '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'
process_plate(sp, dilution_factor=4000,
              plate_folder=f'{data_folder}{experiment_name}2023-10-16_18-37-40_UV-Vis_crude_MeCN.csv',
              well_ids=[0],
              cut_from=cut_from,
              cut_to=200,
              calibrant_shortnames=['phenylacetylene', 'diphenylbutadiyne', 'iodine'],
              calibration_folder=calibration_folder,
              experiment_name=experiment_name,
              do_plot=True,
              use_line=False,
              get_errors_from_fit=True,
              upper_bounds=[1e-10, np.inf, np.inf])