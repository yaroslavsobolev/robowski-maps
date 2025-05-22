from robowski.settings import *
from versatility_examples import process_plate
from robowski.uv_vis_absorption_spectroscopy.calibrator import construct_calibrant
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra

experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Diels-Alder2_WaiShing/'

cut_from = 0

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=False,
    calibration_source_filename='2023-10-25_15-09-00_UV-Vis_dimethyl_acetylenedicarboxylate',
    calibrant_shortnames=['dimethyl_acetylenedicarboxylate'],
    ref_concentrations=[0.0005],
    max_concentrations=[0.1],
    experiment_name=experiment_name,
    no_right_edge_subtraction=True,
    upper_limit_of_absorbance=1e6,
    artefact_generating_upper_limit_of_absorbance=1e6,
    do_smoothing_at_low_absorbance=None
)

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=False,
    calibration_source_filename='2023-10-25_16-51-53_UV-Vis_diphenylisobenzofuran_crude',
    calibrant_shortnames=['diphenylisobenzofuran'],
    ref_concentrations=[0.0002],
    max_concentrations=[1],
    experiment_name=experiment_name,
    no_right_edge_subtraction=True,
    upper_limit_of_absorbance=1e6,
    artefact_generating_upper_limit_of_absorbance=1e6,
    do_smoothing_at_low_absorbance=None
)

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=False,
    calibration_source_filename='2023-10-24_18-46-52_UV-Vis_product',
    calibrant_shortnames=['product'],
    ref_concentrations=[0.0006],
    max_concentrations=[1],
    experiment_name=experiment_name,
    no_right_edge_subtraction=True,
    upper_limit_of_absorbance=1e6,
    artefact_generating_upper_limit_of_absorbance=1e6,
    do_smoothing_at_low_absorbance=None
)


sp = process_wellplate_spectra.SpectraProcessor(
    folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                   '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'
# process_plate(sp, dilution_factor=5000,
#               plate_folder=f'{data_folder}{experiment_name}2023-10-25_16-51-53_UV-Vis_diphenylisobenzofuran_crude.csv',
#               well_ids=[7],
#               cut_from=cut_from,
#               cut_to=300,
#               calibrant_shortnames=['dimethyl_acetylenedicarboxylate', 'diphenylisobenzofuran', 'product'],
#               calibration_folder=calibration_folder,
#               experiment_name=experiment_name,
#               do_plot=True,
#               use_line=False,
#               get_errors_from_fit=True,
#               upper_bounds=[np.inf, np.inf, np.inf])

process_plate(sp, dilution_factor=1000,
              plate_folder=f'{data_folder}{experiment_name}2023-10-30_16-53-59_UV-Vis_crude_2.csv',
              well_ids=[17],
              cut_from=cut_from,
              cut_to=300,
              calibrant_shortnames=['dimethyl_acetylenedicarboxylate', 'diphenylisobenzofuran', 'product'],
              calibration_folder=calibration_folder,
              experiment_name=experiment_name,
              do_plot=True,
              use_line=False,
              get_errors_from_fit=True,
              upper_bounds=[np.inf, np.inf, np.inf])