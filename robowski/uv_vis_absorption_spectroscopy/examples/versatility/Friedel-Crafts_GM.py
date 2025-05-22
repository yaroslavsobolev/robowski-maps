from robowski.settings import *
from versatility_examples import process_plate
from robowski.uv_vis_absorption_spectroscopy.calibrator import construct_calibrant
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra

experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Friedel-Crafts_GM/'

cut_from = 150

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='2023-10-13_15-31-26_UV-Vis_ferrocene',
    calibrant_shortnames=['ferrocene'],
    ref_concentrations=[0.00118],
    max_concentrations=[0.1],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=1e6,
    artefact_generating_upper_limit_of_absorbance=1e6,
    do_smoothing_at_low_absorbance=None
)

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=False,
    calibration_source_filename='2023-10-26_13-31-37_UV-Vis_ferrocene_new',
    calibrant_shortnames=['ferrocene'],
    ref_concentrations=[0.03],
    max_concentrations=[0.1],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=1e6,
    artefact_generating_upper_limit_of_absorbance=1e6,
    do_smoothing_at_low_absorbance=None
)

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=False,
    calibration_source_filename='2023-10-18_18-46-41_UV-Vis_acetylferrocene_new',
    calibrant_shortnames=['acetylferrocene'],
    ref_concentrations=[0.02925],
    max_concentrations=[0.1],
    experiment_name=experiment_name,
    custom_bkg_spectrum_npy_file=data_folder + 'nanodrop-spectrophotometer-measurements/versatility_test/Friedel-Crafts_GM/microspectrometer_data/calibration/references/ferrocene/bkg_spectrum.npy',
    no_right_edge_subtraction=True,
    upper_limit_of_absorbance=1e6,
    artefact_generating_upper_limit_of_absorbance=1e6,
    do_smoothing_at_low_absorbance=None
)

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='2023-11-02_00-07-47_UV-Vis_acetylferrocene_30_20',
    calibrant_shortnames=['acetylferrocene'],
    ref_concentrations=[0.020],
    max_concentrations=[0.1],
    experiment_name=experiment_name,
    custom_bkg_spectrum_npy_file=data_folder + 'nanodrop-spectrophotometer-measurements/versatility_test/Friedel-Crafts_GM/microspectrometer_data/calibration/references/ferrocene/bkg_spectrum.npy',
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
    calibration_source_filename='2023-10-23_15-44-05_UV-Vis_acetylchloride',
    calibrant_shortnames=['acetylchloride'],
    ref_concentrations=[0.01],
    max_concentrations=[0.1],
    experiment_name=experiment_name,
    custom_bkg_spectrum_npy_file=data_folder + 'nanodrop-spectrophotometer-measurements/versatility_test/Friedel-Crafts_GM/microspectrometer_data/calibration/references/ferrocene/bkg_spectrum.npy',
    upper_limit_of_absorbance=1e6,
    artefact_generating_upper_limit_of_absorbance=1e6,
    do_smoothing_at_low_absorbance=None
)


sp = process_wellplate_spectra.SpectraProcessor(
    folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                   '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'
df = process_plate(sp, dilution_factor=1,
              plate_folder=f'{data_folder}{experiment_name}2023-10-18_18-35-39_UV-Vis_crude.csv',
              well_ids=[0],
              cut_from=cut_from,
              cut_to=False,
              calibrant_shortnames=['ferrocene', 'acetylferrocene', 'acetylchloride'],
              calibration_folder=calibration_folder,
              experiment_name=experiment_name,
              do_plot=True,
              use_line=False,
              get_errors_from_fit=True,
              upper_bounds=[np.inf, np.inf, 1e-10])


# df = process_plate(sp, dilution_factor=1,
#               plate_folder=f'{data_folder}{experiment_name}2023-10-27_15-11-19_UV-Vis_centrifuged.csv',
#               well_ids=[44],
#               cut_from=cut_from,
#               cut_to=False,
#               calibrant_shortnames=['ferrocene', 'acetylferrocene', 'acetylchloride'],
#               calibration_folder=calibration_folder,
#               experiment_name=experiment_name,
#               do_plot=True,
#               use_line=False,
#               get_errors_from_fit=True,
#               upper_bounds=[np.inf, np.inf, 1e-10])

# df = process_plate(sp, dilution_factor=1,
#               plate_folder=f'{data_folder}{experiment_name}2023-11-01_11-54-29_UV-Vis_crude_new.csv',
#               well_ids=[51],
#               cut_from=cut_from,
#               cut_to=False,
#               calibrant_shortnames=['ferrocene', 'acetylferrocene', 'acetylchloride'],
#               calibration_folder=calibration_folder,
#               experiment_name=experiment_name,
#               do_plot=True,
#               use_line=False,
#               get_errors_from_fit=True,
#               upper_bounds=[np.inf, np.inf, 1e-10])
