from robowski.settings import *
from versatility_examples import process_plate
from robowski.uv_vis_absorption_spectroscopy.calibrator import construct_calibrant
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra

experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Ullman_BP/'

cut_from = 50

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='2023-10-11_20-33-06_UV-Vis_substrates',
    calibrant_shortnames=['carbazole', '4-bromobenzaldehyde'],
    ref_concentrations=[0.0004, 0.0010],
    max_concentrations=[0.1, 0.1],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=1e6,
    artefactogenic_upper_limit_of_absorbance=1e6,
    do_smoothing_at_low_absorbance=None
)


construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=False,
    calibration_source_filename='2023-10-12_19-58-31_UV-Vis_product',
    calibrant_shortnames=['4-(9H-carbazol-9-yl)benzaldehyde'],
    ref_concentrations=[0.0003],
    max_concentrations=[0.1],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=1e6,
    artefactogenic_upper_limit_of_absorbance=1e6,
    do_smoothing_at_low_absorbance=None
)

sp = process_wellplate_spectra.SpectraProcessor(
    folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                   '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'
df = process_plate(sp, dilution_factor=1914,
              plate_folder=f'{data_folder}{experiment_name}2023-10-14_15-53-35_UV-Vis_crude.csv',
              well_ids=range(6),
              cut_from=cut_from,
              cut_to=200,
              calibrant_shortnames=['carbazole', '4-bromobenzaldehyde', '4-(9H-carbazol-9-yl)benzaldehyde'],
              calibration_folder=calibration_folder,
              experiment_name=experiment_name,
              do_plot=True,
              use_line=False,
              get_errors_from_fit=True)
