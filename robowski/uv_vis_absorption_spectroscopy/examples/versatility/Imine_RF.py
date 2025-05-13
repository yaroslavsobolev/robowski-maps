from versatility_examples import *

experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Imine_RF/'

cut_from = 50

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=False,
    calibration_source_filename='2023-10-10_23-31-55_IM_ref_product',
    calibrant_shortnames=['IM01'],
    ref_concentrations=[0.0005],
    max_concentrations=[0.0011],
    experiment_name=experiment_name,
    custom_bkg_spectrum_npy_file=data_folder + 'nanodrop-spectrophotometer-measurements/versatility_test/Imine_RF/microspectrometer_data/calibration/references/ald01/bkg_spectrum.npy'
)

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=False,
    calibration_source_filename='2023-10-12_21-21-13_IM_ref_substrates',
    calibrant_shortnames=['ald01'],
    ref_concentrations=[0.0002],
    max_concentrations=[0.01],
    experiment_name=experiment_name,
)

sp = process_wellplate_spectra.SpectraProcessor(
    folder_with_correction_dataset='uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                   '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'
process_plate(sp, dilution_factor=200,
              plate_folder=f'{data_folder}{experiment_name}2023-10-13_16-38-07_IM_crude.csv',
              well_ids=[0],
              cut_from=cut_from,
              cut_to=False,
              calibrant_shortnames=['IM01', 'ald01'],
              calibration_folder=calibration_folder,
              experiment_name=experiment_name,
              do_plot=True,
              use_line=False,
              get_errors_from_fit=True)
