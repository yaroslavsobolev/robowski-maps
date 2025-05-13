from robowski.settings import *

import os
import robowski.uv_vis_absorption_spectroscopy.calibrator as calibrator
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

experiment_name = f'BPRF/2024-01-17-run01/'
cut_from = 0

## Without CARY
# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-21_14-53-09_UV-Vis_main_product',
#     calibrant_shortnames=['HRP01'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[0.015],
#     min_concentrations=[0.00004],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.7557,
#     do_reference_stitching=True,
#     bkg_multiplier=0
# )

# # With CARY
# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-21_14-53-09_UV-Vis_main_product',
#     calibrant_shortnames=['HRP01'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[0.015],
#     min_concentrations=[0.00004],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.7557,
#     do_reference_stitching=False,
#     bkg_multiplier=1,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/Hantzsch-ester-HRP01/HRP01_400ug_per_20mL_repeat1.csv',
#     cary_column_name='HRP01_0.4mg_per_20_mL_repeat1',
#     do_record_residuals=True
# )

# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_18-28-35_UV-Vis_starting_materials',
#     calibrant_shortnames=['methoxybenzaldehyde'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[1],
#     min_concentrations=[4e-5],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=True,
#     do_smoothing_at_low_absorbance=None,
#     # forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/pure_benzald.csv',
#     # cary_column_name='pure_benzald_10ul_in_3mL_cuvette_stock_10uL_in_5mL_c1_rep1_1',
#     # forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/methoxybenzald.csv',
#     # cary_column_name='methoxybenzald__c1_rep1_2',
#     # nanodrop_wavelength_shift = -1
#     do_record_residuals=True,
#     do_not_save_data=True
# )

# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_18-28-35_UV-Vis_starting_materials',
#     calibrant_shortnames=['ethyl_acetoacetate'],
#     ref_concentrations=[0.006],
#     max_concentrations=[0.015],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=0.02,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/ethyl_acetoacetate_dilute.csv',
#     cary_column_name='ethyl_acetoacetate__c1_rep1_2',
#     cut_to=None,
#     do_record_residuals=True,
#     do_not_save_data=True
# )

# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_13-38-39_UV-Vis_Knoevenagel',
#     calibrant_shortnames=['dm40_12'],
#     ref_concentrations=[0.01],
#     max_concentrations=[0.02],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm40_12.csv',
#     cary_column_name='dm40_12__c1_rep1',
#     cut_to=None,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/methoxybenzaldehyde/bkg_spectrum.npy',
#     do_record_residuals=True,
#     do_not_save_data=True,
#     bkg_multiplier=0.3,
#     nanodrop_wavelength_shift = -1
# )

# calibrator.construct_calibrant(
#     cut_from=0,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_13-38-39_UV-Vis_Knoevenagel',
#     calibrant_shortnames=['dm40_10'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[0.011],
#     min_concentrations=[1e-12],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     # forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm40_10.csv',
#     # cary_column_name='dm40_10__c1_rep2',
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
#     bkg_multiplier=0.3,
#     do_record_residuals=True,
#     do_not_save_data=True,
#     # nanodrop_wavelength_shift = +2
# )

# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2023-12-20_17-51-07_UV-Vis_side_prod',
#     calibrant_shortnames=['dm35_8'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[1],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     cut_to=None,
#     bkg_multiplier=1,
#     do_record_residuals=True,
#     do_not_save_data=True,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# calibrator.construct_calibrant(
#     cut_from=0,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2023-12-20_17-51-07_UV-Vis_side_prod',
#     calibrant_shortnames=['dm35_9'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[1],
#     min_concentrations=[4.5e-5],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     # forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm35_9.csv',
#     # cary_column_name='dm35_9_c2_rep2',
#     cut_to=None,
#     do_record_residuals=True,
#     do_not_save_data=True,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2023-12-20_17-51-07_UV-Vis_side_prod',
#     calibrant_shortnames=['dm36'],
#     ref_concentrations=[0.0005],
#     max_concentrations=[0.00075],
#     min_concentrations=[4.5e-5],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     cut_to=None,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
#     do_record_residuals=True,
#     do_not_save_data = True,
# )

# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2023-12-20_17-51-07_UV-Vis_side_prod',
#     calibrant_shortnames=['dm37'],
#     ref_concentrations=[0.0005],
#     max_concentrations=[1],
#     min_concentrations=[2.5e-5],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     cut_to=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/Hantzsch_dm37/dm37.csv',
#     cary_column_name='dm_37_SBW1nm_repeat2',
#     do_smoothing_at_low_absorbance=None,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
#     do_record_residuals=True,
#     do_not_save_data=True,
# )

# calibrator.construct_calibrant(
#     cut_from=80,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_14-33-10_UV-Vis_dm053',
#     calibrant_shortnames=['dm053'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[0.0007],
#     min_concentrations=[0.00004e-12],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm053.csv',
#     cary_column_name='dm053__c1_rep2',
#     cut_to=None,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
#     do_record_residuals=False,
#     do_not_save_data=True,
#     bkg_multiplier=0.7
# )


# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-02-02_14-06-53_UV-Vis_dm070',
#     calibrant_shortnames=['dm70'],
#     ref_concentrations=[0.00128],
#     max_concentrations=[1],
#     min_concentrations=[0.000],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     cut_to=None,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
#     do_record_residuals=True,
#     do_not_save_data=True
# )

# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2023-12-26_15-14-38_UV-Vis_ethylaminobutenoate',
#     calibrant_shortnames=['EAB'],
#     ref_concentrations=[0.0005],
#     max_concentrations=[1],#[0.00085],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/EAB.csv',
#     cary_column_name='EAB__c1_rep1',
#     cut_to=None,
#     bkg_multiplier=1.2,
#     do_record_residuals=True,
#     do_not_save_data=True
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-22_12-37-31_UV-Vis_bb017',
#     calibrant_shortnames=['bb017'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[1],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/bb017_lowconc.csv',
#     cary_column_name='bb017_hemiam__c1_rep2',
#     cut_to=None,
#     bkg_multiplier=1,
#     do_record_residuals = True,
#     do_not_save_data = True
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# calibrator.construct_calibrant(
#     cut_from=0,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-02-02_10-24-02_UV-Vis_bb021',
#     calibrant_shortnames=['bb021'],
#     ref_concentrations=[0.0002],
#     max_concentrations=[0.55],
#     min_concentrations=[0.00],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     # forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/bb021.csv',
#     # cary_column_name='bb021__c1_rep3',
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
#     # nanodrop_wavelength_shift = +3
#     do_record_residuals=True,
#     do_not_save_data=True
# )

# plate_folder = data_folder + 'BPRF/2024-01-17-run01/calibrations/2024-03-14_20-52-15_UV-Vis_dm088_4.csv'
# calibrator.take_median_of_nanodrop_spectra(plate_folder)

# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-03-14_20-52-15_UV-Vis_dm088_4_medianned',
#     calibrant_shortnames=['dm088_4'],
#     ref_concentrations=[0.002],
#     max_concentrations=[0.004],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     # do_reference_stitching=True,
#     cut_to=None,
#     bkg_multiplier=1,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm88.csv',
#     cary_column_name='dm88__c2_rep2',
#     do_record_residuals=True,
#     do_not_save_data=True
# )

# ######### DYE
# experiment_name = f'Yaroslav/'
# calibrator.construct_calibrant(
#     cut_from=0,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='2024-01-20_22-50-00_2024_01_20_UV-Vis_53_serial_dilution',
#     calibrant_shortnames=['crystal_violet'],
#     ref_concentrations=[0.60],
#     max_concentrations=[2],
#     min_concentrations=[0.57],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'blackness-noise-calibration-of-nanodrop/crystal-violet-in-ethanol.csv',
#     cary_column_name='crystal-violet-ethanol-1',
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
#     # nanodrop_wavelength_shift = +3
#     do_record_residuals=True,
#     do_not_save_data=True,
#     no_right_edge_subtraction=True,
#     skip_concentrations=[0.55, 0,56, 0.59, 0.60, 0.65, 0.71, 0.73, 0.77, 0.81, 0.83]
# )

# experiment_name = f'Yaroslav/'
# calibrator.construct_calibrant(
#     cut_from=50,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='2024_04_08_UV-Vis__nilered_congored_1',
#     calibrant_shortnames=['nile_red'],
#     ref_concentrations=[0.9],
#     max_concentrations=[2],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'blackness-noise-calibration-of-nanodrop/nilered/NileRed_1.csv',
#     cary_column_name='NileRed_30uL-in-3030uL_rep1_1',
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
#     # nanodrop_wavelength_shift = +3
#     do_record_residuals=True,
#     do_not_save_data=True,
#     no_right_edge_subtraction=True,
#     dont_save_residuals_below_cut_to=True
# )

# experiment_name = f'Yaroslav/'
# calibrator.construct_calibrant(
#     cut_from=1,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='2024_04_08_UV-Vis__nilered_congored_1',
#     calibrant_shortnames=['congo_red'],
#     ref_concentrations=[0.77],
#     max_concentrations=[2],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'blackness-noise-calibration-of-nanodrop/nilered/CongoRed_1.csv',
#     cary_column_name='cONGORed_80uL-in-3020uL_rep1_1',
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
#     bkg_multiplier=0,
#     # nanodrop_wavelength_shift = +3
#     do_record_residuals=True,
#     do_not_save_data=True,
#     no_right_edge_subtraction=True,
#     dont_save_residuals_below_cut_to=True
# )

### VERSATILITY

# experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Friedel-Crafts_GM/'
# cut_from = 150
# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='2023-10-18_18-46-41_UV-Vis_acetylferrocene_new',
#     calibrant_shortnames=['acetylferrocene'],
#     ref_concentrations=[0.02925],
#     max_concentrations=[1],
#     experiment_name=experiment_name,
#     do_record_residuals=True,
#     do_not_save_data=True,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     cut_to=None,
#     no_right_edge_subtraction=True,
#     custom_bkg_spectrum_npy_file=data_folder + 'nanodrop-spectrophotometer-measurements/versatility_test/Friedel-Crafts_GM/microspectrometer_data/calibration/references/ferrocene/bkg_spectrum.npy',
# )


# experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Diels-Alder2_WaiShing/'
# cut_from = 40
# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='2023-10-25_16-51-53_UV-Vis_diphenylisobenzofuran_crude',
#     calibrant_shortnames=['diphenylisobenzofuran'],
#     ref_concentrations=[0.0002],
#     max_concentrations=[1],
#     experiment_name=experiment_name,
#     do_record_residuals=True,
#     do_not_save_data=True,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     cut_to=None,
#     no_right_edge_subtraction=True,
#     upper_limit_of_absorbance=0.95,
#     dont_save_residuals_below_cut_to=True
# )

# experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Diels-Alder2_WaiShing/'
# Validation tests
# experiment_name = 'BPRF/benchmarks/'
# calibrator.construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='2024-04-16_14-30-16_UV-Vis_2_flush_20_dry_12',
#     calibrant_shortnames=['EAB'],
#     ref_concentrations=[0.03],
#     max_concentrations=[1],#[0.00085],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + f'BPRF/2024-01-17-run01/' + 'calibrations/spectrophotometer_data/other-hantzsch/EAB.csv',
#     cary_column_name='EAB__c1_rep1',
#     cut_to=None,
#     bkg_multiplier=1,
#     do_record_residuals=False,
#     do_not_save_data=True
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Diels-Alder2_WaiShing/'
# Validation tests
experiment_name = 'BPRF/benchmarks/'
calibrator.construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='2024-04-16_14-30-16_UV-Vis_2_flush_20_dry_12',
    calibrant_shortnames=['congo_red'],
    ref_concentrations=[0.05],
    max_concentrations=[1],#[0.00085],
    min_concentrations=[0],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=0.95,
    do_reference_stitching=False,
    do_smoothing_at_low_absorbance=None,
    forced_reference_from_agilent_cary_file=data_folder + f'Yaroslav/' + 'blackness-noise-calibration-of-nanodrop/nilered/CongoRed_1.csv',
    cary_column_name='cONGORed_80uL-in-3020uL_rep1_1',
    cut_to=None,
    bkg_multiplier=1,
    do_record_residuals=False,
    do_not_save_data=True
    # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
)