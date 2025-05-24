from robowski.settings import *

import os
import robowski.uv_vis_absorption_spectroscopy.calibrator as calibrator
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
import robowski.uv_vis_absorption_spectroscopy.spectraltools as st

experiment_name = f'BPRF/2024-01-17-run01/'
cut_from = 5

## Without CARY
# calibrator.perform_calibration(
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
# calibrator.perform_calibration(
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

calibrator.perform_calibration(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='calibrations/2024-01-16_18-28-35_UV-Vis_starting_materials',
    calibrant_shortnames=['methoxybenzaldehyde'],
    ref_concentrations=[0.006],
    max_concentrations=[1],
    min_concentrations=[4e-5],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=0.95,
    do_reference_stitching=True,
    do_smoothing_at_low_absorbance=None,
    # forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/pure_benzald.csv',
    # cary_column_name='pure_benzald_10ul_in_3mL_cuvette_stock_10uL_in_5mL_c1_rep1_1',
    forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/methoxybenzald.csv',
    cary_column_name='methoxybenzald__c1_rep1_2',
    nanodrop_wavelength_shift = -1
)

# calibrator.perform_calibration(
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
#     cut_to=None
# )

# calibrator.perform_calibration(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_13-38-39_UV-Vis_Knoevenagel',
#     calibrant_shortnames=['dm40_12'],
#     ref_concentrations=[0.01],
#     max_concentrations=[0.02],
#     min_concentrations=[0.00025],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=True,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm40_12.csv',
#     cary_column_name='dm40_12__c1_rep1',
#     cut_to=None,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/methoxybenzaldehyde/bkg_spectrum.npy',
# )

# calibrator.perform_calibration(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_13-38-39_UV-Vis_Knoevenagel',
#     calibrant_shortnames=['dm40_10'],
#     ref_concentrations=[0.01],
#     max_concentrations=[0.011],
#     min_concentrations=[1e-4],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=True,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm40_10.csv',
#     cary_column_name='dm40_10__c1_rep2',
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
#     bkg_multiplier=0,
#     nanodrop_wavelength_shift = +2
# )

# calibrator.perform_calibration(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2023-12-20_17-51-07_UV-Vis_side_prod',
#     calibrant_shortnames=['dm35_8'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[1],
#     min_concentrations=[4.5e-5],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=True,
#     cut_to=None,
#     bkg_multiplier=1
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# calibrator.perform_calibration(
#     cut_from=80,
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
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm35_9.csv',
#     cary_column_name='dm35_9_c2_rep2',
#     cut_to=None,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# calibrator.perform_calibration(
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
#     do_reference_stitching=True,
#     cut_to=None
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# calibrator.perform_calibration(
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
#     do_smoothing_at_low_absorbance=None
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# # WITH CARY
# calibrator.perform_calibration(
#     cut_from=80,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_14-33-10_UV-Vis_dm053',
#     calibrant_shortnames=['dm053'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[0.00056],
#     min_concentrations=[0.00006],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm053.csv',
#     cary_column_name='dm053__c1_rep2',
#     cut_to=None,
#     bkg_multiplier=0.7
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# # WITHOUT CARY
# calibrator.perform_calibration(
#     cut_from=80,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_14-33-10_UV-Vis_dm053',
#     calibrant_shortnames=['dm053'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[0.0007],
#     min_concentrations=[0.00006],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=True,
#     do_smoothing_at_low_absorbance=None,
#     # forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm053.csv',
#     # cary_column_name='dm053__c1_rep2',
#     cut_to=None,
#     bkg_multiplier=0.7
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )


# calibrator.perform_calibration(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-02-02_14-06-53_UV-Vis_dm070',
#     calibrant_shortnames=['dm70'],
#     ref_concentrations=[0.00128],
#     max_concentrations=[1],
#     min_concentrations=[0.00017],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=True,
#     cut_to=None
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )


############### Stitching EAB spectra from Cary
# data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
# agilent_cary_file = data_folder + f'BPRF/2024-01-17-run01/' + 'calibrations/spectrophotometer_data/other-hantzsch/EAB.csv'
# cary_column_name = 'EAB__c1_rep1'
# data1 = st.read_cary_agilent_csv_spectrum(agilent_cary_file, cary_column_name)
# wavwlengths1, spectrum1 = data1
#
# agilent_cary_file = data_folder + f'BPRF/2024-01-17-run01/' + 'calibrations/spectrophotometer_data/Hantzsch_EAB/EAB_for_stitching.csv'
#
# wavelengths2, spectrum2 = st.read_cary_agilent_csv_spectrum(agilent_cary_file, 'EAB_c2_rep1_1')
# _, spectrum_updated = st.stitch_two_spectra(wavwlengths1, spectrum1, wavelengths2, spectrum2, absorbance_limit=1.0,
#                                          do_plot=True)
#
# wavelengths2, spectrum2 = st.read_cary_agilent_csv_spectrum(agilent_cary_file, 'EAB_c1_rep1_1')
# _, spectrum_updated = st.stitch_two_spectra(wavwlengths1, spectrum_updated, wavelengths2, spectrum2, absorbance_limit=1.0,
#                                          do_plot=True)
#
# wavelengths2, spectrum2 = st.read_cary_agilent_csv_spectrum(agilent_cary_file, 'EAB_c3_rep1_1')
# _, spectrum_updated = st.stitch_two_spectra(wavwlengths1, spectrum_updated, wavelengths2, spectrum2, absorbance_limit=1.0,
#                                          do_plot=True)
#
# wavelengths2, spectrum2 = st.read_cary_agilent_csv_spectrum(agilent_cary_file, 'EAB_c4_rep1_1')
# _, spectrum_updated = st.stitch_two_spectra(wavwlengths1, spectrum_updated, wavelengths2, spectrum2, absorbance_limit=1.0,
#                                          do_plot=True)
#
# st.write_cary_agilent_csv_spectrum(wavwlengths1, spectrum_updated,
#                                 data_folder + f'BPRF/2024-01-17-run01/' + 'calibrations/spectrophotometer_data/Hantzsch_EAB/EAB_stitched.csv',
#                                 'EAB_stitched')


# ##################### Calibrating EAB
# calibrator.perform_calibration(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2023-12-26_15-14-38_UV-Vis_ethylaminobutenoate',
#     calibrant_shortnames=['EAB'],
#     ref_concentrations=[0.0007],
#     max_concentrations=[0.00084],#[0.00085],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/Hantzsch_EAB/EAB_stitched.csv',
#     cary_column_name='EAB_stitched',
#     cut_to=None,
#     bkg_multiplier=1,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# calibrator.perform_calibration(
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
#     bkg_multiplier=1
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# ## with nanodrop
# calibrator.perform_calibration(
#     cut_from=0,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-02-02_10-24-02_UV-Vis_bb021',
#     calibrant_shortnames=['bb021'],
#     ref_concentrations=[0.0002],
#     max_concentrations=[0.00055],
#     min_concentrations=[0.0001],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=True,
#     do_smoothing_at_low_absorbance=None,
#     # forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/bb021.csv',
#     # cary_column_name='bb021__c1_rep3',
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
#     # nanodrop_wavelength_shift = +3
# )

# with CARY and two isomers
# calibrator.perform_calibration(
#     cut_from=0,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-02-02_10-24-02_UV-Vis_bb021',
#     calibrant_shortnames=['bb021'],
#     ref_concentrations=[0.0002],
#     max_concentrations=[0.0004],
#     min_concentrations=[0.0001],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/Hantzsch_other_2/bb021_f1_rep2.csv',
#     cary_column_name='bb021f1_c1_rep2_2',
#     # forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/Hantzsch_other_2/bb021_f2_c2_rep1.csv',
#     # cary_column_name='bb021f2_c2_rep1_2',
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
#     # nanodrop_wavelength_shift = +3
# )

# calibrator.perform_calibration(
#     cut_from=0,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024_04_0_UV-Vis_bb021_F2',
#     calibrant_shortnames=['bb021_f2'],
#     ref_concentrations=[0.0005],
#     max_concentrations=[0.0004],
#     min_concentrations=[0.0001],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/Hantzsch_other_2/bb021_f2_c2_rep1.csv',
#     cary_column_name='bb021f2_c2_rep1_2',
#     cut_to=None
# )


# plate_folder = data_folder + 'BPRF/2024-01-17-run01/calibrations/2024-03-14_20-52-15_UV-Vis_dm088_4.csv'
# calibrator.take_median_of_nanodrop_spectra(plate_folder)

# calibrator.perform_calibration(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-03-14_20-52-15_UV-Vis_dm088_4_medianned',
#     calibrant_shortnames=['dm088_4'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[0.00035],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     cut_to=None,
#     bkg_multiplier=1,
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
#     do_smoothing_at_low_absorbance=None,
#     forced_reference_from_agilent_cary_file=data_folder + experiment_name + 'calibrations/spectrophotometer_data/other-hantzsch/dm88.csv',
#     cary_column_name='dm88__c2_rep2',
# )

# calibrator.perform_calibration(
#     cut_from=0,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/acetic_acid',
#     calibrant_shortnames=['acetic_acid'],
#     ref_concentrations=[0.05978],
#     max_concentrations=[1],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=False,
#     do_smoothing_at_low_absorbance=0.02,
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
#     bkg_multiplier=0
# )