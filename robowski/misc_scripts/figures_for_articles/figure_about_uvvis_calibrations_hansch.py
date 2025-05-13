import matplotlib.pyplot as plt

from calibration_illustration import *

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

experiment_name = f'BPRF/2024-01-17-run01/'

f = plt.figure(figsize=(4.9, 4.5 * 2/3 * 1.1), dpi=300)
cut_from=5

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='calibrations/2024-01-21_14-53-09_UV-Vis_main_product',
    calibrant_shortnames=['HRP01'],
    ref_concentrations=[0.0003],
    max_concentrations=[0.001],
    min_concentrations=[0.0001],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=0.7557,
    do_reference_stitching=True,
    bkg_multiplier=0,
    calibrant_colors=['C0']
)

# construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_18-28-35_UV-Vis_starting_materials',
#     calibrant_shortnames=['methoxybenzaldehyde'],
#     ref_concentrations=[0.0005],
#     max_concentrations=[1],
#     min_concentrations=[4e-5],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=True
# )

# construct_calibrant(
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
#     do_reference_stitching=True,
#     cut_to=None
# )

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='calibrations/2024-01-16_13-38-39_UV-Vis_Knoevenagel',
    calibrant_shortnames=['dm40_12'],
    ref_concentrations=[0.0003],
    max_concentrations=[0.0004],
    min_concentrations=[0.00002],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=0.95,
    do_reference_stitching=True,
    cut_to=None,
    calibrant_colors=['C3'],
    # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
    custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/methoxybenzaldehyde/bkg_spectrum.npy',
)

# construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-16_13-38-39_UV-Vis_Knoevenagel',
#     calibrant_shortnames=['dm40_10'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[1],
#     min_concentrations=[1e-4],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=True,
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
#     bkg_multiplier=0
# )

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='calibrations/2023-12-20_17-51-07_UV-Vis_side_prod',
    calibrant_shortnames=['dm35_8'],
    ref_concentrations=[0.0003],
    max_concentrations=[0.0003],
    min_concentrations=[4.5e-5],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=0.95,
    do_reference_stitching=True,
    cut_to=None,
    bkg_multiplier=1,
    calibrant_colors=['C1']
    # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
)

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='calibrations/2023-12-20_17-51-07_UV-Vis_side_prod',
    calibrant_shortnames=['dm35_9'],
    ref_concentrations=[0.0003],
    max_concentrations=[0.0003],
    min_concentrations=[4.5e-5],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=0.95,
    do_reference_stitching=True,
    cut_to=None,
    calibrant_colors=['C2']
    # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
)

# construct_calibrant(
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

# construct_calibrant(
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
#     do_reference_stitching=True,
#     cut_to=None
#     # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
# )

# construct_calibrant(
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

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='calibrations/2023-12-26_15-14-38_UV-Vis_ethylaminobutenoate',
    calibrant_shortnames=['EAB'],
    ref_concentrations=[0.0005],
    max_concentrations=[0.00045],
    min_concentrations=[0.0001],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=0.95,
    do_reference_stitching=True,
    cut_to=None,
    calibrant_colors=['C4']
    # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
)

construct_calibrant(
    cut_from=cut_from,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='calibrations/2024-01-22_12-37-31_UV-Vis_bb017',
    calibrant_shortnames=['bb017'],
    ref_concentrations=[0.0003],
    max_concentrations=[0.00055],
    min_concentrations=[0],
    experiment_name=experiment_name,
    upper_limit_of_absorbance=0.95,
    do_reference_stitching=True,
    cut_to=None,
    calibrant_colors=['C5']
    # custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy'
)

# construct_calibrant(
#     cut_from=cut_from,
#     lower_limit_of_absorbance=0.007,
#     concentration_column_name='concentration',
#     do_plot=True,
#     calibration_source_filename='calibrations/2024-01-29_15-09-06_UV-Vis_bb021',
#     calibrant_shortnames=['bb021'],
#     ref_concentrations=[0.0003],
#     max_concentrations=[1],
#     min_concentrations=[0],
#     experiment_name=experiment_name,
#     upper_limit_of_absorbance=0.95,
#     do_reference_stitching=True,
#     cut_to=None,
#     custom_bkg_spectrum_npy_file=data_folder + 'BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/bkg_spectrum.npy',
# )

simpleaxis(plt.gca())

# plt.legend(loc='upper right')
plt.xlabel('Wavelength, nm')
plt.ylabel('Absorbance')
plt.xlim(260, 430)
plt.ylim(-0.05, 0.808)
plt.tight_layout()
f.savefig('misc_scripts/figures/calibration-illustration-hansch-1.png', dpi=300)
plt.show()