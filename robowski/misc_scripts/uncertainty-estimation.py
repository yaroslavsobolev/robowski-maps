from robowski.settings import *
'''
Estimate errors (variation) of product concentration measurement at certain experimental conditions
using linear theory of uncertainty propagation. Units are arbitrary, but must be consistent.
'''

import os
from uncertainties import ufloat
from uncertainties.umath import *

import matplotlib.pyplot as plt
import numpy as np

# TODO: For Yankai to fix. This module `calibration_data` does not exist, and the function zeus_uncertainty_from_file() does not exist anywhere in the repository.
calibration_data = importlib.import_module("zeus-pipetter.calibration_data")

# Load and process the volume validation measurements from Zeus
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

def get_error_interpolators_from_all_measurement_files(zeus_volume_validation_measurements_files, tip_types):
    error_interpolators = dict()
    for i, tfile in enumerate(zeus_volume_validation_measurements_files):
        tip_type, statistical_data, error_interpolator = calibration_data.zeus_uncertainty_from_file(tfile, only_tip_type=tip_types[i])
        error_interpolators[tip_type] = error_interpolator
    return error_interpolators

# error_interpolators = get_error_interpolators_from_all_measurement_files([data_folder + \
#     'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_2023_03_26_02_13_50ul.json',
#     data_folder +
#     'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_2023_03_26_15_18_300ul.json'],
#     tip_types=[50, 300])
#
# error_interpolators = get_error_interpolators_from_all_measurement_files([data_folder + \
#     'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_2023_03_27_02_49_50ul_and_300ul.json',
#     data_folder +
#     'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_2023_03_27_02_49_50ul_and_300ul.json'],
#     tip_types=[50, 300])

error_interpolators = get_error_interpolators_from_all_measurement_files([data_folder + \
    'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_2023_03_29_03_08_50ul_recalib.json',
    data_folder +
    'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_2023_03_29_01_13_300ul_recalib.json',
      data_folder +
      'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_03_29_13_50_1000ul_tips_recalib_3rd.json'
      ],
    tip_types=[50, 300, 1000])

def zeus_uncertainty_vs_volume(volume, thresholds_for_tip_choice = [50, 600]):
    '''Assumes that error interpolators were loaded globally into the variable error_interpolators'''
    if volume < thresholds_for_tip_choice[0]:
        error = error_interpolators[50](volume)
    elif volume < thresholds_for_tip_choice[1]:
        error = error_interpolators[300](volume)
    else:
        error = error_interpolators[1000](volume)
    return error

# Preparation of the reaction mixture by pipetting the components into the reaction vessel
reactants = [ufloat(volume, zeus_uncertainty_vs_volume(volume)) for volume in [30, 30, 30, 30]]
additional_solvent_volume = 200 - sum(reactants).n
total_reaction_volume = sum(reactants) + ufloat(additional_solvent_volume, zeus_uncertainty_vs_volume(additional_solvent_volume))
reaction_time = ufloat(16*60*60, 20*60) # in seconds

# Rinetics of the reaction. Here, linear approximation of fourth-order kinetics.
# That is, rate is proportional to product of four concentrations.

# Arrhenius law assuming twofold increase in reaction rate when temperature changes by 10 Kelvin
temperature_flucturation_std = 1 # deg K
temperature_dependence = 2**(ufloat(0, temperature_flucturation_std)/10)

product_concentration = reactants[0] * reactants[1] * reactants[2] * reactants[3] / (total_reaction_volume ** 4) * reaction_time * temperature_dependence
print(f'product concentrations after reaction: {product_concentration:.2eP}, relative error: {product_concentration.s / product_concentration.n * 100:.1f}%')

# Dilution sequence
dilution_transfer_volumes = [ufloat(volume, zeus_uncertainty_vs_volume(volume)) for volume in [1000, 15, 485]]
dilution_transfer_volumes[0] = ufloat(1000, 2)
dilution_transfer_volumes[2] = ufloat(485, error_interpolators[1000](485))

# # single std for second operation -- to account for tip changing at each operation
# dilution_transfer_volumes[1] = ufloat(15, 0.9)

concentration_in_first_well_after_first_operation = product_concentration * total_reaction_volume / \
                                                    (total_reaction_volume + dilution_transfer_volumes[0])
concentration_in_second_well_after_three_operations = concentration_in_first_well_after_first_operation * \
                                                      dilution_transfer_volumes[1] / \
                                                      (dilution_transfer_volumes[1] + dilution_transfer_volumes[2])
dilution_factor = concentration_in_second_well_after_three_operations / product_concentration
print(f'dilution factor: {dilution_factor:.2eP}, relative error: {dilution_factor.s / dilution_factor.n * 100:.1f}%')

# Measurement
# vial_cross_section = ufloat(1, 0.02) # in mm
opticalpathlength = ufloat(1, 0.02) # from vial geometry variation alone
vial_cross_section = 1 # in mm^2
absorbance = concentration_in_second_well_after_three_operations * \
             (dilution_transfer_volumes[1] + dilution_transfer_volumes[2]) / vial_cross_section * opticalpathlength

print(f'absorbance: {absorbance:.2eP}, relative error: {absorbance.s / absorbance.n * 100:.1f}%')

# 50 uL tip, no pre-wet, 15 uL target
tip50ul_target_15ul_data = [14.31, 14.42, 14.39, 14.54, 14.34, 14.49, 14.42, 14.51, 14.56, 14.45, 14.51, 14.46, 14.23,
                            14.31, 14.31, 14.38, 14.48, 14.42, 14.07, 13.97, 14.04, 13.79, 13.88, 13.52, 13.75, 13.91,
                            14.04, 13.6, 13.42, 13.43]
fig1 = plt.figure(1)
plt.errorbar(x=0, y=np.mean(tip50ul_target_15ul_data), yerr=np.std(tip50ul_target_15ul_data), fmt='o', label='no pre-wetting')

tip_type, statistical_data, error_interpolator = calibration_data.zeus_uncertainty_from_file(data_folder + \
    'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_2023_03_27_02_49_50ul_and_300ul.json',
                                                                                             only_tip_type=50, do_plot=False)
plt.errorbar(x=1, y=statistical_data[4, 1], yerr=statistical_data[4, 2], fmt='o', label='with pre-wetting',
             capsize=6)
plt.axhline(y=15)
plt.show()
# plt.errorbar(x=1, y=zeus_uncertainty_vs_volume(volume)