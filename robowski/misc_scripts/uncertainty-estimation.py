'''
Estimate errors (variation) of product concentration measurement at certain experimental conditions
using linear theory of uncertainty propagation. Units are arbitrary, but must be consistent.
'''

import os, json, pandas as pd
import time
from textwrap import wrap
from scipy import interpolate
from uncertainties import ufloat
from uncertainties.umath import *
import importlib
import matplotlib.pyplot as plt
import numpy as np

# Load and process the volume validation measurements from Zeus
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'


def zeus_uncertainty_from_file(tfile, only_tip_type, do_plot=True):
    '''
    Estimate the uncertainty of the volume measurement from the data file
    Parameters
    ----------
    tfile: str
        Path to the file with the data

    only_tip_type: int
        Tip type defined by maximum volume (in microliters). Only this tip type will be loaded from the file.

    do_plot: bool
        Whether to plot the results

    Returns
    -------
    (tip_type, np.array(statistical_data), error_interpolator): tuple
    where
        tip_type: int
            Type of the tip, described by maximum volume in microliters
        statistical_data: numpy.ndarray (3xN)
            The statistical data with columns: target volume, mean measured volume, std (random error) of measured volume
        error_interpolator: scipy.interpolate.interp1d
            The interpolator for the overall error as a function of the target volume
    '''

    # This data structure is a madman's magnum opus. Observe this insanity:
    with open(tfile) as file_handler:
        data = json.load(file_handler)
        data = list(data.values())[0]
    statistical_data = []
    diffs = []
    for entry in data:
        the_only_key = list(entry.keys())[0]
        header = the_only_key.split('_')
        target_volume = int(header[-1][:-2])
        tip_type = int(header[-2][:-2])
        if tip_type != only_tip_type:
            continue
        measured_volumes = entry[the_only_key]['volume']
        # Yankai implored to remove the first point because if the wetting problem/artifact
        # measured_volumes = measured_volumes[:-1]
        measured_volumes = measured_volumes[1:]
        statistical_data.append([target_volume, np.mean(measured_volumes), np.std(measured_volumes)])
        diffs.extend([[target_volume, x - target_volume] for x in measured_volumes])
    statistical_data = np.array(statistical_data)
    df = pd.DataFrame(statistical_data, columns=['target_volume', 'measured_mean', 'measured_std'])
    df.to_csv(tfile.replace('.json', f'_processed_for_{only_tip_type}ul_tiptype.csv'), index=False)
    target_volumes, measured_volumes, measured_std = statistical_data[:, 0], statistical_data[:, 1], statistical_data[:, 2]

    # Optional plotting
    if do_plot:
        fig, ax = plt.subplots(2, 1, sharex=True)
        wrapped_filename = "\n".join(wrap(tfile,60))
        ax[0].set_title(f'Tip type: {tip_type} $\mu$L,\n file: {wrapped_filename}', wrap=True)
        ax[0].errorbar(x=target_volumes, y=measured_volumes, yerr=measured_std,
                       fmt='o-', markersize=3, capsize=5, alpha=0.5)
        ax[0].plot([np.min(target_volumes), np.max(target_volumes)],
                 [np.min(target_volumes), np.max(target_volumes)], color='black')
        ax[1].errorbar(x=target_volumes, y=measured_volumes-target_volumes, yerr=measured_std,
                       fmt='o-', markersize=5, capsize=8, alpha=0.7)
        ax[1].axhline(y=0, color='black')
        ax[1].scatter(np.array(diffs)[:, 0], np.array(diffs)[:, 1], alpha=0.2, color='C1', marker='x')
        plt.xlabel('Intended volume, $\mu$L')
        ax[1].set_ylabel('Measured minus\nintended, $\mu$L')
        ax[0].set_ylabel('Measured volume, $\mu$L')
        plt.tight_layout()
        fig.savefig(tfile.replace('.json', f'_processed_for_{only_tip_type}ul_tiptype.png'), dpi=300)
        plt.show()
    systematic_errors = measured_volumes - target_volumes
    overall_errors = np.sqrt(systematic_errors ** 2 + measured_std ** 2)
    error_interpolator = interpolate.interp1d(target_volumes, overall_errors, fill_value='extrapolate', kind='linear')
    return only_tip_type, np.array(statistical_data), error_interpolator

def get_error_interpolators_from_all_measurement_files(zeus_volume_validation_measurements_files, tip_types):
    error_interpolators = dict()
    for i, tfile in enumerate(zeus_volume_validation_measurements_files):
        tip_type, statistical_data, error_interpolator = zeus_uncertainty_from_file(tfile, only_tip_type=tip_types[i])
        error_interpolators[tip_type] = error_interpolator
    return error_interpolators

error_interpolators = get_error_interpolators_from_all_measurement_files([
    data_folder + 'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_2023_03_29_03_08_50ul_recalib.json',
    data_folder + 'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_2023_03_29_01_13_300ul_recalib.json',
    data_folder + 'multicomp-reactions/pipetter_io/measured_volumes/calibration_results_03_29_14_12_1000ul_tips_recalib.json'
      ],
    tip_types=[50, 300, 1000])

def zeus_uncertainty_vs_volume(volume, thresholds_for_tip_choice = [50, 500]):
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

# Kinetics of the reaction. Here, linear approximation of fourth-order kinetics.
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

