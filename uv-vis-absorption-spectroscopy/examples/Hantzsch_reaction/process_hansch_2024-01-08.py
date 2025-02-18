import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import importlib
process_wellplate_spectra = importlib.import_module("uv-vis-absorption-spectroscopy.process_wellplate_spectra")
organize_run_results = importlib.import_module("misc-scripts.organize_run_results")
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
plt.ioff()

import logging
logging.basicConfig(level=logging.INFO)

def process_run_by_shortname(run_shortname, cut_from=5, dilution_factor=200):
    substrates = ['ethyl_acetoacetate', 'methoxybenzaldehyde', 'ammonium_acetate']
    product_name = 'HRP01'
    run_name = f'BPRF/{run_shortname}/'
    substances_for_fitting = ['methoxybenzaldehyde', 'HRP01']

    sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                                         '2022-12-01/interpolator-dataset/')
    sp.nanodrop_lower_cutoff_of_wavelengths = 220

    df_structure = organize_run_results.load_run_structure(run_name)
    for nanodrop_filepath in df_structure.nanodrop_filepath.unique():
        print('Processing nanodrop filepath', nanodrop_filepath)
        concentrations_here = sp.concentrations_for_one_plate(experiment_folder=data_folder + run_name,
                                                              plate_folder=data_folder + nanodrop_filepath,
                                                              calibration_folder=data_folder + 'BPRF/2023-11-08-run01/' + 'microspectrometer_data/calibration/',
                                                              calibrant_shortnames=substances_for_fitting,
                                                              background_model_folder=data_folder + 'simple-reactions/2023-09-06-run01/microspectrometer_data/background_model/',
                                                              calibrant_upper_bounds=[np.inf, np.inf],
                                                              do_plot=False, return_all_substances=True, cut_from=cut_from,
                                                              ignore_abs_threshold=True, ignore_pca_bkg=True)
        concentrations_here = concentrations_here * dilution_factor
        for vial_id, product_concentrations in enumerate(concentrations_here):
            # index of this vial in the concentrations_df dataframe
            try:
                index_of_this_vial = df_structure.index[(df_structure['container_id'] == vial_id)
                                                        & (df_structure['nanodrop_filepath'] == nanodrop_filepath)][0]
            except IndexError:
                logging.warning(f'Container_id {vial_id} not present in conditions dataframe but was measured in nanodrop file {nanodrop_filepath}')
                continue

            for i, substance_name in enumerate(substances_for_fitting):
                df_structure.loc[index_of_this_vial, f'pc#{substance_name}'] = product_concentrations[i]

    logging.debug(f'global id min {df_structure["global_index"].min()}')
    concentrations_df = df_structure[df_structure[f'pc#{product_name}'].notna()]
    logging.debug(f'cdf global id min {concentrations_df["global_index"].min()}')

    # Calculate yields
    for index, row in concentrations_df.iterrows():
        substrate_concentrations = [concentrations_df.loc[index, f'c#{substrate_name}'] for substrate_name in substrates]
        # multiply the concentration of acetoacetate by 2 because two moles are needed to produce one mole of product
        substrate_concentrations[0] = substrate_concentrations[0] * 2
        substrate_concentration_min = min(substrate_concentrations)
        product_concentration = concentrations_df.loc[index, f'pc#{product_name}']
        # if the product concentration is negative, replace it with zero
        product_concentration = max(product_concentration, 0)
        concentrations_df.loc[index, 'yield'] = product_concentration / substrate_concentration_min

    concentrations_df.to_csv(data_folder + run_name + 'results/' + 'product_concentration.csv', index=False)

def plot_all_spectra_by_shortname(run_shortname):
    run_name = f'BPRF/{run_shortname}/'

    sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                                         '2022-12-01/interpolator-dataset/')
    sp.nanodrop_lower_cutoff_of_wavelengths = 220

    df_structure = organize_run_results.load_run_structure(run_name)
    for nanodrop_filepath in df_structure.nanodrop_filepath.unique():
        sp.show_all_spectra(plate_folder=data_folder + nanodrop_filepath)

if __name__ == '__main__':
    # list_of_runs = tuple([
    #                       '2023-11-08-run01',
    #                       '2023-11-13-run01',
    #                       '2023-11-14-run01',
    #                       '2023-11-21-run01'
    # ])
    list_of_runs = tuple(['2023-12-27-run01_200',
                          '2023-12-27-run02_200',
                          '2023-12-27-run03_200'])

    for i, run_shortname in enumerate(list_of_runs):
        # process_run_by_shortname(run_shortname)
        plot_all_spectra_by_shortname(run_shortname)

    plt.show()