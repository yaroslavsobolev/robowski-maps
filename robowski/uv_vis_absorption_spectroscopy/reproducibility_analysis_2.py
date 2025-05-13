from robowski.settings import *
"""
This script is used to analyze the reproducibility of the concentrations of product RF029A in the plates,
in experiments where the same reaction conditions are repeated in all vials of the same plate.
The concentrations are calculated from the absorbance spectra of the plates.
The absorbance spectra are measured with the CRAIC microspectrometer.
"""

from process_wellplate_spectra import SpectraProcessor
import os
import numpy as np
import matplotlib.pyplot as plt


data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
sp = SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'
dilution_factor = 200
diluted_indices = [i + j for i in [9, 27, 45] for j in range(9)]
undiluted_indices = [i + j for i in [0, 18, 36] for j in range(9)]


def process_or_load_processed(nickname, plate_name, run_name, force_recalculation=False):
    """
    Extract the concentration of product RF029A for one plate and save it to a file. If the file exists, load it.

    Parameters
    ----------
    nickname: str
        The name of the file to save and load the concentrations to/from.
        The file will be saved in data_folder + experiment_name + 'results/validation/{nickname}.npy'
    plate_name: str
        The name of the folder with the plate data.
    run_name: str
        The name of the folder with the experiment run. Must end with '/'.
        Example: 'multicomp-reactions/2023-03-20-run01/'
        This is only needed for storing the cached results of the concentrations and the uf-vis spectral fits.
    force_recalculation: bool
        If True, will recalculate the concentrations from spectra even if the file already exists.

    Returns
    -------
    concentrations: np.ndarray
        The concentrations of product RF029A in mol/L.

    """
    # if 'results' folder does not exist, make it
    if not os.path.isdir(data_folder + run_name + 'results/'):
        os.mkdir(data_folder + run_name + 'results/')

    # if 'results/validation' folder does not exist, make it
    if not os.path.isdir(data_folder + run_name + 'results/validation/'):
        os.mkdir(data_folder + run_name + 'results/validation/')

    filepath = data_folder + run_name + 'results/validation/' + f'{nickname}.npy'
    if force_recalculation or not os.path.isfile(filepath):
        concentrations = sp.concentrations_for_one_plate(experiment_folder=data_folder + run_name,
                                                         plate_folder=craic_folder + plate_name + '/',
                                                         calibration_folder=data_folder + 'multicomp-reactions/2023'
                                                                                          '-01-18-run01/' +
                                                         'microspectrometer_data/calibration/',
                                                         calibrant_shortnames=['IIO029A', 'ald001'],
                                                         background_model_folder=data_folder + 'multicomp-reactions'
                                                         '/2023-03-20-run01/microspectrometer_data/background_model/',
                                                         calibrant_upper_bounds=[np.inf, 1e-10],
                                                         do_plot=False)
        np.save(data_folder + run_name + 'results/validation/' + f'{nickname}.npy',
                concentrations)
    else:
        concentrations = np.load(filepath)

    return concentrations * dilution_factor


def print_stats(input_array, label):
    """
    Print the mean, standard deviation and relative standard deviation of an array.

    Parameters
    ----------
    input_array: np.ndarray
        The array to calculate the statistics for.
    label: str
        The label to print before the statistics.

    Returns
    -------
    rel_std: float
        The relative standard deviation of the array.
    """
    print(
        f'{label}: Mean: {np.mean(input_array)}, std: {np.std(input_array)}, '
        f'rel. std: {np.std(input_array) / np.mean(input_array):.1%}')
    return np.std(input_array) / np.mean(input_array)


def view_reproducibility_for_one_plate(plate_name, run_name, with_dilution=True, title=None, force_recaltulation=False,
                                       correction_factor=1):
    """
    View the reproducibility of the concentrations of product RF029A for one plate. It makes two plots:
    1. The concentrations vs. vial ID.
    2. 2D plate layout with colors indicating the concentrations.

    Also prints the statistics for this plate into the terminal:
    mean, standard deviation and relative standard deviation.

    Parameters
    ----------
    plate_name: str
        The name of the folder with the plate data in the CRAIC absorbance folder.
    run_name: str
        The name of the folder with the experiment run. Must end with '/'.
        Example: 'multicomp-reactions/2023-03-20-run01/'
        This is only needed for storing the cached results of the concentrations and the uf-vis spectral fits.
    with_dilution: bool
        If True, will show the concentrations of the diluted vials only.
        If False, will show concentrations in all vials.
    title: str
        The title to put on the plots.
    force_recaltulation: bool
        If True, will recalculate the concentrations from spectra even if the file already exists.

    Returns
    -------
    None
    """
    if title is None:
        title = plate_name
    if with_dilution:
        concentrations = process_or_load_processed(nickname=plate_name, plate_name=plate_name, run_name=run_name,
                                                   force_recalculation=force_recaltulation)[diluted_indices]
    else:
        concentrations = process_or_load_processed(nickname=plate_name, plate_name=plate_name, run_name=run_name,
                                                   force_recalculation=force_recaltulation)
    concentrations = correction_factor * concentrations
    rel_std = print_stats(concentrations, label=title)
    fig, [ax1, ax2] = plt.subplots(2)
    ax1.scatter(np.arange(concentrations.shape[0]), concentrations, label=title)
    ax1.set_title(f'{title}. Relative error {rel_std * 100:.1f}%')
    ax1.set_xlabel('Vial ID')
    ax1.set_ylabel('Product concentration in mol/L')
    ax1.set_ylim(0, max(concentrations) * 1.4)
    ax1.legend()

    if with_dilution:
        as_plate = concentrations.reshape((3, 9))
    else:
        as_plate = concentrations.reshape((6, 9))
    # plt.figure(2)
    ax2.imshow(as_plate)
    # plt.title(title)
    plt.tight_layout()
    fig.savefig(data_folder + run_name + 'results/validation/' + f'{title}.png', dpi=300)
    plt.show()

    return concentrations


if __name__ == '__main__':
    # view_reproducibility_for_one_plate(plate_name='2023-06-15_17-08-49__plate0000039__one-reaction-2023-05-14-run01',
    #                                    run_name='multicomp-reactions/2023-06-14-run02-dmf/',
    #                                    title='2023-06-14-run01. Plate 39. Same reaction in all vials')

    # view_reproducibility_for_one_plate(plate_name='2023-06-16_20-27-22__plate0000042__multicomp-51room_light_35blacked',
    #                                    run_name='multicomp-reactions/2023-03-20-run01/',
    #                                    title='multicomp-51room_light_35blacked')

    # conc = view_reproducibility_for_one_plate(plate_name='2023-06-15_19-32-20__plate0000041__multicomp-2023-05-15-run01-calibration',
    #                                    run_name='multicomp-reactions/2023-06-15-run01/',
    #                                    title='recalibration',
    #                                    with_dilution=False)
    # print(conc)
    # Conclusion: it shows that the concentrations are now being underestimated by 1.56 times

    # view_reproducibility_for_one_plate(plate_name='2023-06-16_15-14-07__plate0000040__multicomp-2023-06-15-run01-repeatability_9-11_corrupted',
    #                                    run_name='multicomp-reactions/2023-06-15-run02/',
    #                                    title='N2_PA12_9-11-corrupted',
    #                                    correction_factor=1.56)
    # view_reproducibility_for_one_plate(plate_name='2023-06-16_15-58-53__plate0000035__multicomp-2023-06-15-run01-repeatability_N2_AL',
    #                                    run_name='multicomp-reactions/2023-06-15-run02/',
    #                                    title='N2_AL',
    #                                    correction_factor=1.56)
    # view_reproducibility_for_one_plate(plate_name='2023-06-16_16-36-57__plate0000036__multicomp-2023-06-15-run01-repeatability_Air_PA12',
    #                                    run_name='multicomp-reactions/2023-06-15-run02/',
    #                                    title='Air_PA12',
    #                                    correction_factor=1.56)
    # view_reproducibility_for_one_plate(plate_name='2023-06-16_16-51-31__plate0000046__multicomp-2023-06-15-run01-repeatability_Air_AL',
    #                                    run_name='multicomp-reactions/2023-06-15-run02/',
    #                                    title='Air_AL',
    #                                    correction_factor=1.56)

    # view_reproducibility_for_one_plate(plate_name='2023-06-17_18-00-02__plate0000046__multicomp_reactions_2023-06-16-run01-reproduction',
    #                                    run_name='multicomp-reactions/2023-06-16-run01/',
    #                                    title='NW_AL_2023-06-17_noextrapolation',
    #                                    correction_factor=1,
    #                                    force_recaltulation=True)

    view_reproducibility_for_one_plate(plate_name='2023-06-13_11-13-01__plate0000037__multicomponent-reactions-2023-06-13-new',
                                       run_name='multicomp-reactions/2023-06-13-run01/',
                                       title='2023-06-13-new',
                                       correction_factor=1,
                                       force_recaltulation=True)