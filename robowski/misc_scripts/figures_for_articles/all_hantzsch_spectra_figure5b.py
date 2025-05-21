from robowski.settings import *
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra

if __name__ == '__main__':
    sp = process_wellplate_spectra.SpectraProcessor(
        folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                        '2022-12-01/interpolator-dataset/')
    sp.nanodrop_lower_cutoff_of_wavelengths = 220
    sp.use_instrumental_sigmas = True
    substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'dm36', 'dm37', 'dm40_12', 'dm40_10',
                              'ethyl_acetoacetate', 'EAB', 'bb017', 'bb021', 'dm70', 'dm053', 'dm088_4', 'bb021_f2']
    dictionary_for_converting_from_shortnames_to_manuscript_codes = {
        'methoxybenzaldehyde': '19c',
        'HRP01': '19d',
        'dm35_8': '19m',
        'dm35_9': '19k',
        'dm36': '19o',
        'dm37': '19r',
        'dm40_12': '19f',
        'dm40_10': '19g',
        'ethyl_acetoacetate': '19a',
        'EAB': 'EAB',
        'bb017': '19e',
        'bb021': '19h',
        'dm70': '19p',
        'dm053': '19j',
        'dm088_4': '19i',
        'bb021_f2': '(E,E) isomer of 19h'}
    concentrations = sp.plot_calibrant_references(
        calibration_folder=data_folder + 'BPRF/2024-01-17-run01/' + 'microspectrometer_data/calibration/',
        calibrant_shortnames=substances_for_fitting,
        calibrant_labels_on_the_plot=dictionary_for_converting_from_shortnames_to_manuscript_codes)
