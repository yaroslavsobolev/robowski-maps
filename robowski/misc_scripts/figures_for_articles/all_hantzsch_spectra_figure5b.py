from robowski.settings import *
import numpy as np
import matplotlib.pyplot as plt
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra

if __name__ == '__main__':
    sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                         '2022-12-01/interpolator-dataset/')
    sp.nanodrop_lower_cutoff_of_wavelengths = 220
    sp.use_instrumental_sigmas = True

    # print('>>>>>>>>>')
    # print(sp.uncertainty_of_measured_absorbance(221, 0.37))

    # x = sp.load_nanodrop_csv_for_one_plate(plate_folder=data_folder + 'BPRF/2024-01-08-run01/nanodrop_spectra/2024-01-10_12-51-07_UV-Vis_plate_71.csv')

    # well_id = 44
    well_id = 9
    substances_for_fitting = ['methoxybenzaldehyde', 'HRP01', 'dm35_8', 'dm35_9', 'dm36', 'dm37', 'dm40_12', 'dm40_10', 'ethyl_acetoacetate', 'EAB', 'bb017', 'bb021', 'dm70', 'dm053', 'dm088_4', 'bb021_f2']
    # cut_from = 40
    cut_from = 0
    # Condition 154
    # plate_folder = data_folder + 'BPRF/2024-01-08-run01/nanodrop_spectra/2024-01-10_12-51-07_UV-Vis_plate_71.csv'
    # plate_folder = data_folder + 'BPRF/2024-01-08-run02/nanodrop_spectra/2024-01-10_17-10-28_UV-Vis_plate_61.csv'
    # plate_folder = data_folder + 'BPRF/2024-01-17-run01/nanodrop_spectra/2024-01-19_12-22-47_UV-Vis_plate_66.csv'
    plate_folder = data_folder + 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-22-22_UV-Vis_plate_92.csv'
    # plate_folder = data_folder + 'BPRF/2024-02-16-run01/nanodrop_spectra/2024-02-18_17-48-07_UV-Vis_plate74.csv'
    spectrum1 = sp.load_msp_by_id(
        plate_folder=plate_folder,
        well_id=well_id)[:, 1]

    # plate_folder = data_folder + 'BPRF/2024-01-08-run01/nanodrop_spectra/2024-01-10_13-48-13_UV-Vis_plate_73.csv'
    # plate_folder = data_folder + 'BPRF/2024-01-08-run02/nanodrop_spectra/2024-01-10_17-55-20_UV-Vis_plate_66.csv'
    # plate_folder = data_folder + 'BPRF/2024-01-17-run01/nanodrop_spectra/2024-01-19_13-00-17_UV-Vis_plate_67.csv'
    plate_folder = data_folder + 'BPRF/2024-03-06-run01/nanodrop_spectra/2024-03-08_10-44-22_UV-Vis_plate_93.csv'
    # plate_folder = data_folder + 'BPRF/2024-02-16-run01/nanodrop_spectra/2024-02-18_18-02-42_UV-Vis_plate76.csv'
    spectrum2 = sp.load_msp_by_id(
        plate_folder=plate_folder,
        well_id=well_id)[:, 1]

    plt.plot(spectrum1*20)
    plt.plot(spectrum2*200)
    # plt.plot(spectrum1)
    # plt.plot(spectrum2)
    print('len of spectrum1', len(spectrum1))
    plt.show()

    concentrations = sp.multispectrum_to_concentration(target_spectrum_inputs=[spectrum1, spectrum2],
                                                       dilution_factors=[20, 200],
                                                       calibration_folder=data_folder + 'BPRF/2024-01-17-run01/' + 'microspectrometer_data/calibration/',
                                                       calibrant_shortnames=substances_for_fitting,
                                                       background_model_folder=data_folder + 'BPRF/cross_conamination_and_backgound_test/ethanol_background_model/',
                                                       upper_bounds=[np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf],
                                                       do_plot=False, cut_from=cut_from, cut_to=250,
                                                       ignore_abs_threshold=False, ignore_pca_bkg=False,
                                                       plot_calibrant_references=True,
                                                       upper_limit_of_absorbance=0.95,
                                                       obey_stoichiometric_inequalities=False)

    print(concentrations)
    print('Done!')