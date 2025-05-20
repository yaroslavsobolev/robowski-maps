from robowski.settings import *
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra

sp = process_wellplate_spectra.SpectraProcessor(
    folder_with_correction_dataset=repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                    '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 250
sp.use_instrumental_sigmas = False
dilution_factor = 200

##### This constructs the calibration for the product 'IIO029A' and saves for later. Do not rerun unless you know what you do. #######
sp.construct_reference_for_calibrant(calibrant_shortname='IIO029A',
                                     calibration_folder=data_folder + 'multicomp-reactions/2023-01-18-run01/' + 'microspectrometer_data/calibration/',
                                     ref_concentration=0.00011,
                                     do_plot=True, do_reference_refinements=True)

#### This constructs the calibration for the substrate 'ald001' and saves for later. Do not rerun unless you know what you do. #######
sp.construct_reference_for_calibrant(calibrant_shortname='ald001',
                                     calibration_folder=data_folder + 'multicomp-reactions/2023-01-18-run01/' + 'microspectrometer_data/calibration/',
                                     ref_concentration=0.0192096,
                                     do_plot=True, do_reference_refinements=False)