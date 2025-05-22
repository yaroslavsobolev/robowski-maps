from robowski.settings import *
import numpy as np
import pandas as pd
import robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra as process_wellplate_spectra


def process_plate(sp, dilution_factor,
                  plate_folder,
                  well_ids,
                  calibrant_shortnames,
                  calibration_folder,
                  experiment_name,
                  cut_from,
                  cut_to,
                  do_plot=False,
                  use_line=True,
                  ignore_bkg_pca=True,
                  get_errors_from_fit=False,
                  std_calib = 0.0105,
                  upper_bounds='auto'):
    if upper_bounds == 'auto':
        upper_bounds = [np.inf] * len(calibrant_shortnames)
    plate_name = plate_folder.split('/')[-1]
    df = pd.DataFrame(columns=['well_id'] + calibrant_shortnames, dtype=object)
    for well_id in well_ids:
        spectrum = sp.load_msp_by_id(
            plate_folder=plate_folder,
            well_id=well_id)[:, 1]
        process_wellplate_spectra.create_folder_unless_it_exists(data_folder + experiment_name + 'results')
        process_wellplate_spectra.create_folder_unless_it_exists(data_folder + experiment_name + 'results/uv-vis-fits')
        concentrations_here = sp.spectrum_to_concentration(target_spectrum_input=spectrum,
                                                             calibration_folder=calibration_folder,
                                                             calibrant_shortnames=calibrant_shortnames,
                                                             fig_filename=data_folder + experiment_name + f'results/uv-vis-fits/{plate_name}-well{well_id:02d}',
                                                             do_plot=do_plot,
                                                             background_model_folder=data_folder + 'multicomp-reactions/2023-03-20-run01/microspectrometer_data/background_model/',
                                                             upper_bounds=upper_bounds,
                                                             cut_from=cut_from,
                                                             cut_to=cut_to,
                                                             ignore_abs_threshold=True,
                                                             ignore_pca_bkg=ignore_bkg_pca,
                                                             use_line=use_line,
                                                             return_errors=get_errors_from_fit)
        if get_errors_from_fit:
            concentrations_here, concentration_errors = concentrations_here
            concentrations_here = np.array(concentrations_here) * dilution_factor
            concentration_errors = np.array(concentration_errors) * dilution_factor
        else:
            concentrations_here = np.array(concentrations_here) * dilution_factor
        print(f'Well {well_id:02d} concentrations: {[calibrant_shortnames[i] + ": " + str(concentrations_here[i]) for i in range(len(calibrant_shortnames))]}')
        # add to dataframe
        df.loc[len(df)] = [well_id] + list(concentrations_here)
        # print the mean ± std

    for column in df.columns[1:]:
        factor = 1000

        if not get_errors_from_fit:
            std_here = df[column].std()
        else:
            std_here = concentration_errors[calibrant_shortnames.index(column)]
        mean_here = df[column].mean()
        # additional relative error due to uncertainty of calibration solutions
        std_here = mean_here * ( (std_here / mean_here) + std_calib )

        print(f'{column}: ({factor * df[column].mean():.2f} ± {factor * std_here:.2f}) mM, n = {len(df)}')
    print(df)

    return df