import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import importlib
from scipy import interpolate
import matplotlib.text as mtext
from scipy.optimize import curve_fit
import statsmodels.api as sm

process_wellplate_spectra = importlib.import_module("uv-vis-absorption-spectroscopy.process_wellplate_spectra")
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def spectrum_to_concentration_local(sp, target_spectrum_input, calibration_folder, calibrant_shortnames,
                              background_model_folder,
                              lower_limit_of_absorbance=-0.2, fig_filename='temp', do_plot=False,
                              # lower_limit_of_absorbance=0.02
                              upper_bounds=[np.inf, np.inf], use_line=False, cut_from=200, ignore_abs_threshold=False,
                              cut_to=False, ignore_pca_bkg=False, return_errors=False):  # upper_bounds=[np.inf, np.inf]
    calibrants = []
    for calibrant_shortname in calibrant_shortnames:
        dict_here = dict()
        dict_here['coeff_to_concentration_interpolator'], dict_here['reference_interpolator'], dict_here[
            'bkg_spectrum'] = \
            sp.load_calibration_for_one_calibrant(calibrant_shortname, calibration_folder)
        calibrants.append(dict_here.copy())

    bkg_spectrum = calibrants[0]['bkg_spectrum']
    wavelengths = bkg_spectrum[:, 0]
    target_spectrum = target_spectrum_input - bkg_spectrum[:, 1]
    wavelength_indices = np.arange(calibrants[0]['bkg_spectrum'].shape[0])

    thresh_w_indices = [0, 25, 127, 2000]
    thresh_as = [0.67, 0.75, 1.6, 1.6]
    threshold_interpolator = interpolate.interp1d(thresh_w_indices, thresh_as, fill_value='extrapolate')

    if not ignore_abs_threshold:
        mask = np.logical_and(target_spectrum < threshold_interpolator(wavelength_indices),
                              wavelength_indices > cut_from)
    else:
        mask = wavelength_indices > cut_from

    if cut_to:
        mask = np.logical_and(mask, wavelength_indices <= cut_to)

    mask = np.logical_and(mask,
                          target_spectrum > np.min(target_spectrum) + lower_limit_of_absorbance)

    if not ignore_pca_bkg:
        background_interpolators = [interpolate.interp1d(wavelength_indices,
                                                         np.load(background_model_folder + f'component_{i}.npy'),
                                                         fill_value='extrapolate')
                                    for i in range(2)]
    else:
        background_interpolators = [interpolate.interp1d(wavelength_indices,
                                                         np.ones_like(wavelength_indices),
                                                         fill_value='extrapolate')
                                    for i in range(2)]

    if len(wavelength_indices[mask]) == 0:
        print('There is no data that is within mask. Returning zeros.')
        return [0 for i in range(4)]

    if len(calibrant_shortnames) == 2:
        def func(xs, a, b, c, d, e, f):
            return a * calibrants[0]['reference_interpolator'](xs) + b * calibrants[1]['reference_interpolator'](xs) + c \
                   + d * xs + e * background_interpolators[0](xs) + f * background_interpolators[1](xs)
    elif len(calibrant_shortnames) == 3:
        def func(xs, a1, a2, a3, c, d, e, f):
            return a1 * calibrants[0]['reference_interpolator'](xs) + \
                   a2 * calibrants[1]['reference_interpolator'](xs) + \
                   a3 * calibrants[2]['reference_interpolator'](xs) \
                   + c + d * xs + e * background_interpolators[0](xs) + f * background_interpolators[1](xs)
    else:
        raise NotImplementedError

    # p0 = tuple(0.5 if upper_bounds[0] is np.inf else upper_bounds[0],
    #       0.5 if upper_bounds[1] is np.inf else upper_bounds[1],
    #       0,
    #       0,
    #       0,
    #       0)
    p0 = tuple([0.5 if upper_bound is np.inf else upper_bound for upper_bound in upper_bounds] + [0] * 4)
    if use_line:
        linebounds = [-np.inf, np.inf]
    else:
        linebounds = [-1e-15, 1e-15]

    if ignore_pca_bkg:
        bkg_comp_limit = 1e-12
    else:
        bkg_comp_limit = np.inf
    bounds = ([-1e-20] * len(calibrant_shortnames) + [-np.inf, linebounds[0], -1 * bkg_comp_limit, -1 * bkg_comp_limit],
              upper_bounds + [np.inf, linebounds[1], bkg_comp_limit, bkg_comp_limit])
    popt, pcov = curve_fit(func, wavelength_indices[mask], target_spectrum[mask],
                           p0=p0, bounds=bounds)
    perr = np.sqrt(np.diag(pcov))  # errors of the fitted coefficients

    #plot covariance matrix
    fign = plt.figure(figsize=(2, 2.3))
    # convert covariance matrix pcov to correlation matrix
    pcov = pcov / np.outer(perr, perr)
    pcov_to_plot = pcov[:len(calibrant_shortnames), :len(calibrant_shortnames)]
    plt.imshow(pcov_to_plot, vmin=-1*max(np.abs(pcov_to_plot).flatten()), vmax=max(np.abs(pcov_to_plot).flatten()),
               cmap='RdBu_r')
    # make tick labels from calibrant_shortnames
    csnames = ['A', 'B', 'C']
    csnames = ['', '', '']
    plt.yticks(range(len(calibrant_shortnames)), csnames)
    plt.xticks(range(len(calibrant_shortnames)), csnames, rotation=90)
    plt.colorbar(orientation='horizontal', fraction=0.046*2, aspect=10)
    plt.tight_layout()
    fign.savefig('misc-scripts/figures/correlation_matrix.png', dpi=300)
    plt.show()

    concentrations_here = [calibrants[calibrant_index]['coeff_to_concentration_interpolator'](fitted_coeff)
                           for calibrant_index, fitted_coeff in enumerate(popt[:-4])]

    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.9, 4.5), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    ax = ax1
    # ax.plot(wavelengths, target_spectrum_input, label='Raw data', color='grey', alpha=0.2)
    line_experimental = ax.plot(wavelengths, target_spectrum, label='Measured spectrum\nof crude mixture', color='black', linewidth=3, alpha=0.5)
    mask_illustration = np.ones_like(target_spectrum) * np.max(target_spectrum)
    mask_illustration[mask] = 0
    # ax.fill_between(x=wavelengths, y1=0, y2=mask_illustration, color='yellow', alpha=0.3,
    #                 label='Masked data')
    # ax.plot(wavelengths, func(wavelength_indices, *popt), color='r', label='Unmixing model', alpha=0.5)
    label_calibrant_shortnames = ['Methoxychalcone\ncontribution', 'Anisaldehyde\ncontribution', 'Acetophenone\ncontribution']

    lines_components = []
    for calibrant_index in range(len(calibrant_shortnames)):
        cpopt = [x if i == calibrant_index else 0 for i, x in enumerate(popt)]
        lines_components.append(
            ax.plot(wavelengths, func(wavelength_indices, *cpopt), label=label_calibrant_shortnames[calibrant_index], alpha=0.5)
        )
        ax.fill_between(x=wavelengths, y1=0, y2=func(wavelength_indices, *cpopt), alpha=0.1)
    # make a list where only the third from the end item is the same as in popt, while the other ones are zero
    if use_line:
        cpopt = [x if i == len(popt) - 3 else 0 for i, x in enumerate(popt)]
        ax.plot(wavelengths, func(wavelength_indices, *cpopt), label='Line', alpha=0.5)
    if not ignore_pca_bkg:
        cpopt = [x if i == len(popt) - 2 else 0 for i, x in enumerate(popt)]
        ax.plot(wavelengths, func(wavelength_indices, *cpopt), label='Bkg. PC1', alpha=0.5)
        cpopt = [x if i == len(popt) - 1 else 0 for i, x in enumerate(popt)]
        ax.plot(wavelengths, func(wavelength_indices, *cpopt), label='Bkg. PC2', alpha=0.5)
    # plt.ylim(-0.3,
    #          np.max((func(wavelength_indices, *popt)[mask])) * 3)
    # title_str = f'Concentrations:\n'
    # for i in range(len(concentrations_here)):
    #     title_str += f'{np.array(concentrations_here)[i]:.6f} M ({calibrant_shortnames[i]})\n '

    ax.set_ylim(0, 0.4)
    ax.set_xlim(min(wavelengths[mask]), 420)

    # fig1.suptitle(title_str[:-2])
    ax.set_ylabel('Absorbance')
    # ax.legend(loc='upper right')
    simpleaxis(ax)

    # class LegendTitle(object):
    #     def __init__(self, text_props=None):
    #         self.text_props = text_props or {}
    #         super(LegendTitle, self).__init__()
    #
    #     def legend_artist(self, legend, orig_handle, fontsize, handlebox):
    #         x0, y0 = handlebox.xdescent, handlebox.ydescent
    #         title = mtext.Text(x0, y0, r'\underline{' + orig_handle + '}', usetex=True, **self.text_props)
    #         handlebox.add_artist(title)
    #         return title
    #
    # basestring = 'Base'
    # ax.legend([line_experimental[0], 'Individual components:', lines_components[0][0], lines_components[1][0], lines_components[2][0]],
    #            ['Measured spectrum of mixture', '', 'Methoxychalcone', 'Anisaldehyde', 'Acetophenone'],
    #            handler_map={basestring: LegendTitle({'fontsize': 18})})




    # Residuals subplot
    ax = ax2
    ax.plot(wavelengths[mask], target_spectrum[mask] - func(wavelength_indices[mask], *popt), color='black', alpha=0.5,
            label='residuals')
    ax.fill_between(x=wavelengths[mask], y1=0, y2=target_spectrum[mask] - func(wavelength_indices[mask], *popt), color='black', alpha=0.15)
    # ax.legend()
    ax.set_xlabel('Wavelength, nm')
    ax.set_ylabel('Residuals,\nabsorbance units')

    simpleaxis(ax)

    plt.tight_layout()
    fig1.savefig(f"{fig_filename}.png", dpi=300)

    residuals_here = target_spectrum[mask] - func(wavelength_indices[mask], *popt)
    lag = 50
    lb_df = sm.stats.acorr_ljungbox(residuals_here, lags=[lag])
    # print(f'spectrum index {spectrum_index}, LB_pvalue: {lb_df.loc[lag, "lb_pvalue"]}, lag: {lag}')
    # print(lb_df)
    if len(lb_df) == 1:
        # take values from first row of dataframe lb_pvalue
        print(f"LB_pvalue : {lb_df.loc[lag, 'lb_pvalue']}, stat: {lb_df.loc[lag, 'lb_stat']}")
    else:
        print('hmmm')

    if do_plot:
        plt.show()
    else:
        plt.close(fig1)
        plt.close('all')
        plt.clf()

    if return_errors:
        # convert coefficient errors into concentration errors
        upper_confidence_limit = [
            calibrants[calibrant_index]['coeff_to_concentration_interpolator'](fitted_coeff + perr[calibrant_index])
            for calibrant_index, fitted_coeff in enumerate(popt[:-4])]
        concentration_errors = [upper_confidence_limit[i] - concentrations_here[i] for i in
                                range(len(concentrations_here))]
        return concentrations_here, concentration_errors

    return concentrations_here


def process_plate(sp, dilution_factor,
                  plate_folder,
                  well_ids,
                  calibrant_shortnames,
                  calibration_folder,
                  experiment_name,
                  cut_from,
                  cut_to,
                  do_plot=False,
                  use_line=False,
                  ignore_bkg_pca=True,
                  get_errors_from_fit=False,
                  std_calib = 0.0105,
                  upper_bounds='auto'):
    if upper_bounds == 'auto':
        upper_bounds = [np.inf] * len(calibrant_shortnames)
    plate_name = plate_folder.split('/')[-1]
    df = pd.DataFrame(columns=['well_id'] + calibrant_shortnames)
    for well_id in well_ids:
        spectrum = sp.load_msp_by_id(
            plate_folder=plate_folder,
            well_id=well_id)[:, 1]
        process_wellplate_spectra.create_folder_unless_it_exists(data_folder + experiment_name + 'results')
        process_wellplate_spectra.create_folder_unless_it_exists(data_folder + experiment_name + 'results/uv-vis-fits')
        concentrations_here = spectrum_to_concentration_local(sp, target_spectrum_input=spectrum,
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

experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Claisen_WaiShing/'

sp = process_wellplate_spectra.SpectraProcessor(
    folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                   '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220
calibration_folder = data_folder + experiment_name + 'microspectrometer_data/calibration/'
process_plate(sp, dilution_factor=500,
              plate_folder=f'{data_folder}{experiment_name}2023-10-12_14-51-59_UV-Vis_crude.csv',
              well_ids=range(5),
              cut_from=5,
              cut_to=False,
              calibrant_shortnames=['methoxychalcone', 'anisaldehyde', 'acetophenone'],
              calibration_folder=calibration_folder,
              experiment_name=experiment_name,
              do_plot=True)



