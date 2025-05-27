import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit

# we have measured spectra of pure methoxychalcone at different concentrations in ethanol.
# The spectra are organized as the CSV file with numerated columns. Numeration (naming) of the columns does not have to be in order.
# The first column (wavelengths) is unnamed.
# The beginning of the `calibration_of_all_components.csv` looks like this:
# ,0,1,2,3,4,7,5,6,9,8,11,10,13,12,15,14
# 220,0.068,0.102,0.194,0.281,0.015,0.353,0.108,0.187,0.673,0.555,0.049,0.033,0.114,0.078,0.02,0.139
# 221,0.068,0.105,0.2,0.286,0.015,0.338,0.104,0.179,0.631,0.522,0.054,0.035,0.13,0.089,0.017,0.16
# 222,0.066,0.105,0.203,0.291,0.012,0.324,0.101,0.173,0.593,0.494,0.06,0.037,0.149,0.102,0.018,0.182
# 223,0.068,0.108,0.207,0.297,0.011,0.316,0.098,0.168,0.581,0.481,0.064,0.039,0.165,0.112,0.016,0.201
# 224,0.068,0.109,0.209,0.304,0.01,0.3,0.093,0.159,0.555,0.459,0.069,0.041,0.183,0.123,0.016,0.223

# We keep notes of the concentrations in a separate file `calibration_of_all_components.txt`, which is
# a table that looks like this:
#
# nanodrop_col_name,substance,dilution_factor,solvent,concentration
# nanodrop_col_name,substance,dilution_factor,solvent,concentration
# 0,methoxychalcone,1,ethanol,0.00005
# 1,methoxychalcone,1,ethanol,0.0001
# 2,methoxychalcone,1,ethanol,0.0002
# 3,methoxychalcone,1,ethanol,0.0003
# 5,anisaldehyde,1,ethanol,0.00005
# 6,anisaldehyde,1,ethanol,0.0001
# 7,anisaldehyde,1,ethanol,0.0002
# 8,anisaldehyde,1,ethanol,0.0003
# 9,anisaldehyde,1,ethanol,0.0004
# 10,acetophenone,1,ethanol,0.00005
# 11,acetophenone,1,ethanol,0.0001
# 12,acetophenone,1,ethanol,0.0002
# 13,acetophenone,1,ethanol,0.0003
# 14,acetophenone,1,ethanol,0.0004
# 15,solvent_ethanol,1,ethanol,0

# Load the spectral data into dataframe
calibration_source_filename = 'data/calibration_of_all_components'
nanodrop_df = pd.read_csv(calibration_source_filename + '.csv')
# rename first column to "wavelength" and make it float type
nanodrop_df = nanodrop_df.rename(columns={nanodrop_df.columns[0]: "wavelength"})
nanodrop_df["wavelength"] = nanodrop_df["wavelength"].astype(float)

# make an array of wavelengths
wavelengths = nanodrop_df['wavelength'].to_numpy()
wavelength_indices = np.arange(wavelengths.shape[0])

# Load the last file to pandas dataframe and change type of the `nanodrop_col_name` column to string.
all_calibrants_df = pd.read_csv(calibration_source_filename + '.txt')
all_calibrants_df['nanodrop_col_name'] = all_calibrants_df['nanodrop_col_name'].astype(str)
concentration_column_name = 'concentration'

# identifu the background spectrum -- it is one with zero concentration
col_name_with_background = all_calibrants_df.loc[all_calibrants_df[concentration_column_name] == 0].iloc[0][
    'nanodrop_col_name']
absorbances = nanodrop_df[col_name_with_background].to_numpy()
bkg_spectrum = np.array([wavelengths, absorbances]).T

# subtract the absorbances of background spectrum from all the rows of nanodrop_df except the furst one
for col in nanodrop_df.columns[1:]:
    nanodrop_df[col] = nanodrop_df[col] - nanodrop_df[col_name_with_background]

#### now let's do the calibration for methoxychalcone
calibrant_shortname = 'methoxychalcone'
ref_concentration = 0.0002 # reference concentration for methoxychalcone
one_calibrant_df = all_calibrants_df[all_calibrants_df['substance'] == calibrant_shortname]

# identify the spectrum with reference concentrations: column of nanodrop_df with the reference concentration in the
# concentration column
target_concentration = ref_concentration
df_row_with_target_concentration = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == target_concentration].iloc[0]
ref_spectrum = nanodrop_df[df_row_with_target_concentration['nanodrop_col_name']].to_numpy()
# subtract the mean of the last 100 points from the reference spectrum from the reference spectrum
# This makes the plots less confusing, but does not affect the fitting, because vertical offset is one of free parameters anyway,
ref_spectrum -= np.mean(ref_spectrum[-100:])

# make a linear interpolator for the reference spectrum
# this is a function that can be evaluated at any wavelength index
reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')

# sort the concentrations in ascending order
concentrations = sorted([0] + one_calibrant_df[concentration_column_name].to_list())

coeffs = []
spectra = []

fig, axes = plt.subplots(1, 5, figsize=(15, 4))
for concentration in concentrations:
    if concentration == 0:
        coeffs.append(0)
        continue

    df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]

    target_concentration = concentration
    df_row_with_target_concentration = \
    one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == target_concentration].iloc[0]
    target_spectrum = nanodrop_df[df_row_with_target_concentration['nanodrop_col_name']].to_numpy()
    target_spectrum -= np.mean(target_spectrum[-100:])
    spectra.append(np.copy(target_spectrum))

    def model_function(xs, a, b):
        return a * reference_interpolator(xs) + b

    initial_guess_of_scaling_coefficient = concentration / ref_concentration
    initial_guess_of_offset = 0
    # the scaling coefficient is bound by zero from below; other bounds are not set
    popt, pcov = curve_fit(model_function, wavelength_indices, target_spectrum,
                           p0=(initial_guess_of_scaling_coefficient, initial_guess_of_offset),
                           bounds=([-1e-10, -np.inf], [np.inf, np.inf]))

    perr = np.sqrt(np.diag(pcov))
    slope = popt[0]
    slope_error = perr[0]
    coeffs.append(slope)

    ax = axes[concentrations.index(concentration)-1]
    ax.plot(wavelengths, target_spectrum, label=f'Spectrum at\nconcentration {concentration}', color='C0', alpha=0.5)
    ax.plot(wavelengths, model_function(wavelength_indices, *popt), color='r',
            label='Scaled reference\nspectrum,\noffset-corrected', alpha=0.5)
    ax.set_ylim(-0.03, np.max((model_function(wavelength_indices, *popt))) * 1.4)
    # ax.set_title(
    #     f"Scaling the reference spectrum to match the spectrum at concentration {concentration}")
    ax.legend()
    ax.set_title(f'Calibration\nfor {calibrant_shortname}')
    ax.set_xlabel('Wavelength, nm')
    ax.set_ylabel('Absorbance')

plt.tight_layout()
plt.show()

fig3 = plt.figure(figsize=(8, 6))
plt.plot(coeffs, concentrations, 'o-')
plt.xlabel('Best-fit scale coefficient')
plt.ylabel('Concentration, mol/liter')
plt.title(f'Calibration for {calibrant_shortname}')
plt.show()

# we are going to make a linear interpolator for the coefficients, so that we can use it later to calculate
# the concentration for a given value of the scaling coefficient

# For later use in spectral unmixng, let's save the interpolators into a list of dictionaries, one dictionary per substance (calibrant)
calibrants = list()
calibrants.append({'reference_interpolator': interpolate.interp1d(wavelength_indices, ref_spectrum,
                                                                   fill_value='extrapolate'),
                   'coeff_to_concentration_interpolator': interpolate.interp1d(coeffs, concentrations,
                                                                           fill_value='extrapolate')})

# now let's repeat the same steps of calibaration for anisaldehyde and acetophenone

##### CALIBRATION FOR ANYSALDEHYDE
calibrant_shortname = 'anisaldehyde'
ref_concentration = 0.0002 # reference concentration for methoxychalcone
one_calibrant_df = all_calibrants_df[all_calibrants_df['substance'] == calibrant_shortname]

# identify the spectrum with reference concentrations: column of nanodrop_df with the reference concentration in the
# concentration column
target_concentration = ref_concentration
df_row_with_target_concentration = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == target_concentration].iloc[0]
ref_spectrum = nanodrop_df[df_row_with_target_concentration['nanodrop_col_name']].to_numpy()
# subtract the mean of the last 100 points from the reference spectrum from the reference spectrum
# This makes the plots less confusing, but does not affect the fitting, because vertical offset is one of free parameters anyway,
ref_spectrum -= np.mean(ref_spectrum[-100:])

# make a linear interpolator for the reference spectrum
# this is a function that can be evaluated at any wavelength index
reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')

# sort the concentrations in ascending order
concentrations = sorted([0] + one_calibrant_df[concentration_column_name].to_list())

coeffs = []
spectra = []

fig, axes = plt.subplots(1, 5, figsize=(15, 4))
for concentration in concentrations:
    if concentration == 0:
        coeffs.append(0)
        continue

    df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]

    target_concentration = concentration
    df_row_with_target_concentration = \
    one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == target_concentration].iloc[0]
    target_spectrum = nanodrop_df[df_row_with_target_concentration['nanodrop_col_name']].to_numpy()
    target_spectrum -= np.mean(target_spectrum[-100:])
    spectra.append(np.copy(target_spectrum))

    def model_function(xs, a, b):
        return a * reference_interpolator(xs) + b

    initial_guess_of_scaling_coefficient = concentration / ref_concentration
    initial_guess_of_offset = 0
    # the scaling coefficient is bound by zero from below; other bounds are not set
    popt, pcov = curve_fit(model_function, wavelength_indices, target_spectrum,
                           p0=(initial_guess_of_scaling_coefficient, initial_guess_of_offset),
                           bounds=([-1e-10, -np.inf], [np.inf, np.inf]))

    perr = np.sqrt(np.diag(pcov))
    slope = popt[0]
    slope_error = perr[0]
    coeffs.append(slope)

    ax = axes[concentrations.index(concentration)-1]
    ax.plot(wavelengths, target_spectrum, label=f'Spectrum at\nconcentration {concentration}', color='C0', alpha=0.5)
    ax.plot(wavelengths, model_function(wavelength_indices, *popt), color='r',
            label='Scaled reference\nspectrum,\noffset-corrected', alpha=0.5)
    ax.set_ylim(-0.03, np.max((model_function(wavelength_indices, *popt))) * 1.4)
    # ax.set_title(
    #     f"Scaling the reference spectrum to match the spectrum at concentration {concentration}")
    ax.legend()
    ax.set_title(f'Calibration\nfor {calibrant_shortname}')
    ax.set_xlabel('Wavelength, nm')
    ax.set_ylabel('Absorbance')

plt.tight_layout()
plt.show()

fig3 = plt.figure(figsize=(8, 6))
plt.plot(coeffs, concentrations, 'o-')
plt.xlabel('Best-fit scale coefficient')
plt.ylabel('Concentration, mol/liter')
plt.title(f'Calibration for {calibrant_shortname}')
plt.show()

calibrants.append({'reference_interpolator': interpolate.interp1d(wavelength_indices, ref_spectrum,
                                                                   fill_value='extrapolate'),
                   'coeff_to_concentration_interpolator': interpolate.interp1d(coeffs, concentrations,
                                                                           fill_value='extrapolate')})


##### CALIBRATION FOR ANYSALDEHYDE
calibrant_shortname = 'acetophenone'
ref_concentration = 0.0003 # reference concentration for methoxychalcone
one_calibrant_df = all_calibrants_df[all_calibrants_df['substance'] == calibrant_shortname]

# identify the spectrum with reference concentrations: column of nanodrop_df with the reference concentration in the
# concentration column
target_concentration = ref_concentration
df_row_with_target_concentration = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == target_concentration].iloc[0]
ref_spectrum = nanodrop_df[df_row_with_target_concentration['nanodrop_col_name']].to_numpy()
# subtract the mean of the last 100 points from the reference spectrum from the reference spectrum
# This makes the plots less confusing, but does not affect the fitting, because vertical offset is one of free parameters anyway,
ref_spectrum -= np.mean(ref_spectrum[-100:])

# make a linear interpolator for the reference spectrum
# this is a function that can be evaluated at any wavelength index
reference_interpolator = interpolate.interp1d(wavelength_indices, ref_spectrum, fill_value='extrapolate')

# sort the concentrations in ascending order
concentrations = sorted([0] + one_calibrant_df[concentration_column_name].to_list())

coeffs = []
spectra = []

fig, axes = plt.subplots(1, 5, figsize=(15, 4))
for concentration in concentrations:
    if concentration == 0:
        coeffs.append(0)
        continue

    df_row_here = one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == concentration].iloc[0]

    target_concentration = concentration
    df_row_with_target_concentration = \
    one_calibrant_df.loc[one_calibrant_df[concentration_column_name] == target_concentration].iloc[0]
    target_spectrum = nanodrop_df[df_row_with_target_concentration['nanodrop_col_name']].to_numpy()
    target_spectrum -= np.mean(target_spectrum[-100:])
    spectra.append(np.copy(target_spectrum))

    def model_function(xs, a, b):
        return a * reference_interpolator(xs) + b

    initial_guess_of_scaling_coefficient = concentration / ref_concentration
    initial_guess_of_offset = 0
    # the scaling coefficient is bound by zero from below; other bounds are not set
    popt, pcov = curve_fit(model_function, wavelength_indices, target_spectrum,
                           p0=(initial_guess_of_scaling_coefficient, initial_guess_of_offset),
                           bounds=([-1e-10, -np.inf], [np.inf, np.inf]))

    perr = np.sqrt(np.diag(pcov))
    slope = popt[0]
    slope_error = perr[0]
    coeffs.append(slope)

    ax = axes[concentrations.index(concentration)-1]
    ax.plot(wavelengths, target_spectrum, label=f'Spectrum at\nconcentration {concentration}', color='C0', alpha=0.5)
    ax.plot(wavelengths, model_function(wavelength_indices, *popt), color='r',
            label='Scaled reference\nspectrum,\noffset-corrected', alpha=0.5)
    ax.set_ylim(-0.03, np.max((model_function(wavelength_indices, *popt))) * 1.4)
    # ax.set_title(
    #     f"Scaling the reference spectrum to match the spectrum at concentration {concentration}")
    ax.legend()
    ax.set_title(f'Calibration\nfor {calibrant_shortname}')
    ax.set_xlabel('Wavelength, nm')
    ax.set_ylabel('Absorbance')

plt.tight_layout()
plt.show()

fig3 = plt.figure(figsize=(8, 6))
plt.plot(coeffs, concentrations, 'o-')
plt.xlabel('Best-fit scale coefficient')
plt.ylabel('Concentration, mol/liter')
plt.title(f'Calibration for {calibrant_shortname}')
plt.show()

calibrants.append({'reference_interpolator': interpolate.interp1d(wavelength_indices, ref_spectrum,
                                                                   fill_value='extrapolate'),
                   'coeff_to_concentration_interpolator': interpolate.interp1d(coeffs, concentrations,
                                                                           fill_value='extrapolate')})


######### SPECTRAL UNMIXING

ignore_pca_bkg = False
use_line = True
default_limit_of_magnitude_of_weights_for_background_pca_components = np.inf

calibrant_shortnames = ['methoxychalcone', 'anisaldehyde', 'acetophenone']

crude_mixtures_df = pd.read_csv('data/spectra_of_crude_mixtures.csv')
# rename first column to "wavelength" and make it float type
crude_mixtures_df = crude_mixtures_df.rename(columns={crude_mixtures_df.columns[0]: "wavelength"})
crude_mixtures_df["wavelength"] = crude_mixtures_df["wavelength"].astype(float)

sample_id = '0'
# load target spectrum
target_spectrum_raw = crude_mixtures_df[sample_id].to_numpy()
# subtract the background spectrum previously determined
target_spectrum = target_spectrum_raw - bkg_spectrum[:, 1]

# Load the PCA components of the background spectrum, which were previously calculated and saved.
background_model_folder = 'data/ethanol_background_model/'
background_interpolators = [interpolate.interp1d(wavelength_indices,
                                                 np.load(background_model_folder + f'component_{i}.npy'),
                                                 fill_value='extrapolate')
                            for i in range(2)]


# Define the model function to be used by curve_fit
def model_function(*args):
    """
    Spectral model function to be used by scipy.optimize.curve_fit in single-spectrum unmixing.

    Models the measured absorption spectrum as a linear combination of calibrant reference
    spectra plus instrumental baseline and weighted PCA components of background. Since it takes
    wavelength indices as the first argument, the mathematical formulation of this implementation is

    .. math::

        A(i) = \sum_j a_j R_j(i) + b_0 + b_1 \cdot i + \sum_k w_k B_k(i)

    where:

    - :math:`A(i)` is predicted absorbance at wavelength index i
    - :math:`a_j` is scaling coefficient for calibrant j. This coefficient is proportional to concentration.
    - :math:`R_j(i)` is reference spectrum of calibrant j at wavelength index i
    - :math:`b_0` is constant baseline offset (instrumental drift)
    - :math:`b_1` is linear baseline slope (wavelength-dependent drift)
    - :math:`w_k` is weight for PCA background component k
    - :math:`B_k(i)` is PCA background component k at wavelength index i

    The fitted scaling coefficients :math:`a_j` are later converted to concentrations using calibration curves.
    Parameters follow the requirements of scipy.optimize.curve_fit.


    Parameters
    ----------
    args[0] : numpy.ndarray
        Wavelength indices where spectrum should be evaluated
    args[1:-4] : float
        As many scaling coefficients as needed - one for each calibrant's reference spectrum. These will be fitted.
    args[-4:] : float
        [baseline_offset, linear_slope, pca1_weight, pca2_weight]. These will be fitted.

    Returns
    -------
    numpy.ndarray
        Predicted absorption spectrum at the input wavelength indices

    Notes
    -----
    This function accesses calibrants and background_interpolators from the enclosing scope.
    It is designed specifically as the model function to be used by curve_fit for optimization.
    """
    # Extract wavelength indices, at which the model spectrum should be evaluated
    wavelength_indices = args[0]

    # Extract fitted baseline and background parameters (last 4 fitted parameters)
    baseline_offset, linear_slope, pca1_weight, pca2_weight = args[-4:]

    # Extract fitted scaling coefficients for each calibrant (middle fitted parameters)
    calibrant_coefficients = args[1:-4]

    # Build predicted spectrum as linear combination of components
    predicted_spectrum = (
        # Core Beer-Lambert law: sum of scaled reference spectra
            sum([calibrant_coefficients[i] * calibrants[i]['reference_interpolator'](wavelength_indices)
                 for i in range(len(calibrant_coefficients))]) +

            # Instrumental baseline correction: constant offset + linear drift
            baseline_offset + linear_slope * wavelength_indices +

            # Systematic background variations captured by PCA components
            pca1_weight * background_interpolators[0](wavelength_indices) +
            pca2_weight * background_interpolators[1](wavelength_indices)
    )

    return predicted_spectrum


# define the initial guess for the parameters
p0 = tuple([0.0005]*len(calibrant_shortnames) + [0] * 4)

# if use_line is False, will only fit the baseline offset. If True, will also add a baseline component linearly increasing with
# wavelength index. Slope is also a free parameter that is fitted.
if use_line:
    linebounds = [-np.inf, np.inf]
else:
    linebounds = [-1e-15, 1e-15]

# IF ignore_pca_bkg is True, the PCA background components are not fitted at all, and their weights are set to zero.
# In this implementation, this is achieved by setting the bounds for their weights to a very small values around zero.
if ignore_pca_bkg:
    bkg_comp_limit = 1e-12
else:
    bkg_comp_limit = default_limit_of_magnitude_of_weights_for_background_pca_components


# Collect all the bounds for the parameters to be fitted into a form required by the `scipy.optimize.curve_fit` method.
bounds = ([-1e-20] * len(calibrant_shortnames) + [-np.inf, linebounds[0], -1 * bkg_comp_limit, -1 * bkg_comp_limit],
          [np.inf] * len(calibrant_shortnames) + [np.inf, linebounds[1], bkg_comp_limit, bkg_comp_limit])

# perform the curve fitting to find the best values of parameters (returned as `popt` array)
popt, pcov = curve_fit(model_function, wavelength_indices, target_spectrum,
                       p0=p0, bounds=bounds)
perr = np.sqrt(np.diag(pcov))  # errors of the best-fit parameters

# Uses the fitted coefficients to calculate concentrations of the components in the target spectrum
# with the help of previously defined interpolators.
concentrations_here = [calibrants[calibrant_index]['coeff_to_concentration_interpolator'](fitted_coeff)
                       for calibrant_index, fitted_coeff in enumerate(popt[:-4])]

# propagate the uncertainty of the fitted coefficients to the concentrations
concentration_errors = [calibrants[calibrant_index]['coeff_to_concentration_interpolator'](fitted_coeff + perr[calibrant_index])
                          - calibrants[calibrant_index]['coeff_to_concentration_interpolator'](fitted_coeff - perr[calibrant_index])
                          for calibrant_index, fitted_coeff in enumerate(popt[:-4])]

# This is the dilution factor of the crude mixture, which is used to calculate the final concentrations.
dilution_factor_of_the_crude_mixture = 500
# Multiply the concentrations by the dilution factor to get the final concentrations in the crude mixture.
concentrations_here = [conc * dilution_factor_of_the_crude_mixture for conc in concentrations_here]
concentration_errors = [err * dilution_factor_of_the_crude_mixture for err in concentration_errors]


# Print concentration in a nicely formatted way
print("Concentrations of the components in the target spectrum:")
for i, conc in enumerate(concentrations_here):
    print(f"{calibrant_shortnames[i]}: concentration {conc:.6f} mol/L ± {concentration_errors[i]:.6f} mol/L "
          f"(fitted coefficient: {popt[i]:.6f} ± {perr[i]:.6f})")

# Plot the results

fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
ax = ax1
ax.plot(wavelengths, target_spectrum_raw, label='Raw data', color='grey', alpha=0.2)
ax.plot(wavelengths, target_spectrum, label='Data minus bkg.', color='black', alpha=0.5)
ax.plot(wavelengths, model_function(wavelength_indices, *popt), color='r', label='Fit', alpha=0.5)
for calibrant_index in range(len(calibrant_shortnames)):
    cpopt = [x if i == calibrant_index else 0 for i, x in enumerate(popt)]
    ax.plot(wavelengths, model_function(wavelength_indices, *cpopt), label=calibrant_shortnames[calibrant_index], alpha=0.5)
# make a list where only the third from the end item is the same as in popt, while the other ones are zero
if use_line:
    cpopt = [x if i == len(popt) - 3 else 0 for i, x in enumerate(popt)]
    ax.plot(wavelengths, model_function(wavelength_indices, *cpopt), label='Line', alpha=0.5)
if not ignore_pca_bkg:
    cpopt = [x if i == len(popt) - 2 else 0 for i, x in enumerate(popt)]
    ax.plot(wavelengths, model_function(wavelength_indices, *cpopt), label='Bkg. PC1', alpha=0.5)
    cpopt = [x if i == len(popt) - 1 else 0 for i, x in enumerate(popt)]
    ax.plot(wavelengths, model_function(wavelength_indices, *cpopt), label='Bkg. PC2', alpha=0.5)

title_str = f'Concentrations:\n'
for i in range(len(concentrations_here)):
    title_str += f'{np.array(concentrations_here)[i]:.6f} M ({calibrant_shortnames[i]})\n '
fig1.suptitle(title_str[:-2])
ax.set_ylabel('Absorbance')
ax.legend()
# Residuals subplot
ax = ax2
ax.plot(wavelengths, target_spectrum - model_function(wavelength_indices, *popt), color='black', alpha=0.5,
        label='residuals')
ax.legend()
ax.set_xlabel('Wavelength, nm')
ax.set_ylabel('Absorbance')

plt.tight_layout()
plt.show()