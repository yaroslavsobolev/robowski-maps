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

# now let's do the calibration for methoxychalcone
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

    p0 = (concentration / ref_concentration, 0)
    bounds = ([-1e-10, -np.inf], [np.inf, np.inf])
    popt, pcov = curve_fit(model_function, wavelength_indices, target_spectrum, p0=p0, bounds=bounds)

    perr = np.sqrt(np.diag(pcov))
    slope = popt[0]
    slope_error = perr[0]
    coeffs.append(slope)

    fig1 = plt.figure(1)
    plt.plot(wavelengths, target_spectrum, label=f'Spectrum at concentration {concentration}', color='C0', alpha=0.5)
    plt.plot(wavelengths, model_function(wavelength_indices, *popt), color='r', label='Scaled reference spectrum, offset-corrected', alpha=0.5)
    plt.ylim(-0.03,
             np.max((model_function(wavelength_indices, *popt))) * 1.4)
    plt.title(
        f"Scaling the reference spectrum to match the spectrum at concentration {concentration}")
    plt.legend()
    plt.xlabel('Wavelength, nm')
    plt.ylabel('Absorbance')
    plt.show()

fig3, axarr = plt.subplots(1, 2, figsize=(10, 5))
ax1, ax2 = axarr
ax2.plot(coeffs, concentrations, 'o-')
ax2.set_xlabel('Best-fit scale coefficient')
ax2.set_ylabel('Concentration, mol/liter')