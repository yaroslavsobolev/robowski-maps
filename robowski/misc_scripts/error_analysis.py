from uncertainties import ufloat
from uncertainties.unumpy import uarray, nominal_values, std_devs

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

def pipetting_error_interpolation_right_after_calib(volume): # return absolute error

    '''This interpolater is based on measurements right after calibration. Considering calibration
    drift with time and environment changes, this tends to underestimate the pipetting error.
    '''

    csv_error = "_pipetting_error_all_solvents.csv"
    df_error = pd.read_csv(csv_error)

    set_vols = df_error['set_vol'].to_numpy()
    errors = df_error['ethanol_error_%'].to_numpy() / 100 # percentage to numerical

    error_interpolator = interp1d(set_vols, errors, kind='cubic', fill_value="extrapolate")

    return abs(error_interpolator(volume) * volume) # return absolute error


def pipetting_error_interpolation(volumes):
    """
    Calculate the absolute error for given volumes.

    Uses the following rules:
      - 2% error for volumes <= 50
      - 1% error for 50 < volumes <= 300
      - 0.5% error for 300 < volumes <= 1000

    Parameters:
    volumes : int, float, list, or numpy array
        The input volume(s) to calculate the error for.

    Returns:
    list or float:
        Absolute errors corresponding to the input volumes.
    """
    # Ensure volumes is an iterable array-like object
    if isinstance(volumes, (int, float)):  # Handle single number
        volumes = [volumes]
    elif isinstance(volumes, (list, np.ndarray)):  # Handle list or array
        volumes = np.array(volumes)
    else:
        raise TypeError("Input must be a number, list, or numpy array.")

    # Calculate errors
    errors = []
    for volume in volumes:
        if volume <= 50:
            err_abs = volume * 2 / 100
        elif 50 < volume <= 300:
            err_abs = volume * 1 / 100
        elif 300 < volume <= 1000:
            err_abs = volume * 0.5 / 100
        else:
            raise ValueError('Volume must be in the range of 0 to 1000.')
        errors.append(err_abs)

    # Return single value if input was a single number
    if len(errors) == 1:
        return errors[0]
    return errors


# Function to format uarray with percentage uncertainty
def format_uarray_with_percentage(uarray):
    formatted_values = []
    percentage_uncertainties = []
    for val in uarray:
        nominal = val.nominal_value
        std_dev = val.std_dev
        if nominal != 0:  # Avoid division by zero
            percentage_uncertainty = (std_dev / nominal) * 100
            formatted_values.append(f"{nominal:.7f} ± {std_dev:.7f} ({percentage_uncertainty:.2f}%)")
            percentage_uncertainties.append(percentage_uncertainty)
        else:
            formatted_values.append(f"{nominal:.7f} ± {std_dev:.7f} (undefined %)")
    print(formatted_values)
    print(f'Max: {np.max(percentage_uncertainties)}')
    return formatted_values

excel = "20241118_Hantzsch_HPLC_copy.xlsx"
df = pd.read_excel(excel, engine='openpyxl')

vol_cols = ['vol#EtOH',
           'vol#methoxybenzaldehyde',
           'vol#ethyl_acetoacetate',
           'vol#ammonium_acetate']
df = df.sort_values(by=vol_cols, ascending=[True, True, True, True])

vols = np.array([df['vol#EtOH'].to_numpy(),
        df['vol#methoxybenzaldehyde'].to_numpy(),
        df['vol#ethyl_acetoacetate'].to_numpy(),
        df['vol#ammonium_acetate'].to_numpy()])
errors = np.array([pipetting_error_interpolation(vol) for vol in vols])

# get errors for total vol.
v0 = uarray(vols[0], errors[0])
v1 = uarray(vols[1], errors[1])
v2 = uarray(vols[2], errors[2])
v3 = uarray(vols[3], errors[3])
df['EtOH_vol_u'] = v0
df['methoxybenzaldehyde_vol_u'] = v1
df['ethyl_acetoacetate_vol_u'] = v2
df['ammonium_acetate_vol_u'] = v3
vol_tols = v0 + v1 + v2 + v3
df['vol_tol_u'] = vol_tols

# vol_tol_nominals = [vol_tol.nominal_value for vol_tol in vol_tols]
# vol_tol_errors_abs = [vol_tol.std_dev for vol_tol in vol_tols]
# vol_tol_errors_rel = [vol_tol.std_dev / vol_tol.nominal_value for vol_tol in vol_tols]
# vol_tol_errors_rel_percentage = vol_tol_errors_rel * 100

# get error for dilution.
# [529.4, 195.7, 1690.7] are the real volumes used for dilution, specified vols [530, 195, 1690].
# The real total volume in a reaction vial is 540.
# So the dil = (529.4+540) / 540 * (1690.7+195.7) / 195.7 = 19
vol_a1 = [530] * len(v0) # vol. added to crude
vol_a1_error = [pipetting_error_interpolation(vol) for vol in vol_a1]
vol_a1_u = uarray(vol_a1, vol_a1_error)
vol_a2 = [195] * len(v0) # vol. transferred from crude vial to new vial
vol_a2_error = [pipetting_error_interpolation(vol) for vol in vol_a2]
vol_a2_u = uarray(vol_a2, vol_a2_error)
vol_a3 = [1690] * len(v0) # vol. added to new vial
vol_a3_error = [(pipetting_error_interpolation(1000)**2+pipetting_error_interpolation(690)**2)**0.5 \
                for vol in vol_a3]
vol_a3_u = uarray(vol_a3, vol_a3_error)
dil_u = (vol_a1_u + vol_tols) / vol_tols * (vol_a2_u + vol_a3_u) / vol_a2_u
df['dil_u'] = dil_u
# dil_nominals = nominal_values(dil_u)
# dil_std_dev = std_devs(dil_u)

# get errors for conc.
# the stock solutions are made by Hamilton syringes and analytical balance
c_stock_error = 0.011 # Here assume the errors of stock solutions are 2%. This is a rough estimation.
c1_stock_u, c2_stock_u, c3_stock_u \
    = (uarray([0.33333]*len(v1), [0.33333 * c_stock_error]*len(v1)),
       uarray([0.33333]*len(v1), [0.33333 * c_stock_error]*len(v1)),
       uarray([1.33333]*len(v1), [1.33333 * c_stock_error]*len(v1)))

c1 = c1_stock_u * v1 / vol_tols # for methoxybenzaldehyde
c2 = c2_stock_u * v2 / vol_tols # for ethyl_acetoacetate
c3 = c3_stock_u * v3 / vol_tols # for ammonium acetate
df['methoxybenzaldehyde_c_u'] = c1
df['ethyl_acetoacetate_c_u'] = c2
df['ammonium_acetate_c_u'] = c3
# c1_error_abs, c2_error_abs, c3_error_abs \
#     = std_devs(c1), std_devs(c2), std_devs(c3)
# c1_error_rel, c2_error_rel, c3_error_rel\
#     = std_devs(c1) / nominal_values(c1), std_devs(c2) / nominal_values(c2), std_devs(c3) / nominal_values(c3)
# c_error_rel = np.array([c1_error_rel, c2_error_rel, c3_error_rel])
# ##This is error for each substrate.
# c_error_rel_percentage = c_error_rel * 100

HE_c_limit_u = np.min([c1 / 1, c2 / 2, c3 / 1], axis=0)
# HE_c_limit_percentage = std_devs(HE_c_limit_u) / nominal_values(HE_c_limit_u) * 100
THP_c_limit_u = np.min([c1 / 2, c2 / 1, c3 / 2], axis=0)
# THP_c_limit_percentage = std_devs(THP_c_limit_u) / nominal_values(THP_c_limit_u) * 100
HE_c_limit_after_dil_u = HE_c_limit_u / dil_u
THP_c_limit_after_dil_u = THP_c_limit_u / dil_u
df['HE_c_limit_after_dil_u'] = HE_c_limit_after_dil_u
df['THP_c_limit_after_dil_u'] = THP_c_limit_after_dil_u


### error for HPLC measurements
HE_c_HPLC = df['HE_c_HPLC_insert'].to_numpy()
THP_c_HPLC = df['THP_c_HPLC_insert'].to_numpy()
# Next, how to assign the error for the conc measured by HPLC?
# Varification of the HPLC calibration by known compound concentration shows that the HPLC error is small, <5%.
# This error is the real HPLC measurement error.
HE_c_real = [0.001491, 0.000746] # this conc is prepared by balance and hamilton syringe. In mole/liter
HE_c_by_HPLC = [0.0015467, 0.00078381] # this conc measured by HPLC for the two solutions.
THP_c_real = [0.000868, 0.0003, 0.0002]
THP_c_by_HPLC = [9.01E-04,3.04E-04,2.00E-04,]

HE_c_HPLC_error_rel = np.max([abs(i-j)/j for i,j in zip(HE_c_by_HPLC, HE_c_real)]) # use the max as HPLC error 5.06%
THP_c_HPLC_error_rel = np.max([abs(i-j)/j for i,j in zip(THP_c_by_HPLC, THP_c_real)]) # use the max as HPLC error 3.80%

print(f'HE_c_HPLC_error_rel: {HE_c_HPLC_error_rel}')
print(f'THP_c_HPLC_error_rel: {THP_c_HPLC_error_rel}')

HE_c_HPLC_u = uarray(HE_c_HPLC, HE_c_HPLC * HE_c_HPLC_error_rel)
THP_c_HPLC_u = uarray(THP_c_HPLC, THP_c_HPLC * THP_c_HPLC_error_rel)

df['HE_c_HPLC_u'] = HE_c_HPLC_u
df['THP_c_HPLC_u'] = THP_c_HPLC_u


# group the df by the vols
df_group = df.groupby(vol_cols)

# Calculate the mean for the HE column
def ufloat_mean(series):
    # Calculate the mean of ufloat values
    return sum(series) / len(series)

# cal. the median.
HE_c_mean_u = df_group['HE_c_HPLC_u'].apply(np.median).reset_index()
HE_c_mean_u.rename(columns={'HE_c_HPLC_u':'HE_c_HPLC_median_u'}, inplace=True)
# if sample size > 2, median will not reduce error. However, if sample size == 2, median is same as mean and will reduce eror.
# Since the sample size is so small, the error is kept not changed, still with 5.06% and 3.80%
HE_c_mean_u['HE_c_HPLC_median_u'] = uarray(nominal_values(HE_c_mean_u['HE_c_HPLC_median_u'].to_numpy()),
                                           nominal_values(HE_c_mean_u['HE_c_HPLC_median_u'].to_numpy()) * HE_c_HPLC_error_rel)

# cal. the median of THP_c_HPLC
THP_c_mean_u = df_group['THP_c_HPLC_u'].apply(np.median).reset_index()
THP_c_mean_u.rename(columns={'THP_c_HPLC_u':'THP_c_HPLC_median_u'}, inplace=True)
THP_c_mean_u['THP_c_HPLC_median_u'] = uarray(nominal_values(THP_c_mean_u['THP_c_HPLC_median_u'].to_numpy()),
                                           nominal_values(THP_c_mean_u['THP_c_HPLC_median_u'].to_numpy()) * THP_c_HPLC_error_rel)

# merge back with the original df without duplicates
df_no_repeat = df.drop_duplicates(subset=vol_cols).merge(HE_c_mean_u, on=vol_cols).merge(THP_c_mean_u, on=vol_cols)

# get the uncertainties for conc. of HE and THP by HPLC
HE_c_HPLC_mean_u = df_no_repeat['HE_c_HPLC_median_u'].to_numpy() #
THP_c_HPLC_mean_u = df_no_repeat['THP_c_HPLC_median_u'].to_numpy()

# remove repeats
df_no_repeat = df.drop_duplicates(vol_cols).copy()
HE_c_limit_after_dil_u = df_no_repeat['HE_c_limit_after_dil_u'].to_numpy()
THP_c_limit_after_dil_u = df_no_repeat['THP_c_limit_after_dil_u'].to_numpy()

HE_yield_u = HE_c_HPLC_mean_u / HE_c_limit_after_dil_u * 100
THP_yield_u = THP_c_HPLC_mean_u / THP_c_limit_after_dil_u * 100

df_no_repeat['HE_yield_u'] = HE_yield_u
df_no_repeat['THP_yield_u'] = THP_yield_u

# show relative error
df_no_repeat['HE_yield_u_relative_error_%'] = std_devs(HE_yield_u) / nominal_values(HE_yield_u) * 100
df_no_repeat['THP_yield_u_relative_error_%'] = std_devs(THP_yield_u) / nominal_values(THP_yield_u) * 100

# The following row is used as example in the SI error analysis section
row = df_no_repeat[(round(df_no_repeat['vol#EtOH'], 1) == 198.9) & (round(df_no_repeat['vol#methoxybenzaldehyde'], 1) == 162.8)]

for i in range(31):
    print(i)
    a1 = HE_c_HPLC_mean_u[i].s / HE_c_HPLC_mean_u[i].n * 100
    a2 = HE_c_limit_after_dil_u[i].s / HE_c_limit_after_dil_u[i].n * 100
    print(a1)
    print(a2)
    print((a1**2+a2**2)**0.5)
    print(df_no_repeat['HE_yield_u_relative_error_%'].to_numpy()[i])
    print('*****************')


# format_uarray_with_percentage(HE_c_HPLC_mean_u)
# format_uarray_with_percentage(HE_c_limit_after_dil_u)
# format_uarray_with_percentage(HE_yield_u)
#
# format_uarray_with_percentage(THP_c_HPLC_mean_u)
# format_uarray_with_percentage(THP_c_limit_after_dil_u)
# format_uarray_with_percentage(THP_yield_u)
#
df_no_repeat.to_csv('with_uncertainties_final.csv')
df_no_repeat.to_pickle('with_uncertainties_final.pkl')
print('Done')