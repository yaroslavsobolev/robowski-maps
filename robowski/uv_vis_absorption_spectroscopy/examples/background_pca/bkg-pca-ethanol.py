import importlib
import os
from sklearn.decomposition import PCA as pca
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
process_wellplate_spectra = importlib.import_module("uv_vis_absorption_spectroscopy.process_wellplate_spectra")
sp = process_wellplate_spectra.SpectraProcessor(folder_with_correction_dataset='uv_vis_absorption_spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')
sp.nanodrop_lower_cutoff_of_wavelengths = 220

bkg_model_folder = data_folder + 'BPRF/cross_conamination_and_backgound_test/ethanol_background_model/'

def read_peculiar_file_of_yankai(plate_folder):
    nanodrop_df = pd.read_csv(plate_folder)

    # rename first column to "wavelength" and make it float type
    nanodrop_df = nanodrop_df.rename(columns={nanodrop_df.columns[0]: "wavelength"})

    # remove rows where wavelength is lower than nanodrop_lower_cutoff_of_wavelengths
    nanodrop_df = nanodrop_df[nanodrop_df["wavelength"] >= sp.nanodrop_lower_cutoff_of_wavelengths]

    # remove rows where wavelength is higher than nanodrop_upper_cutoff_of_wavelengths
    nanodrop_df = nanodrop_df[nanodrop_df["wavelength"] <= sp.nanodrop_upper_cutoff_of_wavelengths]

    nanodrop_df["wavelength"] = nanodrop_df["wavelength"].astype(float)

    return nanodrop_df

plate_folder = data_folder + 'BPRF/cross_conamination_and_backgound_test/2024-02-17_21-48-29_2024_02_17_UV-Vis_1.csv'
# plate_folder = data_folder + 'BPRF/2024-02-16-run01/nanodrop_spectra/2024-02-18_17-48-07_UV-Vis_plate74.csv'
nanodrop_df = read_peculiar_file_of_yankai(plate_folder)
# iterate over all columns except the 'wavelegth' column
wavelengths = nanodrop_df["wavelength"]
spectra = []
for column in nanodrop_df.columns[1:]:
    first_number_in_column_name_before_underscore = int(column.split('_')[0])
    if (nanodrop_df[column].max() < 0.1) and (first_number_in_column_name_before_underscore in [2, 5, 8, 11,14, 17]):
        spectra.append(nanodrop_df[column])
    plt.plot(nanodrop_df["wavelength"], nanodrop_df[column], label=column)


# Second file
plate_folder = data_folder + 'BPRF/cross_conamination_and_backgound_test/2024-03-10_00-22-33_UV-Vis_pure_ethanol_54_times.csv'
# plate_folder = data_folder + 'BPRF/2024-02-16-run01/nanodrop_spectra/2024-02-18_17-48-07_UV-Vis_plate74.csv'
nanodrop_df = read_peculiar_file_of_yankai(plate_folder)
# iterate over all columns except the 'wavelegth' column
wavelengths = nanodrop_df["wavelength"]
for column in nanodrop_df.columns[1:]:
    if (nanodrop_df[column].max() < 10.1):
        spectra.append(nanodrop_df[column])
    plt.plot(nanodrop_df["wavelength"], nanodrop_df[column], label=column)


plt.legend()
plt.show()

cut_from = 0

# Principal component analysis of the background spectra
nc = 1
pca1 = pca(n_components=nc)
nfeat1 = spectra
nfeat1 = [x - np.mean(x) for x in nfeat1]
X1 = pca1.fit(nfeat1)
expl_var_1 = X1.explained_variance_ratio_
print(f'Explained variance: {expl_var_1}')
sv = X1.components_.T
for spectrum in spectra:
    plt.plot(wavelengths, spectrum, alpha=0.2, color='grey')
colors = ['C0', 'C1', 'C2']
for i in range(nc):
    plt.plot(wavelengths, sv[:, i], label=f'PC {i+1}', color=colors[i], alpha=0.3)

mean_spectrum = np.mean(nfeat1, axis=0)
plt.plot(wavelengths, mean_spectrum, label='mean', color='black', alpha=1)

# Smooth the components and plot them
for i in range(nc):
    smoothed_component = savgol_filter(sv[:, i], polyorder=2, window_length=31)
    plt.plot(wavelengths, smoothed_component, linestyle='--', label=f'smoothed PC {i+1}', color=colors[i])
    # for_saving = np.zeros_like(spectrum)
    for_saving = smoothed_component
    np.save(
        f'{bkg_model_folder}component_{i}.npy',
        for_saving)
plt.title('PCA of background spectra')
plt.xlabel('Wavelength, nm')
plt.ylabel('Absorbance units,\nspectra are shifted to be around zero')
plt.show()

bkg_spectrum = np.array([wavelengths, mean_spectrum]).T
np.save(f'{bkg_model_folder}bkg_spectrum.npy', bkg_spectrum)

# Save the components
# np.save(data_folder + f'simple-reactions/2023-07-05-run01/microspectrometer_data/background_model/component_{i}.npy',
#         for_saving)
