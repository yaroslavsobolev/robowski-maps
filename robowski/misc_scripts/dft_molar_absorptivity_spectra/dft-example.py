from robowski.settings import *
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

experimental_spectrum = pd.read_csv('D:/Docs/Dropbox/robochem/data/nanodrop-spectrophotometer-measurements/versatility_test/Claisen_WaiShing/2023-10-11_16-26-58_UV-Vis_methoxychalcone.csv')
xs = experimental_spectrum['Unnamed: 0']
ys = experimental_spectrum['3'].to_numpy() - experimental_spectrum['4'].to_numpy()
plt.plot(xs, ys, label='experiment')


dft_data = np.loadtxt(repo_data_path + 'misc_scripts/dft_molar_absorptivity_spectra/dft_uv-vis/methoxy-b3lyp-31.txt', skiprows=4, usecols=(0,1))
# compute absorbance for 1 mm light path and 0.0003 mol/L
dft_data[:, 1] *= 0.0003 / 10
plt.plot(dft_data[:, 0], dft_data[:, 1], label='theory, B3LYP', alpha=0.5)

# dft_data = np.loadtxt(repo_data_path + 'misc_scripts/dft_molar_absorptivity_spectra/dft_uv-vis/methoxy-CAMb3lyp-311_t_cq3.txt', skiprows=4, usecols=(0,1))
# # compute absorbance for 1 mm light path and 0.0003 mol/L
# dft_data[:, 1] *= 0.0003 / 10
# plt.plot(dft_data[:, 0], dft_data[:, 1], label='theory, CAM-B3LYP', alpha=0.5)


dft_data = np.loadtxt(repo_data_path + 'misc_scripts/dft_molar_absorptivity_spectra/dft_uv-vis/methoxy_15states_cam-b3lyp__6-311plusg_d-p_.txt', skiprows=4, usecols=(0,1))
# compute absorbance for 1 mm light path and 0.0003 mol/L
dft_data[:, 1] *= 0.0003 / 10
plt.plot(dft_data[:, 0], dft_data[:, 1], label='theory, CAM-B3LYP_15states', alpha=0.5)


plt.legend()
plt.title('Methoxychalcone')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorbance')
plt.show()