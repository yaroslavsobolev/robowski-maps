import os
import numpy as np
import pandas as pd
import importlib
from matplotlib import pyplot as plt
import logging
# set level to info
logging.basicConfig(level=logging.INFO)

st = importlib.import_module('uv-vis-absorption-spectroscopy.spectraltools')

data = np.loadtxt('misc-scripts/dft_molar-absorptivity-spectra/dft-calculation-results/dft-monomer-radical-cation-uv-60-states.txt', skiprows=4)
maxy = np.max(data[:, 1])
plt.plot(data[:, 0], data[:, 1]/maxy*0.25, label='Radical cation, theory', linewidth=5, alpha=0.5)
# for shift in [-0.12, 0.12]:
#     plt.plot(data[:, 0]*(1+shift), data[:, 1]/maxy*0.25, color='C0', alpha=0.5)
plt.errorbar(x=529, y=0.086, xerr=529*0.12, marker='', color='C0', alpha=0.9, capsize=5, capthick=2)

data = np.loadtxt('misc-scripts/dft_molar-absorptivity-spectra/dft-calculation-results/dft-monomer-radical-uv-60-states.txt', skiprows=4)
maxy = np.max(data[:, 1])
plt.plot(data[:, 0], data[:, 1]/maxy*1, label='Radical, theory')

data = np.loadtxt('misc-scripts/dft_molar-absorptivity-spectra/dft-calculation-results/dft-monomer-cation-uv-60-states.txt', skiprows=4)
maxy = np.max(data[:, 1])
plt.plot(data[:, 0], data[:, 1]/maxy*0.1, label='Cation, theory')

cary_file = data_folder + 'Yaroslav/SN1_mystery_product_pi nk/SN1_mystery_product_five-wells.csv'
# cary_column_name = '10mm_plate83_well25_pink_rep1_1'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# plt.plot(wav, spec, label='Crude, 50C, 23 mM alc., 2.7 mM HBr, concentrated HBr stock')

# cary_column_name = '10mm_plate79_well11_pink_rep1_1'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# plt.plot(wav, spec, label='Crude, 50C, 6.5 mM alc., 1.3 mM HBr')

cary_column_name = '10mm_plate72_well11_pink_rep1_1'
wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
plt.plot(wav, spec, label='Crude, 26C, 6.5 mM alc., 1.3 mM HBr')
plt.plot(wav, spec - spec[100], '--', color='C2', label='Crude, 26C, 6.5 mM alc., 1.3 mM HBr, minus the scatter')


cary_file = data_folder + 'Yaroslav/SN1_mystery_product_pi nk/deuterated_TLC_pink_from_rafal_rep1.csv'
cary_column_name = 'deuterated_TLC_pink_from_rafal_rep1'
wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
plt.plot(wav, (spec - spec[100])*0.8, label='Pink product isolated by TLC, in deuterated DMF')

# This sample is misnamed. In truth, it's the well 0, not well 1 on the plate.
# cary_column_name = '10mm_plate83_well1_white_rep1_1'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# plt.plot(wav, spec, label='Crude, 50C, transparent')

# cary_column_name = '10mm_plate83_well3_white_rep1_1'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# plt.plot(wav, spec, label='Crude, 50C, yellow')

plt.xlim(350, 800)
plt.ylim(-0.01, 0.4)
plt.xlabel('Wavelength, nm')
plt.ylabel('Absorbance')
plt.legend()
plt.show()