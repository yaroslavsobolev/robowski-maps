from robowski.settings import *
import os
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import logging
# set level to info
logging.basicConfig(level=logging.INFO)

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'


import robowski.uv_vis_absorption_spectroscopy.spectraltools as st

# make a figure with 6 subplots with sharex = true
fig, axs = plt.subplots(4, 1, figsize=(5, 7), sharex=True)

for ax in axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    # remove y axis ticklabels
    ax.set_yticklabels([])
    # remove ticks
    ax.tick_params(axis='y', which='both', length=0)
    ax.set_ylabel('Absorbance')

ax = axs[0]
ccat_spectrum = np.load(data_folder + 'Yaroslav/mystery_prod/extracted_pink_spectrum.npy')
# remove points below 225 nm
ccat_spectrum = ccat_spectrum[ccat_spectrum[:, 0] > 425]
ax.plot(ccat_spectrum[:, 0], ccat_spectrum[:, 1], label='Unknown product', color='deeppink')
ax.fill_between(x=ccat_spectrum[:, 0], y1=0, y2=ccat_spectrum[:, 1], color='deeppink', alpha=0.4)
ax.set_ylim(0, 0.3)

# ax = axs[1]
# data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
# cary_file = data_folder + 'Yaroslav/r4eal_carbocat/synthetic_carbocat_dcm.csv'
# cary_column_name = 'carbocat_dcm_1uL_in_2500uL_rep1_3'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# ax.plot(wav, spec, label='Carbocation, experiment', color='grey', linewidth=3)
# ax.fill_between(x=wav, y1=0, y2=spec, color='grey', alpha=0.4)
# ax.set_ylim(0, 0.2)
#
# ax = axs[2]
# color = 'C0'
# litdata = np.loadtxt(repo_data_path + 'misc_scripts/literature_spectra/carbocation_Yuichi_Nishimae_2004.txt', delimiter='\t')
# ax.plot(litdata[:, 0], litdata[:, 1], color=color, linewidth=2, alpha=0.9, label='Nishimae et al. (2004)')
# ax.fill_between(x=litdata[:, 0], y1=0, y2=litdata[:, 1], color=color, alpha=0.4)
#
# ax = axs[3]
# color = 'C2'
# litdata = np.loadtxt(repo_data_path + 'misc_scripts/literature_spectra/ammer-sailer-riedle-2012-E12.txt', delimiter='\t')
# ax.plot(litdata[:, 0], litdata[:, 1], color=color, linewidth=2, alpha=0.9, label='Ammer et al. (2012), E3+')
# ax.fill_between(x=litdata[:, 0], y1=0, y2=litdata[:, 1], color=color, alpha=0.4)
#
# ax = axs[4]
# color = 'C3'
# litdata = np.loadtxt(repo_data_path + 'misc_scripts/literature_spectra/ammer-sailer-riedle-2012-E3.txt', delimiter='\t')
# ax.plot(litdata[:, 0], litdata[:, 1], color=color, linewidth=2, alpha=0.9, label='Ammer et al. (2012), E12+')
# ax.fill_between(x=litdata[:, 0], y1=0, y2=litdata[:, 1], color=color, alpha=0.4)


ax = axs[3]
# data = np.loadtxt(repo_data_path + 'misc_scripts/dft_molar-absorptivity-spectra/dft-calculation-results/dft-monomer-radical-cation-uv-60-states.txt', skiprows=4)
nms = np.load(repo_data_path + 'misc_scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_nms.npy')
spectrum = np.load(repo_data_path + 'misc_scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_mean.npy')
perc1 = np.load(repo_data_path + 'misc_scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_perc1.npy')
perc2 = np.load(repo_data_path + 'misc_scripts/figures_for_articles/dft-uv-vis/dimer_OHplus_conf2_1054kjm_uv_perc2.npy')

ax.plot(nms, spectrum, label='Best match, theory', linewidth=2, alpha=0.5)
ax.fill_between(x=nms, y1=0, y2=spectrum, color='C0', alpha=0.4)
ax.set_ylim(0, 0.26)

# for shift in [-0.12, 0.12]:
#     plt.plot(data[:, 0]*(1+shift), data[:, 1]/maxy*0.25, color='C0', alpha=0.5)
# plt.errorbar(x=529, y=0.086, xerr=529*0.12, marker='', color='C0', alpha=0.9, capsize=5, capthick=2)

ax = axs[1]
color='C1'
data = np.loadtxt(repo_data_path + 'misc_scripts/dft_molar-absorptivity-spectra/dft-calculation-results/dft-monomer-radical-uv-60-states.txt', skiprows=4)
maxy = np.max(data[:, 1])
ax.plot(data[:, 0], data[:, 1], color=color, label='Radical, theory')
ax.fill_between(x=data[:, 0], y1=0, y2=data[:, 1], color=color, alpha=0.4)
ax.set_ylim(0, 8500)

ax = axs[2]
color='C2'
data = np.loadtxt(repo_data_path + 'misc_scripts/dft_molar-absorptivity-spectra/dft-calculation-results/dft-monomer-cation-uv-60-states.txt', skiprows=4)
maxy = np.max(data[:, 1])
ax.plot(data[:, 0], data[:, 1], color=color, label='Cation, theory')
ax.fill_between(x=data[:, 0], y1=0, y2=data[:, 1], color=color, alpha=0.4)
ax.set_ylim(0, 22000)


#
# data = np.loadtxt(repo_data_path + 'misc_scripts/dft_molar-absorptivity-spectra/dft-calculation-results/dft-monomer-cation-uv-60-states.txt', skiprows=4)
# maxy = np.max(data[:, 1])
# plt.plot(data[:, 0], data[:, 1]/maxy*0.1, label='Cation, theory')

# cary_file = data_folder + 'Yaroslav/SN1_mystery_product_pi nk/SN1_mystery_product_five-wells.csv'
# cary_column_name = '10mm_plate83_well25_pink_rep1_1'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# plt.plot(wav, spec, label='Crude, 50C, 23 mM alc., 2.7 mM HBr, concentrated HBr stock')

# cary_column_name = '10mm_plate79_well11_pink_rep1_1'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# plt.plot(wav, spec, label='Crude, 50C, 6.5 mM alc., 1.3 mM HBr')

# cary_column_name = '10mm_plate72_well11_pink_rep1_1'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# plt.plot(wav, spec, label='Crude, 26C, 6.5 mM alc., 1.3 mM HBr')
# plt.plot(wav, spec - spec[100], '--', color='C2', label='Crude, 26C, 6.5 mM alc., 1.3 mM HBr, minus the scatter')


# cary_file = data_folder + 'Yaroslav/SN1_mystery_product_pi nk/deuterated_TLC_pink_from_rafal_rep1.csv'
# cary_column_name = 'deuterated_TLC_pink_from_rafal_rep1'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# plt.plot(wav, (spec - spec[100])*0.8, label='Pink product isolated by TLC, in deuterated DMF')

# This sample is misnamed. In truth, it's the well 0, not well 1 on the plate.
# cary_column_name = '10mm_plate83_well1_white_rep1_1'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# plt.plot(wav, spec, label='Crude, 50C, transparent')

# cary_column_name = '10mm_plate83_well3_white_rep1_1'
# wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
# plt.plot(wav, spec, label='Crude, 50C, yellow')

plt.xlim(365, 800)
# plt.ylim(-0.01, 0.4)
plt.xlabel('Wavelength, nm')
# plt.legend()
plt.tight_layout()
plt.show()