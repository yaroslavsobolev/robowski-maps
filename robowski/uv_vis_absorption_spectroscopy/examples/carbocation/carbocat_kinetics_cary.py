from robowski.settings import *
import os
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import logging
# set level to info
logging.basicConfig(level=logging.INFO)

import robowski.uv_vis_absorption_spectroscopy.spectraltools as st


data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
cary_file = data_folder + 'Yaroslav/carbocation_kinetics_cary/mystery_prod_2024-04-19__18-15_start_12min_interval_1.csv'
cary_column_name = 'mystery_prod_run1_18-15start_1_1'
wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
wav_master = np.copy(wav)
specs = []
for i in range(101):
    logging.info(f'Processing column {i+1}')
    cary_column_name = f'mystery_prod_run1_18-15start_2_{i+1}'
    wav, spec = st.read_cary_agilent_csv_spectrum(cary_file, column_name=cary_column_name)
    assert np.allclose(wav, wav_master)
    # spec -= np.mean(spec[:100])
    specs.append(np.copy(spec))

specs = np.array(specs)
# plot all in semilogy scale
# make a figure with two plots one under another, first is normal, second semilogy of the same data
fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
for i, spec in enumerate(specs):
    # color according to colors sampled from viridis colormap
    axs[0].plot(wav, spec, c=plt.cm.viridis(i/100), alpha=0.5)
    axs[1].semilogy(wav, spec, c=plt.cm.viridis(i/100), alpha=0.5)

# add a viridis colorbar, indicating from 0 to 1200 minutes corresponding to 0 and 1 on the colorbar
plt.xlabel('Wavelength, nm')
axs[0].set_ylabel('Absorbance')
axs[1].set_ylabel('Absorbance, log scale')

# add a colorbar to the top and bottom subplots, on the top of the first subplot
cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='viridis'), ax=axs[1], orientation='horizontal')
cbar.set_label('Time, minutes')
# set ticks from 0 to 20, in hours
cbar.set_ticks(np.arange(0, 21)/20)
cbar.set_ticklabels(np.arange(0, 21))
plt.show()