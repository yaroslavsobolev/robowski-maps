import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import corner
# import matplotlib as mpl
# mpl.rcParams["axes.labelpad"] = 20

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

# samples = np.load(f'{data_folder}simple-reactions/2023-09-14-run01/results/kinetics/bootstrapped_keq_fits_size434_2024-12-05.npy')
# samples[:, 5] =  -1 * samples[:, 5]
# samples[:, 6] = 1e5 * samples[:, 6]
# samples[:, 9] = 100 * samples[:, 9]
# labels = ['$\\tilde{\Delta S}_{pKa, COH_2^{+}}$', '$\\tilde{\Delta H}_{pKa, COH_2^{+}}$', '$\Delta S_{net}$', '$\Delta H_{net}$', '$\\tilde{\Delta S}_{dehyd}$/R', '$\\tilde{\Delta H}_{dehyd}/R$', '$10^{5} \kappa_{-}$', '$\Delta S^{\ddag}_{-}$', '$\Delta H^{\ddag}_{-}$', '$10^{2} B_{HBr}$']
# fig = corner.corner(samples, labels=labels, show_titles=True, hist_bin_factor=1, range=[0.95]*10, plot_contours=False, figsize=(12.5, 12.5), bins=14, labelpad=0.06)
# fig.set_size_inches(14.5, 14.5)
# fig.savefig('misc_scripts/figures/bootstrapped_keq_fits_E1.png', dpi=300)
# plt.show()

samples = np.load(f'{data_folder}simple-reactions/2023-09-14-run01/results/kinetics/bootstrapped_keq_fits_size434_2024-12-05.npy')
# reversing signs of enthalpies to match the definitions in the research article
samples[:, 5] =  -1 * samples[:, 5]
samples[:, 1] =  -1 * samples[:, 1]
samples[:, 8] =  -1 * samples[:, 8]

samples[:, 6] = 1e5 * samples[:, 6]
samples[:, 9] = 100 * samples[:, 9]

# for each parameter, print the median, plus and minus 1 sigma (16th and 84th percentiles)
for i in range(samples.shape[1]):
    median = np.median(samples[:, i])
    lower = np.percentile(samples[:, i], 16)
    upper = np.percentile(samples[:, i], 84)
    print(f'{median}^+{upper - median:.3g} _ -{median - lower:.3g}')

truths = np.median(samples, axis=0)
labels = ['$\\tilde{\Delta S}_{pKa}$', '$\\tilde{\Delta H}_{pKa}/R$', '$\Delta S_{net}$', '$\Delta H_{net}$', '$\\tilde{\Delta S}_{dehyd}$/R', '$\\tilde{\Delta H}_{dehyd}/R$', '$10^{5} \kappa_{-} k_{B}/h$', '$\Delta S^{\ddag}_{-}/R$', '$\Delta H^{\ddag}_{-}/R$', '$100 B_{HBr}$']
fig = corner.corner(samples, labels=labels, show_titles=True, hist_bin_factor=1, range=[0.95]*10, plot_contours=False, figsize=(12.5, 12.5), bins=14, labelpad=0.06, truths=truths)
fig.set_size_inches(14.5, 14.5)
fig.savefig('misc_scripts/figures/bootstrapped_keq_fits_E1_2024-12-11b.png', dpi=300)
plt.show()


################## THREE TEMPERATURES ##################
# # samples = np.load(f'{data_folder}simple-reactions/2023-09-14-run01/results/kinetics/bootstrapped_keq_fits_size249_2024-12-14_temps16-21-26.npy')
# samples = np.load(f'{data_folder}simple-reactions/2023-09-14-run01/results/kinetics/bootstrapped_keq_fits_size59_2024-12-17c_temps16-21-26.npy')
# samples_2 = np.load(f'{data_folder}simple-reactions/2023-09-14-run01/results/kinetics/bootstrapped_keq_fits_size305_2024-12-24a_temps16-21-26.npy')
# samples = np.concatenate((samples, samples_2), axis=0)
# print(f'Shape of samples: {samples.shape}')
#
#
# # print medians
# print(list(np.median(samples, axis=0)))
#
# # reversing signs of enthalpies to match the definitions in the research article
# samples[:, 5] =  -1 * samples[:, 5]
# samples[:, 1] =  -1 * samples[:, 1]
# samples[:, 8] =  -1 * samples[:, 8]
#
# samples[:, 6] = 1e5 * samples[:, 6]
# # samples[:, 9] = 100 * samples[:, 9]
#
# # truths = np.load(f'{data_folder}simple-reactions/2023-09-14-run01/results/kinetics/keq_fits_with_errs_tdep_homosc_2024-12-16a_temps16-21-26__from_medians_exact__size249_rep6.npy')[0]
# # truths[5] = -1 * truths[5]
# # truths[6] = 1e5 * truths[6]
# # truths[9] = 100 * truths[9]
# # print(f'truths = {truths}')
#
# truths = np.median(samples, axis=0)
#
# # for each parameter, print the median, plus and minus 1 sigma (16th and 84th percentiles)
# for i in range(samples.shape[1]):
#     median = np.median(samples[:, i])
#     lower = np.percentile(samples[:, i], 16)
#     upper = np.percentile(samples[:, i], 84)
#     print(f'{median}^+{upper - median:.3g} _ -{median - lower:.3g}')
#
# labels = ['$\\tilde{\Delta S}_{pKa}$', '$\\tilde{\Delta H}_{pKa}/R$', '$\Delta S_{net}$', '$\Delta H_{net}$', '$\\tilde{\Delta S}_{dehyd}$/R', '$\\tilde{\Delta H}_{dehyd}/R$', '$10^{5} \kappa_{-} k_{B}/h$', '$\Delta S^{\ddag}_{-}/R$', '$\Delta H^{\ddag}_{-}/R$', '$B_{HBr}$']
# fig = corner.corner(samples, labels=labels, show_titles=True, hist_bin_factor=1, range=[0.95]*10, plot_contours=False, figsize=(12.5, 12.5), bins=14, labelpad=0.06, truths=truths)
# fig.set_size_inches(14.5, 14.5)
# fig.savefig('misc_scripts/figures/bootstrapped_keq_fits_E1_temps16-21-26.png', dpi=300)
# plt.show()