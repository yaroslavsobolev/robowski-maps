from robowski.settings import *
# plot the covariance matrix obtained from the curve_fit
import numpy as np
import matplotlib.pyplot as plt

# load the covariance_matrix from the file

cov_matrix_filepath = repo_data_path + 'kinetics_models/Ugi_kinetics/ugi_v3_REV_outputs/nopka_v3rrr_diff0.2_3point/v3r_REV_theta_pcov_nopka.npy'
# cov_matrix_filepath = repo_data_path + 'kinetics_models/Ugi_kinetics/ugi_v3_REV_outputs/nopka_v3_diff0p2_3points/v3_REV_theta_pcov_nopka_2025-02-15.npy'

cov_matrix = np.load(cov_matrix_filepath)

# convert covariance matrix to correlation matrix
correlation_matrix = np.zeros_like(cov_matrix)
for i in range(cov_matrix.shape[0]):
    for j in range(cov_matrix.shape[1]):
        correlation_matrix[i, j] = cov_matrix[i, j] / np.sqrt(cov_matrix[i, i] * cov_matrix[j, j])

print(correlation_matrix)

# plot the correlation matrix
# parameter_names = ['$k_1$', '$k_{-1}$', '$k_2$', '$k_{-2}$', '$k_3$', '$k_{-3}$', '$k_4$', '$k_{-4}$', '$k_5$', '$k_6$', '$k_{-6}$', '$k_7$', '$k_8$', '$k_9$', '$k_{m1}$', '$k_{-m1}$', '$k_{m3}$']
parameter_names = ['$k_1$', '$k_{-1}$', '$k_2$', '$K_2$', '$K_3$', '$k_{-3}$', '$k_4$', '$K_{4}$', '$k_5$', '$k_6$', '$k_{-6}$', '$k_7$', '$k_8$', '$k_9$', '$k_{m1}$', '$k_{-m1}$',
                   '$K_{m3}$']

fig = plt.figure(figsize=(7.5, 6))
# use parameter names as x and y ticks
plt.imshow(correlation_matrix, cmap='coolwarm', vmin=-1, vmax=1)
plt.xticks(np.arange(len(parameter_names)), parameter_names, rotation=0)
plt.yticks(np.arange(len(parameter_names)), parameter_names)

# add colorbar with axis label (correlation coefficient)
cbar = plt.colorbar()
cbar.set_label('Correlation coefficient')
plt.tight_layout()
fig.savefig(repo_data_path + 'misc_scripts/figures/correlation_matrix_ugi_reformulated_v2.png', dpi=300)
# fig.savefig(repo_data_path + 'misc_scripts/figures/correlation_matrix_ugi.png', dpi=300)
plt.show()

# Load the perr from 'perr_nopka.npy' file, then multiply the errors by (2*5.8/85)**(1/2)
# perr = np.load(repo_data_path + 'kinetics_models/Ugi_kinetics/ugi_v3_REV_outputs/nopka_v3_diff0p2_3points/v3_REV_theta_perr_nopka_2025-02-15.npy')
perr = np.load(repo_data_path + 'kinetics_models/Ugi_kinetics/ugi_v3_REV_outputs/nopka_v3rrr_diff0.2_3point/v3r_REV_theta_perr_nopka.npy')

perr = perr * (2 * 5.8 / 85) ** (0.5)

# print the errors with 3 significant digits ('g')
# assert equal lengths of parameter_names and perr
assert len(parameter_names) == len(perr)
for i, name in enumerate(parameter_names):
    print(f'{name}: {perr[i]:.3g}')