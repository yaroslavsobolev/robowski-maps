import os
import numpy as np
from matplotlib import pyplot as plt

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

# colors = ['C0', 'C1']
# labels = ['batch 1', 'batch 2']
# for i, experiment_name in enumerate(['simple-reactions/2023-11-28-run01/', 'simple-reactions/2023-12-11-run01/']):
#     xs = np.loadtxt(f'{data_folder}{experiment_name}results/kinetics/1000_over_t.txt')
#     ys = np.loadtxt(f'{data_folder}{experiment_name}results/kinetics/logK.txt')
#     plt.scatter(xs, ys, color=colors[i], label=labels[i])
#     # fit line and get slope and intersept
#     z = np.polyfit(xs[:-2], ys[:-2], 1)
#     f = np.poly1d(z)
#     print(f'Intercept {f[0]}, slope {f[1]}')
#     R_gas = 8.31446261815324 # J/(K*mol)
#     print(f'Delta S {f[0] * R_gas} J/(K*mol), delta H {-1 * f[1] * R_gas} kJ/mol')
#     plt.plot(xs, f(xs), '--')
#

colors = ['black']
labels = ['batch 2']
fig = plt.figure(figsize=(4, 3.9))
for i, experiment_name in enumerate(['simple-reactions/2023-12-11-run01/']):
    xs = np.loadtxt(f'{data_folder}{experiment_name}results/kinetics/1000_over_t.txt')
    ys = -1 * np.loadtxt(f'{data_folder}{experiment_name}results/kinetics/logK.txt')
    plt.scatter(xs, ys, color=colors[i], label=labels[i], alpha=0.5)
    # fit line and get slope and intersept
    z = np.polyfit(xs[:-2], ys[:-2], 1)
    f = np.poly1d(z)
    R_gas = 8.31446261815324 # J/(K*mol)
    plt.plot(xs, f(xs), '--', color=colors[i])

    def func(x, a, b):
        return a + b * x

    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(func, xs[:-2], ys[:-2], sigma=0.04*np.ones(len(xs[:-2])), absolute_sigma=True)
    p_sigma = np.sqrt(np.diag(pcov))
    print(f'Delta S {f[0] * R_gas} +- {p_sigma[0] * R_gas} J/(K*mol), delta H {-1 * f[1] * R_gas} +- {-1 * p_sigma[1] * R_gas} kJ/mol')


plt.xlabel('1000/T, K$^{-1}$')
plt.ylabel('ln K')
plt.xlim(3.01, 3.47)
plt.tight_layout()
fig.savefig('misc_scripts/figures/kinetics_of_simpleSN1_2023-12_hoff.png', dpi=300)
# plt.legend()

plt.show()