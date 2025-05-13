import matplotlib.pyplot as plt

from calibration_illustration import *

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

experiment_name = f'nanodrop-spectrophotometer-measurements/versatility_test/Claisen_WaiShing/'

f = plt.figure(figsize=(4.9, 4.5 * 2/3 * 1.1), dpi=300)

construct_calibrant(
    cut_from=5,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='2023-10-11_16-26-58_UV-Vis_methoxychalcone',
    calibrant_shortnames=['methoxychalcone'],
    calibrant_shortnames_labels=['Methoxychalcone'],
    calibrant_colors=['C0'],
    ref_concentrations=[0.0002],
    max_concentrations=[0.0006],
    experiment_name=experiment_name,
)

construct_calibrant(
    cut_from=5,
    lower_limit_of_absorbance=0.007,
    concentration_column_name='concentration',
    do_plot=True,
    calibration_source_filename='2023-10-11_19-50-36_UV-Vis_anisaldehyde',
    calibrant_shortnames=['anisaldehyde', 'acetophenone'],
    calibrant_shortnames_labels=['Anisaldehyde', 'Acetophenone'],
    calibrant_colors=['C1', 'C2'],
    ref_concentrations=[0.0002, 0.0003],
    max_concentrations=[0.001, 0.001],
    experiment_name=experiment_name,
)

simpleaxis(plt.gca())

# plt.legend(loc='upper right')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorbance')
plt.xlim(220, 420)
plt.tight_layout()
f.savefig('misc_scripts/figures/calibration-illustration.png', dpi=300)
plt.show()