from process_wellplate_spectra import SpectraProcessor
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA as pca
from scipy.signal import savgol_filter

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
sp = SpectraProcessor(folder_with_correction_dataset='uv-vis-absorption-spectroscopy/microspectrometer-calibration/'
                                                     '2022-12-01/interpolator-dataset/')


dilution_factor = 200
experiment_name = 'multicomp-reactions/2023-03-20-run01/'
# ##### This constructs the calibration for the product 'IIO029A' and saves for later. Do not rerun unless you know what you do. #######
# sp.construct_reference_for_calibrant(calibrant_shortname='IIO029A',
#                                      calibration_folder=data_folder + 'multicomp-reactions/2023-01-18-run01/' + 'microspectrometer_data/calibration/',
#                                      ref_concentration=0.00011,
#                                      do_plot=True, do_reference_refinements=True)

# #### This constructs the calibration for the substrate 'ald001' and saves for later. Do not rerun unless you know what you do. #######
# sp.construct_reference_for_calibrant(calibrant_shortname='ald001',
#                                      calibration_folder=data_folder + 'multicomp-reactions/2023-01-18-run01/' + 'microspectrometer_data/calibration/',
#                                      ref_concentration=0.0192096,
#                                      do_plot=True, do_reference_refinements=False)

craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'

def process_or_load_processed(nickname, plate_name, force_recalculation=False):
    filepath = data_folder + experiment_name + 'results/mixing_validation/' + f'{nickname}.npy'

    if force_recalculation or not os.path.isfile(filepath):
        concentrations = sp.concentrations_for_one_plate(experiment_folder=data_folder + experiment_name,
                                                         plate_folder=craic_folder + plate_name + '/',
                                                         calibration_folder=data_folder + 'multicomp-reactions/2023-01-18-run01/' + 'microspectrometer_data/calibration/',
                                                         calibrant_shortnames=['IIO029A', 'ald001'],
                                                         background_model_folder=data_folder + 'multicomp-reactions/2023-03-20-run01/microspectrometer_data/background_model/',
                                                         calibrant_upper_bounds=[np.inf, 1e-10],
                                                         do_plot=False)
        np.save(data_folder + experiment_name + 'results/mixing_validation/' + f'{nickname}.npy',
                concentrations)
    else:
        concentrations = np.load(filepath)

    return concentrations * dilution_factor

def print_stats(input_array, label):
    print(
        f'{label}: Mean: {np.mean(input_array)}, std: {np.std(input_array)}, rel. std: {np.std(input_array) / np.mean(input_array):.1%}')
    return np.std(input_array) / np.mean(input_array)

relative_error = 0.146

diluted_indices = [i + j for i in [9, 27, 45] for j in range(9)]
undiluted_indices = [i + j for i in [0, 18, 36] for j in range(9)]
concentrations_ham_stirred = process_or_load_processed('concentrations_ham_stirred', '2023-03-25_14-33-50__plate0000003__mixing-validation')
concentrations_ham_static = process_or_load_processed('concentrations_ham_static', '2023-03-28_11-29-25__plate0000019__mixing-validation')
concentrations_main = process_or_load_processed('concentrations_main', '2023-03-24_14-29-28__plate0000006__-')

print(f'max well id main: {np.argmax(concentrations_main)}')
print(f'max well id main: {np.argmax(concentrations_ham_stirred)}')
fig = plt.figure(1)
# plt.errorbar(x=concentrations_ham_stirred[diluted_indices],
#              xerr=concentrations_ham_stirred[diluted_indices] * relative_error,
#              y=concentrations_ham_static[diluted_indices],
#              yerr=concentrations_ham_static[diluted_indices] * relative_error,
#              marker='o', color='C0', linestyle='None',
#          label='Calculated with dilution', alpha=0.5)
plt.errorbar(x=concentrations_ham_stirred[undiluted_indices] / 2 * 5 * 0.9 / dilution_factor,
             xerr=concentrations_ham_stirred[undiluted_indices] / 2 * 5 * 0.9 / dilution_factor * relative_error,
         y=concentrations_ham_static[undiluted_indices] / 2 * 5 * 0.9 / dilution_factor,
         yerr=concentrations_ham_static[undiluted_indices] / 2 * 5 * 0.9 / dilution_factor * relative_error,
         marker='o', color='C1', linestyle='None',
         label='Calculated without dilution',
         alpha=0.5)
# plt.plot([0, 0.03], [0, 0.03], color='black')
plt.plot([0, 0.03], [0, 0.03], color='black', label='x=y')
plt.xlabel('Stirred during reaction in Eppendorf,\nIIO029A concentration in mol/L')
plt.ylabel('Static during reaction in Eppendorf,\nIIO029A concentration in mol/L')
# plt.legend()
plt.title('Stirring vs. not stirring during the reaction')
plt.tight_layout()
fig.savefig(data_folder + experiment_name + 'results/mixing_validation/comparison_stirred_vs_static_Ep.png', dpi=300)
plt.show()
#
fig2 = plt.figure(2)
plt.plot(concentrations_main[undiluted_indices] / 2 * 5 * 0.9 / dilution_factor,
         concentrations_ham_static[undiluted_indices] / 2 * 5 * 0.9 / dilution_factor, 'o', color='C0',
         label='Calculated without dilution',
         alpha=0.5)
plt.plot(concentrations_main[diluted_indices],
         concentrations_ham_static[diluted_indices], 'o', color='C1',
         label='Calculated with dilution',
         alpha=0.5)
# plt.plot([0, 0.03], [0, 0.03], color='black')
plt.plot([0, 0.03], [0, 0.03], color='black', label='x=y')
plt.xlabel('Reaction in large chamber,\nIIO029A concentration in mol/L')
plt.ylabel('Reaction in small chamber,\nIIO029A concentration in mol/L')
plt.legend()
plt.tight_layout()
fig.savefig(data_folder + experiment_name + 'results/mixing_validation/comparison_from_two_chambers.png', dpi=300)
plt.show()


concentrations_reprod = process_or_load_processed(nickname='concentrations_reprod',
                                                  plate_name='2023-03-28_13-49-08__plate0000021__reprod-test',
                                                  force_recalculation=False)[diluted_indices]
rel_std_reprod = print_stats(concentrations_reprod, label='Same condition, large oven')

fig3 = plt.figure(3)
plt.imshow(concentrations_reprod.reshape((3, 9)))
plt.title('Same condition, large oven')
plt.show()

concentrations_reprod_hotpress = process_or_load_processed(nickname='concentrations_reprod_hotpress',
    plate_name='2023-03-30_10-41-06__plate0000009__multicomponent_0320_27_same_reactions_seal_check')[diluted_indices]
rel_std_reprod_hotpress = print_stats(concentrations_reprod_hotpress, label='Same condition, capped better and with very uniform temperature')
fig31 = plt.figure(31)
plt.imshow(concentrations_reprod_hotpress.reshape((3, 9)))
plt.title('Same condition, capped better and with very uniform temperature')
plt.show()

concentrations_reprod_oven_2 = process_or_load_processed(nickname='concentrations_reprod_oven_2',
    plate_name='2023-03-30_10-58-00__plate0000010__multicomponent_0320_27_same_reactions_seal_check')[diluted_indices]
rel_std_reprod_oven_2 = print_stats(concentrations_reprod_oven_2, label='Same condition, large oven repetition')
fig32 = plt.figure(32)
plt.imshow(concentrations_reprod_oven_2.reshape((3, 9)))
plt.title('Same condition, large oven repetition')
plt.show()

x = [concentrations_reprod_hotpress,
     concentrations_reprod_oven_2,
     concentrations_reprod]
df = pd.DataFrame(x, index=[f'Same conditions, capped better and with very uniform temperature, STD {rel_std_reprod_hotpress:.1%}',
                            f'Same conditions, large oven, repetition, STD {rel_std_reprod_oven_2:.1%}',
                            f'Same conditions, large oven, STD {rel_std_reprod:.1%}'], dtype=object)
figbox = plt.figure(34, figsize=(14,3))
df.T.boxplot(vert=False)
plt.xlabel('Concentration of product, mol/L')
# plt.subplots_adjust(left=0.25)
plt.tight_layout()
plt.show()


concentrations_identical = process_or_load_processed(nickname='concentrations_identical', plate_name='2023-03-28_15-20-56__plate0000022__reprod-test-samesample')
print(f'Identical: Mean: {np.mean(concentrations_identical)}, std: {np.std(concentrations_identical)}, rel. std: {np.std(concentrations_identical) / np.mean(concentrations_identical)}')
as_plate = concentrations_identical.reshape((6, 9))

fig4, ax4 = plt.subplots()
plt.scatter(np.arange(54), concentrations_identical, label='Identical solution in all vials, immediately measured')
plt.title('Identical solution in all vials')
plt.xlabel('Vial ID')
plt.ylabel('Product concentration in mol/L')

fig5 = plt.figure(5)
plt.title('Identical solution in all vials, immediately measured')
plt.imshow(as_plate)
plt.show()


concentrations_identical_evap = process_or_load_processed(nickname='concentrations_identical_evap', plate_name='2023-03-29_13-39-48__plate0000022__evap-check-samesample')
print(f'Identical after evap: Mean: {np.mean(concentrations_identical_evap)}, std: {np.std(concentrations_identical_evap)}, rel. std: {np.std(concentrations_identical_evap) / np.mean(concentrations_identical_evap)}')
as_plate = concentrations_identical_evap.reshape((6, 9))

ax4.scatter(np.arange(54), concentrations_identical_evap, label='Identical solution in all vials, after 16 hours in oven')
ax4.legend()
# plt.title('Identical solution in all wells, after 16 hours in oven')

fig7 = plt.figure(7)
plt.imshow(as_plate)
plt.title('Identical solution in all vials, after 16 hours in oven')
plt.show()


concentrations_identical_diluted = process_or_load_processed(nickname='concentrations_identical_diluted', plate_name='2023-03-29_15-04-17__plate0000019__reprod-check-dilution')
concentrations_identical_diluted = concentrations_identical_diluted[diluted_indices]
concentrations_identical_diluted = concentrations_identical_diluted[:9]
# print(f'Identical, diluted and measured: Mean: {np.mean(concentrations_identical_diluted)}, std: {np.std(concentrations_identical_diluted)}, rel. std: {np.std(concentrations_identical_diluted) / np.mean(concentrations_identical_diluted)}')
print_stats(concentrations_identical_diluted, label='Identical solution in all vials, diluted and measured')
fig8 = plt.figure(8)
plt.scatter(np.arange(9), concentrations_identical_diluted, label='Identical solution in all vials, diluted and measured')
plt.title('Identical solution in all vials, diluted and measured')
plt.xlabel('Vial ID')
plt.ylabel('Product concentration in mol/L')
plt.legend()



concentrations_identical_stirred = process_or_load_processed(nickname='concentrations_identical_stirred', plate_name='2023-04-04_11-31-08__plate0000003__2023-04-03_run01')
concentrations_identical_stirred = concentrations_identical_stirred[diluted_indices]
# print(f'Identical, stirred in Eppendorf during reaction: Mean: {np.mean(concentrations_identical_diluted)}, std: {np.std(concentrations_identical_diluted)}, rel. std: {np.std(concentrations_identical_diluted) / np.mean(concentrations_identical_diluted)}')
print_stats(concentrations_identical_stirred, label='Identical, stirred in Eppendorf during reaction')
fig9 = plt.figure(9)
plt.scatter(np.arange(concentrations_identical_stirred.shape[0]), concentrations_identical_stirred, label='Identical, stirred in Eppendorf during reaction')
plt.title('Identical, stirred in Eppendorf during reaction')
plt.xlabel('Vial ID')
plt.ylabel('Product concentration in mol/L')
plt.legend()


concentrations_identical_capped_in_oven = process_or_load_processed(nickname='concentrations_identical_capped_in_oven', plate_name='2023-04-04_12-15-41__plate0000006__2023-04-03_run01')
concentrations_identical_capped_in_oven = concentrations_identical_capped_in_oven[diluted_indices]
print_stats(concentrations_identical_capped_in_oven, label='Identical, capped between Al plates in large oven')
fig10 = plt.figure(10)
plt.scatter(np.arange(concentrations_identical_capped_in_oven.shape[0]), concentrations_identical_capped_in_oven, label='Identical, capped between Al plates in large oven')
plt.title('Identical, capped between Al plates in large oven')
plt.xlabel('Vial ID')
plt.ylabel('Product concentration in mol/L')
plt.legend()
plt.show()

x = [concentrations_identical_capped_in_oven/(0.0249922896339447/0.21741)*100,
     concentrations_identical/(0.02295312069974119/0.212)*100]
# df = pd.DataFrame(x, index=[f'Repeated experimental condition, N=27, yield (21.7 ± 1.2)%',
#                             f'Repeated measurements of same crude, N=54, yield (21.2 ± 0.4)%'], dtype=object)
df = pd.DataFrame(x, index=[f'Replicates',
                            f'UV-VIS repeats'], dtype=object)
figbox = plt.figure(69, figsize=(4.9*0.9,0.9*4.9/3*0.9), dpi=300)
bp = df.T.boxplot(vert=False, patch_artist=True, widths=(0.6, 0.6), return_type='dict')
for median in bp['medians']:
    median.set_color('black')
    median.set_linewidth(3)
for box in bp['boxes']:
    box.set_facecolor('grey')
    box.set(color='black')
    # change line oclors
    box.set_alpha(0.3)
for whisker in bp['whiskers']:
    whisker.set(color="grey")

plt.gca().grid(False)
def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
simpleaxis(plt.gca())

for cap in bp['caps']:
    cap.set(color="grey")
plt.xlabel('Yield, %')
# plt.subplots_adjust(left=0.25)
plt.tight_layout()
figbox.savefig('misc-scripts/figures/Figure_69.png', dpi=300)
plt.show()

concentrations_aluminum_plate_2023_06_07 = process_or_load_processed(nickname='concentrations_aluminum_plate_2023_06_07', plate_name='2023-06-08_12-41-45__plate0000006__multicomponent-reactions-2023-06-08-dil')
concentrations_aluminum_plate_2023_06_07 = concentrations_aluminum_plate_2023_06_07[diluted_indices]
print_stats(concentrations_aluminum_plate_2023_06_07, label='Identical reactions, in Al vial-plate in blue oven')
fig10 = plt.figure(10)
plt.scatter(np.arange(concentrations_aluminum_plate_2023_06_07.shape[0]), concentrations_aluminum_plate_2023_06_07, label='Identical reactions, in Al vial-plate in blue oven')
plt.title('Identical reactions, in Al vial-plate in blue oven. Rel. STD 20.5')
plt.xlabel('Vial ID')
plt.ylabel('Product concentration in mol/L')
plt.ylim(0, max(concentrations_aluminum_plate_2023_06_07)*1.4)
plt.legend()
plt.show()

as_plate = concentrations_aluminum_plate_2023_06_07.reshape((3, 9))
fig7 = plt.figure(7)
plt.imshow(as_plate)
plt.title('Identical reactions, in Al vial-plate in blue oven. Rel. STD 20.5')
plt.show()