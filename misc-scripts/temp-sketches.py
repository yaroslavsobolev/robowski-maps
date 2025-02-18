# import logging
# import os
# import numpy as np
# import re
# from matplotlib import pyplot as plt
# import pandas as pd


from multiprocessing import Pool

myvar = 1

def f(x):
    return myvar * 2

if __name__ == '__main__':
    myvar = input('Enter a number: ')

    with Pool(5) as p:
        print(p.map(f, [1, 2, 3]))

#
# data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
#
#
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.interpolate import SmoothBivariateSpline
#
# import warnings
# warnings.simplefilter('ignore')
#
# # train_x, train_y = np.meshgrid(np.arange(-5, 5, 0.5), np.arange(-5, 5, 0.5))
# # train_x = train_x.flatten()
# # train_y = train_y.flatten()
#
# train_x = np.random.uniform(-5, 5, 300)
# train_y = np.random.uniform(-5, 5, 300)
#
#
# def z_func(x, y):
#     return np.cos(x) + np.sin(y) ** 2 + 0.05 * x + 0.1 * y
#
# train_z = z_func(train_x, train_y)
# interp_func = SmoothBivariateSpline(train_x, train_y, train_z, s=0.0)
# smth_func = SmoothBivariateSpline(train_x, train_y, train_z)
#
# test_x = np.arange(-9, 9, 0.01)
# test_y = np.arange(-9, 9, 0.01)
# grid_x, grid_y = np.meshgrid(test_x, test_y)
#
# interp_result = interp_func(test_x, test_y).T
# smth_result = smth_func(test_x, test_y).T
# perfect_result = z_func(grid_x, grid_y)
#
# fig, axes = plt.subplots(1, 3, figsize=(16, 8))
# extent = [test_x[0], test_x[-1], test_y[0], test_y[-1]]
# opts = dict(aspect='equal', cmap='nipy_spectral', extent=extent, vmin=-1.5, vmax=2.5)
#
# im = axes[0].imshow(perfect_result, **opts)
# fig.colorbar(im, ax=axes[0], orientation='horizontal')
# axes[0].plot(train_x, train_y, 'w.')
# axes[0].set_title('Perfect result, sampled function', fontsize=21)
#
# im = axes[1].imshow(smth_result, **opts)
# axes[1].plot(train_x, train_y, 'w.')
# fig.colorbar(im, ax=axes[1], orientation='horizontal')
# axes[1].set_title('s=default', fontsize=21)
#
# im = axes[2].imshow(interp_result, **opts)
# fig.colorbar(im, ax=axes[2], orientation='horizontal')
# axes[2].plot(train_x, train_y, 'w.')
# axes[2].set_title('s=0', fontsize=21)
#
# plt.tight_layout()
# plt.show()


# experiment_name = f'BPRF/2024-01-17-run01/'
#
# def read_cary_agilent_csv_spectrum(cary_file, column_name):
#     df = pd.read_csv(cary_file)
#     df = pd.read_csv(cary_file, skiprows=2, names=df.columns)
#     wavelengths = df[column_name]
#     # get next column after the column_name
#     column_index = df.columns.get_loc(column_name)
#
#     next_column_name = df.columns[column_index + 1]
#     ys = df[next_column_name]
#     return wavelengths, ys
#
# cary_file = data_folder + experiment_name + 'calibrations/spectrophotometer_data/Hantzsch-ester-HRP01/HRP01_400ug_per_20mL_repeat1.csv'
# column_name = 'HRP01_0.4mg_per_20_mL_repeat1'
# # cary_file = data_folder + experiment_name + 'calibrations/spectrophotometer_data/Hantzsch_dm37/dm37.csv'
# # column_name = 'dm_37_SBW1nm_repeat2'
# plt.title('HRP01 reference spectrum')
#
# wavelengths, ys = read_cary_agilent_csv_spectrum(cary_file, column_name)
# plt.plot(wavelengths, ys, label='reference from Agilent Cary')
#
# spectrum = np.load(data_folder + '/BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/HRP01/ref_spectrum.npy')
# # spectrum = np.load(data_folder + '/BPRF/2024-01-17-run01/microspectrometer_data/calibration/references/dm37/ref_spectrum.npy')
# wavelengths = 220 + np.arange(spectrum.shape[0])
#
# plt.plot(wavelengths, spectrum*2*0.98, label='reference from nanodrop')
# # plt.plot(wavelengths, spectrum, label='reference from nanodrop')
# plt.legend()
# plt.xlabel('Wavelength, nm')
# plt.ylabel('Absorbance')
#
# plt.show()



# import statsmodels.api as sm
#
#
# xs = np.linspace(0, 2*np.pi*10, 1000)
# ys = 1*np.sin(xs/100) + np.random.normal(0, 0.1, 1000)
# x = sm.stats.acorr_ljungbox(ys, lags=[100])
# print(x)
# plt.plot(xs, ys)
# plt.show()


# from icecream import ic
#
# data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
#
# list_of_runs = tuple([
#     '2023-08-21-run01',
#     '2023-08-22-run01',
#     '2023-08-22-run02',
#     '2023-08-28-run01',
#     '2023-08-29-run01',
#     '2023-08-29-run02'] +
# [
#                           '2023-09-06-run01',
#                           '2023-09-07-run01'
#     ])
#
# def get_missing_for_one_folder(run_shortname):
#     experiment_name = data_folder + 'simple-reactions/' + run_shortname + '/'
#     nanodrop_spectra_folder = experiment_name + 'nanodrop_spectra/'
#     pattern = re.compile(r'(?P<timestamp_string>\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})_plate_(?P<plate_id>\d+).csv')
#     for filename in os.listdir(nanodrop_spectra_folder):
#         if not filename.endswith('.csv'):
#             continue
#         match = pattern.match(filename)
#         if not match:
#             logging.warning(f'Could not parse nanodrop spectrum filename {filename} for timestamp and id')
#             continue
#         nanodrop_spectrum_file = experiment_name + 'nanodrop_spectra/' + filename
#         nanodrop_df = pd.read_csv(nanodrop_spectrum_file)
#         columns = nanodrop_df.columns
#         index_list = [str(x) for x in range(54)]
#         flag_list = [x in columns for x in index_list]
#         while flag_list[-1] == False:
#             flag_list.pop()
#         count_missing = flag_list.count(False)
#         indices_of_missing = [i for i, x in enumerate(flag_list) if not x]
#         missing_elements = [index_list[i] for i in indices_of_missing]
#         if count_missing > 0:
#             print(f'run {run_shortname} has {count_missing} missing column in {filename}, indices: {missing_elements}')
#             # add a zero-filled column with missing label and write dataframe back to csv
#             nanodrop_df = nanodrop_df.rename(columns={nanodrop_df.columns[0]: ""})
#             nanodrop_df[missing_elements[0]] = 0
#             nanodrop_df.to_csv(nanodrop_spectrum_file, index=False)
#     return count_missing
#
# # for run_shortname in list_of_runs:
# #     get_missing_for_one_folder(run_shortname)
#
# # get_missing_for_one_folder('2023-09-06-run01')
#
# d = {'key': {1: 'one'}}
# ic(d['key'][1])
#
# class klass():
#     attr = 'yep'
# ic(klass.attr)

#
# def func(*args):
#     x = args[0]
#     params = args[1:]
#     result = 0
#     for i, param in enumerate(params):
#         result += param * x**i
#     return result
#
# # example of curve_fit
# from scipy.optimize import curve_fit
# xdata = np.linspace(0, 4, 50)
# y = func(xdata, 2.5, 1.3, 0.5, -0.2)
# ydata = y + 0.2 * np.random.normal(size=len(xdata))
# popt, pcov = curve_fit(func, xdata, ydata, p0=[0, 0, 0, 3])
# print(popt)
# print(pcov)
# # plot fitted curve
# import matplotlib.pyplot as plt
# plt.plot(xdata, ydata, 'b-', label='data')
# plt.plot(xdata, func(xdata, *popt), 'r-', label='fit')
# plt.legend()
# plt.show()
#
# import numpy as np
# import matplotlib.pyplot as plt
#
# c=3.e2
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax2 = ax1.twiny()
#
# xvals = np.arange(199.9, 999.9, 0.1)
# data = np.sin(0.03*xvals)
# ax1.plot(xvals, data)
#
# ax1Ticks = ax1.get_xticks()
# ax2Ticks = ax1Ticks
#
# def tick_function(X):
#     V = c/X
#     return ["%.3f" % z for z in V]
#
# ax2.set_xticks(ax2Ticks)
# ax2.set_xbound(ax1.get_xbound())
# ax2.set_xticklabels(tick_function(ax2Ticks))
#
# ax1.set_xlabel("Frequency (GHz)")
# ax2.set_xlabel('Wavelength (mm)')
# ax1.grid(True)
# plt.ylim(ymin=-1.1,ymax=1.1)
# plt.show()