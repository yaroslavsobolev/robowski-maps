import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.interpolate import SmoothBivariateSpline

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

errorbar_folder = f'{data_folder}/nanodrop-spectrophotometer-measurements/nanodrop_errorbar_folder_2024-03-16/'

raw_residuals_folder = errorbar_folder + 'raw_residuals/'

# find all npy files in the folder
npy_files = [f for f in os.listdir(raw_residuals_folder) if f.endswith('.npy')]

# make a dataset for subsequent 2d griddate
residuals_thresh = 0.001
xs = []
ys = []
zs = []
for npy_file in npy_files:
    print(npy_file)
    npy_data = np.load(raw_residuals_folder + npy_file)
    print(npy_data.shape)
    wavelengths = npy_data[:, 0]
    absorbances = npy_data[:, 1]
    residuals = npy_data[:, 2]

    # leave only data where residuals are above thresh
    wavelengths = wavelengths[np.abs(residuals) > residuals_thresh]
    absorbances = absorbances[np.abs(residuals) > residuals_thresh]
    residuals = residuals[np.abs(residuals) > residuals_thresh]

    plt.scatter(wavelengths, absorbances, c='b', label='absorbance')
    xs.extend(wavelengths)
    ys.extend(absorbances)
    zs.extend(residuals)

for i in range(10):
    xs.append(220)
    ys.append(0.12)
    zs.append(0.04)

for i in range(4):
    xs.append(232)
    ys.append(0.2)
    zs.append(0.04)

for i in range(6):
    xs.append(220)
    ys.append(0.21)
    zs.append(0.03)

plt.show()
xs = np.array(xs)
ys = np.array(ys)
zs = np.array(zs)**2
# do griddata() and then plot 2d imshow style (like pcolor or pcolormesh)
from scipy.interpolate import griddata
# define grid
grid_x, grid_y = np.mgrid[np.min(xs):np.max(xs):30j, np.min(ys):np.max(ys):20j]
grid_z = griddata((xs, ys), zs, (grid_x, grid_y), method='linear')
grid_z = np.sqrt(grid_z)
# use pcolormesh
plt.figure()
plt.pcolormesh(grid_x, grid_y, grid_z, vmax=0.1, cmap='viridis', shading='auto')


plt.xlabel('Wavelength, nm')
plt.ylabel('Absorbance')
plt.title('Residuals')
# colorbar with label
plt.colorbar(label='Root-mean square residuals')
plt.show()



def func(x, a, b, c, d, e, f):
    return a + b*x**2 + c*np.exp(d*x) + e*np.exp(-f*x)

# plot zs where x is between 220 and 240
do_show = False

fitted_xs = []
fitted_ys = []
fitted_zs = []

wavelength_step = 20
for wavelength in range(220, 600, wavelength_step):

    ys_here = ys[(xs > wavelength) & (xs <= wavelength + wavelength_step) & (ys<1)]
    zs_here = zs[(xs > wavelength) & (xs <= wavelength + wavelength_step) & (ys<1)]

    # sort by ys_here
    sort_index = np.argsort(ys_here)
    ys_here = ys_here[sort_index]
    zs_here = zs_here[sort_index]

    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(func, ys_here, zs_here, bounds=([0, 0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]))

    # sample the fit
    ys_new = np.linspace(-0.1, 1, 100)
    zs_new = func(ys_new, *popt)

    fitted_xs.extend([wavelength + wavelength_step/2]*len(ys_new))
    fitted_ys.extend(ys_new)
    fitted_zs.extend(np.sqrt(zs_new))

    if do_show:
        plt.plot(ys_here, func(ys_here, *popt), 'r-')
        plt.scatter(ys_here, zs_here, c='b', alpha=0.3)
        plt.xlabel('Absorbance')
        plt.ylabel('Residuals squared')
        plt.title('Residuals where wavelength is between 220 and 240 nm')
        plt.show()

        # plot square roots of z and of fit
        plt.figure()
        plt.scatter(ys_here, np.sqrt(zs_here), c='b', alpha=0.3)
        plt.plot(ys_here, np.sqrt(func(ys_here, *popt)), 'r-')
        plt.xlabel('Absorbance')
        plt.ylabel('Residuals')
        plt.title('Residuals where wavelength is between 220 and 240 nm')
        plt.show()

plt.figure()

grid_x, grid_y = np.mgrid[np.min(fitted_xs):np.max(fitted_xs):400j, 0:np.max(fitted_ys):100j]
grid_z = griddata((fitted_xs, fitted_ys), fitted_zs, (grid_x, grid_y), method='linear')
# grid_z = np.sqrt(grid_z)
# use pcolormesh
plt.figure()
plt.pcolormesh(grid_x, grid_y, grid_z, vmax=0.1, cmap='viridis', shading='linear')


plt.xlabel('Wavelength, nm')
plt.ylabel('Absorbance')
plt.title('Residuals')
# colorbar with label
plt.colorbar(label='Root-mean square residuals')
plt.show()

train_x = xs[(ys < 1)]
train_y = ys[(ys < 1)]
train_z = zs[(ys < 1)]
#
# smth_func = SmoothBivariateSpline(train_x, train_y, train_z, s=len(train_x)/100*7, w=np.ones_like(train_x)*1/(0.03**2))

test_x = np.linspace(np.min(train_x), np.max(train_x), 400)
test_y = np.linspace(0, np.max(train_y), 100)
grid_x, grid_y = np.meshgrid(test_x, test_y)

# pickle the smth_func
interpolator_filename = f'{errorbar_folder}bivariate_spline_interpolator.pkl'
#

import pickle
# with open(interpolator_filename, 'wb') as f:
#     pickle.dump(smth_func, f)

# load from dump
with open(interpolator_filename, 'rb') as f:
    smth_func = pickle.load(f)

grid_z = smth_func(grid_x, grid_y, grid=False)

print(f'Min grid z: {np.min(grid_z)}')

plt.figure()
to_plot_z = np.copy(grid_z)
thresh = 0.005
to_plot_z[to_plot_z < thresh**2] = thresh**2
to_plot_z = np.sqrt(to_plot_z)
plt.pcolormesh(grid_x, grid_y, to_plot_z, cmap='viridis', shading='linear')

plt.xlabel('Wavelength, nm')
plt.ylabel('Absorbance')
plt.title('Residuals')
# colorbar with label
plt.colorbar(label='Root-mean square residuals')
plt.show()

print(np.sqrt(smth_func(225, 0.8)[0][0]))

