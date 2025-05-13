from robowski.settings import *
import os
from matplotlib import pyplot as plt
from scipy import interpolate
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def read_cary_agilent_csv_spectrum(cary_file, column_name):
    # load cary_file as text file and find the index of the first line that is empty
    with open(cary_file, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if line == '\n':
                break
    df = pd.read_csv(cary_file, skipfooter=len(lines)-i, engine='python')
    df = pd.read_csv(cary_file, skiprows=2, names=df.columns, skipfooter=len(lines)-i, engine='python')
    wavelengths = df[column_name]
    # get next column after the column_name
    column_index = df.columns.get_loc(column_name)

    next_column_name = df.columns[column_index + 1]
    ys = df[next_column_name]
    return wavelengths, ys


def write_cary_agilent_csv_spectrum(wavelengths, spectrum, output_file_path, column_name):
    # create a new txt file with the same name as the cary_file
    # write the first two lines of the cary_file
    with open(output_file_path, 'w+') as cary_file:
        # write into the first line column_name followed by a comma
        cary_file.write(column_name + ',\n')
        # write into the second line the 'Wavelength (nm)' followed by a comma and 'Abs'
        cary_file.write('Wavelength (nm),Abs\n')
        # write the rest of the file, comma-separated values of wavelength and spectrum
        for i in range(len(wavelengths)):
            cary_file.write(f'{wavelengths[i]},{spectrum[i]}\n')
        cary_file.write('\n')
        cary_file.write('Written by spectraltools.write_cary_agilent_csv_spectrum()\n')


def stitch_two_spectra(wavelengths1, spectrum1, wavelengths2, spectrum2, absorbance_limit=1.0, do_plot=False):
    spectrum2_interpolator = interpolate.interp1d(wavelengths2, spectrum2, fill_value='extrapolate')

    # find the wavelengths where the second spectrum is below absorbance_limit
    mask = spectrum2_interpolator(wavelengths1) < absorbance_limit

    def func(xs, a, b):
        return a * spectrum2_interpolator(xs) + b

    p0 = (1, 0)
    bounds = ([-1e-20, -np.inf], [np.inf, np.inf])
    popt, pcov = curve_fit(func, wavelengths1[mask], spectrum1[mask],
                           p0=p0, bounds=bounds)


    fig1 = plt.figure(1)
    plt.plot(wavelengths1, spectrum1, label='spectrum1', color='C0', alpha=0.5)
    plt.plot(wavelengths2, spectrum2, label='spectrum2', color='C1', alpha=0.5)
    plt.plot(wavelengths1[mask], func(wavelengths1[mask], *popt), label='fit', color='C3', alpha=0.5)
    # set y to log scale
    plt.yscale('log')
    plt.legend()
    if do_plot:
        plt.show()
    else:
        plt.clf()

    final_spectrum = np.copy(spectrum1)
    final_spectrum[mask] = func(wavelengths1[mask], *popt)
    final_spectrum -= np.min(final_spectrum)

    # plot new ref spectrum in semilog scale
    plt.semilogy(wavelengths1, final_spectrum)
    plt.title(f'Final spectrum')
    plt.show()

    return wavelengths1, final_spectrum

if __name__ == '__main__':
    pass