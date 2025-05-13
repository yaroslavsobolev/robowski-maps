from robowski.settings import *
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

from scipy import signal
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'


def load_and_filter_spectrum(file_name, do_plot=False):
    fs = 1  # Sample frequency (Hz)
    f0 = 3.16761394e-02  # Frequency to be removed from signal (Hz)
    Q = 1  # Quality factor
    # Design notch filter
    b, a = signal.iirnotch(f0, Q, fs)

    # freq, h = signal.freqz(b, a, fs=fs)
    #
    # # Plot
    # # Plot magnitude response of the filter
    # plt.plot(freq * fs / (2 * np.pi), 20 * np.log10(abs(h)),
    #          'r', label='Bandpass filter', linewidth='2')
    #
    # plt.xlabel('Frequency [Hz]', fontsize=20)
    # plt.ylabel('Magnitude [dB]', fontsize=20)
    # plt.title('Notch Filter', fontsize=20)
    # plt.grid()
    # plt.show()

    spectrum = np.loadtxt(file_name, skiprows=2, delimiter=',')
    outputSignal = signal.filtfilt(b, a, spectrum[:, 1])
    if do_plot:
        plt.plot(spectrum[:, 0], spectrum[:, 1], label='Raw spectrum')
        plt.plot(spectrum[:, 0], outputSignal, color='red', label='Filtered spectrum')
        plt.gca().invert_xaxis()
        plt.xlabel('Wavenumber, cm$^{-1}$')
        plt.ylabel('Absorbance in 0.1 mm path')
        plt.legend()
        plt.show()
    spectrum[:, 1] = outputSignal
    return spectrum


if __name__ == '__main__':
    # load_and_filter_spectrum(f'{data_folder}Yaroslav/wet_acetonitrile/ftir/water-in-acn/trans/tr-sc4-water-in-acn-100mMw_rep2.csv', do_plot=True)
    target_frequency = 3630
    spectrum_100mM = load_and_filter_spectrum(
        f'{data_folder}Yaroslav/wet_acetonitrile/ftir/water-in-acn/trans/tr-sc4-water-in-acn-100mMw_rep2.csv')
    spectrum_1000mM = load_and_filter_spectrum(
        f'{data_folder}Yaroslav/wet_acetonitrile/ftir/water-in-acn/trans/tr-sc4-water-in-acn-1000mMw_rep1.csv')
    wavenumbers = spectrum_100mM[:, 0]
    # index where wavenumber is closest to target_frequency
    band_index = np.argmin(np.abs(wavenumbers - target_frequency))
    def band_a(spectrum):
        # index where wavenumber is closest to 3900
        ref_band_index = np.argmin(np.abs(wavenumbers - 3900))
        return spectrum[band_index, 1] - spectrum[ref_band_index, 1]

    xs = np.array([0, 100, 1098.28])
    ys = np.array([0, band_a(spectrum_100mM), band_a(spectrum_1000mM)])
    plt.scatter(xs, ys, color='black')
    # fit a line passing through zero
    def func(x, a):
        return a * x
    popt, pcov = curve_fit(func, xs, ys)
    print(popt)
    perr = np.sqrt(np.diag(pcov))
    print(f'Relative error {perr[0]/popt[0]}')
    rel_err = perr[0]/popt[0] + 0.03
    plt.plot(xs, func(xs, *popt), color='black')
    plt.xlabel('Water concentration, mM')
    plt.ylabel('Absorbance at 3630 cm$^{-1}$ for 0.1 mm path')
    plt.show()

    spectrum_dry_acn = load_and_filter_spectrum(f'{data_folder}Yaroslav/wet_acetonitrile/ftir/water-in-acn/trans/tr-sc16-water-in-acn_reactionDry_rep3.csv')
    concentration_here = band_a(spectrum_dry_acn) / popt[0]
    print(f'Concentration of water in dry acetonitrile: {concentration_here} mM +- {rel_err * concentration_here}')

    spectrum_wet_stock = load_and_filter_spectrum(f'{data_folder}Yaroslav/wet_acetonitrile/ftir/water-in-acn/trans/tr-sc4-water-in-acn_reactionWet_rep1.csv')
    concentration_here = band_a(spectrum_wet_stock) / popt[0]
    print(f'Concentration of water in wet stock: {concentration_here} mM +- {rel_err * concentration_here}')

    spectrum_wet_vial = load_and_filter_spectrum(f'{data_folder}Yaroslav/wet_acetonitrile/ftir/water-in-acn/trans/tr-sc4-water-in-acn_reactionVialWet_rep2.csv')
    concentration_here = band_a(spectrum_wet_vial) / popt[0]
    print(f'Concentration of water in wet vial: {concentration_here} mM +- {rel_err * concentration_here}')

    spectrum_wet_vial_2 = load_and_filter_spectrum(f'{data_folder}Yaroslav/wet_acetonitrile/ftir/water-in-acn/trans/tr-sc4-water-in-acn_reactionVialWet_rep1.csv')
    concentration_here = band_a(spectrum_wet_vial_2) / popt[0]
    print(f'Concentration of water in wet vial 2: {concentration_here} mM +- {rel_err * concentration_here}')

    spectrum_wet_vial_3 = load_and_filter_spectrum(f'{data_folder}Yaroslav/wet_acetonitrile/ftir/water-in-acn/trans/tr-sc4-water-in-acn_reactionVialWet_23m_85RH_nocov_rep1.csv')
    concentration_here = band_a(spectrum_wet_vial_3) / popt[0]
    print(f'Concentration of water in wet vial 3, 23 minutes in 85% RH: {concentration_here} mM +- {rel_err * concentration_here}')

    spectrum_wet_vial_4 = load_and_filter_spectrum(f'{data_folder}Yaroslav/wet_acetonitrile/ftir/water-in-acn/trans/tr-sc4-water-in-acn_reactionVialWet_36m_85RH_cov_rep1.csv')
    concentration_here = band_a(spectrum_wet_vial_4) / popt[0]
    print(f'Concentration of water in wet vial 4, 36 minutes in 85% RH: {concentration_here} mM +- {rel_err * concentration_here}')

    # plt.plot(spectrum_100mM[:, 0], spectrum_100mM[:, 1], label='100 mM water (reference for calibration)')
    plt.plot(spectrum_dry_acn[:, 0], spectrum_dry_acn[:, 1], label='Dry acetonitrile')
    # plt.plot(spectrum_wet_stock[:, 0], spectrum_wet_stock[:, 1], label='Stock solution after 45 minutes in 90% relative humidity')
    # plt.plot(spectrum_wet_vial[:, 0], spectrum_wet_vial[:, 1], label='0.5 mL in vial after 45 minutes in 90% relative humidity')
    # plt.plot(spectrum_wet_vial[:, 0], spectrum_wet_vial[:, 1],
    #          label='0.5 mL in vial after 75 minutes in 99% relative humidity')
    plt.plot(spectrum_wet_vial_3[:, 0], spectrum_wet_vial_3[:, 1], label='0.5 mL in vial after 23 minutes in 85% relative humidity')
    plt.plot(spectrum_wet_vial_4[:, 0], spectrum_wet_vial_4[:, 1], label='0.5 mL in vial after 36 minutes in 85% relative humidity')
    # plt.plot(spectrum_1000mM[:, 0], spectrum_1000mM[:, 1], label='1000 mM water')
    plt.legend()
    plt.xlabel('Wavenumber, cm$^{-1}$')
    plt.ylabel('Absorbance in 0.1 mm path')
    plt.gca().invert_xaxis()
    plt.show()

# plt.show()