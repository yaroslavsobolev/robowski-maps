import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from robowski.settings import *
from robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra import load_msp_file

def plot_absorbance_threshold():
    """
    Plot the absorbance threshold boundary that separates usable from poor quality spectral data.

    This visualization shows the wavelength-dependent absorbance limits used in spectral unmixing
    to exclude regions with instrumental artifacts or non-linear Beer-Lambert behavior.
    """
    # Define threshold parameters
    thresh_w_indices = [0, 25, 127, 2000]
    thresh_as = [0.67, 0.75, 1.6, 1.6]

    # Create interpolator
    threshold_interpolator = interpolate.interp1d(thresh_w_indices, thresh_as, fill_value='extrapolate')

    # Create wavelength arrays
    # Index 0 corresponds to 350 nm, each index represents 1 nm
    wavelength_indices = np.arange(0, 500)  # Cover reasonable range for UV-Vis
    wavelengths = wavelength_indices + 350  # Convert indices to actual wavelengths (nm)

    # Calculate threshold values
    threshold_values = threshold_interpolator(wavelength_indices)

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot the threshold boundary line
    ax.plot(wavelengths, threshold_values, 'k-', linewidth=2, label='Quality threshold')

    # Mark the original data points
    original_wavelengths = np.array(thresh_w_indices) + 350
    ax.plot(original_wavelengths, thresh_as, 'ko', markersize=8)

    # Fill regions
    y_max = max(2.0, np.max(threshold_values) * 1.2)  # Upper limit for plot

    # Fill "Poor quality data" region (above threshold)
    ax.fill_between(wavelengths, threshold_values, y_max,
                    alpha=0.3, color='red', label='Region of poor data quality')

    # Fill "Usable data" region (below threshold)
    ax.fill_between(wavelengths, 0, threshold_values,
                    alpha=0.3, color='green', label='Region of good data quality')

    # Customize the plot
    ax.set_xlabel('Wavelength (nm)', fontsize=12)
    ax.set_ylabel('Absorbance', fontsize=12)
    ax.set_title('Spectral data quality threshold used for CRAIC spectrophotometer', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)

    # Set reasonable axis limits
    ax.set_xlim(350, 850)
    ax.set_ylim(0, y_max)

    # # Add annotations for key points
    # for i, (w_idx, abs_val) in enumerate(zip(thresh_w_indices, thresh_as)):
    #     wavelength = w_idx + 350
    #     if wavelength <= 850:  # Only annotate points within plot range
    #         ax.annotate(f'({wavelength} nm, {abs_val})',
    #                     xy=(wavelength, abs_val),
    #                     xytext=(10, 10),
    #                     textcoords='offset points',
    #                     fontsize=9,
    #                     bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    plt.tight_layout()
    return fig


def plot_threshold_with_typical_spectrum():
    """
    Enhanced plot showing the threshold along with a typical absorption spectrum
    to demonstrate how the threshold works in practice.
    """
    # Define threshold parameters
    thresh_w_indices = [0, 25, 127, 2000]
    thresh_as = [0.67, 0.75, 1.6, 1.6]
    threshold_interpolator = interpolate.interp1d(thresh_w_indices, thresh_as, fill_value='extrapolate')

    # Create wavelength arrays
    wavelength_indices = np.arange(0, 500)
    wavelengths = wavelength_indices + 350
    threshold_values = threshold_interpolator(wavelength_indices)

    filename = data_folder + 'craic_microspectrometer_measurements/absorbance/2023-06-24_09-24-22__plate0000049__multicomp_reactions_2023-06-23-run01/spectrum_-7.msp'
    synthetic_spectrum = load_msp_file(filename)

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 7))

    # Plot the synthetic spectrum
    ax.plot(synthetic_spectrum[:, 0], synthetic_spectrum[:, 1], 'b-', linewidth=2, label='Typical Absorption Spectrum')

    # Plot the threshold boundary
    ax.plot(wavelengths, threshold_values, 'r--', linewidth=2, label='Quality Threshold')

    # Fill regions
    y_max = max(2.0, np.max(threshold_values) * 1.1)

    # Fill "Poor quality data" region (above threshold)
    ax.fill_between(wavelengths, threshold_values, y_max,
                    alpha=0.2, color='red', label='Poor Quality Region')

    # Fill "Usable data" region (below threshold)
    ax.fill_between(wavelengths, 0, threshold_values,
                    alpha=0.2, color='green', label='Usable Data Region')

    # Highlight where spectrum would be rejected
    mask_rejected = synthetic_spectrum[:, 1] > threshold_interpolator(synthetic_spectrum[:,0])
    if np.any(mask_rejected):
        ax.fill_between(synthetic_spectrum[:,0], 0, synthetic_spectrum[:, 1],
                        where=mask_rejected, alpha=0.5, color='orange',
                        label='Spectrum Data Rejected')

    # Customize the plot
    ax.set_xlabel('Wavelength (nm)', fontsize=12)
    ax.set_ylabel('Absorbance', fontsize=12)
    ax.set_title('Spectral data quality threshold used for CRAIC spectrophotometer', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    ax.set_xlim(350, 700)
    ax.set_ylim(0, y_max)

    plt.tight_layout()
    return fig


if __name__ == "__main__":
    # Create both plots
    fig1 = plot_absorbance_threshold()
    # plt.figure()
    # fig2 = plot_threshold_with_typical_spectrum()
    plt.show()