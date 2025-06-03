# UV-Vis Absorption Spectroscopy Analysis

This module provides tools for processing and analyzing UV-Vis absorption spectroscopy data from well plates, with a
focus on spectral unmixing to determine component concentrations in multi-component reaction mixtures.
For tutorial of a worked example of the calibration and spectral unmixing workflow, 
see the Jupyter Notebook `notebooks/spectral_unmixing_tutorial.ipynb` notebook.

## Overview

Spectral unmixing is a mathematical process to decompose a mixed spectrum into its constituent components. The core
assumption is that the absorbance of a mixture is a linear combination of the absorbances of its components (
Beer-Lambert law).

The unmixing process involves:

1. Creating reference spectra for each component
2. Establishing calibration curves relating absorbance to concentration
3. Fitting a model to the measured spectrum that determines the concentrations of each component

### Multi-spectrum analysis

For complex mixtures with components that span a wide range of concentrations, a single measurement may not capture all
components accurately. The `multispectrum_to_concentration` method allows the use of multiple spectra of the same sample
at different dilutions to accurately determine all component concentrations.

The main functionality is in the `process_wellplate_spectra.py` and `calibrator.py` modules, which provide:

1. **Data loading**: Load spectral data from various sources (CRAIC microspectrometer, NanoDrop spectrophotometer).
2. **Spectral processing**: Apply corrections, filtering, and masking to raw spectral data.
3. **Calibration**: Create reference spectra and calibration curves for known compounds.
4. **Spectral unmixing**: Determine component concentrations in complex mixtures based on absorption spectra.

## Data sources

The system supports two main data sources:

- **CRAIC Microspectrometer**: Developed for microspectroscopy. High price (>$200k). Poor performance in UV.
- **NanoDrop Spectrophotometer**: Developed (and commonly used) for UV-Vis spectroscopy of biological samples. Lower
  price (~$10k). High performance in UV.

## Calibration data dile Format

For tutorial of a worked example of the calibration workflow, 
see the Jupyter Notebook `notebooks/spectral_unmixing_tutorial.ipynb` notebook.
The calibration workflow requires two input files with the same base filename but different extensions:

### Spectral data file (.csv)

Contains UV-Vis absorption spectra for all calibration samples.

**Structure:**

- First column (unnamed): Wavelength values in nanometers (e.g., 220, 221, 222, ...). Note that in the specific case of 
 NanoDrop spectrophotometer, the output spectrum may contain wavelengths from 190 nm to 220 nm, but sometimes it does not:
 it varies from sample to sample. The guaranteed wavelength range is from 220 nm upwards. For consistent workflow, 
 the wavelengths below 220 nm are removed upon loading the data from file.

- Subsequent columns (0, 1, 2, 3, ...): Absorbance values for each sample (in absorbance units)
- Each row represents one wavelength point
- Each numbered column represents one sample measurement

**Example:**

```csv
,0,1,2,3,4,5,6,7,8,9,10,11
190,0,0.421,0.705,0.817,0.656,0.868,0.855,0.484,0.248,0.293,0.221,0
191,0,1.274,1.556,1.704,1.576,1.744,1.677,1.247,0.898,0.891,0.858,0
192,0,2.406,2.57,2.637,2.572,2.618,2.565,2.132,1.784,1.759,1.747,0
```

### Metadata file (.txt)

Contains sample information that links CSV column numbers to chemical identities and concentrations.

**Required columns:**

- `nanodrop_col_name`: Column number in the CSV file (0, 1, 2, ...)
- `substance`: Chemical substance name (must match `calibrant_shortnames`)
- `concentration`: Concentration value in mol/L
- `solvent`: Solvent used (for reference)

Optional columns:

- `dilution_factor`: Dilution factor applied. Its value is not currently used in the code.

**Example:**

```
nanodrop_col_name,substance,dilution_factor,solvent,concentration
0,solvent_Ethanol,1,Ethanol,0
1,methoxybenzaldehyde,1,Ethanol,0.0000025
2,methoxybenzaldehyde,1,Ethanol,0.000005
3,methoxybenzaldehyde,1,Ethanol,0.00001
4,methoxybenzaldehyde,1,Ethanol,0.000025
5,ethyl_acetoacetate,1,Ethanol,0.020
6,ethyl_acetoacetate,1,Ethanol,0.016
```

### Requirements for calibration

For each substance you want to calibrate:

1. **Background sample**: At least one sample with `concentration=0` (solvent only)
2. **Reference concentration**: One sample at the specified reference concentration
3. **Concentration series**: Multiple samples spanning the concentration range of interest
4. **Consistent naming**: Substance names in the metadata must match the `calibrant_shortnames` parameter

### Typical calibration workflow

1. Prepare standard solutions of known concentrations
2. Measure UV-Vis spectra using NanoDrop spectrophotometer
3. Export data to CSV format with wavelengths and absorbances
4. Create metadata file linking column numbers to sample information
5. Run `perform_calibration()` with appropriate parameters
6. Use generated calibration files for spectral unmixing of unknown samples