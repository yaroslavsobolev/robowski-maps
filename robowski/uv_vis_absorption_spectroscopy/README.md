# UV-Vis Absorption Spectroscopy Analysis

This module provides tools for processing and analyzing UV-Vis absorption spectroscopy data from well plates, with a
focus on spectral unmixing to determine component concentrations in multi-component reaction mixtures.

## Overview

Spectral unmixing is a mathematical process to decompose a mixed spectrum into its constituent components. The core
assumption is that the absorbance of a mixture is a linear combination of the absorbances of its components (
Beer-Lambert law).

The unmixing process involves:

1. Creating reference spectra for each component
2. Establishing calibration curves relating absorbance to concentration
3. Fitting a model to the measured spectrum that determines the concentrations of each component

### Multi-spectrum Analysis

For complex mixtures with components that span a wide range of concentrations, a single measurement may not capture all
components accurately. The `multispectrum_to_concentration` method allows the use of multiple spectra of the same sample
at different dilutions to accurately determine all component concentrations.

The main functionality is in the `process_wellplate_spectra.py` and `calibrator.py` modules, which provide:

1. **Data Loading**: Load spectral data from various sources (CRAIC microspectrometer, NanoDrop spectrophotometer).
2. **Spectral Processing**: Apply corrections, filtering, and masking to raw spectral data.
3. **Calibration**: Create reference spectra and calibration curves for known compounds.
4. **Spectral Unmixing**: Determine component concentrations in complex mixtures based on absorption spectra.

## Data Sources

The system supports two main data sources:

- **CRAIC Microspectrometer**: Developed for microspectroscopy. High price (>$200k). Poor performance in UV.
- **NanoDrop Spectrophotometer**: Developed (and commonly used) for UV-Vis spectroscopy of biological samples. Lower
  price (~$10k). High performance in UV.

## Calibration Data File Format

The calibration workflow requires two input files with the same base filename but different extensions:

### Spectral Data File (.csv)

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

### Metadata File (.txt)

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

### Relationship between metadata and spectral data

The connection between files works as follows:

- CSV column "0" → TXT row with `nanodrop_col_name=0` (background/solvent)
- CSV column "1" → TXT row with `nanodrop_col_name=1` (methoxybenzaldehyde at 0.0000025 M)
- CSV column "2" → TXT row with `nanodrop_col_name=2` (methoxybenzaldehyde at 0.000005 M)
- And so on...

### Requirements for Calibration

For each substance you want to calibrate:

1. **Background sample**: At least one sample with `concentration=0` (solvent only)
2. **Reference concentration**: One sample at the specified reference concentration
3. **Concentration series**: Multiple samples spanning the concentration range of interest
4. **Consistent naming**: Substance names in the metadata must match the `calibrant_shortnames` parameter

### Typical Workflow

1. Prepare standard solutions of known concentrations
2. Measure UV-Vis spectra using NanoDrop spectrophotometer
3. Export data to CSV format with wavelengths and absorbances
4. Create metadata file linking column numbers to sample information
5. Run `perform_calibration()` with appropriate parameters
6. Use generated calibration files for spectral unmixing of unknown samples

## Main Workflow

### 1. Create a SpectraProcessor instance

```python
from robowski.uv_vis_absorption_spectroscopy.process_wellplate_spectra import SpectraProcessor

sp = SpectraProcessor(
    folder_with_correction_dataset='path/to/correction/dataset',
    spectrum_data_type='craic'  # or 'nanodrop'
)
```

### 2. Create Calibration Reference Spectra

Before analyzing unknown samples, you need to create reference spectra for each component:

```python
sp.construct_reference_for_calibrant(
    calibrant_shortname='component_name',
    calibration_folder='path/to/calibration/folder',
    ref_concentration=0.001,  # Reference concentration in mol/L
    do_plot=True,
    do_reference_refinements=True
)
```

This needs to be done once for each component you want to detect.

### 3. Analyze Unknown Samples

For a single spectrum:

```python
concentrations = sp.spectrum_to_concentration(
    target_spectrum_input=spectrum,
    calibration_folder='path/to/calibration/folder',
    calibrant_shortnames=['component1', 'component2', ...],
    background_model_folder='path/to/background/model',
    do_plot=False
)
```

For analyzing multiple spectra with different dilutions:

```python
concentrations = sp.multispectrum_to_concentration(
    target_spectrum_inputs=[spectrum1, spectrum2],
    dilution_factors=[1, 10],  # Dilution factors for each spectrum
    calibration_folder='path/to/calibration/folder',
    calibrant_shortnames=['component1', 'component2', ...],
    background_model_folder='path/to/background/model',
    do_plot=False
)
```

### 4. Process All Wells in a Plate

To analyze all wells in a plate:

```python
all_concentrations = sp.concentrations_for_one_plate(
    experiment_folder='path/to/experiment',
    plate_folder='path/to/plate/data',
    calibration_folder='path/to/calibration/folder',
    calibrant_shortnames=['component1', 'component2', ...],
    calibrant_upper_bounds=[np.inf, np.inf, ...],
    background_model_folder='path/to/background/model',
    do_plot=False,
    return_all_substances=True
)
```

For plates measured with multiple dilutions:

```python
all_concentrations = sp.multispectrum_concentrations_for_one_plate(
    experiment_folder='path/to/experiment',
    plate_folder_1='path/to/plate1/data',
    plate_folder_2='path/to/plate2/data',
    dilution_factors=[1, 10],
    calibration_folder='path/to/calibration/folder',
    calibrant_shortnames=['component1', 'component2', ...],
    calibrant_upper_bounds=[np.inf, np.inf, ...],
    background_model_folder='path/to/background/model',
    do_plot=False,
    return_all_substances=True,
    return_report=True
)
```

## Advanced Features

### Stoichiometric Constraints

The system can enforce stoichiometric constraints during spectral unmixing:

```python
concentrations = sp.multispectrum_to_concentration(
    # ... other parameters ...
    starting_concentration_dict={'substrate1': 0.1, 'substrate2': 0.2, ...},
    obey_stoichiometric_inequalities=True
)
```

This ensures that the calculated product concentrations don't exceed what's possible given the starting substrate
concentrations.

### Uncertainty Estimation

The system can estimate uncertainties in the calculated concentrations:

```python
concentrations, errors = sp.multispectrum_to_concentration(
    # ... other parameters ...
    return_errors=True
)
```

### Detailed Fitting Reports

For in-depth analysis of the fitting process:

```python
concentrations, report = sp.multispectrum_to_concentration(
    # ... other parameters ...
    return_report=True
)
```

The report contains information about fitting quality, residuals, and other metrics.

## Example Use Cases

- **Reaction Yield Analysis**: Determine product yields in complex reaction mixtures
- **Kinetics Monitoring**: Track changes in component concentrations over time
- **Quality Control**: Ensure consistent composition of products
- **High-Throughput Screening**: Analyze large numbers of reaction conditions efficiently
