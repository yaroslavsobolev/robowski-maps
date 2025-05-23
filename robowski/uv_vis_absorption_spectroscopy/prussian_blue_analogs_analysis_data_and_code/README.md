# Instructions for data and script used in PBA project

## data folder contains all spectra datas
"PBA_reference_calibration_spectra.csv" file contains solvent background spectra; calibration spectra for 3 references, each have 5 different concentrations.
"PBA_sample_spectra.csv" file contains all the sample UV spectra (220-600 nm), each sample spectrum is assigned with specific index.

## result folder contains the fiting results
"Dataset_PBA_composition.csv" file contains the chemical composition of PBAs used for each sample based on their sample index.
"Spectra_unmixing_result.csv" file contains the fitting result for all samples, including the concentration of 3 reference, yield of styrene oxide, yield of benzaldehyde for each sample.
Image files are the comparison between fitting spectrum and experimental spectrum for each sample. Each image contains the results for 20 samples sorted by their sample index.

## "UV_spectra_unmixing_PBA.ipynb" is the demo script used for spectra unmixing for PBA part.

