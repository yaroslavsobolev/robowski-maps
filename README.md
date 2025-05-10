# Code and data supporting the article 'Maxima, anomalies and networks in chemical reaction hyperspaces'

This repository contains the code and data supporting the article 
'Maxima, anomalies and networks in chemical reaction hyperspaces'
by Yankai Jia, Rafał Frydrych, Yaroslav I. Sobolev, Wai-Shing Wong, Bibek Prajapati, Yasemin Bilgi, 
Daniel Matuszczyk, Galymzhan Moldagulov, Juan Carlos Ahumada, Namhun Kim, Eric Larsen, Bartosz A. Grzybowski

To cite this work, use the following citation:

``` 
    BIBTEX citation entry to be added upon publication
```

## If you are only looking for the final datasets of reaction yields

The final reaction yield datasets are tables in CSV format located in the `summary_of_reaction_yield_data` folder.
See the `README.md` file in that folder for more details.

## Obtaining all the raw data

Final yield data and some other datasets are included into this repository. The rest should be downloaded as a ZIP file 
from the [Zenodo repository](https://zenodo.org/uploads/14880579) associated with this work. 
The ZIP file is named `data.zip`, it is about 7 Gb, and 19 Gb when unzipped.
It should be unzipped to any folder of your choosing, and then the path to that folder should be
set as your operating system's ["environment variable"](https://en.wikipedia.org/wiki/Environment_variable) 
named `ROBOCHEM_DATA_PATH`. This environment variable
will be used by the Python scripts to locate the data files.

## Requirements for running the Python scripts

### Setting the Current Working Directory

Most Python scripts in this repository assume that the current working directory (CWD) is the root of the repository.
The exceptions are the scripts in the `zeus-pipetter` folder, and the scripts in the
`visualize-results/examples/kinetics_models` -- these are written to use script's current locations
as CWD, since they are intended to be run on a separate machine or a server for parallel computing. Failure
to set CWD correctly may result in the scripts not finding the necessary data files or folders and throwing a
`FileNotFoundError: No such file or directory:` exception.

### Required packages
For requirements of the automation script used for pipetting robots coupled to the NanoDrop spectrometer,
see the `requirements.txt` file in the `zeus-pipetter` folder and the `README.md` file in the same folder.

For the rest of the repository, the principal requirements are as follows:

`Python 3.7+`

`matplotlib` 3.1.0 or later (if plotting is needed)

`numpy` 1.21.6 or later

`scipy` 1.5.2 or later

`mayavi`4.7.2 or later (for viewing the interactive 3D plots in). 
To work, `mayavi` requires Microsoft Visual C++ 14.0 or greater. If `conda` or `pip` fail to install it automatically,
you can download it from the [Microsoft website](https://visualstudio.microsoft.com/visual-cpp-build-tools/) and
install it manually.

`tqdm` 4.64.1 or later (for progress bars)

`statsmodels` 0.13.5 or later (for Ljung-Box metrics)

`pandas` 1.3.5 or later

Our code was not tested on all combinations of all package versions later than the ones specified above,
but it is likely to work with the latest versions of the libraries.

## Overview of the scripts and functionalities. Reproducing the article's figures.

### Automation of pipetting and NanoDrop spectrophotometer measurement

The `zeus-pipetter` folder contains the scripts for automating the pipetting of the reaction mixtures,
dilution of the final crude mixtures and measuring the UV-Vis absorption spectra with NanoDrop 
spectrometer. See the `README.md` file in that folder for more details.

### Automatic washing of vial plates

The `washing-station` folder contains the scripts for automating the washing of the vial plates.
It interfaces with the FlowEZ pump and the AxiDraw pen plotter that acts as a robotic arm moving the washing
head into and out of each vial, and from one vial to another.

### Automation of CRAIC microspectrometer

The scripts for automating the CRAIC microspectrometer for measuring the UV-Vis absorption spectra of the 
vial plates are in the `uv-vis-absorption-spectroscopy/craic-automation` folder.
They are mostly based on PyAutoGUI for controlling the GUI of the CRAIC software for moving the XY stage and
starting the measurements. The barcodes on the plates are read using a webcam and OpenCV (`cv2`) barcode detection
functionality.

### Spectral unmixing of the UV-Vis absorption spectra

Scripts related to spectral unmixing are in the `uv-vis-absorption-spectroscopy` folder. Most of the functionality
for extracting the component concentrations from the spectra is in the
`process_wellplate_spectra.py` module. This file is supposed to be imported into other scripts that use its functions;
it is not meant to be run as a standalone script. 

For familiarizing oneself the overall workflow, we recommend starting with
the scripts `uv-vis-absorption-spectroscopy/examples/versatility` folder as easy-to-understand, minimal examples of 
spectral processing.
For instance, the `Beckmann_WaiShing.py` is first running the calibration routines `construct_calibrant()` for each of the components
of the reaction mixture, and then executes the `process_plate()` method that iterates over the spectra from separate
vials of the vial plate and applies the unmixing algorithm. The results are returned as `pandas` DataFrame object
containing the concentrations of the components in each vial. Calibration procedures take the spectra of the pure
components and construct the calibration matrix that is used in the unmixing algorithm and the concentration calculation.
Processing for all the other reactions we studied is implemented in the respective scripts in the `uv-vis-absorption-spectroscopy/examples`
folder and follows the same pattern, with the calibration executed first 
(sometimes in a dedicated script usually named like `construct_reference...` or `construct_calibration`), and
the spectral unmixing executed after that for all the vial plates obtained in a given experimental run --
usually in a script whose name begins with `process_...`. The processing results are saved into data folder,
specifically as the `results/product_concentration.csv` table in directory containing the raw data for the given 
experimental run.

### 3D visualization of the yield maps

The `visualize-results` folder contains the scripts for processing the output of the spectral unmixing workflow.
This includes scripts for visualizing the reaction yield data in 3D,
as well as the scripts for fitting the kinetics models to the data.

Interactive viewer of 3D datasets is implemented in `visualize-results/animated_viewer_static.py` script and is
based on `mayavi` library. It allows the user to rotate the 3D plot, zoom in and out, and offers various
options for visualizing the data.

For interactive visualization of 4D datasets, see the scripts whose names starting with `animated-viewer-singlecube0...` in 
the `visualize-results/examples` folder. Some of these scripts also perform fitting of kinetic models to the data
prior to the visualization.

### Outlier filtration and interpolation in Ugi reaction data

For reproducing the Figures S115-S127 from the Supplementary Materials, runs the
`visualize-results/smooth_results.py` script. This script also superimposes the precalculated yields predicted by the
kinetic model of Ugi reaction (Supplementary Materials Section 5.3) over the experimental data points.

### Fitting the kinetic models to the data

Kinetic fitting for the E1 and SN1 reactions is implemented in the `visualize-results/examples/` folder. 

#### E1 reaction
To reproduce
the "corner plot" from precalculated bootstrapping analysis data, run the `visualize-results/examples/kinetics_of_E1_2023_bootstrap.py` script.
This would reproduce the Figure S141 from the Supplementary Materials.
For regenerating the data for the bootstrap analysis from scratch, 
run the `visualize-results/examples/animated-viewer-singlecube-E1-narrowrange_with_perr_v3.py` script.
The same script can be used to plot the model prediction against the experimental data (Figures S138 and S140
from the Supplementary Materials) by uncommenting the appropriate lines in the script.

#### SN1 reactions

For kinetics fitting to the SN1 reaction (Figure 3), 
run the `visualize-results/examples/animated-viewer-singlecube-simpleSN1-2023-08-21.py` script.

#### Ugi reaction

The scripts for kinetics fitting to the Ugi reaction (Figure 4 and Section 5.3 of the Supplementary Materials)
are located in the `visualize-results/examples/kinetics_models` folder. The names of these scripts begin with
`ugi_kinetics_...` followed by the version of the kinetics model and the name of the procedure.
In relation to the description in the Section 5.3 of the Supplementary Materials,
version `v3b` is the "reduced model", `v3` is the full model, and `v3c` is the full model after
the change of variables. The main functionality of each model is contained in the script whose name ends with `_emcee.py`,
whereas scripts whose names end with `_plot.py` and `_optimize.py` are importing the `_emcee.py` script and 
are used for plotting the results of the fitting, or for optimization of the model. The user should pay attention to 
the initial guesses set before running these scripts.
The Markov Chain Monte Carlo (MCMC) sampling is performed by routines in the `_emcee.py` scripts.
Covariance matrices are plotted by the `ugi_kinetics_v3_covplot.py` script based on the precalculated output
of the scripts whose names contain `_perr_`, meaning "parameter errors analysis".


### Calculating proton exchange equilibrium in a complex mixture

The module `visualize-results/examples/kinetics_models/acid_base_equilibrium.py` contains the script for calculating
the proton transfer equilibria in the complex mixtures. The interface is similar to the interface of the 
[`pHcalc` Python package](https://github.com/rnelsonchem/pHcalc). In contrast to [pHcalc](https://github.com/rnelsonchem/pHcalc),
the `acid_base_equilibrium.py` script is more general, as
no assumption is made about the stoichiometry of the acid-base reactions, and the excess of water is not assumed:
water concentration can be comparable to concentrations of other species. The script calculates the equilibrium pH and 
the equilibrium
concentrations of all species in the mixture, given the initial sums of concentrations of both protonated and deprotonated
forms of each species, and the pKa values of the species.
Some examples of the usage of this script are in the `main` function of the `acid_base_equilibrium.py` script itself.


### Numerical exploration of the smoothness of the yield maps of reaction networks

The scripts used as described in the Section 7.5 of the Supplementary Materials are in the
`misc-scripts/smoothness_theory` folder. The `smoothness_analyze.py` can be used to reproduce Figures S154 and S155.
To reproduce the Figure S153 showing the smoothness of experimental yield maps, run the
`misc-scripts/smoothness_theory/smoothness_from_experiments.py` script.


### Calculation of theoretical light absorption spectra from TD-DFT

The respective code and data (Gaussian output files) are in the `misc-scripts/dft_molar-absorptivity-spectra` folder.
For description, see Section 6.2 of the Supplementary Materials.
For calculating the spectra and reproducing the Figures S149, S150, S151 run the 
`misc-scripts/dft_molar-absorptivity-spectra/dft-spectrum-confidence-intervals.py` script.
Gaussian job files and the respective output files are in the `misc-scripts/dft_molar-absorptivity-spectra/dft-calculation-results`
folder.


### Reproducing other figures from the accompanying research article

To reproduce the main-text Figure 5J, run the `visualize-results/visualize_results.py` script.

To reproduce the Figure S17 and recalculating the uncertainty map of the NanoDrop spectrophotometer,
run the `uv-vis-absorption-spectroscopy/absorbance_errorbar_model.py` script.



