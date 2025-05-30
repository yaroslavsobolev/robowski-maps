# Robot-assisted reconstruction of reaction hyperspaces and complex reaction networks they contain

This repository contains the code and data supporting the article

'Robot-assisted reconstruction of reaction hyperspaces and complex reaction networks they contain'
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
will be used by the Python scripts to locate the data files. On Windows, you may have to restart the computer after
changing the environment variable for it to be recognized by the system.

## Requirements for running the Python scripts

For requirements of the automation script used for pipetting robots coupled to the NanoDrop spectrometer,
see the `requirements.txt` file in the `zeus-pipetter` folder and the `README.md` file in the same folder.

For the rest of the repository, the installation procedure is as follows:

### install Microsoft Visual C++ 14.0 or greater (required only for [Mayavi](https://docs.enthought.com/mayavi/mayavi/)).
After downloading `vs_buildtools.exe` from [Microsoft website](https://visualstudio.microsoft.com/visual-cpp-build-tools/): 
When installing, in the section "Desktop & Mobile" click the "Desktop development with C++" option
in the opened window, then click "Install" button in the lower right corner. 

### Creating Anaconda environment with all the required dependencies

The recommended way to install the dependencies is to 
use [Anaconda package manager](https://www.anaconda.com/docs/getting-started/anaconda/install) (version 23.10.0 or greater is recommended).
Launch the terminal with Anaconda activated. The easiest way to do it is by clicking "Start" menu,
finding (e.g. by typing in search) "Anaconda Command Prompt" and running it. It is recommended
that you run it with administrator privileges, so right-click on the "Anaconda Command Prompt" and select "Run as administrator".
In the launched command prompt, navigate to your local folder containing this 
repository and create a new `conda` environment named `robowski` by running this command:
```bash
conda env create -f environment.yml
```

If this command has failed for some reason, you can either attempt the troubleshooting (next section),
or install the environment without Mayavi by running the following command instead:
```bash
conda env create -f environment_without_mayavi.yml
```
The only consequence of not installing Mayavi is that you will not be able to visualize the 3D yield maps of the reactions
using Python scripts. This is not a big issue, because you can instead use the Hyperspace Viewer
distributed with the article (in [Zenodo repository](https://zenodo.org/uploads/14880579) and article's Supplementary 
Files). Hyperspace Viewer is a standalone application that does not even require Python. 
The rest of the code will run normally without Mayavi.

After installation, activate the new environment by the following command:
```bash
conda activate robowski
```
and proceed to the section "Installing this package with `pip`" below.

### Troubleshooting the install

Possible error: during running the above command you get the error message `Intel MKL FATAL ERROR: Cannot load mkl_intel_thread.dll`,
the [official discussion](https://docs.conda.io/projects/conda/en/stable/user-guide/troubleshooting.html#numpy-mkl-library-load-failed) 
(by Anaconda developers) of the causes of this error reads:
> Another software vendor has installed MKL or Intel OpenMP (libiomp5md.dll) files into the C:\Windows\System32 folder. 
> These files are being loaded before Anaconda's and they're not compatible.

Use the [guide about this error on the Anaconda website](https://docs.conda.io/projects/conda/en/stable/user-guide/troubleshooting.html#numpy-mkl-library-load-failed) for troubleshooting,
then retry the `conda env create -f environment.yml` commmand again.

If you get errors saying that the Visual Studio C++ files (e.g. compiler) are not found,
click "Start", type "Developer" and launch the "Developer Command Prompt for VS 2022". Then navigate to 

```
cd C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build
```

or the respective directory of your Visual Studio installation. Run the command

```
vcvarsall.bat amd64_arm
```


### Alternative installation with `pip` instead of Anaconda (not recommended; not tested)

For the alternative method of installation, we provide the `requirements.txt` list of dependencies that you
may try to use to install these dependencies into one of your existing environments by 
```bash
pip install -r requirements.txt
```
This will install all the required dependencies, but it is not guaranteed to work as expected.

### Installing this package with `pip`

Once the `robowski` environment is activated in `conda`, run the following command in the terminal ("Anaconda Command Prompt")
to install this repository's code as "in-place" package (don't miss the dot `.` at the end):

```bash
python -m pip install --no-build-isolation --no-deps -e .
```


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
vial plates are in the `robowski/uv_vis_absorption_spectroscopy/craic-automation` folder.
They are mostly based on PyAutoGUI for controlling the GUI of the CRAIC software for moving the XY stage and
starting the measurements. The barcodes on the plates are read using a webcam and OpenCV (`cv2`) barcode detection
functionality.

### Spectral unmixing of the UV-Vis absorption spectra

Scripts related to spectral unmixing are in the `robowski/uv_vis_absorption_spectroscopy` folder.

For familiarizing oneself the overall workflow, we recommend starting with the tutorial Jupyter Notebook
located at `notebooks/spectral_unmixing_tutorial.ipynb` that demonstrates the calibration procedure and spectral unmixing
for the real data of Claisen-Schmidt condensation reaction and reproduces the Figure 1c,d from the article.
In the second chapter of the notebook, the user is guided through the process of performing the same calibration
and spectral unmixing using the high-level methods from this repository.

Processing for all the reactions we studied is implemented in the respective scripts in 
the `robowski/uv_vis_absorption_spectroscopy/examples`
folder. The calibration script (name begins with `calibration_`) should be executed first, and
the spectral unmixing script (name begins with `process_`) should be executed after that. Unmixing is performed for all 
the vial plates obtained in a given experimental run. 

The processing results are saved into a subdirectory `results/` of the run directory,
specifically as the `results/product_concentration.csv` table relative to the directory containing the raw data for the given 
experimental run.

### 3D visualization of the yield maps

The `visualize_results` folder contains the scripts for processing the output of the spectral unmixing workflow.
This includes scripts for visualizing the reaction yield data in 3D,
as well as the scripts for fitting the kinetics models to the data.

Interactive viewer of 3D datasets is implemented in `visualize_results/animated_viewer_static.py` script and is
based on `mayavi` library. It allows the user to rotate the 3D plot, zoom in and out, and offers various
options for visualizing the data.

For interactive visualization of 3D and 4D datasets, see the scripts whose names starting with `plot_3d_map_for_` 
and `view_4d_data_for_` in 
the `robowski/visualize_results/examples` folder.

### Outlier filtration and interpolation in Ugi reaction data

Ater the initial mapping of the Ugi reaction data was completed, outliers were marked for
repeat measurements. Marking the outliers is implemented in `robowski/outlier_handling/manually_mark_outliers.py`.
The new experimental runs were then requested for remeasuring the outliers -- this is implemented in the
`robowski/outlier_handling/request_new_runs_to_rerun_outliers.py` script. After experiments were completed,
the completed runs were processed and the remeasured outlier data were merged with the original data. 
The script for processing and merging the outlier data is
`robowski/uv_vis_absorption_spectroscopy/examples/Ugi_reaction/process_ugi_remeasured_outliers_2023-07-05_run01.py`.


For reproducing the Figures S115-S127 from the Supplementary Materials, run the
`robowski/visualize_results/smooth_results.py` script. This script also superimposes the precalculated yields predicted by the
kinetic model of Ugi reaction (Supplementary Materials Section 5.3) over the experimental data points.

### Fitting the kinetic models to the data

Kinetic fitting for the E1 and SN1 reactions is implemented in the `robowski/kinetic_models/` folder. 

#### E1 reaction
To reproduce
the "corner plot" from precalculated bootstrapping analysis data, run the `robowski/kinetics_models/E1_kinetics/kinetics_of_E1_2023_bootstrap.py` script.
This would reproduce the Figure S141 from the Supplementary Materials.
For regenerating the data for the bootstrap analysis from scratch, 
run the `robowski/kinetics_models/E1_kinetics/E1_kinetics_v3.py` script.
The same script can be used to plot the model prediction against the experimental data (Figures S138 and S140
from the Supplementary Materials) by uncommenting the appropriate lines in the script.

#### SN1 reactions

For kinetics fitting to the SN1 reaction (Figure 2), 
run the `robowski/kinetics_models/simpleSN1_kinetics/kinetics_of_simpleSN1_2023-11_2023-12.py` script.

#### Ugi reaction

The scripts for kinetics fitting to the Ugi reaction (Figure 4 and Section 5.3 of the Supplementary Materials)
are located in the `robowski/kinetics_models` folder. The names of these scripts begin with
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

#### Hantzsch reaction
1. To reproduce the Figure 5a and Figure S99, run the `robowski/misc_scripts/hntz_discovery_timeline.py` script.
2. To reproduce the main-text figure 5b, run the `robowski/misc_scripts/figures_for_articles/all_hantzsch_spectra_figure5b.py` script.
3. To reproduce the Figure 5d, 5e and 5f, run the `robowski/visualize_results/examples/plot_3d_map_for_hantzsch_by_HPLC.py` script.

#### Reactions with prussian blue analogs

The code for spectral unmixing of the UV-Vis absorption spectra of the crude mixture for reactions with prussian 
blue analogs is organized as a Jupyter notebook located at 
the `robowski/uv_vis_absorption_spectroscopy/prussian_blue_analogs_analysis_data_and_code/UV_spectra_unmixing_PBA.ipynb`
file. The respective experimental data is located at 
the `robowski/uv_vis_absorption_spectroscopy/prussian_blue_analogs_analysis_data_and_code/data`.


### Calculating proton exchange equilibrium in a complex mixture

The module `robowski/kinetics_models/acid_base_equilibrium.py` contains the script for calculating
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
`robowski/smoothness_theory` folder. The `smoothness_analyze.py` can be used to reproduce Figures S154 and S155.
To reproduce the Figure S153 showing the smoothness of experimental yield maps, run the
`robowski/smoothness_theory/smoothness_from_experiments.py` script. 

Implementation of the simulation of
the reaction network is implemented in the `robowski/smoothness_theory/smoothness_v2.py` script, and for running
the calculations on the multi-CPU server the `robowski/smoothness_theory/run_simulation.py` script is used,
followed by `robowski/smoothness_theory/construct_short_db.py` script for constructing the concise
database of the simulation results. Be warned that the simulartions are computationally expensive and may take
a long time to run, depending on the number of CPU cores available on your machine.


### Calculation of theoretical light absorption spectra from TD-DFT

The respective code and data (Gaussian output files) are in the `robowski/misc_scripts/dft_molar-absorptivity-spectra` folder.
For description, see Section 6.2 of the Supplementary Materials.
For calculating the spectra and reproducing the Figures S149, S150, S151 run the 
`robowski/misc_scripts/dft_molar-absorptivity-spectra/dft-spectrum-confidence-intervals.py` script.
Gaussian job files and the respective output files are in the `robowski/misc_scripts/dft_molar-absorptivity-spectra/dft-calculation-results`
folder.


### Reproducing other figures from the accompanying research article

To reproduce the main-text Figure 4j, run the `robowski/visualize_results/visualize_results.py` script.
To reproduce the Figure S17 and recalculating the uncertainty map of the NanoDrop spectrophotometer,
run the `robowski/uv-vis-absorption-spectroscopy/absorbance_errorbar_model.py` script.



