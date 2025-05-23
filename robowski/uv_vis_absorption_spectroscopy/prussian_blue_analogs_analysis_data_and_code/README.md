# Description of code and data files for the Prussian Blue Analogs (PBA) analysis

## 📁 `data/` — Experimental data

- Table `PBA_reference_calibration_spectra.csv`  
  Contains:
  - Spectra wavelength in nanometers (column `Wavelength`).
  - Solvent background spectrum (column `bkg`: ).
  - Calibration spectra for 3 references (each at 5 concentrations). Column names refer to the reference type and their 
concentration, for example, `ref2_c2` column contains the calibration spectrum for reference substance #2 of its second concentration
  - (see concentration details in the script).

- `PBA_sample_spectra.csv`  
  Contents:
  - Spectra wavelength in nanometers (Column `Wavelength`).
  - Each column corresponds spectrum to a sample, column name is a unique index assigned to each sample (from 0 to 755).

---

## 📁 `result/` — Results of the analysis

- Table `Dataset_PBA_composition.csv`  
  Contents:
  - Sample index (column `index`).
  - Ratio of different metals in PBA used for the sample. Column `Mn` to column `Cu` stand for the ratio on metal 
  B site, which add up to 1; column `M1-Fe` and column `M1-Co` stand for the ratio on metal A site, which also add up to 1.

- Table `Spectra_unmixing_result.csv`  
  Contents:
  - Column `Index` is the sample's index.
  - Column `ref1` to column `ref3`: fitted concentration values for the 3 reference compounds for each sample.
  - Column `yield_epox` and column `yield_alde`: yield of styrene oxide and benzaldehyde for each sample.

- `*.png` images  
  - Each image shows the comparison between experimental and fitted UV spectra.
  - Each image includes results for 20 samples, sorted by sample index.

---

## 📜 `UV_spectra_unmixing_PBA.ipynb` — Jupyter Notebook used for the analysis

This Jupyter Notebook demonstrates the full workflow of:
- Spectra preprocessing
- SVR-based spectral unmixing
- Plotting fits and residuals
- Exporting and interpreting the results
