# Summary of reaction yield data

This folder contains tables in comma-delimited (CSV) format with the yields of reactions at different conditions - one condition per row.
Yields are with respect to the least abundant reactant. Concentrations are given in mol/L,
unless otherwise specified. Concentrations of substances before the reaction are given in the column starting with
`c#`: for example, `c#HBr`. Concentrations of substances measured in the crude reaction mixture after the reaction
are given in the column starting with `pc#` followed by the substance identifier: for example, `pc#E1OH02`.
Volumes of stock solution pipetted into the reaction mixture prior to the reaction are given in the column starting with
`vol#` followed by the substance identifier: for example, `vol#Acetonitrile`

In the case of Hantzsch reaction, the concentrations of components as measured by HPLC after the reaction are named like
`19d_c_HPLC_[M]`, where `19d` is the substance code (same as in the research article, consult Figure 5 there),
`c_HPLC_[M]` means concentration by HPLC in mol/L. The columns containing the yields are named similarly, but with `_y_`
instead of `_c_`. For example: `19d_y_HPLC_[%]`.

## E1 reaction

The file `E1/raw_yields.csv` contains the raw yields of the E1 reaction (scheme in the Figure 2a of the research paper).

## First SN1 reaction 

The file `SN1_simple/raw_yields.csv` contains the data for the SN1 reaction shown in Figure 2b of the research paper.

## Second SN1 reaction

THe file `SN1/raw_yields.csv` contains the data for the SN1 reactions shown in Figure 3a of the research paper.

## Ugi reaction

The file `Ugi/raw_yields.csv` contains the raw yields of the Ugi reaction at different conditions.

The file `Ugi/postprocessed_yields.csv` contains the yields after filtering out the outliers, performing smoothing
and interpolation for more intermediate values of p-TSA. See the Supplementary Materials Section 3.9 for
more details.

The file `Ugi/raw_yields_together_with_model_predictions_2025-01-31.csv` contains the raw yields of the Ugi reaction at different conditions
together with the model predictions. The model predictions are based on the kinetic model of the Ugi reaction described
in the Supplementary Materials Section 5.3.

## Hantzsch reaction

The files `Hantzsch/Hantzsch_26deg_results_HPLC.csv` and `Hantzsch/Hantzsch_80deg_results_HPLC.csv` contain
the Hantzsch reaction yield data for 26 and 80 degrees Celcius, respectively. 
Hantzsch reaction scheme is shown Figure 5 in the research article.

The columns containing concentrations of components as measured by HPLC after the reaction are named
like `19d_c_HPLC_[M]`, where `19d` is the substance code (same as in the research article, consult Figure 6 there),
`c_HPLC_[M]` means concentration by HPLC in mol/L. 

The columns containing the yields are named similarly, but with `_y_` instead
of `_c_`. For example: `19d_y_HPLC_[%]`.

## Reaction with Prussian Blue Analogs (PBA)

- Table `PBA/Dataset_PBA_composition.csv`  
  Contents:
  - Sample index (column `index`).
  - Ratio of different metals in PBA used for the sample. Column `Mn` to column `Cu` stand for the ratio on metal 
  B site, which add up to 1; column `M1-Fe` and column `M1-Co` stand for the ratio on metal A site, which also add up to 1.


- Table `PBA/Spectra_unmixing_result.csv`  
  Contents:
  - Column `Index` is the sample's index.
  - Column `ref1` to column `ref3`: fitted concentration values for the 3 reference compounds for each sample.
  - Column `yield_epox` and column `yield_alde`: yield of styrene oxide and benzaldehyde for each sample.

# Preprocessing the datasets for the Hyperspace Viewer

## 3D datasets

The preprocessing (converting the 3D data to the format suitable for HyperspaceViewer) is done by
the `convert_dataset_to_3dviewer_hyvu_format.py` script. The commands we used to do this preprocessing were
executed in this folder (`summary_of_reaction_yield_data`) and are:

```bash
python convert_3d_dataset_to_hyvu_format.py --input_csv=SN1/raw_yields.csv --output_csv=SN1/SN1_raw_yields_hyvu.csv --X='c#SN1OH01' --Y='c#HBr' --Z='Temperature' --V='yield' --Xscale=1000 --Yscale=1000 --Zscale=1 --Xrename='[SN10H01](mM)' --Yrename='[HBr](mM)'
python convert_3d_dataset_to_hyvu_format.py --input_csv=SN1/resampled_SN1_yield.csv --output_csv=SN1/resampled_SN1_yield_hyvu.csv --X='Alcohol(mM)' --Y='HBr(mM)' --Z='Temperature(°C)' --V='Yield' --Xscale=1 --Yscale=1 --Zscale=1 --Xrename='[Alcohol](mM)' --Yrename='[HBr](mM)'
python convert_3d_dataset_to_hyvu_format.py --input_csv=SN1/resampled_SN1_yield_15d.csv --output_csv=SN1/resampled_SN1_yield_15d_hyvu.csv --X='Alcohol(mM)' --Y='HBr(mM)' --Z='Temperature(°C)' --V='yield of 15d' --Xscale=1 --Yscale=1 --Zscale=1 --Xrename='[Alcohol](mM)' --Yrename='[HBr](mM)'

python convert_3d_dataset_to_hyvu_format.py --input_csv=SN1_simple/raw_yields.csv --output_csv=SN1_simple/simpleSN1_raw_yields_hyvu.csv --X='c#SN1OH03' --Y='temperature' --Z='c#HBr' --V='yield' --Xscale=1000 --Yscale=1 --Zscale=1000 --Xrename='[Alcohol](mM)' --Yrename='Temperature(°C)' --Zrename='[HBr](mM)'
python convert_3d_dataset_to_hyvu_format.py --input_csv=SN1_simple/simpleSN1_resampled_yields.csv --output_csv=SN1_simple/simpleSN1_resampled_yields_hyvu.csv --X='Alcohol(mM)' --Y='Temperature(°C)' --Z='HBr(mM)' --V='yield' --Xscale=1 --Yscale=1 --Zscale=1 --Xrename='[Alcohol](mM)' --Yrename='Temperature(°C)' --Zrename='[HBr](mM)'

python convert_3d_dataset_to_hyvu_format.py --input_csv=E1/E1_resampled_yields.csv --output_csv=E1/E1_resampled_yields_hyvu.csv --X='Alcohol(mM)' --Y='Temperature(°C)' --Z='HBr(mM)' --V='yield' --Xscale=1 --Yscale=1 --Zscale=1 --Xrename='[Alcohol](mM)' --Yrename='Temperature(°C)' --Zrename='[HBr](mM)'
```

## 4D datasets

The Ugi reaction yields were converted to the format suitable for HyperspaceViewer using the `convert4d.py` script:

```bash
python convert_4d_dataset_to_hyvu_format.py --input_csv=Ugi/postprocessed_yields.csv --output_name=Ugi/ugi_postprocessed_yields_hyvu/ugi_hyvu --prefix_of_relative_dir='../CSV/ugi_postprocessed_yields_hyvu/' --X='am001' --Y='ald001' --Z='ic001' --T='ptsa' --V='yield' --Xrename='[amine](mM)' --Yrename='[aldehyde](mM)' --Zrename='[isocyanide](mM)' --Xscale=1000 --Yscale=1000 --Zscale=1000 --Tlabel_prefix='[pTSA]=' --Tlabel_suffix=' M'
```

