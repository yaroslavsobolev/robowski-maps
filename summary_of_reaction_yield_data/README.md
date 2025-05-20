# Summary of reaction yield data

This folder contains tables in comma-delimited (CSV) format with the yields of reactions at different conditions - one condition per row.
Yields are with respect to the least abundant reactant. Concentrations are given in mol/L,
unless otherwise specified. Concentrations of substances before the reaction are given in the column starting with
`c#`: for example, `c#HBr`. Concentrations of substances measured in the crude reaction mixture after the reaction
are given in the column starting with `pc#` followed by the substance identifier: for example, `pc#E1OH02`.
Volumes of stock solution pipetted into the reaction mixture prior to the reaction are given in the column starting with
`vol#` followed by the substance identifier: for example, `vol#Acetonitrile`

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
Hantzsch reaction scheme is shown Figure 6 in the research article.

The columns containing concentrations of components as measured by HPLC after the reaction are named
like `19d_c_HPLC_[M]`, where `19d` is the substance code (same as in the research article, consult Figure 6 there),
`c_HPLC_[M]` means concentration by HPLC in mol/L. 

The columns containing the yields are named similarly, but with `_y_` instead
of `_c_`. For example: `19d_y_HPLC_[%]`.

# Preprocessing the datasets for the interactive viewer (YipMan)

## 3D datasets

The preprocessing (converting the 3D data to the format suitable for YipMan) is done by
the `convert_dataset_to_3dviewer_yipman_format.py` script. The commands we used to do this preprocessing were
executed in this folder (`summary_of_reaction_yield_data`) and are:

```bash
python convert_3d_dataset_to_yipman_format.py --input_csv=SN1/raw_yields.csv --output_csv=SN1/SN1_raw_yields_yipman.csv --X='c#SN1OH01' --Y='c#HBr' --Z='Temperature' --V='yield' --Xscale=1000 --Yscale=1000 --Zscale=1 --Xrename='[SN10H01](mM)' --Yrename='[HBr](mM)'
python convert_3d_dataset_to_yipman_format.py --input_csv=SN1/resampled_SN1_yield.csv --output_csv=SN1/resampled_SN1_yield_yipman.csv --X='Alcohol(mM)' --Y='HBr(mM)' --Z='Temperature(°C)' --V='Yield' --Xscale=1 --Yscale=1 --Zscale=1 --Xrename='[Alcohol](mM)' --Yrename='[HBr](mM)'
python convert_3d_dataset_to_yipman_format.py --input_csv=SN1/resampled_SN1_yield_15d.csv --output_csv=SN1/resampled_SN1_yield_15d_yipman.csv --X='Alcohol(mM)' --Y='HBr(mM)' --Z='Temperature(°C)' --V='yield of 15d' --Xscale=1 --Yscale=1 --Zscale=1 --Xrename='[Alcohol](mM)' --Yrename='[HBr](mM)'
```

## 4D datasets

The Ugi reaction yields were converted to the format suitable for YipMan using the `convert4d.py` script:

```bash
python convert_4d_dataset_to_yipman_format.py --input_csv=Ugi/postprocessed_yields.csv --output_name=Ugi/ugi_postprocessed_yields_yipman/ugi_yipman --prefix_of_relative_dir='../CSV/ugi_postprocessed_yields_yipman/' --X='am001' --Y='ald001' --Z='ic001' --T='ptsa' --V='yield' --Xrename='[amine](mM)' --Yrename='[aldehyde](mM)' --Zrename='[isocyanide](mM)' --Xscale=1000 --Yscale=1000 --Zscale=1000 --Tlabel_prefix='[pTSA]=' --Tlabel_suffix=' M'
```