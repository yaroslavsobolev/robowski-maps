# Plotting the multidimensional experimental yield maps for various reactions

## Scripts used to produce the yield map figures in the research paper:

Table of the reaction, figure in the research article, and the corresponding script:

| Reaction        | Figure    | Script                                  |
|-----------------|-----------|-----------------------------------------|
| E1              | Figure 2a | `plot_3d_map_for_E1_reaction.py`        |
| Simple SN1      | Figure 2b | `plot_3d_map_for_simpleSN1_reaction.py` |
| Less simple SN1 | Figure 3  | `plot_3d_map_for_SN1_reaction.py`       |
| Ugi             | Figure 4  | `view_4d_data_for_ugi_reaction.py`      |

## Auxiliary scripts

### Script `controls_of_unmixing_precision_for_hantzsch_reaction.py`

This script is used to evaluate the precision of determining the concentrations in the
Hantzsch reaction by comparing the results obtained from different methods:
 1. UV-VIS (by spectral unmixing)
 2. NMR (by standard method involving the integration of the peaks in the spectrum)

The conclusion is that the UV-VIS method is unreliable for determining the concentrations of products in this
reaction. The results for substances other than Hantzsch ester and product 19e (in the notation of the research article)
are entirely hopeless. The results for Hantzsch ester and piperidone don't look so hopeless on this plot, but subsequent
comparisons (UV-VIS vs. NMR) we performed at other conditions show that the concentrations by the UV-VIS method are
not reliable at all for the Hantzsch reaction: too many components in the unmixing model make it unstable,
causing large uncertainties in the concentration values.

### Script `plot_3d_map_for_hantzsch_reaction_from_uv_spectra.py`

This script is used to visualize the results of a Hantzsch reaction based on the spectral unmixing of UV spectra.
These concentration values are extremely unreliable: too many components in the unmixing model make it unstable,
causing large uncertainties in the concentration values. Therefore, these plots of concentrations were never presented
in the accompanying research article. You can plot it out of curiosity, but don't expect it to be useful for anything.
The maps of yield plotted in Figure 5 of the article are based on the HPLC measurements, not on unmixing of UV-Vis
spectra.

### Script `plot_onion_map_for_ugi_reaction.py`

This script is an implementation of a visualization approach suggested by one of the Reviewers
of the research article:
> (four reactants and one product): use three orthogonal axes for isocyanide, amine, and aldehyde concentrations; 
> indicate the concentration of p-toluenesulfonic acid by node size; and yield by node color.

This visualization method is not used in the research article, but we leave it here for 
entertainment of a curious reader.