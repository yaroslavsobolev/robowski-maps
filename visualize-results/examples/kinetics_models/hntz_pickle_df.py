from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import importlib
import logging

logging.basicConfig(level=logging.ERROR)

NCORES = 77
organize_run_results = importlib.import_module("misc-scripts.organize_run_results")


list_of_runs = tuple(['2024-03-04-run01',
                      '2024-03-04-run02'])

substances = ['c#ethyl_acetoacetate',  'c#methoxybenzaldehyde', 'c#ammonium_acetate']
substance_titles = ['Acetoacetate', 'Methoxy', 'Ammonium acetate']
# substrates = ['c#SN1OH03', 'c#HBr']

df_results = organize_run_results.join_data_from_runs([f'BPRF/{x}/' for x in list_of_runs],
                                 round_on_columns=None)
column_to_plot = 'yield#bb017'
df_results.dropna(subset=[column_to_plot], inplace=True)
df_results = df_results[~df_results[column_to_plot].isin([np.inf, -np.inf])]

# pickle df_results
df_results.to_pickle('hntz_df_results.pkl')
# unpickle df_results
df_results = pd.read_pickle('hntz_df_results.pkl')

# pickle df_results