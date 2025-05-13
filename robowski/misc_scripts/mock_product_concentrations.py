from robowski.settings import *
import os
import numpy as np
import pandas as pd

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

experiment_name = 'multicomp-reactions/2023-06-19-run01/'
concentrations_df = pd.read_csv(data_folder + experiment_name + 'outVandC/' + 'outC.csv')
# add a column for the product concentration, fill it with ones
concentrations_df['IIO029A'] = concentrations_df['ptsa'] *0 + np.random.normal(0.5, 0.1, len(concentrations_df['ptsa']))
concentrations_df['yield'] = concentrations_df['ptsa'] * 0 + np.random.normal(0.5, 0.1, len(concentrations_df['ptsa']))
concentrations_df.to_csv(data_folder + experiment_name + 'results/' + 'product_concentration.csv', index=False)