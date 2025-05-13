from robowski.settings import *
import json
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


# logging.basicConfig(level=logging.INFO)

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

target_file = f'{data_folder}simple-reactions/2023-09-14-run01/2023-09-14-run01.xlsx'

df = pd.read_excel(target_file, sheet_name='reactions_with_run_info')

timing_dicts = []
target_plate_barcode = 61

# iterate over all rows in the dataframe and if column 'plate_barcode' is equal to target_plate_barcode, load JSON from 'full_status' column and append to timing_dicts
for index, row in df.iterrows():
    if row['plate_barcode'] == target_plate_barcode:
        json_here = (row['full_status'])
        # convert json string to dict python object
        dict_here = json.loads(json_here)

        # keep only the timestamps
        for key, value in dict_here.items():
            dict_here[key] = value[-1]

        timing_dicts.append(dict_here)

# make dataframe from timing_dicts
df_timings = pd.DataFrame(timing_dicts, dtype=object)

latest_timestamp_of_the_entire_dataframe = df_timings.max().max()
print(f'Latest timestamp of the entire dataframe: {latest_timestamp_of_the_entire_dataframe}')

earlierst_timestamp_of_the_entire_dataframe = df_timings.min().min()
print(f'Earliest timestamp of the entire dataframe: {earlierst_timestamp_of_the_entire_dataframe}')

print(f'Latest minus earliest, in minutes {(latest_timestamp_of_the_entire_dataframe - earlierst_timestamp_of_the_entire_dataframe) / 60}')

latest_acetonitrile = df_timings['Acetonitrile'].max()

print('Latest global minus the latest acetonitrile', (latest_timestamp_of_the_entire_dataframe - latest_acetonitrile) / 60)