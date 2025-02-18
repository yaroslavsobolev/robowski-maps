import numpy as np
import os
import pandas as pd
import importlib
# data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
# craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'
organize_run_results = importlib.import_module("misc-scripts.organize_run_results")

list_of_runs = tuple(['2023-06-20-run01',
                '2023-06-21-run01',
                '2023-06-21-run02',
                '2023-06-22-run01',
                '2023-06-22-run02',
                '2023-06-22-run03',
                '2023-06-23-run01',
                '2023-06-23-run02',
                '2023-06-26-run01',
                '2023-06-26-run02',
                '2023-06-27-run01',
                '2023-06-27-run02',
                '2023-06-27-run03',
                '2023-06-28-run01',
                '2023-06-28-run02'])

for i, run_name in enumerate(list_of_runs):
    print(run_name)
    if i<13:
        print('skipping')
        continue
    organize_run_results.organize_run_structure(f'multicomp-reactions/{run_name}/')
    organize_run_results.outV_to_outC_by_lookup(experiment_name=f'multicomp-reactions/{run_name}/',
                           lookup_run='multicomp-reactions/2023-06-19-run01/')