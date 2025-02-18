import importlib
import os

import matplotlib.pyplot as plt
import numpy as np

organize_run_results = importlib.import_module("misc-scripts.organize_run_results")
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
craic_folder = data_folder + 'craic_microspectrometer_measurements/absorbance/'


# make a bar plot with compositions from dictionary:
composition = {'water': 22, 'dioxane': 10, 'acetic acid': 8, 'HBr':1}
# make a stacked bar plot with compositions from dictionary:
net_width = 0
for key in ['water', 'dioxane', 'acetic acid', 'HBr']:
    net_width += composition[key]
    plt.barh(y=0, width=composition[key], height=0.5, left=net_width-composition[key])

# annotate the plot with compositions:
# net_width = 0
# for key in ['water', 'dioxane', 'acetic acid', 'HBr']:
#     net_width += composition[key]
#     plt.annotate(key, (net_width-composition[key]/2, 0), ha='center', va='center', fontsize=12)
plt.xlabel('Molar composition')
# remove all axes
plt.gca().axes.get_yaxis().set_visible(False)
plt.gca().axes.get_xaxis().set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.show()




# experiment_name = 'multicomp-reactions/2023-07-04-run01/'
# df = organize_run_results.load_df_from_run_info(data_folder + experiment_name + 'pipetter_io/run_info.csv')
# df.to_pickle('tests/test_organize_run_results/expected_outputs/run_info_v1.00.pkl')
# df = organize_run_results.load_df_from_dilution_info(experiment_name)

# df = organize_run_results.join_data_from_runs(['multicomp-reactions/2023-06-20-run01/',
#                                                'multicomp-reactions/2023-06-21-run01/',
#                                                'multicomp-reactions/2023-06-21-run02/'
#                                               ])
# df.to_pickle('tests/test_organize_run_results/expected_outputs/joined_data_from_runs.pkl')

# import matplotlib.pyplot as plt
# import numpy as np
# from mayavi import mlab
#
# mlab.clf()
# x, y, z = np.mgrid[-5:5:30j, -5:5:30j, -5:5:30j]
# values = x*x*0.5 + y*y + z*z*2.0
# print(np.max(values))
# mlab.contour3d(values, contours=[0.05, 30, 66], opacity=0.3, vmin=0, vmax=70.0, extent=[-5, 5, -5, 5, -5, 5],
#                resolution=32)
# mlab.show()

# data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
# list_of_runs = tuple(['2023-06-20-run01',
#                       '2023-06-21-run01',
#                       '2023-06-21-run02',
#                       '2023-06-22-run01',
#                       '2023-06-22-run02',
#                       '2023-06-22-run03',
#                       '2023-06-23-run01',
#                       '2023-06-23-run02',
#                       '2023-06-26-run01',
#                       '2023-06-26-run02',
#                       '2023-06-27-run01',
#                       '2023-06-27-run02',
#                       '2023-06-27-run03',
#                       '2023-06-28-run01',
#                       '2023-06-28-run02',
#                       '2023-06-28-run03'])
# for run in list_of_runs:
#     dilution_filename = data_folder + 'multicomp-reactions/' + run + '/dilution/dilution_info.csv'
#     with open(dilution_filename, 'r') as original:
#         data = original.read()
#         if '#version' not in data:
#             print(f'run {run}, prepending version number to dilution_info.txt')
#             with open(dilution_filename, 'w') as modified:
#                 modified.write("#version: 1.00\n" + data)

