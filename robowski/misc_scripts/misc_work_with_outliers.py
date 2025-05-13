from robowski.settings import *


from matplotlib import pyplot as plt

from organize_run_results import *

df_output, yield_lists = merge_repeated_outliers(original_run = 'multicomp-reactions/2023-06-19-run01/',
                        outlier_runs=['multicomp-reactions/2023-06-30-run01/',
                                      'multicomp-reactions/2023-07-04-run01/'])
for j, yields in enumerate(yield_lists):
    for i, yield_here in enumerate(yields):
        plt.scatter(j, yield_here, color=f'C{i}')
plt.legend()
plt.show()

run_name = '2023-07-04-run01'
organize_run_structure(f'multicomp-reactions/{run_name}/')
outV_to_outC_by_lookup(experiment_name=f'multicomp-reactions/{run_name}/',
                       lookup_run='multicomp-reactions/2023-06-19-run01/')