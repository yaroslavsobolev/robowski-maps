from robowski.settings import *
from matplotlib import pyplot as plt
from robowski.misc_scripts.organize_run_results import *
from robowski.uv_vis_absorption_spectroscopy.examples.Ugi_reaction.process_2023_06_19_run01 import process_run_by_shortname

# this should be executed after experiments for this run have been completed

# organize run structure for the run that was used to remeasure the outliers
run_name = '2023-07-04-run01'
organize_run_structure(f'multicomp-reactions/{run_name}/')
outV_to_outC_by_lookup(experiment_name=f'multicomp-reactions/{run_name}/',
                       lookup_run='multicomp-reactions/2023-06-19-run01/')

# processing the run dedicated to remeasuring the conditions that were marked as outliers in the main runs
run_shortname = '2023-07-04-run01'
process_run_by_shortname(run_shortname)

# After the data of the 2023-07-04-run01 has been processed (unmixing performed),
# merge the data of the repeated outliers run (2023-07-04-run01) with the data from the original run
df_output, yield_lists = merge_repeated_outliers(original_run = 'multicomp-reactions/2023-06-19-run01/',
                        outlier_runs=['multicomp-reactions/2023-06-30-run01/',
                                      'multicomp-reactions/2023-07-04-run01/'])

for j, yields in enumerate(yield_lists):
    for i, yield_here in enumerate(yields):
        plt.scatter(j, yield_here, color=f'C{i}')
plt.legend()
plt.show()