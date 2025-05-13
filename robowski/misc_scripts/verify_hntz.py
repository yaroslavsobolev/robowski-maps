from robowski.settings import *
from hntz_discovery_timeline import *
from assert_identity_of_two_csv_tables import compare_csv_tables

if __name__ == '__main__':
    list_of_runs = tuple(['2024-02-16-run01'])
    ## Uncomment to recalibrate the analysis
    for set_index in [0]:
        print(f'Processing calibrant set {set_index}...')
        combined_set = []
        for i in range(set_index+1):
            combined_set += calibrant_sets[i]['calibrant_set']
        print(f'Combined set: {combined_set}')
        folder_for_prod_conc = f'results/historical/calibrantset_{set_index}'

        for i, run_shortname in enumerate(list_of_runs):
            process_run_by_shortname(run_shortname, combined_set, folder_for_prod_conc)

    datafolder1 = 'E:/robowski-maps-data-upload/'
    datafolder2 = 'E:/hntz_expected/'
    run_name = 'BPRF/2024-02-16-run01/'
    csv_filepath_1 = datafolder1 + run_name + 'results/historical/calibrantset_0/' + 'product_concentration.csv'
    csv_filepath_2 = datafolder2 + run_name + 'results/' + 'product_concentration.csv'
    are_identical, result = compare_csv_tables(csv_filepath_1, csv_filepath_2)

    if are_identical:
        print("The tables are identical.")
    else:
        print("The tables are not identical.")
        print("Differences:")
        print(result)