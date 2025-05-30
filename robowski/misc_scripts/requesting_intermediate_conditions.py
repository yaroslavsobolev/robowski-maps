"""
Experimental Design Script: Generate Intermediate Reaction Conditions

This script processes experimental data from multiple chemical reaction runs and generates
new intermediate experimental conditions through linear interpolation. The goal is to
systematically explore the parameter space between existing experimental points to
identify optimal reaction conditions.

It is intended to be used interactively in an interpreter or as part of a larger
experiment planning workflow.

Workflow:
1. Load and consolidate data from multiple experimental runs of Ugi reaction
2. Filter out outliers and padding rows (zero-concentration conditions)
3. Extract substrate concentrations and catalyst variations for Ugi reaction
4. Generate intermediate conditions by interpolating between existing condition points
5. Prepare new experimental conditions for running them on the automated experimentation platform

The script works with multi-component reactions involving:
- Substrates: ic001, am001, ald001 (isocyanide, amine, aldehyde components)
- Catalyst: ptsa (p-toluenesulfonic acid)
- Product: IIO029A (target reaction product)

The interpolation generates both concentration and volume specifications for new
experiments, enabling systematic exploration of reaction space between known conditions.

Input Data:
- Multiple experimental run directories containing reaction results
- Each run contains spectral data, concentrations, and yield measurements
- Data includes outlier flags and experimental metadata

Output:
- Global dataframes `new_concentrations` and `new_volumes` populated with
  interpolated experimental conditions ready for laboratory execution

"""

from robowski.settings import *
import pandas as pd
import numpy as np
import os
import robowski.visualize_results.visualize_results as visualize_results

run_name = 'multicomp-reactions/2023-06-19-run01/'

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
                    '2023-06-28-run02',
                    '2023-06-28-run03'])

# Load experimental data from multiple reaction runs and combine into single dataset
df_results = visualize_results.join_data_from_runs([f'multicomp-reactions/{run}/' for run in list_of_runs])

# Clean data: remove outliers and rows with zero substrate concentrations (padding rows)
print(f"There are {df_results[df_results['is_outlier'] == 1].shape[0]} outliers.")
df_results = df_results[df_results['is_outlier'] == 0]

### Uncomment to load interpolated product concentration data instead of raw results
# df_results = pd.read_csv(data_folder + experiment_name + f'results/interpolated_product_concentration.csv')

# Convert negative yields to zero (physically meaningful constraint)
df_results['yield'] = df_results['yield'].apply(lambda x: 0 if x < 0 else x)

# Define chemical components: three substrates + one catalyst
substances = ['ic001','am001','ald001','ptsa']
product = 'IIO029A'

# Remove experimental conditions where all substrates are zero (control/blank conditions)
padding_rows_count = (df_results[substances] == 0).all(axis=1).sum()
print(f"There are {padding_rows_count} padding rows (with zero concentrations of substrates).")
df_results = df_results[(df_results[substances] != 0).any(axis=1)]

# Extract substrate concentration arrays for analysis
substrate_cs = []
for substance in substances:
    substrate_cs.append(df_results[substance].to_numpy())

xs0, ys0, zs0, cats = substrate_cs # Unpack: ic001, am001, ald001, ptsa concentrations

print('Max concentrations of substrates: ')
for x in [xs0, ys0, zs0]:
    print(max(x))

minimal_concentration_of_substrates = np.min(np.array([xs0, ys0, zs0]))

# Find unique catalyst concentrations for systematic interpolation
unique_cats = sorted(list(set(list(cats))))
print(f'Unique cats: {unique_cats}')


# make interpolation between two conditions
for substance in substances:
    print(substance)
    print(df_results[substance].unique())

# linear interpolation between two conditions
def linear_interp_between_df_rows(number_of_points_between, df_original, targets_ids,
                                  only_use_these=None):
    """
    Create linearly interpolated rows between two existing dataframe rows.

    Generates intermediate experimental conditions by linear interpolation between
    two specified rows in the dataframe. Each column is interpolated independently,
    creating a smooth transition between the target experimental conditions.

    Parameters
    ----------
    number_of_points_between : int
        Total number of interpolation points to generate between the two target rows.
        Excludes the endpoints (target rows themselves).
    df_original : pandas.DataFrame
        Source dataframe containing experimental conditions to interpolate between.
        Must contain numeric columns for substrate concentrations and experimental parameters.
    targets_ids : list of int
        Two-element list containing row indices [start_row, end_row] to interpolate between.
        These rows define the endpoints of the interpolation.
    only_use_these : list of int, optional
        Subset indices to select specific interpolated points. If None, returns all
        interpolated points. Useful for selecting every nth point or custom spacing.

    Returns
    -------
    pandas.DataFrame
        New dataframe containing interpolated experimental conditions with same column
        structure as input. Each row represents an intermediate experimental condition.
        Missing values are filled with zeros.

    """
    interp_conditions = df_original.iloc[:0,:].copy()
    for column in interp_conditions.columns:
        if only_use_these is not None:
            interp_conditions[column] = np.linspace(df_original[column].iloc[targets_ids[0]],
                                                    df_original[column].iloc[targets_ids[1]],
                                                    number_of_points_between)[1:-1][only_use_these]
        else:
            interp_conditions[column] = np.linspace(df_original[column].iloc[targets_ids[0]],
                                                   df_original[column].iloc[targets_ids[1]],
                                                   number_of_points_between)[1:-1]
    interp_conditions = interp_conditions.fillna(0)
    return interp_conditions

# Initialize global dataframes to collect interpolated experimental conditions
# These will be populated by the interpolation functions
new_concentrations = df_results.iloc[:0, :].copy()
new_volumes = df_volumes.iloc[:0, :].copy()


def append_interpolation(targets, number_of_points_between, only_use_these=None):
    """
    Find experimental conditions matching target parameters and interpolate between them.

    Locates dataframe rows matching specified substrate and catalyst concentrations,
    then generates interpolated experimental conditions between these target points.
    Results are appended to global dataframes for concentration and volume specifications.

    Parameters
    ----------
    targets : list of dict
        Two-element list containing target experimental conditions as dictionaries.
        Each dict must contain keys: 'ic001', 'am001', 'ald001', 'ptsa' with
        corresponding concentration values to match in the experimental data.
    number_of_points_between : int
        Number of interpolation points to generate between the target conditions.
    only_use_these : list of int, optional
        Indices to select subset of interpolated points. If None, uses all points.

    Global Variables Modified
    -------------------------
    new_concentrations : pandas.DataFrame
        Appended with interpolated concentration specifications.
    new_volumes : pandas.DataFrame
        Appended with interpolated volume specifications.

    """
    global new_concentrations, new_volumes
    targets_ids = [
        df_results[(df_results['ic001'] == target['ic001']) & (df_results['am001'] == target['am001'])
                   & (df_results['ald001'] == target['ald001']) & (df_results['ptsa'] == target['ptsa'])].index[0]
        for target in targets
    ]
    interpolated_points_concentration = linear_interp_between_df_rows(number_of_points_between=number_of_points_between,
                                                                      df_original=df_results.copy(),
                                                                      targets_ids=targets_ids,
                                                                      only_use_these=only_use_these)
    new_concentrations = new_concentrations.append(interpolated_points_concentration, ignore_index=True)
    interpolated_points_volumes = linear_interp_between_df_rows(number_of_points_between=number_of_points_between,
                                                                df_original=df_volumes.copy(),
                                                                targets_ids=targets_ids,
                                                                only_use_these=only_use_these)
    new_volumes = new_volumes.append(interpolated_points_volumes, ignore_index=True)

def append_interpolation_by_indices(param_values_by_index_start, param_values_by_index_stop,
                                    number_of_points_between, cats_range_indices,
                                    column_names=('ic001', 'am001', 'ald001'),
                                    catalyst_name='ptsa', only_use_these=None):
    """
    Generate systematic interpolations across multiple catalyst concentrations using parameter indices.

    Creates interpolated experimental conditions by mapping parameter indices to actual
    concentration values, then interpolating between start and stop conditions across
    a range of catalyst concentrations. This enables systematic exploration of the
    experimental parameter space.

    Parameters
    ----------
    param_values_by_index_start : list of int
        Parameter indices defining the starting experimental condition.
        Each index maps to sorted unique values for the corresponding substrate.
    param_values_by_index_stop : list of int
        Parameter indices defining the ending experimental condition.
        Must have same length as param_values_by_index_start.
    number_of_points_between : int
        Number of interpolation points to generate between start and stop conditions.
    cats_range_indices : list of int
        Catalyst concentration indices to include in the interpolation.
        Each index maps to sorted unique catalyst concentrations.
    column_names : tuple of str, optional
        Substrate column names corresponding to parameter indices.
        Default: ('ic001', 'am001', 'ald001') for isocyanide, amine, aldehyde.
    catalyst_name : str, optional
        Column name for catalyst concentration. Default: 'ptsa'.
    only_use_these : list of int, optional
        Subset indices for selecting specific interpolated points.

    Global Variables Modified
    -------------------------
    new_concentrations : pandas.DataFrame
        Extended with systematically interpolated concentration conditions.
    new_volumes : pandas.DataFrame
        Extended with corresponding volume specifications.

    """
    global new_concentrations, new_volumes
    for cats_id in cats_range_indices:
        targets = []
        for param_value_indices in [param_values_by_index_start, param_values_by_index_stop]:
            target = {column_name: sorted(df_results[column_name].unique())[param_value_indices[i]]
                      for i, column_name in enumerate(column_names)}
            target[catalyst_name] = unique_cats[cats_id]
            targets.append(target)
        append_interpolation(targets=targets, number_of_points_between=number_of_points_between,
                             only_use_these=only_use_these)