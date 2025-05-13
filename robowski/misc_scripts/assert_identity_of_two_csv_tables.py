from robowski.settings import *
import pandas as pd
import numpy as np


def compare_csv_tables(csv_filepath_1, csv_filepath_2):
    """
    Compare two CSV tables for equality.

    Args:
        csv_filepath_1 (str): Path to the first CSV file.
        csv_filepath_2 (str): Path to the second CSV file.

    Returns:
        bool: True if the tables are identical, False otherwise.
        pd.DataFrame: DataFrame highlighting differences (if any).
    """
    # Load the CSV files into DataFrames
    df1 = pd.read_csv(csv_filepath_1)
    df2 = pd.read_csv(csv_filepath_2)

    # Check if shapes are the same
    if df1.shape != df2.shape:
        return False, f"Shape mismatch: {df1.shape} vs {df2.shape}"

    # Initialize a DataFrame to store differences
    differences = pd.DataFrame(index=df1.index, columns=df1.columns, dtype=object)

    # Compare each element
    for col in df1.columns:
        if col not in df2.columns:
            return False, f"Column '{col}' not found in the second table."

        if pd.api.types.is_numeric_dtype(df1[col]) and pd.api.types.is_numeric_dtype(df2[col]):
            # Compare numerical values using np.isclose
            differences[col] = ~np.isclose(df1[col], df2[col], equal_nan=True)
        else:
            # Compare non-numerical values directly
            differences[col] = df1[col] != df2[col]

    # Check if there are any differences
    if differences.any().any():
        return False, differences
    return True, None


def compare_results_of_same_run(run_name):
    datafolder1 = 'E:/robowski-maps-data-upload/'
    datafolder2 = 'D:/Docs/Dropbox/robochem/data/'
    csv_filepath_1 = datafolder1 + run_name + 'results/' + 'product_concentration.csv'
    csv_filepath_2 = datafolder2 + run_name + 'results/' + 'product_concentration.csv'
    are_identical, result = compare_csv_tables(csv_filepath_1, csv_filepath_2)

    if are_identical:
        print("The tables are identical.")
    else:
        print("The tables are not identical.")
        print("Differences:")
        print(result)


# Example usage
if __name__ == "__main__":
    # compare_results_of_same_run('simple-reactions/2023-07-05-run01/')

    compare_results_of_same_run('simple-reactions/2023-09-06-run01/')


    # csv_filepath_1 = "table1.csv"
    # csv_filepath_2 = "table2.csv"
    # are_identical, result = compare_csv_tables(csv_filepath_1, csv_filepath_2)
    #
    # if are_identical:
    #     print("The tables are identical.")
    # else:
    #     print("The tables are not identical.")
    #     print("Differences:")
    #     print(result)