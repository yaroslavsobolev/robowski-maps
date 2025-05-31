#!/usr/bin/env python3
"""
CSV Decimator Script
Decimates a CSV dataset by the 'ptsa' column using a user-specified factor.
"""

import pandas as pd
import os
import sys


def decimate_csv(input_filename, decimate_by_factor):
    """
    Decimate a CSV file by the ptsa column.

    Args:
        input_filename (str): Path to input CSV file
        decimate_by_factor (int): Decimation factor (keep every Nth row)

    Returns:
        str: Path to output file
    """

    # Read the CSV file
    try:
        df = pd.read_csv(input_filename)
        print(f"Loaded CSV with {len(df)} rows")
    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        return None
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return None

    # Check if ptsa column exists
    if 'ptsa' not in df.columns:
        print("Error: 'ptsa' column not found in the CSV file.")
        print(f"Available columns: {list(df.columns)}")
        return None

    # Sort by ptsa column to ensure proper ordering for decimation
    df_sorted = df.sort_values('ptsa').reset_index(drop=True)
    print(f"Sorted data by 'ptsa' column")
    print(f"ptsa range: {df_sorted['ptsa'].min():.6f} to {df_sorted['ptsa'].max():.6f}")
    print(f"Unique ptsa values: {df_sorted['ptsa'].nunique()}")

    # Apply decimation (keep every Nth row)
    decimated_df = df_sorted.iloc[::decimate_by_factor].reset_index(drop=True)
    print(f"After decimation by factor {decimate_by_factor}: {len(decimated_df)} rows")

    # Generate output filename
    base_name, ext = os.path.splitext(input_filename)
    output_filename = f"{base_name}_decimated{ext}"

    # Save decimated data
    try:
        decimated_df.to_csv(output_filename, index=False)
        print(f"Decimated data saved to: {output_filename}")
        return output_filename
    except Exception as e:
        print(f"Error saving file: {e}")
        return None


if __name__ == "__main__":
    decimate_csv(input_filename='postprocessed_yields.csv',
                 decimate_by_factor=5)