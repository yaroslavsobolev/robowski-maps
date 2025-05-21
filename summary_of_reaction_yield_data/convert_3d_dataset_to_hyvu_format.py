# this is a command line script that loads the given target .CSV table, reorders the column according to command-line
# arguments (--X, --Y, --Z, --V), and saves the result to a new .CSV file. First column remains the first column. Second column is X.
# Third column is Y. Fourth column is Z. Fifth column is V.
# If the table does not have 'spectrum_dir' column, it will be added
# to the table and filled with 'C:\' string in all rows.
#
# Example command:
# python convert_3d_dataset_to_hyvu_format.py --input_csv=SN1/raw_yields.csv --output_csv=SN1/SN1_raw_yields_hyvu.csv --X='c#SN1OH01' --Y='c#HBr' --Z='Temperature' --V='yield'
#
# Sequence of columns in the original .CSV:
# reactions,DMF,SN1OH01,SN1OH01_dil_x_100,HBr,HBr_dil_x_10,HBr_dil_x_10000,vial_id,reaction_plate_id,diluted_plate_id,craic_folder,is_outlier,craic_folder_undil,c#SN1OH01,c#HBr,Temperature,pc#SN1OH01,pc#SN1Br01s1,pc#carbocat,yield,yield#SN1Br01,yield#SN1Br01s1,run_name
#
# Sequence of columns in the new .CSV:
# reactions,c#SN1OH01,c#HBr,Temperature,yield,DMF,SN1OH01,SN1OH01_dil_x_100,HBr,HBr_dil_x_10,HBr_dil_x_10000,vial_id,reaction_plate_id,diluted_plate_id,craic_folder,is_outlier,craic_folder_undil,pc#SN1OH01,pc#SN1Br01s1,pc#carbocat,yield#SN1Br01,yield#SN1Br01s1,run_name,spectrum_dir

import pandas as pd
import argparse

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Reorder columns in a CSV file for compliance with the format of interactive 3D viewer (HyperspaceViewer).")
    parser.add_argument("--input_csv", required=True, help="Path to the input CSV file.")
    parser.add_argument("--output_csv", required=True, help="Path to the output CSV file.")
    parser.add_argument("--X", required=True, help="Column name for X coordinate.")
    parser.add_argument("--Y", required=True, help="Column name for Y coordinate.")
    parser.add_argument("--Z", required=True, help="Column name for Z coordinate.")
    parser.add_argument("--V", required=True, help="Column name for default scalar field to be plotted.")
    parser.add_argument("--Xscale", required=False, help="Scale factor for X coordinate.", default=1.0)
    parser.add_argument("--Yscale", required=False, help="Scale factor for Y coordinate.", default=1.0)
    parser.add_argument("--Zscale", required=False, help="Scale factor for Z coordinate.", default=1.0)
    parser.add_argument("--Xrename", required=False, help="New name for X coordinate column.", default=None)
    parser.add_argument("--Yrename", required=False, help="New name for Y coordinate column.", default=None)
    parser.add_argument("--Zrename", required=False, help="New name for Z coordinate column.", default=None)
    args = parser.parse_args()

    # Load the input CSV file
    df = pd.read_csv(args.input_csv)

    # Ensure the specified columns exist in the input CSV
    for col in [args.X, args.Y, args.Z, args.V]:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in the input CSV file.")

    # Reorder columns
    new_order = [args.V, args.X, args.Y, args.Z] + \
                [col for col in df.columns if col not in [args.X, args.Y, args.Z, args.V]]

    # Reorder the DataFrame
    df = df[new_order]

    # Add 'spectrum_dir' column if it does not exist
    if 'spectrum_dir' not in df.columns:
        df['spectrum_dir'] = 'C:\\'

    # scale the values in X, Y, Z columns
    df[args.X] = df[args.X].astype(float) * float(args.Xscale)
    df[args.Y] = df[args.Y].astype(float) * float(args.Yscale)
    df[args.Z] = df[args.Z].astype(float) * float(args.Zscale)

    # Rename the columns if new names are provided
    if args.Xrename:
        df.rename(columns={args.X: args.Xrename}, inplace=True)
    if args.Yrename:
        df.rename(columns={args.Y: args.Yrename}, inplace=True)
    if args.Zrename:
        df.rename(columns={args.Z: args.Zrename}, inplace=True)

    # Save the reordered DataFrame to the output CSV file
    df.to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()