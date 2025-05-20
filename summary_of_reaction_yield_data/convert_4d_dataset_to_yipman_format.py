import pandas as pd
import argparse

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert 4D dataset into the format of interactive 3D viewer (YipMan).")
    parser.add_argument("--input_csv", required=True, help="Path to the input CSV file.")
    parser.add_argument("--output_name", required=True, help="Path to the output CSV file.")
    parser.add_argument("--prefix_of_relative_dir", required=True, help="Path to the output CSV file.")
    parser.add_argument("--X", required=True, help="Column name for X coordinate.")
    parser.add_argument("--Y", required=True, help="Column name for Y coordinate.")
    parser.add_argument("--Z", required=True, help="Column name for Z coordinate.")
    parser.add_argument("--T", required=True, help="Column name for time (fourth dimension) coordinate.")
    parser.add_argument("--V", required=True, help="Column name for default scalar field to be plotted.")
    parser.add_argument("--Xscale", required=False, help="Scale factor for X coordinate.", default=1.0)
    parser.add_argument("--Yscale", required=False, help="Scale factor for Y coordinate.", default=1.0)
    parser.add_argument("--Zscale", required=False, help="Scale factor for Z coordinate.", default=1.0)
    parser.add_argument("--Xrename", required=False, help="New name for X coordinate column.", default=None)
    parser.add_argument("--Yrename", required=False, help="New name for Y coordinate column.", default=None)
    parser.add_argument("--Zrename", required=False, help="New name for Z coordinate column.", default=None)
    parser.add_argument("--Tlabel_prefix", required=False, help="Label prefix for the fourth dimension coordinate.", default=None)
    parser.add_argument("--Tlabel_suffix", required=False, help="Label sufffix for the fourth dimension coordinate.",
                        default=None)
    args = parser.parse_args()

    # Load the input CSV file
    df_input = pd.read_csv(args.input_csv)

    # make a dataframe with fields 'filename', 'path', 't', 't_label'
    df_summary = pd.DataFrame(columns=['filename', 'path', 't', 't_label'], dtype=object)

    # make a list of unique values of the T column
    unique_t_values = df_input[args.T].unique()
    print(f'Number of unique time points: {len(unique_t_values)}')

    for i, unique_t in enumerate(unique_t_values):
        df = df_input[df_input[args.T] == unique_t]

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
        filepath = f'{args.output_name}_{i:04d}.csv'
        df.to_csv(filepath, index=False)

        # since filepath can contain slashes, we need to get the last part of the path
        filename = filepath.split('/')[-1]
        df_summary = pd.concat([df_summary,
                                pd.DataFrame([[filename,
                                               f'{args.prefix_of_relative_dir}/{filename}',
                                               unique_t,
                                               f'{args.Tlabel_prefix}{unique_t:.3f}{args.Tlabel_suffix}']],
                                             columns=['filename', 'path', 't', 't_label'])], ignore_index=True)

    # Save the summary DataFrame to the output CSV file, but without the column names
    df_summary.to_csv(f'{args.output_name}_tesseract.csv', index=False, header=False)

if __name__ == "__main__":
    main()