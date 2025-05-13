import shortuuid
import pandas as pd
import numpy as np


def make_hantzsch_plane():
    excel_template = 'D:/Docs/Dropbox/robochem/data/BPRF/2024-03-04-run01/misc/2024-03-04-template.xlsx'

    # # For 400 mM ammonium acetate plane
    # output_excel_file = 'D:/Docs/Dropbox/robochem/data/BPRF/2024-03-04-run01/2024-03-04-run01.xlsx'
    # fixed_ammonium_acetate = 150 # corresponds to 400 mM

    # # For 100 mM ammonium acetate plane
    # output_excel_file = 'D:/Docs/Dropbox/robochem/data/BPRF/2024-09-30-run01/2024-09-30-run01.xlsx'
    # fixed_ammonium_acetate = 37.5 # corresponds to 100 mM

    # For 400 mM ammonium acetate plane
    output_excel_file = 'D:/Docs/Dropbox/robochem/data/BPRF/2024-10-01-run01/2024-10-01-run01.xlsx'
    fixed_ammonium_acetate = 150 # corresponds to 400 mM

    # open first worksheet from excel, foll columns with numbers and save back to the same excel file
    df = pd.read_excel(excel_template, sheet_name=0)

    # add 54*6 new rows with same values in all columns
    for i in range(54*6-1):
        df = df.append(df.iloc[0], ignore_index=True)

    # change the values of the 'plate_barcode' column to i + 20
    for i in range(6):
        for j in range(54):
            df.loc[i*54+j, 'plate_barcode'] = i + 20
            df.loc[i * 54 + j, 'plate_barcodes_for_dilution'] = i + 30
            df.loc[i * 54 + j, 'plate_barcodes_for_dilution_2'] = i + 40

    # copy index values to "local_index' and 'global_index' column
    df['local_index'] = df.index
    df['global_index'] = df.index

    # fill column 'uuid' with shortuuid.uuid()
    df['uuid'] = df['uuid'].apply(lambda x: shortuuid.uuid())


    colname_x = 'vol#methoxybenzaldehyde'
    colname_y = 'vol#ethyl_acetoacetate'
    colname_z = 'vol#ammonium_acetate'
    colname_ethanol = 'vol#Ethanol'
    row_indexer = 0
    randomized_list_of_rows = np.random.permutation(df.index)
    for xs in np.linspace(15, 150, 18):
        for ys in np.linspace(15, 150, 18):
            row_index_here = randomized_list_of_rows[row_indexer]
            df.loc[row_index_here, colname_x] = xs
            df.loc[row_index_here, colname_y] = ys
            df.loc[row_index_here, colname_z] = fixed_ammonium_acetate
            df.loc[row_index_here, 'global_index'] = row_indexer
            ethanol = 500 - xs - ys - fixed_ammonium_acetate
            df.loc[row_index_here, colname_ethanol] = ethanol
            row_indexer += 1


    # save to the same excel file, to the
    df.to_excel(output_excel_file, index=False)


def hantzsch_cube():
    excel_template = 'D:/Docs/Dropbox/robochem/data/BPRF/2024-03-04-run01/misc/2024-03-04-template.xlsx'

    nplates = 1

    # For 400 mM ammonium acetate plane
    output_excel_file = 'D:/Docs/Dropbox/robochem/data/BPRF/2024-10-10-run01/2024-10-10-run01.xlsx'
    ammonium_acetate = 150  # corresponds to 400 mM

    # open first worksheet from excel, foll columns with numbers and save back to the same excel file
    df = pd.read_excel(excel_template, sheet_name=0)

    # add 54*6 new rows with same values in all columns
    for i in range(54 * nplates - 1):
        df = df.append(df.iloc[0], ignore_index=True)

    # change the values of the 'plate_barcode' column to i + 20
    for i in range(nplates):
        for j in range(54):
            df.loc[i * 54 + j, 'plate_barcode'] = i + 20
            df.loc[i * 54 + j, 'plate_barcodes_for_dilution'] = i + 30
            df.loc[i * 54 + j, 'plate_barcodes_for_dilution_2'] = i + 40

    # copy index values to "local_index' and 'global_index' column
    df['local_index'] = df.index
    df['global_index'] = df.index

    # fill column 'uuid' with shortuuid.uuid()
    df['uuid'] = df['uuid'].apply(lambda x: shortuuid.uuid())

    colname_x = 'vol#methoxybenzaldehyde'
    colname_y = 'vol#ethyl_acetoacetate'
    colname_z = 'vol#ammonium_acetate'
    colname_ethanol = 'vol#Ethanol'
    row_indexer = 0
    randomized_list_of_rows = np.random.permutation(df.index.tolist()[:27])
    for ammonium_acetate in np.linspace(37.5, 150, 3):
        for aldehyde in np.linspace(15, 150, 3):
            for eaa in np.linspace(15, 150, 3):
                row_index_here = randomized_list_of_rows[row_indexer]
                df.loc[row_index_here, colname_x] = aldehyde
                df.loc[row_index_here, colname_y] = eaa
                df.loc[row_index_here, colname_z] = ammonium_acetate
                df.loc[row_index_here, 'global_index'] = row_indexer
                ethanol = 500 - aldehyde - eaa - ammonium_acetate
                df.loc[row_index_here, colname_ethanol] = ethanol
                row_indexer += 1

    # save to the same excel file, to the
    df.to_excel(output_excel_file, index=False)


if __name__ == '__main__':
    for i in range(20):
        print(shortuuid.uuid())