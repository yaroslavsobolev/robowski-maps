import os
import shutil

import numpy as np
import pandas as pd

import config


def clean_excel_folder(excel_path):
    for excel in excel_path:
        excel_dir = os.path.dirname(excel)
        ## remove pipetter_io folder if exists
        if os.path.exists(excel_dir + '/pipetter_io'):
            input_usr = input(f"Remove pipetter_io folder in {excel_dir}? (y/n)")
            if input_usr == 'y' or input_usr == 'Y':
                shutil.rmtree(excel_dir + '/pipetter_io')
            print(f"Removed pipetter_io folder in {excel_dir}. You can proceed.")

        ## remove dilution folder if exists
        if os.path.exists(excel_dir + '/dilution'):
            input_usr = input(f"Remove dilution folder in {excel_dir}? (y/n)")
            if input_usr == 'y' or input_usr == 'Y':
                shutil.rmtree(excel_dir + '/dilution')
            print(f"Removed dilution folder in {excel_dir}. You can proceed.")


def process_all_excel_files(excel_path):

    for excel in excel_path:
        ## check if excel name is the sama as the excel folder name
        excel_name = os.path.basename(excel).split('.')[0]
        excel_dir = os.path.dirname(excel).split('/')[-1]
        assert excel_name == excel_dir, f"Excel file name: {excel_name} not same as Excel folder name: {excel_dir}"

        ## process excel files one by one
        df = pd.read_excel(excel, sheet_name= config.sheet_name_for_run_info, engine='openpyxl')
        process_one_excel_file(df=df, excel_path=excel)


def get_run_info(df_in_one_plate, excel_name):
    plate_barcode = df_in_one_plate.iloc[0]['plate_barcode']
    start_time_unix = df_in_one_plate.iloc[0]['timestamp']
    start_time_datetime = pd.to_datetime(start_time_unix, unit='s')
    end_time_unix = df_in_one_plate.iloc[-1]['timestamp']
    end_time_datetime = pd.to_datetime(end_time_unix, unit='s')
    return f'{plate_barcode},{excel_name},{start_time_unix},{start_time_datetime},{end_time_unix},{end_time_datetime},,\n'


def get_dilution_info(df_in_one_plate):
    for i in range(1, len(df_in_one_plate)):
        assert df_in_one_plate.iloc[i]['plate_barcode'] == df_in_one_plate.iloc[i-1]['plate_barcode'],\
            f"plate_barcode: {df_in_one_plate.iloc[i]['plate_barcode']} not same " \
            f"as previous plate_barcode: {df_in_one_plate.iloc[i-1]['plate_barcode']}"
        assert df_in_one_plate.iloc[i]['plate_barcodes_for_dilution'] == df_in_one_plate.iloc[i-1]['plate_barcodes_for_dilution'],\
            f"plate_barcodes_for_dilution: {df_in_one_plate.iloc[i]['plate_barcodes_for_dilution']} not same " \
            f"as previous plate_barcodes_for_dilution: {df_in_one_plate.iloc[i-1]['plate_barcodes_for_dilution']}"

    dilution_from = df_in_one_plate.iloc[0]['plate_barcode']
    dilution_to = df_in_one_plate.iloc[0]['plate_barcodes_for_dilution']

    return f'{dilution_from},{dilution_to}\n'


def process_one_excel_file(df, excel_path):
    ## get split index
    split_index = []
    for index in range(1, len(df)):
        if df.iloc[index]['plate_barcode'] != df.iloc[index-1]['plate_barcode']:
            split_index.append(index)
    ## group by split index
    df_grouped_by_plate_id = np.split(df, split_index)
    # print(df_grouped_by_plate_id)

    ## get excel directory from excel path
    excel_directory = os.path.dirname(excel_path)
    excel_name = os.path.basename(excel_path).split('.')[0]
    ## make a folder named pipetter_io if not exist in the excel directory
    if not os.path.exists(excel_directory + '/pipetter_io'):
        os.mkdir(excel_directory + '/pipetter_io')
        ## make a file called run_info.csv in the folder if not exist
        if not os.path.exists(excel_directory + '/pipetter_io/run_info.csv'):
            with open(excel_directory + '/pipetter_io/run_info.csv', 'w') as f:
                f.write(f'{config.version_of_run_info}\n')
    ## make a folder named dilution if not exist
    if not os.path.exists(excel_directory + '/dilution'):
        os.mkdir(excel_directory + '/dilution')
        col_name = 'from,to'
        ## make a file called dilution_info.csv in the folder if not exist
        if not os.path.exists(excel_directory + '/dilution/dilution_info.csv'):
            with open(excel_directory + '/dilution/dilution_info.csv', 'w') as f:
                f.write(f'{config.version_of_dilution_info}\n')
                f.write(f'{col_name}\n')

    ## creacte files for run info and dilution info
    for df_in_one_plate in df_grouped_by_plate_id:
        run_info = get_run_info(df_in_one_plate = df_in_one_plate, excel_name = excel_name)
        dilution_info = get_dilution_info(df_in_one_plate = df_in_one_plate)
        with open(excel_directory + '/pipetter_io/run_info.csv', 'a') as f:
            f.write(run_info)
        with open(excel_directory + '/dilution/dilution_info.csv', 'a') as f:
            f.write(dilution_info)

    ## log
    print(f'run_info.csv and dilution_info.csv created in {excel_directory}')


if __name__ == '__main__':

    data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
    reaction_name = 'simple-reactions'
    run_name = ['2023-07-26-run01',
                '2023-07-26-run02',
                '2023-07-27-run01',
                '2023-07-27-run02',
                '2023-07-28-run01'
                ]
    excel_path = [data_folder + reaction_name + '/' + run + '/' + f'{run}.xlsx' for run in run_name]

    clean_excel_folder(excel_path)

    process_all_excel_files(excel_path)

