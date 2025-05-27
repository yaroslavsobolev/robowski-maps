"""
This module analyzes an Excel file and tells the schedule
of actions needs to be done.

"""

import PySimpleGUI as sg, pandas as pd, os, numpy as np
from datetime import datetime
# from pandasgui import show
import config


## get Exce file path by pysimplegui
def get_excel_path():
    example_path = 'C:/Users/Chemiluminescence/Dropbox/robochem/data/simple-reactions/2023-08-07-run01/2023-08-07-run01.xlsx'.replace('\\', '/')

    sg.theme('DarkAmber')
    layout = [[sg.Text('Please select the Excel file to be processed.', font='Any 15')],
              [sg.Input(example_path,font='Any 15'), sg.FileBrowse( font='Any 15')],
              [sg.OK(font='Any 15'), sg.Cancel(font='Any 15')]]

    window = sg.Window('Excel file selection', layout, size=(700, 150))
    event, values = window.read()
    window.close()
    return values[0]

def get_info_from_excel(excel_path):

    df = pd.read_excel(excel_path, sheet_name='reactions_with_run_info', engine='openpyxl')
    ## add one column "recation_duration" to df if not exist
    if 'reaction_duration' not in df.columns:
        duration_from_user = input('Please enter the duration of the reaction in hours: ')
        df['reaction_duration'] = duration_from_user
        ## save the df to excel
        with pd.ExcelWriter(excel_path, engine='openpyxl', mode='a', if_sheet_exists="replace") as writer:
            df.to_excel(writer, sheet_name=config.sheet_name_for_run_info, index=False)

    ## split df according to 'plate_barcode'
    split_index = [index for index in range(1, len(df)) if df.iloc[index]['plate_barcode'] != df.iloc[index-1]['plate_barcode']]
    split_list = np.split(df, split_index)

    ## get info for each plate
    info_list = []
    for plate_df in split_list:
        info = {}
        reaction_duration_in_seconds = plate_df.iloc[-1]['reaction_duration'] * 3600
        info['plate_pipetting'] = plate_df.iloc[-1]['plate_barcode']
        info['plate_dilution'] = plate_df.iloc[-1]['plate_barcodes_for_dilution']
        info['time_pipetting'] =  datetime.fromtimestamp (plate_df.iloc[-1]['timestamp']).strftime('%Y-%m-%d %H:%M')
        info['reaction_duration'] = plate_df.iloc[-1]['reaction_duration']
        info['times_UV-Vis/dilution'] = datetime.fromtimestamp (plate_df.iloc[-1]['timestamp'] + reaction_duration_in_seconds).strftime('%Y-%m-%d %H:%M')
        info_list.append(info)

    return pd.DataFrame(info_list)

def show_dataframe_by_pysimplegui(df, path: str):
    run_name = path.split('/')[-1]
    sg.theme('DarkAmber')
    layout = [[sg.Text("run_name:", font='Any 20'), sg.Text(f"{run_name}", font='Any 20 bold', text_color='lime')],
              [sg.Table(values=df.values.tolist(),
                        headings=df.columns.tolist(),
                        display_row_numbers=True,
                        auto_size_columns=True,
                        num_rows=min(25, len(df)),
                        size=(1000, 600),
                        row_height=50,
                        justification='center',
                        text_color='lime',
                        font='Any 15')
               ]
             ]

    window = sg.Window('Table', layout, grab_anywhere=False, size=(1000, 300))
    event, values = window.read()
    print(os.path.dirname(path) + '/schedule.png')
    # window.save_window_screenshot_to_disk(os.path.dirname(path) + '/schedule.png')
    window.close()

if __name__ == '__main__':
    excel_path = get_excel_path()
    df_plate = get_info_from_excel(excel_path = excel_path)
    show_dataframe_by_pysimplegui(df = df_plate, path = excel_path)

