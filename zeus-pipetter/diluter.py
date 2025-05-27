"""
Serial Dilution Controller

This module automates serial dilution workflows using the liquid handling system.
It reads dilution parameters from an Excel file and executes pipetting steps via robotic hardware.

Main Features:
--------------
- Initializes hardware: Zeus module, gantry, and pipetter
- GUI-driven input for selecting Excel files, dilution parameters, and containers
- Supports one- and two-step dilutions (A→B, A→C, B→C)
- Dynamically generates pipetting events from Excel input
- Logs all actions to console and file (`main_roboski2.log`)
- Allows manual volume override and event re-execution

Usage:
------
Run as a script to interactively plan and execute dilutions

Author: Yankai Jia
"""

import logging, winsound, os
import copy
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
logger_path = data_folder[:-5] + 'pipetter_files/main_roboski2.log'
def setup_logger():
    # better logging format in console
    class CustomFormatter(logging.Formatter):
        grey = "\x1b[38;20m"
        yellow = "\x1b[33;20m"
        red = "\x1b[31;20m"
        bold_red = "\x1b[31;1m"
        reset = "\x1b[0m"
        format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"

        FORMATS = {
            logging.DEBUG: grey + format + reset,
            logging.INFO: yellow + format + reset,
            logging.WARNING: yellow + format + reset,
            logging.ERROR: red + format + reset,
            logging.CRITICAL: bold_red + format + reset
        }

        def format(self, record):
            log_fmt = self.FORMATS.get(record.levelno)
            formatter = logging.Formatter(log_fmt)
            return formatter.format(record)

    # create logger with 'main'
    logger = logging.getLogger('diluter')
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(logger_path)
    fh.setLevel(logging.INFO)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(CustomFormatter())
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

module_logger = setup_logger()

import time, pickle, re, importlib, json, os, sys, PySimpleGUI as sg, pandas as pd
import config
import zeus, pipetter, planner as pln, breadboard as brb, config

sg.theme('DarkAmber')

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

events_A_to_B = {0: [], 1: [], 2: []}  # {dilution_step: events}
events_A_to_C = {0: [], 1: [], 2: []}
events_B_to_C = {0: [], 1: [], 2: []}


def initiate_hardware() -> (zeus.ZeusModule, pipetter.Gantry, pipetter.Pipetter):
    # initiate zeus
    zm = zeus.ZeusModule(id=1)
    time.sleep(3)
    module_logger.info("zeus is loaded as: zm")

    # initiate gantry
    gt = pipetter.Gantry(zeus=zm)
    time.sleep(3)
    module_logger.info("gantry is loaded as: gt")
    # gt.configure_grbl() # This only need to be done once.
    gt.home_xy()
    if gt.xy_position == (0, 0):
        module_logger.info("gantry is homed")

    # initiate pipetter
    pt = pipetter.Pipetter(zeus=zm, gantry=gt)
    time.sleep(2)
    module_logger.info("pipetter is loaded as: pt")

    return zm, gt, pt


## Function to load the Excel file with run info by pysimplegui
def load_excel_path_by_pysimplegui():
    # change the theme of the GUI
    sg.theme('DarkAmber')
    layout = [
        [sg.Text('Please select the Excel file with run info for dilution')],
        [sg.Text('Excel file'), sg.Input(size=(30, 10)), sg.FileBrowse()],
        [sg.Submit(), sg.Cancel()]
    ]
    window = sg.Window('Excel file', layout, size=(1000, 300), text_justification='c', font='Helvetica 25')
    event, values = window.read()
    window.close()
    excel_path = values[0]
    module_logger.info(f'Excel file path: {excel_path}')
    print(f'Excel file path: {excel_path}')

    return excel_path


def get_liquid_class_table_index(solvent: str, tip_type: str, mode: str = 'empty') -> dict:
    liquid_class_dict = {
        'water_empty_50ul_clld': 21,
        'water_empty_300ul_clld': 1,
        'water_empty_1000ul_clld': 2,
        'water_part_50ul_clld': 3,
        'water_part_300ul_clld': 4,
        'water_part_1000ul_clld': 5,
        'serum_empty_50ul_clld': 6,
        'serum_empty_300ul_clld': 7,
        'serum_empty_1000ul_clld': 8,
        'serum_part_50ul_clld': 9,
        'serum_part_300ul_clld': 10,
        'serum_part_1000ul_clld': 11,
        'ethanol_empty_50ul_plld': 12,
        'ethanol_empty_300ul_plld': 13,
        'ethanol_empty_1000ul_plld': 14,
        'glycerin_empty_50ul_plld': 15,
        'glycerin_empty_300ul_plld': 16,
        'glycerin_empty_1000ul_plld': 17,
        'DMF_empty_300ul_clld': 22,
        'DMF_empty_1000ul_clld': 23,
        'DMF_empty_50ul_clld': 24,
        'Dioxane_empty_50ul_plld': 27,
        'Dioxane_empty_300ul_plld': 28,
        'Dioxane_empty_1000ul_plld': 29,
        'DCE_empty_50ul_clld': 36,
        'DCE_empty_300ul_clld': 37,
        'DCE_empty_1000ul_clld': 38,
    }

    solvent_para = {solvent, mode, tip_type}
    # define a set of paras
    for liquid_class, index in liquid_class_dict.items():
        solvent_para_here = set(liquid_class.split('_'))
        # print(f'solvent_para_here: {solvent_para_here}')
        # print(f'solvent_para: {solvent_para}')
        if solvent_para.issubset(solvent_para_here):
            return index


def generate_events_for_one_container(excel_path: str,
                                      dilution_model:str,
                                      df_for_one_container: object,
                                      solvent: str,
                                      diluting_solvent_container: object,
                                      source_container: object,
                                      source_container_volume: float,
                                      destination_container: object,
                                      dilution_factor: float,
                                      step_id: int,
                                      final_volume: float = 1000,
                                      source_container_diluted_factor=1, ):
    with open(data_folder + 'pipetter_files\\event_template.pickle', 'rb') as f:
        event = pickle.load(f)
    ## this is the content of the event_template.pickle
    event_long_str = """{
    'condition_uuid': '466QJxA8QwRJt39GCyLbER', 
    'plate_barcode': 11, 
    'excel_path_for_conditions': 'C:/Users/Chemiluminescence/Dropbox/robochem/data/simple-reactions/2023-07-17-run01/2023-07-17-run01.xlsx', 
    'substance': 'Dioxane', 
    'event_label': '466QJxA8QwRJt39GCyLbER_Dioxane',
    'source_container': Container(name='jar_100ml', container_id='', containerGeometryTableIndex=2, container_shape='cylindrical', diameter=520, bottomHeight=0, bottomSection=10000, bottomPosition=2120, immersionDepth=30, leavingHeight=40, jetHeight=130, startOfHeightBottomSearch=50, dispenseHeightAfterBottomSearch=50, liquid_volume=6519.7, volume_max=100000, area=2123.7, min_z=15.0, top_z=70, safety_margin_for_lldsearch_position=40, solvent='Dioxane', liquid_surface_height=1813, xy=[-322, -165], substance='Dioxane', substance_density=1.033), 
    'destination_container': Container(name='vial_2ml', container_id='', containerGeometryTableIndex=0, container_shape='cylindrical', diameter=118, bottomHeight=0, bottomSection=10000, bottomPosition=2210, immersionDepth=20, leavingHeight=20, jetHeight=130, startOfHeightBottomSearch=0, dispenseHeightAfterBottomSearch=80, liquid_volume=0, volume_max=2000, area=75.7, min_z=48, top_z=32, safety_margin_for_lldsearch_position=40, solvent='', liquid_surface_height=2110, xy=(0.0, -235.5), substance='', substance_density=1.0),
    'transfer_volume': 448.2040188, 
    'tip_type': '1000ul', 
    'asp_containerGeometryTableIndex': 2,
     'asp_deckGeometryTableIndex': 1, 
    'liquidClassTableIndex': 29, 
    'disp_containerGeometryTableIndex': 0, 
    'disp_deckGeometryTableIndex': 1, 
    'asp_qpm': 1, 
    'asp_lld': 0, 
    'disp_lld': 0, 
    'asp_mixVolume': 0, 
    'asp_mixFlowRate': 0, 
    'asp_mixCycles': 0, 
    'disp_mixVolume': 0, 
    'disp_mixFlowRate': 0, 
    'disp_mixCycles': 0, 
    'searchBottomMode': 0, 
    'is_event_conducted': False, 
    'event_finish_time': None}
      """

    event.condition_uuid = df_for_one_container['uuid']
    event.excel_path_for_conditions = excel_path
    event.plate_barcode = None
    event.substance = None
    event.event_label = None
    event.solvent = solvent

    transfer_volume_step0 = 0
    transfer_volume_step1 = 0
    transfer_volume_step2 = 0

    assert dilution_factor > 1, f'dilution_factor is {dilution_factor}, which is not valid.'

    if dilution_factor <= 100:  # A-> B
        transfer_volume_step0 = 100
        transfer_volume_step1 = final_volume / dilution_factor
        transfer_volume_step2 = final_volume - transfer_volume_step1

    elif dilution_factor > 100 and dilution_factor <= 300 and dilution_model == 'A->B':  # A-> B
        dilution_factor_step0 = 1500 // source_container_volume
        dilution_factor_step1 = dilution_factor / dilution_factor_step0
        transfer_volume_step0 = source_container_volume * (dilution_factor_step0 - 1)
        transfer_volume_step1 = final_volume / dilution_factor_step1
        transfer_volume_step2 = final_volume - transfer_volume_step1

    elif dilution_factor > 100 and dilution_factor <= 300 and dilution_model == 'A->C':  # A-> B
        dilution_factor_here = dilution_factor / source_container_diluted_factor
        transfer_volume_step0 = 100 #摆烂版本，这个数值需要后续手动修改
        transfer_volume_step1 = final_volume / dilution_factor_here
        transfer_volume_step2 = final_volume - transfer_volume_step1

    elif dilution_factor > 300 and dilution_model == 'B->C':  # serial dilution: B->C

        dilution_factor_here = dilution_factor / source_container_diluted_factor

        if dilution_factor_here <= 100:
            transfer_volume_step0 = 0
            transfer_volume_step1 = final_volume / dilution_factor_here
            transfer_volume_step2 = final_volume - transfer_volume_step1
        elif dilution_factor_here > 100 and dilution_factor_here <= 300:
            dilution_factor_step0 = 1500 // source_container_volume
            dilution_factor_step1 = dilution_factor_here / dilution_factor_step0
            transfer_volume_step0 = source_container_volume * (dilution_factor_step0 - 1)
            transfer_volume_step1 = final_volume / dilution_factor_step1
            transfer_volume_step2 = final_volume - transfer_volume_step1

    if step_id == 0:
        event.source_container = diluting_solvent_container
        event.destination_container = source_container
        event.transfer_volume = transfer_volume_step0

    elif step_id == 1:
        event.source_container = source_container
        event.destination_container = destination_container
        event.transfer_volume = transfer_volume_step1

    elif step_id == 2:
        event.source_container = diluting_solvent_container
        event.destination_container = destination_container
        event.transfer_volume = transfer_volume_step2

    if event.transfer_volume <= 50:
        event.tip_type = '50ul'
        lc_index = get_liquid_class_table_index(solvent=solvent, tip_type=event.tip_type, mode='empty')
        event.liquidClassTableIndex = lc_index
    elif event.transfer_volume <= 300 and event.transfer_volume > 50:
        event.tip_type = '300ul'
        lc_index = get_liquid_class_table_index(solvent=solvent, tip_type=event.tip_type, mode='empty')
        event.liquidClassTableIndex = lc_index
    else:
        event.tip_type = '1000ul'
        lc_index = get_liquid_class_table_index(solvent=solvent, tip_type=event.tip_type, mode='empty')
        event.liquidClassTableIndex = lc_index

    event.asp_containerGeometryTableIndex = event.source_container.containerGeometryTableIndex
    event.disp_containerGeometryTableIndex = event.destination_container.containerGeometryTableIndex
    event.asp_deckGeometryTableIndex = brb.deckGeometryTableIndex[event.tip_type]
    event.disp_deckGeometryTableIndex = event.asp_deckGeometryTableIndex

    return event

def define_tip_and_liquid_class(event):

    solvent = event.source_container.solvent

    if event.transfer_volume <= 50:
        event.tip_type = '50ul'
        lc_index = get_liquid_class_table_index(solvent=solvent, tip_type=event.tip_type, mode='empty')
        event.liquidClassTableIndex = lc_index
    elif event.transfer_volume <= 300 and event.transfer_volume > 50:
        event.tip_type = '300ul'
        lc_index = get_liquid_class_table_index(solvent=solvent, tip_type=event.tip_type, mode='empty')
        event.liquidClassTableIndex = lc_index
    else:
        event.tip_type = '1000ul'
        lc_index = get_liquid_class_table_index(solvent=solvent, tip_type=event.tip_type, mode='empty')
        event.liquidClassTableIndex = lc_index

    return event


def generate_events_for_one_plate(excel_path,
                                  plate_barcode,
                                  final_volume: float,
                                  solvent: dict,
                                  diluting_solvent_container: tuple):
    global events_A_to_B
    global events_A_to_C
    global events_B_to_C

    events_A_to_B = {0: [], 1: [], 2: []}  # {dilution_step: events}
    events_A_to_C = {0: [], 1: [], 2: []}
    events_B_to_C = {0: [], 1: [], 2: []}

    df = pd.read_excel(excel_path, sheet_name=config.sheet_name_for_run_info, engine='openpyxl')
    df_for_one_plate = df[df['plate_barcode'] == plate_barcode]

    # plate A : the plate to be diluted. It can be the plate with reaction crude or the plate with diluted products.
    # plate B or C : the target plate for dilution.
    # dilution mode: A->B, A->C, or B->C

    stock_solution_columns = [col for col in df_for_one_plate.columns if '#' in col]
    plate_A_container_volume = round(df_for_one_plate.iloc[1][stock_solution_columns].sum())

    num_dilutions = len([i for i in df_for_one_plate.columns if 'dilution_factor' in i])

    plate_A_brb_id = brb.plate_list[0]
    plate_B_brb_id = brb.plate_list[1]
    plate_C_brb_id = brb.plate_list[2]

    try:
        dilution_factor_first = df_for_one_plate.iloc[0]['dilution_factor']
    except:
        print('dilution_factor is not found in the Excel.')
        sys.exit()

    try:
        dilution_factor_second = df_for_one_plate.iloc[0]['dilution_factor_2']
    except:
        dilution_factor_second = 0

    if num_dilutions > 2:
        print(f'num_dilutions is {num_dilutions}, which is not supported yet.')
        sys.exit()

    # first dilution: A->B
    for i in range(3):
        solvent_here = solvent['reaction_crude'] if i == 1 else solvent['diluting_solvent_0']
        for index, row in df_for_one_plate.iterrows():
            index_real = index % 54
            # print(index_real)
            # print(plate_A_brb_id.containers[index_real])
            event_here = generate_events_for_one_container(excel_path=excel_path,
                                                           dilution_model = 'A->B',
                                                           df_for_one_container=row,
                                                           solvent=solvent_here,
                                                           diluting_solvent_container=diluting_solvent_container[0],
                                                           source_container=plate_A_brb_id.containers[index_real],
                                                           source_container_volume=plate_A_container_volume,
                                                           source_container_diluted_factor=1,
                                                           destination_container=plate_B_brb_id.containers[index_real],
                                                           dilution_factor=dilution_factor_first,
                                                           final_volume=final_volume,
                                                           step_id=i)

            if event_here.transfer_volume > 0:
                events_A_to_B[i].append(event_here)

    if events_A_to_B[0] == []:
        dilution_factor_for_A_so_far = 1
    else:
        dilution_factor_for_A_so_far = events_A_to_B[0][0].transfer_volume / plate_A_container_volume + 1


    if num_dilutions == 1:
        print('All events are generated.')
        return 0

    elif num_dilutions == 2:
        # second dilution: A->C or B->C
        dilution_mode = 'A->C' if dilution_factor_second <= 300 else 'B->C'

        if dilution_mode == 'A->C':
            for i in range(3):
                solvent_here = solvent['reaction_crude'] if i == 1 else solvent['diluting_solvent_1']
                for index, row in df_for_one_plate.iterrows():
                    index_real = index % 54
                    event_here = generate_events_for_one_container(excel_path=excel_path,
                                                                   dilution_model='A->C',
                                                                   df_for_one_container=row,
                                                                   solvent=solvent_here,
                                                                   diluting_solvent_container=diluting_solvent_container[1],
                                                                   source_container=plate_A_brb_id.containers[index_real],
                                                                   source_container_volume=plate_A_container_volume,
                                                                   source_container_diluted_factor=dilution_factor_for_A_so_far,
                                                                   destination_container=plate_C_brb_id.containers[index_real],
                                                                   dilution_factor=dilution_factor_second,
                                                                   final_volume=final_volume,
                                                                   step_id=i)
                    if event_here.transfer_volume > 0:
                        events_A_to_C[i].append(event_here)

        elif dilution_mode == 'B->C':
            for i in range(3):
                solvent_here = solvent['reaction_crude'] if i == 1 else solvent['diluting_solvent_1']
                for index, row in df_for_one_plate.iterrows():
                    index_real = index % 54
                    event_here = generate_events_for_one_container(excel_path=excel_path,
                                                                   dilution_model='B->C',
                                                                   df_for_one_container=row,
                                                                   solvent=solvent_here,
                                                                   diluting_solvent_container=diluting_solvent_container,
                                                                   source_container=plate_B_brb_id.containers[index_real],
                                                                   source_container_volume=plate_A_container_volume,
                                                                   source_container_diluted_factor=dilution_factor_first,
                                                                   destination_container=plate_C_brb_id.containers[index_real],
                                                                   dilution_factor=dilution_factor_second,
                                                                   final_volume=final_volume,
                                                                   step_id=i)
                    if event_here.transfer_volume > 0:
                        events_B_to_C[i].append(event_here)

        else:
            raise ValueError(f'dilution_mode is {dilution_mode}, which is not supported yet.')

        print('All events are generated.')
        return 0


def beep_n():
    duration = 600  # milliseconds
    freq = 1000  # Hz
    # time.sleep(0.2)
    for i in range(10):
        winsound.Beep(freq, duration)


def load_dilution_info_by_pysimplegui():
    source_plate_barcode = None
    destination_barcode_0 = None
    dilution_factor_0 = None
    destination_barcode_1 = None
    dilution_factor_1 = None
    num_dilutions = 0
    # make a pysimplegui window to ask which plate to dilute
    layout = [[sg.Text('Which plate do you want to dilute?')],
              [sg.InputText()],
              [sg.Submit(), sg.Cancel()]]
    window = sg.Window('Dilution', layout, size=(1000, 250), text_justification='c', font='Helvetica 25')
    event, values = window.read()

    if event == 'Submit':
        source_plate_barcode = int(values[0])
        window.close()
    elif event == 'Cancel':
        sys.exit()

    # make a window to ask how many times to dilute, choose from 1 or 2.
    layout = [[sg.Text(f'How many times to dilute for plate {source_plate_barcode}?')],
              [sg.InputCombo(['1', '2'], default_value='1')],
              [sg.Submit(), sg.Cancel()]]

    window = sg.Window('Dilution', layout, size=(1000, 250), text_justification='c', font='Helvetica 25')
    event, values = window.read()
    if event == 'Submit':
        num_dilutions = int(values[0])
        window.close()
    elif event == 'Cancel':
        sys.exit()

    if num_dilutions == 1:
        #  make a window to ask about the desination plate barcode and dilution factor
        layout = [[sg.Text(f'Dilute 1 time. What is the destination plate for plate {source_plate_barcode}?')],
                  [sg.InputText()],
                  [sg.Text('What is the dilution factor?')],
                  [sg.InputText()],
                  [sg.Submit(), sg.Cancel()]]
        window = sg.Window('Dilution', layout, size=(1000, 250), text_justification='c', font='Helvetica 25')
        event, values = window.read()
        if event == 'Submit':
            destination_barcode_0 = int(values[0])
            dilution_factor_0 = float(values[1])
            window.close()
        elif event == 'Cancel':
            sys.exit()
    if num_dilutions == 2:
        # make a window to ask the two destination plate barcodes and dilution factors
        layout = [
            [sg.Text(f'Dilute 2 times  for plate {source_plate_barcode}.')],
            [sg.Text(f'What is the first destination plate?')],
            [sg.InputText()],
            [sg.Text('What is the first dilution factor?')],
            [sg.InputText()],
            [sg.Text(f'What is the second destination plate?')],
            [sg.InputText()],
            [sg.Text('What is the second dilution factor?')],
            [sg.InputText()],
            [sg.Submit(), sg.Cancel()]
        ]
        window = sg.Window('Dilution', layout, size=(1000, 600), text_justification='c', font='Helvetica 25')
        event, values = window.read()
        if event == 'Submit':
            destination_barcode_0 = int(values[0])
            dilution_factor_0 = float(values[1])
            destination_barcode_1 = int(values[2])
            dilution_factor_1 = float(values[3])
            window.close()
        elif event == 'Cancel':
            sys.exit()

    return source_plate_barcode, destination_barcode_0, dilution_factor_0, destination_barcode_1, dilution_factor_1


def load_start_end_event_id_by_pysimplegui():
    start_event_id, end_event_id = 0, 53
    # make a window to say that all dilution events are generated. Ask for the start and end event id.
    section1 = [[sg.Text('A_to_B:')],
        [sg.Text('start_alpha(add to crude):'), sg.InputText(default_text='0', size=(5, 1))],
        [sg.Text('end_alpha:'), sg.InputText(default_text='54', size=(5, 1))],
        [sg.Text('start_beta(crude to new):'), sg.InputText(default_text='0', size=(5, 1))],
        [sg.Text('end_beta:'), sg.InputText(default_text='54', size=(5, 1))],
        [sg.Text('start_gamma(add to new):'), sg.InputText(default_text='0', size=(5, 1))],
        [sg.Text('end_gamma:'), sg.InputText(default_text='54', size=(5, 1))],
    ]
    section2 = [[sg.Text('A_to_C:')],
        [sg.Text('start_alpha(add to crude):'), sg.InputText(default_text='0', size=(5, 1))],
        [sg.Text('end_alpha:'), sg.InputText(default_text='54', size=(5, 1))],
        [sg.Text('start_beta(crude to new):'), sg.InputText(default_text='0', size=(5, 1))],
        [sg.Text('end_beta:'), sg.InputText(default_text='54', size=(5, 1))],
        [sg.Text('start_gamma(add to new):'), sg.InputText(default_text='0', size=(5, 1))],
        [sg.Text('end_gamma:'), sg.InputText(default_text='54', size=(5, 1))],
    ]
    section3 = [[sg.Text('B_to_C:')],
        [sg.Text('start_alpha(add to crude):'), sg.InputText(default_text='0', size=(5, 1))],
        [sg.Text('end_alpha:'), sg.InputText(default_text='54', size=(5, 1))],
        [sg.Text('start_beta(crude to new):'), sg.InputText(default_text='0', size=(5, 1))],
        [sg.Text('end_beta:'), sg.InputText(default_text='54', size=(5, 1))],
        [sg.Text('start_gamma(add to new):'), sg.InputText(default_text='0', size=(5, 1))],
        [sg.Text('end_gamma:'), sg.InputText(default_text='54', size=(5, 1))],
    ]

    layout = [[sg.T('Please enter the starting and ending ids.', font='Helvetica 28')],
              [sg.Column(section1), sg.VSeperator(), sg.Column(section2),  sg.VSeperator(), sg.Column(section3)],
              [sg.B('Submit', button_color=('yellow', 'firebrick4')),
               sg.B('Cancel', button_color=('yellow', 'firebrick4'))],
             ]

    window = sg.Window('Visible / Invisible Element Demo', layout, font='Helvetica 18', finalize=True, size=(1500, 800),)

    while True:  # Event Loop
        event, values = window.read()
        print(event, values)
        if event == sg.WIN_CLOSED or event == 'Submit':

            break
        elif event == 'Cancel':
            window.close()
            raise Exception('User stopped actions!')

    window.close()

    return [int(i) for i in values.values()]

def load_plate_container_id_by_pysimplegui():
    plate_id = None
    container_id = None
    # make a window to ask for the plate id and the container id
    layout = [[sg.Text('Input the diluting solvent position?', font='Helvetica 28')],
              [sg.Text('What is the plate id?')],
              [sg.InputText(default_text='6')],
              [sg.Text('What is the container id?')],
              [sg.InputText(default_text='0')],
              [sg.Submit(), sg.Cancel()]
              ]
    window = sg.Window('Dilution', layout, size=(1000, 500), text_justification='c', font='Helvetica 25')
    event, values = window.read()
    if event == 'Submit':
        plate_id = int(values[0])
        container_id = int(values[1])
        window.close()
    elif event == 'Cancel':
        sys.exit()
    return plate_id, container_id


def display_dilution_summary_by_pysimplegui(
        all_events: tuple,
        final_volume: float,
        source_container_volume: float,

        dilution_factor_0: float,
        source_plate_barcode_0:int,
        destination_plate_barcode_0: int,

        dilution_factor_1: float,
        source_plate_barcode_1:int,
        destination_plate_barcode_1:int,
):

    dilution_mode_0 = 'A->B'
    events_A_to_B, events_A_to_C, events_B_to_C = all_events[0], all_events[1], all_events[2]
    events_for_second_dilution = events_A_to_C if events_A_to_C[1]!=[] else events_B_to_C
    dilution_mode_1 = 'A->C' if events_A_to_C[1]!=[] else 'B->C'

    assert len(set([i.transfer_volume for i in events_A_to_B[0]])) <= 1, 'transfer_volume is not the same for one operation.'
    assert len(set([i.transfer_volume for i in events_A_to_B[1]])) <= 1, 'transfer_volume is not the same for one operation.'
    assert len(set([i.transfer_volume for i in events_A_to_B[2]])) <= 1, 'transfer_volume is not the same for one operation.'
    assert len(set([i.transfer_volume for i in events_A_to_C[0]])) <= 1, 'transfer_volume is not the same for one operation.'
    assert len(set([i.transfer_volume for i in events_A_to_C[1]])) <= 1, 'transfer_volume is not the same for one operation.'
    assert len(set([i.transfer_volume for i in events_A_to_C[2]])) <= 1, 'transfer_volume is not the same for one operation.'
    assert len(set([i.transfer_volume for i in events_B_to_C[0]])) <= 1, 'transfer_volume is not the same for one operation.'
    assert len(set([i.transfer_volume for i in events_B_to_C[1]])) <= 1, 'transfer_volume is not the same for one operation.'
    assert len(set([i.transfer_volume for i in events_B_to_C[2]])) <= 1, 'transfer_volume is not the same for one operation.'

    try: first_dilution_first_step_volume = events_A_to_B[0][0].transfer_volume
    except: first_dilution_first_step_volume = 0

    try: second_dilution_first_step_volume = events_for_second_dilution[0][0].transfer_volume
    except: second_dilution_first_step_volume = 0

    section1 = [
        [sg.Text('First Dilution', font='Helvetica 28')],
        [sg.Text(f'Dilution factor: {dilution_factor_0}')],
        [sg.Text(f'Dilution model: {dilution_mode_0}')],
        [sg.Text(f'Source plate barcode: {source_plate_barcode_0}')],
        [sg.Text(f'Destination plate barcode: {destination_plate_barcode_0}')],
        [sg.Text(f'Initial crude volume: {source_container_volume}')],
        [sg.Text(f'Final dilution volume: {final_volume}')],
        [sg.Text(f'First step volume: {first_dilution_first_step_volume}')],
        [sg.Text(f'Second step volume: {events_A_to_B[1][0].transfer_volume}')],
        [sg.Text(f'Third step volume: {events_A_to_B[2][0].transfer_volume}')],
    ]

    section2 = [
        [sg.Text('Second Dilution', font='Helvetica 28')],
        [sg.Text(f'Dilution factor: {dilution_factor_1}')],
        [sg.Text(f'Dilution model: {dilution_mode_1}')],
        [sg.Text(f'Source plate barcode: {source_plate_barcode_1}')],
        [sg.Text(f'Destination plate barcode: {destination_plate_barcode_1}')],
        [sg.Text(f'Initial crude volume: {source_container_volume}')],
        [sg.Text(f'Final dilution volume: {final_volume}')],
        [sg.Text(f'First step volume: {second_dilution_first_step_volume}')],
        [sg.Text(f'Second step volume: {events_for_second_dilution[1][0].transfer_volume}')],
        [sg.Text(f'Third step volume: {events_for_second_dilution[2][0].transfer_volume}')],
    ]

    layout = [[sg.T('Dilution summary, check carefully!', font='Helvetica 28')],
              [sg.Column(section1), sg.VSeperator(), sg.Column(section2)],
              [sg.B('Correct, proceed!', button_color=('yellow', 'firebrick4')),
               sg.B('Wrong, exit!', button_color=('yellow', 'firebrick4'))],
             ]

    window = sg.Window('Visible / Invisible Element Demo', layout, font='Helvetica 18', finalize=True, size=(1000, 800),)

    while True:  # Event Loop
        event, values = window.read()
        print(event, values)
        if event == sg.WIN_CLOSED or event == 'Correct, proceed!':
            break
        elif event == 'Wrong, exit!':
            window.close()
            raise Exception('User stopped actions!')

    window.close()


## yes or no to check the volume in vials
def check_volume_pysimplegui():

    font_size = 20
    custom_font = ('Any', font_size)
    layout = [
        [sg.Text('Check liquid surface height in vials?', font=custom_font)],
        [sg.Radio('Yes', 'YN', key='yes', default=True, font=custom_font),
         sg.Radio('No', 'YN', key='no', font=custom_font)],
        [sg.Button('OK', font=custom_font)],
    ]

    window = sg.Window('Yes or No Input', layout, finalize=True, size= (1000, 300))

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break
        elif event == 'OK':
            if values['yes']:
                result = True
            elif values['no']:
                result = False
            else:
                result = None
            break

    window.close()

    return result

def input_volumes_for_dilution():
    sg.theme('DarkAmber')
    section1 = [[sg.Text('First Dilution')],
                [sg.Text('Volume 0'), sg.InputText(key='VOLUME1', default_text='0')],
                [sg.Text('Volume 1'), sg.InputText(key='VOLUME2', default_text='0')],
                [sg.Text('Volume 2'), sg.InputText(key='VOLUME3', default_text='0')],]
    section2 = [[sg.Text('Second Dilution')],
                [sg.Text('Volume 0'), sg.InputText(key='VOLUME4', default_text='0')],
                [sg.Text('Volume 1'), sg.InputText(key='VOLUME5', default_text='0')],
                [sg.Text('Volume 2'), sg.InputText(key='VOLUME6', default_text='0')],]
    layout = [[sg.T('Addition volumes, check carefully!', font='Helvetica 28')],
              [sg.Column(section1), sg.VSeperator(), sg.Column(section2)],
              [sg.Button('Submit'), sg.Button('Cancel')]]
    window = sg.Window('Volume Input', layout)
    while True:
        event, values = window.read()
        if event in (None, 'Cancel'):   # if user closes window or clicks cancel
            raise TypeError("The volumes needs to be specified!")
        if event == 'Submit':
            break
    window.close()
    return [int(i) for i in values.values()]

def specify_volumes_manually(volumes:list):
    for event in events_A_to_B[0]:
        event.transfer_volume = volumes[0]
        event = define_tip_and_liquid_class(event)
    for event in events_A_to_B[1]:
        event.transfer_volume = volumes[1]
        event = define_tip_and_liquid_class(event)
    for event in events_A_to_B[2]:
        event.transfer_volume = volumes[2]
        event = define_tip_and_liquid_class(event)

    for event in events_A_to_C[0]:
        event.transfer_volume = volumes[3]
        event = define_tip_and_liquid_class(event)
    for event in events_A_to_C[1]:
        event.transfer_volume = volumes[4]
        event = define_tip_and_liquid_class(event)
    for event in events_A_to_C[2]:
        event.transfer_volume = volumes[5]
        event = define_tip_and_liquid_class(event)
if __name__ == '__main__':

    # initiate hardware
    zm, gt, pt = initiate_hardware()
    time.sleep(1)

    final_volume = 1000
    source_container_volume = 500
    # source_container_volume = 400

    # load excel path
    run_info_path = load_excel_path_by_pysimplegui()

    # load dataframe from excel
    df = pd.read_excel(run_info_path, sheet_name=config.sheet_name_for_run_info, engine='openpyxl')
    reaction_plate_barcode = set(df['plate_barcode'].to_list())

    if_check_volume_bool = check_volume_pysimplegui()

    # assign diluting solutions and reacction crude to the breadboard
    stock_solution_containers = \
        pln.assign_stock_solutions_to_containers_and_check_volume_new(excel_path=run_info_path,
                                                                  sheet_name='stock_solutions_dilution',
                                                                  check_volume_by_pipetter=if_check_volume_bool,
                                                                  pt=pt)

    reaction_crude_container, diluting_solvent_container_0, diluting_solvent_container_1 = None, None, None

    for container in stock_solution_containers:
        if container.substance == 'reaction_crude':
            reaction_crude_container = container
        elif container.substance == 'diluting_solvent_0':
            diluting_solvent_container_0 = container
        elif container.substance == 'diluting_solvent_1':
            diluting_solvent_container_1 = container

    for i in range(54):
        brb.plate0.containers[i].id = {'container_id': 0, 'plate_id': i}
        brb.plate0.containers[i].liquid_surface_height = reaction_crude_container.liquid_surface_height+10
        brb.plate0.containers[i].liquid_surface_height_real_value_float = float(brb.plate0.containers[i].liquid_surface_height)
        brb.plate0.containers[i].liquid_volume = reaction_crude_container.liquid_volume
        brb.plate0.containers[i].solvent = reaction_crude_container.solvent
        brb.plate0.containers[i].substance = reaction_crude_container.substance
        brb.plate0.containers[i].container_id = i

    # ask which plate to dilute
    source_plate_barcode, \
    destination_barcode_0, \
    dilution_factor_0, \
    destination_barcode_1, \
    dilution_factor_1 = load_dilution_info_by_pysimplegui()

    if dilution_factor_0 <= 50:
        final_volume = 500

    assert source_plate_barcode in reaction_plate_barcode, \
        f'plate {source_plate_barcode} is not found in the excel file.'
    assert dilution_factor_0 > 1, \
        f'dilution_factor_0 is {dilution_factor_0}, which is not valid.'


    # update the dilution info to the Excel file
    df.loc[df['plate_barcode'] == source_plate_barcode, 'plate_barcodes_for_dilution'] = destination_barcode_0
    df.loc[df['plate_barcode'] == source_plate_barcode, 'dilution_factor'] = dilution_factor_0
    df.loc[df['plate_barcode'] == source_plate_barcode, 'timestamp_dilution'] = 0

    if destination_barcode_1 != None:
        df.loc[df['plate_barcode'] == source_plate_barcode, 'plate_barcodes_for_dilution_2'] = destination_barcode_1
        df.loc[df['plate_barcode'] == source_plate_barcode, 'dilution_factor_2'] = dilution_factor_1
        df.loc[df['plate_barcode'] == source_plate_barcode, 'timestamp_dilution_2'] = 0

    # save the dilution info to Excel
    with pd.ExcelWriter(run_info_path, engine='openpyxl', mode='a', if_sheet_exists="replace") as writer:
        df.to_excel(writer, sheet_name=config.sheet_name_for_run_info, index=False, index_label='uuid')

    solvent_dict = {'diluting_solvent_0': diluting_solvent_container_0.solvent,
                    'diluting_solvent_1': diluting_solvent_container_1.solvent,
                    'reaction_crude': reaction_crude_container.solvent}

    generate_events_for_one_plate(excel_path=run_info_path,
                                  plate_barcode=source_plate_barcode,
                                  final_volume=final_volume,
                                  solvent=solvent_dict,
                                  diluting_solvent_container=(diluting_solvent_container_0, diluting_solvent_container_1))

    # specify volumes for each steps
    volumes = [] # there are 6 volumes, each one corresponding to one step of the three steps for the two dilutions
    # volumes = input_volumes_for_dilution()

    ## this is the actual volume that will be used.
    volumes = [500, 20, 480, 0, 100, 400]

    specify_volumes_manually(volumes)

    # display a summary of the dilution info
    display_dilution_summary_by_pysimplegui(all_events=(events_A_to_B, events_A_to_C, events_B_to_C),
                                            final_volume = final_volume,
                                            source_container_volume = source_container_volume,
                                            dilution_factor_0=dilution_factor_0,
                                            source_plate_barcode_0=source_plate_barcode,
                                            destination_plate_barcode_0=destination_barcode_0,
                                            dilution_factor_1=dilution_factor_1,
                                            source_plate_barcode_1 = source_plate_barcode,
                                            destination_plate_barcode_1=destination_barcode_1,)

    ####### RESET EVENTS FOR B TO C######
    # this is to make the second dilution from B to C. deepcopy it from events_A_to_B.
    events_B_to_C = copy.deepcopy(events_A_to_B)
    events_B_to_C[0] = [] # make the first step empty, because no dilution is needed for the first step.
    for index, event in enumerate(events_B_to_C[1]): # construct the second step for B to C
        event.source_container = brb.plate_list[1].containers[index]
        event.destination_container = brb.plate_list[2].containers[index]
        event.source_container.liquid_surface_height = 2130
        event.transfer_volume = 100
        event.liquidClassTableIndex = 13
        event.tip_type = '300ul'

    for index, event in enumerate(events_B_to_C[2]): # construct the third step for B to C
        event.source_container = brb.plate_list[6].containers[1]
        event.destination_container = brb.plate_list[2].containers[index]
        event.transfer_volume = 400
        event.tip_type = '1000ul' # 300ul tip
        event.liquidClassTableIndex = 14 # 13 is the index for 300ul tip
    ####### RESET EVENTS FOR B TO C######

    # # make events_A_to_C empty
    events_A_to_C = {0: [], 1: [], 2: []} # when diluting from B to C for the second dilution, ensure that events_A_to_C is empty by running this line.

    ##############if stop in the middle, copy paste and run the following##################################
    ids = load_start_end_event_id_by_pysimplegui()
    A_to_B_ids = {0:[ids[0],ids[1]],1:[ids[2],ids[3]],2:[ids[4],ids[5]]}
    A_to_C_ids = {0:[ids[6],ids[7]],1:[ids[8],ids[9]],2:[ids[10],ids[11]]}
    B_to_C_ids = {0:[ids[12],ids[13]],1:[ids[14],ids[15]],2:[ids[16],ids[17]]}
    print(f'ids for dilutions: {A_to_B_ids}, {A_to_B_ids}, {B_to_C_ids}')

    # first dilution: A->B
    for i in range(3):
        assert A_to_B_ids[i][0] <= A_to_B_ids[i][1], f"The input id is wrong {A_to_B_ids[i][0]}---{A_to_B_ids[i][1]}!"
        if A_to_B_ids[i][0]==0 and A_to_B_ids[i][1]==0:
            print(f"A->B step {i} is passed!")
        else:
            pln.run_events_chem_dilution_new(zm=zm, pt=pt, logger=module_logger,
                                         event_list=events_A_to_B[i],
                                         dilution_mode='A->B',
                                         step_id=i,
                                         start_event_id=A_to_B_ids[i][0],
                                         end_event_id=A_to_B_ids[i][1],
                                         log_to_excel=True,
                                         excel_path=run_info_path,
                                         barcode_of_plate_for_reactions=source_plate_barcode,
                                         prewet_tip=True)
        module_logger.info(f'First dilution step {i} is done.')

    beep_n()
    module_logger.info('The first dilution is done.')

    # second dilution: A->C or B->C
    if events_A_to_C[1] != []:
        for i in range(3):
            assert A_to_C_ids[i][0] <= A_to_C_ids[i][1], f"The input id is wrong {A_to_C_ids[i][0]}---{A_to_C_ids[i][1]}!"
            if A_to_C_ids[i][0] == 0 and A_to_C_ids[i][1] == 0:
                print(f"A->C step {i} is passed!")
            else:
                pln.run_events_chem_dilution_new(zm=zm, pt=pt, logger=module_logger,
                                         event_list=events_A_to_C[i],
                                         dilution_mode='A->C',
                                         step_id=i,
                                         start_event_id=A_to_C_ids[i][0],
                                         end_event_id=A_to_C_ids[i][1],
                                         log_to_excel=True,
                                         excel_path=run_info_path,
                                         barcode_of_plate_for_reactions=source_plate_barcode,
                                         prewet_tip=True )
            module_logger.info(f'Second dilution step {i} is done.')
        beep_n()

    elif events_B_to_C[1] != []:
        ####### RUN ONLY B TO C######
        for i in range(3):
            assert B_to_C_ids[i][0] <= B_to_C_ids[i][1], f"The input id is wrong {B_to_C_ids[i][0]}---{B_to_C_ids[i][1]}!"
            if B_to_C_ids[i][0] == 0 and B_to_C_ids[i][1] == 0:
                print(f"B->C step {i} is passed!")
            else:
                pln.run_events_chem_dilution_new(zm=zm, pt=pt, logger=module_logger,
                                         event_list=events_B_to_C[i],
                                         dilution_mode='B->C',
                                         step_id=i,
                                         start_event_id=B_to_C_ids[i][0],
                                         end_event_id=B_to_C_ids[i][1],
                                         log_to_excel=True,
                                         excel_path=run_info_path,
                                         barcode_of_plate_for_reactions=source_plate_barcode,
                                         purpose = 'diluting',prewet_tip=True)
            module_logger.info(f'Second dilution step {i} is done.')
        beep_n()
        module_logger.info('All dilution(s) are done.')
        ####### RUN ONLY B TO C######

    else:
        beep_n()
        module_logger.info('No second dilution is needed. All dilutions are done.')
    ############################################################################################################
