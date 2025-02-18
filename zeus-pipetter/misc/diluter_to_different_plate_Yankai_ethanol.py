import logging

import winsound

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
    logger = logging.getLogger('main')
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler('C:\\Yankai\\Dropbox\\robochem\\pipetter_files\\main_roboski2.log')
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

import copy, time, pickle, re, importlib, json, os, sys, PySimpleGUI as sg, pandas as pd
import zeus, pipetter, planner as pln, breadboard as brb, config

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

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
    window = sg.Window('Excel file', layout, size=(1000, 300), text_justification='c', font = 'Helvetica 25')
    event, values = window.read()
    window.close()
    excel_path = values[0]
    module_logger.info(f'Excel file path: {excel_path}')
    print(f'Excel file path: {excel_path}')

    return excel_path

def check_plate_barcodes_for_dilution(run_info_path: str):

    df_run_info = pd.read_excel(run_info_path, sheet_name='reactions_with_run_info', engine='openpyxl')

    # make a pysimlegui window to input the two barcodes of plates
    layout = [
            [sg.Text('Input barcode of plate with reactions: slot 1')],
            [sg.Text('Barcode'), sg.InputText(size=(10, 10))],
            [sg.Text('Input barcode of plate for dilution: slot 2')],
            [sg.Text('Barcode'), sg.InputText(size=(10, 10))],
            [sg.Submit(), sg.Cancel()]
            ]
    window = sg.Window('Plate barcodes', layout,  size=(1000, 300), text_justification='c', font = 'Helvetica 25')
    event, values = window.read()
    window.close()

    barcode_of_plate_for_reactions = int(values[0])
    barcode_of_plate_for_dilution= int(values[1])

    print(f'barcode_of_plate_for_reactions: {barcode_of_plate_for_reactions}')
    print(f'barcode_of_plate_for_dilution you just specified: {barcode_of_plate_for_dilution}')

    barcode_of_plate_for_reactions_from_df = 0
    # get the barcodes of plates from the dataframe
    for index, row in df_run_info.iterrows():
        if row['plate_barcode'] == barcode_of_plate_for_reactions:
            barcode_of_plate_for_reactions_from_df = row['plate_barcodes_for_dilution']
            print(f'barcode_of_plate_for_reactions_from_df: {barcode_of_plate_for_reactions_from_df}')
            break
        ## if the barcode is not found in the dataframe, exit the program
        elif index == len(df_run_info) - 1:
            print(f'barcode_of_plate_for_reactions: {barcode_of_plate_for_reactions} is not found in the dataframe.\n'
                  f'Please check the Excel file.')
            module_logger.info(f'barcode_of_plate_for_reactions: {barcode_of_plate_for_reactions} is not found in the dataframe')
            raise Exception(f'barcode_of_plate_for_reactions: {barcode_of_plate_for_reactions} is not found in the dataframe')

    if barcode_of_plate_for_reactions_from_df == barcode_of_plate_for_dilution:
        ## make a pysimplegui window to tell the user the barcodes are correct
        layout = [[sg.Text('The barcodes are correct. Please continue')],[sg.Submit(), sg.Cancel()]]
        window = sg.Window('Plate barcodes', layout, size=(1000, 300), text_justification='c', font='Helvetica 25')
        event, values = window.read()
        window.close()

        return barcode_of_plate_for_reactions, barcode_of_plate_for_dilution

    else:
        ## make a pysimplegui window to tell the user the barcodes are incorrect ana ask for options from three buttons
        layout = [
            [sg.Text('The barcodes do NOT match. Please proceed with one of the following options')],
            [sg.Button('1. I will change the breadboard plate to match the Excel.'),],
            [sg.Button('2. I insist using the breadboard plate. Overwrite the plate barcode in the Excel.')],
            [sg.Button('3. Something is wrong, I will exit the program')]
        ]

        window = sg.Window('Plate barcodes', layout, size=(1200, 400), text_justification='c', font='Helvetica 25')
        event, values = window.read()
        window.close()
        ## check what the user choose
        if event == '1. I will change the breadboard plate to match the Excel.':
            ## make a pysimplegui window to input OK for changing the plate on the breadboard
            layout = [[sg.Text('Please change the plate on the breadboard and click OK')], [sg.Submit(), sg.Cancel()]]
            window = sg.Window('Plate barcodes', layout, size=(1000, 250), text_justification='c', font='Helvetica 25')
            event, values = window.read()
            window.close()

            return barcode_of_plate_for_reactions, barcode_of_plate_for_reactions_from_df

        elif event == '2. I insist using the breadboard plate. Overwrite the plate barcode in the Excel.':
            ## overwrite the barcode in the dataframe and save back to the Excel
            df_run_info.loc[df_run_info['plate_barcodes_for_dilution'] == barcode_of_plate_for_reactions_from_df, 'plate_barcodes_for_dilution'] = barcode_of_plate_for_dilution

            with pd.ExcelWriter(run_info_path, engine='openpyxl', mode='a', if_sheet_exists="replace") as writer:
                df_run_info.to_excel(writer, sheet_name= config.sheet_name_for_run_info, index=False)

            print(f'barcode of the plate for dilution is overwritten by user: from {barcode_of_plate_for_reactions_from_df} to {barcode_of_plate_for_dilution}')
            module_logger.info(f'barcode of the plate for dilution is overwritten by user: from {barcode_of_plate_for_reactions_from_df} to {barcode_of_plate_for_dilution}')

            return barcode_of_plate_for_reactions, barcode_of_plate_for_dilution

        elif event == '3. Something is wrong, stop the program.':
            raise Exception('Something is wrong, stop the program.')

##  generate_dilution_events()
# step1: dilution original reactions, adding volume: 1400ul
# step2: transfer liquid from original reaction to new vial, transfer volume: 20ul
# step3: dilution new vial, adding volume: 480ul

# volume setting versions:
# 1. 200ul, 1600ul, 22.5ul, 477.5ul
# 2. 200ul, 1400ul, 20ul, 480ul
# 3. 200ul, 1000ul, 15ul, 485ul ## this one is now used, 2023-03-22 14:23

def get_liquid_class_table_index(solvent: str,  tip_type: str, mode: str = 'empty') -> dict:
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

def generate_dilution_event(source_container: object,
                            destination_container: object,
                            volume: float,
                            solvent:str):

    ## load a Event object as template
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


    # asign new parameters to the event
    event.source_container = source_container
    event.destination_container = destination_container
    event.transfer_volume = volume
    event.event_label = f'Dilution: transfer from {source_container.container_id} ' \
                        f'to {destination_container.container_id}'

    if volume <= 50:
        event.tip_type = '50ul'
        lc_index = get_liquid_class_table_index(solvent=solvent, tip_type=event.tip_type, mode='empty')
        event.liquidClassTableIndex = lc_index
    elif volume <= 300 and volume > 50:
        event.tip_type = '300ul'
        lc_index = get_liquid_class_table_index(solvent=solvent, tip_type=event.tip_type, mode='empty')
        event.liquidClassTableIndex = lc_index
    else:
        event.tip_type = '1000ul'
        lc_index = get_liquid_class_table_index(solvent=solvent, tip_type=event.tip_type, mode='empty')
        event.liquidClassTableIndex = lc_index

    event.asp_lld = 0
    event.disp_lld = 0

    event.asp_containerGeometryTableIndex = source_container.containerGeometryTableIndex
    event.disp_containerGeometryTableIndex = destination_container.containerGeometryTableIndex

    event.condition_uuid = "Dilution"

    return event


def generate_events_for_diluting_old_vial(solvent: str, volume_added_to_old_vial: float, rows_to_dilute: tuple=(0, 9, 18, 27, 36, 45)): # diluting volume 1400ul
    event_list_dilute_old_vial = []
    # generate dilution events
    source_container = brb.plate_list[6].containers[0]
    source_container.liquid_surface_height, source_container.liquid_volume \
    = pt.check_volume_in_container(container=source_container)

    for i in rows_to_dilute:
        for vial_index in range(i, i+9):
            destination_container = brb.plate_list[0].containers[vial_index]
            destination_container.liquid_surface_height = 2100
            event_temp = generate_dilution_event(source_container=source_container,
                                                destination_container=destination_container,
                                                volume=volume_added_to_old_vial,
                                                solvent=solvent)
            event_list_dilute_old_vial.append(event_temp)
    # time.sleep(2)


    return event_list_dilute_old_vial


# step2: transfer liquid from original reaction to new vial, transfer volume: 15ul
def generate_events_for_transferring_liquid_from_old_vials_to_new(solvent: str, volume_transfered_from_old_to_new_vial: float): # transfer volume 20ul
    event_list_dilution_old_to_new = []
    for vial_index in range(54):
        source_container = brb.plate_list[0].containers[vial_index]
        source_container.liquid_surface_height = 2150
        destination_container = brb.plate_list[1].containers[vial_index]
        destination_container.liquid_surface_height = 2100
        event_temp = generate_dilution_event(source_container=source_container,
                                             destination_container=destination_container,
                                             volume=volume_transfered_from_old_to_new_vial,
                                             solvent=solvent)
        event_list_dilution_old_to_new.append(event_temp)

    return event_list_dilution_old_to_new

# step3: dilution new vial, adding volume: 485ul
def generate_events_for_diluting_new_vial(solvent:str, volume_added_to_new_vial: float): # diluting volume 485ul

    event_list_dilute_new_vial = []
    source_container = brb.plate_list[6].containers[0]
    source_container.liquid_surface_height, source_container.liquid_volume \
    = pt.check_volume_in_container(container=source_container)
    for vial_index in range(54):
        destination_container = brb.plate_list[1].containers[vial_index]
        destination_container.liquid_surface_height = 2100
        event_temp = generate_dilution_event(source_container=source_container,
                                             destination_container=destination_container,
                                             volume=volume_added_to_new_vial,
                                             solvent=solvent)
        event_list_dilute_new_vial.append(event_temp)

    return event_list_dilute_new_vial
def dilute(zm, pt, dilution_factor: int, volume_of_reaction: int, final_volume_of_dilution: int, first_step_dilution_factor: int = 3):

    volume_added_to_old_vial = volume_of_reaction * (first_step_dilution_factor-1)
    volume_transfered_from_old_to_new_vial = final_volume_of_dilution / ((dilution_factor/first_step_dilution_factor)-1)
    volume_added_to_new_vial = final_volume_of_dilution - volume_transfered_from_old_to_new_vial

    event_list_dilute_old_vial = generate_events_for_diluting_old_vial(solvent=solvent,
                                                                       volume_added_to_old_vial=volume_added_to_old_vial)
    event_list_dilution_old_to_new = generate_events_for_transferring_liquid_from_old_vials_to_new(solvent=solvent,
                                                                                                   volume_transfered_from_old_to_new_vial=volume_transfered_from_old_to_new_vial)
    event_list_dilute_new_vial = generate_events_for_diluting_new_vial(solvent=solvent,
                                                                       volume_added_to_new_vial=volume_added_to_new_vial)
    #step1
    pln.run_events_chem_dilution(zm=zm, pt=pt, logger=module_logger, event_list= event_list_dilute_old_vial,
                                 start_event_id= 0)
    #step2
    pln.run_events_chem_dilution(zm=zm, pt=pt, logger=module_logger, event_list=event_list_dilution_old_to_new,
                                 start_event_id=0, change_tip_after_every_pipetting=True)
    # step3
    pln.run_events_chem_dilution(zm=zm, pt=pt, logger=module_logger, event_list=event_list_dilute_new_vial,
                                    start_event_id=0, log_to_excel=True, excel_path=run_info_path,
                                    barcode_of_plate_for_reactions = barcode_of_plate_for_reactions,)

    return event_list_dilute_old_vial, event_list_dilution_old_to_new, event_list_dilute_new_vial




def beep_n():
    duration = 600  # milliseconds
    freq = 1000  # Hz
    # time.sleep(0.2)
    for i in range(10):
        winsound.Beep(freq, duration)

if __name__ == '__main__':

    # specify volumes for dilution

    ## for multicomponent reaction, volume of reaction was 200 ul, and the final dilution factor was 200X.
    # volume_added_to_old_vial = 1000
    # volume_transfered_from_old_to_new_vial = 15
    # volume_added_to_new_vial = 485

    # for SN1,WSWC001 volume of reaction was 500 ul, and the final dilution factor was 100X.
    # volume_added_to_old_vial = 1000
    # volume_transfered_from_old_to_new_vial = 15
    # volume_added_to_new_vial = 485

    # for SN1,WSWC001 volume of reaction was 500 ul, and the final dilution factor was 100X.
    volume_added_to_old_vial = 1000
    volume_transfered_from_old_to_new_vial = 40
    volume_added_to_new_vial = 960

    solvent = 'ethanol'

    ## for Hantzsch volume of reaction was 500 ul, and the final dilution factor was 200X.
    # volume_added_to_old_vial = 1000
    # volume_transfered_from_old_to_new_vial = 15
    # volume_added_to_new_vial = 985

    ## for E1 , volume of reaction is 500 ul, and the final dilution factor is 30X
    # volume_added_to_old_vial = 1000
    # volume_transfered_from_old_to_new_vial = 40
    # volume_added_to_new_vial = 360

    # specify solvent
    ## for SN1, Dioxane was the solvent
    # solvent = 'Dioxane'

    ## for E1, MeCN was the solvent. However, ethanol is used in the code because their LC paras are almost the same.
    # solvent = 'ethanol'

    # load excel path
    run_info_path = load_excel_path_by_pysimplegui()

    # load run info by dataframe
    df_run_info = pd.read_excel(run_info_path, sheet_name='reactions_with_run_info', engine='openpyxl')

    # check if proper plates are placed in the breadboard
    barcode_of_plate_for_reactions, barcode_of_plate_for_dilution = \
    check_plate_barcodes_for_dilution(run_info_path = run_info_path)
    #
    ## initiate hardware
    zm, gt, pt = initiate_hardware()
    time.sleep(1)


    event_list_dilute_old_vial = \
        generate_events_for_diluting_old_vial(solvent=solvent,volume_added_to_old_vial=volume_added_to_old_vial)

    event_list_dilution_old_to_new = \
        generate_events_for_transferring_liquid_from_old_vials_to_new(solvent=solvent,volume_transfered_from_old_to_new_vial=volume_transfered_from_old_to_new_vial)

    event_list_dilute_new_vial = \
        generate_events_for_diluting_new_vial(solvent=solvent,volume_added_to_new_vial=volume_added_to_new_vial)

    # step1
    pln.run_events_chem_dilution(zm=zm, pt=pt, logger=module_logger, event_list= event_list_dilute_old_vial,
                                 start_event_id= 0, prewet_tip= False)
    # # step2
    pln.run_events_chem_dilution(zm=zm, pt=pt, logger=module_logger, event_list=event_list_dilution_old_to_new,
                                 start_event_id=0, change_tip_after_every_pipetting=True)
    # step3
    pln.run_events_chem_dilution(zm=zm, pt=pt, logger=module_logger, event_list=event_list_dilute_new_vial,
                                    start_event_id=0, log_to_excel=True, excel_path=run_info_path,
                                    barcode_of_plate_for_reactions = barcode_of_plate_for_reactions,)