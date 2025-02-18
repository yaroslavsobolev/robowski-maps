import csv
import logging, copy, asyncio, PySimpleGUI as sg
import sys

import pyautogui, json
pyautogui.FAILSAFE = False
import config
import pandas as pd

# create logger
module_logger = logging.getLogger('main.measure_spectrum')

import zeus, pipetter, planner as pln, sound, breadboard as brb, prepare_reaction as prep, nanodrop
import time, os, pickle

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
time_for_measuring_one_spectrum = 5

def initiate_hardware() -> (zeus.ZeusModule, pipetter.Gantry, pipetter.Pipetter):
    # initiate zeus
    zm = zeus.ZeusModule(id=1)
    time.sleep(1)
    module_logger.info("zeus is loaded as: zm")

    # initiate gantry
    gt = pipetter.Gantry(zeus=zm)
    time.sleep(1)

    module_logger.info("gantry is loaded as: gt")
    # gt.configure_grbl() # This only need to be done once.
    gt.home_xy()
    if gt.xy_position == (0, 0):
        module_logger.info("gantry is homed")

    # initiate pipetter
    pt = pipetter.Pipetter(zeus=zm, gantry=gt)
    time.sleep(1)
    module_logger.info("pipetter is loaded as: pt")

    return zm, gt, pt

def construct_liquid_transfer_events_for_measurement_old():
    # mn.data_folder
    with open(data_folder + 'pipetter_files\\event_template.pickle', 'rb') as f:
        event = pickle.load(f)

    event_list = []
    source_plate = brb.plate2

    for index, container in enumerate(source_plate.containers):
        container.liquid_surface_height = 2160
        container.liquid_volume = 1000

        event_here = copy.deepcopy(event)
        event_here.source_container = container
        event_here.destination_container = brb.nanodrop_pedestal

        event_here.transfer_volume = 6 # volume 4 ul generally,including ethanol, water
        # event_here.transfer_volume = 6 # volume 6 ul for DCM
        # event_here.transfer_volume = 2 # volume 3 ul for 1,2-DCE
        # event_here.transfer_volume = 30 # volume 30 ul for 1,2-DCE with O-ring

        event_here.lld = 0
        event_here.tip_type = '50ul'
        event_here.liquidClassTableIndex = 40 ## LC 40 is only for nanodrop, it is based on 27 (dioxane)
        event_here.asp_containerGeometryTableIndex = 0 # this is for 2-ml vial
        event_here.disp_containerGeometryTableIndex = 6 # this is for nanodrop pedestal
        event_here.event_label = f'{index}'
        event_list.append(event_here)

    return event_list

def construct_liquid_transfer_events_for_measurement(df: pd.DataFrame) -> list:
    # mn.data_folder
    with open(data_folder + 'pipetter_files\\event_template.pickle', 'rb') as f:
        event = pickle.load(f)

    event_list = []
    source_plate = brb.plate2


    for index, container in enumerate(source_plate.containers):
        if index <= len(df)-1: # only work on non-empty rows
            # print(f'index #{index}')
            row_for_this_container = df[df['local_index'] % 54 == index]
            # print(f'row for this container: {row_for_this_container}')
            container.liquid_surface_height = 2160 # in 0.1 mm
            container.liquid_volume = 1000

            event_here = copy.deepcopy(event)
            event_here.source_container = container
            event_here.destination_container = brb.nanodrop_pedestal

            event_here.transfer_volume = 6 # volume 4 ul generally,including ethanol, water
            # event_here.transfer_volume = 6 # volume 6 ul for DCM
            # event_here.transfer_volume = 2 # volume 3 ul for 1,2-DCE
            # event_here.transfer_volume = 30 # volume 30 ul for 1,2-DCE with O-ring

            event_here.lld = 0
            event_here.tip_type = '50ul'
            event_here.liquidClassTableIndex = 40 ## LC 40 is only for nanodrop, it is based on 27 (dioxane)
            event_here.asp_containerGeometryTableIndex = 0 # this is for 2-ml vial
            event_here.disp_containerGeometryTableIndex = 6 # this is for nanodrop pedestal
            print(f'row for this container: {row_for_this_container}')
            event_here.condition_uuid = f'{row_for_this_container["uuid"].values[0]}' # this is for the event id
            event_here.event_label = f'{index}_{row_for_this_container["uuid"].values[0]}'
            event_list.append(event_here)
        else:
            print(f'Vial #{index} is empty!')

    return event_list

def run_one_event_chem(pt: object, event=None):

    nd.open_lid()
    time.sleep(0.5)

    if zm.tip_on_zeus=='':
        pt.pick_tip('50ul')
    else:
        pt.discard_tip()
        pt.pick_tip('50ul')

    pt.transfer_liquid(event)
    gt.move_to_idle_position()
    time.sleep(0.2)
    nd.close_lid()
    time.sleep(0.1)
    pt.discard_tip()
    gt.move_to_idle_position()
    time.sleep(0.2)

class ZeusError(Exception):
    pass

async def pick_tip():

    print("picking up tips...")

    if zm.tip_on_zeus=='':
        pt.pick_tip('50ul')
    elif zm.tip_on_zeus != '':
        pt.discard_tip()
        time.sleep(0.1)
        pt.pick_tip('50ul')

async def aspirate_next_sample(event=None):

    assert event != None, "Event for asp is incorrect!!"

    zm.move_z(880)
    gt.move_xy(event.source_container.xy)

    await asyncio.sleep(0.1)

    print(f'aspirating liquid... {event.source_container.id}')

    try:
        # print(f'Aspiration volume for zeus: {int(round(transfer_event.aspirationVolume * 10))}')
        zm.aspiration(aspirationVolume=int(round(event.transfer_volume*5*10)),  # volume in 0.1 ul,## transfer_volume is doubled, this is to avoid bubbles.
                     containerGeometryTableIndex=event.asp_containerGeometryTableIndex,
                     deckGeometryTableIndex=event.asp_deckGeometryTableIndex,
                     liquidClassTableIndex=event.liquidClassTableIndex,
                     qpm=event.asp_qpm,
                     lld=event.asp_lld,
                     liquidSurface=event.source_container.liquid_surface_height,
                     lldSearchPosition=event.source_container.liquid_surface_height - 50,
                     mixVolume=event.asp_mixVolume,
                     mixFlowRate=event.asp_mixFlowRate,
                     mixCycles=event.asp_mixCycles)

        zm.wait_until_zeus_responds_with_string('GAid')

    except ZeusError:
        if zm.zeus_error_code(zm.r.received_msg) == '81':
            # Empty tube detected during aspiration
            gt.logger.info('ZEUS ERROR: Empty tube during aspiration. Dispensing and trying again.')
            time.sleep(2)
            exit()

    await asyncio.sleep(0.1)

    gt.move_to_idle_position()

    await asyncio.sleep(0.1)

def dispense_sample(event=None):
    try:
        pt.dispense_liquid(event)
    except:
        gt.move_to_idle_position()
        print("Pipetting error happened!")
        sound.beep_for_tip_changing()
        sound.beep_for_tip_changing()
        time.sleep(1)

async def measure_one_spectrum_by_pyautogui():

    ## activate nanodrop software window
    pyautogui.moveTo(2700, 510)  # Find where button.png appears on the screen and click it.
    pyautogui.click()

    # start measurement
    pyautogui.moveTo(2604, 111)  # Find where button.png appears on the screen and click it.
    pyautogui.click()
    print('measuring spectrum')
    print(f'timestampe here: {time.strftime("%Y-%m-%d %H:%M:%S")}')

async def input_sample_name(sample_name: str):

    ## activate nanodrop software window
    pyautogui.moveTo(2700, 510)  # Find where button.png appears on the screen and click it.
    pyautogui.click()
    # input sample name
    pyautogui.moveTo(4943, 168)
    pyautogui.click()
    pyautogui.press('backspace', presses=26)
    pyautogui.press('delete', presses=10)
    pyautogui.write(sample_name, interval=0.0)

async def main(events_for_measurement = None, only_do_ids = ()):

    assert events_for_measurement != None, "event list is empty!"

    # if only_do_ids is not empty, only do those events
    if len(only_do_ids) > 0: events_for_measurement= [events_for_measurement[i] for i in only_do_ids]

    # start with cleaning the pedestal
    print('cleaning pedestal...')
    # nd.flush_pedestal_not_async()
    # nd.dry_pedestal_not_async()

    await asyncio.gather(pick_tip())
    await asyncio.gather(aspirate_next_sample(events_for_measurement[0]),
                         input_sample_name(events_for_measurement[0].event_label))

    for num, event in enumerate(events_for_measurement):

        nd.open_lid()
        time.sleep(1)
        nd.close_air()
        dispense_sample(event)
        gt.move_to_idle_position()
        nd.close_lid()
        time.sleep(4) # settle down time for solvent

        if num != len(events_for_measurement)-1:
            await asyncio.gather(measure_one_spectrum_by_pyautogui())
            pt.discard_tip()
            time.sleep(5)
            await asyncio.gather(nd.flush_pedestal(),)
            await asyncio.gather(input_sample_name(events_for_measurement[num+1].event_label),
                                 pick_tip())
            await asyncio.gather(nd.dry_pedestal(),
                                 aspirate_next_sample(events_for_measurement[num+1]))
        elif num == len(events_for_measurement)-1:
            await asyncio.gather(measure_one_spectrum_by_pyautogui())

        sound.beep_for_measurement()

    time.sleep(5)
    print('cleaning pedestal...')
    wash_nanodrop_pedestal()

    sound.beep_for_success()
    print("all measurements are done!!")
    nd.close_all()

def gt_move_to_pedestal():
    pd_xy = brb.nanodrop_pedestal.xy
    nd.open_lid()
    time.sleep(2)
    gt.move_xy(pd_xy)

def load_excel_path_by_pysimplegui():
    # change the theme of the GUI
    sg.theme('DarkAmber')
    layout = [
        [sg.Text('Please select the Excel file with run info for dilution')],
        [sg.Text('Excel file'), sg.Input(size=(30, 10)), sg.FileBrowse()],
        [sg.Submit(), sg.Cancel()]
    ]
    window = sg.Window('Excel file', layout, size=(1400, 300), text_justification='c', font='Helvetica 25')
    event, values = window.read()
    window.close()
    excel_path = values[0]
    module_logger.info(f'Excel file path: {excel_path}')
    print(f'Excel file path: {excel_path}')

    return excel_path

def load_plate_by_pysimplegui():
    # change the theme of the GUI
    sg.theme('DarkAmber')
    layout = [
        [sg.Text('Please input the plate barcode to be measured')],
        [sg.Input(size=(30, 10))],
        [sg.Submit(), sg.Cancel()]
    ]
    window = sg.Window('Plate barcode', layout, size=(1400, 300), text_justification='c', font='Helvetica 25')
    event, values = window.read()
    window.close()
    plate_barcode = values[0]
    module_logger.info(f'Plate to be measured: {plate_barcode}')
    print(f'Plate to be measured: {plate_barcode}')

    return int(plate_barcode)

def if_wash_nanodrop_pedestal():
    sg.theme('DarkAmber')

    result = None
    custom_font = ('Any', 20)
    layout = [[sg.Text('Wash Nanodrop pedestal?', font=custom_font)],
            [sg.Radio('Yes', 'YN', key='yes', default=True, font=custom_font),sg.Radio('No', 'YN', key='no', font=custom_font)],
            [sg.Button('OK', font=custom_font)]]

    window = sg.Window('Yes or No Input', layout, finalize=True, size=(1000, 300))

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED: break
        elif event == 'OK':
            if values['yes']: result = True
            elif values['no']:result = False
            else:result = None
            break

    window.close()

    return result
    nd.flush_pedestal_not_async()

def wash_nanodrop_pedestal():

    time.sleep(0.5)
    nd.flush_pedestal_not_async()
    nd.dry_pedestal_not_async()
    time.sleep(0.5)
    nd.dry_pedestal_not_async()

def yes_or_no(title:str):
    sg.theme('DarkAmber')
    layout = [[sg.Text(title)],[sg.Button('Yes'), sg.Button('No')]]

    window = sg.Window('Yes or No', layout, size=(1000, 300), font='Helvetica 20')

    while True:
        event, values = window.read()

        if event == sg.WIN_CLOSED or event == 'No':
            choice = False
            break
        elif event == 'Yes':
            choice = True
            break

    window.close()

    return choice

def load_start_end_event_id_by_pysimplegui():

    # make a window to say that all dilution events are generated. Ask for the start and end event id.
    layout = [
        [sg.Text('Please enter the starting and ending id.')],
        [sg.Text('Starting id:'), sg.InputText(default_text='0')],
        [sg.Text('Ending id:'), sg.InputText(default_text='53')],
        [sg.Text('Only do ids:'), sg.InputText(default_text='')],
        [sg.Submit(), sg.Cancel()]
    ]

    window = sg.Window('Dilution', layout, size=(1000, 500), text_justification='c', font='Helvetica 25')
    event, values = window.read()
    if event == 'Submit':
        if values[2] != '':
            only_do_ids = tuple([int(x) for x in values[2].split(',')])
            start_event_id = int(values[0])
            end_event_id = int(values[1])

        elif values[2] == '':
            only_do_ids = tuple()
            start_event_id = int(values[0])
            end_event_id = int(values[1])

        window.close()

    elif event == 'Cancel':
        window.close()

    return start_event_id, end_event_id, only_do_ids

if __name__ == '__main__':

    # ## this is for nd_2000c
    # for event in events_for_measurement:
    #     event.destination_container.xy = (-600, -150)

    nd = nanodrop.Nanodrop(nd_id='nd_2000c')
    nd.close_lid()
    time.sleep(2)

    zm, gt, pt = initiate_hardware()

    if yes_or_no(title='Do you need spectrum match to UUID?'):

        reaction_excel = load_excel_path_by_pysimplegui()
        df = pd.read_excel(reaction_excel, sheet_name=config.sheet_name_for_run_info, engine='openpyxl')

        plate_to_be_measured = load_plate_by_pysimplegui()

        if plate_to_be_measured in [int(x) for x in df['plate_barcodes_for_dilution'].values if str(x) != 'nan']:
            dilution_mode = 'plate_barcodes_for_dilution'
        elif plate_to_be_measured in [int(x) for x  in df['plate_barcodes_for_dilution_2'].values if str(x) != 'nan']:
            dilution_mode = 'plate_barcodes_for_dilution_2'
        else:
            raise ValueError('The plate barcode is not in the Excel file!')

        df_this_plate = df[df[dilution_mode] == plate_to_be_measured]

        events_for_measurement = \
            construct_liquid_transfer_events_for_measurement(df=df_this_plate)

    else:
        events_for_measurement = construct_liquid_transfer_events_for_measurement_old()

    if if_wash_nanodrop_pedestal():
        wash_nanodrop_pedestal()

    if yes_or_no(title="Have you prepared the software and done blanking?"):
        start_event_id, end_event_id, only_do_ids = load_start_end_event_id_by_pysimplegui()

        only_do_ids_here = []

        if len(only_do_ids) == 0:
            only_do_ids_here = tuple(range(start_event_id, end_event_id + 1))
        elif len(only_do_ids) > 0:
            only_do_ids_here = only_do_ids

        asyncio.run(main(events_for_measurement,
                         only_do_ids=only_do_ids_here))
    else:
        raise ValueError("You didn't do blanking!")

