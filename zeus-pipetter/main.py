"""
workflow:
1. initiate hardware
2. generate event list for detecting surface height of stock solutions
3. run events for surface detection, get liquid surface heights and write to excel
4. generate event list for pipetting
"""
import logging, copy, time, pickle, re, importlib, json, os, PySimpleGUI as sg, pandas as pd, numpy as np

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
    logger = logging.getLogger('main')
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

logger = setup_logger()


import zeus, pipetter, planner as pln, breadboard as brb, prepare_reaction as prep

# import matplotlib
# matplotlib.use('TKAgg')


def initiate_hardware() -> (zeus.ZeusModule, pipetter.Gantry, pipetter.Pipetter):
    # initiate zeus
    zm = zeus.ZeusModule(id=1)
    time.sleep(3)
    logger.info("zeus is loaded as: zm")

    # initiate gantry
    gt = pipetter.Gantry(zeus=zm)
    time.sleep(3)
    logger.info("gantry is loaded as: gt")
    # gt.configure_grbl() # This only need to be done once.
    gt.home_xy()
    if gt.xy_position == (0, 0):
        logger.info("gantry is homed")

    # initiate pipetter
    pt = pipetter.Pipetter(zeus=zm, gantry=gt)
    time.sleep(2)
    logger.info("pipetter is loaded as: pt")

    return zm, gt, pt


def assign_stock_solutions_to_containers_and_check_volume(excel_path:str, check_volume_by_pipetter: bool = True):
    sheet_name_for_stock_solutions = 'stock_solutions'
    stock_solution_containers = []
    ## load stock solutions from excel
    df_stock_solutions = pd.read_excel(excel_path, sheet_name=sheet_name_for_stock_solutions, engine='openpyxl')

    for index, row in df_stock_solutions.iterrows():
    # ## print out the column names
    #     columns = df_stock_solutions.columns
    #     print(f"columns: {columns}")
    ## assign stock solutions to containers
        substance_name, \
        plate_id,\
        container_id, \
        solvent, \
        density, \
        volume_ml, \
        liquid_surface_height,\
        mode= row['substance'], row['breadboard_plate_id'], row['container_id'], row['solvent'], row['density'], row['volume'], row['liquid_surface_height'], row['pipetting_mode']
        container_for_this_stock_solution = brb.plate_list[plate_id].containers[container_id]
        container_for_this_stock_solution.substance = substance_name
        container_for_this_stock_solution.substance_density = density
        container_for_this_stock_solution.solvent = solvent
        container_for_this_stock_solution.pipetting_mode = mode
        if not check_volume_by_pipetter:
            container_for_this_stock_solution.liquid_surface_height = liquid_surface_height
            container_for_this_stock_solution.liquid_volume = volume_ml
        else:
            container_for_this_stock_solution.liquid_surface_height, \
            container_for_this_stock_solution.liquid_volume = \
                pt.check_volume_in_container(container=container_for_this_stock_solution,
                                             liquidClassTableIndex=12,change_tip_after_each_check=True)

        stock_solution_containers.append(container_for_this_stock_solution)

        ## update the excel after checking the volume
        if check_volume_by_pipetter:
            df_stock_solutions.loc[index, 'volume'] = container_for_this_stock_solution.liquid_volume
            df_stock_solutions.loc[index, 'liquid_surface_height'] = container_for_this_stock_solution.liquid_surface_height
            with pd.ExcelWriter(excel_path, engine='openpyxl', mode='a', if_sheet_exists="replace") as writer:
                df_stock_solutions.to_excel(writer, sheet_name=sheet_name_for_stock_solutions, index=False)

    return stock_solution_containers


def sort_events_according_to_aspiration_volume(event_list_chem):

        output_list = []

        # # group event by plate id
        # get mask for change of plate id
        split_index = []
        for index in range(1, len(event_list_chem)):
            if event_list_chem[index].destination_container.id['plate_id'] != \
                event_list_chem[index-1].destination_container.id['plate_id']:
                split_index.append(index)
        # group event by mask
        events_grouped_by_plate_id = np.split(event_list_chem, split_index)

        # group and sort events in each plate according to substance
        for index, plate_event in enumerate(events_grouped_by_plate_id):
            # get mask for change of substance
            split_index = []
            for index in range(1, len(plate_event)):
                if plate_event[index].substance != \
                    plate_event[index-1].substance:
                    split_index.append(index)
            # group event by mask
            split_list = np.split(plate_event, split_index)

            for event_of_one_substance in split_list:
                # sort events in each group according to aspiration volume
                sorted_list = sorted(event_of_one_substance, key=lambda x: x.aspirationVolume, reverse=True)
                output_list.extend(sorted_list)

        return output_list


def turn_off_lld(event_list):
    for event in event_list:
        event.asp_lld = 0


def sort_events_by_substance_volume(event_list):
    event_list_sorted = []

    # split the event list by substance
    split_index = []
    for index in range(1, len(event_list)):
        if event_list[index].substance != \
                event_list[index-1].substance:
            split_index.append(index)
    split_list = np.split(event_list, split_index)

    for event_of_one_substance in split_list:
        # sort events in each group according to aspiration volume
        sorted_list = sorted(event_of_one_substance, key=lambda x: x.transfer_volume, reverse=True)
        event_list_sorted.extend(sorted_list)

    return event_list_sorted

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


if __name__ == '__main__':

    ## initiate hardware
    zm, gt, pt = initiate_hardware()

    excel_path_before_treatment, \
    plate_barcodes, \
    reaction_temperature,\
    plate_barcode_for_dilution = prep.GUI_get_excel_path_plate_barcodes_temperature_etc()

    logger.info(f"excel_path_before_treatment: {excel_path_before_treatment}\n" \
    f"plate_barcodes: {plate_barcodes}\n" \
    f"reaction_temperature: {reaction_temperature}\n" \
    f"plate_barcode_for_dilution: {plate_barcode_for_dilution}")

    excel_path_for_conditions, _ = prep.prepare_excel_file_for_reaction(reaction_temperature=reaction_temperature,
                                                                        excel_path=excel_path_before_treatment,
                                                                        plate_barcodes=plate_barcodes,
                                                                        plate_barcodes_for_dilution=plate_barcode_for_dilution)

    check_volume_by_pipetter = check_volume_pysimplegui()

    stock_solution_containers = \
        pln.assign_stock_solutions_to_containers_and_check_volume(excel_path = excel_path_for_conditions,
                                                                  check_volume_by_pipetter = check_volume_by_pipetter, pt=pt)

    df_reactions_grouped_by_plate_id,  substance_addition_sequence = prep.extract_reactions_df_to_run(excel_path_for_conditions)


    event_list_to_run = pln.generate_event_list_new(excel_path_for_conditions = excel_path_for_conditions,
                        df_reactions_grouped_by_plate_id = df_reactions_grouped_by_plate_id,
                        substance_addition_sequence = substance_addition_sequence,
                        stock_solution_containers = stock_solution_containers,
                        asp_lld = 0)

    event_list_to_run_sorted = sort_events_by_substance_volume(event_list_to_run)

    assert len(event_list_to_run_sorted) > 0, "event_list_to_run_sorted is empty"

    # for i in range(9):
    #     event_list_to_run[i].destination_container.xy = brb.plate_list[2].containers[i].xy

    # pln.run_events_chem(zm=zm, pt=pt,
    #                     event_list=event_list_to_run_sorted,
    #                     prewet_tip=False,
    #                     pause_after_every_plate_min=0,
    #                     test_mode=False,
    #                     pause_after_every_addition_sec=0,)

    # for i in range(5):
    #     event_list_to_run[i].destination_container.id['plate_id'] = 1
    #     event_list_to_run[i].destination_container.id['container_id'] = event_list_to_run[i].destination_container.id['container_id']+45
    # for index, event in enumerate(event_list_to_run):
    #     print(event.destination_container.id)

# c = brb.plate2.containers[-6].xy
# event_list_to_run_sorted[1].liquidClassTableIndex = 23
# event_list_to_run_sorted[0].liquidClassTableIndex = 24
# event_list_to_run_sorted[0].destination_container.xy = c
# event_list_to_run_sorted[1].destination_container.xy = c
# event_list_to_run_sorted[0].source_container.liquid_surface_height += 50
# event_list_to_run_sorted[1].source_container.liquid_surface_height += 50
#
# event_list_to_run_sorted[0].transfer_volume = 1000
# event_list_to_run_sorted[1].transfer_volume = 80
#
# event_list_to_run_sorted[0].destination_container.liquid_surface_height = 2000
# event_list_to_run_sorted[1].destination_container.liquid_surface_height = 2000
#
# event_list_to_run_sorted[0].destination_container.liquid_surface_height_real_value_float = 1900
# event_list_to_run_sorted[1].destination_container.liquid_surface_height_real_value_float = 1900


    pln.run_events_chem(zm=zm, pt=pt,
                        event_list=event_list_to_run_sorted,
                        prewet_tip=False,
                        pause_after_every_plate_min=0,
                        test_mode=False,
                        pause_after_every_addition_sec=0)