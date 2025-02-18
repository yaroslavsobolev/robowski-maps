
import logging


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
    fh = logging.FileHandler('C:\\Users\\Chemiluminescence\\Dropbox\\robochem\\pipetter_files\\main.log')
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


import copy, time, pickle, re, importlib, json, os
from datetime import datetime
from typing import List
from openpyxl import Workbook
from openpyxl import load_workbook
import PySimpleGUI as sg
import pandas as pd

import zeus
import pipetter
import planner as pln
import breadboard as brb


## TODO 2023-03-21:
# 1. module_logger is not working;
# 2. dispensing height is not adjustable, no idea what is wrong.
# 3. coordinates need to be further optimized.


data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

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

## initiate hardware
zm, gt, pt = initiate_hardware()
time.sleep(1)

## note: sometimes the arduino port fails to initiate, restart the computer to fix it.
## better way need to be found.

##  generate_dilution_events()
# step1: dilution original reactions, adding volume: 1400ul
# step2: transfer liquid from original reaction to new vial, transfer volume: 20ul
# step3: dilution new vial, adding volume: 480ul

# volume setting versions:
# 1. 200ul, 1600ul, 22.5ul, 477.5ul
# 2. 200ul, 1400ul, 20ul, 480ul
# 3. 200ul, 1000ul, 15ul, 485ul ## this one is now used, 2023-03-22 14:23

# specify volumes for dilution
volume_added_to_old_vial = 1000
volume_transfered_from_old_to_new_vial = 15
volume_added_to_new_vial = 485

event_list_dilute_old_vial = []
event_list_dilution_old_to_new = []
event_list_dilute_new_vial = []

event_list_transfer_to_54_vials = []

def generate_dilution_event(source_container: object = None,
                            destination_container: object = None,
                            volume: float = 0,
                            asp_liquid_surface: int = 0,
                            disp_liquid_surface: int = 0):
    # load template event
    with open('multicomponent_reaction\\template\\dilution_template.pickle', 'rb') as f:
        event_template = pickle.load(f)
    # always creat a copy of the template event
    event = copy.deepcopy(event_template)
    # asign new parameters to the event
    event.source_container = source_container
    event.destination_container = destination_container
    event.aspirationVolume = volume
    event.dispensingVolume = volume
    event.event_label = f'transfer from {source_container.container_id} to {destination_container.container_id}'

    print(f'volume: {volume}')
    if volume <= 50:
        event.tip_type = '50ul'
        event.asp_liquidClassTableIndex = 24
        event.disp_liquidClassTableIndex = 24

    elif volume <= 300 and volume > 50:
        event.tip_type = '300ul'
        event.asp_liquidClassTableIndex = 22
        event.disp_liquidClassTableIndex = 22

    else:
        event.tip_type = '1000ul'
        event.asp_liquidClassTableIndex = 23
        event.disp_liquidClassTableIndex = 23

    event.asp_liquidSurface = asp_liquid_surface
    event.asp_lld = 1
    event.asp_lldSearchPosition = asp_liquid_surface - 50

    event.disp_liquidSurface = disp_liquid_surface
    event.disp_lld = 0

    event.asp_containerGeometryTableIndex = source_container.containerGeometryTableIndex
    event.disp_containerGeometryTableIndex = destination_container.containerGeometryTableIndex

    return event


# step1: dilution original reactions, adding volume: 1400ul
def dilute_old_vial(skip_vials=(), rows_to_dilute=(0, 18, 36)): # diluting volume 1400ul
    global event_list_dilute_old_vial
    event_list_dilute_old_vial = []

    # generate dilution events
    for i in rows_to_dilute:
        for vial_index in range(i, i+9):
            if vial_index in skip_vials:
                print(f'skipping vial with index {vial_index}')
                continue
            source_container = copy.deepcopy(brb.plate_list[6].containers[0])
            destination_container = copy.deepcopy(brb.plate_list[2].containers[vial_index])
            event_temp = generate_dilution_event(source_container=source_container,
                                                destination_container=destination_container,
                                                volume=volume_added_to_old_vial,
                                                asp_liquid_surface = 1600,
                                                disp_liquid_surface = 2100)
            event_list_dilute_old_vial.append(event_temp)
    # time.sleep(2)

    ## run dilution events
    pln.run_events_chem_dilution(zm=zm, pt=pt, logger=logger,
                        event_list= event_list_dilute_old_vial, start_event_id=0)


## this is for transfering liquid from jar to 54 containers,
# for testing spectrophotometer repeatability, 2021-03-22 14:23
def transfer_to_54_vials(volume_added_to_vial = 500): # diluting volume 1400ul
    global event_list_transfer_to_54_vials
    event_list_transfer_to_54_vials = []

    # generate dilution events
    for vial_index in range(54):
            source_container = copy.deepcopy(brb.plate_list[6].containers[0])
            destination_container = copy.deepcopy(brb.plate_list[2].containers[vial_index])
            event_temp = generate_dilution_event(source_container=source_container,
                                                destination_container=destination_container,
                                                volume=volume_added_to_vial,
                                                asp_liquid_surface = 1600,
                                                disp_liquid_surface = 2000)
            event_list_transfer_to_54_vials.append(event_temp)
    # time.sleep(2)

    ## run dilution events
    pln.run_events_chem_dilution(zm=zm, pt=pt, logger=logger,
                        event_list= event_list_transfer_to_54_vials, start_event_id=0)


# step2: transfer liquid from original reaction to new vial, transfer volume: 15ul
def transfer_liquid_from_old_vial_to_new( skip_vials=(), rows_to_dilute=(0, 18, 36)): # transfer volume 20ul
    global event_list_dilution_old_to_new
    event_list_dilution_old_to_new = []

    for i in rows_to_dilute:
        for vial_index in range(i, i + 9):
            if vial_index in skip_vials:
                print(f'skipping vial with index {vial_index}')
                continue
            source_container = copy.deepcopy(brb.plate_list[2].containers[vial_index])
            destination_container = copy.deepcopy(brb.plate_list[2].containers[vial_index+9])
            event_temp = generate_dilution_event(source_container=source_container,
                                                 destination_container=destination_container,
                                                 volume=volume_transfered_from_old_to_new_vial,
                                                 asp_liquid_surface=1900,
                                                 disp_liquid_surface=2100)
            event_list_dilution_old_to_new.append(event_temp)

    # time.sleep(1)
    pln.run_events_chem_dilution(zm=zm, pt=pt, logger=logger,
                        event_list=event_list_dilution_old_to_new, start_event_id=0,
                        change_tip_after_every_pipetting = True)


# step3: dilution new vial, adding volume: 485ul

def dilute_new_vial(skip_vials=(), rows_to_dilute=(0, 18, 36)): # diluting volume 485ul
    global event_list_dilute_new_vial
    # TODO: These two added nines in two difference places of these loops are confusing. Logic here should be more transparent.
    for i in [x + 9 for x in rows_to_dilute]:
        for vial_index in range(i, i + 9):
            if vial_index in skip_vials:
                print(f'skipping vial with index {vial_index}')
                continue
            source_container = copy.deepcopy(brb.plate_list[6].containers[0])
            destination_container = copy.deepcopy(brb.plate_list[2].containers[vial_index])
            event_temp = generate_dilution_event(source_container=source_container,
                                                 destination_container=destination_container,
                                                 volume=volume_added_to_new_vial,
                                                 asp_liquid_surface=1600,
                                                 disp_liquid_surface=2100)
            event_list_dilute_new_vial.append(event_temp)

    # time.sleep(2)
    pln.run_events_chem_dilution(zm=zm, pt=pt, logger=logger,
                        event_list=event_list_dilute_new_vial, start_event_id=0)

if __name__ == '__main__':
    mins_to_wait = 0
    print(f'Waiting for {mins_to_wait} minutes')
    time.sleep(60*mins_to_wait)

    print('Starting dilution')
    # dilute_old_vial(skip_vials = ())
    transfer_liquid_from_old_vial_to_new(skip_vials = ())
    dilute_new_vial(skip_vials = ())

# transfer one product to 54 vails
#     transfer_to_54_vials()
