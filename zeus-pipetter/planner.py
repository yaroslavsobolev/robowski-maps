
import logging, os, pickle, json, numpy as np, winsound, copy, math, re, pandas as pd, time, csv
from math import ceil
from datetime import datetime
from typing import Dict, Any
import PySimpleGUI as sg

module_logger = logging.getLogger('main.planner')


import breadboard as brb, config


class EventInterpreter:
    '''This class is used to interpret the MS Excel file to a list of pipetting df.
    **How to use:
        1, upon calling __init__, the MS Excel file will be read into a
    pandas dataframe called self.reaction_df. Also, an empty dataframe called self.pd_output
    will be created.
        2, call add_events_to_df() to add events to self.pd_output.
    '''
    def __init__(self,
                 dataframe_filename: str,
                 is_for_bio=False
                 ):
        self.dataframe_filename = dataframe_filename

        self.logger = logging.getLogger('main.planner.EventInterpreter')
        self.logger.info(f'Interpretating dataframe_filename from {self.dataframe_filename}')

        self.reaction_df = pd.read_excel(io=self.dataframe_filename,
                                         sheet_name= 'reactions_with_run_info')

        # remove the first column from df
        # only choose the columns stating with 'vol#'1
        self.reaction_df = self.reaction_df.loc[:, [i for i in self.reaction_df.columns if 'vol#' in i]]

        self.pd_output = pd.DataFrame(columns=['reaction_id',
                                               'event_id',
                                               'plate_number',
                                               'plate_number_barcode',
                                               'plate_id_on_breadboard',
                                               'container_id',
                                               'substance',
                                               'transfer_volume',
                                               ], dtype=object)
        self.reaction_plates = pd.DataFrame(dtype=object)
        self.is_for_bio = is_for_bio
        if self.is_for_bio:
            self.containers_per_plate: int = 96
        else:
            self.containers_per_plate: int = 54

    def correct_order(self):  # TODO: use this function to correct the order of the adding substances
        pass

    def add_events_to_df(self):
        # print('plate_id')
        event_id = 0

        # generate a list of plates to accommodate all the reactions
        if self.reaction_df.shape[0] // self.containers_per_plate == 0:
            plate_list = [0]  # note: list(range(0)) = []
        else:
            plate_list = list(range(self.reaction_df.shape[0] // self.containers_per_plate))
        #
        for plate_id in plate_list:
            # print(plate_id)
            df_here = self.reaction_df.copy()
            dataframes_for_this_plate = df_here[plate_id * self.containers_per_plate: (plate_id + 1) * self.containers_per_plate]
            # print(this_plate)
            for enum, substance in enumerate(dataframes_for_this_plate.columns):
                # print(this_plate[0].columns)
                for reaction_id in dataframes_for_this_plate[substance].index:
                    transfer_volume = dataframes_for_this_plate[substance][reaction_id]
                    if transfer_volume < 0.5 or math.isnan(transfer_volume):
                        continue
                    plate_id_on_breadboard = plate_id % 3  # The breadboard has 3 plates for chemical reactions using
                    # 2_mL vials, so the plate_id_on_breadboard is the remainder of plate_id divided by 3
                    if self.is_for_bio:
                        plate_id_on_breadboard = 3  # for bio, the 96 wells are only on plate_id = 3
                    dict_here = {
                        # event_id is the index of pipetting events
                        'event_id': event_id,
                        # reaction_id is the index of reactions, this is also the index of the MS excel file
                        "reaction_id": reaction_id,
                        # plate_number is the index of the plate indicated by the barcode on the physical plate
                        "plate_number": plate_id,
                        # plate bar code
                        'plate_number_barcode': str(plate_id).zfill(2),
                        # this is the index of the plate on the breadboard
                        'plate_id_on_breadboard': plate_id_on_breadboard,
                        # container_id is the index of the container on the plate
                        'container_id': reaction_id % self.containers_per_plate,  # output: 0-53
                        # what substance is going to be pipetted in this event
                        'substance': substance[4:], ## change from 'vol#Dioxane' to 'Dioxane'
                        # how much volume is going to be pipetted in this event
                        'transfer_volume': transfer_volume
                    }
                    # self.volume_update(volume = transfer_volume, source_container = source_container, destination_container = destination_container)
                    # print(_dict)
                    self.pd_output = pd.concat([self.pd_output, pd.DataFrame(dict_here, index=[event_id], dtype=object)],
                                               ignore_index=True)
                    event_id += 1


class TransferEventConstructor:

    def __init__(self, event_dataframe, containers_for_stock, pipeting_to_balance: bool = False):


        self.substance_name: str = event_dataframe['substance']

        self.event_label: str = ' event_id:' + str(event_dataframe['event_id']) + '   ' + \
                                'substance:' + str(event_dataframe['substance']) + '   ' + \
                                'transfer_volume:' + str(event_dataframe['transfer_volume'])
                                # + " " + 'plate_number_barcode:' + str(event_dataframe['plate_number_barcode'])

        self.source_container: object = [container for container in containers_for_stock
                                         if container.substance == self.substance_name][0]

        if pipeting_to_balance:
            self.destination_container: object = brb.balance_cuvette
        else:
            self.destination_container: object = brb.plate_list[event_dataframe['plate_id_on_breadboard']].containers[
                event_dataframe['container_id']]

        # for aspiration
        self.aspirationVolume: int = event_dataframe['transfer_volume']

        self.tip_type: str = self.choose_tip_type(self.aspirationVolume)

        self.asp_containerGeometryTableIndex = self.source_container.containerGeometryTableIndex

        self.asp_deckGeometryTableIndex: int = self.get_deck_index(self.tip_type)

        self.asp_liquidClassTableIndex: int = \
            self.get_liquid_class_index(solvent=self.source_container.solvent,
                                        mode=self.source_container.solvent,
                                        tip_type=self.tip_type)

        # self.asp_liquidSurface: int = self.get_liquid_surface(self.source_container)
        self.asp_liquidSurface: int = self.source_container.liquid_surface_height
        self.asp_lldSearchPosition: int = self.asp_liquidSurface - 50

        self.dispensingVolume: int = self.aspirationVolume

        self.disp_containerGeometryTableIndex: int = self.destination_container.containerGeometryTableIndex
        self.disp_deckGeometryTableIndex: int = self.asp_deckGeometryTableIndex
        self.disp_liquidClassTableIndex: int = self.asp_liquidClassTableIndex

        self.disp_liquidSurface: int = self.get_liquid_surface(self.destination_container)
        self.disp_lldSearchPosition: int = self.disp_liquidSurface - 50

        # default values
        self.asp_qpm: int = 1
        self.asp_lld: int = 1
        self.disp_lld: int = 0
        self.asp_mixVolume: int = 0
        self.asp_mixFlowRate: int = 0
        self.asp_mixCycles: int = 0
        self.disp_mixVolume: int = 0
        self.disp_mixFlowRate: int = 0
        self.disp_mixCycles: int = 0
        self.searchBottomMode: int = 0

        # event conducting status
        self.is_event_conducted: bool = False
        self.event_finish_time: str = None  # time: (UTC time, local time)

    def get_source_container(self, substance_name: str, source_containers=None):

        """
        source_substance_containers exp:
        {'DMF': {'plate_id': 4, 'container_id': 2}, 'amine': {'plate_id': 7, 'container_id': 15}}
        """
        # print(substance_name)
        source_containers = brb.source_substance_containers  # global variable, mutable object (list),
        # so it should not be used directly as a default argument
        for container in source_containers:  # iterate through keys
            # print(this_substance)
            # print(substance_name)
            # print(list(this_substance.keys())[0])
            if substance_name == container.substance:
                return container
                # plate_id = this_substance[substance_name]['plate_id']
                # container_id = this_substance[substance_name]['container_id']
                # # print(f'substance is found in container: brb.plate_list[{plate_id}].containers[{container_id}]')
                # return brb.plate_list[plate_id].containers[container_id]

        print(f'Subtance:: {substance_name} is not found in the stock container!')

    def get_deck_index(self, tip_type: str):
        # print(f'tipe_type: {tip_type}')

        tip_index_dict = {'50ul': 3, '300ul': 0, '1000ul': 1}

        return tip_index_dict[tip_type]

    def get_liquid_class_index(self, solvent: str, mode: str, tip_type: str):
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
        solvent_para = {solvent, mode, tip_type}  # define a set of paras
        for liquid_class, index in liquid_class_dict.items():
            solvent_para_here = set(liquid_class.split('_'))
            # print(f'solvent_para_here: {solvent_para_here}')
            # print(f'solvent_para: {solvent_para}')
            if solvent_para.issubset(solvent_para_here):
                return index

    def choose_tip_type(self, transfer_volume: int):
        if transfer_volume <= 50:
            return "50ul"
        if transfer_volume > 50 and transfer_volume <= 300:
            return "300ul"
        if transfer_volume > 300 and transfer_volume < 2000:
            return "1000ul"

    def get_liquid_surface(self, container: object) -> int:

        # container.liquid_volume. This is in uL
        # container.area This is in mm^2
        # container.liquid_volume / container.area This is in mm
        if container.container_shape == 'cylindrical':
            liquid_height = ((container.liquid_volume / container.area) * 10)

        elif container.container_shape == 'conical_1500ul':
            if container.liquid_volume <= 500:
                liquid_height = 17 * 10
            else:
                liquid_height = 17 * 10 + ((container.liquid_volume - 500) / container.area) * 10
        else:
            print('Container shape is not defined!')

        return round(container.bottomPosition - liquid_height)

class Event:

    def __init__(self, event_dataframe: pd.DataFrame = None,
                 column_to_generate_event:str = None,
                 stock_solution_containers: list = None,
                 excel_path_for_conditions: str = None,
                 pipeting_to_balance: bool = False,
                 asp_lld : int = 0,):

        self.condition_uuid: str = event_dataframe['uuid']
        self.plate_barcode: str = event_dataframe['plate_barcode']
        self.excel_path_for_conditions = excel_path_for_conditions
        self.substance: str = column_to_generate_event[4:] # remove "vol#"
        self.event_label: str = self.condition_uuid + '_' + self.substance

        # print(f'self.substance: {self.substance}')
        # TODO: assert

        assert len([container for container in stock_solution_containers if container.substance == self.substance])>0,\
            'Check substance and stock solution name! No match found!'

        self.source_container: object = [container for container in stock_solution_containers
                                         if container.substance == self.substance][0]
        assert self.source_container.substance == self.substance, "substance name not match"

        if pipeting_to_balance:
            self.destination_container: object = brb.balance_cuvette
        else:
            breadboard_plate_id = event_dataframe['slot_id']
            container_id = event_dataframe['container_id']
            self.destination_container: object = brb.plate_list[breadboard_plate_id].containers[container_id]


        # for aspiration
        self.transfer_volume: int = event_dataframe[column_to_generate_event]
        # assert self.transfer_volume >=10 and self.transfer_volume <= 1000, f"transfer volume is invalid: {self.transfer_volume}"

        self.tip_type: str = self.choose_tip_type(self.transfer_volume)
        assert self.tip_type in ['50ul', '300ul', '1000ul'], f"tip type is invalid: {self.tip_type}"

        self.asp_containerGeometryTableIndex = self.source_container.containerGeometryTableIndex

        self.asp_deckGeometryTableIndex: int =  {'50ul': 3, '300ul': 0, '1000ul': 1}[self.tip_type]
        assert self.asp_deckGeometryTableIndex in [0, 1, 3], "tip type is invalid!"

        self.solvent: str = self.source_container.solvent

        self.liquidClassTableIndex: int = self.get_liquid_class_index(solvent=self.source_container.solvent,
                                                                      mode=self.source_container.mode,tip_type=self.tip_type)

        assert self.liquidClassTableIndex in range(100), f"liquid class is invalid: {self.liquidClassTableIndex}"

        self.disp_containerGeometryTableIndex: int = self.destination_container.containerGeometryTableIndex
        self.disp_deckGeometryTableIndex: int = self.asp_deckGeometryTableIndex

        # default values
        self.asp_qpm: int = 1
        self.asp_lld: int = asp_lld
        self.disp_lld: int = 0
        self.asp_mixVolume: int = 0
        self.asp_mixFlowRate: int = 0
        self.asp_mixCycles: int = 0
        self.disp_mixVolume: int = 0
        self.disp_mixFlowRate: int = 0
        self.disp_mixCycles: int = 0
        self.searchBottomMode: int = 0

        # event conducting status
        self.is_event_conducted: bool = False
        self.event_finish_time: str = None  # time: (UTC time, local time)

    def get_liquid_class_index(self, solvent: str, mode: str, tip_type: str):
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
            'hbrhac1v1_empty_300ul_clld': 31,
            'hbrhac1v1_empty_50ul_clld': 31,
            'hbrhac1v1_empty_1000ul_clld': 31,
            'DCE_empty_50ul_clld': 36,
            'DCE_empty_300ul_clld': 37,
            'DCE_empty_1000ul_clld': 38,
            'MeCN_empty_50ul_clld': 39,
            'MeCN_empty_300ul_clld': 40,
            'MeCN_empty_1000ul_clld': 41,
        }

        solvent_para = {solvent, mode, tip_type}  # define a set of paras

        for liquid_class, index in liquid_class_dict.items():
            solvent_para_here = set(liquid_class.split('_'))
            # print(f'solvent_para_here: {solvent_para_here}')
            # print(f'solvent_para: {solvent_para}')
            if solvent_para.issubset(solvent_para_here):
                return index

    def choose_tip_type(self, transfer_volume: int):
        if transfer_volume <= 50:
            return "50ul"
        elif transfer_volume > 50 and transfer_volume <= 300:
            return "300ul"
        elif transfer_volume > 300 and transfer_volume < 2000:
            return "1000ul"
        else:
            raise ValueError(f'Transfer volume is out of range! at {self.condition_uuid}, {self.substance}')

    def execute_event(self, purpose: str = 'reaction'):

        self.source_container.liquid_volume -= self.transfer_volume
        self.destination_container.liquid_volume += self.transfer_volume
        print(f'self.source_container.liquid_surface_height_real_value_float before: {self.source_container.liquid_surface_height_real_value_float}')
        self.source_container.liquid_surface_height_real_value_float += (self.transfer_volume / self.source_container.area)*10
        print(f'self.source_container.liquid_surface_height_real_value_float after: {self.source_container.liquid_surface_height_real_value_float}')
        print(f'self.tansfer_volume: {self.transfer_volume}')
        print(f'self.source_container.area: {self.source_container.area}')
        print(f'self.tarnsfer_volume / self.source_container.area: {(self.transfer_volume / self.source_container.area)*10}')

        self.source_container.liquid_surface_height = ceil(self.source_container.liquid_surface_height_real_value_float)

        self.destination_container.liquid_surface_height_real_value_float -= (self.transfer_volume / self.destination_container.area)*10
        self.destination_container.liquid_surface_height = ceil(self.destination_container.liquid_surface_height_real_value_float)

        self.is_event_conducted = True
        self.event_finish_time = round(time.time())

        if purpose == 'reaction':
            ## change the reaction/substance status in the excel file.
            df_reactions = pd.read_excel(self.excel_path_for_conditions, sheet_name=config.sheet_name_for_run_info)
            # change the status of the substance
            # print(f'condition_uuid: {self.condition_uuid}')
            # print('full_status: ', df_reactions[df_reactions['uuid']==self.condition_uuid]['full_status'].item())
            dict_substance = json.loads(df_reactions[df_reactions['uuid']==self.condition_uuid]['full_status'].item())
            dict_substance[self.substance] = ['completed', self.event_finish_time]
            dict_substance_after = json.dumps(dict_substance)
            df_reactions.loc[df_reactions['uuid'] == self.condition_uuid, 'full_status']\
                = dict_substance_after

            # change the status of the reaction if necessary
            substances_for_this_condition = list(dict_substance.keys())  ## after python 3.6, dict is ordered.
            if self.substance in substances_for_this_condition[:-1]:
                df_reactions.loc[df_reactions['uuid'] == self.condition_uuid, 'status']\
                    = 'not_completed'
            elif self.substance == substances_for_this_condition[-1]:
                df_reactions.loc[df_reactions['uuid'] == self.condition_uuid, 'status']\
                    = 'completed'
                df_reactions.loc[df_reactions['uuid'] == self.condition_uuid, 'timestamp']\
                    = self.event_finish_time
            # save the dataframe to excel file
            with  pd.ExcelWriter(self.excel_path_for_conditions, engine='openpyxl', mode='a', if_sheet_exists="replace") as writer:
                df_reactions.to_excel(writer, sheet_name=config.sheet_name_for_run_info, index=False)

def interprete_events_from_excel_to_dataframe(dataframe_filename: str,
                                              is_for_bio: bool) -> pd.DataFrame:
    # generate empty dataframes
    event_dataframes = EventInterpreter(dataframe_filename=dataframe_filename,
                                        is_for_bio=is_for_bio)

    # add all events to the dataframe
    event_dataframes.add_events_to_df()

    module_logger.info(f'event_dataframe is generated with {len(event_dataframes.pd_output.index)} events.')

    event_dataframes.pd_output.to_json(
        f'C:\\Yankai\\event_dataframes\\event_dataframe_{datetime.now().strftime("%Y_%m_%d_%H_%M")}.json',
        orient='records', lines=True)
    module_logger.info(f'event_dataframe is saved as event_dataframe_{datetime.now().strftime("%Y_%m_%d_%H_%M")}.json')

    return event_dataframes.pd_output


def generate_event_object(logger: object, excel_to_generate_dataframe: str,containers_for_stock,
                          is_pipeting_to_balance: bool = False,
                          is_for_bio: bool = False) -> tuple:

    # generate event dataframes from excel
    event_dataframe = \
        interprete_events_from_excel_to_dataframe(dataframe_filename=excel_to_generate_dataframe,
                                                  is_for_bio=is_for_bio)

    logger.info(f"All events are generated to dataframes from excel here: {excel_to_generate_dataframe}")

    # generate event list
    event_list = generate_event_list_new(event_dataframe=event_dataframe,
                                     pipeting_to_balance=is_pipeting_to_balance,
                                     containers_for_stock=containers_for_stock)

    logger.info("All event objects are generated from dataframes.")

    return event_dataframe, event_list


def generate_event_list_new(df_reactions_grouped_by_plate_id,
                        substance_addition_sequence,
                        stock_solution_containers,
                        excel_path_for_conditions: str, asp_lld):
    pipetting_to_balance = False
    ## create a list of events
    event_list = []
    for df_reactions in df_reactions_grouped_by_plate_id:
        for substance in substance_addition_sequence:
            ## substance: vol#Dioxane
            for index, df_row in df_reactions.iterrows():
                ## pass the row to the Event class only if the substance volume is not 0
                ## AND IMPORTANTLY: the full_status of this substance should be "not_started"
                # full_status_dict = json.loads(df_reactions[df_reactions['uuid'] == df_reactions.condition_uuid]['full_status'].item())
                full_status_dict = json.loads(df_row['full_status'])
                # print(f'substance: {substance}, full_status_dict: {full_status_dict}')
                # print(f'substance[4:]: {substance}')
                # print(f'full_status_dict[substance[4:]]: {full_status_dict[substance[4:]]}')
                # print(f'full_status_dict[substance[4:]][0]: {full_status_dict[substance[4:]][0]}')
                # print(f'event.condition_uuid: {df_row["uuid"]}')
                # if df_row[substance] != 0 and full_status_dict[substance[4:]][0] == 'not_started':
                # print(f'substance: {substance}, substance to be added: {full_status_dict.keys()}')
                if  substance[4:] in full_status_dict.keys():
                    # print(f'substance: {substance}, full_status_dict: {full_status_dict[substance[4:]][0]}, event.condition_uuid: {df_row["uuid"]}')
                    if full_status_dict[substance[4:]][0] == 'not_started':
                        event = Event(event_dataframe=df_row,
                                      column_to_generate_event=substance,
                                      pipeting_to_balance=pipetting_to_balance,
                                      stock_solution_containers=stock_solution_containers,
                                      excel_path_for_conditions=excel_path_for_conditions,
                                      asp_lld=asp_lld)
                        # print(f'event: {event.event_label}. volume: {event.transfer_volume}')
                        event_list.append(event)
                    elif full_status_dict[substance[4:]][0] == 'completed':
                        module_logger.info(f'{substance} for {df_row["uuid"]} was done in previous run. skip this event.')

                    else:
                        raise ValueError(f'{substance} for {event.condition_uuid} is not "not_started" or "completed".\n'
                                         f'check the EXCEL file.')

    return event_list


def do_calibration_on_events(zm: object, pt: object, logger: object,
                             calibration_event_list: list[object],
                             change_tip_after_every_pipetting: bool,
                             repeat_n_times: int, starting_index: int,
                             prewet_tip:bool, use_calibration):
    '''This function is used to calibrate the pipetting of substances.'''
    results_for_calibration = []

    assert len(calibration_event_list) > 0, "calibration_event_list is empty!"

    if zm.tip_on_zeus:
        pt.discard_tip()

    starting_index = starting_index
    ending_index = len(calibration_event_list)

    for event_index in range(starting_index, ending_index):

        if zm.tip_on_zeus != calibration_event_list[event_index].tip_type:
            pt.change_tip(calibration_event_list[event_index].tip_type)
            logger.info(f'The tip is changed to : {calibration_event_list[event_index].tip_type}')
            time.sleep(0.5)

        # if prewet_tip:
        #     prewet_new_tip(zm=zm, pt=pt, pipetting_event=calibration_event_list[event_index], use_calibration  = use_calibration)

        result = pt.pipetting_to_balance_and_weight_n_times(zm= zm, pt= pt, transfer_event=calibration_event_list[event_index],
                                                            n_times=repeat_n_times,
                                                            change_tip_after_every_pipetting=change_tip_after_every_pipetting,
                                                            prewet_tip=prewet_tip, use_calibration=use_calibration)
        pt.discard_tip()

        results_for_calibration.append(result)


        # pt.home_xy()

        time.sleep(1)
        logger.info(f"Performed one measurement: {calibration_event_list[event_index].event_label}")
        logger.info(f'Result: {result}')

        # if change_tip_after_every_pipetting and event_index != 0:
        #     pt.discard_tip()
        #     time.sleep(0.5)

        time.sleep(0.5)

        result_dict = {f'{calibration_event_list[event_index].substance}_{calibration_event_list[event_index].transfer_volume}': results_for_calibration}
        logger.info(f'Result_here: result_dict')

    # save result_dict to jason file
    folder_for_result = '//'.join((calibration_event_list[0].excel_path_for_conditions.split('/')[:-1]))
    with open(folder_for_result+f'//calibration_results_{datetime.now().strftime("%Y_%m_%d")}_{calibration_event_list[0].substance}.json', 'a') as f:
            json.dump(result_dict, f, indent=4)

    # pt.discard_tip()

    return results_for_calibration


def prewet_new_tip( zm: object, pt: object , pipetting_event: object, use_calibration: bool = False):

    print('Prewetting new tip...')

    event_adjusted = copy.deepcopy(pipetting_event)

    max_volume = int(re.findall(r'\d+', zm.tip_on_zeus)[0])
    event_adjusted.transfer_volume = max_volume
    event_adjusted.destination_container = event_adjusted.source_container
    print(f'lqiuid_height_prewetting: {event_adjusted.source_container.liquid_surface_height}')
    event_adjusted.disp_liquidSurface = 1800

    module_logger.info(f'Prewetting tip with {max_volume}ul of {event_adjusted.substance}')
    pt.transfer_liquid(event_adjusted, use_calibration=use_calibration)
    module_logger.info('Prewet done! Continue with pipetting...')


def beep():
    duration = 600  # milliseconds
    freq = 1000  # Hz
    # time.sleep(0.2)
    winsound.Beep(freq, duration)

def beep_n():
    duration = 600  # milliseconds
    freq = 1000  # Hz
    # time.sleep(0.2)
    for i in range(5):
        winsound.Beep(freq, duration)


def assign_stock_solutions_to_containers_and_check_volume(excel_path:str, pt, check_volume_by_pipetter: bool = True):
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
        mode\
        = row['substance'], row['breadboard_plate_id'], row['container_id'], row['solvent'], row['density'], row['volume'], row['liquid_surface_height'], row['pipetting_mode']

        container_for_this_stock_solution = brb.plate_list[plate_id].containers[container_id]
        container_for_this_stock_solution.substance = substance_name
        container_for_this_stock_solution.substance_density = density
        container_for_this_stock_solution.solvent = solvent
        container_for_this_stock_solution.pipetting_mode = mode

        if not check_volume_by_pipetter:
            container_for_this_stock_solution.liquid_surface_height = liquid_surface_height+40 # 2mm deeper under the liquid surface
            container_for_this_stock_solution.liquid_surface_height_real_value_float = float(liquid_surface_height)
            container_for_this_stock_solution.liquid_volume = volume_ml

        else:

            height, volume = pt.check_volume_in_container(container=container_for_this_stock_solution,
                                             liquidClassTableIndex=13,change_tip_after_each_check=True)
            container_for_this_stock_solution.liquid_surface_height = height+20# 2mm deeper under the liquid surface
            container_for_this_stock_solution.liquid_volume = volume
            print(f'the liquid surface height detected: {height}')
            print(f'the volume detected: {volume}')
            container_for_this_stock_solution.liquid_surface_height_real_value_float=float(height)
            print(f'the liquid surface height in float: {float(height)}')

        stock_solution_containers.append(container_for_this_stock_solution)

        ## update the excel after checking the volume
        if check_volume_by_pipetter:
            df_stock_solutions.loc[index, 'volume'] = container_for_this_stock_solution.liquid_volume
            df_stock_solutions.loc[index, 'liquid_surface_height'] = container_for_this_stock_solution.liquid_surface_height
            with pd.ExcelWriter(excel_path, engine='openpyxl', mode='a', if_sheet_exists="replace") as writer:
                df_stock_solutions.to_excel(writer, sheet_name=sheet_name_for_stock_solutions, index=False)

    return stock_solution_containers


def assign_stock_solutions_to_containers_and_check_volume_new(excel_path:str, pt: object,
                                                          sheet_name: str = None,
                                                          check_volume_by_pipetter: bool = True):
    if sheet_name is None:
        sheet_name_for_stock_solutions = 'stock_solutions'
    else:
        sheet_name_for_stock_solutions = sheet_name

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

            height, volume = pt.check_volume_in_container(container=container_for_this_stock_solution,
                                             liquidClassTableIndex=13,change_tip_after_each_check=True)

            container_for_this_stock_solution.liquid_surface_height = height
            container_for_this_stock_solution.liquid_surface_height_real_value_float = float(height)
            container_for_this_stock_solution.liquid_volume = volume

        stock_solution_containers.append(container_for_this_stock_solution)

        ## update the excel after checking the volume
        if check_volume_by_pipetter:
            df_stock_solutions.loc[index, 'volume'] = container_for_this_stock_solution.liquid_volume
            df_stock_solutions.loc[index, 'liquid_surface_height'] = container_for_this_stock_solution.liquid_surface_height
            with pd.ExcelWriter(excel_path, engine='openpyxl', mode='a', if_sheet_exists="replace") as writer:
                df_stock_solutions.to_excel(writer, sheet_name=sheet_name_for_stock_solutions, index=False)

    return stock_solution_containers

def run_one_event_chem(pt: object, purpose='mixing', event=None):
    pt.transfer_liquid(event)
    time.sleep(0.1)
    event.execute_event()
    beep()

def run_events_chem(zm: object, pt: object,
                    event_list=None,
                    change_tip_after_every_pipetting: bool = False,
                    prewet_tip: bool = True,
                    pause_after_every_plate_min: int = 0,
                    test_mode: bool = False,
                    pause_after_every_addition_sec = 0):

    # discard tip if there is one on the pipetter
    if zm.tip_on_zeus != '': # tip_on_zeus: '' or '50ul' or '300ul' or '1000ul'
        pt.discard_tip()

    # group all events by plate_id
    split_index = []
    for index in range(1, len(event_list)):
        if event_list[index].plate_barcode != event_list[index - 1].plate_barcode:
            split_index.append(index)

    events_grouped_by_plate_id = np.split(event_list, split_index)

    # run events in each plate
    for plate_index, events_in_one_plate in enumerate(events_grouped_by_plate_id):

        pt.home_xy()
        start_time_datetime = datetime.now()
        module_logger.info(f'Gantry is homed at {start_time_datetime}, before a new plate is processed.')

        for index, event in enumerate(events_in_one_plate):
            # change tip if needed
            if zm.tip_on_zeus != event.tip_type:
                pt.change_tip(event.tip_type) if not test_mode else print("this is test mode.")
                module_logger.info(f'The tip is changed to : {event.tip_type}')
                time.sleep(0.5)

                # do prewet every time a new tip is taken
                if prewet_tip:
                    prewet_new_tip(zm=zm, pt=pt, pipetting_event=event)

            try:
                if not test_mode:
                    run_one_event_chem(pt= pt, event= event)
                    time.sleep(pause_after_every_addition_sec)
                else:
                    print(f"Test mode one done: {event.event_label}")

            except Exception as e:
                module_logger.error(f'Error in transfer liquid.\n')
                module_logger.error(f"The run is stopped at event {event}.\n")
                pt.discard_tip()
                raise e

            module_logger.info(f"Event has been performed: {event.event_label}, "
                        f"{event.transfer_volume}uL,"
                        f'slot: {event.destination_container.id["plate_id"]},'
                        f'plate_code: {event.plate_barcode}')

            ## change tip when the substance changed
            if index < len(events_in_one_plate)-1:
                if events_in_one_plate[index].substance != events_in_one_plate[index + 1].substance:
                        pt.discard_tip()

            # time.sleep(0.5)

            #  change tip after every pipetting if specified
            if change_tip_after_every_pipetting:
                pt.discard_tip()
                time.sleep(0.5)

            # check if this is the last event in the plate
            if index == len(events_in_one_plate)-1:

                ## play some sound to notify the user
                beep_n()
                time.sleep(2)
                beep_n()
                time.sleep(2)

                module_logger.info(f'This is the last event in this plate: {event.event_label}. plate_barcode: {event.plate_barcode}.')

                ## check if this is the last plate
                if plate_index == len(events_grouped_by_plate_id)-1:
                    ## play some sound to notify the user
                    beep_n()
                    time.sleep(2)
                    ## prommp a pysimplegui window for completion
                    sg.popup_ok('All events are done. Please remove the plate and click OK to finish.')
                else:
                    continue

        pt.discard_tip()

        ## pause after each plate
        if pause_after_every_plate_min > 0:
            ## make a pysimplegui window to show the pause time
            time.sleep(pause_after_every_plate_min * 60)

        #####the following setction is commented out because the malfuction of the pipettor pressure sensor.#### 20240226Yankai
        # # check if this plate is the last plate and re-measure the liquid surface height if not the last plate
        # if plate_index != len(events_grouped_by_plate_id)-1:
        #     # measure the liquid surface height again and update the value to the stock container
        #     assign_stock_solutions_to_containers_and_check_volume(excel_path=event.excel_path_for_conditions,
        #                                                          check_volume_by_pipetter=True, pt=pt)


def run_events_chem_nps(zm: object, pt: object, logger: object, start_event_id: int,
                    event_list_path=None, event_list=None,
                    change_tip_after_every_pipetting: bool = False,
                    prewet_tip: bool = True) -> dict[Any, Any]:

    # for event list either specify a path or a list. Only speficify one of them.
    if event_list_path is not None:
        with open(event_list_path, 'rb') as f:
            event_list = pickle.load(f)

    liquid_surface_height_from_zeus = {}

    if zm.tip_on_zeus:
        pt.discard_tip()
    #
    # for event in event_list:
    #     event.asp_lld = 1

    for event_index in range(start_event_id, len(event_list)):

        if zm.tip_on_zeus != event_list[event_index].tip_type:
            pt.change_tip(event_list[event_index].tip_type)
            logger.info(f'The tip is changed to : {event_list[event_index].tip_type}')
            # do prewet every time a new tip is taken
            time.sleep(0.5)
            if prewet_tip:
                prewet_new_tip(zm=zm, pt=pt, pipetting_event=event_list[event_index])

        # record start time
        event_start_time = int(time.time())  # unix time
        event_start_time_datetime = datetime.fromtimestamp(event_start_time)
        event_list[event_index].event_start_time_utc = event_start_time
        event_list[event_index].event_start_time_datetime = str(event_start_time_datetime)
        try:
            liquid_surface_height_from_zeus_here = pt.transfer_liquid(event_list[event_index])
            beep()

        except Exception as e:
            logger.error(f'Error in transfer liquid.\n '
                         f'\t\t\t\t\t\tConsider adding more liquid to source container: '
                         f'{event_list[event_index].source_container.container_id}\n'
                         f'\t\t\t\t\t\tNext, proceed with: event_number{event_index+1}, {event_list[event_index].event_label}')

            with open(f'multicomponent_reaction\\event_list_chem_id_finished_at_{event_index-1}.pickle', 'wb') as f:
                pickle.dump(event_list, f)

            pt.discard_tip()

            return liquid_surface_height_from_zeus

        # record finish time
        event_finish_time = int(time.time())  # UTC time
        event_finish_time_datetime = datetime.fromtimestamp(event_finish_time)
        event_list[event_index].event_finish_time = event_finish_time
        event_list[event_index].event_finish_time_datetime = str(event_finish_time_datetime)
        event_list[event_index].is_event_conducted = True

        # calculate volume and liquid height
        liquid_surface_height_from_zeus[event_list[event_index].substance + '_height'] = \
            liquid_surface_height_from_zeus_here
        volume_here = (-liquid_surface_height_from_zeus_here + event_list[event_index].source_container.bottomPosition) \
                      * event_list[event_index].source_container.area / 10  # in uL
        liquid_surface_height_from_zeus[event_list[event_index].substance + '_volume'] = round(volume_here, 1)

        time.sleep(0.05)
        logger.info(f"Performed one event: event_number {event_index}, "
                    f"{event_list[event_index].event_label}")

        # check tip type and change tip if needed
        if event_index != len(event_list) - 1:  # check if this is the last event.
            if event_list[event_index].substance != event_list[event_index + 1].substance:
                # pt.change_tip(event_list[event_index + 1].tip_type)
                pt.discard_tip()

        time.sleep(0.5)

        # with open(f'multicomponent_reaction\\event_list_chem_{datetime.now().strftime("%Y_%m_%d_%H_%M")}.pickle', 'wb') as f:
        #     pickle.dump(event_list, f)

        if change_tip_after_every_pipetting:
            pt.discard_tip()
            time.sleep(0.5)

    pt.discard_tip()

    return liquid_surface_height_from_zeus


def run_events_chem_dilution(zm: object, pt: object, logger: object,
                             start_event_id: int,
                             skip_vial_id: tuple = (),
                             event_list=None,
                             change_tip_after_every_pipetting: bool = False,
                             log_to_excel: bool=False, excel_path: str=None,
                             barcode_of_plate_for_reactions: int = None,purpose = 'diluting',
                             end_event_id = None, prewet_tip=True):



    if zm.tip_on_zeus:
        pt.discard_tip()

    end_event_id = len(event_list) if end_event_id is None else end_event_id  # assign the end_event_id if it is not specified

    for event_index in range(start_event_id, end_event_id):

        print(f'event_index: {event_index}')
        print(f'event: {event_list[event_index]}')

        if event_index in skip_vial_id:
            print(f"Skip event {event_index}.\n")
            logger.info(f"Skip event {event_index}.\n")
            continue

        if zm.tip_on_zeus != event_list[event_index].tip_type:
            pt.change_tip(event_list[event_index].tip_type)
            logger.info(f'The tip is changed to : {event_list[event_index].tip_type}')

            # do prewet every time a new tip is taken
            if prewet_tip:
                prewet_new_tip(zm=zm, pt=pt, pipetting_event=event_list[event_index])

        try:
            pt.transfer_liquid(event_list[event_index], purpose=purpose)
            event_list[event_index].source_container.liquid_surface_height += \
                ceil((event_list[event_index].transfer_volume / event_list[event_index].source_container.area) * 10)
            beep()

        except Exception as e:
            print(f"Error in transfer liquid.\n")

            # with open(f'multicomponent_reaction\\event_list_chem_id_{event_index}.pickle', 'wb') as f:
            #     pickle.dump(event_list, f)

            pt.discard_tip() # discard the tip to trash bin if there is an error
            raise e # raise the error to stop the program


        logger.info(f"Performed one pipetting for dilution: event_number {event_index}, "
                    f"volume: {event_list[event_index].transfer_volume},"
                    f"from {event_list[event_index].source_container.id},"
                    f"to {event_list[event_index].destination_container.id},")

        if log_to_excel:
            df = pd.read_excel(excel_path, sheet_name=config.sheet_name_for_run_info, engine='openpyxl')
            ## add a column called "timestamp_dilution" to the dataframe
            if "timestamp_dilution" not in df.columns:
                df["timestamp_dilution"] = "0"
            ## write the timestamp to the dataframe
            ## locat the cell by plate barcode and the container id
            print(f'barcode_of_plate_for_reactions: {barcode_of_plate_for_reactions}')
            print(f'event_list[event_index].source_container.id[container_id]: {event_list[event_index].source_container.id["container_id"]}')

            df.loc[(df['plate_barcode'] == barcode_of_plate_for_reactions)&(df['container_id']== event_list[event_index].destination_container.id['container_id']), 'timestamp_dilution'] \
                = int(time.time())

            ## write the dataframe to the excel file
            with pd.ExcelWriter(excel_path, engine='openpyxl', mode='a', if_sheet_exists="replace") as writer:
                df.to_excel(writer, sheet_name=config.sheet_name_for_run_info, index= False)


        # check tip type and change tip if needed
        if event_index != len(event_list) - 1:  # check if this is the last event.
            if event_list[event_index].substance != event_list[event_index + 1].substance:
                pt.change_tip(event_list[event_index + 1].tip_type)

        # time.sleep(0.1)

        # with open(f'multicomponent_reaction\\event_list_chem_dilution_\\pickle_output\\'
        #           f'{datetime.now().strftime("%Y_%m_%d_%H_%M")}.pickle', 'wb') as f:
        #     pickle.dump(event_list, f)

        if change_tip_after_every_pipetting:
            pt.discard_tip()
            # time.sleep(0.1)

    if zm.tip_on_zeus:
        pt.discard_tip()


def run_events_chem_dilution_new(zm: object, pt: object, logger: object,
                             start_event_id: int,
                             end_event_id: int,
                             dilution_mode: str,
                             skip_vial_id: tuple = (),
                             event_list=None,
                             step_id:int = 0,
                             change_tip_after_every_pipetting: bool = False,
                             log_to_excel: bool=False,
                             excel_path: str=None,
                             barcode_of_plate_for_reactions: int = None,
                             purpose = 'diluting',
                             prewet_tip=True):

    if event_list == [] or event_list is None:
        return 0

    if zm.tip_on_zeus: # discard previous tip if there is one on the pipetter
        pt.discard_tip()

    for event_index in range(start_event_id, end_event_id):

        if event_index in skip_vial_id: # skip the vial if it is in the skip_vial_id list
            print(f"Skip event {event_index}.\n")
            logger.info(f"Skip event {event_index}.\n")
            continue

        if zm.tip_on_zeus != event_list[event_index].tip_type: # change tip if needed
            pt.change_tip(event_list[event_index].tip_type)
            logger.info(f'The tip is changed to : {event_list[event_index].tip_type}')

            # do prewet every time a new tip is taken
            if prewet_tip:
                prewet_new_tip(zm=zm, pt=pt, pipetting_event=event_list[event_index])

        try:
            pt.transfer_liquid(event_list[event_index])
            event_list[event_index].execute_event(purpose= 'dilution')
            beep()

        except Exception as e:
            print(f"Error in transfer liquid.\n")
            pt.discard_tip() # discard the tip to trash bin if there is an error
            raise e # raise the error to stop the program


        logger.info(f"Performed one pipetting for dilution: event_number {event_index}, "
                    f"volume: {event_list[event_index].transfer_volume},"
                    f"from {event_list[event_index].source_container.id},"
                    f"to {event_list[event_index].destination_container.id},")

        if log_to_excel and step_id == 2:
            df = pd.read_excel(excel_path, sheet_name=config.sheet_name_for_run_info, engine='openpyxl')
            if dilution_mode == "A->B":
                ## write the timestamp to the dataframe
                df.loc[df['uuid'] == event_list[event_index].condition_uuid, 'timestamp_dilution'] \
                    = int(time.time())
            elif dilution_mode == "A->C" or dilution_mode == "B->C":
                ## write the timestamp to the dataframe
                df.loc[df['uuid'] == event_list[event_index].condition_uuid, 'timestamp_dilution_2'] \
                    = int(time.time())

            ## write the dataframe to the excel file
            with pd.ExcelWriter(excel_path, engine='openpyxl', mode='a', if_sheet_exists="replace") as writer:
                df.to_excel(writer, sheet_name=config.sheet_name_for_run_info, index=False, index_label='uuid')


        # check tip type and change tip if needed
        if event_index != len(event_list) - 1:  # check if this is the last event.
            if event_list[event_index].substance != event_list[event_index + 1].substance:
                pt.change_tip(event_list[event_index + 1].tip_type)

        # change tip after every 5 pipetting, also change if step_id == 1 (This means the transfer is from A to B, or A to C, or B to C)
        if change_tip_after_every_pipetting or step_id ==1:
            pt.discard_tip()
            # time.sleep(0.1)

    if zm.tip_on_zeus:
        pt.discard_tip()


if __name__ == "__main__":
    print('You are runing the module: ', __name__)
