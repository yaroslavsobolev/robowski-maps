"""
breadboard.py

This module defines the physical and logical structure of the breadboard used in
the automated liquid handling platform. It provides classes and utilities to configure
and manage lab containers, plates, decks, and pipette tips for robotic liquid transfers.

Core Features:
--------------
1. **Container Definitions**:
   Predefined container types such as:
   - vial_2ml
   - well_bio (microplate wells)
   - bottle_20ml
   - jar_100ml
   - tube_1.5ml (tube_1500ul)
   - balance_cuvette
   - nanodrop_pedestal

2. **Plate Setup**:
   Eight plates are defined on the breadboard:
   - plate0 to plate2: 2 mL vials (54 vials each)
   - plate3 to plate5: 20 mL bottles (8 bottles each)
   - plate6: 100 mL jars (2 jars)
   - plate7: 1.5 mL tubes (20 tubes)

3. **Coordinate Generation**:
   Generates Cartesian coordinates for all wells/containers
   on each plate using interpolation from corner coordinates
   loaded from a configuration JSON (`brb.json`).

4. **Container Assignment & Substance Tracking**:
   Each container can be assigned:
   - Coordinates (xy)
   - Substance identity, volume, solvent, and liquid surface height
   - Safety and physical parameters for liquid handling (e.g. jet height, immersion depth)

5. **Deck Geometry for Pipette Tips**:
   Support for different pipette tip geometries (50 µL, 300 µL, 1000 µL, balance tips).
   Tip racks can be initialized and reloaded dynamically.

6. **Tip Rack Management**:
   - Load or reload tip racks with updated positions and availability.
   - Mark tips as used after pipetting steps.
   - Manage tip status persistently via `tip_rack.json`.

7. **Robot Selection**:
   Robot setup can be auto-detected via environment variable `PIPETTER_NAME`
   or selected manually via a PySimpleGUI prompt.
   Supported robots: `Robowski#1`, `Robowski#2`.

8. **Data Structures**:
   - `Container`: stores the metadata of a single container.
   - `Plate`: holds a collection of containers and provides substance management.
   - `Deck_para`: defines geometry and movement parameters for tip decks.
   - `Liquid`: stores physical properties of liquids used in experiments.

Environment Variables:
----------------------
- `PIPETTER_NAME`: Robot identifier (`Robowski#1` or `Robowski#2`)
- `ROBOCHEM_DATA_PATH`: Base directory for experiment and robot configuration files

External Dependencies:
----------------------
- numpy
- json
- os
- copy
- logging
- PySimpleGUI (for GUI-based robot selection)

Usage:
------
- Run this module as a script to reload tip racks manually.
- Use `plate_on_breadboard()` to instantiate all plate objects with assigned containers and positions.
- Call `load_new_tip_rack('300ul')` or similar to update tip availability when racks are changed.

Note:
-----
Make sure that `ROBOCHEM_DATA_PATH` and `PIPETTER_NAME` are correctly set before using this module.

Author: Yankai Jia
"""

import logging
# create logger
module_logger = logging.getLogger('pipette_calibration.breadboard')

from dataclasses import dataclass
import numpy as np, copy, json, os

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

## choose a robot from pysimplegui
def choose_robot_by_pysimplegui():
    import PySimpleGUI as sg
    sg.theme('DarkAmber')  # Add a touch of color
    # All the stuff inside your window.
    layout = [[sg.Text('Choose a robot', size=(50, 1), font=("Helvetica", 30), justification='center')],
              [sg.Button('Roboski#1', size=(20, 1), font=("Helvetica", 25, "bold")),sg.Button('Roboski#2', size=(20, 1), font=("Helvetica", 25, "bold"))],
              [sg.Button('Cancel', size=(10, 1), font=("Helvetica", 25, "bold"))]
             ]
    # Create the Window
    window = sg.Window('Choose a robot', layout, auto_size_text= True, auto_size_buttons= True,size = (800, 300), element_justification='c')
    # Event Loop to process "events" and get the "values" of the inputs
    while True:
        event, values = window.read()
        # print(event, values)
        if event == sg.WIN_CLOSED or event == 'Cancel':  # if user closes window or clicks cancel
            window.close()
            return None
        elif event == 'Roboski#1':
            window.close()
            return 'Roboski#1'
        elif event == 'Roboski#2':
            window.close()
            return 'Roboski#2'

# the variables robot_name, CONFIG_PATH  and STATUS_PATH are global and are used in other modules.
try:
    robot_name = os.environ['PIPETTER_NAME']
except KeyError:
    raise ValueError('Environment variable PIPETTER_NAME is not set. Choose robowski by GUI.')
    robot_name = choose_robot_by_pysimplegui()

if robot_name == 'Roboski#1':
    CONFIG_PATH = 'config//roboski1//'
    STATUS_PATH = data_folder + "\\pipetter_files\\roboski1\\"
elif robot_name == 'Robowski#2':
    CONFIG_PATH = 'config//roboski2//'
    STATUS_PATH = data_folder + "\\pipetter_files\\roboski2\\"
else:
    raise ValueError('No robot is chosen. Please choose a robot.')

print(f'CONFIG_PATH = {CONFIG_PATH}')
assert CONFIG_PATH in ['config//roboski1//','config//roboski2//'], 'CONFIG_PATH is not correct.'

# load config file from json
with open(CONFIG_PATH + 'brb.json', 'r') as config_file:
    config_brb = json.load(config_file)

floor_z = 2220
ZeusTraversePosition = 880
balance_traverse_height = ZeusTraversePosition

bottom_z_of_vial_2ml = 2210
bottom_z_of_well_bio = 2190
bottom_z_of_bottle_20ml = 2175
bottom_z_of_jar_100ml = 2120
bottom_z_of_tube_1500ul = 2190
bottom_z_of_balance_cuvette = 1550
bottom_z_of_nanodrop_pedestal = 1300

source_substance_containers: list = []

@dataclass
class Container:
    name: str
    container_id: str
    # contaniner geometry table for zeus
    containerGeometryTableIndex: int
    container_shape: str
    diameter: int
    bottomHeight: int
    bottomSection: int
    bottomPosition: int
    immersionDepth: int
    leavingHeight: int
    jetHeight: int
    startOfHeightBottomSearch: int
    dispenseHeightAfterBottomSearch: int
    # for liquid transfer
    liquid_volume: float
    volume_max: float
    area: float  # container's horizontal cross-section area is in square mm
    min_z: float  # location of container's bottom above the floor in mm
    top_z: float  # height of container
    safety_margin_for_lldsearch_position: int
    solvent: str
    liquid_surface_height: int
    liquid_surface_height_real_value_float: float


    # coordinate
    xy: tuple = (0, 0)

    substance: str = ' '
    substance_density: float = 1.0
    mode = 'empty'


vial_2ml = Container(
    name='vial_2ml',
    containerGeometryTableIndex=0,
    container_shape="cylindrical",
    # diameter = 98, ## old value 20230322
    diameter= 118, ## new value, 20230322
    bottomHeight=0,
    bottomSection=10000,
    bottomPosition=bottom_z_of_vial_2ml,
    immersionDepth=20,
    leavingHeight= 20,
    jetHeight=50*2,
    startOfHeightBottomSearch=0,
    dispenseHeightAfterBottomSearch=80,

    liquid_volume=0,
    volume_max=2000,
    area=75.7,  # mm^2, wall thickness=0.88mm, OD=11.58, ID = 11.58-0.88*2
    min_z=floor_z - 2172,
    top_z=32,
    safety_margin_for_lldsearch_position=40,
    solvent='',
    xy=(0, 0),  # coordinate
    substance='',
    substance_density=1.0,
    container_id='',
    liquid_surface_height = bottom_z_of_vial_2ml - 100,
    liquid_surface_height_real_value_float=bottom_z_of_vial_2ml - 100
)

well_bio = Container(
    name='well_bio',
    containerGeometryTableIndex=4,
    container_shape='cylindrical',
    diameter=68,
    bottomHeight=0,
    bottomSection=10000,
    bottomPosition=bottom_z_of_well_bio,
    immersionDepth=20,
    leavingHeight=20,
    jetHeight=60,
    startOfHeightBottomSearch=30,
    dispenseHeightAfterBottomSearch=80,

    liquid_volume=0,
    volume_max=200,
    area=6.8 * 6.8 * 3.14 / 4,  # container's horizontal cross-section area is in square mm
    min_z=3.4,  # location of container's bottom above the floor in mm
    top_z=12,
    safety_margin_for_lldsearch_position=40,
    solvent='',
    xy=(0, 0),  # coordinate
    substance='',
    substance_density=1.0,
    container_id='',
    liquid_surface_height = bottom_z_of_well_bio - 100,
    liquid_surface_height_real_value_float=bottom_z_of_well_bio - 100

)

bottle_20ml = Container(
    name='bottle_20ml',
    containerGeometryTableIndex=1,
    container_shape='cylindrical',
    diameter=255,
    bottomHeight=0,
    bottomSection=10000,
    bottomPosition=bottom_z_of_bottle_20ml,
    immersionDepth=30,
    leavingHeight=30,
    jetHeight=130,
    startOfHeightBottomSearch=20,
    dispenseHeightAfterBottomSearch=50,

    liquid_volume=0,
    volume_max=20000,
    area=510.7,
    min_z=(floor_z - 2165) / 10,
    top_z=62,
    safety_margin_for_lldsearch_position=40,
    solvent='',
    xy=(0, 0),  # coordinate
    substance='',
    substance_density=1.0,
    container_id='',
    liquid_surface_height = bottom_z_of_bottle_20ml - 100,
    liquid_surface_height_real_value_float=bottom_z_of_bottle_20ml - 100)

jar_100ml = Container(
    name='jar_100ml',
    containerGeometryTableIndex=2,
    container_shape='cylindrical',
    diameter=520,  # ID of tube
    bottomHeight=0,
    bottomSection=10000,
    bottomPosition=bottom_z_of_jar_100ml,
    immersionDepth=30,
    leavingHeight=40,
    jetHeight=130,
    startOfHeightBottomSearch=50,
    dispenseHeightAfterBottomSearch=50,

    liquid_volume=0,
    volume_max=100000,
    area=2123.7,
    min_z=(floor_z - 2070) / 10,
    top_z=70,
    safety_margin_for_lldsearch_position=40,
    solvent='',
    xy=(0, 0),  # coordinate
    substance='',
    substance_density=1.0,
    container_id='',
    liquid_surface_height = bottom_z_of_jar_100ml - 100,
    liquid_surface_height_real_value_float=bottom_z_of_jar_100ml - 100
    )

tube_1500ul = Container(
    name='tube_1500ul',
    containerGeometryTableIndex=5,
    container_shape='conical_1500ul',
    diameter=88,
    bottomHeight=195,
    bottomSection=0,
    bottomPosition=bottom_z_of_tube_1500ul,
    immersionDepth=20,
    leavingHeight=20,
    jetHeight=80,
    startOfHeightBottomSearch=30,
    dispenseHeightAfterBottomSearch=80,

    liquid_volume=0,
    volume_max=2000,
    area=8.8 * 8.8 * 3.14 / 4,
    min_z=floor_z - 2172,
    top_z=32,
    safety_margin_for_lldsearch_position=40,
    solvent='',
    xy=(0, 0),  # coordinate
    substance='',
    substance_density= 1.0,
    container_id='',
    liquid_surface_height = bottom_z_of_tube_1500ul - 100,
    liquid_surface_height_real_value_float=bottom_z_of_tube_1500ul - 100
    )

balance_cuvette = Container(
    name='balance_cuvette',
    containerGeometryTableIndex=3,
    container_shape='cylindrical',
    diameter=400,  # ID of tube
    bottomHeight=0,
    bottomSection=10000,
    bottomPosition=bottom_z_of_balance_cuvette,
    immersionDepth=20,
    leavingHeight=20,
    jetHeight=50,
    startOfHeightBottomSearch=30,
    dispenseHeightAfterBottomSearch=100,

    liquid_volume=0,
    volume_max=80000,
    area=40*40, # in mm^2
    min_z=(floor_z - 1590) / 10,
    top_z=70,
    safety_margin_for_lldsearch_position=40,
    solvent='water',
    xy=config_brb['balance_cuvette_xy'],# for balance: xy=(-820, -240),  # coordinate
    substance='water',
    substance_density= 1.0,
    container_id='balance_cuvette',
    liquid_surface_height = bottom_z_of_balance_cuvette - 100,
    liquid_surface_height_real_value_float=bottom_z_of_balance_cuvette - 100
    )

nanodrop_pedestal = Container(
    name='nanodrop_pedestal',
    containerGeometryTableIndex=6,
    container_shape='cylindrical',
    diameter=400,  # ID of tube
    bottomHeight=0,
    bottomSection=10000,
    bottomPosition=bottom_z_of_nanodrop_pedestal,
    immersionDepth=20,
    leavingHeight=20,
    jetHeight=0,
    startOfHeightBottomSearch=30,
    dispenseHeightAfterBottomSearch=100,

    liquid_volume=0,
    volume_max=80000,
    area=40*40, # in mm^2
    min_z=(floor_z - 1590) / 10,
    top_z=70,
    safety_margin_for_lldsearch_position=40,
    solvent='water',
    xy=config_brb['nanodrop_pedestal_xy'],# for balance: xy=(-820, -240),  # coordinate
    substance='water',
    substance_density= 1.0,
    container_id='nanodrop_pedestal',
    liquid_surface_height = 1230,
    liquid_surface_height_real_value_float=1230
    )

def generate_container_coordinates(Nwells, topleft, topright, bottomleft, bottomright):
    '''generate coordinates for all wells of a  plate from coordinates of corner wells.'''
    # left_side_wells
    coordinates = []
    xs = np.linspace(topleft[0], bottomleft[0], Nwells[0])
    ys = np.linspace(topleft[1], bottomleft[1], Nwells[0])
    left_side_wells = np.stack((xs, ys)).T

    # right side wells
    xs = np.linspace(topright[0], bottomright[0], Nwells[0])
    ys = np.linspace(topright[1], bottomright[1], Nwells[0])
    right_side_wells = np.stack((xs, ys)).T

    wells = []
    for i in range(Nwells[0]):
        xs = np.linspace(left_side_wells[i, 0], right_side_wells[i, 0], Nwells[1])
        ys = np.linspace(left_side_wells[i, 1], right_side_wells[i, 1], Nwells[1])
        wells.append(np.stack((xs, ys)).T)
    well_positions = np.vstack(wells)

    for well_index in range(well_positions.shape[0]):
        coordinates.append(tuple(well_positions[well_index, :]))

    return coordinates


plate0_vial_2mL_coordinates = generate_container_coordinates(Nwells=(6, 9),
                                                             topleft=config_brb['plate0'][0],
                                                             topright=config_brb['plate0'][1],
                                                             bottomleft=config_brb['plate0'][2],
                                                             bottomright=config_brb['plate0'][3])


plate1_vial_2mL_coordinates = generate_container_coordinates(Nwells=(6, 9),
                                                             topleft=config_brb['plate1'][0],
                                                             topright=config_brb['plate1'][1],
                                                             bottomleft=config_brb['plate1'][2],
                                                             bottomright=config_brb['plate1'][3])

plate2_vial_2mL_coordinates = generate_container_coordinates(Nwells=(6, 9),
                                                             topleft=config_brb['plate2'][0],
                                                             topright=config_brb['plate2'][1],
                                                             bottomleft=config_brb['plate2'][2],
                                                             bottomright=config_brb['plate2'][3])

plate3_bottle_20ml_coordinates = generate_container_coordinates(Nwells=(2, 4),
                                                             topleft=config_brb['plate3'][0],
                                                             topright=config_brb['plate3'][1],
                                                             bottomleft=config_brb['plate3'][2],
                                                             bottomright=config_brb['plate3'][3])

plate4_bottle_20ml_coordinates = generate_container_coordinates(Nwells=(2, 4),
                                                             topleft=config_brb['plate4'][0],
                                                             topright=config_brb['plate4'][1],
                                                             bottomleft=config_brb['plate4'][2],
                                                             bottomright=config_brb['plate4'][3])


plate5_bottle_20ml_coordinates = generate_container_coordinates(Nwells=(2, 4),
                                                                topleft=config_brb['plate5'][0],
                                                                topright=config_brb['plate5'][1],
                                                                bottomleft=config_brb['plate5'][2],
                                                                bottomright=config_brb['plate5'][3])

plate6_jar_100ml_coordinates = [config_brb['plate6'][0],config_brb['plate6'][1]]

plate7_tube_1500ul_coordinates = [config_brb['plate7'][0],config_brb['plate7'][1]]

# return a list of container(object) in one plate.
# this function puts geometry and coordinate of containers together into one specific plate
def container_list(container_geom: object, container_coordinates) -> list:
    # input exp: vial_2ml, plate0_vial_2mL_coordinates
    container_list = []
    for container_index in range(len(container_coordinates)):
        container_geom.xy = container_coordinates[container_index]
        container_temp = copy.copy(container_geom)
        container_list.append(container_temp)
    return container_list

class Plate:
    def __init__(self, plate_id: str, containers=None):

        self.containers = containers if containers is not None else []
        self.plate_id = plate_id

        self.logger = logging.getLogger('pipette_calibration.breadboard.Plate')
        self.logger.debug(f'A Plate object is created with plate_id: {plate_id}.')

    def add_container(self, container_list: list):
        for container in container_list:
            if container not in self.containers:
                self.containers.append(container)

    def add_substance_to_container(self,
                                   substance_name: str,
                                   container_id: int,
                                   # liquid_volume: int, # calculate from liquid_surface_height not from user input
                                   solvent: str,
                                   substance_density: float,
                                   liquid_surface_height: int):
        liquid_volume_in_container = (self.containers[container_id].bottomPosition - liquid_surface_height)/10 \
                                     * self.containers[container_id].area # in ul. 1 mm^3 is 1 ul.
        self.containers[container_id].substance = substance_name
        self.containers[container_id].liquid_surface_height = liquid_surface_height
        self.containers[container_id].liquid_volume = round(liquid_volume_in_container, 1)
        self.containers[container_id].solvent = solvent
        self.containers[container_id].substance_density = float(substance_density)
        self.logger.info(f'Container {container_id} is filled with {substance_name} in {solvent} solvent. ')


    def assign_container_id(self, plate_id: int):
        for container_index in range(len(self.containers)):
            # print(f'brb.plate_list[{plate_id}].containers[{container_index}]')
            self.containers[container_index].id = {'plate_id':plate_id, 'container_id': container_index}

def plate_on_breadboard():
    plate0_containers = container_list(vial_2ml, plate0_vial_2mL_coordinates)
    plate0 = Plate(plate_id='plate0', containers=plate0_containers)
    plate0.assign_container_id(plate_id=0)

    plate1_containers = container_list(vial_2ml, plate1_vial_2mL_coordinates)
    plate1 = Plate(plate_id='plate1', containers=plate1_containers)
    plate1.assign_container_id(plate_id=1)

    plate2_containers = container_list(vial_2ml, plate2_vial_2mL_coordinates)
    plate2 = Plate(plate_id='plate2', containers=plate2_containers)
    plate2.assign_container_id(plate_id=2)

    plate3_containers = container_list(bottle_20ml, plate3_bottle_20ml_coordinates)
    plate3 = Plate(plate_id='plate3', containers=plate3_containers)
    plate3.assign_container_id(plate_id=3)

    plate4_containers = container_list(bottle_20ml, plate4_bottle_20ml_coordinates)
    plate4 = Plate(plate_id='plate4', containers=plate4_containers)
    plate4.assign_container_id(plate_id=4)

    plate5_containers = container_list(bottle_20ml, plate5_bottle_20ml_coordinates)
    plate5 = Plate(plate_id='plate5', containers=plate5_containers)
    plate5.assign_container_id(plate_id=5)

    plate6_containers = container_list(jar_100ml, plate6_jar_100ml_coordinates)
    plate6 = Plate(plate_id='plate6', containers=plate6_containers)
    plate6.assign_container_id(plate_id=6)

    plate7_containers = container_list(tube_1500ul, plate7_tube_1500ul_coordinates)
    plate7 = Plate(plate_id='plate7', containers=plate7_containers)
    plate7.assign_container_id(plate_id=7)

    module_logger.info('All plates in breadboard are created.')
    print('All plates in breadboard are created/updated.')

    return plate0, plate1, plate2, plate3, plate4, plate5, plate6, plate7

plate0, plate1, plate2, plate3, plate4, plate5, plate6, plate7 = plate_on_breadboard()
plate_list = [plate0, plate1, plate2, plate3, plate4, plate5, plate6, plate7]




@dataclass
class Deck_para:
    deckGeometryTableIndex: int
    endTraversePosition: int
    beginningofTipPickingPosition: int
    positionofTipDepositProcess: int

deckgeom_300ul = Deck_para(deckGeometryTableIndex=0, endTraversePosition=ZeusTraversePosition,
                      beginningofTipPickingPosition=config_brb['beginningofTipPickingPosition'], positionofTipDepositProcess=config_brb['positionofTipDepositProcess'])

deckgeom_1000ul = Deck_para(deckGeometryTableIndex=1, endTraversePosition=ZeusTraversePosition,
                       beginningofTipPickingPosition=config_brb['beginningofTipPickingPosition'], positionofTipDepositProcess=config_brb['positionofTipDepositProcess'])

deckgeom_balance = Deck_para(deckGeometryTableIndex=2, endTraversePosition=balance_traverse_height,
                        beginningofTipPickingPosition=config_brb['beginningofTipPickingPosition_for_balance'], positionofTipDepositProcess=config_brb['positionofTipDepositProcess_for_balance'])

deckgeom_50ul = Deck_para(deckGeometryTableIndex=3, endTraversePosition=ZeusTraversePosition,
                     beginningofTipPickingPosition=config_brb['beginningofTipPickingPosition'],positionofTipDepositProcess=config_brb['positionofTipDepositProcess'])  # this is the same as 300ul tips


deckGeometryTableIndex = {'300ul': deckgeom_300ul.deckGeometryTableIndex,
                          '1000ul': deckgeom_1000ul.deckGeometryTableIndex,
                          'balance': deckgeom_balance.deckGeometryTableIndex,
                          '50ul': deckgeom_50ul.deckGeometryTableIndex}


# decks are for pipetting tips
def generate_deck_coordinates(Nwells, topleft, topright, bottomleft, bottomright):

    # left_side_wells
    xs = np.linspace(topleft[0], bottomleft[0], Nwells[0])
    ys = np.linspace(topleft[1], bottomleft[1], Nwells[0])
    left_side_wells = np.stack((xs, ys)).T

    # right side wells
    xs = np.linspace(topright[0], bottomright[0], Nwells[0])
    ys = np.linspace(topright[1], bottomright[1], Nwells[0])
    right_side_wells = np.stack((xs, ys)).T

    wells = []
    for i in range(Nwells[0]):
        xs = np.linspace(left_side_wells[i, 0], right_side_wells[i, 0], Nwells[1])
        ys = np.linspace(left_side_wells[i, 1], right_side_wells[i, 1], Nwells[1])
        wells.append(np.stack((xs, ys)).T)
    return np.vstack(wells)


def create_deck(template_well, Nwells, topleft, topright, bottomleft, bottomright):
    well_positions = generate_deck_coordinates(Nwells, topleft, topright, bottomleft, bottomright)
    deck = {'tips': list()}
    for well_index in range(well_positions.shape[0]):
        deck['tips'].append(template_well.copy())
        deck['tips'][-1]['xy'] = list(well_positions[well_index, :])
    return deck


def load_new_tip_rack(rack_reload):
    # tip_rack = {}
    with open(STATUS_PATH + 'tip_rack.json') as json_file:
        tip_rack = json.load(json_file)

    tip = {'300ul': {'tip_vol': 300,
                     'xy': (300, -100), # dummie coordinates
                     'tipTypeTableIndex': 4,
                     'deckGeometryTableIndex': 0,
                     'ZeusTraversePosition': 880,
                     'exists': True,
                     'substance': 'None'
                     },
           '1000ul': {'tip_vol': 1000,
                      'xy': (-296.5, -32.5), # dummie coordinates
                      'tipTypeTableIndex': 6,
                      'deckGeometryTableIndex': 1,
                      'ZeusTraversePosition': 880,
                      'exists': True,
                      'substance': 'None'
                      },
           '50ul': {'tip_vol': 50,
                    'xy': (300, -100), # dummie coordinates
                    'tipTypeTableIndex': 2,
                    'deckGeometryTableIndex': 0,
                    'ZeusTraversePosition': 880,
                    'exists': True,
                    'substance': 'None'
                    },
           }
    if rack_reload == '50ul':
        print(f"the 50ul rack corner coords: {config_brb['rack_50ul'][0]},{config_brb['rack_50ul'][1]},{config_brb['rack_50ul'][2]},{config_brb['rack_50ul'][3]}")

        tip_rack['50ul'] = create_deck(template_well=tip['50ul'],
                                       Nwells=(8, 12),
                                       topleft=config_brb['rack_50ul'][0],
                                       topright=config_brb['rack_50ul'][1],
                                       bottomleft=config_brb['rack_50ul'][2],
                                       bottomright=config_brb['rack_50ul'][3],
                                       )

    if rack_reload == '300ul':

        print(f"the 300ul rack corner coords: {config_brb['rack_300ul'][0]},{config_brb['rack_300ul'][1]},{config_brb['rack_300ul'][2]},{config_brb['rack_300ul'][3]}")

        tip_rack['300ul'] = create_deck(template_well=tip['300ul'],
                                        Nwells=(8, 12),
                                        topleft=config_brb['rack_300ul'][0],
                                        topright=config_brb['rack_300ul'][1],
                                        bottomleft=config_brb['rack_300ul'][2],
                                        bottomright=config_brb['rack_300ul'][3],
                                        )
    if rack_reload == '1000ul':
        print(f"the1000ul rack corner coords: {config_brb['rack_1000ul'][0]},{config_brb['rack_1000ul'][1]},{config_brb['rack_1000ul'][2]},{config_brb['rack_1000ul'][3]}")

        tip_rack['1000ul'] = create_deck(template_well=tip['1000ul'],
                                         Nwells=(8, 12),
                                         topleft=config_brb['rack_1000ul'][0],
                                         topright= config_brb['rack_1000ul'][1],
                                         bottomleft=config_brb['rack_1000ul'][2],
                                         bottomright=config_brb['rack_1000ul'][3],
                                         )

    with open(STATUS_PATH + 'tip_rack.json', 'w', encoding='utf-8') as f:
        json.dump(tip_rack, f, ensure_ascii=False, indent=4)
        print("tip_rack.json is updated.")

    return tip_rack

def mark_next_n_tip_as_used(tip_type, n):
    with open(STATUS_PATH + 'tip_rack.json') as json_file:
        tip_rack = json.load(json_file)
    i = 0

    while i < n:
        for tip in tip_rack[tip_type]['tips']:
            if tip['exists'] == True:
                tip['exists'] = False
                break
        i+=1

    # save the revised tip_rack to jason
    with open(STATUS_PATH + 'tip_rack.json', 'w', encoding='utf-8') as f:
        json.dump(tip_rack, f, ensure_ascii=False, indent=4)


with open(STATUS_PATH + 'tip_rack.json') as json_file:
    tip_rack = json.load(json_file)

tip_rack_50ul = tip_rack['50ul']
tip_rack_300ul = tip_rack['300ul']
tip_rack_1000ul = tip_rack['1000ul']

@dataclass
class Liquid:
    name: str
    liquid_class_index: int
    density: float

    def __post_init__(self):
      pass

    def __repr__(self):
        return f'liquid_name: {self.name}, liquid_class_index: {self.liquid_class_index}'

if __name__ == "__main__":

    print('This is main.')

    ## run this ONLY when changing new tip rack.
    # load_new_tip_rack(rack_reload ='300ul')
    # module_logger.info('New tip rack: 300ul is loaded.')
    # #
    # load_new_tip_rack(rack_reload ='1000ul')
    # module_logger.info('New tip rack: 1000ul is loaded.')
    #
    # load_new_tip_rack(rack_reload ='50ul')
    # module_logger.info('New tip rack: 50ul is loaded.')

    # plate_on_breadboard()
