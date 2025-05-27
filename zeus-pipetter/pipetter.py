"""
Pipetter Control System

This module defines the control logic for automated liquid handling operations using a pipetting robot
composed of a Zeus vertical axis and a 2D XY gantry system. It provides full functionality for pipetting,
tip handling, liquid level detection, and integration with an analytical balance.

Core Components:
----------------
- **Gantry**: Controls XY-stage motion for accurate positioning under the robot head (Zeus).
- **Pipetter**: Orchestrates actions such as pick_tip, discard_tip, transfer_liquid, and interactions with a balance.

Key Features:
-------------
- Move to coordinates with collision-aware logic (e.g., transitions between balance and nanodrop zones)
- Tip handling with JSON-based state persistence (`tip_rack.json`)
- Multi-step pipetting for volumes exceeding tip capacity, with optional calibration
- Automatic volume detection using Zeus sensors
- Integration with serial-connected balance (tare, read value, internal adjustment)
- Transfer tracking and error logging
- Optional prewetting for improved accuracy

Typical Workflow:
-----------------
1. Initialize hardware (Zeus + Gantry + Pipetter)
2. Pick or change tips
3. Detect or input liquid levels
4. Transfer or aspirate/dispense liquids
5. Measure and record volumes using the balance for pipetting calibration
"""

import logging

pipetter_logger = logging.getLogger('main.pipetter')

import copy, json,time, numpy as np, serial,re, winsound, sound, asyncio, os


import breadboard as brb
from calibration import calibrations as calib

## get config file path from breadboard
CONFIG_PATH = brb.CONFIG_PATH
STATUS_PATH = brb.STATUS_PATH

# load config file from json
with open(CONFIG_PATH + 'pipetter.json', 'r') as config_file:
    config_pt = json.load(config_file)

robot_name = os.environ['PIPETTER_NAME']


class Gantry():
    """
    The xy gantry moves with Zeus on it.

    Gantry() take zeus object as argument, which is used to request position of Z drive.
    Only when the Z drive position is in safe traverse height will the gantry be able to move.
    """

    xy_position = ((0, 0)) # this is to store the gantry position after every move.

    def __init__(self,
                 zeus: object, # pass the zeus module to gantry, this is for checking traverse height,
                 max_x: int = config_pt['gantry']['max_x'], # -820,
                 max_y: int = config_pt['gantry']['max_y'], #-360,
                 horiz_speed: int = config_pt['gantry']['horiz_speed'],#200*60,# horizontal speed in mm / min
                 xy_offset: tuple = config_pt['gantry']['xy_offset'], #(-2.5, 0),# offsets in x and y, negative to right, closer; positive, to left, further
                 zeus_traverse_position: int = config_pt['gantry']['zeus_traverse_position'], #880,
                 trash_xy: tuple = config_pt['gantry']['trash_xy'], #(-500, -70),
                 idle_xy: tuple = config_pt['gantry']['idle_xy'], #(-500, -220),
                 idle_xy_balance:tuple = config_pt['gantry']['idle_xy_balance'], #(-500, -220),

                 ):
        self.logger = logging.getLogger('main.gantry.Gantry')
        self.logger.info('gantry is initiating...')
        print('gantry is initiating...')
        self.serial = serial.Serial('COM5', 115200, timeout=0.2)
        self.horiz_speed = horiz_speed # horizontal speed in mm/min
        self.xy_offset = xy_offset
        self.max_x = max_x
        self.max_y = max_y

        self.trash_xy = trash_xy
        self.idle_xy = idle_xy
        self.idle_xy_balance = idle_xy_balance


        self.zm = zeus
        self.zeus_traverse_position = zeus_traverse_position
        # self.home_xy()
        self.left_balance = config_pt['gantry']['balance_left_boundary']
        self.upper_balance = config_pt['gantry']['balance_chamber_upper_boundary']
        self.lower_balance = config_pt['gantry']['balance_chamber_lower_boundary']
        print('gantry is initiated!')

    def send_to_xy_stage(self, command, wait_for_ok=True, verbose=False, read_all=False,
                         ensure_traverse_height=True) -> None:
        ser = self.serial
        ser.write(str.encode(command + '\r\n'))
        # ser.write(str.encode(command))
        if verbose:
            print('SENT: {0}'.format(command))
        # time.sleep(1)

        if wait_for_ok:
            if verbose:
                print('Waiting for ok...')
            while True:
                line = ser.readline()
                if b'Alarm' in line:
                    print('GRBL ALARM: GRBL wend into alarm. Overrode it with $X.')
                    self.send_to_xy_stage('$X')
                    break
                if verbose:
                    print(line)
                if b'ok' in line:
                    break

        if read_all:
            if verbose:
                print('Reading all...')
            while True:
                line = ser.readline()
                if verbose:
                    print(line)
                if line == b'':
                    break

    def configure_grbl(self) -> None:
        with open(CONFIG_PATH + "\\grbl_settings.txt", 'r') as grbl_config_file:
            for line in grbl_config_file:
               self.send_to_xy_stage(command = line.split('    (')[0], read_all= True, verbose= True)
        print("XY stage configured!")

    def xy_pos(self) -> None:
        self.send_to_xy_stage(command= '?', read_all= True, verbose= False)

    def time_that_xy_motion_takes(self, dx: int, dy: int, acceleration=2000, max_speed=333.33333):

        travel_times = []
        for distance in [abs(dx), abs(dy)]:
            halfdistance = distance / 2
            # constant acceleration scenario
            constant_acceleration_halftime = np.sqrt(halfdistance * 2 / acceleration)
            speed_at_midpoint = constant_acceleration_halftime * acceleration
            if speed_at_midpoint <= max_speed:
                time_here = constant_acceleration_halftime * 2
            else:
                # this means that the stage reaches max speed before midpoint and then
                #   continues at his max speed
                constant_acceleration_halftime = max_speed / acceleration
                dist_traveled_at_constant_acceleration = acceleration * (constant_acceleration_halftime ** 2) / 2
                distance_traveled_at_constant_speed = halfdistance - dist_traveled_at_constant_acceleration
                const_speed_halftime = distance_traveled_at_constant_speed / max_speed
                time_here = 2 * (constant_acceleration_halftime + const_speed_halftime)
            travel_times.append(time_here)
        # print(max(travel_times))
        return max(travel_times)

    def get_location(self, x, y, location=None):
        if x >= -460:
            location = 'out'
        else:
            if y >= -255:
                location = 'in_balance'
            else:
                location = 'in_nanodrop'

        return location

    def move_to_one_position(self, target_xy,
                             ensure_traverse_height,
                             block_until_motion_is_completed,
                             use_time_estimate,
                             verbose = False):

        zeus_at_traverse_height = self.zm.pos <= self.zeus_traverse_position
        if ensure_traverse_height and not zeus_at_traverse_height:
            print(f'ERROR: ZEUS was not in traverse height before motion, but instead at {self.zm.pos}.\n'
                  f'Motion aborted!')

            raise ValueError('ZEUS was not in traverse height before motion')

        self.send_to_xy_stage(
            command='G0 X{0:.3f} Y{1:.3f}'.format(target_xy[0], target_xy[1]),
            read_all=False, ensure_traverse_height=ensure_traverse_height)

        if block_until_motion_is_completed:
            if use_time_estimate:
                time.sleep(self.time_that_xy_motion_takes(dx=target_xy[0] - self.xy_position[0],
                                                          dy=target_xy[1] - self.xy_position[1]))
            else:
                t0 = time.time()
                time.sleep(0.1)
                finished_moving = False
                for i in range(100):
                    if finished_moving:
                        break
                    if verbose:
                        print(f'Status read {i}')
                    self.serial.write(str.encode('?' + '\r\n'))
                    while True:
                        line = self.serial.readline()
                        if verbose:
                            print(line)
                        if b'Idle' in line:
                            finished_moving = True
                        if line == b'':
                            break
                # print(f'{time.time() - t0}')
                if verbose:
                    print('Finished moving xy stage')

    def move_xy(self, xy: tuple,
                verbose=False,
                ensure_traverse_height=True,
                block_until_motion_is_completed=True,
                use_time_estimate=True) -> None:

        # print(f'Moving to {xy}...')
        # soft limit for x
        if xy[0] < self.max_x or xy[0] > 0:
            self.logger.error(f'XY STAGE ERROR: target X is beyond the limit ({self.max_x}, 0). Motion aborted.')
            raise ValueError(f'XY STAGE ERROR: target X is beyond the limit ({self.max_x}, 0). Motion aborted.')
        # soft limit for y
        if xy[1] < self.max_y or xy[1] > 0:
            self.logger.error(f'XY STAGE ERROR: target Y is beyond the limit ({self.max_y}, 0). Motion aborted.')
            raise ValueError(f'XY STAGE ERROR: target Y is beyond the limit ({self.max_y}, 0). Motion aborted.')

        current_location = self.get_location(self.xy_position[0], self.xy_position[1])
        target_location = self.get_location(xy[0], xy[1])
        transition_sequence = None

        # moving between: out <-> in_balance
        if {current_location, target_location} == {'out', 'in_balance'}:
            transition_location = self.idle_xy_balance
        # moving between: out <-> in_nanodrop
        elif {current_location, target_location} == {'out', 'in_nanodrop'}:
            transition_location = self.idle_xy
        # moving between: in_balance <-> in_nanodrop
        elif {current_location, target_location} == {'in_balance','in_nanodrop'}:
            transition_location = 'warning'
            if current_location == 'in_balance':
                transition_sequence = (self.idle_xy_balance, self.idle_xy)
            elif current_location == 'in_nanodrop':
                transition_sequence = (self.idle_xy, self.idle_xy_balance)
        else:
            # no need to go to transition location.
            transition_location = None

        if transition_location == self.idle_xy_balance or transition_location == self.idle_xy:
            # move to transition_location first
            self.move_to_one_position(target_xy=transition_location,
                                      ensure_traverse_height=ensure_traverse_height,
                                      block_until_motion_is_completed=block_until_motion_is_completed,
                                      use_time_estimate=use_time_estimate,
                                      verbose=verbose)

        elif transition_location == 'warning':
            x = input('Moving between Balance and Nanodrop, continue?')
            if x in ['y', 'Y']:
                    self.move_to_one_position(target_xy=transition_sequence[0],
                                              ensure_traverse_height=ensure_traverse_height,
                                              block_until_motion_is_completed=block_until_motion_is_completed,
                                              use_time_estimate=use_time_estimate,
                                              verbose=verbose)
                    self.move_to_one_position(target_xy=transition_sequence[1],
                                              ensure_traverse_height=ensure_traverse_height,
                                              block_until_motion_is_completed=block_until_motion_is_completed,
                                              use_time_estimate=use_time_estimate,
                                              verbose=verbose)
            else:
                raise ValueError('Moving between Balance and Nanodrop stopped due to possible collision.')
        elif transition_location is None:
            pass

        else:
            raise ValueError('transition_location not defined.')

        self.move_to_one_position(target_xy = (xy[0] + self.xy_offset[0], xy[1] + self.xy_offset[1]),
                                  ensure_traverse_height=ensure_traverse_height,
                                  block_until_motion_is_completed=block_until_motion_is_completed,
                                  use_time_estimate=use_time_estimate,
                                  verbose=verbose)

        self.xy_position = xy

    def home_xy(self, ensure_traverse_height=True) -> None:
        self.logger.info('The gantry is homing...')
        self.zm.move_z(self.zm.ZeusTraversePosition)
        self.send_to_xy_stage(command = '$H', read_all=True, verbose=False,
                              ensure_traverse_height=ensure_traverse_height)
        self.xy_pos()

    def kill_alarm(self) -> None:
        self.send_to_xy_stage("$X", read_all= True, verbose= True)

    def close_gantry(self)-> None:
        time.sleep(2)
        self.serial.close()

    def view_grbl_settings(self)-> None:
        self.send_to_xy_stage('$$', read_all=True, verbose=True)
        self.xy_pos()

    def move_through_wells(self, plate: object, dwell_time=0.1, ensure_traverse_height=True):
        for container in plate.containers:
            print(f'This is well index: {container}')
            self.move_xy(container.xy, ensure_traverse_height=ensure_traverse_height)
            time.sleep(dwell_time)

    def move_to_trash_bin(self, ensure_traverse_height=True):
        self.move_xy(self.trash_xy, ensure_traverse_height=ensure_traverse_height)

    def move_to_idle_position(self, ensure_traverse_height=True):
        self.move_xy(self.idle_xy, ensure_traverse_height=ensure_traverse_height)


class Pipetter():

    def __init__(self,
                 zeus: object,
                 gantry: object,
                 robot_name:str = robot_name,
                 is_balance_involved: bool = False,
                 ):
        self.zeus = zeus
        self.gantry = gantry
        self.robot_name = robot_name
        if is_balance_involved:
            self.balance = serial.Serial(port=config_pt["balance_port"]["port"],
                                         baudrate=config_pt["balance_port"]["baudrate"],
                                         stopbits=serial.STOPBITS_ONE,
                                         parity=serial.PARITY_NONE,
                                         timeout=config_pt["balance_port"]["timeout"])
            print('A balance is installed.')
        else:
            print('No balance is installed.')

        self.logger = logging.getLogger('main.pipetter.Pipetter')
        self.logger.info('creating an instance of Pepetter')

    def home_xy(self):
        self.gantry.home_xy()

    def beep_n(self):
        duration = 600  # milliseconds
        freq = 1000  # Hz
        # time.sleep(0.2)
        for i in range(10):
            winsound.Beep(freq, duration)

    def pick_tip(self, tip_type: str):

        with open(STATUS_PATH + 'tip_rack.json') as json_file:
            tip_rack = json.load(json_file)

        self.zeus.move_z(config_pt['gantry']['zeus_traverse_position'])
        self.zeus.wait_until_zeus_reaches_traverse_height()

        # if the rack is empty then ask user to reload
        if not any(item['exists'] for item in tip_rack[tip_type]['tips']):
            sound.beep_for_tip_changing()
            input(f'ERROR: The tip rack is empty. Please reload the tip rack and hit enter.')
            tip_rack = brb.load_new_tip_rack(rack_reload=tip_type)

        # In the rack, find the first tip that exists
        for item in tip_rack[tip_type]['tips']:
            if item['exists']:
                # pick up tip
                self.gantry.move_xy(item['xy'], ensure_traverse_height=True)
                time.sleep(0.5)
                self.zeus.pickUpTip(tipTypeTableIndex=item['tipTypeTableIndex'],
                                    deckGeometryTableIndex=item['deckGeometryTableIndex'])

                self.zeus.wait_until_zeus_responds_with_string('GTid')
                tip_is_picked = self.zeus.getTipPresenceStatus()

                if tip_is_picked:
                    # print(f'tip_status: {self.zeus.getTipPresenceStatus()}')
                    self.zeus.tip_on_zeus = tip_type
                    logging.info(f'Now the tip on zeus is : {tip_type}')
                    item['exists'] = False
                    # wait_until_zeus_reaches_traverse_height()
                    # update json file
                    with open(STATUS_PATH + 'tip_rack.json', 'w', encoding='utf-8') as f:
                        json.dump(tip_rack, f, ensure_ascii=False, indent=4)
                else:
                    raise ValueError('No tip is picked up.')
                return True
        self.logger.info('ERROR: No tips in rack.')
        raise Exception

    def discard_tip(self):
        self.zeus.move_z(config_pt['gantry']['zeus_traverse_position'])
        self.zeus.wait_until_zeus_reaches_traverse_height()
        self.gantry.move_xy(self.gantry.trash_xy)
        self.zeus.discardTip(deckGeometryTableIndex=1)
        self.zeus.tip_on_zeus = ''
        # time.sleep(0.25)
        # self.zeus.move_z(self.zeus.ZeusTraversePosition)
        self.zeus.wait_until_zeus_responds_with_string('GUid')

    def change_tip(self, tip_rack: str):
        if self.zeus.tip_on_zeus != '':
            self.discard_tip()
        self.pick_tip(tip_rack)
        # self.zeus.wait_until_zeus_responds_with_string('GUid')

    def check_volume_in_container(self, container: object,
                                  liquidClassTableIndex: int = 13, lld: int = 1,
                                  lldSearchPosition: int = 1300, liquidSurface: int = 1300,
                                  tip_for_volume_check: str = '300ul',
                                  change_tip_after_each_check: bool = True):

        if change_tip_after_each_check:
            self.change_tip(tip_for_volume_check)
        else:
            if self.zeus.tip_on_zeus != tip_for_volume_check:
                self.change_tip(tip_for_volume_check)

        self.zeus.move_z(self.zeus.ZeusTraversePosition, raise_exception=False)
        time.sleep(1)
        self.gantry.move_xy(container.xy)

        self.zeus.volumeCheck(containerGeometryTableIndex=container.containerGeometryTableIndex,
                              deckGeometryTableIndex=brb.deckGeometryTableIndex[tip_for_volume_check],
                              liquidClassTableIndex=liquidClassTableIndex,
                              lld=lld,
                              lldSearchPosition=lldSearchPosition,
                              liquidSurface=liquidSurface)

        received_msg = self.zeus.r.received_msg
        while 'yl' not in received_msg:
            time.sleep(1)
            received_msg = self.zeus.r.received_msg

        if not self.zeus.zeus_had_error(received_msg):
            # print(received_msg)
            liquid_surface = received_msg[received_msg.find('yl') + 2:received_msg.find('yl') + 6]
            # volume = received_msg[received_msg.find('aw') + 2:received_msg.find('aw') + 8] # This value is from Zeus and not precise.

            ## calculate the volume manually
            volume = ((container.bottomPosition-int(liquid_surface)) / 10)  * container.area # this is in mm^3, uL
            print(f'Volume Check done! liquid_surface: {liquid_surface}, volume: {volume}')

        else:
            print(f'Liquid level not detected')
            liquid_surface = 0
            volume = 0
            self.zeus.move_z(config_pt['gantry']['zeus_traverse_position'], raise_exception=False)
            # self.change_tip(tip_for_volume_check, raise_exception=False)
            self.change_tip(tip_for_volume_check)

        print(f'Volume Check done! liquid_surface: {liquid_surface}, volume: {volume}')
        return (int(liquid_surface), float(int(volume) / 10))  # after / 10, volume is in ul


    def draw_liquid(self, transfer_event: object, n_retries=3) -> bool:

        self.zeus.move_zeus_to_traverse_height()
        self.gantry.move_xy(transfer_event.source_container.xy)

        for retry in range(n_retries):
            try:
                # print(f'Aspiration volume for zeus: {int(round(transfer_event.aspirationVolume * 10))}')
                self.zeus.aspiration(aspirationVolume=int(round(transfer_event.transfer_volume * 10)), # volume in 0.1 ul
                                     containerGeometryTableIndex=transfer_event.asp_containerGeometryTableIndex,
                                     deckGeometryTableIndex=transfer_event.asp_deckGeometryTableIndex,
                                     liquidClassTableIndex=transfer_event.liquidClassTableIndex,
                                     qpm=transfer_event.asp_qpm,
                                     lld=transfer_event.asp_lld,
                                     liquidSurface = transfer_event.source_container.liquid_surface_height,
                                     lldSearchPosition= transfer_event.source_container.liquid_surface_height - 50,
                                     mixVolume=transfer_event.asp_mixVolume,
                                     mixFlowRate=transfer_event.asp_mixFlowRate,
                                     mixCycles=transfer_event.asp_mixCycles)
                print(f'liquid surface height in source container: {transfer_event.source_container.liquid_surface_height} ')
                self.logger.info(f'liquid surface height in source container: {transfer_event.source_container.liquid_surface_height} ')

                # Replace this sleep with a proper check. Why do you even need a sleep if next function is zeus.wait_until...???
                # time.sleep(2)
                self.zeus.wait_until_zeus_responds_with_string('GAid')

                # time.sleep(0.5)
                self.zeus.move_z(self.zeus.ZeusTraversePosition)

                return True

            except ZeusError:
                if self.zeus.zeus_error_code(self.zeus.r.received_msg) == '81':
                    # Empty tube detected during aspiration
                    self.logger.info('ZEUS ERROR: Empty tube during aspiration. Dispensing and trying again.')
                    time.sleep(2)
                    self.zeus.move_z(self.zeus.ZeusTraversePosition)
                    time.sleep(2)
                    self.dispense_liquid(transfer_event)
                    time.sleep(2)
                    continue

        self.logger.info(f'Tried {n_retries} but zeus error is still there')
        raise Exception

    def dispense_liquid(self, transfer_event: object) -> None:

        self.zeus.move_zeus_to_traverse_height()
        self.gantry.move_xy(transfer_event.destination_container.xy)


        self.zeus.dispensing(dispensingVolume=int(round(transfer_event.transfer_volume * 10)),
                         containerGeometryTableIndex=transfer_event.disp_containerGeometryTableIndex,
                         deckGeometryTableIndex=transfer_event.disp_deckGeometryTableIndex,
                         liquidClassTableIndex=transfer_event.liquidClassTableIndex,
                         lld=transfer_event.disp_lld,
                         lldSearchPosition=transfer_event.destination_container.liquid_surface_height - 50,
                         liquidSurface=transfer_event.destination_container.liquid_surface_height,
                         searchBottomMode=transfer_event.searchBottomMode,
                         mixVolume=transfer_event.disp_mixVolume,
                         mixFlowRate=transfer_event.disp_mixVolume,
                         mixCycles=transfer_event.disp_mixCycles)
        # print(f'DEBUG::dispense_liquid():: disp_liquidSurface: {transfer_event.disp_liquidSurface} ')

        # time.sleep(0.25)
        # wait_until_zeus_reaches_traverse_height()
        self.zeus.wait_until_zeus_responds_with_string('GDid')
        # time.sleep(0.5)
        self.zeus.move_z(self.zeus.ZeusTraversePosition)
        return True

    def calib(self, transfer_event:object):

        if transfer_event.solvent == 'water':
            print('The solvent is water. No calibration is needed.')
            return transfer_event

        transfer_event =copy.deepcopy(transfer_event)

        target_volume = transfer_event.transfer_volume
        setting_volume = calib.Interpolation(transfer_event.solvent).interp(target_volume,
                                                                            verbose= True)

        transfer_event.transfer_volume = setting_volume
        print(f'@@@@@@@@Volume valibrated. Target volume: {target_volume}. Setting volume: {setting_volume}')

        if transfer_event.transfer_volume <=50:
            transfer_event.tip_type = '50ul'
        elif transfer_event.transfer_volume <=300:
            transfer_event.tip_type = '300ul'
        elif transfer_event.transfer_volume<=1000:
            transfer_event.tip_type = '1000ul'


        return transfer_event

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

    def transfer_liquid(self, transfer_event: object, max_volume: int=None, use_calibration: bool = True) -> None:

        # if transfer_event.solvent == 'water': ## water is the built-in solvent that needs no software calibration
        #     use_calibration = False

        if max_volume is None:

            print(transfer_event.tip_type, transfer_event.substance)
            max_volume = int(transfer_event.tip_type[:-2]) # extract volume from tip type

        print(f'transfer volume is set to : {transfer_event.transfer_volume}ul')

        # if transfer volume exceeds max_volume, do several pipettings
        N_max_vol_pipettings = int(transfer_event.transfer_volume // max_volume)
        # print(f'N_max_vol_pipettings: {N_max_vol_pipettings}')

        for i in range(N_max_vol_pipettings):
            # print(f'Pipetting {i+1} of {N_max_vol_pipettings}')
            _split_event_1 = copy.deepcopy(transfer_event)
            _split_event_1.transfer_volume = max_volume # volume in ul

            if use_calibration:
                # this is for calibration
                _split_event_1 = self.calib(_split_event_1)

                # change tip if necessary
                if self.zeus.tip_on_zeus != _split_event_1.tip_type:
                    self.change_tip(_split_event_1.tip_type)
                    # liquid class parameter is also needed to be updated.
                    _split_event_1.liquidClassTableIndex = \
                        self.get_liquid_class_index(solvent=_split_event_1.solvent,
                                                    mode='empty',
                                                    tip_type=_split_event_1.tip_type)

            # print(f'Aspiration volume: {_split_event_1.aspirationVolume}ul ')
            # print(f'Dispensing volume: {_split_event_1.dispensingVolume}ul ')
            self.draw_liquid(_split_event_1)
            # print(f'DEBUG: transfer_liquid()::disp_height: {transfer_event.disp_liquidSurface}')
            liquid_surface_height_from_zeus = self.detect_liquid_surface()
            self.dispense_liquid(_split_event_1)

        volume_of_last_pipetting = transfer_event.transfer_volume % max_volume

        if volume_of_last_pipetting:
            _split_event_2 = copy.deepcopy(transfer_event)
            _split_event_2.transfer_volume = volume_of_last_pipetting
            _split_event_2.transfer_volume = volume_of_last_pipetting

            if use_calibration:
                _split_event_2 = self.calib(_split_event_2)

            # change tip if necessary
            if self.zeus.tip_on_zeus != _split_event_2.tip_type:
                self.change_tip(_split_event_2.tip_type)
                # liquid class parameter is also needed to be updated.
                _split_event_2.liquidClassTableIndex = self.get_liquid_class_index(solvent=_split_event_2.solvent,
                                                                                   mode='empty',
                                                                                   tip_type=_split_event_2.tip_type)

            self.draw_liquid(_split_event_2)
            liquid_surface_height_from_zeus = self.detect_liquid_surface()
            self.dispense_liquid(_split_event_2)
            time.sleep(0.5)

        # self.logger.info(f'Aspiration volume: {transfer_event.aspirationVolume}ul '
        #                  f'Dispensing volume: {transfer_event.dispensingVolume}ul')


    def detect_liquid_surface(self) -> int:
        self.zeus.sendCommand('GNid0001')
        time.sleep(0.25)
        liquid_surface_height_detected = re.findall('[0-9]+', self.zeus.r.received_msg)[1]

        # self.logger.info(f'liquid_surface_height_detected: {liquid_surface_height_detected}')
        return int(liquid_surface_height_detected)

    def send_command_to_balance(self, command, read_all=True, verbose=False):

        self.balance.write(str.encode(command + '\n'))
        while True:
            line = self.balance.readline()
            if verbose:
                print(f'Response from balance. {line}')
            if line == b'':
                break

    def open_balance_door(self):
        self.send_command_to_balance('WS 1')

    def close_balance_door(self):
        self.send_command_to_balance('WS 0')
        # time.sleep(1)
        # self.send_command_to_balance('C3\n')


    def balance_tare(self, verbose = True):
        self.balance.write(str.encode('T\n'))
        taring_complete = False
        while True:
            line = self.balance.readline()
            if b'T' in line:
                taring_complete = True
            if verbose:
                pass
                # print(f'Tare in progress. Response from balance. {line}')
            if line == b'':
                if taring_complete:
                    print('The balance is tared.')
                    break

    def balance_cancel(self,verbose = False):

        self.balance.write(str.encode('@\n'))

        cancelling_complete = False

        while True:
            line = self.balance.readline()

            if b'A' in line: cancelling_complete = True

            if verbose: print(f'Cancelling in progress. Response from balance. {line}')

            if line == b'':
                if cancelling_complete:
                    break

    def balance_cancel_all(self,verbose = False):

            self.balance.write(str.encode('C\n'))

            cancelling_complete = False

            while True:
                line = self.balance.readline()

                if b'C A' in line: cancelling_complete = True

                if verbose: print(f'Cancelling in progress. Response from balance. {line}')

                if line == b'':
                    if cancelling_complete:
                        break

    def balance_internal_adjustment(self, verbose=False):

        self.send_command_to_balance('C3\n')
        adjustment_complete = False
        time_here = time.time()
        while True:
            line = self.balance.readline()

            if b'C3 B' in line:
                print('Balance internal adjustment has started.')

            if b'C3 A' in line:
                adjustment_complete = True
                print('Balance internal adjustment has completed.')

            if b'C3 I' in line or b'C3 L' in line:
                print('Balance internal adjustment has ERROR. No respons will be further given.')
                break

            if verbose:
                print(f'Adjustment in progress. Response from balance. {line}')
                # pass

            if line == b'':
                if adjustment_complete:
                    break

            if time.time() - time_here > 120:
                print(f'Balance internal adjustment took more than 120 seconds. Aborting.')
                break

    def balance_zero(self, verbose=True):
        self.balance.write(str.encode('Z\n'))
        zeroing_complete = False
        time_stamp = time.time()
        while True:
            line = self.balance.readline()
            if b'Z' in line:
                zeroing_complete = True
            if verbose:
                print(f'Tare in progress. Response from balance. {line}')
                # pass
            if line == b'':
                if zeroing_complete:
                    print('The balanced is zeroed.')
                    break
            if time.time() - time_stamp > 10:
                self.logger.error(f'Balance zeroing took more than 10 seconds. Aborting.')
                zeroing_complete = False
                break
        return zeroing_complete

    def balance_value(self, verbose=False):
        value = 0
        self.balance.write(str.encode('SI\n'))
        measurement_successful = False
        while True:
            line = self.balance.readline()
            # print(f'balance line: {line}')
            if (b'S D' in line) or (b'S S' in line):
                raw_parsed = line.split(b' g\r\n')[0][-8:] # this is for balance XPE205
                raw_parsed = re.findall(r"[-+]?(?:\d*\.*\d+)", line.decode("utf-8") )[0] # this is for balance ME204

                if verbose:
                    print(f"Raw parsed: {raw_parsed}")
                value: float = float(raw_parsed)
                measurement_successful = True
            if b'S I' in line:
                'Balance command not executed.'
                time.sleep(5)
                self.balance.write(str.encode('SI\n'))
            if verbose:
                print(line)
            if line == b'':
                if measurement_successful:
                    break
        return value

    def balance_stable_value(self, verbose=False):
        value = 0
        self.balance.write(str.encode('S\n'))
        measurement_successful = False
        while True:
            line = self.balance.readline()
            # print(f'balance line: {line}')
            if (b'S D' in line) or (b'S S' in line):
                raw_parsed = line.split(b' g\r\n')[0][-8:] # this is for balance XPE205
                # raw_parsed = re.findall(r"[-+]?(?:\d*\.*\d+)", line.decode("utf-8") )[0] # this is for balance ME204
                # print(f'responds from balance: {line}')
                if verbose:
                    print(f"Raw parsed: {raw_parsed}")
                value: float = float(raw_parsed)
                measurement_successful = True

            if b'S I' in line:
                print('Balance command not executed.')
                time.sleep(2)
                self.balance.write(str.encode('SI\n'))

            if verbose:
                print(line)

            if line == b'':
                if measurement_successful:
                    break
        return value


    def move_to_balance(self, xy: tuple = brb.balance_cuvette.xy):
        self.open_balance_door()
        self.zeus.move_z(config_pt['gantry']['zeus_traverse_position'])
        self.gantry.move_xy(xy)

    def prewet_new_tip(self,zm: object, pt: object, pipetting_event: object, use_calibration = False):

        print('Prewetting new tip...')

        event_adjusted = copy.deepcopy(pipetting_event)

        max_volume = int(re.findall(r'\d+', zm.tip_on_zeus)[0])
        event_adjusted.transfer_volume = max_volume
        event_adjusted.destination_container = event_adjusted.source_container
        event_adjusted.disp_liquidSurface = 1800

        print(f'Prewetting tip with {max_volume}ul of {event_adjusted.substance}')
        pt.transfer_liquid(event_adjusted, use_calibration=use_calibration)
        print('Prewet done! Continue with pipetting...')

    def pipetting_to_balance_and_weight(self, zm, pt, transfer_event, use_calibration, timedelay= 2, prewet_tip = True,):
        global xy_position
        global weighted_values

        # do prewet every time a new tip is taken
        if prewet_tip:
            self.prewet_new_tip(zm=zm, pt=pt, pipetting_event=transfer_event)


        # self.close_balance_door()
        print('Waiting for balance to settle...')
        time.sleep(timedelay)

        # balance_tare()
        # self.balance_zero(verbose=False)
        # print('Balance zeroed.')

        weight_before = self.balance_value()
        print(f'weight_before: {weight_before} g')

        print('Waiting for liquid transfer...')
        self.open_balance_door()
        self.transfer_liquid(transfer_event=transfer_event,use_calibration=use_calibration)
        transfer_event.execute_event()

        time.sleep(0.5)
        self.gantry.move_to_trash_bin()
        self.close_balance_door()
        time.sleep(timedelay)
        weight_after = self.balance_value()
        print(f'weight_after: {weight_after} g')

        pipetting_weight = round((weight_after - weight_before) * 1000, 6)  # mg
        pipetting_volume = round(pipetting_weight / transfer_event.source_container.substance_density, 2)
        self.logger.info(f'Weight of liquid transferred: {pipetting_weight} mg')
        self.logger.info(f'Volume of aliquottransferred: {pipetting_volume} ul')

        sound.beep_for_each_pipetting()

        return pipetting_weight, pipetting_volume

    def pipetting_to_balance_and_weight_n_times(self, zm, pt, transfer_event, n_times=3,
                                                change_tip_after_every_pipetting:bool = False, prewet_tip = True, use_calibration: bool = True):

        # print(f'this is transfer_event: {transfer_event}')
        dict_for_one_event = {}
        transfer_volume = transfer_event.transfer_volume
        dict_for_one_event[f'{transfer_event.substance}_{transfer_volume}ul'] = \
            {'weight': [], 'volume': [], 'liquid_class_index': [], 'tip_type': []}
        temp_dict = dict_for_one_event[f'{transfer_event.substance}_{transfer_volume}ul']

        for i in range(n_times):

            print(f'this is n_times: {i}/{n_times}')
            # print(transfer_event.event_label)

            if i != 0:
                prewet_tip = False

            weight, volume = self.pipetting_to_balance_and_weight(zm=zm, pt=pt, transfer_event=transfer_event,
                                                                  prewet_tip=prewet_tip, use_calibration=use_calibration)

            if change_tip_after_every_pipetting:
                self.change_tip(transfer_event.tip_type)
                time.sleep(0.5)
                print(f'Changed tip to {transfer_event.tip_type} after {i}th pipetting.')

            temp_dict['weight'].append(weight)
            temp_dict['volume'].append(volume)
            temp_dict['liquid_class_index'].append(transfer_event.liquidClassTableIndex)
            temp_dict['tip_type'].append(transfer_event.tip_type)

            # input(f'Press enter to continue...')

        print(dict_for_one_event)

        return dict_for_one_event

    def close_ports_and_zeus(self):
        self.balance.close()
        self.gantry.serial.close()
        self.zeus.switchOff()

class ZeusError(Exception):
    pass


def initiate_hardware():
    # initiate zeus
    zm = zeus.ZeusModule(id=1)
    time.sleep(3)
    print("zeus is loaded as: zm")

    # initiate gantry
    gt = Gantry(zeus=zm)
    time.sleep(3)
    print("gantry is loaded as: gt")
    # gt.configure_grbl() # This only need to be done once.
    gt.home_xy()
    if gt.xy_position == (0, 0):
        print("gantry is homed")

    # initiate pipetter
    pt = Pipetter(zeus=zm, gantry=gt)
    time.sleep(2)
    print("pipetter is loaded as: pt")

    return zm, gt, pt


if __name__ == '__main__':
    import zeus, breadboard as brb

    # initiate zeus
    zm = zeus.ZeusModule(id=1)
    time.sleep(3)
    print("zeus is loaded as: zm")

    # initiate gantry
    gt = Gantry(zeus=zm)
    time.sleep(3)
    print("gantry is loaded as: gt")
    # gt.configure_grbl() # This only need to be done once.
    gt.home_xy()
    if gt.xy_position == (0, 0):
        print("gantry is homed")

    # initiate pipetter
    pt = Pipetter(zeus=zm, gantry=gt)
    time.sleep(2)
    print("pipetter is loaded as: pt")


    # pt.check_volume_in_container(container = brb.plate5.containers[0],
    #                              containerGeometryTableIndex = brb.bottle_20ml.containerGeometryTableIndex,
    #                              deckGeometryTableIndex = brb.deckgeom_50ul.index,
    #                              liquidClassTableIndex = '21',
    #                              lld = 1,
    #                              lldSearchPosition = '1700',
    #                              liquidSurface='1700',
    #                              tip_for_volume_check='50ul')
    # time.sleep(2)
    # print(zm.r.received_msg)
    #
    # time.sleep(2)
    # print(zm.r.received_msg)



