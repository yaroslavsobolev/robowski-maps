import logging
import breadboard as brb
module_logger = logging.getLogger('main.zeus')
import can
import signal
import time
from time import sleep
from colorama import init, Fore, Back, Style
from threading import Thread, Lock
import codecs
import json
import matplotlib.pyplot as plt
import matplotlib
# set the backend to use for matplotlib
matplotlib.use('TkAgg')
import pandas as pd
import sound
import matplotlib.colors as mcolors
import numpy as np
import csv
import sys
import pprint

DEBUG = 0
INFO = 1
WARNING = 1
ERROR = 1
# KICK_MASK = 0x0400
KICK_MASK = 0b10000000000
SENDER_ID_MASK = 0x03E0
RECEIVER_ID_MASK = 0x001F
EOM_MASK = 0b10000000

class Unbuffered(object):

    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)

# sys.stdout = Unbuffered(sys.stdout)
def printMSG(type, msg):
    ts = time.time()  # Timestamp
    if (type == 'info' and INFO == 1):
        pass
        # print('OK!')
        # print(
        #     Fore.WHITE + "(" + "{0:f}".format(ts) + ") INFO: " + msg + Style.RESET_ALL)

    elif (type == 'debug' and DEBUG == 1):
        print(Fore.MAGENTA +
              "(" + "{0:f}".format(ts) + ") DEBUG: " + msg + Style.RESET_ALL)

    elif (type == 'warning' and WARNING == 1):
        print(Fore.YELLOW + "(" + "{0:f}".format(ts) + ") WARNING: " +
              msg + Style.RESET_ALL)

    elif (type == 'error' and ERROR == 1):
        print(Fore.RED + "(" + "{0:f}".format(ts)
              + ") ERROR: " + msg + Style.RESET_ALL)
        # raise Exception(msg)

class DeckGeometry(object):

    def __init__(self, index=None, endTraversePosition=0,
                 beginningofTipPickingPosition=0, positionofTipDepositProcess=0):
        if index is None:
            raise ValueError(
                "Cannot initialize DeckGeometry instance with unspecified index.")
        self.index = index
        self.endTraversePosition = endTraversePosition
        self.beginningofTipPickingPosition = beginningofTipPickingPosition
        self.positionofTipDepositProcess = positionofTipDepositProcess


class LiquidClass(object):

    def __init__(self, id=0, index=None, liquidClassForFilterTips=0,
                 aspirationMode=0, aspirationFlowRate=0, overAspiratedVolume=0,
                 aspirationTransportVolume=0, blowoutAirVolume=0, aspirationSwapSpeed=0,
                 aspirationSettlingTime=0, lld=0, clldSensitivity=0, plldSensitivity=0,
                 adc=0, dispensingMode=0, dispensingFlowRate=0, stopFlowRate=0,
                 stopBackVolume=0, dispensingTransportVolume=0, acceleration=0,
                 dispensingSwapSpeed=0, dispensingSettlingTime=0, flowRateTransportVolume=0):
        if index is None:
            raise ValueError(
                "Cannot initialize LiquidClass instance with unspecified index.")
        self.id = id
        self.index = index
        self.liquidClassForFilterTips = liquidClassForFilterTips
        self.aspirationMode = aspirationMode
        self.aspirationFlowRate = aspirationFlowRate
        self.overAspiratedVolume = overAspiratedVolume
        self.aspirationTransportVolume = aspirationTransportVolume
        self.blowoutAirVolume = blowoutAirVolume
        self.aspirationSwapSpeed = aspirationSwapSpeed
        self.aspirationSettlingTime = aspirationSettlingTime
        self.lld = lld
        self.clldSensitivity = clldSensitivity
        self.plldSensitivity = plldSensitivity
        self.adc = adc
        self.dispensingMode = dispensingMode
        self.dispensingFlowRate = dispensingFlowRate
        self.stopFlowRate = stopFlowRate
        self.stopBackVolume = stopBackVolume
        self.dispensingTransportVolume = dispensingTransportVolume
        self.acceleration = acceleration
        self.dispensingSwapSpeed = dispensingSwapSpeed
        self.dispensingSettlingTime = dispensingSettlingTime
        self.flowRateTransportVolume = flowRateTransportVolume


class remoteFrameListener(can.Listener):

    def __init__(self, p):
        self.parent = p
        self.remote_flag = 0
        self.waiting_for_remote_flag = 0
        self.kick_flag = 0
        self.waiting_for_kick_flag = 0
        self.data_flag = 0
        self.received_msg = ""
        self.msg_complete_flag = 0
        self.last_transmitted = ""
        self.ltMutex = Lock()
        self.rfMutex = Lock()
        self.kfMutex = Lock()
        self.msg_last = can.Message()

    def on_message_received(self, msg):

        """This function is called whenever a message is received on the bus. It
        is the main handler for incoming messages. It is responsible for
        parsing the incoming message and taking appropriate action based on the
        message type. It also sets flags to indicate that a remote frame or a
        kick frame has been received. It also assembles the received message"""

        printMSG("debug", f'Received message: arbitration id = {msg.arbitration_id:X}')

        # REMOTE FRAME ACTION
        if (msg.is_remote_frame == True):
            #  if((msg.arbitration_id == 0x0000) or (msg.arbitration_id == 0x0020)):
            #  return
            if (self.getRemoteFlag() == 0):
                self.setRemoteFlag(1)
                printMSG("debug", "Received remote frame with ID = {}, DLC = {}.".format(
                    self.parseMsgID(msg.arbitration_id, "s"), msg.dlc))
                #  print(Fore.BLUE + "{}".format(msg) + Style.RESET_ALL)
            return

        # Prevent handler from responding multiple times to same message
        if (msg == self.msg_last):
            printMSG("debug", 'Same message as last. Skipping.')
            return
        self.msg_last = msg

        # KICK FRAME ACTION
        #  elif(msg.arbitration_id == 0x0420):
        if ((msg.arbitration_id & KICK_MASK) or (msg.arbitration_id ==0x0420)):
            if msg.arbitration_id == 0x0420:
                printMSG("debug", 'msg.arbitration_id == 0x0420')
            if (msg.arbitration_id == 0x0401):
                return
            if (self.getKickFlag() == 0):
                self.setKickFlag(1)
                printMSG("debug", "Received Kick.")
                #  print(Fore.BLUE + "{}".format(msg) + Style.RESET_ALL)
                if self.parent.auto_response:
                    self.parent.sendRemoteFrame(1)
            return

        else:
            # Ignore the messages that we've sent out
            if self.parseMsgID(msg.arbitration_id, "r") == 1:
                return

            # print(Fore.BLUE + "{}".format(msg) + Style.RESET_ALL)
            if (self.msg_complete_flag == 1):
                self.received_msg = ""
                self.msg_complete_flag = 0

            # self.received_msg += msg.data.replace(" ", "")[:-1]
            # self.received_msg += msg.data.replace(" ", "")
            self.received_msg += msg.data[:-1].decode('iso-8859-1')  # .replace(" ", "")

            if (self.msg_is_last(msg) == 0):
                if self.parent.auto_response:
                    self.parent.sendRemoteFrame(8)

            else:
                self.msg_complete_flag = 1
                # printMSG(
                #     "debug", "Assembled message {}".format(self.received_msg))
                #  if self.parent.auto_response:
                #  self.parent.sendRemoteFrame(8)

                # I'm really not sure whether the kick flag should be reset here
                self.setKickFlag(0)
                # printMSG("debug", 'Kick flag is set to 0.')

                ret = self.parent.parseErrors(self.received_msg) # This is where the error message is parsed

                # print(f'ret before = {ret}')

                if ret =='Movement error on pipetting drive.': # This is to ignore the pipetting drive error
                    ret = 'NONE'

                # print(f'ret after = {ret}')

                # TODO: Figure out why sometimes it's "NONE" (capitals) and sometimes "None" and unify the format
                if (str(ret) != "None") and (str(ret) != "NONE"):
                    printMSG("error", "{}".format(ret))
                    return ret

    def remote_received(self):
        if (self.remote_flag == 1):
            self.remote_flag = 0
            return 1
        else:
            return 0

    def setLastTransmitted(self, lt):
        self.ltMutex.acquire()
        lt = lt.replace("", "")
        self.last_transmitted = lt
        printMSG("warning", "Setting last transmitted to\
                {}".format(lt))
        self.ltMutex.release()

    def getLastTransmitted(self):
        self.ltMutex.acquire()
        ret = self.last_transmitted
        self.ltMutex.release()
        return ret

    def kick_received(self):
        if (self.kick_flag == 1):
            self.kick_flag = 0
            return 1
        else:
            return 0

    def data_received(self):
        if (self.data_flag == 1):
            self.data_flag = 0
            return 1
        else:
            return 0

    def msg_is_last(self, msg):
        size = len(msg.data)
        if (size > 0):
            printMSG(
             "debug", "control byte length = {}".format(len(msg.data)))
            control_byte = msg.data[(len(msg.data) - 1)]
            if ((control_byte & EOM_MASK) > 0):
                return 1

        return 0

    def parseMsgID(self, id, field):
        if field == "r":
            return id & RECEIVER_ID_MASK
        elif field == "s":
            return (id & SENDER_ID_MASK) >> 5
        else:
            return 0

    def setRemoteFlag(self, val):
        self.rfMutex.acquire()
        self.remote_flag = val
        self.rfMutex.release()

    def setKickFlag(self, val):
        self.kfMutex.acquire()
        self.kick_flag = val
        self.kfMutex.release()

    def getRemoteFlag(self):
        self.rfMutex.acquire()
        val = self.remote_flag
        self.rfMutex.release()
        return val

    def getKickFlag(self):
        self.kfMutex.acquire()
        val = self.kick_flag
        self.kfMutex.release()
        return val


def split_by_n(seq, n):
    """A generator to divide a sequence into chunks of n units."""
    while seq:
        yield seq[:n]
        seq = seq[n:]


class ZeusError(Exception):
    pass


class ZeusModule:
    CANBus = None
    transmission_retries = 10
    remote_timeout = 1
    errorTable = {
        "20": "No communication to EEPROM.",
        "30": "Undefined command.",
        "31": "Undefined parameter.",
        "32": "Parameter out of range.",
        "35": "Voltage outside the permitted range.",
        "36": "Emergency stop is active or was sent during action.",
        "38": "Empty liquid class.",
        "39": "Liquid class write protected.",
        "40": "Parallel processes not permitted.",
        "50": "Initialization failed.",
        "51": "Pipetting drive not initialized.",
        "52": "Movement error on pipetting drive.",
        "53": "Maximum volume of the tip reached.",
        "54": "Maximum volume in pipetting drive reached.",
        "55": "Volume check failed.",
        "56": "Conductivity check failed.",
        "57": "Filter check failed.",
        "60": "Initialization failed.",
        "61": "Z-drive is not initialized.",
        "62": "Movement error on the z-drive.",
        "63": "Container bottom search failed.",
        "64": "Z-position not possible.",
        "65": "Z-position not possible.",
        "66": "Z-position not possible.",
        "67": "Z-position not possible.",
        "68": "Z-position not possible.",
        "69": "Z-position not possible.",
        "70": "Liquid level not detected.",
        "71": "Not enough liquid present.",
        "72": "Auto calibration of the pressure sensor not possible.",
        "74": "Early liquid level detection.",
        "75": "No tip picked up or no tip present.",
        "76": "Tip already picked up.",
        "77": "Tip not discarded.",
        "80": "Clot detected during aspiration.",
        "81": "Empty tube detected during aspiration.",
        "82": "Foam detected during aspiration.",
        "83": "Clot detected during dispensing.",
        "84": "Foam detected during dispensing.",
        "85": "No communication to the digital potentiometer.",
    }

    ZeusTraversePosition = 880

    def __init__(self, id=None, tip_on_zeus='', init_module=True, auto_response=True):
        # colorama.init()
        init()
        self.logger = logging.getLogger("main.zeus.ZeusModule")
        self.id = id
        self.tip_on_zeus = tip_on_zeus
        self.auto_response = auto_response
        if id is None:
            raise ValueError(
                "Cannot initialize ZeusModule instance with unspecified id.")
        elif id not in range(1, 32):
            raise ValueError(
                "Cannot initialize ZeusModule instance with out of range id."
                " Valid id range is [1-31]")
        self.initCANBus()
        self.pos = 0
        self.minZPosition = 0
        self.maxZPosition = 2340
        self.r = remoteFrameListener(self)
        self.remoteFrameNotifier = can.Notifier(self.CANBus, [self.r])

        with open('config/liquid_class_table_para_ALL.json') as json_file:
            liquid_class_table_para = json.load(json_file)

        self.liquid_class_table_para = liquid_class_table_para

        self.logger.info(f"ZeusModule {self.id} is initializing...")

        if init_module:
            # self.getFirmwareVersion()
            self.initZDrive()
            printMSG("debug", 'sleeping before initDosingDrive')
            sleep(3)
            printMSG("debug", f'Kick flag = {self.r.getKickFlag()}')
            self.initDosingDrive()
            printMSG("debug", 'sleeping after initDosingDrive')
            sleep(3)

    def setAutoResponse(self, auto):
        self.auto_response = auto

    def cmdHeader(self, command):
        # return command + "id" + str(self.id).zfill(4)
        return command + "id0000"

    def assembleIdentifier(self, msg_type, master_id=0):
        identifier = 0
        identifier |= self.id
        if (master_id > 0):
            identifier |= (master_id << 5)
        if msg_type == 'kick':
            identifier |= 1 << 10
        # print(identifier)
        return identifier

    def sendRemoteFrame(self, dlc):
        # SEND REMOTE FRAME
        msg = can.Message(
            is_extended_id=False,
            is_remote_frame=True,
            arbitration_id=0x0020)
        msg.dlc = dlc
        printMSG(
            "debug", "ZeusModule {}: sending remote frame with dlc = {}...".format(self.id, msg.dlc))
        # print(Fore.GREEN + "{}".format(msg) + Style.RESET_ALL)
        try:
            self.CANBus.send(msg)
        except can.canError:
            print(
                Fore.RED + "ERROR: Remote frame not sent!" + Style.RESET_ALL)

    def waitForRemoteFrame(self):
        # WAIT FOR REMOTE RESPONSE
        # sleep(self.remote_timeout)
        s = time.time()
        c = time.time()
        # LOOP HERE UNTIL TIMEOUT EXPIRES
        while ((c - s) < self.remote_timeout):
            c = time.time()
            if (self.r.remote_received() == 1):
                printMSG("debug", f'waitForRemoteFrame: Received remote frame after {time.time() - s} s')
                # printMSG("debug", "ACK Received.")
                return 1
        return 0

    def sendKickFrame(self):
        # print('sending kick')
        n = 0
        msg = can.Message(
            is_extended_id=False,
            arbitration_id=self.assembleIdentifier('kick'), data=0)
        printMSG("debug",
                 "ZeusModule, {}: sending kick frame...".format(self.id))
        # print(Fore.GREEN + "{}".format(msg) + Style.RESET_ALL)
        while (n < self.transmission_retries):
            try:
                self.CANBus.send(msg)
            except can.canError:
                printMSG("error", "Kick not sent!")

            # WAIT FOR REMOTE RESPONSE
            if (self.waitForRemoteFrame() == 1):
                return
            if (n < self.transmission_retries):
                printMSG("warning", "Timeout waiting for kick response. Issuing retry {} of {}".format(
                    n + 1, self.transmission_retries))
            n += 1

        # Should not get here unless no response is received.
        printMSG(
            "error", "No remote frame received. Check connection to Zeus module.")
        exit(1)

    def waitForKickFrame(self):
        printMSG("debug", 'waitForKickFrame')
        # WAIT FOR REMOTE RESPONSE
        # sleep(self.remote_timeout)
        # sleep(0.001)
        s = time.time()
        c = time.time()
        # LOOP HERE UNTIL TIMEOUT EXPIRES
        while ((c - s) < self.remote_timeout):
            c = time.time()
            if (self.r.kick_received() == 1):
                printMSG("debug", "waitforkickframwe - ACK Received.")
                return 1

        return 0

    def sendDataObject(self, i, cmd_len, data):
        byte = 0
        printMSG(
            "debug", "ZeusModule {}: sending data frame {} of {}...".format(self.id, i + 1, cmd_len))
        printMSG('debug', "Outstring = {}".format(data))
        printMSG(
            # "debug", "data pre append = {}".format(data.encode('hex')))
            "debug", "data pre append (hex) = {}".format(codecs.encode(data.encode(), 'hex')))
        # Assemble the 8th (status) byte
        # Add EOM bit if this is the last frame of the message.
        if (i == (cmd_len - 1)):
            printMSG("debug", "appending EOM bit to byte")
            byte |= 1 << 7
        # Add number of data bytes
        byte |= (len(data) << 4)
        printMSG("debug", "num data bytes = {}".format(len(data)))
        byte |= ((i + 1) % 31)
        printMSG("debug", "frame counter = {}".format(((i + 1) % 31)))
        printMSG("debug", "control byte (binary) = {0:b}".format(byte))
        printMSG("debug", "control byte (hex) = {0:X}".format(byte))
        # PAD FRAME WITH ZEROES
        while (len(data) < 7):
            data += " "

        # # APPEND CONTROL BYTE
        # data += chr(byte)

        # APPEND CONTROL BYTE (encodings are devil's work)
        data = bytearray(data.encode('iso-8859-1'))
        # print(data)
        data.append(byte)
        # print(data)

        printMSG(
            # "debug", "data post append = {}".format(codecs.encode(data.encode(), 'hex')))
            "debug", "data post append = {}".format(data))
        printMSG(
            # "debug", "data pre append = {}".format(data.encode('hex')))
            "debug", "data post append (hex) = {0}".format(['{0:X}'.format(x) for x in data]))
        msg = can.Message(
            is_extended_id=False,
            arbitration_id=self.assembleIdentifier('data'),
            data=data)
        # self.r.setLastTransmitted(msg.data)
        try:
            self.CANBus.send(msg)
            # print(Fore.GREEN + "{}".format(msg) + Style.RESET_ALL)

        except can.CanError as err:
            printMSG("error", "Message not sent!")
            # raise can.CanError(err)

    def sendCommand(self, cmd):
        data = list(split_by_n(cmd, 7))
        # print(f'The split list sent is : {data}')
        cmd_len = len(data)
        # printMSG(
        #     "debug", "ZeusModule {}: sending packet {} in {} data frame(s)...".format(self.id, cmd, cmd_len))
        printMSG(
            "debug", "ZeusModule {}: sending packet {} in {} data frame(s)...".format(self.id, cmd, cmd_len))

        # Send kick frame and wait for remote response
        self.sendKickFrame()
        for i in range(0, cmd_len):
            #        n = 0
            # outstring = bytearray(data[i])
            outstring = data[i]
            for n in range(0, self.transmission_retries):
                # SEND DATA FRAME
                self.sendDataObject(i, cmd_len, outstring)

                # WAIT FOR REMOTE RESPONSE UNLESS FRAME IS LAST IN COMMAND
                # if (int(outstring[-1]) & EOM_MASK):
                #     return

                # Yaroslav's version
                printMSG("debug", f'outstring is {outstring}')
                if i == cmd_len - 1:
                    printMSG("debug", 'This frame is last in command.')
                    return
                if (self.waitForRemoteFrame() == 1):
                    printMSG("debug", 'Received remote frame')
                    break
                else:
                    printMSG("warning", "Timeout waiting for remote response. Issuing retry {} of {}".format(
                        n + 1, self.transmission_retries))
        print(f'cmd sent to zeus is : {cmd}')
        self.waitForKickFrame()

    def initCANBus(self):
        printMSG(
            "info", "ZeusModule {}: initializing CANBus...".format(self.id))
        can.rc['interface'] = 'kvaser'
        can.rc['channel'] = '0'
        can.rc['bitrate'] = 500000
        self.CANBus = can.interface.Bus()
        #  self.CANBus = can.interface.Bus(can_filters=[{"can_id": 0x01,
        #  "can_mask": 0xFF}])

    def initDosingDrive(self):
        cmd = self.cmdHeader('DI')
        printMSG(
            "info", "ZeusModule {}: initializing dosing drive...".format(self.id))
        self.sendCommand(cmd)

    def initZDrive(self):
        cmd = self.cmdHeader('ZI')
        printMSG(
            "info", "ZeusModule {}: initializing z-drive...".format(self.id))
        # self.pos = self.maxZPosition # This bug is commented out by Yankai.The following line is used
        self.pos = self.minZPosition
        self.sendCommand(cmd)

    def moveZDrive(self, pos, speed):
        cmd = self.cmdHeader('GZ')
        if (pos > self.maxZPosition) or (pos < self.minZPosition):
            raise ValueError(
                "ZeusModule {}: requested z-position out of range. "
                " Valid range for z-position is between {} and {}"
                .format(self.id, self.minZPosition, self.maxZPosition))
        if (speed == "slow"):
            speed = 0
        elif (speed == "fast"):
            speed = 1
        else:
            raise ValueError(
                "ZeusModule {}: invalid z-axis drive speed specified."
                " Accepted values for z-axis drive speed are \'slow\'"
                " and \'fast\'.")
        # print(
        #     "ZeusModule {}: moving z-drive from position {} to position {}."
        #     .format(self.id, self.pos, pos))
        cmd = cmd + 'gy' + str(pos).zfill(4) + 'gw' + str(speed)
        self.pos = pos
        self.sendCommand(cmd)

    def pickUpTip(self, tipTypeTableIndex, deckGeometryTableIndex):
        cmd = self.cmdHeader('GT')
        cmd = cmd + 'tt' + str(tipTypeTableIndex).zfill(
            2) + 'go' + str(deckGeometryTableIndex).zfill(2)

        self.sendCommand(cmd)

    def discardTip(self, deckGeometryTableIndex):
        cmd = self.cmdHeader('GU')
        cmd = cmd + 'go' + str(deckGeometryTableIndex).zfill(2)
        self.sendCommand(cmd)

    def sendString(self, string):
        # cmd = self.cmdHeader(string)
        self.sendCommand(string)

    def volumeCheck(self, containerGeometryTableIndex, deckGeometryTableIndex,
                    liquidClassTableIndex, lld, lldSearchPosition, liquidSurface):
        # print(f'containerGeometryTableIndex is {containerGeometryTableIndex}')
        # print(f'deckGeometryTableIndex is {deckGeometryTableIndex}')
        # print(f'liquidClassTableIndex is {liquidClassTableIndex}')
        # print(f'lld is {lld}')
        # print(f'lldSearchPosition is {lldSearchPosition}')
        # print(f'liquidSurface is {liquidSurface}')

        cmd = self.cmdHeader('GJ')
        cmd = cmd +\
                'ge' + str(containerGeometryTableIndex).zfill(2) + \
                'go' + str(deckGeometryTableIndex).zfill(2) + \
                'lq' + str(liquidClassTableIndex).zfill(2) + \
                'lb' + str(lld) + \
                'zp' + str(lldSearchPosition).zfill(4) + \
                'cf' + str(liquidSurface).zfill(4)
        self.sendCommand(cmd)

    def aspiration(self, aspirationVolume=0, containerGeometryTableIndex=0,
                   deckGeometryTableIndex=0, liquidClassTableIndex=0, qpm=0,
                   lld=0, lldSearchPosition=0, liquidSurface=0, mixVolume=0,
                   mixFlowRate=0, mixCycles=0):
        cmd = self.cmdHeader('GA')
        cmd = cmd + 'ai' + str(aspirationVolume).zfill(5) + \
              'ge' + str(containerGeometryTableIndex).zfill(2) + \
              'go' + str(deckGeometryTableIndex).zfill(2) + \
              'lq' + str(liquidClassTableIndex).zfill(2) + \
              'gq' + str(qpm) + \
              'lb' + str(lld) + \
              'zp' + str(lldSearchPosition).zfill(4) + \
              'cf' + str(liquidSurface).zfill(4) + \
              'ma' + str(mixVolume).zfill(5) + \
              'mb' + str(mixFlowRate).zfill(5) + \
              'dn' + str(mixCycles).zfill(2)
        # print(f"YANKAI_note: The command sent to Zeus is : {cmd}")
        self.sendCommand(cmd)

    def dispensing(self, dispensingVolume=0, containerGeometryTableIndex=0,
                   deckGeometryTableIndex=0, qpm=0, liquidClassTableIndex=0,
                   lld=0, lldSearchPosition=0, liquidSurface=0,
                   searchBottomMode=0, mixVolume=0, mixFlowRate=0, mixCycles=0):
        cmd = self.cmdHeader('GD')
        cmd = cmd + 'di' + str(dispensingVolume).zfill(5) + \
              'ge' + str(containerGeometryTableIndex).zfill(2) + \
              'go' + str(deckGeometryTableIndex).zfill(2) + \
              'gq' + str(qpm) + \
              'lq' + str(liquidClassTableIndex).zfill(2) + \
              'lb' + str(lld) + \
              'zp' + str(lldSearchPosition).zfill(4) + \
              'cf' + str(liquidSurface).zfill(4) + \
              'zm' + str(searchBottomMode) + \
              'ma' + str(mixVolume).zfill(5) + \
              'mb' + str(mixFlowRate).zfill(5) + \
              'dn' + str(mixCycles).zfill(2)
        # print(f"DEBUG: zeus.dispensing():: The command sent to Zeus is : {cmd}")
        # print(f"DEBUG: zeus.dispensing():: liquidSurface is : {liquidSurface}")
        self.sendCommand(cmd)

    def switchOff(self):
        cmd = self.cmdHeader('AV')
        self.sendCommand(cmd)

    def calculateContainerVolume(self, containerGeometryTableIndex=0,
                                 deckGeometryTableIndex=0,
                                 liquidClassTableIndex=0, lld=0,
                                 lldSearchPosition=0, liquidSurface=0):
        cmd = self.cmdHeader('GJ')
        cmd = cmd + 'ge' + str(containerGeometryTableIndex).zfill(2) + \
              'go' + str(deckGeometryTableIndex).zfill(2) + \
              'lq' + str(liquidClassTableIndex).zfill(2) + \
              'lb' + str(lld) + \
              'zp' + str(lldSearchPosition).zfill(4) + \
              'cf' + str(liquidSurface).zfill(3)
        self.sendCommand(cmd)

    def calculateContainerVolumeAfterPipetting(self):
        cmd = self.cmdHeader('GN')
        self.sendCommand(cmd)

    def getErrorCode(self):
        cmd = self.cmdHeader('RE')
        self.sendCommand(cmd)

    def getFirmwareVersion(self):
        # cmd = self.cmdHeader('RF')
        self.sendCommand("RF")
        # self.sendCommand(cmd)

    def getParameterValue(self, parameterName):
        if (len(parameterName) > 2):
            raise ValueError(
                "ZeusModule {}: Invalid parameter \'{}\' requested. "
                " Parameter format must be two lower-case letters."
                .format(parameterName))
        cmd = self.cmdHeader('RA')
        cmd = cmd + 'ra' + parameterName
        self.sendCommand(cmd)

    def getInstrumentInitializationStatus(self):
        cmd = self.cmdHeader('QW')
        self.sendCommand(cmd)

    def getNameofLastFaultyParameter(self):
        cmd = self.cmdHeader('VP')
        self.sendCommand(cmd)

    def getTipPresenceStatus(self):
        cmd = self.cmdHeader('RT')
        self.sendCommand(cmd)
        time.sleep(0.3)
        if 'rt1' in self.r.received_msg:
            return True
        else:
            return False

    def getTechnicalStatus(self):
        cmd = self.cmdHeader('QT')
        self.sendCommand(cmd)

    def getAbsoluteZPosition(self):
        cmd = self.cmdHeader('RZ')
        self.sendCommand(cmd)
        time.sleep(0.2)
        self.sendCommand(cmd)
        output = int(self.r.received_msg[-4:])
        return output

    def getCycleCounter(self):
        cmd = self.cmdHeader('RV')
        self.sendCommand(cmd)

    def getLifetimeCounter(self):
        cmd = self.cmdHeader('RY')
        self.sendCommand(cmd)

    def getInstrumentStatus(self):
        cmd = self.cmdHeader('RQ')
        self.sendCommand(cmd)
        pass

    def getLiquidClassRevision(self):
        cmd = self.cmdHeader('XB')
        self.sendCommand(cmd)
        pass

    def setEmergencyStop(self, state):
        if state in {1, 'on', 'ON', 'True', 'true'}:
            cmd = self.cmdHeader('AB')
        if state in {0, 'off', 'OFF', 'False', 'false'}:
            cmd = self.cmdHeader('AW')
        self.sendCommand(cmd)

    def switchDigitalOutput(self, out1State, out2State):
        cmd = self.cmdHeader('OU')
        cmd += 'ou'
        if out1State in {1, 'on', 'ON', 'True', 'true'}:
            cmd += str(1)
        if out2State in {0, 'off', 'OFF', 'False', 'false'}:
            cmd += str(0)

        cmd += 'oy'
        if out2State in {1, 'on', 'ON', 'True', 'true'}:
            cmd += str(1)
        if out2State in {0, 'off', 'OFF', 'False', 'false'}:
            cmd += str(0)

    def switchLEDStatus(self, blueState, redState):
        if (blueState not in set([0, 1])) or (redState not in set([0, 1])):
            raise ValueError(
                "ZeusModule {}: requested LED state out of range. "
                " Valid range for LED state is [0,1]".format(self.id))
        cmd = self.cmdHeader('SL')
        cmd += 'sl' + str(blueState) + \
               'sk' + str(redState)
        print("Switching status of blue LED to {} and red LED to {}"
              .format(blueState, redState))
        self.sendCommand(cmd)

    def testModeCommand(self, status):
        if (status not in set([0, 1])):
            raise ValueError(
                "ZeusModule {}: requested LED state out of range. "
                " Valid range for LED state is [0,1]".format(self.id))
        cmd = self.cmdHeader('AT')
        cmd += 'at' + str(status)
        self.sendCommand(cmd)

    def setDosingDriveInCleaningPosition(self):
        cmd = self.cmdHeader('GX')
        self.sendCommand(cmd)

    def setContainerGeometryParameters(self, container:object):
        # pprint(vars(container))
        print(f'Container name: {container}')
        cmd = self.cmdHeader('GC')
        cmd = cmd + 'ge' + str(container.containerGeometryTableIndex).zfill(2) + \
              'cb' + str(container.diameter).zfill(3) + \
              'bg' + str(container.bottomHeight).zfill(4) + \
              'gx' + str(container.bottomSection).zfill(5) + \
              'ce' + str(container.bottomPosition).zfill(4) + \
              'ie' + str(container.immersionDepth).zfill(4) + \
              'yq' + str(container.leavingHeight).zfill(4) + \
              'yr' + str(container.jetHeight).zfill(4) + \
              'ch' + \
              str(container.startOfHeightBottomSearch).zfill(4) + \
              'ci' + \
              str(container.dispenseHeightAfterBottomSearch).zfill(
                  4)
        print(f'cmd sent is : {cmd}')
        self.sendCommand(cmd)

    def getContainerGeometryParameters(self, index):
        if index is None:
            raise ValueError(
                "Please specify a valid container geometry table index.")
        cmd = self.cmdHeader('GB')
        cmd = cmd + 'ge' + str(index).zfill(2)
        self.sendCommand(cmd)
        # time.sleep(3)
        # paras = self.r.received_msg
        time.sleep(3)
        paras = self.r.received_msg

        return paras

    def setDeckGeometryParameters(self, deckGeometryParameters):
        cmd = self.cmdHeader('GO')
        cmd = cmd + 'go' + str(deckGeometryParameters.deckGeometryTableIndex).zfill(2) + \
              'te' + str(deckGeometryParameters.endTraversePosition).zfill(4) + \
              'tm' + \
              str(deckGeometryParameters.beginningofTipPickingPosition).zfill(4) + \
              'tr' + \
              str(deckGeometryParameters.positionofTipDepositProcess).zfill(
                  4)
        self.sendCommand(cmd)

    def getDeckGeometryParameters(self, index):
        if index is None:
            raise ValueError(
                "Please specify a valid deck geometry table index.")
        cmd = self.cmdHeader('GR')
        cmd = cmd + 'go' + str(index).zfill(2)
        self.sendCommmand(cmd)

    def setLiquidClassParameters(self, liquidClassParameters):
        if liquidClassParameters.index is None:
            raise ValueError(
                "Please specify a valid deck geometry table index.")
        # cmd = self.cmdHeader('GL')
        cmd = 'GL' + 'id' + str(liquidClassParameters.id).zfill(4) + \
              'lq' + str(liquidClassParameters.index).zfill(2) + \
              'uu' + str(liquidClassParameters.liquidClassForFilterTips) + \
              ' ' + str(liquidClassParameters.aspirationMode) + \
              ' ' + str(liquidClassParameters.aspirationFlowRate).zfill(5) + \
              ' ' + str(liquidClassParameters.overAspiratedVolume).zfill(4) + \
              ' ' + str(liquidClassParameters.aspirationTransportVolume).zfill(5) + \
              ' ' + str(liquidClassParameters.blowoutAirVolume).zfill(5) + \
              ' ' + str(liquidClassParameters.aspirationSwapSpeed).zfill(4) + \
              ' ' + str(liquidClassParameters.aspirationSettlingTime).zfill(3) + \
              ' ' + str(liquidClassParameters.lld) + \
              ' ' + str(liquidClassParameters.clldSensitivity) + \
              ' ' + str(liquidClassParameters.plldSensitivity) + \
              ' ' + str(liquidClassParameters.adc) + \
              ' ' + str(liquidClassParameters.dispensingMode) + \
              ' ' + str(liquidClassParameters.dispensingFlowRate).zfill(5) + \
              ' ' + str(liquidClassParameters.stopFlowRate).zfill(5) + \
              ' ' + str(liquidClassParameters.stopBackVolume).zfill(3) + \
              ' ' + str(liquidClassParameters.dispensingTransportVolume).zfill(5) + \
              ' ' + str(liquidClassParameters.acceleration).zfill(3) + \
              ' ' + str(liquidClassParameters.dispensingSwapSpeed).zfill(4) + \
              ' ' + str(liquidClassParameters.dispensingSettlingTime).zfill(3) + \
              ' ' + str(liquidClassParameters.flowRateTransportVolume).zfill(5)
        print(f'cmd send by GL: {cmd}')
        self.sendCommand(cmd)

    def getLiquidClassParameters(self, id, index):
        cmd = self.cmdHeader('GM')
        cmd = cmd + 'lq' + str(index).zfill(2)  # 'lq' was revised from 'iq'. iq is a typo. Yankai_20230106
        print(f'cmd send is : {cmd}')
        self.sendCommand(cmd)

    def firmwareUpdate(self, filename):
        pass

    def parseErrors(self, errorString):
        if len(errorString.replace(" ", "")) == 0:
            return "NONE"
        if errorString == "":
            return "NONE"

        printMSG("info", "Message from Zeus: {}".format(errorString))

        #  else:
        # print(Fore.MAGENTA + "DEBUG: Received error string '{}' with
        # length: {}".format(errorString, len(errorString.replace(" ",
        # ""))) + Style.RESET_ALL)

        cmd = str(errorString[:2])
        # print("cmd = {}".format(cmd))
        eidx = errorString.find("er")
        if (eidx == -1):
            return "NONE"

        ec = str(errorString[(eidx + 2): (eidx + 4)])
        if ec == '00':
            # NO ERROR
            return
        #  print(Fore.MAGENTA + "DEBUG: ec = {}".format(ec))

        defaultError = "Unknown error code returned."

        if cmd == 'DI':
            if ec in {'00', '30', '35', '36', '40', '50', '52'}:
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'ZI':
            if ec in {'00', '30', '35', '36', '40', '60', '62'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GZ':
            if ec in {'00', '31', '32', '35', '36', '40', '61', '62', '64'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GT':
            if ec in {'00', '31', '32', '35', '36', '40', '51', '52', '61', '62', '65', '75', '76'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GU':
            if ec in {'00', '30', '31', '32', '35', '36', '40', '51', '52', '61', '62', '65', '69', '75', '77'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GA':
            if ec in {'00', '30', '31', '32', '35', '36', '38', '40', '51', '52', '53', '54', '55', '56', '57', '61',
                      '62', '65', '66', '67', '68', '70', '71', '72', '74', '75', '80', '81', '82', '85'}:
                self.logger.error(self.errorTable[ec])

                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GD':
            if ec in {'00', '30', '31', '32', '35', '36', '38', '40', '51', '52', '54', '55', '57', '61', '62', '63',
                      '65', '66', '67', '68', '70', '72', '74', '75', '83', '84', '85'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError
        elif cmd == 'GD':
            if ec in {'00', '30', '31', '32', '35', '36', '38', '40', '51', '52', '54', '55', '57', '61', '62', '63',
                      '65', '66', '67', '68', '70', '72', '74', '75', '83', '84', '85'}:
                self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GJ':
            if ec in {'00', '30', '31', '32', '35', '36', '38', '40', '51', '52', '56', '57', '61', '62', '65', '66',
                      '67', '68', '70', '72', '74', '85'}:
                self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]

            else:
                return defaultError

        elif cmd == 'AB':
            if ec in {'00', '30'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError
        elif cmd == 'AW':
            if ec in {'00', '30'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError
        elif cmd == 'XA':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GK':
            if ec in {'00', '30', '31', '32', '35', '51', '52', '54'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GC':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GO':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GB':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GR':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GL':
            if ec in {'00', '20', '30', '31', '32', '39'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GM':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GQ':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GS':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GV':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GW':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GG':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GE':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GH':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError

        elif cmd == 'GI':
            if ec in {'00', '20', '30', '31', '32'}:
                # self.logger.error(self.errorTable[ec])
                return self.errorTable[ec]
            else:
                return defaultError
        else:
            return "Error code returned '{}' corresponds to unknown command.".format(errorString)

    def wait_until_zeus_reaches_traverse_height(self, n_retries=70):
        traverse_height = self.ZeusTraversePosition
        # time.sleep(0.5)
        for i in range(n_retries):
            # print(f'Waiting for Zeus to get back to traverse height: attempt {i}')
            self.getAbsoluteZPosition()

            time.sleep(0.15)

            idx = self.r.received_msg.find("gy")

            if idx == -1:
                # this means that there is an error. Retry
                self.parseErrors(self.r.received_msg)
                time.sleep(0.5)
                continue
            else:
                position = int(self.r.received_msg[idx + 2:])
                # print(f'Current position (true): {position}')
            if position <= traverse_height:
                # print('Traverse height is reached.')
                return True
        print(f'Traverse height was not reached after {n_retries} retries. This is dangerous, so we do emergency stop')
        raise Exception
        return False

    def zeus_had_error(self, errorString):
        cmd = str(errorString[:2])
        eidx = errorString.find("er")

        if (eidx == -1):
            return False

        ec = str(errorString[(eidx + 2): (eidx + 4)])

        if ec == '00':
            return False

        return True

    def zeus_error_code(self, errorString):
        cmd = str(errorString[:2])
        eidx = errorString.find("er")
        if (eidx == -1):
            return False
        return str(errorString[(eidx + 2): (eidx + 4)])

    def wait_until_zeus_responds_with_string(self, search_pattern, n_retries=200):
        # module_logger.info(f'Waiting for search_pattern :  {search_pattern}')
        for i in range(n_retries):
            # print(f'Waiting for Zeus to respond: attempt {i}')
            # zm.getAbsoluteZPosition()
            # time.sleep(0.6)
            if self.zeus_had_error(self.r.received_msg) and 'er52' not in self.r.received_msg:
                # print('Zeus responded with error message. Aborting all operations.')
                self.logger.error('Zeus responded with error message. Aborting all operations.')
                self.logger.error(f'error message: {self.r.received_msg}')
                time.sleep(0.2)
                self.move_z(ZeusModule.ZeusTraversePosition)
                sound.beep_for_error()
                raise ZeusError
            elif self.zeus_had_error(self.r.received_msg) and 'er52' in self.r.received_msg: # this error will be ignored by Zeus.
                module_logger.warning(f'Warning: Zeus has Movement error on pipetting drive (code: 52). This is ignored.')
                time.sleep(0.5)
                return True

            idx = self.r.received_msg.find(search_pattern)

            if idx == -1:
                # this means that there is no pattern. Retry
                # zm.parseErrors(zm.r.received_msg)
                time.sleep(0.1)
                continue
            else:
                # print(f'Competion response received after {i} attempts.')
                return True

        print(f'Response not received after {n_retries} retries. This is dangerous, so we do emergency stop')
        self.logger.error(f'Response not received after {n_retries} retries. This is dangerous, so we do emergency stop')
        raise Exception
        return False

    def move_z(self, z, raise_exception=True):
        self.moveZDrive(z, 'fast')

        if raise_exception:
            self.wait_until_zeus_responds_with_string('GZid')

    def zeus_is_at_traverese_height(self):
        if self.pos <= self.ZeusTraversePosition:
            return True
        else:
            print(f'ERROR: ZEUS was not in traverse height before motion, but instead at {self.pos}')
            return False

    def move_zeus_to_traverse_height(self):
        if self.zeus_is_at_traverese_height():
            return
        else:
            self.move_z(self.ZeusTraversePosition)

    def check_last_faulty_para(self):
        self.sendCommand('VPid0001')

    """
        The following methods are for reading and setting liquid class parametes to Zeus.

        Working steps:
            1 extract all the parameters for the built-in liquid classes.
            2 store all the extracted parameters to a dictionary. step 1 and 2 need to be done only once.
                ## extract_all_built_in_liquid_class_parameters_to_a_dict()
            3 name your new liquid class from index 21, because the first 20 index are read only. 
              Copy one of the liquid class parameters to your new liquid class.
                ## copy_para_from_to(index_from, index_to)            
            4 modify your new liquid class by setting new parameters
                a. set the calibration paras to perfect pipetting mode, e.g., 
                        00100 00100 00200 00200 00500 00500 01000 01000 02000 02000 05000 05000 07500 07500 10000 10000
                b, measure calibration curve and set the calibration paras to the measured values, e.g., 
                        00100 00109 00200 00214 00500 00538 01000 01069 02000 02095 05000 05187 07500 07763 10000 10340
            5 update your local dictionary, i.e., Json file.
                ## update_liquid_dict()
            6 set your new liquid class to zeus
                ## set_liquid_class_to_zeus( liquid_index)
            7 check your liquid class parameters from zeus
                ## request_parameters_from_zeus(liquid_index):

        Notes
            1 do not forget the space between parameters when sending it the zeus
            2 there is one glitch with 'GG' and 'GH' commands. See inside function: set_liquid_class_to_zeus( liquid_index)

        Yankai Jia 2023/01/23
        """

    def import_from_json(self):
        return self.liquid_class_table_para

    def extract_liquid_class_parameter(self, liquid_index, id='0001'):
        cmd = 'GMid' + id + 'lq' + str(liquid_index).zfill(2)
        print(f'cmd send is : {cmd}')
        self.sendCommand(cmd)
        time.sleep(1)  # This delay is IMPORTANT. Without this delay, the msg will return the previous data.
        # This value should be larger than 0.1s.
        msg_received_from_Zeus = self.r.received_msg
        print(f'msg_received_from_Zeus len is : {len(msg_received_from_Zeus)}')
        return msg_received_from_Zeus

    def extract_calibration_aspiration(self, liquid_index, id='0001'):

        cmd = 'GEid' + id + 'gg' + str(liquid_index).zfill(2)
        print(f'cmd send is : {cmd}')
        self.sendCommand(cmd)
        time.sleep(0.5)  # This delay is IMPORTANT. Without this delay, the msg will return the previous data.
        # This value should be larger than 0.1s.
        msg_received_from_Zeus = self.r.received_msg
        print(f'msg_received_from_Zeus for calibration_aspiration is : {msg_received_from_Zeus}')
        return msg_received_from_Zeus

    def extract_calibration_dispensing(self, liquid_index, id='0001'):
        cmd = 'GIid' + id + 'gh' + str(liquid_index).zfill(2)
        print(f'cmd send is : {cmd}')
        self.sendCommand(cmd)
        time.sleep(0.5)
        msg_received_from_Zeus = self.r.received_msg
        # print(f'msg_received_from_Zeus for calibration_dispensing is : {msg_received_from_Zeus}')
        return msg_received_from_Zeus

    def extract_qpm_aspiration(self, liquid_index, id='0001'):
        cmd = 'GSid' + id + 'gv' + str(liquid_index).zfill(2)
        print(f'cmd send is : {cmd}')
        self.sendCommand(cmd)
        time.sleep(0.5)
        msg_received_from_Zeus = self.r.received_msg
        print(f'msg_received_from_Zeus for qpm_aspiration is : {msg_received_from_Zeus}')
        return msg_received_from_Zeus

    def extract_qpm_dispensing(self, liquid_index, id='0001'):
        cmd = 'GWid' + id + 'gp' + str(liquid_index).zfill(2)
        print(f'cmd send is : {cmd}')
        self.sendCommand(cmd)
        time.sleep(1)
        msg_received_from_Zeus = self.r.received_msg
        print(f'msg_received_from_Zeus for qpm_dispensing is : {msg_received_from_Zeus}')
        return msg_received_from_Zeus

    def fill_one_liquid_class_parameter(self, liquid_index, id='0001'):
        msg = self.extract_liquid_class_parameter(id=id, liquid_index=liquid_index)
        para_container = ''.join(i for i in msg if i.isdigit())
        # print(para_container)
        # global liquid_class_table_para
        var_dict = {'id': 4,
                    'index': 2,
                    'liquidClassForFilterTips': 1,
                    'aspirationMode': 1,
                    'aspirationFlowRate': 5,
                    'overAspiratedVolume': 4,
                    'aspirationTransportVolume': 5,
                    'blowoutAirVolume': 5,
                    'aspirationSwapSpeed': 4,
                    'aspirationSettlingTime': 3,
                    'lld': 1,
                    'clldSensitivity': 1,
                    'plldSensitivity': 1,
                    'adc': 1,
                    'dispensingMode': 1,
                    'dispensingFlowRate': 5,
                    'stopFlowRate': 5,
                    'stopBackVolume': 3,
                    'dispensingTransportVolume': 5,
                    'acceleration': 3,
                    'dispensingSwapSpeed': 4,
                    'dispensingSettlingTime': 3,
                    'flowRateTransportVolume': 5}  # store variables and its len
        # if not sum(var_dict.values()) +  ==  len (para_container):
        #     print("Msg string length does not match the needed length!")
        #     return
        n = 0
        for i in var_dict:
            # print(f'liquid index is : {liquid_index}')
            self.liquid_class_table_para['liquid_class_para'][liquid_index][i] = int(para_container[n:n + var_dict[i]])
            n += var_dict[i]
        return self.liquid_class_table_para

    def extract_all_built_in_liquid_class_parameters_to_a_dict(self):
        # global liquid_class_table_para
        for i in range(40):
            liquid_index = str(i).zfill(2)

            self.liquid_class_table_para['liquid_class_para'][liquid_index] = {}
            self.fill_one_liquid_class_parameter(liquid_index, id='0001')

            self.liquid_class_table_para['calibration']['aspiration'][liquid_index] = {}
            string1 = self.extract_calibration_aspiration(liquid_index, id='0001')
            self.liquid_class_table_para['calibration']['aspiration'][liquid_index] = string1[:-3]

            self.liquid_class_table_para['calibration']['dispensing'][liquid_index] = {}
            string2 = self.extract_calibration_dispensing(liquid_index, id='0001')
            self.liquid_class_table_para['calibration']['dispensing'][liquid_index] = string2[:-3]

            self.liquid_class_table_para['qpm']['aspiration'][liquid_index] = {}
            string3 = self.extract_qpm_aspiration(liquid_index, id='0001')
            self.liquid_class_table_para['qpm']['aspiration'][liquid_index] = string3[:-4]

            self.liquid_class_table_para['qpm']['dispensing'][liquid_index] = {}
            string4 = self.extract_qpm_dispensing(liquid_index, id='0001')
            self.liquid_class_table_para['qpm']['dispensing'][liquid_index] = string4

    # extract_all_built_in_liquid_class_parameters_to_a_dict() # Do this only when needed

    def copy_para_from_to(self, index_from, index_to):
        self.liquid_class_table_para['liquid_class_para'][str(index_to).zfill(2)] = \
            self.liquid_class_table_para['liquid_class_para'][str(index_from).zfill(2)]
        self.liquid_class_table_para['liquid_class_para'][str(index_to).zfill(2)][
            'index'] = index_to  # Update new liquid index

        string1 = self.liquid_class_table_para['calibration']['aspiration'][str(index_from).zfill(2)]
        new_string1 = 'GGid0001gg' + str(index_to).zfill(2) + string1[
                                                              12:]  # Update new liquid index and change GE (reequest) to GG (set)
        self.liquid_class_table_para['calibration']['aspiration'][str(index_to).zfill(2)] = new_string1

        string2 = self.liquid_class_table_para['calibration']['dispensing'][str(index_from).zfill(2)]
        new_string2 = 'GHid0001gh' + str(index_to).zfill(2) + string2[
                                                              12:]  # Update new liquid index and change GI(reequest) to GH (set)
        self.liquid_class_table_para['calibration']['dispensing'][str(index_to).zfill(2)] = new_string2

        string3 = self.liquid_class_table_para['qpm']['aspiration'][str(index_from).zfill(2)]
        new_string3 = 'GQid0001gv' + str(index_to).zfill(2) + string3[
                                                              12:]  # Update new liquid index and change GS (reequest) to GQ (set)
        self.liquid_class_table_para['qpm']['aspiration'][str(index_to).zfill(2)] = new_string3

        string4 = self.liquid_class_table_para['qpm']['dispensing'][str(index_from).zfill(2)]
        new_string4 = 'GVid0001gp' + str(index_to).zfill(2) + string4[
                                                              12:]  # Update new liquid index and change GW (reequest) to Gv (set)
        self.liquid_class_table_para['qpm']['dispensing'][str(index_to).zfill(2)] = new_string4

        self.update_liquid_dict_from_classvar_to_json()

    def update_liquid_dict_from_classvar_to_json(self):
        # update json file
        with open('config/liquid_class_table_para_ALL.json', 'w', encoding='utf-8') as f:
            json.dump(self.liquid_class_table_para, f, ensure_ascii=False, indent=4)

    def update_liquid_dict_from_json_to_classvar(self):
        # update json file
        with open('config/liquid_class_table_para_ALL.json', 'r', encoding='utf-8') as f:
            self.liquid_class_table_para = json.load(f)

    def set_liquid_class_to_zeus(self, liquid_index):
        # write liquid class parameters
        para1 = self.liquid_class_table_para['liquid_class_para'][str(liquid_index).zfill(2)]
        print(f'liquid class para: {para1}')
        lc_param = LiquidClass(**para1)
        self.setLiquidClassParameters(lc_param)
        time.sleep(5)
        print('Liquid class parameters set')

        ## write calibration curve
        # aspiration
        para2 = self.liquid_class_table_para['calibration']['aspiration'][str(liquid_index).zfill(2)]
        print(f'calibration_asp: {para2}')
        """
        There is a firmware malfunction here. Instead of send the string: GGid0001gg21ck + parameters. You should remove
        gg21 from the string and send this:GGid0001ck + parameters. But before send this, you should do this:
        zm.sendCommand('GHid0001gh21') and zm.sendCommand('RAid0000ragh'). Yaroslav figured this out. There was a lot of frustration
        before this was figured out.
        """
        self.sendCommand('GGid0001gg' + str(liquid_index).zfill(2))
        time.sleep(5)
        self.sendCommand('RAid0000ragg')
        time.sleep(5)
        para2_new = para2[:8] + para2[12:]
        self.sendCommand(para2_new)
        print(f'pra2_new: {para2_new}')
        time.sleep(5)
        print('Calibration for asp curve set')

        # dispensing
        para3 = self.liquid_class_table_para['calibration']['dispensing'][str(liquid_index).zfill(2)]
        print(f'calibration_disp: {para3}')
        self.sendCommand('GHid0001gh' + str(liquid_index).zfill(2))
        # time.sleep(5)
        self.sendCommand('RAid0000ragh')
        time.sleep(5)
        para3_new = para3[:8] + para3[12:]
        self.sendCommand(para3_new)
        time.sleep(5)
        print('Calibration for disp curve set')

        #
        # set_liquid_class_to_zeus(liquid_index = 21)
        # set_liquid_class_to_zeus(liquid_index = 22)
        ## write qpm
        # aspiration
        para4 = self.liquid_class_table_para['qpm']['aspiration'][str(liquid_index).zfill(2)]
        print(f'qpm for asp curve: {para4}')
        self.sendCommand(para4)
        time.sleep(5)
        print('qpm for asp curve set')

        # dispensing
        para5 = self.liquid_class_table_para['qpm']['dispensing'][str(liquid_index).zfill(2)]
        print(f'qpm for disp curve: {para5}')
        self.sendCommand(para5)
        time.sleep(5)
        print('qpm for disp curve set')


    def request_parameters_from_zeus(self, liquid_index):

        # liquid parameters
        self.sendCommand('GMid0001lq' + str(liquid_index).zfill(2))
        time.sleep(0.5)
        self.logger.info(f"liquid parameters: {self.r.received_msg}")
        print(f"liquid parameters: {self.r.received_msg}")

        # calibrations
        self.sendCommand('GEid0001gg' + str(liquid_index).zfill(2))
        time.sleep(0.5)
        self.logger.info(f"calibration_asp {self.r.received_msg}")
        print(f"calibration_asp {self.r.received_msg}")

        self.sendCommand('GIid0001gh' + str(liquid_index).zfill(2))
        time.sleep(0.5)
        self.logger.info(f"calibration_disp {self.r.received_msg}")
        print(f"calibration_disp {self.r.received_msg}")

        # qpm
        self.sendCommand('GSid0001gv' + str(liquid_index).zfill(2))
        time.sleep(0.5)
        self.logger.info(f"qpm_asp {self.r.received_msg}")
        time.sleep(0.5)
        print(f"qpm_disp {self.r.received_msg}")

        self.sendCommand('GWid0001gp' + str(liquid_index).zfill(2))
        time.sleep(0.5)
        self.logger.info(f"qpm_disp {self.r.received_msg}")
        print(f"qpm_disp {self.r.received_msg}")

    # set_liquid_class_to_zeus( liquid_index =23 )

    def get_pressure_data(self):
        import csv
        self.sendCommand('QHid0000')
        msg_here = self.r.received_msg
        time.sleep(1)
        msg_here = self.r.received_msg
        num_of_point = int(msg_here.split(' ')[1])
        start_id = 5000 - num_of_point

        data = []

        num_of_extract_needed =  num_of_point // 50 +1  # +1 is for the remainder of points
        remainder = num_of_point % 50

        for i in range(num_of_extract_needed):

            if i < num_of_extract_needed-1:
                self.sendCommand('QIid0000li'+str(start_id)+'ln50')
            else:
                self.sendCommand('QIid0000li'+str(start_id)+'ln'+str(remainder))

            time.sleep(0.2)
            str_here = self.r.received_msg
            time.sleep(0.2)
            str_here  = self.r.received_msg

            str_list = str_here.split(' ')
            str_list[0] = str_list[0][-4:]

            print(str_list)

            print(str_list)

            data_here = [int(i) for i in str_list if i.isdigit()]

            for data_point in data_here:
                data.append(data_point)

            start_id += 50

        print(data)

        ## save the pressure data to csv
        try:
            with open(f'pressure_data/pressure_data_{time.time()}.csv', 'a', newline='\n') as myfile:
                wr = csv.writer(myfile)
                wr.writerow(data)
        except: print('Error in saving pressure data to csv')

        return data

    def plot_pressure_curve(self):
        # read data from csv to dataframe
        df = pd.read_csv('pressure_data/pressure_data.csv', header= None, on_bad_lines='warn')
        df = df.T
        # pressure_data = np.genfromtxt('pressure_data/pressure_data.csv', delimiter=',')
        plt.plot(df.iloc[:,-1])
        plt.show()

    def get_pressure_data(self):
        zm.sendCommand('QHid0000')
        msg_here = zm.r.received_msg
        time.sleep(1)
        msg_here = self.r.received_msg
        num_of_point = int(msg_here.split(' ')[1])
        start_id = 5000 - num_of_point

        data = []

        num_of_extract_needed = num_of_point // 50 + 1  # +1 is for the remainder of points
        remainder = num_of_point % 50

        for i in range(num_of_extract_needed):

            if i < num_of_extract_needed - 1:
                self.sendCommand('QIid0000li' + str(start_id) + 'ln50')
            else:
                self.sendCommand('QIid0000li' + str(start_id) + 'ln' + str(remainder))

            time.sleep(0.2)
            str_here = self.r.received_msg
            time.sleep(0.2)
            str_here = self.r.received_msg

            str_list = str_here.split(' ')
            str_list[0] = str_list[0][-4:]

            print(str_list)

            print(str_list)

            data_here = [int(i) for i in str_list if i.isdigit()]

            for data_point in data_here:
                data.append(data_point)

            start_id += 50

        print(data)

        ## save the pressure data to csv

        with open('pressure_data/pressure_data.csv', 'a', newline='\n') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(data)


if __name__ == '__main__':

    print('This is main of zeus.py')

    # load container parameters
    zm = ZeusModule(id=1)


