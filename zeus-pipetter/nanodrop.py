##  connect to arduino and send commands "1" and "0"
import serial, time
import asyncio

# make a class for nanodrop
class Nanodrop:
    def __init__(self, nd_id):

        self.lid_closed = True
        self.nd_id = nd_id

        if self.nd_id == 'nd_2000':
            self.com_id = 'COM4'
            self.command_open_lid = b'1'
            self.command_close_lid = b'0'
            self.command_open_liquid = b'40'
            self.command_close_liquid = b'41'
            self.command_open_air = b'30'
            self.command_close_air = b'31'
            self.command_open_vacuum = b'20'
            self.command_close_vacuum = b'21'
            self.flush_time = 20
            self.dry_time = 12

        elif self.nd_id == 'nd_2000c':
            self.com_id = 'COM6'
            self.command_open_lid = b'1'
            self.command_close_lid = b'0'
            self.command_open_liquid = b'30'
            self.command_close_liquid = b'31'
            self.command_open_air = b'40'
            self.command_close_air = b'41'
            self.command_open_vacuum = b'20'
            self.command_close_vacuum = b'21'
            self.flush_time = 12
            self.dry_time = 18

        else:
            raise ValueError('The Nanodrop id is wrong!')

        try:
            self.serial = serial.Serial(self.com_id, 9600, timeout=1)
        except:
            print("Arduino not connected")

    def open_lid(self):
        if self.lid_closed:
            self.close_vacumm()
            time.sleep(0.2)
            # self.open_air()
            # time.sleep(0.2)
            self.serial.write(self.command_open_lid)
            self.lid_closed = False
        else:
            print("lid already open")

    def close_lid(self):
        if not self.lid_closed:
            # self.open_air()
            # time.sleep(0.2)
            self.serial.write(self.command_close_lid)
            time.sleep(0.2)
            # self.close_air()
            # time.sleep(0.2)
            self.lid_closed = True
        else:
            print("lid already closed")

    def open_liquid(self):
        self.serial.write(self.command_open_liquid)
    def close_liquid(self):
        self.serial.write(self.command_close_liquid)
    def open_air(self):
        self.close_air()
        time.sleep(0.2)
        self.serial.write(self.command_open_air)
    def close_air(self):
        self.serial.write(self.command_close_air)

    def open_vacumm(self):
        self.serial.write(self.command_open_vacuum)

    def close_vacumm(self):
        self.serial.write(self.command_close_vacuum )

    async def flush_pedestal(self):
        time_stamp = time.time()
        self.open_vacumm()
        time.sleep(0.2)
        self.close_air()
        time.sleep(0.2)
        self.open_liquid()
        print("flushing...")
        await asyncio.sleep(self.flush_time)
        self.close_liquid()
        time.sleep(0.1)
        self.close_vacumm()
        time.sleep(0.1)
    def flush_pedestal_not_async(self):
        self.open_vacumm()
        time.sleep(0.2)
        self.close_air()
        time.sleep(0.2)
        self.open_liquid()
        print("flushing...")
        time.sleep(self.flush_time)
        self.close_liquid()
        time.sleep(0.1)
        self.close_vacumm()
        time.sleep(0.1)

    async def dry_pedestal(self):
        time_stamp = time.time()
        self.open_vacumm()
        time.sleep(0.2)
        self.close_liquid()
        time.sleep(0.2)
        self.open_air()
        print('drying pedestal...')
        await asyncio.sleep(self.dry_time)
        self.close_vacumm()
        time.sleep(0.3)
        self.close_air()

    def dry_pedestal_not_async(self):
        self.open_vacumm()
        time.sleep(0.5)
        self.close_liquid()
        time.sleep(0.5)
        self.open_air()
        print('drying pedestal...')
        time.sleep(self.dry_time)
        self.close_vacumm()
        time.sleep(0.5)
        self.close_air()

    async def flush_then_dry_pedestal(self):

        self.close_air()
        time.sleep(0.1)
        self.open_vacumm()
        time.sleep(0.1)
        self.open_liquid()
        time.sleep(self.dry_time)
        self.close_liquid()
        time.sleep(0.1)
        self.close_vacumm()
        time.sleep(0.1)

        self.open_vacumm()
        time.sleep(0.1)
        self.close_liquid()
        time.sleep(0.1)
        self.open_air()
        await asyncio.sleep(self.dry_time)
        self.close_vacumm()
        time.sleep(0.1)
        self.open_air()
        time.sleep(0.1)

    def close_all(self):
        self.close_lid()
        time.sleep(0.5)
        self.close_vacumm()
        time.sleep(0.5)
        self.close_liquid()
        time.sleep(0.5)
        self.close_air()
        time.sleep(0.5)

    def close_serial(self):
        self.serial.close()

    def test50(self):
        for i in range(50):
            nd.open_lid()
            time.sleep(3)
            nd.close_lid()
            time.sleep(3)

if __name__ ==  '__main__':
    nd = Nanodrop('nd_2000c')
    print(1)