import numpy as np
import time
from Fluigent.SDK import fgt_init, fgt_close
from Fluigent.SDK import fgt_set_pressure, fgt_get_pressure, fgt_get_pressureRange
from pyaxidraw import axidraw


def generate_well_coordinates(Nwells, topleft, topright, bottomleft, bottomright):
    '''generate coordinates for all wells of a well plate from coordinates of corner wells.'''
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

xy_offset = (2.5, -5)

# #Aluminum
# vials = generate_well_coordinates(Nwells=(6, 9),
#                                  topleft =   (xy_offset[0] + 28.5  ,   xy_offset[1] + 10),
#                                  topright=   (xy_offset[0] + 132.8 ,   xy_offset[1] + 8.5),
#                                  bottomleft= (xy_offset[0] + 29.8,   xy_offset[1] + 74.5),
#                                  bottomright=(xy_offset[0] + 133.0 , xy_offset[1] + 73.5))

vials = generate_well_coordinates(Nwells=(6, 9),
                                 topleft =   (xy_offset[0] + 26.4  ,   xy_offset[1] + 7.4),
                                 topright=   (xy_offset[0] + 130.7 ,   xy_offset[1] + 5.7),
                                 bottomleft= (xy_offset[0] + 27.8,   xy_offset[1] + 71.8),
                                 bottomright=(xy_offset[0] + 131.5 , xy_offset[1] + 70.4))


ad = axidraw.AxiDraw()
ad.interactive()
ad.options.pen_pos_up = 100
ad.options.pen_pos_down = 0
ad.options.pen_rate_raise = 75
ad.options.pen_rate_lower = 12
ad.options.units = 2
ad.options.pen_pos_up = 100
ad.options.pen_pos_down = 0
print(ad.connect())
ad.update()
ad.penup()

fgt_init()
pressureMeasurement = fgt_get_pressure(0)
print('Current pressure: {}'.format(pressureMeasurement))


def start_flow():
    fgt_set_pressure(0, 55) # 45 for acetone, 70 for ethanol. # 45 is changed to 30 on 20240201 because of overflow.


def stop_flow():
    fgt_set_pressure(0, 0)


def wash_vial():
    ad.options.pen_rate_lower = 3
    ad.update()
    start_flow()
    ad.pendown()
    ad.options.pen_rate_lower = 12
    ad.update()
    time.sleep(2)

    ad.options.pen_pos_up = 35
    ad.update()
    ad.penup()
    time.sleep(1)

    # for i in range(1):
    #     ad.pendown()
    #     time.sleep(0.5)
    #     ad.penup()
    #     time.sleep(1)

    penposups = [50, 60, 75]
    delay_times = [0.7, 1.5, 1.5]

    for i, penposup in enumerate(penposups):
        ad.pendown()
        time.sleep(0.5)
        ad.options.pen_pos_up = penposup
        ad.update()
        ad.penup()
        time.sleep(delay_times[i])

    ad.options.pen_rate_lower = 3
    ad.update()
    time.sleep(1)
    stop_flow()
    ad.pendown()
    x0, y0 = ad.current_pos()

    # this move brings the sucking needle closer to the wall and therefore to the outer rim of the bottom,
    #   thus allowing to extract more liquid from the meniscus near the wall
    ad.lineto(x0-1, y0)

    time.sleep(3)

    ad.lineto(x0, y0)

    ad.options.pen_pos_up = 100
    ad.update()
    ad.penup()


def wash_vials(vial_ids=list(range(54)), idle_run=False, delay=0.5):
    ad.options.pen_pos_up = 100
    ad.update()
    ad.penup()
    time.sleep(1)
    for i, vial in enumerate(vials):
        # print(f'vial: {i}-{vial}')
        if i not in vial_ids:
            continue
        # print(f'Washing vial number {i}')
        assert vial[0]>=0, f'X coordinate of vial {i} is samller than 0!'
        assert vial[1] >= 0, f'Y coordinate of vial {i} is samller than 0!'
        ad.moveto(vial[0], vial[1])
        # print(f'vial location: {vial[0]},{vial[1]}')
        if idle_run:
            time.sleep(delay)
        else:
            time.sleep(0.2)
            wash_vial()
    ad.moveto(0, 0)
    time.sleep(1)
    ad.pendown()
    print('Wash completed')

def corners_check(delay=3):
    wash_vials(idle_run=True, vial_ids=[0, 8, 45, 53], delay=delay)

#wash_vials(idle_run=0)
#corners_check(delay=3)

if __name__ == '__main__':
    assert ad.options.pen_pos_up == 100
    assert ad.options.units == 2

    #corners_check()
    wash_vials(idle_run=False,  vial_ids=list(range(54)))


