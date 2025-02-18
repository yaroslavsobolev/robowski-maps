import pyautogui, time, os, shutil
import numpy as np
from datetime import datetime
import pandas as pd
import cv2

# cap = cv2.VideoCapture(0, cv2.CAP_DSHOW)
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
screenshots_folder = 'screenshots_raman/'
absolute_path_to_temp_folder = 'C:\\robochem\\uv-vis-absorption-spectroscopy\\craic-automation\\temp_folder'
darkpos = (0, 0)

desktop_icons = {'lambdafire': ( 645,   37),
                 'stage_control': ( 567,   32)}

goto_window = {'ref_default': ( 764,  191),
                     'ref_now': ( 764,  191),
                     'x': ( 882,  169),
                     'y': ( 881,  210),
                     'goto': ( 764,  191)}

main_window = {'ref_default': (1751,  486),
               'ref_now': (1751,  486),
               'collect_dark_and_ref': ( 971, 1107),
               'collect_spectrum': (1260, 1106),
               'beamsplitter_one': (1425,  989),
               'clear_all': ( 726,  489),
               'my_folder': ( 539, 1103),
               'spectrometer_params': ( 905,  490),
               'camera_params': ( 987,  486),
               'Goto': (0, 0)}

popup_window_buttons = {'ref_default': (0, 0),
                        'ref_now': (0, 0),
                        'sampling_time': (1082,  643),
                        'spec_param_ok': ( 994,  770),
                        'scans_to_average': (1083,  571),
                        'absorbance_radiobutton': ( 763,  495),
                        'enter_file_name': (1054,  630),
                        'scan_and_take_image': ( 837,  674),
                        'stage_set_up': (1064,   35),
                        'stage_set_up_2': (1107,   64),
                        'stage_speed' : ( 244,  277),
                        'stage_speed_ok': (( 243,  389)),
                        # 'x': (0, 0),
                        # 'x': (0, 0),
                        # 'x': (0, 0),
                        # 'x': (0, 0),
                        # 'x': (0, 0),
                        # 'x': (0, 0),
                        'x': (0, 0)}

error_window_buttons = {'ref_default': (0, 0),
                        'ref_now': (0, 0),
                        'error_indicator' : (1204,  574),
                        'error_continue_button': (1080,  687),
                        'error_continue_button2': (1080,  720)}

error_indicator_when_error = (140, 221, 43)

stage_control_window = {'motion_indicator': (0, 0)}
indicator_color_when_stage_is_idle = (255, 164, 137)

# set_camera_params_clicks = [
#     (1085,  399),
#     (1087,  598),
#     (1021,  678),
#     ( 835,  400),
#     ( 835,  600),
#     (1115,  858)]

# # if the CRAIC software is not running, start it
# if pyautogui.locateCenterOnScreen(screenshots_folder + 'craictech.png') is None:
#     x, y = desktop_icons['lambdafire']
#     pyautogui.doubleClick(x, y)
#     print('LambdaFire is launching...')
#     time.sleep(20)
#     while pyautogui.locateCenterOnScreen(screenshots_folder + 'craictech.png') is None:
#         print('...')
#         time.sleep(0.2)
#     print('...done.')
#
# if pyautogui.locateCenterOnScreen(screenshots_folder + 'motion_indicator.png') is None:
#     x, y = desktop_icons['stage_control']
#     pyautogui.doubleClick(x, y)
#     print('Stage Control is launching...')
#     while pyautogui.locateCenterOnScreen(screenshots_folder + 'motion_indicator.png') is None:
#         print('...')
#         time.sleep(0.2)
#     print('...done.')
#     x,y = pyautogui.locateCenterOnScreen(screenshots_folder + 'goto_abs.png')
#     pyautogui.click(x, y)
#     time.sleep(0.5)
#     pyautogui.click(1502,   40)
#     time.sleep(0.5)
#     pyautogui.mouseDown()
#     pyautogui.dragTo(1240,   16, 0.5, button='left')

# update the reference locations of windows by finding certain images on the screen
# main_window['ref_now'] = pyautogui.locateCenterOnScreen(screenshots_folder + 'craictech.png')
goto_window['ref_now'] = pyautogui.locateCenterOnScreen(screenshots_folder + 'goto.png')
stage_control_window['motion_indicator'] = pyautogui.locateCenterOnScreen(screenshots_folder + 'motion_indicator.png')

def get_zero_point_by_plate_id(plate_id):
    # if plate_id>34:
    #     return (-50, -693)
    # else:
    #     return (0, 0)
    return (0, 0)


def click_button(window_dictionary, button_name, delay=0.5):
    x = window_dictionary[button_name][0] - window_dictionary['ref_default'][0] + window_dictionary['ref_now'][0]
    y = window_dictionary[button_name][1] - window_dictionary['ref_default'][1] + window_dictionary['ref_now'][1]
    pyautogui.click(x, y)
    time.sleep(delay)

# def set_camera_parameters():
#     click_button(main_window, 'camera_params')
#     time.sleep(0.5)
#     for location in set_camera_params_clicks:
#         pyautogui.click(location[0], location[1])
#     time.sleep(3)

def set_folder_to_temp_folder():
    time.sleep(2)
    click_button(main_window, 'my_folder')
    pyautogui.press('enter')
    pyautogui.write(absolute_path_to_temp_folder)
    pyautogui.press('enter')
    pyautogui.press('tab', presses=4)
    pyautogui.press('space')
    pyautogui.press('tab')
    pyautogui.press('space')

def erase_text_field():
    pyautogui.press('end')
    with pyautogui.hold('ctrl'):
        pyautogui.press('backspace')

def set_spectrometer_parameters(scans_to_average, sampling_time=50):
    click_button(main_window, 'spectrometer_params', delay=1)
    click_button(popup_window_buttons, 'sampling_time')
    erase_text_field()
    pyautogui.write(f'{sampling_time}')
    click_button(popup_window_buttons, 'scans_to_average')
    erase_text_field()
    pyautogui.write(f'{scans_to_average}')
    click_button(popup_window_buttons, 'spec_param_ok')
    time.sleep(3)


def initial_setup():
    # set_spectrometer_parameters(100, 200)
    # set_spectrometer_parameters(10, 50)
    # set_camera_parameters()
    set_folder_to_temp_folder()
    # click_button(main_window, 'beamsplitter_one')


def stage_is_moving():
    return not pyautogui.pixelMatchesColor(int(round(stage_control_window['motion_indicator'][0])),
                                           int(round(stage_control_window['motion_indicator'][1])),
                                           indicator_color_when_stage_is_idle)


def move_to(position):
    x,y = position
    click_button(goto_window, 'x', delay=0.01)
    pyautogui.press('end')
    pyautogui.press('backspace', presses=10)
    pyautogui.write(f'{int(round(x))}')
    click_button(goto_window, 'y', delay=0.01)
    pyautogui.press('end')
    pyautogui.press('backspace', presses=10)
    pyautogui.write(f'{int(round(y))}')
    click_button(goto_window, 'goto')
    # waits until the stage stops by checking the color of the indicator on the Stage Control window
    print(f'Stage is moving to {position}..')
    while stage_is_moving():
        time.sleep(0.1)
        print('...stage is moving...')
    print('...stage motion finished.')



def measure_dark_and_background(scans_to_average=100, delay=60):
    move_to((0, 0))
    set_spectrometer_parameters(scans_to_average, 50)
    click_button(main_window, 'collect_dark_and_ref')
    print(f'Dark & ref, waiting {delay} seconds...')
    time.sleep(delay)
    set_spectrometer_parameters(10, 50)
    move_to(darkpos)


def collecting_spectra():
    return not (pyautogui.locateCenterOnScreen(screenshots_folder + 'collecting_spectra.png') is None)


def check_for_file_absence_error_and_deal_with_it():
    if not pyautogui.pixelMatchesColor(int(round(error_window_buttons['error_indicator'][0])),
                                           int(round(error_window_buttons['error_indicator'][1])),
                                           (255, 255, 255)):
        print('File absence error. Continuing.')
        click_button(error_window_buttons, 'error_continue_button')
        click_button(error_window_buttons, 'error_continue_button2')
        # pyautogui.click('screenshots/continue.PNG')
        time.sleep(0.5)
        check_for_file_absence_error_and_deal_with_it()
    else:
        print(f"Color of error-detection pixel: {pyautogui.pixel(int(round(error_window_buttons['error_indicator'][0])), int(round(error_window_buttons['error_indicator'][1])))}")

# if pyautogui.pixelMatchesColor(int(round(error_window_buttons['error_indicator'][0])),
    #                                        int(round(error_window_buttons['error_indicator'][1])),
    #                                        error_indicator_when_error):
    #     print('File absence error. Continuing.')
    #     click_button(error_window_buttons, 'error_continue_button')
    #     time.sleep(0.5)
    # else:
    #     print(f"Color of error-detection pixel: {pyautogui.pixel(int(round(error_window_buttons['error_indicator'][0])), int(round(error_window_buttons['error_indicator'][1])))}")

def collect_spectrum(id, suffix='', position=None, deltax=400):
    for file_extension in ['msp', 'jpg']:
        file_path = absolute_path_to_temp_folder + f'\\spectrum_-{id}{suffix}.{file_extension}'
        if os.path.exists(file_path):
            os.remove(file_path)

    # set stage speed 100 um/s
    click_button(popup_window_buttons, 'stage_set_up', delay=0.1)
    click_button(popup_window_buttons, 'stage_set_up_2', delay=0.5)
    click_button(popup_window_buttons, 'stage_speed', delay=0.5)
    pyautogui.press('end')
    pyautogui.press('backspace', presses=10)
    pyautogui.write('100')
    click_button(popup_window_buttons, 'stage_speed_ok', delay=0.5)

    click_button(main_window, 'collect_spectrum')
    time.sleep(1)
    check_for_file_absence_error_and_deal_with_it()
    # click_button(popup_window_buttons, 'absorbance_radiobutton', delay=0.01)
    click_button(popup_window_buttons, 'enter_file_name', delay=0.01)
    erase_text_field()
    pyautogui.write(f'spectrum_-{id}{suffix}')
    click_button(popup_window_buttons, 'scan_and_take_image', delay=0.01)
    time_when_exposure_started = time.time()
    time.sleep(0.5)
    print('Collecting spectrum...')
    time.sleep(2)
    current_deltax_sign = 1
    while collecting_spectra():
        current_deltax_sign = -1 * current_deltax_sign
        print(f'...shifting to {current_deltax_sign} times {deltax}')
        move_to((position[0] + current_deltax_sign * deltax, position[1]))
        time.sleep(0.1)
    print(f'Collected spectrum with id {id}{suffix}')

    # set stage speed 10000 um/s
    click_button(popup_window_buttons, 'stage_set_up', delay=0.1)
    click_button(popup_window_buttons, 'stage_set_up_2', delay=0.5)
    click_button(popup_window_buttons, 'stage_speed', delay=0.5)
    pyautogui.press('end')
    pyautogui.press('backspace', presses=10)
    pyautogui.write('10000')
    click_button(popup_window_buttons, 'stage_speed_ok', delay=0.5)
    return time_when_exposure_started

def copy_temp_folder_to_archive(plate_id, exp_name):
    now = datetime.now()
    datetime_string_now = now.strftime("%Y-%m-%d_%H-%M-%S")
    unixtime = time.mktime(now.timetuple())
    target_folder_name = f'{datetime_string_now}__plate{plate_id}__{exp_name}'
    shutil.copytree(absolute_path_to_temp_folder, data_folder + 'craic_microspectrometer_measurements/raman/' + target_folder_name)
    df = pd.read_csv(data_folder + 'craic_microspectrometer_measurements/raman/database_about_these_folders.csv')
    df = df.append({'timestamp':unixtime, 'datetime':datetime_string_now, 'plate_id':plate_id, 'exp_name':exp_name, 'folder':target_folder_name},
                   ignore_index=True)
    df.to_csv(data_folder + 'craic_microspectrometer_measurements/raman/database_about_these_folders.csv',
              index=False)
    # # This initialized the archive CSV. Never run.
    # df = pd.DataFrame([[unixtime, datetime_string_now, plate_id, exp_name, target_folder_name]],
    #                   columns=['timestamp', 'datetime', 'plate_id', 'exp_name', 'folder'])


# x_range = 103581,
# y_range = 65166,
def measure_one_plate(id, experiment_name='test',
                      zero_pos=(0, 0),
                      darkpos=(0, 0),
                      x_range=103993,
                      y_range=64928,
                      nwells_x=9, nwells_y=6):
    click_button(main_window, 'clear_all')
    click_button(main_window, 'clear_all')
    move_to(zero_pos)
    print('Initial position reached.')
    xs = np.linspace(0, x_range, nwells_x)
    ys = np.linspace(y_range, 0, nwells_y) # scanning in reverse because the bottles are flipped, hence the plate pattern is mirrored
    well_id = 1
    endtimes = []
    starttimes = []
    for x in xs:
        for y in ys:
            move_to((x, y))
            t0 = time.time()
            starttimes.append(collect_spectrum(id=well_id, suffix='', position=(x, y)) - t0)
            endtimes.append(time.time() - t0)
            well_id += 1
        # # this next move is to avoid diagonal motion and instead do a "knight move"
        # # The goal is to avoid illuminating wells that were not yet measured and thus to avoid photobleaching.
        # print('Knight move...')
        # move_to((x, y_range))
    np.savetxt(absolute_path_to_temp_folder + '\\starttimes.txt', np.array(starttimes))
    np.savetxt(absolute_path_to_temp_folder + '\\endtimes.txt', np.array(endtimes))
    print('Measuring complete. Copying data to archive.')
    copy_temp_folder_to_archive(plate_id=id, exp_name=experiment_name)
    print('Data copied to archive.')
    move_to(darkpos)
    print('Returned to zero position.')


def scan_barcode():
    ret, img = cap.read()
    # img = np.flip(img, axis=0)
    # img = np.flip(img, axis=1)
    img = img[0:135, 60:600, :]
    # img = cv2.resize(img, dsize=(img.shape[1], img.shape[0]*3), interpolation=cv2.INTER_CUBIC)
    img_for_scanner = cv2.bitwise_not(cv2.resize(img, dsize=(img.shape[1], img.shape[0]*6), interpolation=cv2.INTER_CUBIC))
    img_for_scanner = cv2.normalize(img_for_scanner, None, 0, 255, cv2.NORM_MINMAX)
    cv2.imwrite("barcode_frame.jpg", np.flip(np.flip(img, axis=0), axis=1))
    cv2.imwrite("barcode_frame_for_scanner.jpg", img_for_scanner)
    bd = cv2.barcode.BarcodeDetector()
    retval, decoded_info, decoded_type, points = bd.detectAndDecode(img_for_scanner)
    if retval:
        # Remove the last digit because the last digit in EAN-8 protocol is the checksum digit ("check digit")
        goodinfo = [x[:-1] for x in decoded_info if x != '']
        # print(f'Barcode found: {decoded_info}, good one is {goodinfo}')
        print(f'Barcode found: {goodinfo}')
        return retval, goodinfo
    else:
        print('No barcode found')
        return retval, decoded_info
    # cap.release()

def n_last_codes_are_the_same(x):
    return (x.count(x[0]) == len(x)) and (x[0] != '')

def repeat_single_well(ntimes):
    t0 = time.time()
    print('Repeating one well')
    endtimes = []
    starttimes = []
    for i in range(ntimes):
        print(f'Repetition {i}')
        time_of_exposure_start = collect_spectrum(i+1)
        starttimes.append(time_of_exposure_start - t0)
        endtimes.append(time.time() - t0)
    times = np.array(endtimes)
    np.savetxt(absolute_path_to_temp_folder + '\\starttimes.txt', starttimes)
    np.savetxt(absolute_path_to_temp_folder + '\\endtimes.txt', endtimes)



if __name__ == '__main__':

    for repid in range(19):
        # print(pyautogui.locateCenterOnScreen(screenshots_folder + 'craictech.png'))
        experiment_name = \
    f'12testArep{repid+20}_2024-10-22-run01_'
        plate_barcode = 67
        click_button(main_window, 'clear_all')
        click_button(main_window, 'clear_all')

        # repeat_single_well(54)

        # initial_setup()
        # measure_dark_and_background()
        move_to(darkpos)

        yrange = 65763
        measure_one_plate(id=plate_barcode,
                          experiment_name=experiment_name,
                          zero_pos=(0, yrange),
                          darkpos=(0, yrange),
                          x_range=104305,
                          y_range=yrange,
                          nwells_x=9, nwells_y=6)

        # # automatic operation with barcode scanner
        # scanned_barcodes = ['']*10
        # last_measured_barcode = ''
        # while True:
        #     time.sleep(1)
        #     retval, decoded_info = scan_barcode()
        #     if not retval:
        #         scanned_barcodes.pop(0)
        #         scanned_barcodes.append('')
        #         continue
        #     else:
        #         scanned_barcodes.pop(0)
        #         scanned_barcodes.append(decoded_info[0])
        #
        #     if n_last_codes_are_the_same(scanned_barcodes):
        #         if (scanned_barcodes[-1] == last_measured_barcode):
        #             print('This plate has been just measured.')
        #         else:
        #             measure_one_plate(id=scanned_barcodes[-1],
        #                               experiment_name=experiment_name)
        #             last_measured_barcode = scanned_barcodes[-1]
        #             scanned_barcodes = [''] * 10