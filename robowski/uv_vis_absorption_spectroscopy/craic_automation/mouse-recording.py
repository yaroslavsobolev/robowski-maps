import pyautogui, sys, time
delay = 1

print(pyautogui.locateCenterOnScreen('screenshots/continue.PNG'))

# pyautogui.click('screenshots/continue.PNG')

print('Press Ctrl-C to quit.')
try:
    while True:
        x, y = pyautogui.position()
        positionStr = '(' + str(x).rjust(4) + ', ' + str(y).rjust(4) + ')\n'
        # positionStr = str(pyautogui.position())
        print(positionStr, end='')
        time.sleep(delay)
        print(pyautogui.pixel(1204,  574))
        # print('\b' * len(positionStr), end='', flush=True)
except KeyboardInterrupt:
    print('\n')