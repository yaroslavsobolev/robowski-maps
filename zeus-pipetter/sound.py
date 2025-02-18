import winsound

def beep_for_measurement():
    frequency = 600  # Frequency of the beep sound in Hertz (1000 Hz is a typical value for a beep)
    duration = 600  # Duration of the beep sound in milliseconds (1 second)
    winsound.Beep(frequency, duration)

def beep_for_tip_changing():
    frequency = 800  # Frequency of the beep sound in Hertz (1000 Hz is a typical value for a beep)
    duration = 500  # Duration of the beep sound in milliseconds (1 second)
    for i in range(10):
     winsound.Beep(frequency, duration)

def beep_for_each_pipetting():
    frequency = 1000  # Frequency of the beep sound in Hertz (1000 Hz is a typical value for a beep)
    duration = 600  # Duration of the beep sound in milliseconds (1 second)
    winsound.Beep(frequency, duration)

def beep_for_error():
    for i in range(10):
        frequency = 1500  # Frequency of the beep sound in Hertz (1000 Hz is a typical value for a beep)
        duration = 200  # Duration of the beep sound in milliseconds (1 second)
        winsound.Beep(frequency, duration)

def beep_for_success():
    for i in range(50):
        frequency = 600  # Frequency of the beep sound in Hertz (1000 Hz is a typical value for a beep)
        duration = 300  # Duration of the beep sound in milliseconds (1 second)
        winsound.Beep(frequency, duration)


if __name__ == '__main__':
    # beep_for_each_pipetting()
    beep_for_tip_changing()
    # beep_for_measurement()