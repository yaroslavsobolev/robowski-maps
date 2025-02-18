import time

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

pd1 = pd.read_csv(data_folder + 'simple-reactions\\2023-12-02-run01\\'
                                 'nanodrop_spectra\\2023-12-05_00-09-38_2023_12_04_UV-Vis_plate_66.csv')
pd2 = pd.read_csv(data_folder + 'simple-reactions\\2023-12-02-run01\\'
                                    'nanodrop_spectra\\2023-12-05_00-30-45_2023_12_04_UV-Vis_plate_66_autolength.csv')

## take the first 10 columns of pd1
pd1 = pd1.iloc[:, :10]

## devide each column other than the column 0 of pd2 by 33
pd2.iloc[:, 1:] = pd2.iloc[:, 1:] / 33

## plot the column 1 against column 0 of pd1 and pd2
plt.plot(pd1.iloc[:, 0], pd1.iloc[:, 9], label='pd1')
plt.plot(pd2.iloc[:, 0], pd2.iloc[:, 9], label='pd2')
plt.legend(loc="upper right")
plt.xlabel('wavelength (nm)')
plt.ylabel('abs')
plt.show()
