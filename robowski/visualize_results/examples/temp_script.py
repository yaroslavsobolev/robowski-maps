from robowski.settings import *
import numpy as np
import robowski.uv_vis_absorption_spectroscopy.spectraltools as spt
import matplotlib.pyplot as plt

import importlib_resources
repo_data_path = str(importlib_resources.files("robowski")) + '/'

if __name__ == '__main__':
    my_resources = importlib_resources.files("robowski")
    filename = my_resources / repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/2022-12-01/interpolator-dataset/agilent_data.npy'
    filename2 = repo_data_path + repo_data_path + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/2022-12-01/interpolator-dataset/agilent_data.npy'
    data = np.load(filename2)
    print(data)

    plt.plot([0, 1], [1,3])
    spt.simpleaxis(plt.gca())

    plt.show()

