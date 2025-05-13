import numpy as np
import importlib_resources
import robowski.uv_vis_absorption_spectroscopy.spectraltools as spt
import matplotlib.pyplot as plt

if __name__ == '__main__':
    my_resources = importlib_resources.files("robowski")
    filename = my_resources / 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/2022-12-01/interpolator-dataset/agilent_data.npy'
    f = str(my_resources)
    print(f)
    filename2 = f + '/' + 'uv_vis_absorption_spectroscopy/microspectrometer-calibration/2022-12-01/interpolator-dataset/agilent_data.npy'
    data = np.load(filename2)
    print(data)

    plt.plot([0, 1], [1,3])
    spt.simpleaxis(plt.gca())

    plt.show()

