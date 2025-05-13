from robowski.settings import *
from ugi_kinetics_v3_emcee import *

if __name__ == '__main__':

    p0 = np.load('theta_v3_REV_opt_85c_2025-02-15.npy')

    plot_yields(filename=None, theta=p0, sizefactor=3, show_legends=False, same_ylims=True)


    # plot_yields(filename=None, theta=theta, do_load=True, do_show=True)