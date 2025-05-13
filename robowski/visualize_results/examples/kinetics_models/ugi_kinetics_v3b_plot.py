from robowski.settings import *
from ugi_kinetics_v3b_emcee import *

if __name__ == '__main__':
    # theta = np.array([  2e+02, 1.73297381e+00,  # k1, k_1, \
    #                     1.11069998e+02, 2.42879249e+00,  # k2, k_2, \
    #                     1.36625198e-05, 5.25799912e-02,  # k3, k_3, \
    #                     2.01185205e-01, 1.44369384e-05,  # k4, k_4, \
    #                     8.792766806237957, # k5, \
    #                     0.07, 0.00014145287902037243,      # km1, k_m1, \
    #                     -1.9257310211345142, # pKa_ptsa, \
    #                     12.842106004719119, # pKa_BuNH3, \ # 12.842106004719119
    #                     7.8645547494439295, # pKa_p1, \
    #                   ])

    # plot_yields(filename=None, theta=theta)
    # plot_yields(filename=None, theta=theta, do_load=True, do_show=True)
    # plot_yields(filename='theta_intermediate_during_v3bREV_opt.npy', sizefactor=3, show_legends=False, same_ylims=True)
    plot_yields(filename='theta_v3bREV_opt.npy', sizefactor=3, show_legends=False, same_ylims=True)
    # plot_yields(filename='theta_opt.npy', sizefactor=3, show_legends=False, same_ylims=True)