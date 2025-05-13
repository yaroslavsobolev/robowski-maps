from ugi_kinetics_v4_emcee import *

if __name__ == '__main__':

    theta = np.load('ugi_v4_REV_outputs/ugi_v4_theta_final.npy')

    print(theta)
    assert len(theta) == 36
    plot_yields(filename=None, theta=theta, sizefactor=3, show_legends=False, same_ylims=True)
