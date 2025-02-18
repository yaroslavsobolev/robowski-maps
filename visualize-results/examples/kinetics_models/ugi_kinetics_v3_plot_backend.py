from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import emcee
from _decimal import Decimal
from tqdm import tqdm
from multiprocessing import Pool
import corner

if __name__ == '__main__':

    # # server settings
    filename = 'visualize-results/examples/kinetics_models/ugi_v3_outputs/ugi_v2_emcee_backend_gauss_2025_01_14a.h5'

    sampler = emcee.backends.HDFBackend(filename)

    # nplots = 5
    # fig, axes = plt.subplots(nplots, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    # # labels = range(len(p0))
    # for i in range(nplots):
    #     ax = axes[i]
    #     ax.plot(samples[:, :, i], "k", alpha=0.3)
    #     ax.set_xlim(0, len(samples))
    #     # ax.set_ylabel(labels[i])
    #     ax.yaxis.set_label_coords(-0.1, 0.5)
    #
    # axes[-1].set_xlabel("step number")
    # plt.show()

    tau = sampler.get_autocorr_time()
    print(f'tau: {tau}')

    reader = emcee.backends.HDFBackend(filename)

    samples = reader.get_chain()

    # tau = reader.get_autocorr_time()
    # burnin = int(2 * np.max(tau))
    # thin = int(0.5 * np.min(tau))
    # samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
    # log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)
    # log_prior_samples = reader.get_blobs(discard=burnin, flat=True, thin=thin)

    # print("burn-in: {0}".format(burnin))
    # print("thin: {0}".format(thin))
    print("flat chain shape: {0}".format(samples.shape))
    # print("flat log prob shape: {0}".format(log_prob_samples.shape))
    # print("flat log prior shape: {0}".format(log_prior_samples.shape))

    # flat_samples = sampler.get_chain(discard=30, thin=5, flat=True)
    # print(flat_samples.shape)
    #
    # best_index = np.argmax(sampler.flatlnprobability)
    # theta_max = sampler.flatchain[best_index]
    # np.save('theta_best.npy', theta_max)
    #
    # print(f'Best theta: {theta_max}')


    labels = [str(i) for i in range(31)]
    fig = corner.corner(
        samples, labels=labels
    )

    fig.savefig('visualize-results/examples/kinetics_models/ugi_v3_outputs/tempcorner.png', dpi=300)
    plt.show()

