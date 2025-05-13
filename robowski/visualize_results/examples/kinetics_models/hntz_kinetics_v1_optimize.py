from hntz_kinetics_v1_emcee import *
from scipy.optimize import minimize
from scipy.optimize import least_squares
import numpy as np
import logging

if __name__ == '__main__':
    # set logging level to ERROR
    # logging.basicConfig(level=logging.ERROR)

    finite_diff_rel_step = np.array([0.001]*8)

    # factor = 500
    # p0 = np.array([factor * 0.01, factor * 1e-7,
    #                factor * 0.1, factor * 1e-7,
    #                factor * 0.1 / 1000, factor * 0.01 / 1000,
    #                factor * 0.1 / 1000, factor * 0.01 / 1000])

    # p0 = np.array([8.10448175e+00, 1.36644571e-03,
    #                2.45759559e+02, 1.84622923e-02,
    #                1.05885445e+00, 5.63083827e+01,
    #                5.02497025e-02, 7.90050804e+00])
    #
    # # after diagonal opt
    # p0 = np.array([1.65313301e+02, 3.17122692e-02,
    #                2.95588818e+02, 2.05623993e-03,
    #                6.75413472e-03, 4.08117605e+01,
    #                8.71656953e+00, 8.59488696e+01])
    #
    # p0 = np.array([1.94178004e+02, 3.38318064e-02,
    #                4.87277423e+02, 9.73005996e-03,
    #                3.37912085e-03, 4.08440879e+01,
    #                6.59385331e+00, 8.78474051e+01])

    p0 = np.load('theta_opt_2024-09-25.npy')

    # TRF VERSION
    res = least_squares(residuals_parallel, p0, bounds=(lower_bounds, upper_bounds), verbose=2, jac='2-point',
                        diff_step=finite_diff_rel_step, xtol=None)

    # # DOGBOX version
    # res = least_squares(residuals_parallel, p0, bounds=(lower_bounds, upper_bounds), verbose=2, jac='2-point',
    #                     diff_step=finite_diff_rel_step, xtol=1e-12, method='dogbox')

    # save xoln.x to a file
    np.save('theta_opt.npy', res.x)
    print(f'res success: {res.success}\nres message: {res.message}\nsolution:{res.x}')

    # J = res.jac
    # pcov = np.linalg.inv(J.T.dot(J))
    # perr = np.sqrt(np.diag(pcov))
    # print(f'perr: {perr}')
    #
    # # save perr to a file
    # np.save('theta_perr_opt.npy', perr)

    y_sigmas = 0.02
    U, s, Vh = np.linalg.svd(res.jac, full_matrices=False)
    tol = np.finfo(float).eps * s[0] * max(res.jac.shape)
    w = s > tol
    cov = (Vh[w].T / s[w] ** 2) @ Vh[w]  # robust covariance matrix
    # scale the cov matrix by the sigmas of the y values
    cov = cov * y_sigmas ** 2
    perr = np.sqrt(np.diag(cov))  # 1sigma uncertainty on fitted parameters

    print(f'perr: {perr}')

    # save perr to a file
    np.save('theta_perr_opt.npy', perr)