from ugi_kinetics_v4_emcee import *
from scipy.optimize import minimize
from scipy.optimize import least_squares
import numpy as np

if __name__ == '__main__':
    nll = lambda *args: -log_likelihood_parallel(*args)

    finite_diff_rel_step = np.array([0.02]*17 + [0.002] * 14 + [0.002] * 2 + [0.05] * 3)
    finite_diff_rel_step[8] = 0.1

    theta = np.load('theta_v3_REV_intermediate_during_opt.npy')
    theta[8] = theta[8] / 6
    new_params = [0.1, 12, 1000, 30, 3]
    # append the new params to the theta
    theta = np.append(theta, new_params)
    p0 = theta

    x_scale = p0 * 1
    x_scale[17:33] = 1

    # TRF VERSION
    res = least_squares(residuals_parallel, p0, bounds=(lower_bounds, upper_bounds), verbose=2, jac='2-point',
                        diff_step=finite_diff_rel_step, xtol=None, x_scale=x_scale)

    # # DOGBOX version
    # res = least_squares(residuals_parallel, p0, bounds=(lower_bounds, upper_bounds), verbose=2, jac='2-point',
    #                     diff_step=finite_diff_rel_step, xtol=1e-12, method='dogbox')

    # save xoln.x to a file
    np.save('ugi_v4_theta_opt.npy', res.x)
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
    np.save('ugi_v4_theta_perr_opt.npy', perr)