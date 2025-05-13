from robowski.settings import *
from ugi_kinetics_v3r_emcee import *
from scipy.optimize import least_squares
import numpy as np

if __name__ == '__main__':
    nll = lambda *args: -log_likelihood_parallel(*args)

    old_parameter_names = ['k1', 'k_1', 'k2', 'k_2', 'k3', 'k_3', 'k4', 'k_4', 'k5', 'k6', 'k_6', 'k7', 'k8', 'k9', 'km1', 'k_m1', 'km3']
    parameter_names     = ['k1', 'k_1', 'k2',  'K2', 'K3', 'k_3', 'k4', 'K4', 'k5', 'k6', 'k_6', 'k7', 'k8', 'k9', 'km1', 'k_m1', 'Km3']

    p0 = np.load('theta_v3_REV_intermediate_during_opt.npy')

    print('initial guess before variable change:')
    for i, name in enumerate(old_parameter_names):
        print(f'{name}: {p0[i]}')

    # make a new function that only uses the first 17 parameters and fizes pKa values to the ones in p0
    # Then only do optimization of the 17 parameters, which are the rate constants
    def residuals_parallel_nokpa(p, *args):
        p_full = np.copy(p0)
        p_full[:17] = p
        return residuals_parallel(p_full, *args)

    def residuals_parallel_nokpa_with_timeout(p, *args):
        p_full = np.copy(p0)
        p_full[:17] = p
        return residuals_parallel_with_timeout(p_full, *args)

    p0_reform = np.copy(p0)

    # express the new parameters through old kinetic parameters by solving the following expressions:

    # p0_reform[1] = p0[0] / p0[1] # K1 = k1 / k_1
    p0_reform[3] = p0[2] / p0[3] # K2 = k2 / k_2
    # Km3 = km3 * k_m1
    p0_reform[16] = p0[16] * p0[15]

    # # K8 = K2 * k8
    # p0_reform[12] = p0_reform[3] * p0[12]
    # K3 = K1 / k3
    p0_reform[4] = p0[0] / p0[4]
    p0_reform[7] = p0[6] / p0[7]

    print('initial guess after variable change:')
    for i, name in enumerate(parameter_names):
        print(f'{name}: {p0_reform[i]}')

    x_scale = p0_reform[0:17]

    # TRF VERSION FOR COVARIANCE MATRIX
    finite_diff_rel_step = np.array([0.2] * 17)
    finite_diff_rel_step[12] = 0.0001
    # bounds for 17 parameters are not important anymore, because this execution is for uncertainty estimation
    lower_bounds = [0] * 17
    upper_bounds = [1e20] * 17
    res = least_squares(residuals_parallel_nokpa_with_timeout, p0_reform[:17], bounds=(lower_bounds, upper_bounds), verbose=2, jac='3-point',
                        diff_step=finite_diff_rel_step, x_scale=x_scale, xtol=10, gtol=10, ftol=10)

    print(f'res success: {res.success}\nres message: {res.message}\nsolution:{res.x}')

    # In this situation, and assuming one can trust the input sigma uncertainties, one can use the output
    # Jacobian matrix jac returned by least_squares to estimate the covariance matrix.
    # Moreover, assuming the covariance matrix is diagonal, or simply ignoring non-diagonal terms, one can also obtain
    # the 1σ uncertainty perr in the model parameters (often called "formal errors")
    # as follows (see Section 15.4.2 of Numerical Recipes 3rd ed.)
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
    np.save('v3r_REV_theta_perr_nopka.npy', perr)
    np.save('v3r_REV_theta_pcov_nopka.npy', cov)

    J = res.jac
    pcov = np.linalg.inv(J.T.dot(J))
    pcov = pcov * (y_sigmas ** 2)
    perr = np.sqrt(np.diag(pcov))
    print(f'default perr: {perr}')

    # save perr to a file
    np.save('default_v3r_REV_theta_perr_nokpa.npy', perr)
    np.save('default_v3r_REV_theta_pcov_nopka.npy', pcov)
