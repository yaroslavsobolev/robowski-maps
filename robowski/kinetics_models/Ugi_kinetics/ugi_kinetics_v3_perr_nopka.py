from robowski.settings import *
from ugi_kinetics_v3_emcee import *
from scipy.optimize import minimize
from scipy.optimize import least_squares
import numpy as np

if __name__ == '__main__':
    nll = lambda *args: -log_likelihood_parallel(*args)
    parameter_names = ['k1', 'k_1', 'k2', 'k_2', 'k3', 'k_3', 'k4', 'k_4', 'k5', 'k6', 'k_6', 'k7', 'k8', 'k9', 'km1', 'k_m1', 'km3', 'pKa_ptsa', 'pKa_BuNH3', 'pKa_DMFH', 'pKa_p1', 'pKa_p2', 'pKa_p3', 'pKa_c', 'pKa_ch', 'pKa_q1', 'pKa_q2h', 'pKa_m1h', 'pKa_m4h', 'pKa_me2nh2', 'pKa_q4']

    p0 = np.load('ugi_v3_REV_outputs/theta_v3_REV_opt_85c_2025-02-15.npy')

    print('initial guess:')
    print_parameters_with_names(p0)

    # # make a new function that only uses the first 17 parameters and fizes pKa values to the ones in p0
    # # Then only do optimization of the 17 parameters, which are the rate constants
    # def residuals_parallel_nokpa(p, *args):
    #     p_full = np.copy(p0)
    #     p_full[:17] = p
    #     return residuals_parallel(p_full, *args)

    def residuals_parallel_nokpa_with_timeout(p):
        p_full = np.copy(p0)
        p_full[:17] = p
        return residuals_parallel_with_timeout(p_full)

    x_scale = p0[0:17]

    # bounds for 17 parameters
    # TRF VERSION FOR COVARIANCE MATRIX
    # finite_diff_rel_step = np.array([0.001] * 17)

    # finite_diff_rel_step = np.array([0.5] * 17)
    # res = least_squares(residuals_parallel_nokpa, p0[:17], bounds=(lower_bounds[:17], upper_bounds[:17]), verbose=2, jac='2-point',
    #                     diff_step=finite_diff_rel_step, x_scale=x_scale, xtol=10, gtol=10, ftol=10)

    finite_diff_rel_step = np.array([0.2] * 17)
    finite_diff_rel_step[12] = 0.0001
    # bounds for 17 parameters are not important anymore, because this execution is for uncertainty estimation
    lower_bounds = [0] * 17
    upper_bounds = [1e20] * 17
    res = least_squares(residuals_parallel_nokpa_with_timeout, p0[:17], bounds=(lower_bounds[:17], upper_bounds[:17]), verbose=2, jac='3-point',
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
    np.save('v3_REV_theta_perr_nopka_2025-02-15.npy', perr)
    np.save('v3_REV_theta_pcov_nopka_2025-02-15.npy', cov)

    J = res.jac
    pcov = np.linalg.inv(J.T.dot(J))
    pcov = pcov * (y_sigmas ** 2)
    perr = np.sqrt(np.diag(pcov))
    print(f'default perr: {perr}')

    # save perr to a file
    np.save('default_v3_REV_theta_perr_nokpa_2025-02-15.npy', perr)
    np.save('default_v3_REV_theta_pcov_nopka_2025-02-15.npy', pcov)
