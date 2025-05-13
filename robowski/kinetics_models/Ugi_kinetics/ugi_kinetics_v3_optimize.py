from robowski.settings import *
from ugi_kinetics_v3_emcee import *
from scipy.optimize import minimize
from scipy.optimize import least_squares
import numpy as np

if __name__ == '__main__':
    nll = lambda *args: -log_likelihood_parallel(*args)

    finite_diff_rel_step = np.array([0.02]*17 + [0.002] * 14)

    parameter_names = ['k1', 'k_1', 'k2', 'k_2', 'k3', 'k_3', 'k4', 'k_4', 'k5', 'k6', 'k_6', 'k7', 'k8', 'k9', 'km1', 'k_m1', 'km3', 'pKa_ptsa', 'pKa_BuNH3', 'pKa_DMFH', 'pKa_p1', 'pKa_p2', 'pKa_p3', 'pKa_c', 'pKa_ch', 'pKa_q1', 'pKa_q2h', 'pKa_m1h', 'pKa_m4h', 'pKa_me2nh2', 'pKa_q4']

    p0 = np.array([3.19062740e+03, 4.23573978e+01, 1.16954112e+04, 2.25787748e+05,
        1.68951419e+04, 5.19849862e-02, 2.90540599e+01, 3.14087236e+02,
         3.98937158e+00, 2.90783462e+04, 3.24156477e+01, 4.41237515e-03,
         1.62517291e+12, 1.96772914e+01, 4.39400954e-01, 4.98340123e+03,
         4.06512128e+03, -8.03814143e-01, 1.13990883e+01, -7.63072776e+00,
         7.95663838e+00, 8.91516973e+00, 9.00000706e+00, 1.51030182e+01,
         -2.96035411e+00, 1.57570445e+01, -2.32678186e-01, 1.45038405e+01,
         -5.83141924e-01, 7.57789143e+00, 1.28345461e+01])
    p0[21] = 7.3
    p0[22] = 12
    p0[30] = 5.9
    p0[29] = 11.97


    print('initial guess:')
    print_parameters_with_names(p0)

    # print each parameters, comparing with lower_bounds and upper_bounds
    for i in range(len(p0)):
        print(f'{parameter_names[i]}: {p0[i]}')
        print(f'lower_bound: {lower_bounds[i]}, p0 above lower bound: {p0[i] > lower_bounds[i]}')
        print(f'upper_bound: {upper_bounds[i]}, p0 below upper bound: {p0[i] < upper_bounds[i]}')

    x_scale = p0 * 1
    x_scale[17:] = 1

    # # TRF VERSION
    # res = least_squares(residuals_parallel, p0, bounds=(lower_bounds, upper_bounds), verbose=2, jac='2-point',
    #                     diff_step=finite_diff_rel_step, xtol=None, x_scale=x_scale)

    # TRF VERSION FOR COVARIANCE MATRIX
    finite_diff_rel_step = np.array([0.001] * 17 + [0.0001] * 14)
    res = least_squares(residuals_parallel, p0, bounds=(lower_bounds, upper_bounds), verbose=2, jac='2-point',
                        diff_step=finite_diff_rel_step, x_scale=x_scale, xtol=10, gtol=10, ftol=10)

    # save xoln.x to a file
    np.save('v3_theta_opt.npy', res.x)
    print(f'res success: {res.success}\nres message: {res.message}\nsolution:{res.x}')

    # J = res.jac
    # pcov = np.linalg.inv(J.T.dot(J))
    # perr = np.sqrt(np.diag(pcov))
    # print(f'perr: {perr}')
    #
    # # save perr to a file
    # np.save('theta_perr_opt.npy', perr)

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
    np.save('v3_REV_theta_perr_opt.npy', perr)
    np.save('v3_REV_theta_pcov_opt.npy', cov)

    J = res.jac
    pcov = np.linalg.inv(J.T.dot(J))
    pcov = pcov * (y_sigmas ** 2)
    perr = np.sqrt(np.diag(pcov))
    print(f'default perr: {perr}')

    # save perr to a file
    np.save('default_v3_REV_theta_perr_opt.npy', perr)
    np.save('default_v3_REV_theta_pcov_opt.npy', pcov)
