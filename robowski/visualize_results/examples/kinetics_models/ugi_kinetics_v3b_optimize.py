from ugi_kinetics_v3b_emcee import *
from scipy.optimize import minimize
from scipy.optimize import least_squares

if __name__ == '__main__':
    # nll = lambda *args: -log_likelihood_parallel(*args)

    finite_diff_rel_step = [0.02]*11 + [0.005] + [0.0001] * 2

    p0 = np.load('theta_intermediate_during_v3bREV_opt.npy')

    # for the kinetic rates, the x_scale should be the same as the initial guess
    x_scale = p0 * 1
    # for the pKa values, the x_scale should be 1
    x_scale[11:] = 1


    nll = lambda *args: -log_likelihood_parallel(*args)

    # TRF VERSION
    res = least_squares(residuals_parallel, p0, bounds=(lower_bounds, upper_bounds), verbose=2, jac='2-point',
                        diff_step=finite_diff_rel_step, xtol=None, x_scale=x_scale)

    # save xoln.x to a file
    np.save('theta_v3bREV_opt.npy', res.x)
    print(f'res success: {res.success}\nres message: {res.message}\nsolution:{res.x}')