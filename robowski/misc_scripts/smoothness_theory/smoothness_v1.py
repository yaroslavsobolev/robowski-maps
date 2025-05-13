from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from tqdm import tqdm
from multiprocessing import Pool

def right_hand_side_of_ode(t, y):
    deriv_scale = 2
    derivs = np.zeros_like(y)
    derivs[0] = -y[0] * y[0]
    for i in range(1, len(y) - 1):
        derivs[i] = y[i - 1] * y[i - 1] - y[i] * y[i]
    derivs[-1] = y[-2] * y[-2]

    return derivs * deriv_scale

def solve_ode(t, y0):
    return solve_ivp(right_hand_side_of_ode, (t[0], t[-1]), y0, t_eval=t)

if __name__ == '__main__':
    # t = np.linspace(0, 100, 100)
    # N = 10
    # y0 = np.zeros(N)
    # C0 = 2
    # y0[0] = C0
    #
    # sol = solve_ode(t, y0)
    #
    # plt.plot(sol.t, sol.y.T, label='Solution')
    # plt.legend()
    # plt.show()
    #
    # y0 = np.zeros(N)
    # delta = 0.01
    # y0[0] = C0 + delta
    # sol2 = solve_ode(t, y0)
    # y2_deriv = (sol2.y.T - sol.y.T)/delta
    # plt.plot(sol.t, y2_deriv, label=range(N))
    # plt.ylabel('Partial derivative of last product by initial concentration of substrate')
    # plt.xlabel('Time')
    # plt.legend()
    # plt.show()

    xs = np.linspace(0, 1, 100)
    kt = 2
    ys = 1 / (1/xs + kt)
    plt.plot(xs, ys, color='C1', linewidth=3)
    plt.xlabel('Starting concentration of substrate')
    plt.ylabel('Final concentration of product')
    # remove ticklabels and ticks
    plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    # plot point and tangent using the derivative of ys at n=25 point
    n = 25
    plt.scatter(xs[n], ys[n], color='red')
    deltax = 0.3
    derivative = (ys[n+1] - ys[n-1])/(xs[n+1] - xs[n-1])
    plt.plot([xs[n] - deltax, xs[n] + deltax], [ys[n] - derivative*deltax, ys[n] + derivative*deltax], '--', color='black')

    plt.show()



