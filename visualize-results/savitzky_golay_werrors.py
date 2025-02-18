import numpy as np
from scipy.linalg import svd, inv
from numpy.polynomial import polynomial as P
from numpy.polynomial.polynomial import polyvander

def solve_polyfit(xarr, yarr, degree, weight, deriv=None):
    # Obtain solution using polyfit
    z = P.polyfit(xarr, yarr, deg=degree, w=weight)

    if deriv is not None:
        z = P.polyder(z, m=deriv)
    return z

def solve_leastsq(yarr, ycov, vander, vanderT, deriv=None):
    ycovinv = inv(ycov)
    prefactor = inv(np.dot(np.dot(vanderT, ycovinv), vander))
    postfactor = np.dot(np.dot(vanderT, ycovinv), yarr)
    z = np.dot(prefactor, postfactor)

    if deriv is not None:
        z = P.polyder(z, m=deriv)

    return z


'''
Simple Savitzky Golay interpolator, assuming the data has errorbars.
Utilizes numpy.polyfit
Author: Surhud More (Kavli IPMU)
Email bug reports to surhudkicp [at] gmail.com

TODO:
    Python 3 compatibility
    Check that window_length and degree are integers
    Checks that y and error are arrays
'''
def savgol_filter_werror(y, window_length, degree, error=None, cov=None,
        deriv=None):
    ynew = y * 0.0

    # Check that window_length is odd
    if window_length % 2 == 0:
        print("Window length must be odd\n")
        exit(11)

    # Take care that the window does not spill out of our array
    margin = int(window_length/2)
    xarr = np.arange(-margin, margin+1)

    if cov is not None:
        vander = polyvander(xarr, deg=degree)
        vanderT = np.transpose(vander)
    else:
        weight = 1./error


    for i in range(margin, y.size-margin):
        if cov is None:
            z = solve_polyfit(xarr,
                              y[i-margin:i+margin+1],
                              degree,
                              weight[i- margin:i+margin+1],
                              deriv=deriv)
        else:
            z = solve_leastsq(y[i-margin:i+margin+1],
                              cov[i-margin:i+margin+1,i-margin:i+margin+1],
                              vander,
                              vanderT,
                              deriv=deriv)
        ynew[i] = P.polyval(0.0, z)


    # Now fit the left boundary, by fitting the first window_length points with
    # a degree order polynomial
    if cov is None:
        z = solve_polyfit(xarr,
                          y[:window_length],
                          degree,
                          weight[:window_length],
                          deriv=deriv)
    else:
        z = solve_leastsq(y[:window_length],
                          cov[:window_length, :window_length],
                          vander,
                          vanderT,
                          deriv=deriv)
    for i in range(margin):
        ynew[i] = P.polyval(xarr[i], z)

    # Now fit the right boundary, by fitting the last window_length points with
    # a degree order polynomial
    if cov is None:
        z = solve_polyfit(xarr,
                          y[-window_length:],
                          degree,
                          weight[-window_length:],
                          deriv=deriv)
    else:
        z = solve_leastsq(y[-window_length:],
                          cov[-window_length:, -window_length:],
                          vander,
                          vanderT,
                          deriv=deriv)

    for i in range(margin):
        ynew[y.size-margin+i] = P.polyval(xarr[i+margin+1], z)

    return ynew
