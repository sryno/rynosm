from __future__ import print_function
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import argparse


__author__ = 'Sean M. Ryno'
__copyright__ = 'Copyright 2017, Sean M. Ryno'
__credits__ = 'Sean M. Ryno'
__license__ = 'GPL v3.0'
__version__ = '0.2'
__maintainer__ = 'Sean M. Ryno'
__email__ = 'sean.m.ryno@gmail.com'
__status__ = 'Development'


def parseTXT(inFile):
    fin = open(inFile, 'r')
    angles, values = [], []
    for line in fin:
        if line.strip() != '':
            angles.append(float(line.split()[0]))
            values.append(float(line.split()[1]))
        else:
            pass

    return angles, values


def ang2rad(angle):
    return np.radians(angle)

def changePhase(angle):
    return angle + 180.0


def rbFunc(x, c0, c1, c2, c3, c4, c5):
    return c0 + (c1 * np.cos(x)) + (c2 * np.cos(x)**2) + (c3 * np.cos(x)**3) + (c4 * np.cos(x)**4) + (c5 * np.cos(x)**5)


def plotData(xdata, ydata):
    plt.plot(xdata, ydata, 'bo', label='data')
    popt, pcov = curve_fit(rbFunc, xdata, ydata, ftol=1e-15, xtol=1e-15, gtol=1e-15)
    # plt.plot(xdata, rbFunc(xdata, *popt), 'r-', label='fit')
    xnew = np.linspace(min(xdata), max(xdata), 300)
    xSmooth = interp1d(xdata, rbFunc(xdata, *popt), kind='cubic')
    ynew = xSmooth(xnew)
    residuals = ydata - rbFunc(xdata, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('')
    print('C0 = {:>8.3f}'.format(popt[0]))
    print('C1 = {:>8.3f}'.format(popt[1]))
    print('C2 = {:>8.3f}'.format(popt[2]))
    print('C3 = {:>8.3f}'.format(popt[3]))
    print('C4 = {:>8.3f}'.format(popt[4]))
    print('C5 = {:>8.3f}'.format(popt[5]))
    print('')
    print('R^2 = {:>7.4f}'.format(r_squared))
    plt.plot(xnew, ynew, 'r-', label='fit')
    plt.xlabel('Angle (radians)')
    plt.ylabel('Energy')
    plt.legend()
    plt.show()

    return


if __name__ == '__main__':
    """
    Fits a Ryckaert-Bellemans Dihedral to a set of points.
    """

    parser = argparse.ArgumentParser(description='Fits a Ryckaert-Bellemans Dihedral to a set of points. Angles should start at 0 or 180.')
    parser.add_argument('INPUT_FILE', nargs=1, help='Input text file with ANGLE VALUE on each line. Angle should be in degrees.')
    args = parser.parse_args()

    ang, val = parseTXT(vars(args)['INPUT_FILE'][0])

    if ang[0] < 180.0:
        for i in range(len(ang)):
            ang[i] = changePhase(ang[i])

    for i in range(len(ang)):
        ang[i] = ang2rad(ang[i])

    plotData(ang, val)
