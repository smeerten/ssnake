#!/usr/bin/env python

# Copyright 2016 - 2017 Bas van Meerten and Wouter Franssen

# This file is part of ssNake.
#
# ssNake is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ssNake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ssNake. If not, see <http://www.gnu.org/licenses/>.


import numpy as np
from scipy.special import wofz



def voigtLine(x, pos, lor, gau, integral, Type = 0):
    lor = np.abs(lor)
    gau = np.abs(gau)
    axis = x - pos
    if Type == 0: #Exact: Freq domain simulation via Faddeeva function
        if gau == 0.0: #If no gauss, just take lorentz
           lor = 1.0 / (np.pi * 0.5 * lor * (1 + (axis /(0.5 * lor))**2) )
           return integral * lor
        else:
            sigma = gau / (2 * np.sqrt(2 * np.log(2)))
            z = (axis + 1j * lor / 2) / (sigma * np.sqrt(2))
            return integral * wofz(z).real / (sigma * np.sqrt(2 * np.pi))
    elif Type == 1: #Approximation: THOMPSON et al (doi: 10.1107/S0021889887087090 )
        sigma = gau / (2 * np.sqrt(2 * np.log(2)))
        lb = lor / 2
        f = (sigma**5 + 2.69269 * sigma**4 * lb + 2.42843 * sigma**3 * lb**2 + 4.47163 * sigma**2 * lb**3 + 0.07842* sigma * lb**4 + lb**5) ** 0.2
        eta = 1.36603 * (lb/f) - 0.47719 * (lb/f)**2 + 0.11116 * (lb/f)**3
        lor = f / (np.pi * (axis**2 + f**2))
        gauss = np.exp( -axis**2 / (2 * f**2)) / (f * np.sqrt(2 * np.pi))
        return integral * (eta * lor + (1 - eta) * gauss)


def apodize(t,shift,sw,axLen,lor,gauss,cos2,hamming,wholeEcho = False):
    t2 = t - shift
    x = np.ones(axLen)
    if lor is not None:
        x = x * np.exp(-np.pi * lor * abs(t2))
    if gauss is not None:
        x = x * np.exp(-((np.pi * gauss * t2)**2) / (4 * np.log(2)))
    if cos2 is not None:
        x = x * (np.cos(cos2 * (-0.5 * shift * np.pi * sw / axLen + np.linspace(0, 0.5 * np.pi, axLen)))**2)
    if hamming is not None:
        alpha = 0.53836  # constant for hamming window
        x = x * (alpha + (1 - alpha) * np.cos(hamming * (-0.5 * shift * np.pi * sw / axLen + np.linspace(0, np.pi, axLen))))
    if wholeEcho:
        x[-1:-(int(len(x) / 2) + 1):-1] = x[:int(len(x) / 2)]

    return x

def fib(n):
    start = np.array([[1, 1], [1, 0]], dtype='int64')
    temp = start[:]
    for i in range(n):
        temp = np.dot(start, temp)
    return temp[0, 0], temp[0, 1], temp[1, 1]


def zcw_angles(m, symm=0):
    samples, fib_1, fib_2 = fib(m)
    js = np.arange(samples, dtype='Float64') / samples
    if symm == 0:
        # full
        c = (1., 2., 1.)
    elif symm == 1:
        # hemi
        c = (-1., 1., 1.)
    elif symm == 2:
        # oct
        c = (-1., 1., 4.)
    j_samples = fib_2 * js
    phi = 2 * np.pi / c[2] * np.mod(j_samples, 1.0)
    theta = np.arccos(c[0] * (c[1] * np.mod(js, 1.0) - 1))
    weight = np.ones(samples) / samples
    return phi, theta, weight

def euro(val, num):
    firstDigit = '%.0e' % val
    firstDigit = int(firstDigit[0])
    order = int(np.floor(np.log10(val)))
    if num < 1:
        return
    numStep = int(num) // 3 + 1
    if firstDigit == 1:
        subset = [1, 2, 5]
    elif firstDigit == 2:
        subset = [2, 5, 10]
    elif firstDigit == 5:
        subset = [5, 10, 20]
    else:
        return
    returnVal = np.tile(subset, numStep)
    orderArray = np.repeat(range(numStep), 3) + order
    returnVal = returnVal * 10.0**orderArray
    return returnVal[:num]
