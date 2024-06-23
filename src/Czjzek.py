#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2016 - 2024 Bas van Meerten and Wouter Franssen

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

import math
import multiprocessing
import numpy as np
import scipy.integrate as SI
import scipy.special as SP
try: #If numba exists, use jit, otherwise make a mock decorator
    from numba import jit
except ImportError:
    def jit(func):
        return func


@jit
def gammaFunc(gamma, a11pa15, a51pre, a55part):
    """
    Calculates the gamma angle part of an extended Czjzek distribution. This
    Is the last integral of the calculation, so all terms are bundled here.

    Parameters
    ----------
    gamma: float
        The gamma angle in radians
    a11pa15: float
        The a11 term with part of the a15 mixed in
    a51pre: float
        Part of the a51 term
    a55part: float
        Part of the a55 term

    Returns
    -------
    float
        Exponent of all relevant terms and the gamma angle.
    """

    cos2G = math.cos(gamma)
    sin2G = math.sin(gamma)
    a51 = a51pre * cos2G
    a55 = a55part * sin2G
    return math.exp(a11pa15 + a51 + a55)

def alphaFunc(alpha, preCalc, a11, a15pre, a55pre, a51prepre1, a51prepre2):
    """
    Calculates the alpha angle part of an extended Czjzek distribution.

    Parameters
    ----------
    alpha: float
        The alpha angle in radians
    preCalc: float
        Some pre-calculated constants
    a11: float
        The a11 term
    a15pre:
        Precalc of the a15 term
    a55pre: float
        Precalc of the a55 term
    a51prepre1: float
        Precalc of a part of the a51 term
    a51prepre2: float
        Precalc of a part of the a51 term

    Returns
    -------
    float
        Integral over the gamma angle
    """
    cos2A = math.cos(alpha)
    sin2A = math.sin(alpha)
    a11pa15 = a11 - a15pre * cos2A
    a55part = sin2A * a55pre
    a51pre = a51prepre1 + a51prepre2 * cos2A
    Amp = SI.quad(gammaFunc, 0, 2 * np.pi, args=(a11pa15, a51pre, a55part), epsrel=0.1, epsabs=0)
    #Scale with 0.5, as integral is done to 2pi instead of pi (and sin(gamma) is used not sin(2 * gamma)
    return Amp[0] * 0.5

def betaFunc(beta, eta, eta0, preCalc):
    """
    Calculates the alpha angle part of an extended Czjzek distribution.

    Parameters
    ----------
    beta: float
        The beta angle in radians
    eta: float
        The eta value for which the Czjzek intensity needs to be found
    eta0: float
        The base eta value of the extended Czjzek distribution
    preCalc: float
        Some pre-calculated constants

    Returns
    -------
    float
        Integral over the alpha and gamma angle
    """
    cosB = math.cos(beta)
    cosBS = cosB**2
    sinB = math.sin(beta)
    sinBS = sinB**2
    preVal = preCalc[0] * sinB
    a11 = - 0.5 * (3 * cosBS - 1) * preCalc[1] * preCalc[2] + preCalc[3]
    a15pre = eta0 * preCalc[1] / 2 * sinBS * preCalc[2]
    a55pre = -cosB * preCalc[4] * preCalc[2]
    a51prepre1 = preCalc[1] / 2 * sinBS * eta * preCalc[2]
    a51prepre2 = 0.5 * (1 + cosBS) * preCalc[4] * preCalc[2]
    Amp = SI.quad(alphaFunc, 0, np.pi, args=(preCalc, a11, a15pre, a55pre, a51prepre1, a51prepre2), epsrel=0.1, epsabs=0)
    #Scale with 0.5, as integral is done to pi instead of pi/2 (and sin(alpha) is used not sin(2 * alpha)
    return Amp[0] * preVal * 0.5

def extendedCzjzek(inp):
    """
    Calculates the intensity of a extended Czjzek distribution for a particular Cq
    and eta value.

    Parameters
    ----------
    inp: list[floats]
        List with the inputs [cq, eta, cq0, eta0, sigma, d]

    Returns
    -------
    float
        Extended Czjzek intensity for the Cq and eta value of the input
    """
    cq, eta, cq0, eta0, sigma, d = inp
    #Prevent overflow error (if this condition is true, the point is far outside the distribution, and is 0 anyway)
    if abs(cq**2*(1+eta**2/3) - cq0**2)/(2 * sigma**2) > 1000:
        return 0.0
    #Prevent overflow error on eta range
    if cq0 / sigma * abs(eta0 - eta) > 10:
        return 0.0
    N1 = cq ** (d - 1) / (sigma ** d) * eta * (1 - eta**2 / 9)
    afterfact = -1.0 / (2 * sigma ** 2)
    preCalc = [N1,np.sqrt(3),  2.0/np.sqrt(3) * cq * cq0 * afterfact, (cq0**2 * (1 + eta0**2 / 3) + cq**2 * (1 + eta**2/3)) * afterfact]
    eta0deveta = eta0 / preCalc[1] * eta # eta0 divided by sqrt3 multiply by eta
    preCalc.append(eta0deveta)
    Amp = SI.quad(betaFunc, 0, np.pi, args=(eta, eta0, preCalc), epsrel=0.1, epsabs=0)
    return Amp[0]

@jit
def tFunc(t, pre2, fact, eta):
    """
    Function used for integrating over 't' in in extended Czjzek distribution
    with eta0 == 0.

    Parameters
    ----------
    t: float
        Value of t
    pre: float
        A pre-calculated constant (see the function 'extendedCzjzekNoEta0')
    pre2: float
        A pre-calculated constant (see the function 'extendedCzjzekNoEta0')
    fact: float
        A pre-calculated factor (see the function 'extendedCzjzekNoEta0')
    eta: float
        The eta value for which the Czjzek intensity is to be found

    Returns
    -------
    float
        Term that should be integrated over
    """
    expo = math.exp(fact * (3*t**2-1) + pre2)
    z = eta * abs(fact) * (1-t**2)
    bessel = SP.iv(0, z)
    return expo * bessel

def extendedCzjzekNoEta0(inp):
    """
    Function used to calculate the extended Czjzek distribution intensity for a
    particular value of Cq and eta. The function is optimized for the case
    eta0 == 0.

    Parameters
    ----------
    inp: list of floats
        Value of [cq, eta, cq0, sigma, d]

    Returns
    -------
    float
        Extended Czjzek intensity for the Cq and eta value of the input
    """
    cq, eta, cq0, sigma, d = inp
    #Prevent overflow error (if this condition is true, the point is far outside the distribution, and is 0 anyway)
    if abs(cq**2*(1+eta**2/3) - cq0**2)/(2 * sigma**2) > 1000:
        return 0.0
    #Prevent overflow error on eta range
    if cq0 / sigma * eta > 10:
        return 0.0
    pre = cq**(d-1) / sigma**d * eta * (1 - eta**2 / 9.0)
    pre2 = -(cq0**2 + cq**2 * (1 + eta**2 / 3.0)) / (2 * sigma**2)
    fact = cq * cq0 / (2*sigma**2)
    Amp = SI.quad(tFunc, 0, 1, args=(pre2, fact, eta), epsrel=0.0001, epsabs=0)[0]
    return Amp * pre

def normalCzjzekFunc(cq, eta, sigma, d):
    """
    Function used to calculate the normal Czjzek distribution intensity for a
    Cq and eta value.

    Parameters
    ----------
    cq: float
        Cq value of the quadrupole interaction
    eta: float
        eta value of the quadrupole interaction
    sigma: float
        Sigma value (i.e. width) of the distribution
    d: float
        D value of the distribution

    Returns
    -------
    float
        Normal Czjzek intensity for the Cq and eta value of the input
    """
    return cq**(d - 1) * eta / (np.sqrt(2 * np.pi) * sigma**d) * (1 - eta**2 / 9.0) * np.exp(-cq**2 / (2.0 * sigma**2) * (1 + eta**2 / 3.0))

def czjzekIntensities(sigma, d, cq, eta, cq0=0, eta0=0):
    """
    Function used to calculate the (extended) Czjzek distribution intensity for a
    Cq and eta grid. Based on the optional input of cq0 and eta0 values, it either
    calculates a normal Czjzek, or an extended Czjzek.

    Parameters
    ----------
    sigma: float
        Sigma value (i.e. width) of the distribution
    d: float
        D value of the distribution
    cq: ndarray
        A 1-D array with cq values
    eta: ndarray
        A 1-D array with eta values
    cq0: float, optional
        Base cq value of the distribution
    eta0: float, optional
        Base eta value of the distribution

    Returns
    -------
    ndarray
        1-D array of the normalized Czjzek intensity distribution
    """
    if cq0 == 0.0 and eta0 == 0.0:
        if sigma == 0.0:  # protect against divide by zero
            czjzek = np.zeros_like(cq)
        else:
            czjzek = normalCzjzekFunc(cq, eta, sigma, d)
    else:
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        if eta0 != 0.0:
            eta0 = 1 - abs(abs(eta0)%2 - 1) #scale continuously between 0--1
            fit = pool.map_async(extendedCzjzek, [(cq[i], eta[i], cq0, eta0, sigma, d) for i in range(len(cq))])
        else:
            fit = pool.map_async(extendedCzjzekNoEta0, [(cq[i], eta[i], cq0, sigma, d) for i in range(len(cq))])
        pool.close()
        pool.join()
        czjzek = np.array(fit.get())
    pos = np.isnan(czjzek)
    czjzek[pos] = 0.0 #Convert any issues to 0
    if np.sum(czjzek) == 0.0: #Protect against divide by zero
        czjzek = 0 * czjzek
    else:
        czjzek = czjzek / np.sum(czjzek)
    return czjzek
