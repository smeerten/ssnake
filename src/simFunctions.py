#!/usr/bin/env python3

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

import tempfile
import os
import shutil
import subprocess
import numpy as np
from safeEval import safeEval
import functions as func
import specIO as io
import Czjzek

class SimException(Exception):
    pass

def d2tens(b):
    """
    Calculates the Wigner (small) d-matrices of rank 2 for an array of beta angles.

    Parameters
    ----------
    b : ndarray
        A 1-D array with beta angles in radians.

    Returns
    -------
    ndarray
        A 3-D matrix with the wigner d-matrices as a function of the beta angles.
    """
    cosb = np.cos(b)
    cosb2 = np.cos(b/2.0)
    sinb2 = np.sin(b/2.0)
    # -2
    dm22 = sinb2**4
    # -1
    dm11 = 0.5 * sinb2**2 *(4 * cosb + 2)
    dm12 = -2 * sinb2**3 * cosb2
    # 0
    d00 = 0.5 * (3 * cosb**2 - 1)
    d01 = -np.sqrt(6) * sinb2 * cosb2 * cosb
    d02 = np.sqrt(6) * sinb2**2 * cosb2**2
    # 1
    d11 = 0.5 * cosb2**2 *(4 * cosb - 2)
    d12 = -2 * sinb2 * cosb2**3
    # 2
    d22 = cosb2**4
    d = np.array([[ d22,   d12,   d02,  dm12, dm22],
                  [-d12,   d11,   d01,  dm11, dm12],
                  [ d02,  -d01,   d00,  d01,  d02 ],
                  [-dm12,  dm11, -d01,  d11,  d12 ],
                  [ dm22, -dm12,  d02, -d12,  d22 ]])
    return np.rollaxis(d, 2, 0)

def D2tens(a, b, g):
    """
    Calculates the Wigner D-matrices of rank 2 for an array of alpha, beta, gamma angles.

    Parameters
    ----------
    a : ndarray
        A 1-D array with alpha angles in radians.
    b : ndarray
        A 1-D array with beta angles in radians.
        Should have the same length as a.
    g : ndarray
        A 1-D array with gamma angles in radians.
        Should have the same length as a.

    Returns
    -------
    ndarray
        A 3-D matrix with the wigner d-matrices as a function of the beta angles.
    """
    d = d2tens(b)
    m = np.arange(-2, 3)
    atmp = np.exp(1j * m[:, np.newaxis] * a[:, np.newaxis, np.newaxis])
    gtmp = np.exp(1j * m * g[:, np.newaxis, np.newaxis])
    return atmp * gtmp * d

def d4tens(b):
    """
    Calculates the Wigner (small) d-matrices of rank 4 for an array of beta angles.

    Parameters
    ----------
    b : ndarray
        A 1-D array with beta angles in radians.

    Returns
    -------
    ndarray
        A 3-D matrix with the wigner d-matrices as a function of the beta angles.
    """
    cosb = np.cos(b)
    sinb2 = np.sin(b/2.0)
    cosb2 = np.cos(b/2.0)
    # -4
    dm44 = sinb2**8
    # -3
    dm33 = 0.5 * sinb2**6 * (8 * cosb + 6)
    dm34 = -2 * np.sqrt(2) * sinb2**7 * cosb2
    # -2
    dm22 = sinb2**4 * (7 * (cosb - 1)**2 + 21 * (cosb - 1) + 15)
    dm23 = -0.5 * np.sqrt(7.0/2) * sinb2**5 * cosb2 * (8 * cosb + 4)
    dm24 = 2 * np.sqrt(7) * sinb2**6 * cosb2**2
    # -1
    dm11 = sinb2**2 * (7 * (cosb - 1)**3 + 105.0/4 * (cosb - 1)**2 + 30 * (cosb - 1) + 10)
    dm12 = -np.sqrt(2) * sinb2**3 * cosb2 * (7 * (cosb - 1)**2 + 35.0/2 * (cosb - 1) + 10)
    dm13 = 0.5 * np.sqrt(7) * sinb2**4 * cosb2**2 * (8 * cosb + 2)
    dm14 = -2 * np.sqrt(14) * sinb2**5 * cosb2**3
    # 0
    d00 = 1.0/8 * (35 * cosb**4 - 30 * cosb**2 + 3)
    d01 = -0.5 *np.sqrt(5) * sinb2 * cosb2 * (7 * (cosb - 1)**3 + 21 * (cosb - 1)**2 + 18 * (cosb - 1) + 4)
    d02 = np.sqrt(5.0/2) * sinb2**2 * cosb2**2 * (7 * (cosb - 1)**2 + 14 * (cosb - 1) + 6)
    d03 = -2 * np.sqrt(35) * sinb2**3 * cosb2**3 * cosb
    d04 = np.sqrt(70) * sinb2**4 * cosb2**4
    # 1
    d11 = cosb2**2 * (7 * (cosb - 1)**3 + 63.0/4 * (cosb - 1)**2 + 9 * (cosb - 1) + 1)
    d12 = -np.sqrt(2) * sinb2 * cosb2**3 * (7 * (cosb - 1)**2 + 21.0/2 * (cosb - 1) + 3)
    d13 = 0.5 * np.sqrt(7) * sinb2**2 * cosb2**4 * (8 * cosb - 2)
    d14 = -2 * np.sqrt(14) * sinb2**3 * cosb2**5
    # 2
    d22 = cosb2**4 * (7 * (cosb - 1)**2 + 7 * (cosb - 1) + 1)
    d23 = -0.5 * np.sqrt(7.0/2) * sinb2 * cosb2**5 * (8 * cosb - 4)
    d24 = 2 * np.sqrt(7) * sinb2**2 * cosb2**6
    # 3
    d33 = 0.5 * cosb2**6 * (8 * cosb - 6)
    d34 = -2 * np.sqrt(2) * sinb2 * cosb2**7
    # 4
    d44 = cosb2**8
    d = np.array([[d44,   d34,  d24,  d14,  d04, dm14, dm24, dm34, dm44],
                  [-d34,  d33,  d23,  d13,  d03, dm13, dm23, dm33, dm34],
                  [ d24, -d23,  d22,  d12,  d02, dm12, dm22, dm23, dm24],
                  [-d14,  d13, -d12,  d11,  d01, dm11, dm12, dm13, dm14],
                  [ d04, -d03,  d02, -d01,  d00, d01,  d02,  d03,  d04 ],
                  [-dm14, dm13,-dm12, dm11,-d01, d11,  d12,  d13,  d14 ],
                  [ dm24,-dm23, dm22,-dm12, d02,-d12,  d22,  d23,  d24 ],
                  [-dm34, dm33,-dm23, dm13,-d03, d13, -d23,  d33,  d34 ],
                  [ dm44,-dm34, dm24,-dm14, d04,-d14,  d24, -d34,  d44 ]])
    return np.rollaxis(d, 2, 0) # make the b values lie along the first dim

def D4tens(a, b, g):
    """
    Calculates the Wigner D-matrices of rank 4 for an array of alpha, beta, gamma angles.

    Parameters
    ----------
    a : ndarray
        A 1-D array with alpha angles in radians.
    b : ndarray
        A 1-D array with beta angles in radians.
        Should have the same length as a.
    g : ndarray
        A 1-D array with gamma angles in radians.
        Should have the same length as a.

    Returns
    -------
    ndarray
        A 3-D matrix with the wigner d-matrices as a function of the beta angles.
    """
    d = d4tens(b)
    m = np.arange(-4, 5)
    atmp = np.exp(1j * m[:, np.newaxis] * a[:, np.newaxis, np.newaxis])
    gtmp = np.exp(1j * m * g[:, np.newaxis, np.newaxis])
    return atmp * gtmp * d

def csaSpace(delta):
    """
    The space part of the CSA Hamiltonian defined in irreducible spherical tensor operator format.

    Parameters
    ----------
    delta : ndarray
        A 1-D array (of length 3) with the CSA tensor values (xyz definition) in Hz.

    Returns
    -------
    float
        The V00 element.
    ndarray
        The V2n elements ordered as [V2-2, V2-1, V20, V21, V22].
    """
    V00 = -np.sqrt(1.0/3.0) * (delta[0] + delta[1] + delta[2])
    V20 = np.sqrt(1.0/6.0) * (2*delta[2] - delta[0] - delta[1])
    V2pm2 = 0.5 * (delta[1] - delta[0])
    return V00, np.array([V2pm2, 0, V20, 0, V2pm2])

def csaSpin():
    """
    The spin part of the CSA Hamiltonian defined in irreducible spherical tensor operator format.
    Does not include the spin operator Iz and the magnetic field B0.

    Returns
    -------
    float
        The T00 element.
    float
        The T20 element.
    """
    T00 = -np.sqrt(1.0 / 3)
    T20 = np.sqrt(1.0 / 6.0) * 2
    return T00, T20

def firstQuadSpace(eta):
    """
    The space part of the first order quadrupole Hamiltonian defined in irreducible spherical tensor operator format.

    Parameters
    ----------
    eta : float
        The eta value of the quadrupole interaction.

    Returns
    -------
    ndarray
        The V2n elements ordered as [V2-2, V2-1, V20, V21, V22].
    """
    V20 = 1.0
    V2pm2 = np.sqrt(1 / 6.0) * eta
    return np.array([V2pm2, 0, V20, 0, V2pm2])

def firstQuadSpin(I, m1, m2):
    """
    The spin part of the first order quadrupole Hamiltonian defined in irreducible spherical tensor operator format.

    Parameters
    ----------
    I : float
        The spin quantum number.
    m1, m2: float
        The quantum numbers of the energy levels which are involved.

    Returns
    -------
    float
        The T20 element.
    """
    spin20_1 = 3 * m1**2 - I * (I + 1)
    spin20_2 = 3 * m2**2 - I * (I + 1)
    return spin20_1 - spin20_2

def secQuadSpace(eta):
    """
    The space part of the second order quadrupole Hamiltonian defined in irreducible spherical tensor operator format.

    Parameters
    ----------
    eta : float
        The eta value of the quadrupole interaction.

    Returns
    -------
    float
        The V00 element.
    ndarray
        The V2n elements ordered as [V2-2, V2-1, V20, V21, V22].
    ndarray
        The V4n elements ordered as [V4-4, V4-3, V4-2, V4-1, V40, V41, V42, V43, V44].
    """
    V00 = -1.0 / 5 * (3 + eta**2)
    V20 = 1.0 / 14 * (eta**2 - 3)
    V2pm2 = 1.0 / 7 * np.sqrt(3.0 / 2) * eta
    V40 = 1.0 / 140 * (18  + eta**2)
    V4pm2 = 3.0 / 70 * np.sqrt(5.0 / 2) * eta
    V4pm4 = 1 / (4 * np.sqrt(70)) * eta**2
    return V00, np.array([V2pm2, 0, V20, 0, V2pm2]), np.array([V4pm4, 0, V4pm2, 0, V40, 0, V4pm2, 0, V4pm4])

def secQuadSpin(I, m1, m2):
    """
    The spin part of the second order quadrupole Hamiltonian defined in irreducible spherical tensor operator format.

    Parameters
    ----------
    I : float
        The spin quantum number.
    m1, m2: float
        The quantum numbers of the energy levels which are involved.

    Returns
    -------
    float
        The T00 element.
    float
        The T20 element.
    float
        The T40 element.
    """
    spin00_1 = m1 * (I * (I + 1) - 3 * m1**2)
    spin20_1 = m1 * (8 * I * (I + 1) - 12 * m1**2 - 3)
    spin40_1 = m1 * (18 * I * (I + 1) - 34 * m1**2 - 5)
    spin00_2 = m2 * (I * (I + 1) - 3 * m2**2)
    spin20_2 = m2 * (8 * I * (I + 1) - 12 * m2**2 - 3)
    spin40_2 = m2 * (18 * I * (I + 1) - 34 * m2**2 - 5)
    return spin00_1 - spin00_2, spin20_1 - spin20_2, spin40_1 - spin40_2

def relaxationFunc(x, freq, sw, axMult, extra, amp, const, coeff, T):
    """
    Simulation function used for fitting relaxation curves.

    Parameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 1-D method only the last array in the list is used.
    freq : list of float
        The list of frequency per dimension in Hz (not used).
    sw : list of float
        The list of spectral width per dimension in Hz (not used).
    axMult : float
        The multiplier of the x-axis (not used).
    extra : list
        The extra parameters of the function (not used).
    amp : float
        The amplitude of the curve.
    const : float
        The constant.
    coeff : float
        The coefficient.
    T : float
        The relaxation time. Has the same units as x.

    Returns
    -------
    ndarray
        The relaxation curve. Has the same length as x[-1].
    """
    x = x[-1]
    return amp * (const + coeff * np.exp(-x / abs(T)))

def diffusionFunc(x, freq, sw, axMult, extra, amp, const, coeff, D):
    """
    Simulation function used for fitting diffusion curves.

    Parameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 1-D method only the last array in the list is used.
    freq : list of float
        The list of frequency per dimension in Hz (not used).
    sw : list of float
        The list of spectral width per dimension in Hz (not used).
    axMult : float
        The multiplier of the x-axis (not used).
    extra : list
        The extra parameters of the function [gamma, delta, triangle].
        Where gamma is the gyromagnetic ratio in MHz/T, delta is the duration of the gradient in seconds, and triangle (capital delta) is the time between the start of gradients in seconds.
    amp : float
        The amplitude of the curve.
    const : float
        The constant.
    coeff : float
        The coefficient.
    D : float
        The diffusion constant in m^2/s.

    Returns
    -------
    ndarray
        The diffusion curve. Has the same length as x[-1].
    """
    x = x[-1]
    gamma, delta, triangle = extra
    return amp * (const + coeff * np.exp(-(abs(gamma) *1e6 * 2 * np.pi * abs(delta) * x)**2 * abs(D) * (abs(triangle) - abs(delta) / 3.0)))

def functionRun(x, freq, sw, axMult, extra, *parameters):
    """
    Simulation function used for function fitting.
    The words between @ symbols are replaced by the fit values and the resulting string is evaluated.

    Parameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 1-D method only the last array in the list is used.
    freq : list of float
        The list of frequency per dimension in Hz (not used).
    sw : list of float
        The list of spectral width per dimension in Hz (not used).
    axMult : float
        The multiplier of the x-axis (not used).
    extra : list
        The extra parameters of the function [names, function].
        Where names is a list of strings with the names of the parameters and function is a string with the function for fitting.
    *parameters
        The parameters used in the fit. Should have the same length as names.

    Returns
    -------
    ndarray
        The curve resulting from the evaluation of the function.
    """
    names, function = extra
    x = x[-1]
    for i, elem in enumerate(names):
        function = function.replace('@' + elem + '@', str(parameters[i]))
    return safeEval(function, length=len(x), x=x)

def externalFitRunScript(x, freq, sw, axMult, extra, bgrnd, mult, *parameters):
    """
    Simulation function used for external fitting.
    The words between @ symbols are replaced by the fit values and the resulting string is used as an input script for the given command.
    The output of the command is processed (apodization, Fourier, regrid) as specified.

    Parameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 1-D method only the last array in the list is used.
    freq : list of float
        The list of frequency per dimension in Hz (not used).
    sw : list of float
        The list of spectral width per dimension in Hz (not used).
    axMult : float
        The multiplier of the x-axis (not used).
    extra : list
        The extra parameters of the function [names, command, script, output, spec].
        Where names is a list of strings with the names of the parameters, command is a string with the command for fitting, script is a string with the script to be modified, output is a list of two strings to which the stdout and sterr are written, and spec is a boolean which is True when the output should be a spectrum.
    bgrnd : float
        The offset value added to the output curve.
    mult : float
        The value by which the output curve is multiplied.
    *parameters
        The parameters used in the fit.
        The last three values are interpreted as [amp, lor, gauss] and are used as the amplitude, the Lorentzian apodization (Hz) and Gaussian apodization (chemical shift distribution) in axMult unit.
        Should have be three longer than names.

    Returns
    -------
    ndarray
        The curve result from the command.
    """
    names, command, script, output, spec = extra
    amp, lor, gauss = parameters[-3:]
    gauss /= axMult
    x = x[-1]
    if script is None:
        return None
    for i, elem in enumerate(names):
        script = script.replace('@' + elem + '@', str(parameters[i]))
    directory_name = tempfile.mkdtemp()
    inputFileName = "script.in"
    fullPath = os.path.join(directory_name, inputFileName)
    with open(fullPath, "w") as text_file:
        text_file.write(script)
    process = subprocess.Popen(command + ' ' + fullPath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=directory_name)
    if output:
        output[0], output[1] = process.communicate()
    else:
        process.wait()
    fileList = os.listdir(directory_name)
    fileList.remove(inputFileName)
    if not fileList:
        shutil.rmtree(directory_name, ignore_errors=True)
        return None
    outputFileName = fileList[0]
    masterData = io.autoLoad(os.path.join(directory_name, outputFileName))
    masterData.noUndo = True
    masterData.apodize(lor, gauss, [None, None], 0, 0, 0, 0, 0)
    if masterData.spec[0] != spec:
        masterData.complexFourier(0)
    masterData.regrid([x[0], x[-1]], len(x), 0)
    shutil.rmtree(directory_name, ignore_errors=True)
    return mult * amp * np.real(masterData.getHyperData(0))

def fib(n):
    """
    Calculates three Fibonacci numbers starting from a given point in the series.

    Parameters
    ----------
    n : int
        Start at the nth Fibonacci number.

    Returns
    -------
    int
        The n Fibonacci number.
    int
        The n+1 Fibonacci number.
    int
        The n+2 Fibonacci number.
    """
    start = np.array([[1, 1], [1, 0]], dtype=np.int64)
    temp = start[:]
    for i in range(n):
        temp = np.dot(start, temp)
    return temp[0, 0], temp[0, 1], temp[1, 1]

def zcw_angles(m, symm=0):
    """
    Calculates two angle sets for powder averaging based on the Zaremba, Conroy, and Wolfsberg (ZCW) method.
    The number of orientations depends on the given Cheng number.

    Parameters
    ----------
    m : int
        The Cheng number.
        The number of orientations is equal to the m+2 Fibonacci number.
    symm : {0, 1, 2}, optional
        The symmetry of the problem. When 0 the orientations run over the entire sphere, when 1 the orientations run over a hemispere, and when 2 the orientations run over and octant.

    Returns
    -------
    ndarray
        The phi angles.
    ndarray
        The theta angles.
    ndarray
        The weights of the different orientations.
    """
    samples, fib_1, fib_2 = fib(m)
    js = np.arange(samples, dtype=np.float64) / samples
    if symm == 0:
        c = (1., 2., 1.)
    elif symm == 1:
        c = (-1., 1., 1.)
    elif symm == 2:
        c = (-1., 1., 4.)
    j_samples = fib_2 * js
    phi = 2 * np.pi / c[2] * np.mod(j_samples, 1.0)
    theta = np.arccos(c[0] * (c[1] * np.mod(js, 1.0) - 1))
    weight = np.ones(samples) / samples
    return phi, theta, weight

def peakSim(x, freq, sw, axMult, extra, bgrnd, mult, pos, amp, lor, gauss):
    """
    Simulates an FID with Lorentzian and Gaussian broadening.

    Parameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 1-D method only the last array in the list is used.
    freq : list of float
        The list of frequency per dimension in Hz (not used).
    sw : list of float
        The list of spectral width per dimension in Hz.
        The last value is used to determine the dwell time.
    axMult : float
        The multiplier of the x-axis.
    extra : list
        The extra parameters of the function (not used).
    bgrnd : float
        The offset value added to the FID.
    mult : float
        The value by which the FID is multiplied.
    pos : float
        The frequency of the peak (in Hz*axMult).
    amp : float
        The amplitude of the peak.
    lor : float
        The Lorentzian broadening of the peak.
    gauss : float
        The Gaussian broadening of the peak (in Hz*axMult), corresponding to CS distribution.

    Returns
    -------
    ndarray
        The simulated FID
    """
    x = x[-1]
    pos /= axMult
    gauss /= axMult
    if pos < np.min(x) or pos > np.max(x):
        return np.zeros_like(x)
    lor = np.abs(lor)
    gauss = np.abs(gauss)
    length = len(x)
    t = np.fft.fftfreq(length, sw[-1]/float(length))
    return float(mult) * float(amp) / abs(sw[-1]) * np.exp(2j * np.pi * (pos - x[length//2]) * t - np.pi * np.abs(lor * t) - ((np.pi * np.abs(gauss) * t)**2) / (4 * np.log(2)))

def makeSpectrum(x, sw, v, gauss, lor, weight):
    """
    Creates an FID from a list of frequencies with corresponding weights.
    Also applies Lorentzian and Gaussian broadening.

    Parameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 1-D method only the last array in the list is used.
    sw : list of float
        The list of spectral width per dimension in Hz.
    v : ndarray
        The list of frequencies to be used in the FID.
    gauss : float
        Gaussian broadening in Hz.
    lor : float
        Lorentzian broadening in Hz.
    weight : ndarray
        The weights corresponding to the frequencies. Should have the same length as v.

    Returns
    -------
    ndarray
        The FID. Has the same length as x[-1].
    """
    length = len(x)
    t = np.abs(np.fft.fftfreq(length, sw / float(length)))
    diff = (x[1] - x[0]) * 0.5
    final, _ = np.histogram(v, length, range=[x[0]-diff, x[-1]+diff], weights=weight)
    apod = np.exp(-np.pi * np.abs(lor) * t - ((np.pi * np.abs(gauss) * t)**2) / (4 * np.log(2)))
    inten = np.fft.ifft(final) * apod
    inten *= len(inten)  / abs(sw)
    return inten

def makeMQMASSpectrum(x, sw, v, gauss, lor, weight, slope, fold_D1=True):
    """
    Creates an 2D FID from a list of frequencies with corresponding weights.
    Also applies Lorentzian and Gaussian broadening.

    Parameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 2-D method only the last two arrays in the list are used for the first and second dimension respectively.
    sw : list of float
        The list of spectral width per dimension in Hz.
        The second last value is used for the first dimension and the last value is used for the second dimension.
    v : list of ndarray
        The list of frequencies to be used in the 2D FID.
        The first ndarray contains in frequencies in the first (indirect) dimension and the second the frequencies in the second (direct) dimension.
        Both ndarray should have equal length.
    gauss : array_like
        Gaussian broadening in Hz for the first and second dimension.
    lor : array_like
        Lorentzian broadening in Hz for the first and second dimension.
    weight : ndarray
        The weights corresponding to the frequencies. Should have the same length as the arrays in v.
    slope : slope (t2/t1) along which the CS and CS distribution is refocused 
            to simulate CS distribution with gaussian g=broadening 
    fold_D1 : bool
        if fold_D1 is True, calculates a folded spectrum: any frequency that falls outside the 
        spectrum window is folded in.
        if fold_D1 is False, only retains frequencies that are within spetrum window.

    Returns
    -------
    ndarray
        The 2D FID. Has a shape of (len(x[-2]), len(x[-1])).
    """
    length_final1 = len(x[-2])
    step1 = x[-2][1] - x[-2][0]
    diff1 = step1/2
    sw1=abs(sw[-2])
    if fold_D1:
        sp_min = x[-2][0]  # current min freq
        sp_max = x[-2][-1] # current max freq
        if sp_min > sp_max: 
        # case where sw < 0 (sw was scaled with negative factor)
        # xax is reversed that is min is in x[-2][-1]
#            print("reversal!!!")
            sp_min, sp_max = sp_max, sp_min 
            v[0] *= -1  # frequencies used in histogram are reversed as well ?
            step1 *= -1 # make step > 0
            diff1 *= -1 # makes diff1 > 0
#        print(f"sp_min/max={sp_min}, {sp_max}  and sw1={sw1}, step1={step1}")
        maxf = np.max(v[0]) # max freq to reach
        minf = np.min(v[0])
#        print(f"f_min/max={minf}, {maxf}")
        # enlarge the histogram search to min/max frequencies in width multiple of initial sw1
        maxD1 = sp_max + (int(np.ceil((maxf-sp_max)/sw1)))*sw1 + diff1
        minD1 = sp_min - (int(np.ceil((sp_min-minf)/sw1)))*sw1 - diff1
        length1 = int(np.round((maxD1 - minD1 )/step1))
#        print(f"maxD1= {maxD1}, minD1={minD1}")
#        print(f"supposed steps number: {(maxD1 - minD1 )/step1}")
#        print(f"length1={length1}")
    else:
        sp_min = x[-2][0]  # current min freq
        sp_max = x[-2][-1] # current max freq
        if sp_min > sp_max:
#            print("reversal!!!")
            sp_min, sp_max = sp_max, sp_min 
            v[0] *= -1
            step1 *= -1
            diff1 *= -1
        maxD1 = sp_max+diff1
        minD1 = sp_min-diff1
        length1 = length_final1 

    split_len = int(np.round(length1/length_final1))
#    print(f"D1min/max={minD1}, {maxD1}    split_len={split_len}")
    t1 = np.fft.fftfreq(length_final1, sw[-2]/length_final1)
    t1 = t1[:, np.newaxis]

    length2 = len(x[-1])
    t2 = np.fft.fftfreq(length2, sw[-1]/length2)
    diff2 = (x[-1][1] - x[-1][0])*0.5
    minD2 = x[-1][0] - diff2
    maxD2 = x[-1][-1] + diff2
    if minD2 > maxD2:
        minD2, maxD2 = maxD2, minD2
        v[1] *= -1

    final, _, _ = np.histogram2d(v[0], v[1], [length1, length2], range=[[minD1, maxD1], [minD2, maxD2]], weights=weight)
#    print(f"histogram shape is {final.shape}")
    final = final.reshape(split_len, length_final1, length2)
#    print(f"histogram reshape is {final.shape}")
    final = final.sum(axis=0)
#    print(f"histogram final shape is {final.shape}")
    
    final = np.fft.ifftn(final)
    apod2 = np.exp(-np.pi * np.abs(lor[1] * t2) - 
                   ((np.pi * np.abs(gauss[1]) * (t2 + t1*slope))**2) / (4 * np.log(2)))
    apod1 = np.exp(-np.pi * np.abs(lor[0] * t1) - 
                   ((np.pi * np.abs(gauss[0]) * t1)**2) / (4 * np.log(2)))
    final *= apod1 * apod2 * length_final1 / sw1 * length2 / abs(sw[-1])
    return final

def carouselAveraging(spinspeed, v, weight, vConstant):
    """
    Performs carousel averaging for finite spinning samples.

    Parameters
    ----------
    spinspeed : float
        The spinning speed in Hz.
    v : 2-D ndarray
        The anisotropic part of the frequency.
        The first dimension contains the contributions of different alpha and beta angles.
        The second dimension contains a full rotation over gamma.
    weight : array_like
        Gaussian broadening in Hz for the first and second dimension.
    vConstant : float
        The offset frequency (isotropic value).

    Returns
    -------
    ndarray
        The 2D FID. Has a shape of (len(x[-2]), len(x[-1])).
    """
    numssb = v.shape[1]
    dt = 1.0 / spinspeed / numssb
    prod = np.exp(1j * np.cumsum(v * dt * 2 * np.pi, axis=1))
    tot = np.fft.fft(prod, axis=1)
    tot *= np.conj(tot)
    weight2 = weight[:, np.newaxis] / numssb**2
    tot *= weight2
    v = np.fft.fftfreq(numssb, 1.0 / numssb) * spinspeed
    return v + vConstant[:, np.newaxis], tot

def csaFreqBase(angle, tensor, D2, spinspeed, numssb):
    """
    Calculates the CSA frequencies for given alpha and beta angles and over a full circle over gamma.

    Parameters
    ----------
    angle : float
        The spinning angle in radians.
    tensor : array_like
        The three CSA tensor values in the xyz format (in Hz).
    D2 : array_like
        The second rank wigner rotation matrices used to calculate the CSA frequencies.
        The first dimension should contain the angular dependence.
    spinspeed : float
        The spinning frequency in Hz.
    numssb : int
        The number of alpha angles to calculate.

    Returns
    -------
    ndarray
        The array with frequencies. Has the same length as D2.
    float
        The isotropic frequency.
    """
    d2 = d2tens(np.array([angle]))[0, :, 2]
    A0, A2 = csaSpace(tensor)
    T0, T2 = csaSpin()
    dat0 = A0 * T0
    dat2 = A2 * T2
    dat2 = np.matmul(dat2, D2)
    factor2 = d2[2]
    vConstant = 0.0
    if spinspeed == np.inf:
        v = np.real(dat2[:, 2] * factor2 + dat0)
    elif spinspeed == 0.0:
        v = np.real(dat2[:, 2] + dat0)
    else:
        gammastep = 2 * np.pi / numssb
        gval = np.arange(numssb) * gammastep
        spinD2 = np.exp(1j * np.arange(-2, 3)[:, np.newaxis] * gval) * d2[:, np.newaxis]
        vConstant = np.real(dat0 + dat2[:, 2] * factor2) 
        dat2[:, 2] = 0
        v = np.matmul(dat2, spinD2)
    return v, vConstant

def csaFunc(x, freq, sw, axMult, extra, bgrnd, mult, spinspeed, t11, t22, t33, amp, lor, gauss):
    """
    Uses the quadCSAFunc function for the specific case where the quadrupole interaction is zero.
    """
    shiftdef, numssb, angle, D2, weight, MAStype = extra
    extra = [False, 0.5, numssb, angle, D2, None, weight, MAStype, shiftdef]
    return quadCSAFunc(x, freq, sw, axMult, extra, bgrnd, mult, spinspeed, t11, t22, t33, 0.0, 0.0, 0.0, 0.0, 0.0, amp, lor, gauss, 0)

def quadFunc(x, freq, sw, axMult, extra, bgrnd, mult, spinspeed, pos, cq, eta, amp, lor, gauss, lorST):
    """
    Uses the quadCSAFunc function for the specific case where the CSA interaction is zero.
    """
    satBool, I, numssb, angle, D2, D4, weight, MAStype = extra
    extra = [satBool, I, numssb, angle, D2, D4, weight, MAStype, 0]
    return quadCSAFunc(x, freq, sw, axMult, extra, bgrnd, mult, spinspeed, pos, pos, pos, cq, eta, 0.0, 0.0, 0.0, amp, lor, gauss, lorST)

def quadFreqBase(I, m1, m2, cq, eta, freq, angle, D2, D4, numssb, spinspeed):
    """
    Calculates the quadrupole frequencies for given alpha and beta angles and over a full circle over gamma.

    Parameters
    ----------
    I : float
        The spin quantum number.
    m1, m2 : float
        The levels between which to calculate the frequency.
    cq : float
        The Cq value of the quadrupole coupling in Hz.
    eta : float
        The asymmetry parameter of the quadrupole coupling.
    freq : float
        The Larmor frequency of the nucleus in Hz.
    angle : float
        The spinning angle in radians.
    D2 : array_like
        The second rank wigner rotation matrices used to calculate the quadrupole frequencies.
        The first dimension should contain the angular dependence.
    D4 : array_like
        The fourth rank wigner rotation matrices used to calculate the quadrupole frequencies.
        The first dimension should contain the angular dependence.
        Should have the same length as D2.
    numssb : int
        The number of alpha angles to calculate.
    spinspeed : float
        The spinning frequency in Hz.

    Returns
    -------
    ndarray
        The array with frequencies. Has the same length as D2 and D4.
    float
        The isotropic frequency.
    """
    if freq == 0.0:
        raise SimException("Sim: Frequency cannot be zero")
    pre2 = -cq**2 / (4 * I *(2 * I - 1))**2 * 2 / freq
    pre1 = cq / (4 * I *(2 * I - 1))
    firstA2 = pre1 * firstQuadSpace(eta)
    secA0, secA2, secA4 = secQuadSpace(eta)
    secA0 *= pre2
    secA2 *= pre2
    secA4 *= pre2
    d2 = d2tens(np.array([angle]))[0, :, 2]
    d4 = d4tens(np.array([angle]))[0, :, 4]
    factor2 = d2[2]
    factor4 = d4[4]
    firstspin2 = firstQuadSpin(I, m1, m2)
    secspin0, secspin2, secspin4 = secQuadSpin(I, m1, m2)
    dat0 = secA0 * secspin0
    dat2 = secA2 * secspin2 + firstA2 * firstspin2
    dat4 = secA4 * secspin4
    dat2 = np.matmul(dat2, D2)
    dat4 = np.matmul(dat4, D4)
    vConstant = 0
    if spinspeed == np.inf:
        v = np.real(dat4[:, 4] * factor4  + dat2[:, 2] * factor2 + dat0)
    elif spinspeed == 0.0:
        v = np.real(dat4[:, 4]  + dat2[:, 2] + dat0)
    else:
        gammastep = 2 * np.pi / numssb
        gval = np.arange(numssb) * gammastep
        spinD2 = np.exp(1j * np.arange(-2, 3)[:, np.newaxis] * gval) * d2[:, np.newaxis]
        spinD4 = np.exp(1j * np.arange(-4, 5)[:, np.newaxis] * gval) * d4[:, np.newaxis]
        vConstant = np.real(dat0 + dat2[:, 2] * factor2 + dat4[:, 4] * factor4)
        dat4[:, 4] = 0
        dat2[:, 2] = 0
        v = np.matmul(dat2, spinD2) + np.matmul(dat4, spinD4)
    return v, vConstant

def quadCSAFunc(x, freq, sw, axMult, extra, bgrnd, mult, spinspeed, t11, t22, t33, cq, eta, alphaCSA, betaCSA, gammaCSA, amp, lor, gauss, lorST):
    """
    Calculates an FID of a powder averaged site under influence of CSA and a quadrupole interaction.
    This function works for static, finite, and infinite spinnning.

    Parameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 1-D method only the last array in the list is used.
    freq : list of float
        The list of frequency per dimension in Hz. Only the last value is used.
    sw : list of float
        The list of spectral width per dimension in Hz. Only the last value is used.
    axMult : float
        The multiplier of the x-axis.
    extra : list
        The extra parameters defined as [satBool, I, numssb, angle, D2, D4, weight, MAStype, shiftdef].
        satBool is a boolean when True will include the quadrupole satellites.
        I is the spin quantum number.
        numssb is the number of sidebands to be simulated.
        angle is the spinning angle in radians.
        D2 is the second rank wigner rotation matrices.
        D4 is the fourth rank wigner rotation matrices.
        weight are the weights corresponding to the orientations of D2 and D4.
        MAStype=0 performs a static simulation, MAStype=1 performs a finite spinning simulation, and MAStype=2 performs an infinite spinning simulation.
        shiftdef is the definition in which t11, t22, and t33 are given (see shiftConversion).
    bgrnd : float
        The offset value added to the FID.
    mult : float
        The value by which the FID is multiplied.
    spinspeed : float
        The spinning speed in kHz.
    t11, t22, t33 : float
        The tensor values given in the definition specified by shiftdef.
    cq : float
        The quadrupole coupling constant Cq given in MHz.
    eta : float
        The asymmetry parameter of the quadrupole coupling.
    alphaC, betaC, alphaC : float
        The angles of the relative orientation of the CSA tensor with respect to the quadrupole tensor given in degrees.
    amp : float
        The amplitude of the peak.
    lor : float
        The Lorentzian broadening of the peak.
    gauss : float
        The Gaussian broadening of the peak.
    lorST : float
        The Lorentzian broadening of STs of the peak.

    Returns
    -------
    ndarray
        The simulated FID
    """
    alphaCSA *= np.pi / 180.0   # Degrees to radians
    betaCSA *= np.pi / 180.0    # Degrees to radians
    gammaCSA *= np.pi / 180.0   # Degrees to radians
    x = x[-1]
    satBool, I, numssb, angle, D2, D4, weight, MAStype, shiftdef = extra
    if MAStype == 0:
        spinspeed = 0.0
    elif MAStype == 2:
        spinspeed = np.inf
    if not satBool and (I % 1) == 0.0:
        return np.zeros_like(x)          # Integer spins have no central transition
    if shiftdef == 2:                    # If heaberlen, make eta continuous, and between 0--1
        t33 = 1 - abs(abs(t33) % 2 - 1)
    elif shiftdef == 3:                  # For Hertzfeld-Berger
        t33 = 1 - abs(abs(t33 + 1)%4 - 2)
    tensor = np.array(func.shiftConversion([t11, t22, t33], shiftdef)[1], dtype=float)
    tensor /= float(axMult)
    gauss /= axMult
    cq *= 1e6
    eta = 1 - abs(abs(eta) % 2 - 1)      # Force eta to 0--1 in a continuous way: 0.9 == 1.1, 0 == 2
    spinspeed *= 1e3
    freq = freq[-1]
    sw = sw[-1]
    mList = np.arange(-I, I)
    totalEff = len(mList) * (I**2 + I) - np.sum(mList * (mList + 1))
    if not satBool:
        mList = [-0.5]
    spectrum = np.zeros(len(x), dtype=complex)
    relativeD2 = D2tens(np.array([alphaCSA]), np.array([betaCSA]), np.array([gammaCSA]))
    vCSA, vConstantCSA = csaFreqBase(angle, tensor, np.matmul(relativeD2, D2), spinspeed, numssb)
    for m in mList:
        eff = I**2 + I - m * (m + 1)
        eff /= totalEff
        if I == 0.5:
            v = vCSA
            vConstant = vConstantCSA
        else:
            v, vConstant = quadFreqBase(I, m, m+1, cq, eta, freq, angle, D2, D4, numssb, spinspeed)
            v += vCSA
            vConstant += vConstantCSA
        tot = weight
        if spinspeed not in (0.0, np.inf):
            v, tot = carouselAveraging(spinspeed, v, weight, vConstant)
        if m == -0.5:
            lb = lor
        else:
            lb = lorST
        spectrum += eff * makeSpectrum(x, sw, v, gauss, lb, tot)
    return mult * amp * spectrum

def quadCzjzekFunc(x, freq, sw, axMult, extra, bgrnd, mult, pos, sigma, cq0, eta0, amp, lor, gauss):
    """
    Calculates an FID of a quadrupole spectrum with an (extended) Czjzek distribution using a library of spectra.

    Parameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 1-D method only the last array in the list is used.
    freq : list of float
        The list of frequency per dimension (not used).
    sw : list of float
        The list of spectral width per dimension in Hz. Only the last value is used.
    axMult : float
        The multiplier of the x-axis.
    extra : list
        The extra parameters defined as [method, d, lib, cq, eta].
        method is a boolean when True the extended version of the Czjzek distribution is used.
        d is the order parameter of the distribution.
        lib is the library (array) of FIDs.
        cq (in Hz) and eta are the quadrupole parameters corresponding to the FIDs in lib and both should have the same length as lib.
    bgrnd : float
        The offset value added to the FID.
    mult : float
        The value by which the FID is multiplied.
    pos : float
        The isotropic chemical shift in the units defined by axMult.
    sigma : float
        The width of the Czjzek distribution in MHz.
    cq0 : float
        The extended Czjzek quadrupole coupling constant Cq0 given in MHz (only used when method=True).
    eta0 : float
        The extended Czjzek asymmetry parameter of the quadrupole coupling (only used when method=True).
    amp : float
        The amplitude of the peak.
    lor : float
        The Lorentzian broadening of the peak.
    gauss : float
        The Gaussian broadening of the peak in the units defined by axMult.

    Returns
    -------
    ndarray
        The simulated FID
    """
    x = x[-1]
    sw = sw[-1]
    method, d, lib, cq, eta = extra
    if method == 0:
        cq0 = 0
        eta0 = 0
    pos /= axMult
    gauss /= axMult
    sigma = abs(sigma) * 1e6
    cq0 *= 1e6
    czjzek = Czjzek.czjzekIntensities(sigma, d, cq, eta, cq0, eta0)
    fid = np.dot(czjzek, lib)
    length = len(x)
    t = np.fft.fftfreq(length, sw/float(length))
    pos -= x[len(x)//2]
    apod = np.exp(2j * np.pi * pos * t - np.pi * np.abs(lor) * np.abs(t) - ((np.pi * np.abs(gauss) * t)**2) / (4 * np.log(2)))
    return mult * amp * fid * apod

def mqmasFunc(x, freq, sw, axMult, extra, bgrnd, mult, spinspeed, pos, sigmaCS, cq, eta, amp, lor2, lor1, gauss2=0, gauss1=0):
    """
    Calculates a 2-D FID of an MQMAS spectrum.
    This function works for static, finite, and infinite spinnning.

    Parameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 2-D method only the last two arrays in the list are used for the first and second dimension respectively.
    freq : list of float
        The list of frequency per dimension in Hz.
        The second last value is used for the first dimension and the last value is used for the second dimension.
    sw : list of float
        The list of spectral width per dimension in Hz.
        The second last value is used for the first dimension and the last value is used for the second dimension.
    axMult : float
        The multiplier of the x-axis.
    extra : list
        The extra parameters defined as [I, mq, numssb, angle, D2, D4, weight, shear, scale, MAStype].
        I is the spin quantum number.
        mq is the multiple quantum transition that is used in the indirect dimension.
        numssb is the number of sidebands to be simulated.
        angle is the spinning angle in radians.
        D2 is the second rank wigner rotation matrices.
        D4 is the fourth rank wigner rotation matrices.
        weight are the weights corresponding to the orientations of D2 and D4.
        shear is the shearing factor.
        scale is the scaling factor of the indirect axis.
        MAStype=0 performs a static simulation, MAStype=1 performs a finite spinning simulation, and MAStype=2 performs an infinite spinning simulation.
        foldF1 (bool): calculates a folded spectrum in F1 if True
    bgrnd : float
        The offset value added to the FID.
    mult : float
        The value by which the FID is multiplied.
    spinspeed : float
        The spinning speed in kHz.
    pos : float
        The isotropic chemical shift value defined in terms of axMult.
    sigmaCS : float
        The Gaussian broadening along the chemical shift axis in units defined by axMult.
    cq : float
        The quadrupole coupling constant Cq given in MHz.
    eta : float
        The asymmetry parameter of the quadrupole coupling.
    amp : float
        The amplitude of the peak.
    lor2 : float
        The Lorentzian broadening in the direct dimension.
    gauss2 : float
        The Gaussian broadening in the direct dimension defined (not used).
    lor1 : float
        The Lorentzian broadening in the indirect dimension.
    gauss1 : float
        The Gaussian broadening in the indirect dimension (not used).

    Returns
    -------
    ndarray
        The simulated 2-D FID
    """
    freq1 = freq[-2]
    freq2 = freq[-1]
    I, mq, numssb, angle, D2, D4, weight, shear, scale, MAStype, foldF1 = extra
    if MAStype == 0:
        spinspeed = 0.0
    elif MAStype == 2:
        spinspeed = np.inf
    pos /= axMult
    sigmaCS /= axMult
    spinspeed *= 1e3
    cq *= 1e6
    eta = 1 - abs(abs(eta)%2 - 1)
    tot = weight
    v2, vConstant2 = quadFreqBase(I, -0.5, 0.5, cq, eta, freq2, angle, D2, D4, numssb, spinspeed)
    v1, vConstant1 = quadFreqBase(I, -0.5*mq, 0.5*mq, cq, eta, freq2, angle, D2, D4, numssb, spinspeed)
    if spinspeed not in (0.0, np.inf):
        v2, tot2 = carouselAveraging(spinspeed, v2, weight, vConstant2)
        v1, tot1 = carouselAveraging(spinspeed, v1, np.ones_like(weight), vConstant1)
        tot = tot1*tot2
    v2 += pos
    v1 += mq*pos - v2 * shear
    v1 *= scale
    slope = (mq-shear)*scale  # t1/t2 slope along which to apply CS gaussian distribution
    return mult * amp * makeMQMASSpectrum(x, sw, [np.real(v1.flatten()), np.real(v2.flatten())], [0, sigmaCS], [lor1, lor2], np.real(tot).flatten(), slope, foldF1)

def mqmasCzjzekFunc(x, freq, sw, axMult, extra, bgrnd, mult, pos, sigmaCS, sigma, cq0, eta0, amp, lor2, lor1, gauss2=0, gauss1=0):
    """
    Calculates a 2-D FID of an MQMAS spectrum with an (extended) Czjzek distribution using a library of 1-D spectra.

    , shear_factorParameters
    ----------
    x : list of ndarray
        A list of axes values for the simulation.
        As this is a 1-D method only the last array in the list is used.
    freq : list of float
        The list of frequency per dimension in Hz.
        The second last value is used for the first dimension and the last value is used for the second dimension.
    sw : list of float
        The list of spectral width per dimension in Hz.
        The second last value is used for the first dimension and the last value is used for the second dimension.
    axMult : float
        The multiplier of the x-axis.
    extra : list
        The extra parameters defined as [I, mq, cq, eta, lib, shear, scale, method, d].
        I is the spin quantum number.
        mq is the multiple quantum transition that is used in the indirect dimension.
        cq (in Hz) and eta are the quadrupole parameters corresponding to the FIDs in lib and both should have the same length as lib.
        lib is the library (array) of FIDs.
        shear is the shearing factor.
        scale is the scaling factor of the indirect axis.
        method is a boolean when True the extended version of the Czjzek distribution is used.
        d is the order parameter of the distribution.
    bgrnd : float
        The offset value added to the FID.
    mult : float
        The value by which the FID is multiplied.
    pos : float
        The isotropic chemical shift in the units defined by axMult.
    sigma : float
        The width of the Czjzek distribution in MHz.
    sigmaCS : float
        The Gaussian broadening along the chemical shift axis in units defined by axMult.
    cq0 : float
        The extended Czjzek quadrupole coupling constant Cq0 given in MHz (only used when method=True).
    eta0 : float
        The extended Czjzek asymmetry parameter of the quadrupole coupling (only used when method=True).
    amp : float
        The amplitude of the peak.
    lor2 : float
        The Lorentzian broadening in the direct dimension.
    gauss2 : float
        The Gaussian broadening in the direct dimension (not used).
    lor1 : float
        The Lorentzian broadening in the indirect dimension.
    gauss1 : float
        The Gaussian broadening in the indirect dimension (not used).

    Returns
    -------
    ndarray
        The simulated 2-D FID
    """
    if freq[-1] == 0.0 or freq[-2] == 0.0:
        raise SimException("Sim: Frequency cannot be zero")
    I, mq, cq, eta, lib, shear, scale, method, d = extra
    if method == 1:
        cq0 *= 1e6
    else:
        cq0 = 0
        eta0 = 0
    pos /= axMult
    sigmaCS /= axMult
    sigma *= 1e6
    czjzek = Czjzek.czjzekIntensities(sigma, d, cq, eta, cq0, eta0)
    length2 = len(x[-1])
    czjzek *= length2 / abs(sw[-2])
    newLib = czjzek[..., np.newaxis]*lib
    length1 = len(x[-2])
    t1 = np.fft.fftfreq(length1, sw[-2]/float(length1))
    t1 = t1[:, np.newaxis]
    t2 = np.fft.fftfreq(length2, sw[-1]/float(length2))
    V40 = 1.0 / 140 * (18  + eta**2)
    T40_m = mq * (18 * I * (I + 1) - 34 * (mq/2.0)**2 - 5)
    T40_1 = (18 * I * (I + 1) - 34 * (0.5)**2 - 5)
    pre = -cq**2 / (4 * I *(2 * I - 1))**2 * 2 / freq[-2]
    cosb = np.cos(np.arctan(np.sqrt(2)))
    factor = 1.0/8 * (35 * cosb**4 - 30 * cosb**2 + 3)
    offset = pre * V40 * T40_m * factor
    shearFactor = T40_m * freq[-1] / (T40_1 * freq[-2])
    offset *= scale
    ind = np.digitize(offset, x[-2]-(x[-2][1]-x[-2][0])/2.0)
    fid = np.zeros((length1, length2), dtype=complex)
    for i, _ in enumerate(ind):
        fid[ind[i]-1] += newLib[i]
    fid = np.fft.ifft(fid, axis=0)
    posIndirect = pos * (mq - shearFactor) * scale
    offsetMat = np.exp(2j * np.pi * (posIndirect * t1 + (pos - x[-1][length2//2])*t2))
    shiftGauss = np.exp(-((np.pi * np.abs(sigmaCS) * (t2 + t1*(mq-shearFactor)*scale))**2) / (4 * np.log(2)))
#    apod2 = np.exp(-np.pi * np.abs(lor2 * t2) - ((np.pi * np.abs(gauss2) * t2)**2) / (4 * np.log(2)))
#    apod1 = np.exp(-np.pi * np.abs(lor1 * t1) - ((np.pi * np.abs(gauss1) * t1)**2) / (4 * np.log(2)))
    apod2 = np.exp(-np.pi * np.abs(lor2 * t2) )
    apod1 = np.exp(-np.pi * np.abs(lor1 * t1) )
    fid *= offsetMat * apod1 * apod2 * shiftGauss
    shearMat = np.exp((shearFactor-shear) * 2j * np.pi * t1 * x[-1])
    fid = np.fft.fft(fid, axis=1) * shearMat
    return mult * amp * fid * length1 / length2

def genLib(length, minCq, maxCq, minEta, maxEta, numCq, numEta, extra, freq, sw, spinspeed):
    """
    Generate a library of FIDs for Czjzek distribution fitting.

    Parameters
    ----------
    length : int
        The length of the FIDs to generate.
    minCq, maxCq, minEta, maxEta : float
        The Cq (in MHz) and eta values between which to generate the FIDs.
    numCq, numEta : int
        The number of FIDs to generate along Cq and eta respectively.
    extra : list
        The extra parameters as used by quadFunc.
    freq : float
        The Larmor frequency of the nucleus in Hz.
    sw : float
        The spectral width in Hz.
    spinspeed : float
        The spinning frequency in Hz.

    Returns
    -------
    ndarray
        The library of FIDs with shape (numCq*numEta, length).
    ndarray
        The Cq values corresponding to the FIDs in library (in Hz).
    ndarray
        The eta values corresponding to the FIDs in library.
    """
    cq, eta = np.meshgrid(np.linspace(minCq, maxCq, numCq), np.linspace(minEta, maxEta, numEta))
    cq = cq.flatten()
    eta = eta.flatten()
    x = np.fft.fftshift(np.fft.fftfreq(length, 1/float(sw)))
    lib = np.zeros((len(cq), length), dtype=complex)
    for i, (cqi, etai) in enumerate(zip(cq, eta)):
        lib[i] = quadFunc([x], [freq], [sw], 1.0, extra, 0.0, 1.0, spinspeed, 0.0, cqi, etai, 1.0, 0.0, 0.0, 0.0)
    return lib, cq*1e6, eta
