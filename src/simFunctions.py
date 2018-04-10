#!/usr/bin/env python

# Copyright 2016 - 2018 Bas van Meerten and Wouter Franssen

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
import scipy.ndimage
import functions as func
import extendedCzjzek as eCzjz

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

def voigtLine(x, pos, lor, gau, integral, Type=0):
    lor = np.abs(lor)
    gau = np.abs(gau)
    axis = x - pos
    if Type == 0:  # Exact: Freq domain simulation via Faddeeva function
        if gau == 0.0:  # If no gauss, just take lorentz
            f = 1.0 / (np.pi * 0.5 * lor * (1 + (axis / (0.5 * lor))**2))
            return integral * f
        elif lor == 0.0: # If no lorentz, just take gauss
            sigma = gau / (2 * np.sqrt(2 * np.log(2)))
            f = np.exp(-axis**2/(2 * sigma**2)) / np.sqrt(2 * np.pi * sigma**2)
            return integral * f 
        else:
            sigma = gau / (2 * np.sqrt(2 * np.log(2)))
            z = (axis + 1j * lor / 2) / (sigma * np.sqrt(2))
            return integral * wofz(z).real / (sigma * np.sqrt(2 * np.pi))
    elif Type == 1:  # Approximation: THOMPSON et al (doi: 10.1107/S0021889887087090 )
        sigma = gau / (2 * np.sqrt(2 * np.log(2)))
        lb = lor / 2
        f = (sigma**5 + 2.69269 * sigma**4 * lb + 2.42843 * sigma**3 * lb**2 + 4.47163 * sigma**2 * lb**3 + 0.07842 * sigma * lb**4 + lb**5) ** 0.2
        eta = 1.36603 * (lb / f) - 0.47719 * (lb / f)**2 + 0.11116 * (lb / f)**3
        lor = f / (np.pi * (axis**2 + f**2))
        gauss = np.exp(-axis**2 / (2 * f**2)) / (f * np.sqrt(2 * np.pi))
        return integral * (eta * lor + (1 - eta) * gauss)

def makeSpectrum(x, sw, v, gauss, lor, weight):
    # Takes axis, frequencies and intensities and makes a spectrum with lorentz and gaussian broadening
    length = len(x)
    t = np.abs(np.fft.fftfreq(length, sw/float(length)))
    diff = (x[1]-x[0])*0.5
    final, junk = np.histogram(v, length, range=[x[0]-diff, x[-1]+diff], weights=weight)
    apod = np.exp(-np.pi * np.abs(lor) * t -((np.pi * np.abs(gauss) * t)**2) / (4 * np.log(2)))
    inten = np.fft.ifft(final) * apod
    inten = inten / sw * len(inten) 
    return inten

def makeMQMASSpectrum(x, sw, v, gauss, lor, weight, offset, shear):
    length1 = len(x[0])
    t1 = np.fft.fftfreq(length1, sw[0]/float(length1))
    diff1 = (x[0][1] - x[0][0])*0.5
    length2 = len(x[1])
    t2 = np.fft.fftfreq(length2, sw[1]/float(length2))
    diff2 = (x[1][1] - x[1][0])*0.5
    t1 = t1[:, np.newaxis]
    final, junk = np.histogram(v, length2, range=[x[1][0]-diff2, x[1][-1]+diff2], weights=weight)
    apod2 = np.exp(-np.pi * np.abs(lor[1] * t2)-((np.pi * np.abs(gauss[1]) * t2)**2) / (4 * np.log(2)))
    final = np.fft.fft(np.fft.ifft(final) * apod2) /sw[0] * length2/sw[1]
    if np.abs(shear) > 1e-7:
        offset += shear*x[1]
    apod1 = np.exp(2j * np.pi * offset * t1)
    apod1 *= np.exp(- np.pi * np.abs(lor[0] * t1) -((np.pi * np.abs(gauss[0]) * t1)**2) / (4 * np.log(2)))
    return final * apod1 # time domain in dim 0 and freq domain in dim 1

def csaAngleStuff(cheng):
    phi, theta, weight = zcw_angles(cheng, symm=2)
    sinT2 = np.sin(theta)**2
    return weight, [[sinT2 * np.cos(phi)**2, sinT2 * np.sin(phi)**2, np.cos(theta)**2]]

def tensorDeconvtensorFunc(x, t11, t22, t33, lor, gauss, multt, sw, weight):
    v = np.dot([t11, t22, t33], multt)
    return makeSpectrum(x, sw, v, gauss, lor, weight)

def tensorMASDeconvtensorFunc(x, pos, delta, eta, lor, gauss, sw, spinspeed, cheng, numssb):
    numssb = float(numssb)
    omegar = 2 * np.pi * 1e3 * spinspeed
    phi, theta, weight = zcw_angles(cheng, symm=2)
    sinPhi = np.sin(phi)
    cosPhi = np.cos(phi)
    sin2Theta = np.sin(2 * theta)
    cos2Theta = np.cos(2 * theta)
    tresolution = 2 * np.pi / omegar / numssb
    t = np.linspace(0, tresolution * (numssb - 1), numssb)
    cosOmegarT = np.cos(omegar * t)
    cos2OmegarT = np.cos(2 * omegar * t)
    angleStuff = [np.array([np.sqrt(2) / 3 * sinPhi * cosPhi * 3]).transpose() * cosOmegarT,
                  np.array([-1.0 / 3 * 3 / 2 * sinPhi**2]).transpose() * cos2OmegarT,
                  np.transpose([cos2Theta / 3.0]) * (np.array([np.sqrt(2) / 3 * sinPhi * cosPhi * 3]).transpose() * cosOmegarT),
                  np.array([1.0 / 3 / 2 * (1 + cosPhi**2) * cos2Theta]).transpose() * cos2OmegarT,
                  np.array([np.sqrt(2) / 3 * sinPhi * sin2Theta]).transpose() * np.sin(omegar * t),
                  np.array([cosPhi * sin2Theta / 3]).transpose() * np.sin(2 * omegar * t)]
    omegars = 2 * np.pi * delta * (angleStuff[0] + angleStuff[1] + eta * (angleStuff[2] + angleStuff[3] + angleStuff[4] + angleStuff[5]))
    numssb = angleStuff[0].shape[1]
    QTrs = np.concatenate([np.ones([angleStuff[0].shape[0], 1]), np.exp(-1j * np.cumsum(omegars, axis=1) * tresolution)[:, :-1]], 1)
    for j in range(1, numssb):
        QTrs[:, j] = np.exp(-1j * np.sum(omegars[:, 0:j] * tresolution, 1))
    rhoT0sr = np.conj(QTrs)
    # calculate the gamma-averaged FID over 1 rotor period for all crystallites
    favrs = np.zeros(numssb, dtype=complex)
    for j in range(numssb):
        favrs[j] += np.sum(weight * np.sum(rhoT0sr * np.roll(QTrs, -j, axis=1), 1) / numssb**2)
    # calculate the sideband intensities by doing an FT and pick the ones that are needed further
    inten = np.real(np.fft.fft(favrs))
    posList = np.array(np.fft.fftfreq(numssb, 1.0 / numssb)) * spinspeed * 1e3 + pos
    return makeSpectrum(x, sw, posList, gauss, lor, inten)
   
def quad1DeconvtensorFunc(x, I, pos, cq, eta, lor, gauss, angleStuff, freq, sw, weight):
    m = np.arange(-I, I)
    v = []
    cq *= 1e6
    weights = []
    for i in m:
        tmp = (cq / (4 * I * (2 * I - 1)) * (I * (I + 1) - 3 * (i + 1)**2)) - (cq / (4 * I * (2 * I - 1)) * (I * (I + 1) - 3 * (i)**2))
        v = np.append(v, tmp * (angleStuff[0] - eta * angleStuff[1]) + pos)
        weights = np.append(weights, weight)
    return makeSpectrum(x, sw, v, gauss, lor, weights) / (2 * I)


def quad1DeconvsetAngleStuff(cheng):
    phi, theta, weight = zcw_angles(cheng, symm=2)
    angleStuff = [0.5 * (3 * np.cos(theta)**2 - 1), 0.5 * np.cos(2 * phi) * (np.sin(theta)**2)]
    return weight, angleStuff


def quad1MASFunc(x, pos, cq, eta, lor, gauss, sw, spinspeed, cheng, I, numssb):
    numssb = float(numssb)
    omegar = 2 * np.pi * 1e3 * spinspeed
    phi, theta, weight = zcw_angles(cheng, symm=2)
    sinPhi = np.sin(phi)
    cosPhi = np.cos(phi)
    sin2Theta = np.sin(2 * theta)
    cos2Theta = np.cos(2 * theta)
    tresolution = 2 * np.pi / omegar / numssb
    t = np.linspace(0, tresolution * (numssb - 1), int(numssb))
    cosOmegarT = np.cos(omegar * t)
    cos2OmegarT = np.cos(2 * omegar * t)
    angleStuff = [np.array([np.sqrt(2) / 3 * sinPhi * cosPhi * 3]).transpose() * cosOmegarT,
                  np.array([-1.0 / 3 * 3 / 2 * sinPhi**2]).transpose() * cos2OmegarT,
                  np.transpose([cos2Theta / 3.0]) * (np.array([np.sqrt(2) / 3 * sinPhi * cosPhi * 3]).transpose() * cosOmegarT),
                  np.array([1.0 / 3 / 2 * (1 + cosPhi**2) * cos2Theta]).transpose() * cos2OmegarT,
                  np.array([np.sqrt(2) / 3 * sinPhi * sin2Theta]).transpose() * np.sin(omegar * t),
                  np.array([cosPhi * sin2Theta / 3]).transpose() * np.sin(2 * omegar * t)]
    m = np.arange(-I, 0)  # Only half the transitions have to be calculated, as the others are mirror images (sidebands inverted)
    eff = I**2 + I - m * (m + 1)  # The detection efficiencies of the top half transitions
    #Scale the intensities to sum to 1
    if np.floor(I) != I: #If not even:
        scale = np.sum(eff[0:-1] * 2) + eff[-1]
    else:
        scale = np.sum(eff) * 2
    eff = eff / scale
    splitting = np.arange(I - 0.5, -0.1, -1)  # The quadrupolar couplings of the top half transitions
    sidebands = np.zeros(int(numssb))
    for transition in range(len(eff)):  # For all transitions
        if splitting[transition] != 0:  # If quad coupling not zero: calculate sideban pattern
            delta = splitting[transition] * 2 * np.pi * 3 / (2 * I * (2 * I - 1)) * cq * 1e6  # Calc delta based on Cq [MHz] and spin quantum
            omegars = delta * (angleStuff[0] + angleStuff[1] + eta * (angleStuff[2] + angleStuff[3] + angleStuff[4] + angleStuff[5]))
            QTrs = np.concatenate([np.ones([angleStuff[0].shape[0], 1]), np.exp(-1j * np.cumsum(omegars, axis=1) * tresolution)[:, :-1]], 1)
            for j in range(1, int(numssb)):
                QTrs[:, j] = np.exp(-1j * np.sum(omegars[:, 0:j] * tresolution, 1))
            rhoT0sr = np.conj(QTrs)
            # calculate the gamma-averaged FID over 1 rotor period for all crystallites
            favrs = np.zeros(int(numssb), dtype=complex)
            for j in range(int(numssb)):
                favrs[j] += np.sum(weight * np.sum(rhoT0sr * np.roll(QTrs, -j, axis=1), 1) / numssb**2)
            # calculate the sideband intensities by doing an FT and pick the ones that are needed further
            partbands = np.real(np.fft.fft(favrs))
            sidebands += eff[transition] * (partbands + np.roll(np.flipud(partbands), 1))
        else:  # If zero: add all the intensity to the centreband
            sidebands[0] += eff[transition]
    posList = np.array(np.fft.fftfreq(int(numssb), 1.0 / numssb)) * spinspeed * 1e3 + pos
    return makeSpectrum(x, sw, posList, gauss, lor, sidebands)


def quad2StaticsetAngleStuff(cheng):
    phi, theta, weight = zcw_angles(cheng, symm=2)
    cosT2 = np.cos(theta)**2
    cosT4 = cosT2**2
    cos2P = np.cos(2 * phi)
    angleStuff = [-27 / 8.0 * cosT4 + 15 / 4.0 * cosT2 - 3 / 8.0,
                  (-9 / 4.0 * cosT4 + 2 * cosT2 + 1 / 4.0) * cos2P,
                  -1 / 2.0 * cosT2 + 1 / 3.0 + (-3 / 8.0 * cosT4 + 3 / 4.0 * cosT2 - 3 / 8.0) * cos2P**2]
    return weight, angleStuff


def quad2MASsetAngleStuff(cheng):
    phi, theta, weight = zcw_angles(cheng, symm=2)
    cosT2 = np.cos(theta)**2
    cosT4 = cosT2**2
    cos2P = np.cos(2 * phi)
    angleStuff = [21 / 16.0 * cosT4 - 9 / 8.0 * cosT2 + 5 / 16.0,
                  (-7 / 8.0 * cosT4 + cosT2 - 1 / 8.0) * cos2P,
                  1 / 12.0 * cosT2 + (+7 / 48.0 * cosT4 - 7 / 24.0 * cosT2 + 7 / 48.0) * cos2P**2]
    return weight, angleStuff


def quad2tensorFunc(x, I, pos, cq, eta, lor, gauss, angleStuff, freq, sw, weight):
    cq *= 1e6
    v = -1 / (6 * freq) * (3 * cq / (2 * I * (2 * I - 1)))**2 * (I * (I + 1) - 3.0 / 4) * (angleStuff[0] + angleStuff[1] * eta + angleStuff[2] * eta**2) + pos
    return makeSpectrum(x, sw, v, gauss, lor, weight)

def czjzekIntensities(sigma, d, wq, eta, wq0 = 0, eta0 = 0):
    #Calculates an intensity distribution for a Czjzek library
    #wq: omega_q grid (2D flattened)
    #eta: eta grid (2D flattened)
    if sigma == 0.0:  # protect against divide by zero
        czjzek = np.zeros_like(wq)
        czjzek[:, 0] = 1
    else:
        czjzek = wq**(d - 1) * eta / (np.sqrt(2 * np.pi) * sigma**d) * (1 - eta**2 / 9.0) * np.exp(-wq**2 / (2.0 * sigma**2) * (1 + eta**2 / 3.0))
    if np.sum(czjzek) == 0.0: #Protect against divide by zero
        czjzek = 0 * czjzek
    else:
        czjzek = czjzek / np.sum(czjzek)
    return czjzek

def quad2CzjzektensorFunc(x, sigma, d, pos, width, gauss, wq, eta, lib, freq, sw, wq0 = 0, eta0 = 0):
    sigma = sigma * 1e6
    if wq0 == 0 and eta0 == 0: #If normal Czjzek
        czjzek = czjzekIntensities(sigma, d, wq, eta, wq0, eta0)
    else:
        czjzek = eCzjz.getInts(sigma, d, eta0, wq0, wq, eta)
    fid = np.dot(czjzek, lib)
    length = len(x)
    t = np.fft.fftfreq(length, sw/float(length))
    pos -= x[int(len(x)/2)]
    apod = np.exp(2j*np.pi*pos*t-np.pi * np.abs(width) * np.abs(t) -((np.pi * np.abs(gauss) * t)**2) / (4 * np.log(2)))
    spectrum = np.real(np.fft.fft(fid * apod))
    spectrum = spectrum / sw * len(spectrum)
    return spectrum

def mqmasAngleStuff(cheng):
    phi, theta, weight = zcw_angles(cheng, symm=2)
    sinT2 = np.sin(theta)**2
    sinT4 = sinT2**2
    cos2P = np.cos(2 * phi)
    angleStuff = [9 * (35 * sinT4 - 40 * sinT2 + 8),
                  30 * cos2P * (6 * sinT2 - 7 * sinT4),
                  35 * cos2P**2 * sinT4 -20 * sinT2 + 4]
    return weight, angleStuff

def mqmasFreq(I, Cq, eta, angleStuff):
    v0Q = - Cq**2 * (3 + eta**2) / (40 * (I*(2*I-1))**2)
    v4Q = Cq**2 / (1120.0 * (2*I*(2*I-1))**2) * (angleStuff[0] + angleStuff[1]*eta + angleStuff[2]*eta**2)
    return v0Q, v4Q

def mqmasFunc(x, I, p, shear, scale, pos, cq, eta, lor, gauss, angleStuff, freq, sw, weight):
    cq *= 1e6
    v0Q, v4Q = mqmasFreq(I, cq, eta, angleStuff)
    C10 = I*(I+1) - 3/4.0
    C14 = -7/18.0 * (18*I*(I+1) - 17/2.0 - 5)
    v2 = pos + (C10 * v0Q + C14 * v4Q) / freq[-1]
    Cp0 = p * (I*(I+1) - 3/4.0 * p**2)
    Cp4 = -7/18.0 * p * (18*I*(I+1) - 17/2.0 * p**2 - 5)
    #v1 = p * pos + (Cp0 * v0Q + Cp4 * v4Q) / freq[-2]
    shearFactor = Cp4 / C14 * freq[-1] / freq[-2]
    offset = p * pos + Cp0 * v0Q / freq[-2] - shearFactor * (pos + C10 * v0Q / freq[-1])
    offset *= scale
    return makeMQMASSpectrum(x, sw, v2, gauss, lor, weight, offset, (shearFactor-shear))
