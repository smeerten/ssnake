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
import tempfile
import os
import shutil
import subprocess
from safeEval import safeEval
import functions as func
import specIO as io
import Czjzek as Czjzek

def d2tens(b):
    cosb = np.cos(b)
    sinb = np.sin(b)
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
    d = np.array([[d22, d12, d02, dm12, dm22],
                  [-d12, d11, d01, dm11, dm12],
                  [d02, -d01, d00, d01, d02],
                  [-dm12, dm11, -d01, d11, d12],
                  [dm22, -dm12, d02, -d12, d22]])
    return np.rollaxis(d, 2, 0)

def D2tens(a, b, g):
    d = d2tens(b)
    m = np.arange(-2, 3)
    atmp = np.exp(1j * m[:,np.newaxis] * a[:,np.newaxis,np.newaxis])
    gtmp = np.exp(1j * m * g[:,np.newaxis,np.newaxis])
    return atmp * gtmp * d

def d4tens(b):
    cosb = np.cos(b)
    sinb = np.sin(b)
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
    d = np.array([[d44,  d34, d24, d14, d04, dm14, dm24, dm34, dm44],
                  [-d34, d33, d23, d13, d03, dm13, dm23, dm33, dm34],
                  [d24, -d23, d22, d12, d02, dm12, dm22, dm23, dm24],
                  [-d14, d13, -d12, d11,  d01, dm11, dm12, dm13, dm14],
                  [d04, -d03, d02, -d01, d00, d01, d02, d03, d04],
                  [-dm14, dm13, -dm12, dm11, -d01, d11, d12, d13, d14],
                  [dm24, -dm23, dm22, -dm12, d02, -d12, d22, d23, d24],
                  [-dm34, dm33, -dm23, dm13, -d03, d13, -d23, d33, d34],
                  [dm44, -dm34, dm24, -dm14, d04, -d14, d24, -d34, d44]])
    return np.rollaxis(d, 2, 0) # make the b values lie along the first dim

def D4tens(a, b, g):
    d = d4tens(b)
    m = np.arange(-4, 5)
    atmp = np.exp(1j * m[:,np.newaxis] * a[:,np.newaxis,np.newaxis])
    gtmp = np.exp(1j * m * g[:,np.newaxis,np.newaxis])
    return atmp * gtmp * d

def csaSpace(delta): # CSA spherical tensor cartesian space functions
    # delta in Hz
    V00 = - np.sqrt(1.0/3.0) * (delta[0] + delta[1] + delta[2])
    V20 =   np.sqrt(1.0/6.0) * (2*delta[2] - delta[0] - delta[1])
    V2pm2 = 0.5 * (delta[0] - delta[1])
    return V00, np.array([V2pm2, 0, V20, 0, V2pm2])

def csaSpin(): # CSA spherical tensor spin space functions
    T00 = -np.sqrt(1.0 / 3)  # times Iz and B0, but is 1
    T20 = np.sqrt(1.0 / 6.0) * 2  # times Iz and B0, but is 1
    return T00, T20

def firstQuadSpace(eta):
    V20 = 1.0 
    V2pm2 = np.sqrt(1 / 6.0) * eta
    return np.array([V2pm2, 0, V20, 0, V2pm2])

def firstQuadSpin(I, m1, m2): #T values
    spin20_1 = 3 * m1**2 - I * (I + 1)
    spin20_2 = 3 * m2**2 - I * (I + 1)
    return spin20_1 - spin20_2

def secQuadSpace(eta):
    V00 = -1.0 / 5 * (3 + eta**2)
    V20 = 1.0 / 14 * (eta**2 - 3)
    V2pm2 = 1.0 / 7 * np.sqrt(3.0 / 2) * eta
    V40 = 1.0 / 140 * (18  + eta**2)
    V4pm2 = 3.0 / 70 * np.sqrt(5.0 / 2) * eta
    V4pm4 = 1 / (4 * np.sqrt(70)) * eta**2
    return V00, np.array([V2pm2, 0, V20, 0, V2pm2]), np.array([V4pm4, 0, V4pm2, 0, V40, 0, V4pm2, 0, V4pm4])

def secQuadSpin(I, m1, m2): #T values
    spin00_1 = m1 * (I * (I + 1) - 3 * m1**2)
    spin20_1 = m1 * (8 * I * (I + 1) - 12 * m1**2 - 3)
    spin40_1 = m1 * (18 * I * (I + 1) - 34 * m1**2 - 5)
    spin00_2 = m2 * (I * (I + 1) - 3 * m2**2)
    spin20_2 = m2 * (8 * I * (I + 1) - 12 * m2**2 - 3)
    spin40_2 = m2 * (18 * I * (I + 1) - 34 * m2**2 - 5)
    return spin00_1 - spin00_2, spin20_1 - spin20_2, spin40_1 - spin40_2

def relaxationFunc(x, freq, sw, axMult, extra, amp, const, coeff, T):
    x = x[-1]
    return amp * (const + coeff * np.exp(-x / T))

def diffusionFunc(x, freq, sw, axMult, extra, amp, const, coeff, D):
    x = x[-1]
    gamma, delta, triangle = extra
    return amp * (const + coeff * np.exp(-(gamma * delta * x)**2 * D * (triangle - delta / 3.0)))

def functionRun(x, freq, sw, axMult, extra, *parameters):
    names, function = extra
    x = x[-1]
    for i, elem in enumerate(names):
        function = function.replace('@' + elem, str(parameters[i]))
    return safeEval(function, length=len(x), x=x)

def SIMPSONRunScript(x, freq, sw, axMult, extra, bgrnd, *parameters):
    names, command, script, output, spec = extra
    amp, lor, gauss = parameters[-3:]
    x = x[-1]
    if script is None:
        return None
    for i, elem in enumerate(names):
        script = script.replace('@' + elem, str(parameters[i]))
    directory_name = tempfile.mkdtemp()
    inputFileName = "simpsonScript.in"
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
    masterData.apodize(lor, gauss, 0, 0, 0, 0, 0, 0)
    if masterData.spec[0] != spec:
        masterData.fourier(0)
    masterData.regrid([x[0], x[-1]], len(x), 0)
    shutil.rmtree(directory_name, ignore_errors=True)
    return amp * np.real(masterData.getHyperData(0))

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
    
def peakSim(x, freq, sw, axMult, extra, bgrnd, pos, amp, lor, gauss):
    x = x[-1]
    pos /= axMult
    lor = np.abs(lor)
    gauss = np.abs(gauss)
    length = len(x)
    t = np.fft.fftfreq(length, sw[-1]/float(length))
    return float(amp) / sw[-1] * np.exp(2j * np.pi * (pos - x[length//2]) * t - np.pi * np.abs(lor * t) - ((np.pi * np.abs(gauss) * t)**2) / (4 * np.log(2)))

def makeSpectrum(x, sw, v, gauss, lor, weight):
    # Takes axis, frequencies and intensities and makes a spectrum with lorentz and gaussian broadening
    length = len(x)
    t = np.abs(np.fft.fftfreq(length, sw/float(length)))
    diff = (x[1] - x[0]) * 0.5
    final, junk = np.histogram(v, length, range=[x[0]-diff, x[-1]+diff], weights=weight)
    apod = np.exp(-np.pi * np.abs(lor) * t - ((np.pi * np.abs(gauss) * t)**2) / (4 * np.log(2)))
    inten = np.fft.ifft(final) * apod
    inten *= len(inten)  / sw
    return inten

def makeMQMASSpectrum(x, sw, v, gauss, lor, weight):
    length1 = len(x[-2])
    t1 = np.fft.fftfreq(length1, sw[-2]/float(length1))
    diff1 = (x[-2][1] - x[-2][0])*0.5
    length2 = len(x[-1])
    t2 = np.fft.fftfreq(length2, sw[-1]/float(length2))
    diff2 = (x[-1][1] - x[-1][0])*0.5
    t1 = t1[:, np.newaxis]
    final, junk, junk = np.histogram2d(v[0], v[1], [length1, length2], range=[[x[-2][0]-diff1, x[-2][-1]+diff1], [x[-1][0]-diff2, x[-1][-1]+diff2]], weights=weight)
    final = np.fft.ifftn(final)
    apod2 = np.exp(-np.pi * np.abs(lor[1] * t2) - ((np.pi * np.abs(gauss[1]) * t2)**2) / (4 * np.log(2)))
    apod1 = np.exp(-np.pi * np.abs(lor[0] * t1) - ((np.pi * np.abs(gauss[0]) * t1)**2) / (4 * np.log(2)))
    final *= apod1 * apod2 * length1 / sw[-2] * length2 / sw[-1]
    return final

def csaFunc(x, freq, sw, axMult, extra, bgrnd, spinspeed, t11, t22, t33, amp, lor, gauss):
    x = x[-1]
    shiftdef, numssb, angle, D2, weight = extra
    freq = freq[-1]
    sw = sw[-1]
    spinspeed *= 1e3
    tensor = np.array(func.shiftConversion([t11, t22, t33], shiftdef)[1])
    tensor /= axMult
    A0, A2 = csaSpace(tensor)
    T0, T2 = csaSpin()
    d2 = d2tens(np.array([angle]))[0,:,2]
    if spinspeed != np.inf and spinspeed != 0.0:
        gammastep = 2 * np.pi / numssb
        gval = np.arange(numssb) * gammastep
        dt = 1.0 / spinspeed / numssb
        spinD2 = np.exp(1j * np.arange(-2, 3)[:,np.newaxis] * gval) * d2[:,np.newaxis]
    factor2 = d2[2]
    dat0 = A0 * T0
    dat2 = A2 * T2
    dat2 = np.matmul(dat2, D2)
    if spinspeed == np.inf:
        v = np.real(dat2[:,2] * factor2 + dat0)
        tot = weight
    elif spinspeed == 0.0:
        v = np.real(dat2[:,2] + dat0)
        tot = weight
    else:
        vConstant = np.real(dat0 + dat2[:,2] * factor2) 
        dat2[:,2] = 0
        v = np.matmul(dat2, spinD2)
        prod = np.exp(1j * np.cumsum(v * dt * 2 * np.pi, axis=1))
        tot = np.fft.fft(prod, axis=1)
        tot *= np.conj(tot)
        weight2 = weight[:,np.newaxis] / numssb**2
        tot *= weight2
        v = np.fft.fftfreq(numssb, 1.0 / numssb) * spinspeed
        v = v + vConstant[:,np.newaxis]
    return amp * makeSpectrum(x, sw, v, lor, gauss, tot)
    
def quadFunc(x, freq, sw, axMult, extra, bgrnd, spinspeed, pos, cq, eta, amp, lor, gauss):
    x = x[-1]
    satBool, I, numssb, angle, D2, D4, weight = extra
    if not satBool and (I % 1) == 0.0:
        # Integer spins have no central transition
        return np.zeros_like(x)
    cq *= 1e6
    if eta < 0.0:
        eta = 0.0
    elif eta > 1.0:
        eta = 1.0
    spinspeed *= 1e3
    freq = freq[-1]
    sw = sw[-1]
    pos /= axMult
    mList = np.arange(-I, I)
    totalEff = len(mList) * (I**2 + I) - np.sum(mList * (mList + 1))
    if not satBool:
        mList = [-0.5]
    spectrum = np.zeros(len(x), dtype=complex)
    for m in mList:
        eff = I**2 + I - m * (m + 1)
        eff /= totalEff
        v, tot = quadFreq(I, m, m+1, spinspeed, numssb, angle, D2, D4, weight, freq, pos, cq, eta)
        spectrum += eff * makeSpectrum(x, sw, v, lor, gauss, tot)
    return amp * spectrum

def quadFreq(I, m1, m2, spinspeed, numssb, angle, D2, D4, weight, freq, pos, cq, eta):
    pre2 = -cq**2 / (4 * I *(2 * I - 1))**2 * 2 / freq
    pre1 = cq / (4 * I *(2 * I - 1))
    firstA2 = pre1 * firstQuadSpace(eta)
    secA0, secA2, secA4 = secQuadSpace(eta)
    secA0 *= pre2
    secA2 *= pre2
    secA4 *= pre2
    d2 = d2tens(np.array([angle]))[0,:,2]
    d4 = d4tens(np.array([angle]))[0,:,4]
    if spinspeed != np.inf and spinspeed != 0.0:
        gammastep = 2 * np.pi / numssb
        gval = np.arange(numssb) * gammastep
        dt = 1.0 / spinspeed / numssb
        spinD2 = np.exp(1j * np.arange(-2, 3)[:,np.newaxis] * gval) * d2[:,np.newaxis]
        spinD4 = np.exp(1j * np.arange(-4, 5)[:,np.newaxis] * gval) * d4[:,np.newaxis]
    factor2 = d2[2]
    factor4 = d4[4]
    firstspin2 = firstQuadSpin(I, m1, m2)
    secspin0, secspin2, secspin4 = secQuadSpin(I, m1, m2)
    dat0 = secA0 * secspin0 + pos
    dat2 = secA2 * secspin2 + firstA2 * firstspin2
    dat4 = secA4 * secspin4
    dat2 = np.matmul(dat2, D2)
    dat4 = np.matmul(dat4, D4)
    if spinspeed == np.inf:
        v = np.real(dat4[:,4] * factor4  + dat2[:,2] * factor2 + dat0)
        tot = weight
    elif spinspeed == 0.0:
        v = np.real(dat4[:,4]  + dat2[:,2] + dat0)
        tot = weight
    else:
        vConstant = np.real(dat0 + dat2[:,2] * factor2 + dat4[:,4] * factor4) 
        dat4[:,4] = 0
        dat2[:,2] = 0
        v = np.matmul(dat2, spinD2) + np.matmul(dat4, spinD4)
        prod = np.exp(1j * np.cumsum(v * dt * 2 * np.pi, axis=1))
        tot = np.fft.fft(prod, axis=1)
        tot *= np.conj(tot)
        weight2 = weight[:,np.newaxis] / numssb**2
        tot *= weight2
        v = np.fft.fftfreq(numssb, 1.0 / numssb) * spinspeed
        v = v + vConstant[:,np.newaxis]
    return v, tot

def quadCzjzekFunc(x, freq, sw, axMult, extra, bgrnd, d, pos, sigma, wq0, eta0, amp, lor, gauss):
    x = x[-1]
    freq = freq[-1]
    sw = sw[-1]
    method, lib, wq, eta = extra
    if method == 0:
        wq0 = 0
        eta0 = 0
    pos /= axMult
    sigma *= 1e6
    wq0 *= 1e6
    czjzek = Czjzek.czjzekIntensities(sigma, d, wq, eta, wq0, eta0)
    fid = np.dot(czjzek, lib)
    length = len(x)
    t = np.fft.fftfreq(length, sw/float(length))
    pos -= x[len(x)//2]
    apod = np.exp(2j * np.pi * pos * t - np.pi * np.abs(lor) * np.abs(t) - ((np.pi * np.abs(gauss) * t)**2) / (4 * np.log(2)))
    return fid * apod * amp

def mqmasFunc(x, freq, sw, axMult, extra, bgrnd, spinspeed, pos, cq, eta, amp, lor2, gauss2, lor1, gauss1):
    x1 = x[-2]
    x2 = x[-1]
    freq1 = freq[-2]
    freq2 = freq[-1]
    I, mq, numssb, angle, D2, D4, weight, shear, scale = extra
    pos /= axMult
    spinspeed *= 1e3
    cq *= 1e6
    v2, tot2 = quadFreq(I, -0.5, 0.5, spinspeed, numssb, angle, D2, D4, weight, freq2, pos, cq, eta)
    v1, tot1 = quadFreq(I, -mq/2.0, mq/2.0, spinspeed, numssb, angle, D2, D4, np.ones_like(weight), freq1, mq*pos, cq, eta)
    v1 -= v2 * shear
    v1 *= scale
    return amp * makeMQMASSpectrum(x, sw, [np.real(v1.flatten()), np.real(v2.flatten())], [gauss1, gauss2], [lor1, lor2], np.real(tot1*tot2).flatten())

def genLib(length, minWq, maxWq, minEta, maxEta, numWq, numEta, extra, freq, sw, spinspeed):
    I = extra[1]
    wq, eta = np.meshgrid(np.linspace(minWq, maxWq, numWq), np.linspace(minEta, maxEta, numEta))
    wq = wq.flatten()
    eta = eta.flatten()
    cq = wq * (4 * I * (2 * I - 1) / (2 * np.pi))
    x = np.fft.fftshift(np.fft.fftfreq(length, 1/float(sw)))
    lib = np.zeros((len(wq), length), dtype=complex)
    for i, (cqi, etai) in enumerate(zip(cq, eta)):
        lib[i] = quadFunc([x], [freq], [sw], 1.0, extra, 0.0, spinspeed, 0.0, cqi, etai, 1.0, 0.0, 0.0)
    return lib, wq*1e6, eta

def mqmasCzjzekFunc(x, freq, sw, axMult, extra, bgrnd, d, pos, sigma, sigmaCS, wq0, eta0, amp, lor2, gauss2, lor1, gauss1):
    I, mq, wq, eta, lib, shear, scale, method = extra
    if method == 1:
        wq0 *= 1e6
    else:
        wq0 = 0
        eta0 = 0
    pos /= axMult
    sigma *= 1e6
    wq0 *= 1e6
    czjzek = Czjzek.czjzekIntensities(sigma, d, wq, eta, wq0, eta0)
    length2 = len(x[-1])
    czjzek *= length2 / sw[-1] / sw[-2]
    newLib = czjzek[...,np.newaxis]*lib
    length1 = len(x[-2])
    t1 = np.fft.fftfreq(length1, sw[-2]/float(length1))
    t1 = t1[:, np.newaxis]
    diff1 = (x[-2][1] - x[-2][0])*0.5
    t2 = np.fft.fftfreq(length2, sw[-1]/float(length2))
    diff2 = (x[-1][1] - x[-1][0])*0.5
    apod2 = np.exp(-np.pi * np.abs(lor2 * t2) - ((np.pi * np.abs(gauss2) * t2)**2) / (4 * np.log(2)))
    apod1 = np.exp(-np.pi * np.abs(lor1 * t1) - ((np.pi * np.abs(gauss1) * t1)**2) / (4 * np.log(2)))
    cq = wq*2*I*(2*I-1)/(2*np.pi)
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
    for i in range(len(ind)):
        fid[ind[i]-1] += newLib[i]
    fid = np.fft.ifft(fid, axis=0)
    posIndirect = pos * (mq - shearFactor) * scale
    offsetMat = np.exp(2j * np.pi * ((posIndirect - x[-2][length1//2])*t1 + (pos - x[-1][length2//2])*t2))
    shiftGauss = np.exp(-((np.pi * np.abs(sigmaCS) * (t2 + t1*(mq-shearFactor)*scale))**2) / (4 * np.log(2)))
    fid *= offsetMat * apod1 * apod2 * shiftGauss
    shearMat = np.exp((shearFactor-shear) * 2j * np.pi * t1 * x[-1])
    fid = np.fft.fft(fid, axis=1) * shearMat
    return amp * fid
