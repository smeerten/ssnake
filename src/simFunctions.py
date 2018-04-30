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

def d2tens(b):
    cosb = np.cos(b)
    sinb = np.sin(b)
    sin2b = 2*sinb*cosb
    tens = np.zeros((len(b),5,5))
    tens[:,0,0] = 0.25 * (1 + cosb)**2
    tens[:,1,0] = 0.5  *sinb * (1 +cosb)
    tens[:,2,0] = np.sqrt(3.0/8) * sinb**2
    tens[:,3,0] = 0.5 * sinb * (1 - cosb)
    tens[:,4,0] = 0.25 * (1-cosb)**2
    tens[:,0,1] = -tens[:,1,0]
    tens[:,1,1] = 0.5 * (2 * cosb**2 + cosb - 1)
    tens[:,2,1] = np.sqrt(3.0/8) * sin2b
    tens[:,3,1] = -0.5 * (2 * cosb**2 - cosb - 1)
    tens[:,4,1] = tens[:,3,0]
    tens[:,0,2] = tens[:,2,0]
    tens[:,1,2] = -tens[:,2,1]
    tens[:,2,2] = 0.5 * (3 * cosb**2 - 1)
    tens[:,3,2] = tens[:,2,1]
    tens[:,4,2] = tens[:,0,2]
    tens[:,0,3] = -tens[:,4,1]
    tens[:,1,3] = -tens[:,3,1]
    tens[:,2,3] = -tens[:,2,1]
    tens[:,3,3] = tens[:,1,1]
    tens[:,4,3] = -tens[:,0,1]
    tens[:,0,4] =  tens[:,4,0]
    tens[:,1,4] = -tens[:,3,0]
    tens[:,2,4] = tens[:,2,0]
    tens[:,3,4] = -tens[:,1,0]
    tens[:,4,4] = tens[:,0,0]
    return tens

def D2tens(a,b,g):
    d = d2tens(b)
    tens = np.zeros((len(a),5,5),dtype=np.complex128)
    part = [2,1,0,-1,-2]
    for i, elem in enumerate(part):
        for k, elem2 in enumerate(part):
            tens[:,i,k] = np.exp(-1j * (a * elem2 + g *elem))
    tens *= d 
    return np.transpose(tens,(0,2,1)) 

def D2tensFinal(a,b, use = np.array([0,1,2,3,4]), type = 'Complex', d = None):
    part = np.array([2,1,0,-1,-2])[use]
    if type == 'Complex':
        arot = np.zeros((len(a),len(part)),dtype=np.complex128)
        for i, elem in enumerate(part):
            arot[:,i] = np.exp(-1j * a * elem)
    elif type == 'Real':
        arot = np.zeros((len(a),len(part)))
        for i, elem in enumerate(part):
            if elem == 0:
                arot[:,i] = 1
            else:
                arot[:,i] = np.cos(a * elem)
    if d is None:
        d = d2tensFinal(b, use)
    return arot * d

def d2tensFinal(b, use=np.array([0,1,2,3,4])):
    tens = np.zeros((len(b),len(use)))
    if np.any(np.in1d(use,[0,2,4])): #If needed
        cos2b = np.cos(2 * b)
    if np.any(np.in1d(use,[1,3])): #If needed
        sin2b = np.sin(2 * b)
    for i, val in enumerate(use):
        if val == 0 or val == 4:
            tens[:,i] = np.sqrt(3.0/8)/2 * (1 - cos2b)
        elif val == 1:
            tens[:,i] = np.sqrt(3.0/8) * sin2b
        elif val ==2:
            tens[:,i] = 0.25 + 0.75 * cos2b
        elif val==3:
            tens[:,i] = -np.sqrt(3.0/8) * sin2b
    return tens
    
def voigtLine(x, freq, sw, axMult, extra, bgrnd, pos, amp, lor, gauss):
    x = x[-1]
    pos /= axMult
    lor = np.abs(lor)
    gauss = np.abs(gauss)
    axis = x - pos
    if extra[0] == 0:  # Exact: Freq domain simulation via Faddeeva function
        if gauss == 0.0:  # If no gauss, just take lorentz
            f = 1.0 / (np.pi * 0.5 * lor * (1 + (axis / (0.5 * lor))**2))
            return amp * f
        elif lor == 0.0: # If no lorentz, just take gauss
            sigma = gauss / (2 * np.sqrt(2 * np.log(2)))
            f = np.exp(-axis**2/(2 * sigma**2)) / np.sqrt(2 * np.pi * sigma**2)
            return amp * f 
        else:
            sigma = gauss / (2 * np.sqrt(2 * np.log(2)))
            z = (axis + 1j * lor / 2) / (sigma * np.sqrt(2))
            return amp * wofz(z).real / (sigma * np.sqrt(2 * np.pi))
    else:  # Approximation: THOMPSON et al (doi: 10.1107/S0021889887087090 )
        sigma = gauss / (2 * np.sqrt(2 * np.log(2)))
        lb = lor / 2
        f = (sigma**5 + 2.69269 * sigma**4 * lb + 2.42843 * sigma**3 * lb**2 + 4.47163 * sigma**2 * lb**3 + 0.07842 * sigma * lb**4 + lb**5) ** 0.2
        eta = 1.36603 * (lb / f) - 0.47719 * (lb / f)**2 + 0.11116 * (lb / f)**3
        lor = f / (np.pi * (axis**2 + f**2))
        gauss = np.exp(-axis**2 / (2 * f**2)) / (f * np.sqrt(2 * np.pi))
        return amp * (eta * lor + (1 - eta) * gauss)

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

def csaSpace(delta): #CSA spherical tensor cartesian space functions
    #delta in Hz
    v00 = - np.sqrt(1.0/3.0) * (delta[0] + delta[1] + delta[2])
    v20 =   np.sqrt(1.0/6.0) * (2*delta[2] - delta[0] - delta[1])
    v2pm2 = 0.5 * (delta[0] - delta[1])
    return v00, v20, v2pm2

def csaSpin(): #CSA spherical tensor spin space functions
    T00 = -np.sqrt(1.0/3)  #times Iz and B0, but is 1
    T20 = np.sqrt(1.0/6.0) * 2  #times Iz and B0, but is 1
    return T00, T20

def csaAngleStuff(cheng):
    alpha, beta, weight = zcw_angles(cheng, symm=2)
    D = D2tensFinal(alpha,beta, use = [0,2], type='Real')
    return weight, D

def csaFunc(x, freq, sw, axMult, extra, bgrnd, spinspeed, t11, t22, t33, amp, lor, gauss):
    mas = extra[0]
    shiftdef = extra[1]
    extra = extra[2:]
    freq = freq[-1]
    sw = sw[-1]
    tensor = np.array(func.shiftConversion([t11, t22, t33], shiftdef)[1])
    tensor /= axMult
    if mas:
        return csaMASFunc(x, freq, sw, extra, spinspeed, tensor, amp, lor, gauss)
    else:
        return csaStaticFunc(x, freq, sw, extra, tensor, amp, lor, gauss)

def csaStaticFunc(x, freq, sw, extra, tensor, amp, lor, gauss):
    x = x[-1]
    weight, D = extra
    v00, v20, v2pm2 = csaSpace(tensor)
    T00, T20 = csaSpin()
    vecttmp = np.array([2 * v2pm2, v20]) * T20 # already multiply with T20, before the data gets larger
    vect = np.tile(vecttmp,(len(weight),1))
    v = np.sum(D * vect, 1) + v00 * T00
    return amp * makeSpectrum(x, sw, v, gauss, lor, weight)

def csaMASFunc(x, freq, sw, extra, spinspeed, tensor, amp, lor, gauss):
    x = x[-1]
    cheng, numssb = extra
    spinspeed *= 1000
    alpha, beta, weight = zcw_angles(cheng, symm=2)
    D = D2tens(alpha,beta,np.zeros(len(alpha))) #The wigner rotation
    gammastep = 2 * np.pi / numssb
    gval = np.arange(numssb) * gammastep
    dt = 1.0/spinspeed/numssb
    v00, v20, v2pm2 = csaSpace(tensor)
    T00, T20 = csaSpin()
    vecttmp = np.array([v2pm2, 0, v20, 0, v2pm2]) * T20
    vect = np.tile(vecttmp, (len(alpha), 1, 1))
    vect = np.matmul(vect, D)[:,0,:] #Rotate powder
    #Get the frequency for each gamma angle (for all crystallites)
    v = np.zeros((len(alpha), numssb))
    d = d2tensFinal(np.ones(1) * np.arctan(np.sqrt(2)))[0,:] # get the magic angle rotation
    for i in range(numssb):
        rot = np.exp(1j* np.array([2, 1, 0, -1, -2]) * gval[i]) * d #total rotation
        v[:,i] = np.real(np.sum(vect * rot,1))
    v = np.exp(1j * v * dt * 2 *np.pi) #get phase evolution for each gamma point
    prod = np.cumprod(v, 1) #cumulative phase rotation
    tot = np.array(prod) #for the first gamma angle (0) no rotation, make copy of prod
    ph = np.zeros((len(alpha), 1), dtype=np.complex128)
    for i in range(numssb - 1): #efficient cumprod(roll(v,i))
        ph[:,0] = v[:,-(i+1)]
        prod[:,1:] = ph * prod[:,:-1]
        prod[:,0] = ph[:,0]
        tot += prod
    weight = np.array([weight]).transpose() / numssb **2
    tot = np.sum(weight * tot,0)
    tot = np.roll(tot, 1) #Roll to make last point first (is t=0 or 2pi)
    tot = np.real(np.fft.fft(tot))
    posList = np.array(np.fft.fftfreq(numssb, 1.0 / numssb)) * spinspeed + v00 * T00
    return amp * makeSpectrum(x, sw, posList, gauss, lor, tot)

def quad1Func(x, freq, sw, axMult, extra, bgrnd, spinspeed, pos, cq, eta, amp, lor, gauss):
    mas = extra[0]
    extra = extra[1:]
    freq = freq[-1]
    sw = sw[-1]
    pos /= axMult
    if mas:
        return quad1MASFunc(x, freq, sw, extra, spinspeed, pos, cq, eta, amp, lor, gauss)
    else:
        return quad1StaticFunc(x, freq, sw, extra, pos, cq, eta, amp, lor, gauss)

def quad1StaticFunc(x, freq, sw, extra, pos, cq, eta, amp, lor, gauss):
    x = x[-1]
    I, weight, angleStuff = extra
    m = np.arange(-I, I)
    v = []
    cq *= 1e6
    weights = []
    for i in m:
        tmp = (cq / (4 * I * (2 * I - 1)) * (I * (I + 1) - 3 * (i + 1)**2)) - (cq / (4 * I * (2 * I - 1)) * (I * (I + 1) - 3 * (i)**2))
        v = np.append(v, tmp * (angleStuff[0] - eta * angleStuff[1]) + pos)
        weights = np.append(weights, weight)
    return amp * makeSpectrum(x, sw, v, gauss, lor, weights) / (2 * I)

def quad1StaticsetAngleStuff(cheng):
    phi, theta, weight = zcw_angles(cheng, symm=2)
    angleStuff = [0.5 * (3 * np.cos(theta)**2 - 1), 0.5 * np.cos(2 * phi) * (np.sin(theta)**2)]
    return weight, angleStuff

def quad1MASFunc(x, freq, sw, extra, spinspeed, pos, cq, eta, amp, lor, gauss):
    x = x[-1]
    I, cheng, numssb = extra
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
    return amp * makeSpectrum(x, sw, posList, gauss, lor, sidebands)

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

def quad2Func(x, freq, sw, axMult, extra, bgrnd, spinspeed, pos, cq, eta, amp, lor, gauss):
    x = x[-1]
    I, weight, angleStuff = extra
    freq = freq[-1]
    sw = sw[-1]
    pos /= axMult
    cq *= 1e6
    v = -1 / (6 * freq) * (3 * cq / (2 * I * (2 * I - 1)))**2 * (I * (I + 1) - 3.0 / 4) * (angleStuff[0] + angleStuff[1] * eta + angleStuff[2] * eta**2) + pos
    return amp * makeSpectrum(x, sw, v, gauss, lor, weight)

def quad2CzjzekFunc(x, freq, sw, axMult, extra, bgrnd, d, pos, sigma, wq0, eta0, amp, lor, gauss):
    x = x[-1]
    method, lib, wq, eta = extra
    if method == 0:
        wq0 = 0
        eta0 = 0
    freq = freq[-1]
    sw = sw[-1]
    pos /= axMult
    sigma = sigma * 1e6
    czjzek = Czjzek.czjzekIntensities(sigma, d, wq, eta, wq0, eta0)
    fid = np.dot(czjzek, lib)
    length = len(x)
    t = np.fft.fftfreq(length, sw/float(length))
    pos -= x[int(len(x)/2)]
    apod = np.exp(2j*np.pi*pos*t-np.pi * np.abs(lor) * np.abs(t) -((np.pi * np.abs(gauss) * t)**2) / (4 * np.log(2)))
    return fid * apod * amp / sw * len(fid)

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

def mqmasFunc(x, freq, sw, axMult, extra, bgrnd, pos, cq, eta, amp, lor2, gauss2, lor1, gauss1):
    I, mq, weight, angleStuff, shear, scale = extra
    pos /= axMult
    cq *= 1e6
    v0Q, v4Q = mqmasFreq(I, cq, eta, angleStuff)
    C10 = I*(I+1) - 3/4.0
    C14 = -7/18.0 * (18*I*(I+1) - 17/2.0 - 5)
    v2 = pos + (C10 * v0Q + C14 * v4Q) / freq[-1]
    Cp0 = mq * (I*(I+1) - 3/4.0 * mq**2)
    Cp4 = -7/18.0 * mq * (18*I*(I+1) - 17/2.0 * mq**2 - 5)
    #v1 = mq * pos + (Cp0 * v0Q + Cp4 * v4Q) / freq[-2]
    shearFactor = Cp4 / C14 * freq[-1] / freq[-2]
    offset = mq * pos + Cp0 * v0Q / freq[-2] - shearFactor * (pos + C10 * v0Q / freq[-1])
    offset *= scale
    return amp * makeMQMASSpectrum(x, sw, v2, [gauss1, gauss2], [lor1, lor2], weight, offset, (shearFactor-shear))

def mqmasCzjzekFunc(x, freq, sw, axMult, extra, bgrnd, d, pos, sigma, sigmaCS, wq0, eta0, amp, lor2, gauss2, lor1, gauss1):
    I, mq, wq, eta, lib, shear, scale, method = extra
    if method == 1:
        wq0 *= 1e6
    else:
        wq0 = 0
        eta0 = 0
    pos /= axMult
    sigma = sigma * 1e6
    czjzek = Czjzek.czjzekIntensities(sigma, d, wq, eta, wq0, eta0)
    length2 = len(x[1])
    czjzek *= length2 / sw[-1] / sw[-2]
    newLib = czjzek[...,np.newaxis]*lib
    length1 = len(x[0])
    t1 = np.fft.fftfreq(length1, sw[0]/float(length1))
    t1 = t1[:, np.newaxis]
    diff1 = (x[0][1] - x[0][0])*0.5
    t2 = np.fft.fftfreq(length2, sw[1]/float(length2))
    diff2 = (x[1][1] - x[1][0])*0.5
    apod2 = np.exp(-np.pi * np.abs(lor2 * t2) -((np.pi * np.abs(gauss2) * t2)**2) / (4 * np.log(2)))
    apod1 = np.exp(-np.pi * np.abs(lor1 * t1) -((np.pi * np.abs(gauss1) * t1)**2) / (4 * np.log(2)))
    cq = wq*2*I*(2*I-1)/(2*np.pi)
    v0Q = - cq**2 * (3 + eta**2) / (40 * (I*(2*I-1))**2)
    C10 = I*(I+1) - 3/4.0
    C14 = -7/18.0 * (18*I*(I+1) - 17/2.0 - 5)
    Cp0 = mq * (I*(I+1) - 3/4.0 * mq**2)
    Cp4 = -7/18.0 * mq * (18*I*(I+1) - 17/2.0 * mq**2 - 5)
    shearFactor = Cp4 / C14 * freq[-1] / freq[-2]
    offset = Cp0 * v0Q / freq[-2] - shearFactor * (C10 * v0Q / freq[-1])
    offset *= scale
    ind = np.digitize(offset, x[0]-(x[0][1]-x[0][0])/2.0)
    fid = np.zeros((length1, length2), dtype=complex)
    for i in range(len(ind)):
        fid[ind[i]-1] += newLib[i]
    fid = np.fft.ifft(fid, axis=0)
    posIndirect = pos * (mq - shearFactor) * scale
    offsetMat = np.exp(2j*np.pi*((posIndirect - x[0][length1//2])*t1 + (pos - x[1][length2//2]) *t2))
    shiftGauss = np.exp(-((np.pi * np.abs(sigmaCS) * (t2 + t1*(mq-shearFactor)*scale))**2) / (4 * np.log(2)))
    fid *= offsetMat * apod1 * apod2 * shiftGauss
    shearMat = np.exp((shearFactor-shear) * 2j * np.pi * t1 * x[1])
    fid = np.fft.fft(fid, axis=1) * shearMat
    return amp * fid
