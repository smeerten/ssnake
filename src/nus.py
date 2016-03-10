#!/usr/bin/env python

# Copyright 2015 Bas van Meerten and Wouter Franssen

#This file is part of ssNake.
#
#ssNake is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ssNake is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with ssNake. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import scipy.optimize

def ent_ffm(missingPoints,fid,posArray):
    fid[posArray] = missingPoints[:len(posArray)] + 1j*missingPoints[len(posArray):]
    spec = np.fft.fft(fid)
    zn = np.fft.fft((np.imag(spec)+1j*np.real(spec))/np.abs(spec))
    return (np.sum(np.abs(spec)), np.append(np.imag(zn[posArray]),np.real(zn[posArray]))) #Current entropy is the norm 

def ffm(inp):
    l=len(inp[1])
    res = scipy.optimize.minimize(ent_ffm, np.zeros(l*2), method='L-BFGS-B', args=inp,jac=True)
    inp[0][inp[1]] = res['x'][:l] + 1j*res['x'][l:]
    return inp[0]

def clean(inp):
    residuals = inp[0]
    mask = inp[1]
    gamma = inp[2]
    stopLevel = inp[3]
    maxIter = inp[4]
    replica = np.zeros(len(residuals))
    for i in range(maxIter):
        findMax = np.argmax(residuals)
        maxAmp = residuals[findMax]
        if maxAmp < stopLevel:
            break
        replica[findMax] += maxAmp * gamma
        residuals -= maxAmp * gamma * np.roll(mask, findMax)  
    replica += residuals
    return replica
