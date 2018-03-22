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
import scipy.constants as SC


def apodize(t, shift, sw, axLen, lor, gauss, cos2, hamming, wholeEcho=False):
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


def shiftConversion(Values, Type):
    # Calculates the chemical shift tensor based on:
    # Values: a list with three numbers
    # Type: an integer defining the input shift convention
    # Returns a list of list with all calculated values

    if Type == 0:  # If from standard
        deltaArray = Values
    if Type == 1:  # If from xyz
        deltaArray = Values  # Treat xyz as 123, as it reorders them anyway
    if Type == 2:  # From haeberlen
        iso = Values[0]
        delta = Values[1]
        eta = Values[2]
        delta11 = delta + iso  # Treat xyz as 123, as it reorders them anyway
        delta22 = (eta * delta + iso * 3 - delta11) / 2.0
        delta33 = iso * 3 - delta11 - delta22
        deltaArray = [delta11, delta22, delta33]
    if Type == 3:  # From Hertzfeld-Berger
        iso = Values[0]
        span = Values[1]
        skew = Values[2]
        delta22 = iso + skew * span / 3.0
        delta33 = (3 * iso - delta22 - span) / 2.0
        delta11 = 3 * iso - delta22 - delta33
        deltaArray = [delta11, delta22, delta33]
    Results = []  # List of list with the different definitions
    # Force right order
    deltaSorted = np.sort(deltaArray)
    D11 = deltaSorted[2]
    D22 = deltaSorted[1]
    D33 = deltaSorted[0]
    Results.append([D11, D22, D33])
    # Convert to haeberlen convention and xxyyzz
    iso = (D11 + D22 + D33) / 3.0
    xyzIndex = np.argsort(np.abs(deltaArray - iso))
    zz = deltaArray[xyzIndex[2]]
    yy = deltaArray[xyzIndex[0]]
    xx = deltaArray[xyzIndex[1]]
    Results.append([xx, yy, zz])
    aniso = zz - iso
    if aniso != 0.0:  # Only is not zero
        eta = (yy - xx) / aniso
    else:
        eta = 'ND'
    Results.append([iso, aniso, eta])
    # Convert to Herzfeld-Berger Convention
    span = D11 - D33
    if span != 0.0:  # Only if not zero
        skew = 3.0 * (D22 - iso) / span
    else:
        skew = 'ND'
    Results.append([iso, span, skew])
    return Results


def quadConversion(Values,I, Type, Q = None):
    # Calculates the chemical shift tensor based on:
    # Values: a list with two or three numbers (Cq/Wq and Eta, Or Vxx Vyy Vzz)
    # Type: an integer defining the input shift convention
    # I: spin quntum number
    # Q: Quad moment in fm^2
    # Returns a list of list with all calculated values
        if Type == 0:  # Cq as input
            # Cq, eta
            # Czz is equal to Cq, via same definition (scale) Cxx and Cyy can be found
            Czz = Values[0]
            Eta = Values[1]
            Cxx = Czz * (Eta - 1) / 2
            Cyy = -Cxx - Czz
            Values = [ Cxx, Cyy, Czz]
        if Type == 1:
            #Wq, eta
            Vmax = Values[0]
            Eta = Values[1]
            Czz = Vmax * (2.0 * I * (2 * I - 1)) / 3.0
            Cxx = Czz * (Eta - 1) / 2
            Cyy = -Cxx - Czz
            Values = [ Cxx, Cyy, Czz]
        if Type == 2:
            #Vxx, Vyy, Vzz
            Vxx = Values[0]
            Vyy = Values[1]
            Vzz = Values[2]
            # Force traceless
            if not np.isclose(Vxx + Vyy + Vzz, 0.0):
                Diff = (Vxx + Vyy + Vzz) / 3.0
                Vxx = Vxx - Diff
                Vyy = Vyy - Diff
                Vzz = Vzz - Diff
            Scaling = SC.elementary_charge * Q / SC.Planck
            Czz = Vzz * Scaling / 1e6  # scale for Cq definition in MHz
            Cxx = Vxx * Scaling / 1e6
            Cyy = Vyy * Scaling / 1e6
            Values = [ Cxx, Cyy, Czz]

        #Conversion
        CArray = np.array(Values)
        Cindex = np.argsort(np.abs(CArray))
        Csort = CArray[Cindex]
        if Csort[2] < 0:  # If Czz negative due to weird input, make it positive
            Csort = -Csort
        CqNew = Csort[2]
        if CqNew == 0.0:
            EtaNew = None
        else:
            EtaNew = np.abs((Csort[0] - Csort[1]) / Csort[2])  # Abs to avoid -0.0 rounding error
        WqNew = CqNew * 3.0 / (2.0 * I * (2 * I - 1))
        if Q != None:
            Scaling = SC.elementary_charge * Q / SC.Planck
            Vxx = Csort[0] / Scaling * 1e6
            Vyy = Csort[1] / Scaling * 1e6
            Vzz = Csort[2] / Scaling * 1e6
        else:
            Vxx = None
            Vyy = None
            Vzz = None

        return [[CqNew, EtaNew], [WqNew, EtaNew], [Vxx, Vyy, Vzz]]

