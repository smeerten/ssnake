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


import numpy as np
from scipy.special import wofz
import scipy.constants as SC
import scipy.linalg
from fractions import Fraction

def apodize(t, shift=0.0, lor=None, gauss=None, cos2=[None, None], hamming=None, wholeEcho=False):
    """
    Calculates the window function for apodization.

    Parameters
    ----------
    t : array_like
        The time values at which to calculate the apodization function.
    shift : float, optional
        A shift in time of the function.
        A positive value shift the curve to later times.
        By default a shift of 0.0 is used.
    lor : float, optional
        The Lorentzian component of the apodization window.
        By default Lorentzian apodization is not applied.
    gauss : float, optional
        The Gaussian component of the apodization window.
        By default Gaussian apodization is not applied.
    cos2 : array_like, optional
        Defines the squared cosine apodization component.
        Should have a length of at least two.
        The first value is the frequency (two times the number of periods in the time domain).
        The second value is the phase shift in degrees.
        By default squared cosine apodization is not applied.
    hamming : float, optional
        The Hamming window component.
        By default Hamming apodization is not applied.
    wholeEcho : bool, optional
        When True the apodization window is made symmetric in time to match echo experiments.
        By default wholeEcho is False.

    Returns
    -------
    ndarray
        The apodization curve.
    """
    t2 = t - shift
    axLen = len(t2)
    sw = 1.0/(t[1]-t[0])
    x = np.ones(axLen)
    if lor is not None:
        x *= np.exp(-np.pi * lor * abs(t2))
    if gauss is not None:
        x *= np.exp(-((np.pi * gauss * t2)**2) / (4 * np.log(2)))
    if cos2[0] is not None and cos2[1] is not None:
        x *= (np.cos(np.radians(cos2[1]) + cos2[0] * (-0.5 * shift * np.pi * sw / axLen + np.linspace(0, 0.5 * np.pi, axLen)))**2)
    if hamming is not None:
        alpha = 0.53836  # constant for hamming window
        x *= (alpha + (1 - alpha) * np.cos(hamming * (-0.5 * shift * np.pi * sw / axLen + np.linspace(0, np.pi, axLen))))
    if wholeEcho:
        x[-1:-(int(len(x) / 2) + 1):-1] = x[:int(len(x) / 2)]
    return x

def lpsvd(fullFid, nPredict, maxFreq, forward=False, L=None):
    """
    Performs a LPSVD (Linear Predictive Singular Value Decomposition) on a given FID.
    Both forward and backward prediction are available.
    The method follows: R. Kumaresan, D. W. Tufts IEEE Trans. Acoust. Speech Signal Processing 
    vol. ASSP-30, 837-840, 1982.  

    Parameters
    ----------
    fullFid : array_like
        The FID to use for linear prediction.
    nPredict : int
        The number of points to predict.
        For backward prediction, the points are added to the start of the FID.
        For forward prediction, the points are added at the end of the FID.
    maxFreq : int
        Maximum number of frequencies to include in the prediction.
    forward : bool, optional
        True for forward prediction, False for backward prediction.
        By default backward prediction is used.
    L : int, optional
        Number of datapoints from the beginning of the FID to use in the prediction of the frequency components.
        By default all datapoints from fullFid are used.

    Returns
    -------
    ndarray
        The data including the LPSVD predicted values.
    """
    fid = fullFid[:L]
    N = len(fid)
    M = int(np.floor(N * 3 / 4.0))
    H = scipy.linalg.hankel(fid[1:N-M+1], fid[N-M:])
    U, S, Vh = np.linalg.svd(H, full_matrices=1)
    sigVal = len(S[S > (S[0]*1e-6)]) # Number of significant singular values
    if sigVal > maxFreq:
        sigVal = maxFreq
    bias = np.mean(S[sigVal:])
    S = S[:sigVal] - bias
    U = U[:, :sigVal]
    Sinv = np.diag(1.0/S)
    Vh = Vh[:sigVal]
    q = np.dot(np.dot(np.conj(Vh.T), np.dot(Sinv, np.conj(U.T))), fid[:(N-M)])
    s = np.roots(np.append(-q[::-1], 1)) # Find the roots of the polynomial
    s = s[np.abs(s) < 1.0] # Accept only values within the unit circle
    sLog = np.log(s)
    freq = np.imag(sLog)
    lb = np.real(sLog)
    Nfull = len(fullFid)
    Z = s**np.arange(Nfull)[:, np.newaxis]
    a = np.linalg.lstsq(Z, fullFid)[0]
    amp = np.abs(a)
    phase = np.angle(a)
    if forward:
        xpredict = np.arange(Nfull, Nfull + nPredict)
    else:
        xpredict = np.arange(-nPredict, 0)
    if len(amp) > 0:
        reconstructed = amp * np.exp(xpredict[:, np.newaxis] * (1j * freq + lb) + 1j * phase)
        reconstructed = np.sum(reconstructed, axis=1)
    else:
        reconstructed = np.zeros_like(xpredict)
    if forward:
        data = np.concatenate((fullFid, reconstructed))
    else:
        data = np.concatenate((reconstructed, fullFid))
    return data

def euro(val, num):
    """
    Generates a given number of elements from the Euro series (e.g., 0.1, 0.2, 0.5, 1.0, 2.0, ...) starting at a given value.

    Parameters
    ----------
    val : float
        The value at which to start the Euro series.
        When this value is not part of the euro series an exception is thrown.
    num : int
        Number of elements from the Euro series should be generated.

    Returns
    -------
    ndarray
        The Euro series.

    Raise
    -----
    ValueError
        If val is not in the Euro series
    """
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
        raise ValueError(str(val) + " is not part of the Euro series")
    returnVal = np.tile(subset, numStep)
    orderArray = np.repeat(range(numStep), 3) + order
    returnVal = returnVal * 10.0**orderArray
    return returnVal[:num]

def shiftConversion(Values, Type):
    """
    Generates a list of lists of the various chemical shift tensor definitions in solid-state NMR from a specific set of tensor values.

    Parameters
    ----------
    Values : array_like
        The chemical shift values defined in a specific definition.
        Should have a length of 3.
    Type : int
        The definition in which Values is defined.
        0=standard, 1=xyz, 2=Haeberlen, 3=Hertzfeld-Berger.

    Returns
    -------
    list of lists
        The chemical shift tensor definitions, which are (in order) standard, xyz, Haeberlen, and Herzfeld-Berger.
    """
    if Type == 0:  # If from standard
        deltaArray = Values
    if Type == 1:  # If from xyz
        deltaArray = Values  # Treat xyz as 123, as it reorders them anyway
    if Type == 2:  # From haeberlen
        iso = Values[0]
        delta = Values[1]
        eta = Values[2]
        if eta < 0.0:
            eta = 0.0
        elif eta > 1.0:
            eta = 1.0
        delta11 = delta + iso  # Treat xyz as 123, as it reorders them anyway
        delta22 = (eta * delta + iso * 3 - delta11) / 2.0
        delta33 = iso * 3 - delta11 - delta22
        deltaArray = [delta11, delta22, delta33]
    if Type == 3:  # From Hertzfeld-Berger
        iso = Values[0]
        span = Values[1]
        skew = Values[2]
        if skew < -1.0:
            skew = -1.0
        elif skew > 1.0:
            skew = 1.0
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

def quadConversion(Values, I, Type, Q=None):
    """
    Generates a list of lists of the various quadrupole tensor definitions in solid-state NMR from a specific set of tensor values.

    Parameters
    ----------
    Values : array_like
        The chemical shift values defined in a specific definition.
        Should have a length of 2 for the Cq/eta and Wq/eta definitions and 3 for the Vxx/Vyy/Vzz definition.
    I : float
        Spin quantum number.
    Type : int
        The definition in which Values is defined.
        0=Cq/eta, 1=Wq/eta, 2=Vxx/Vyy/Vzz.
    Q : float, optional
        Quadrupole moment in fm^2.
        Has to be defined when Type is 2.
        By default None.

    Returns
    -------
    list of lists
        The quadrupole tensor definitions, which are (in order) Cq/eta, Wq/eta, and Vxx/Vyy/Vzz.
        When Q is None, Vxx=None, Vyy=None, Vzz=None.

    Raises
    ------
    ValueError
        When Q is undefined when Type is 2.
    """
    if Type == 0:  # Cq as input
        # Cq, eta
        # Czz is equal to Cq, via same definition (scale) Cxx and Cyy can be found
        Czz = Values[0]
        Eta = Values[1]
        if Eta > 1.0:
            Eta = 1.0
        elif Eta < 0.0:
            Eta = 0.0
        Cxx = Czz * (Eta - 1) / 2
        Cyy = -Cxx - Czz
        Values = [Cxx, Cyy, Czz]
    if Type == 1:
        # Wq, eta
        Vmax = Values[0]
        Eta = Values[1]
        if Eta > 1.0:
            Eta = 1.0
        elif Eta < 0.0:
            Eta = 0.0
        Czz = Vmax * (2.0 * I * (2 * I - 1)) / 3.0
        Cxx = Czz * (Eta - 1) / 2
        Cyy = -Cxx - Czz
        Values = [Cxx, Cyy, Czz]
    if Type == 2:
        #Vxx, Vyy, Vzz
        if Q is None:
            raise ValueError("Q cannot be None, when quadrupole type is Vxx/Vyy/Vzz")
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
        Values = [Cxx, Cyy, Czz]
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
    if Q is not None:
        Scaling = SC.elementary_charge * Q / SC.Planck
        Vxx = Csort[0] / Scaling * 1e6
        Vyy = Csort[1] / Scaling * 1e6
        Vzz = Csort[2] / Scaling * 1e6
    else:
        Vxx = None
        Vyy = None
        Vzz = None
    return [[CqNew, EtaNew], [WqNew, EtaNew], [Vxx, Vyy, Vzz]]

def ACMEentropy(phaseIn, data, x, firstOrder=True):
    """
    Calculates the cost value for autophasing.

    Parameters
    ----------
    phaseIn : list of float
        Should contain two values, the first one is the zero order phase, and the second one the first order phase.
    data : ndarray
        The data to be phased.
    x : ndarray
        The x-axis corresponding to the first order phasing.
    firstOrder : bool, optional
        If True, the first order phasing is included.
        True by default.

    Returns
    -------
        The cost value.
    """
    phase0 = phaseIn[0]
    if firstOrder:
        phase1 = phaseIn[1]
    else:
        phase1 = 0.0
    L = len(data)
    s0 = data * np.exp(1j * (phase0 + phase1 * x))
    s2 = np.real(s0)
    ds1 = np.abs((s2[3:L] - s2[1:L - 2]) / 2.0)
    p1 = ds1 / sum(ds1)
    p1[np.where(p1 == 0)] = 1
    h1 = -p1 * np.log(p1)
    H1 = sum(h1)
    Pfun = 0.0
    as1 = s2 - np.abs(s2)
    sumas = sum(as1)
    if np.real(sumas) < 0:
        Pfun = Pfun + sum(as1**2) / 4 / L**2
    return H1 + 1000 * Pfun

# functions for calculating quadrupolar constants and factors

def C0(I, m):
    """Returns the C0 term for second order quadrupolar interaction for energy level m of spin I as a Fraction"""
    mf = Fraction(m).limit_denominator(2)
    If = Fraction(I).limit_denominator(2)
    return mf*(If*(If + 1) - 3*mf**2)

def C2(I, m):
    """Returns the C2 term for second order quadrupolar interaction for energy level m of spin I as a Fraction"""
    mf = Fraction(m).limit_denominator(2)
    If = Fraction(I).limit_denominator(2)
    return mf*(8*If*(If + 1) - 12*mf**2 - 3)

def C4(I, m):
    """Returns the C4 term for second order quadrupolar interaction for energy level m of spin I as a Fraction"""
    mf = Fraction(m).limit_denominator(2)
    If = Fraction(I).limit_denominator(2)
    return mf*(18*If*(If + 1) - 34*mf**2 - 5)

def D0(I, m, n):
    """Returns the D0 term for second order quadrupolar interaction for coherence m->n of spin I as a Fraction
       Sign may not be the conventional one but ratios are not affected as long as order is consistent"""
    return C0(I, m) - C0(I, n)

def D2(I, m, n):
    """Returns the D2 term for second order quadrupolar interaction for coherence m->n of spin I as a Fraction"""
    return C2(I, m) - C2(I, n)

def D4(I, m, n):
    """Returns the D4 term for second order quadrupolar interaction for coherence m->n of spin I as a Fraction"""
    return C4(I, m) - C4(I, n)

def m_n_order(m, n):
    """Returns the coherence order m->n transition for spin I"""
    mf = Fraction(m).limit_denominator(2)
    nf = Fraction(n).limit_denominator(2)
    return (mf - nf)

def R(I, m, n):
    """Returns the slope for MQMAS/STMAS correlation of <m|n> coherence with CT(-1Q) (-1/2|1/2) of spin I as a Fraction
    The shearing ratio will then be -R to align the correlation parallel to D2 axis"""
    return D4(I, m, n)/D4(I, -1/2, 1/2)

def scale_CarRef_ratio(I, m, n):
    """Returns the scaling ratio for MQMAS/STMAS of Em-En transition with CT(-1Q) (-1/2->1/2) of spin I as a Fraction"""
    return - R(I, m, n) - m_n_order(m, n) 
    
def scale_SW_ratio(I, m, n):
    """Returns the shearing ratio for MQMAS/STMAS correlation of m->n transition with CT(-1Q) (-1/2->1/2) of spin I as a Fraction"""
    return 1 / scale_CarRef_ratio(I, m, n) 



