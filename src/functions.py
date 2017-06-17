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

