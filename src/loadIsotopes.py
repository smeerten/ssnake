#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2016 - 2020 Bas van Meerten and Wouter Franssen

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

import sys

def fOrNone(inp):
    """Converts a string to a float and dashes to None"""
    if inp == '-':
        return None
    return float(inp)

def getIsotopeInfo(isoPath):
    """
    Loads the isotope table from a given path.

    Parameters
    ----------
    isoPath : str
        The path to the file with the isotope properties.

    Returns
    -------
    dict
        A dictionary with the isotope properties.
        Unknown or undefined values are set to None.
    """
    if sys.version_info < (3,):
        with open(isoPath) as isoFile:
            isoList = [line.strip().split('\t') for line in isoFile]
    else:
        with open(isoPath, encoding='UTF-8') as isoFile:
            isoList = [line.strip().split('\t') for line in isoFile]
    isoList = isoList[1:] #Cut off header
    nameList = []
    fullNameList = []
    formatNameList = []
    atomNumList = []
    atomMassList = []
    spinList = []
    abundanceList = []
    gammaList = []
    qList = []
    freqRatioList = []
    refSampleList = []
    sampleConditionList = []
    linewidthFactorList = []
    lifetimeList = []
    sensList = []
    for i, _ in enumerate(isoList):
        isoN = isoList[i]
        atomNumList.append(int(isoN[0]))
        nameList.append(isoN[1])
        fullNameList.append(isoN[2])
        atomMassList.append(fOrNone(isoN[3]))
        formatNameList.append(nameList[-1])
        if atomMassList[-1] is not None:
            formatNameList[-1] = '%d' % (atomMassList[i]) + formatNameList[-1]
        spinList.append(fOrNone(isoN[4]))
        abundanceList.append(fOrNone(isoN[5]))
        gammaList.append(fOrNone(isoN[6]))
        qList.append(fOrNone(isoN[7]))
        freqRatioList.append(fOrNone(isoN[8]))
        refSampleList.append(isoN[9])
        sampleConditionList.append(isoN[10])
        if isoN[4] == '0.5' or spinList[i] is None or qList[i] is None:
            linewidthFactorList.append(None)
        else:
            linewidthFactorList.append((2 * spinList[i] + 3) * qList[i]**2 / (spinList[i]**2 * (2 * spinList[i] - 1)))  # Linewidth due to quadrupolar broadening: (2I + 3) * Q /(I^2 * (2I - 1))
        if gammaList[-1] is not None and abundanceList[-1] is not None and spinList[-1] is not None:
            sensList.append(abundanceList[-1] * abs(gammaList[-1])**3 * spinList[-1] * (spinList[-1] + 1))              # Sensitivity: chi * gamma**3 * I * (I + 1)
        else:
            sensList.append(None)
        lifetimeList.append(isoN[11])
    isotopes = {'atomNum':atomNumList, 'name':nameList, 'fullName':fullNameList, 'atomMass':atomMassList,
                'formatName':formatNameList, 'spin':spinList, 'abundance':abundanceList, 'q':qList, 'freqRatio':freqRatioList,
                'refSample':refSampleList, 'sampleCondition':sampleConditionList, 'linewidthFactor':linewidthFactorList,
                'sensitivity':sensList, 'lifetime':lifetimeList, 'gamma':gammaList}
    return isotopes
