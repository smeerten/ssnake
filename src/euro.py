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

def euro(val, num):
    firstDigit = '%.0e' % val
    firstDigit = int(firstDigit[0])
    order = int(np.log10(val))
    if num < 1:
        return
    numStep = int(num)//3 + 1
    if firstDigit == 1:
        subset = [1,2,5]
    elif firstDigit == 2:
        subset = [2,5,10]
    elif firstDigit==5:
        subset = [5,10,20]
    else:
        return
    returnVal = np.tile(subset,numStep)
    orderArray = np.repeat(range(numStep),3)+order
    returnVal = returnVal*10.0**orderArray
    return returnVal[:num]
