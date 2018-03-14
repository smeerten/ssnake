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

#########################################################################
# the hyper complex data class

class HComplexData(object):

    def __init__(self, data=None, hyper=None):
        if data is None:
            self.data = np.array([], dtype=complex)
            self.hyper = np.array([])
        else:
            if hyper is None:
                # Data is not hypercomplex
                self.data = np.array([data], dtype=complex)
                self.hyper = np.array([0]) # The data is only complex in the last dimension
            else:
                if len(hyper) != len(data):
                    raise RuntimeError('Length of hyper and data mismatch')
                self.data = np.array(data, dtype=complex)
                self.hyper = np.array(hyper)
                
    def __len__(self):
        if self.data:
            return self.data[0]

    def __add__(self, other):
        tmpData = np.copy(self.data)
        tmpHyper = np.copy(self.hyper)
        if isinstance(other, HComplexData):
            tmpData[np.isin(self.hyper, other.hyper, assume_unique=True)] += other.data[np.isin(other.hyper, tmpHyper, assume_unique=True)]
            diffList = np.setdiff1d(other.hyper, tmpHyper, assume_unique=True)
            insertOrder = np.searchsorted(tmpHyper, diffList)
            tmpData = np.insert(tmpData, insertOrder, other.data[np.isin(other.hyper, diffList, assume_unique=True)], axis=0)
            tmpHyper = np.insert(tmpHyper, insertOrder, diffList)
            return HComplexData(tmpData, tmpHyper)
        else:
            tmpData[0] += other # Real values should only be added to the real data
            return HComplexData(tmpData, tmpHyper)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        tmpData = np.copy(self.data)
        tmpHyper = np.copy(self.hyper)
        if isinstance(other, HComplexData):
            tmpData[np.isin(tmpHyper, other.hyper, assume_unique=True)] -= other.data[np.isin(other.hyper, tmpHyper, assume_unique=True)]
            diffList = np.setdiff1d(other.hyper, tmpHyper, assume_unique=True)
            insertOrder = np.searchsorted(tmpHyper, diffList)
            tmpData = np.insert(tmpData, insertOrder, -1*other.data[np.isin(other.hyper, diffList, assume_unique=True)], axis=0)
            tmpHyper = np.insert(tmpHyper, insertOrder, diffList)
            return HComplexData(tmpData, tmpHyper)
        else:
            tmpData[0] -= other
            return HComplexData(tmpData, tmpHyper)

    def __rsub__(self, other):
        tmpData = -1 * self.data
        tmpHyper = np.copy(self.hyper)
        tmpData[0] += other
        return HComplexData(tmpData, tmpHyper)
        
    def ndim(self):
        return self.data.ndim - 1 # One extra dimension to contain the hypercomplex information

    def isComplex(self, axis):
        if axis == (self.ndim()-1):
            return True
        return bool(np.max(self.hyper) & (2**axis))
        
a = HComplexData([[1]], [0])
b = HComplexData([[1], [2]], [0,1])
print (b+10).data
