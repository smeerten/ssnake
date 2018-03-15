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

def parity(x):
    # Find the parity of an integer
    parity = False
    if x<0:
        raise RuntimeError('Parity does not work for negative values')
    while x:
        parity = not parity
        x = x & (x - 1)
    return int(parity)

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

    def ndim(self):
        return self.data.ndim - 1 # One extra dimension to contain the hypercomplex information

    def shape(self):
        return self.data[0].shape
    
    def __len__(self):
        if self.data:
            return self.data[0]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return np.all(other.data == self.data) and np.all(other.hyper == self.hyper)
        return False
        
    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return np.any(other.data != self.data) and np.any(other.hyper != self.hyper)
        return True
        
    def __neg__(self):
        return HComplexData(-self.data, np.copy(self.hyper))

    def __pos__(self):
        return HComplexData(+self.data, np.copy(self.hyper))

    def __abs__(self):
        return HComplexData(np.abs(self.data), np.copy(self.hyper))

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
        return self.__add__(-other)
    
    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if isinstance(other, HComplexData):
            tmpHyper = np.unique(np.concatenate((self.hyper, other.hyper)))
            tmpHyper.sort()
            tmpData = np.zeros((len(tmpHyper),) + self.shape(), dtype=complex)
            for i, idim in enumerate(self.hyper):
                for j, jdim in enumerate(other.hyper):
                    if parity(idim & jdim):
                        tmpData[(idim^jdim) == tmpHyper] -= self.data[i] * other.data[j]
                    else:
                        tmpData[(idim^jdim) == tmpHyper] += self.data[i] * other.data[j]
            return HComplexData(tmpData, tmpHyper)
        else:
            return HComplexData(self.data*other, np.copy(self.hyper))

    def __rmul__(self, other):
        return self.__mul__(other)

    def conj(self):
        self.data.conj()
        self.data[1:] = -self.data[1:]

    def isAllReal(self):
        tmp = 0
        if len(self.hyper) > 1:
            tmp = np.count_nonzero(self.data[1:])
        tmp += np.count_nonzero(tmp.imag)
        return not bool(tmp)

    def __div__(self, other):
        if isinstance(other, HComplexData):
            if len(other.hyper) > 1:
                # Recursive calculation of the multicomplex division
                tmpOther = HComplexData(np.copy(other.data), np.copy(other.hyper))
                tmpSelf = HComplexData(np.copy(self.data), np.copy(self.hyper))
                while not tmpOther.isAllReal():
                    tmpObj = HComplexData(np.copy(tmpOther.data), np.copy(tmpOther.hyper))
                    tmpObj.conj()
                    tmpOther *= tmpObj
                    tmpSelf *= tmpObj
                tmpSelf.data /= tmpOther.data[0]
                # Zero divisors might introduce incorrect division by zero errors
                return tmpSelf
            return HComplexData(self.data/other.data, np.copy(self.hyper))
        else:
            return HComplexData(self.data/other, np.copy(self.hyper))

    def __truediv__(self, other):
        return self.__div__(other)

    def __pow__(self, other):
        if isinstance(other, HComplexData):
            if len(self.hyper) > 1 or len(other.hyper) > 1:
                raise RuntimeError('Division of data with more than one complex axis is not permitted')
            return HComplexData(self.data**other.data, np.copy(self.hyper))
        else:
            if len(self.hyper) > 1:
                raise RuntimeError('Power with more than one complex axis is not permitted')
            return HComplexData(self.data**other, np.copy(self.hyper))

    def __rpow__(self, other):
        if len(self.hyper) > 1:
            raise RuntimeError('Power with more than one complex axis is not permitted')
        return HComplexData(other**self.data, np.copy(self.hyper))

    def __getitem__(self, key):
        return HComplexData(self.data[:, key], self.hyper)

    def __setitem__(self, key, value):
        if isinstance(value, HComplexData):
            self.data[:, key] = 0
            diffList = np.setdiff1d(value.hyper, self.hyper, assume_unique=True)
            insertOrder = np.searchsorted(self.hyper, diffList)
            self.data = np.insert(self.data, insertOrder, 0, axis=0)
            self.hyper = np.insert(self.hyper, insertOrder, diffList)
            self.data[np.isin(self.hyper, value.hyper), key] = value.data
        else:
            self.data[0, key] = value
            self.data[1:, key] = 0
    
    def isComplex(self, axis):
        if axis == (self.ndim()-1):
            return True
        return bool(np.max(self.hyper) & (2**axis))
        
a = HComplexData([5, 2.5], [0,3])
b = HComplexData([[2,6], [2.5,4]], [0,3])
c = HComplexData([[10,20,30], [100, 200, 300]], [0,1])
