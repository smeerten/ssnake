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
import warnings

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
                self.hyper = np.array([0])
            else:
                if len(hyper) != len(data):
                    raise RuntimeError('Length of hyper and data mismatch')
                self.data = np.array(data, dtype=complex)
                self.hyper = np.array(hyper)

    def ndim(self):
        return self.data.ndim - 1 # One extra dimension to contain the hypercomplex information

    def shape(self):
        return self.data[0].shape

    def getHyperData(self, hyperVal):
        return self.data[hyperVal == self.hyper][0]
    
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
        if isinstance(other, HComplexData):
            tmpHyper = np.unique(np.concatenate((self.hyper, other.hyper)))
            tmpHyper.sort()
            tmpData = np.zeros((len(tmpHyper),) + np.broadcast(self.data[0], other.data[0]).shape, dtype=complex)
            for i in self.hyper:
                tmpData[i==tmpHyper] = self.data[i==self.hyper]
            for i in other.hyper:
                tmpData[i==tmpHyper] += other.data[i==other.hyper]
            return HComplexData(tmpData, tmpHyper)
        else:
            tmpData[0] += other # Real values should only be added to the real data
            return HComplexData(tmpData, self.hyper)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-np.asarray(other))
    
    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if isinstance(other, HComplexData):
            tmpHyper = np.concatenate((self.hyper, other.hyper))
            for i in other.hyper:
                xorHyper = i^self.hyper
                tmpHyper = np.concatenate((tmpHyper, xorHyper))
            tmpHyper = np.unique(tmpHyper)
            tmpHyper.sort()
            tmpData = np.zeros((len(tmpHyper),) + np.broadcast(self.data[0], other.data[0]).shape, dtype=complex)
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
        tmpData = np.conj(self.data)
        tmpData[1:] = -tmpData[1:]
        return HComplexData(tmpData, np.copy(self.hyper))

    def isAllReal(self):
        tmp = 0
        if len(self.hyper) > 1:
            tmp = np.count_nonzero(self.data[1:])
        tmp += np.count_nonzero(tmp.imag)
        return not bool(tmp)

    def insertDim(self, axis):
        watershedBits = 2**axis - 1
        lowBits = self.hyper & watershedBits
        self.hyper = (self.hyper - lowBits) * 2 + lowBits
    
    def removeDim(self, axis):
        if axis < 1:
            raise RuntimeError("Cannot remove axis below 1 from hyper dimensions")
        if self.isComplex(axis):
            self.data = self.real(axis)
        watershedBits = 2**axis - 1
        lowBits = self.hyper & watershedBits
        self.hyper = (self.hyper - lowBits) / 2 + lowBits
    
    def __div__(self, other):
        if isinstance(other, HComplexData):
            if len(other.hyper) > 1:
                # Recursive calculation of the multicomplex division
                warnings.warn("Calculation of multicomplex data may not result in the correct value")
                tmpOther = HComplexData(np.copy(other.data), np.copy(other.hyper))
                tmpSelf = HComplexData(np.copy(self.data), np.copy(self.hyper))
                while not tmpOther.isAllReal():
                    tmpObj = HComplexData(np.copy(tmpOther.data), np.copy(tmpOther.hyper))
                    tmpObj.iconj()
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
        if not isinstance(key, tuple):
            try:
                key = tuple(key)
            except TypeError:
                key = (key, )
        return HComplexData(self.data[(slice(None), ) + key], self.hyper)

    def __setitem__(self, key, value):
        if not isinstance(key, tuple):
            try:
                key = tuple(key)
            except TypeError:
                key = (key, )
        if isinstance(value, HComplexData):
            self.data[(slice(None), ) + key] = 0
            diffList = np.setdiff1d(value.hyper, self.hyper, assume_unique=True)
            insertOrder = np.searchsorted(self.hyper, diffList)
            self.data = np.insert(self.data, insertOrder, 0, axis=0)
            self.hyper = np.insert(self.hyper, insertOrder, diffList)
            self.data[(np.isin(self.hyper, value.hyper), ) + key] = value.data
        else:
            self.data[(slice(0,1), ) + key] = value
            self.data[(slice(1,None), ) + key] = 0
    
    def isComplex(self, axis):
        # Axis ndim-1 are the regular complex numbers, which are always complex
        if axis == (self.ndim()-1):
            return True
        return self.isHyperComplex(axis)

    def isHyperComplex(self, axis):
        return bool(np.max(self.hyper) & (2**axis))

    def real(self, axis=-1):
        if axis < 0:
            axis = self.ndim() + axis
        if not self.isHyperComplex(axis):
            return HComplexData(np.real(self.data), np.copy(self.hyper))
        bit = 2**axis
        select = np.logical_not(self.hyper & bit)
        return HComplexData(self.data[select], self.hyper[select])
        
    def imag(self, axis=-1):
        if axis < 0:
            axis = self.ndim() + axis
        if not self.isHyperComplex(axis):
            return HComplexData(np.imag(self.data), np.copy(self.hyper))
        bit = 2**axis
        select = np.array(self.hyper & bit, dtype=bool)
        return HComplexData(self.data[select], self.hyper[select]-bit)

    def abs(self, axis=-1):
        if axis < 0:
            axis = self.ndim() + axis
        if not self.isHyperComplex(axis):
            return HComplexData(np.abs(self.data), np.copy(self.hyper))
        bit = 2**axis
        bArray = np.array(self.hyper & bit, dtype=bool)
        tmpHyper = np.concatenate((self.hyper[np.logical_not(bArray)], self.hyper[bArray] - bit))
        tmpHyper = np.unique(tmpHyper)
        tmpHyper.sort()
        tmpData = np.zeros((len(tmpHyper),) + self.data[0].shape, dtype=complex)
        for i, idim in enumerate(tmpHyper):
            if idim in self.hyper and (idim+bit) in self.hyper:
                tmpData[i] += np.sqrt(np.real(self.data[idim==self.hyper][0])**2 + np.real(self.data[(idim+bit)==self.hyper][0])**2)
                tmpData[i] += 1j * np.sqrt(np.imag(self.data[idim==self.hyper][0])**2 + np.imag(self.data[(idim+bit)==self.hyper][0])**2)
            elif idim in self.hyper:
                tmpData[i] = self.data[idim==self.hyper]
            elif (idim+bit) in self.hyper:
                tmpData[i] = self.data[idim==self.hyper]
        return HComplexData(tmpData, tmpHyper)

    def complexReorder(self, axis=0):
        if not self.isHyperComplex(axis):
            # If the data is not complex along that axis return the data unchanged
            return HComplexData(np.copy(self.data), np.copy(self.hyper))
        bit = 2**axis
        bArray = np.array(self.hyper & bit, dtype=bool)
        tmpHyper = np.concatenate((self.hyper, self.hyper[bArray] - bit, self.hyper[np.logical_not(bArray)] + bit))
        tmpHyper = np.unique(tmpHyper)
        tmpHyper.sort()
        tmpData = np.zeros((len(tmpHyper),) + self.data[0].shape, dtype=complex)
        tmpBArray = np.array(self.hyper & bit, dtype=bool)
        tmpData[np.logical_not(tmpBArray)] = np.real(self.data[np.logical_not(bArray)]) + 1j*np.real(self.data[bArray])
        tmpData[tmpBArray] = np.imag(self.data[np.logical_not(bArray)]) + 1j*np.imag(self.data[bArray])
        return HComplexData(tmpData, tmpHyper)
        
    def insert(self, pos, other, axis=-1):
        if axis < 0:
            axis = self.ndim() - axis
        if not isinstance(other, HComplexData):
            other = HComplexData(other)
        tmpHyper = np.unique(np.concatenate((self.hyper, other.hyper)))
        tmpHyper.sort()
        tmpData = []
        for i, idim in enumerate(tmpHyper):
            if idim in self.hyper and idim in other.hyper:
                tmpData.append(np.insert(self.data[idim==self.hyper][0], pos, other.data[idim==other.hyper][0], axis=axis))
            elif idim in self.hyper:
                tmpData.append(np.insert(self.data[idim==self.hyper][0], pos, np.zeros_like(other.data[0]), axis=axis))
            elif idim in other.hyper:
                tmpData.append(np.insert(np.zeros_like(self.data[0]), pos, other.data[idim==other.hyper][0], axis=axis))
            else:
                tmpData.append(np.insert(np.zeros_like(self.data[0]), pos, np.zeros_like(other.data[0]), axis=axis))
        return HComplexData(np.array(tmpData), tmpHyper)

    def delete(self, pos, axis):
        if axis >= 0:
            axis += 1
        return HComplexData(np.delete(self.data, pos, axis), np.copy(self.hyper))
    
    def concatenate(self, axis):
        if axis >= 0:
            axis += 1
        tmpData = np.swapaxes(self.data, 0, 1)
        tmpData = np.concatenate(tmpData, axis)
        return HComplexData(tmpData, np.copy(self.hyper))

    def split(self, sections, axis):
        if axis >= 0:
            axis += 1
        return HComplexData(np.swapaxes(np.split(self.data, sections, axis), 0, 1), np.copy(self.hyper))

    def states(self, axis, TPPI=False):
        if self.shape()[axis] % 2 != 0:
            raise RuntimeError("For conversion the data has to have even datapoints in dimension axis")
        slicing1 = (slice(None), ) * (axis+1) + (slice(None, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axis)
        slicing2 = (slice(None), ) * (axis+1) + (slice(1, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axis)
        addHyper = self.hyper + 2**axis
        insertOrder = np.searchsorted(self.hyper, addHyper)
        tmp1 = self.data[slicing1]
        tmp2 = self.data[slicing2]
        if TPPI:
            tmp1[slicing2] *= -1
            tmp2[slicing2] *= -1
        self.data = np.insert(tmp1, insertOrder, tmp2, axis=0)
        self.hyper = np.insert(self.hyper, insertOrder, addHyper)

    def echoAntiEcho(self, axis):
        if self.shape()[axis] % 2 != 0:
            raise RuntimeError("For conversion the data has to have even datapoints in dimension axis")
        slicing1 = (slice(None), ) * (axis+1) + (slice(None, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axis)
        slicing2 = (slice(None), ) * (axis+1) + (slice(1, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axis)
        addHyper = self.hyper + 2**axis
        insertOrder = np.searchsorted(self.hyper, addHyper)
        tmp1 = self.data[slicing1] + self.data[slicing2]
        tmp2 = 1j * (self.data[slicing1] - self.data[slicing2])
        self.data = np.insert(tmp1, insertOrder, tmp2, axis=0)
        self.hyper = np.insert(self.hyper, insertOrder, addHyper)

    def mean(self, axis=-1, **kwargs):
        if axis >= 0:
            axis += 1
        return HComplexData(np.mean(self.data, axis=axis, **kwargs), np.copy(self.hyper))

    def sum(self, axis=-1, **kwargs):
        if axis >= 0:
            axis += 1
        return HComplexData(np.sum(self.data, axis=axis, **kwargs), np.copy(self.hyper))

    def max(self, axis=-1):
        argVals = np.argmax(self.data[0], axis=axis)
        ind = list(np.indices(argVals.shape))
        ind.insert(axis, argVals)
        return HComplexData(self.data[(slice(None), ) + tuple(ind)], np.copy(self.hyper))

    def min(self, axis=-1):
        argVals = np.argmin(self.data[0], axis=axis)
        ind = list(np.indices(argVals.shape))
        ind.insert(axis, argVals)
        return HComplexData(self.data[(slice(None), ) + tuple(ind)], np.copy(self.hyper))

    def argmax(self, axis=-1):
        return HComplexData(np.argmax(self.data[0], axis=axis))

    def argmin(self, axis=-1):
        return HComplexData(np.argmin(self.data[0], axis=axis))

    def expand_dims(self, axis=-1):
        if axis >= 0:
            axis += 1
        return HComplexData(np.expand_dims(self.data, axis), np.copy(self.hyper))

    def append(self, values, axis=-1):
        if axis >= 0:
            axis += 1
        if isinstance(values, HComplexData):
            # Fix for unequal hyper
            return HComplexData(np.append(self.data, values.data, axis=axis), np.copy(self.hyper))
        else:
            return HComplexData(np.append(self.data, values, axis=axis), np.copy(self.hyper))

    def reshape(self, shape):
        newShape = tuple(len(self.data)) + shape
        return HComplexData(self.data.reshape(newShape), np.copy(self.hyper))
        
    def diff(self, axis=-1):
        if axis >= 0:
            axis += 1
        return HComplexData(np.diff(self.data, axis), np.copy(self.hyper))

    def cumsum(self, axis=-1):
        if axis >= 0:
            axis += 1
        return HComplexData(np.cumsum(self.data, axis), np.copy(self.hyper))

    def hilbert(self, axis=-1):
        import scipy.signal
        if axis >= 0:
            axis += 1
        tmpData = scipy.signal.hilbert(np.real(self.data), axis=axis)
        return HComplexData(tmpData, np.copy(self.hyper))

    def regrid(self, newX, oldX, axis=-1):
        from scipy import interpolate as intp
        if axis >= 0:
            axis += 1
        tmpData = np.apply_along_axis(lambda data, newX, oldX: intp.interp1d(oldX, data, fill_value=0, bounds_error=False)(newX), axis, self.data, newX, oldX)
        return HComplexData(tmpData, np.copy(self.hyper))

    def resize(self, size, pos, axis):
        if axis >= 0:
            axis += 1
        oldSize = self.data.shape[axis]
        if size > oldSize:
            slicing1 = (slice(None), ) * axis + (slice(None, pos), )
            slicing2 = (slice(None), ) * axis + (slice(pos, None), )
            zeroShape = np.array(self.data.shape)
            zeroShape[axis] = size - oldSize
            tmpData = np.concatenate((self.data[slicing1], np.zeros(zeroShape), self.data[slicing2]), axis=axis)
        else:
            difference = oldSize - size
            removeBegin = int(np.floor(difference / 2))
            removeEnd = difference - removeBegin
            if pos < removeBegin:
                slicing = (slice(None), ) * axis + (slice(difference, None), )
                tmpData = self.data[slicing]
            elif oldSize - pos < removeEnd:
                slicing = (slice(None), ) * axis + (slice(None, size), )
                tmpData = self.data[slicing]
            else:
                slicing1 = (slice(None), ) * axis + (slice(None, pos - removeBegin), )
                slicing2 = (slice(None), ) * axis + (slice(pos + removeEnd, None), )
                tmpData = np.append(self.data[slicing1], self.data[slicing2], axis=axis)
        return HComplexData(tmpData, np.copy(self.hyper))

    def reorder(self, pos, newLength=None, axis=-1):
        if axis >= 0:
            axis += 1
        if newLength is None:
            newLength = max(pos) + 1
        if (max(pos) >= newLength) or (min(pos) < 0):
            raise RuntimeError("Positions out of bounds in reorder")
        newShape = np.array(self.data.shape)
        newShape[axis] = newLength
        tmpData = np.zeros(newShape, dtype=complex)
        slicing = (slice(None), ) * axis + (pos, )
        tmpData = np.zeros(newShape, dtype=complex)
        tmpData[slicing] = self.data
        return HComplexData(tmpData, np.copy(self.hyper))

    def apply_along_axis(self, func, axis, *args, **kwargs):
        if axis >= 0:
            axis += 1
        tmpData = np.apply_along_axis(func, axis, self.data, *args, **kwargs)
        return HComplexData(tmpData, np.copy(self.hyper))

    def roll(self, shift, axis):
        if axis >= 0:
            axis += 1
        return HComplexData(np.roll(self.data, shift, axis=axis), np.copy(self.hyper))
    
    def fft(self, axis):
        if axis >= 0:
            axis += 1
        return HComplexData(np.fft.fft(self.data, axis=axis), np.copy(self.hyper))

    def ifft(self, axis):
        if axis >= 0:
            axis += 1
        return HComplexData(np.fft.ifft(self.data, axis=axis), np.copy(self.hyper))

    def fftshift(self, axis):
        if axis >= 0:
            axis += 1
        return HComplexData(np.fft.fftshift(self.data, axes=axis), np.copy(self.hyper))

    def ifftshift(self, axis):
        if axis >= 0:
            axis += 1
        return HComplexData(np.fft.ifftshift(self.data, axes=axis), np.copy(self.hyper))

