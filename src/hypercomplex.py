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

import warnings
import numpy as np

def parity(x):
    # Find the parity of an integer
    parity = False
    if x < 0:
        raise RuntimeError('Parity does not work for negative values')
    while x:
        parity = not parity
        x = x & (x - 1)
    return int(parity)


class HComplexException(Exception):
    pass

#########################################################################
# the hyper complex data class


class HComplexData(object):
    """
    Data object for hypercomplex data.
    Methods which depend on the hypercomplex nature of the data are part of this class.

    The first dimension of self.data contains the hypercomplex matrices.
    self.hyper has the same length as self.data and represents whether a matrix is imaginary along a certain dimension.
    self.hyper is encoded in a binary way (0 means real, 1 means imaginary) with the least significant bit representing the first dimension.
    For example, when self.hyper = [0,1,2,3] this means self.data[0] is real in all dimensions,
    self.data[1] is imaginary in the first dimension and real in the second,
    self.data[2] is real in the first dimension and imaginary in the second,
    self.data[3] is imaginary in both first and second dimension.
    self.data contains complex values, the imaginary values are from the last dimension which is always complex and is not listed in self.hyper.
    """

    def __init__(self, data=None, hyper=None):
        """
        Initializes the HComplexData

        Parameters
        ----------
        data : array_like, optional
            Data to be used as hypercomplex data.
            If data is None, an empty hypercomplex array is created.
        hyper : array of ints, optional
            should have the same length as data.
            Is used as the hyper list of the hypercomplex data
            If hyper is None, data is assumed to be regular complex data and an additional dimension
            is added to hold the hypercomplex information.
        """
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
                    raise HComplexException('Length of hyper and data mismatch')
                self.data = np.array(data, dtype=complex)
                self.hyper = np.array(hyper)

    def ndim(self):
        """
        Number of dimensions of the hypercomplex data, which is one less than the dimension of self.data.

        Returns
        -------
        int
            Number of dimensions.
        """
        return self.data.ndim - 1

    def shape(self):
        """
        Shape of the hypercomplex data, which does not include the first dimension of self.data.

        Returns
        -------
        tuple of ints
            Shape of the data.
        """
        return self.data[0].shape

    def getHyperData(self, hyperVal):
        """
        Returns the complex data corresponding to a specific hyper value.

        Parameters
        ----------
        hyperVal : int
            The value of hyper for which the complex data should be returned.

        Returns
        -------
        ndarray
            Complex data corresponding to hyperValue.
        """
        return self.data[hyperVal == self.hyper][0]

    def __repr__(self, *args):
        return self.__class__.__name__ + '(' + repr(self.data) + ', ' + repr(self.hyper) + ')'

    def __len__(self):
        if len(self.data) > 0:
            return len(self.data[0])

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
        return HComplexData(self.data, np.copy(self.hyper))

    def __abs__(self):
        return HComplexData(np.abs(self.data), np.copy(self.hyper))

    def __add__(self, other):
        tmpData = self.copy()
        return tmpData.__iadd__(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        if isinstance(other, HComplexData):
            tmpHyper = np.unique(np.concatenate((self.hyper, other.hyper)))
            tmpHyper.sort()
            tmpData = np.zeros((len(tmpHyper),) + np.broadcast(self.data[0], other.data[0]).shape, dtype=complex)
            for i in self.hyper:
                tmpData[i == tmpHyper] = self.data[i == self.hyper]
            for i in other.hyper:
                tmpData[i == tmpHyper] += other.data[i == other.hyper]
            self.data = tmpData
            self.hyper = tmpHyper
        else:
            self.data[0] += other # Real values should only be added to the real data
        return self

    def __sub__(self, other):
        if isinstance(other, list):
            other = np.asarray(other)
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __isub__(self, other):
        if isinstance(other, list):
            other = np.asarray(other)
        return self.__iadd__(-other)

    def __mul__(self, other):
        tmpData = self.copy()
        return tmpData.__imul__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __imul__(self, other):
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
            self.data = tmpData
            self.hyper = tmpHyper
        else:
            self.data *= other
        return self

    def __div__(self, other):
        tmpData = self.copy()
        return tmpData.__idiv__(other)

    # TODO: implement inverse division

    def __idiv__(self, other):
        if isinstance(other, HComplexData):
            if len(other.hyper) > 1:
                # Recursive calculation of the multicomplex division
                warnings.warn("Calculation of multicomplex data may not result in the correct value")
                tmpOther = HComplexData(np.copy(other.data), np.copy(other.hyper))
                while not tmpOther.isAllReal():
                    tmpObj = HComplexData(np.copy(tmpOther.data), np.copy(tmpOther.hyper))
                    tmpObj = tmpObj.conjAll()
                    tmpOther *= tmpObj
                    self *= tmpObj
                self.data /= tmpOther.data[0]
                # Zero divisors might introduce incorrect division by zero errors
            else:
                self.data /= other.data
        else:
            self.data /= other
        return self

    def __truediv__(self, other):
        return self.__div__(other)

    def __pow__(self, other):
        tmpData = self.copy()
        return tmpData.__ipow__(other)

    def __rpow__(self, other):
        if len(self.hyper) > 1:
            raise HComplexException('Power with more than one complex axis is not permitted')
        return HComplexData(other**self.data, np.copy(self.hyper))

    def __ipow__(self, other):
        if isinstance(other, HComplexData):
            if len(self.hyper) > 1 or len(other.hyper) > 1:
                raise HComplexException('Division of data with more than one complex axis is not permitted')
            self.data **= other.data
        else:
            if len(self.hyper) > 1:
                raise HComplexException('Power with more than one complex axis is not permitted')
            self.data **= other

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
            self.data[(np.in1d(self.hyper, value.hyper), ) + key] = value.data
        else:
            self.data[(slice(0, 1),) + key] = value
            self.data[(slice(1, None),) + key] = 0

    def conj(self, axis):
        """
        Compute the complex conjugate along a specific axis.

        Parameters
        ----------
        axis : int
            The axis along which to calculate the complex conjugate.
            If the axis is hypercomplex, the hypercomplex conjugate of that axis is returned.
            Otherwise the complex conjugate of the entire data is returned.

        Returns
        -------
        HComplexData
            The complex conjugated data.
        """
        if axis < 0:
            axis = self.ndim() + axis
        if not self.isHyperComplex(axis) or axis == (self.ndim()-1):
            return HComplexData(np.conj(self.data), np.copy(self.hyper))
        tmpData = np.copy(self.data)
        imagBool = np.array(self.hyper & (2**axis), dtype=bool)
        tmpData[imagBool] = -tmpData[imagBool]
        return HComplexData(tmpData, np.copy(self.hyper))

    def conjAll(self):
        """
        Compute the complex conjugate of all dimensions.

        Returns
        -------
        HComplexData
            The complex conjugated data.
        """
        tmpData = np.conj(self.data)
        tmpData[1:] = -tmpData[1:]
        return HComplexData(tmpData, np.copy(self.hyper))

    def isAllReal(self):
        """
        Test if the data contains only real values

        The data is allowed to have hypercomplex dimensions as long as they only contain zeros.

        Returns
        -------
        bool
            True if all imaginary dimensions contain only zeros.
        """
        tmp = 0
        if len(self.hyper) > 1:
            tmp = np.count_nonzero(self.data[1:])
        tmp += np.count_nonzero(tmp.imag)
        return not bool(tmp)

    def insertDim(self, axis):
        """
        Add a dimension in self.hyper.
        This should always be accompanied with a modification of self.data.

        Parameters
        ----------
        axis : int
            The axis along which a dimension was added.
        """
        watershedBits = 2**axis - 1
        lowBits = self.hyper & watershedBits
        self.hyper = (self.hyper - lowBits) * 2 + lowBits

    def removeDim(self, axis):
        """
        Remove a dimension from self.hyper.
        This should always be accompanied with a modification of self.data.

        Parameters
        ----------
        axis : int
            The axis along which a dimension was removed.
            If the axis had hypercomplex data associated with it, the imaginary parts are dropped.
            If axis was the last dimension, the second last dimension is promoted to last dimension.
        """
        if axis < 0:
            axis = self.data.ndim - axis
        if self.isHyperComplex(axis):
            tmpdata = self.real(axis)
            self.data = tmpdata.data
            self.hyper = tmpdata.hyper
        if axis == self.ndim() and self.isHyperComplex(axis-1):
            tmpdata = self.real(axis-1)
            self.data = tmpdata.data
            self.hyper = tmpdata.hyper
        watershedBits = 2**axis - 1
        lowBits = self.hyper & watershedBits
        self.hyper = (self.hyper - lowBits) // 2 + lowBits

    def isComplex(self, axis):
        """
        Test whether an axis is complex.

        Parameters
        ----------
        axis : int
            The axis to test.

        Returns
        -------
        bool
            True means the data along axis is complex.
            This includes the last dimension.
        """
        if axis == (self.ndim()-1):
            return True
        return self.isHyperComplex(axis)

    def isHyperComplex(self, axis):
        """
        Test whether an axis is hypercomplex.

        Parameters
        ----------
        axis : int
            The axis to test.

        Returns
        -------
        bool
            True means the data along axis is hypercomplex.
            This excludes the last dimension.
        """
        return bool(np.max(self.hyper) & (2**axis))

    def real(self, axis=-1):
        """
        Return the real part along axis.

        Parameters
        ----------
        axis : int, optional
            The axis along which to get the real data.
            If the axis is hypercomplex the imaginary data along that axis is removed.
            Otherwise the imaginary data (the last dimension) is set to zero.
            By default the last dimension is used.

        Returns
        -------
        HComplexData
            Real data along axis.
        """
        if axis < 0:
            axis = self.ndim() + axis
        if not self.isHyperComplex(axis):
            return HComplexData(np.real(self.data), np.copy(self.hyper))
        bit = 2**axis
        select = np.logical_not(self.hyper & bit)
        return HComplexData(self.data[select], self.hyper[select])

    def imag(self, axis=-1):
        """
        Return the imaginary part along axis.

        Parameters
        ----------
        axis : int, optional
            The axis along which to get the imaginary data.
            If the axis is hypercomplex the imaginary data is made the new real data.
            The old real data along that dimension is removed.
            By default the last dimension is used.

        Returns
        -------
        HComplexData
            Imaginary data along axis.
        """
        if axis < 0:
            axis = self.ndim() + axis
        if not self.isHyperComplex(axis):
            return HComplexData(np.imag(self.data), np.copy(self.hyper))
        bit = 2**axis
        select = np.array(self.hyper & bit, dtype=bool)
        return HComplexData(self.data[select], self.hyper[select]-bit)

    def abs(self, axis=-1):
        """
        Return the absolute data along axis.

        Parameters
        ----------
        axis : int, optional
            The axis along which to get the absolute data.
            If the axis is hypercomplex the absolute data is calculated from the hypercomplex pair.
            Otherwise the regular complex values are used (the last dimension).
            By default the last dimension is used.

        Returns
        -------
        HComplexData
            Absolute data along axis.
        """
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
                tmpData[i] += np.sqrt(np.real(self.data[idim == self.hyper][0])**2 + np.real(self.data[(idim+bit) == self.hyper][0])**2)
                tmpData[i] += 1j * np.sqrt(np.imag(self.data[idim == self.hyper][0])**2 + np.imag(self.data[(idim+bit) == self.hyper][0])**2)
            elif idim in self.hyper:
                tmpData[i] = self.data[idim == self.hyper]
            elif idim + bit in self.hyper:
                tmpData[i] = self.data[idim == self.hyper]
        return HComplexData(tmpData, tmpHyper)

    def complexReorder(self, axis=0):
        """
        Return data where the regular imaginary values are exchanged with those of axis.

        Parameters
        ----------
        axis : int, optional
            The axis along which to exchange with the regular imaginary values.
            When the data along this axis is not hypercomplex, the data is returned unchanged.
            By default the first dimension is used.

        Returns
        -------
        HComplexData
            Reordered hypercomplex data.
        """
        tmpData = self.copy()
        return tmpData.icomplexReorder(axis)

    def icomplexReorder(self, axis=0):
        """
        Reorders data so that the regular imaginary values are exchanged with those of axis.
        The operation is applied in place for efficiency.

        Parameters
        ----------
        axis : int, optional
            The axis along which to exchange with the regular imaginary values.
            When the data along this axis is not hypercomplex, the data is returned unchanged.
            By default the first dimension is used.

        Returns
        -------
        HComplexData
            A pointer to self
        """
        if not self.isHyperComplex(axis):
            # If the data is not complex along that axis return the data unchanged
            return self
        bit = 2**axis
        bArray = np.array(self.hyper & bit, dtype=bool)
        tmpHyper = np.concatenate((self.hyper, self.hyper[bArray] - bit, self.hyper[np.logical_not(bArray)] + bit))
        tmpHyper = np.unique(tmpHyper)
        tmpHyper.sort()
        tmpData = np.zeros((len(tmpHyper),) + self.data[0].shape, dtype=complex)
        tmpBArray = np.array(self.hyper & bit, dtype=bool)
        tmpData[np.logical_not(tmpBArray)] = np.real(self.data[np.logical_not(bArray)]) + 1j*np.real(self.data[bArray])
        tmpData[tmpBArray] = np.imag(self.data[np.logical_not(bArray)]) + 1j*np.imag(self.data[bArray])
        self.data = tmpData
        self.hyper = tmpHyper
        return self

    def moveaxis(self, axis1, axis2):
        """
        Move axes of an array to new positions.
        Other axes remain in their original order.

        Parameters
        ----------
        axis1 : int or sequence of int
            Original positions of the axes to move. These must be unique.
        axis2 : int or sequence of int
            Destination positions for each of the original axes. These must also be unique.

        Returns
        -------
        HComplexData
            A copy of the data with moved axes.
        """
        if isinstance(axis1, (int, float)):
            axis1 = [axis1]
        if isinstance(axis2, (int, float)):
            axis2 = [axis2]
        axis1 = np.array(axis1)
        axis2 = np.array(axis2)
        axis1[axis1 >= 0] += 1
        axis2[axis2 >= 0] += 1
        tmpData = np.moveaxis(self.data, axis1, axis2)
        return HComplexData(tmpData, np.copy(self.hyper))

    def insert(self, pos, other, axis=-1):
        """
        Insert data along the given axis before the given position.

        The resulting data will have hypercomplex dimensions of both input matrices.

        Parameters
        ----------
        pos : int
            Position where to insert the data.
        other : HComplexData or ndarray
            Data to insert.
            Should have the same dimensions as self.data, except for dimension axis.
        axis : int, optional
            The axis along which to insert the data.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the data with the inserted data.
        """
        if axis < 0:
            axis = self.ndim() - axis
        if not isinstance(other, HComplexData):
            other = HComplexData(other)
        tmpHyper = np.unique(np.concatenate((self.hyper, other.hyper)))
        tmpHyper.sort()
        tmpData = []
        for idim in tmpHyper:
            if idim in self.hyper and idim in other.hyper:
                tmpData.append(np.insert(self.data[idim == self.hyper][0], [pos], other.data[idim == other.hyper][0], axis=axis))
            elif idim in self.hyper:
                tmpData.append(np.insert(self.data[idim == self.hyper][0], [pos], np.zeros_like(other.data[0]), axis=axis))
            elif idim in other.hyper:
                tmpData.append(np.insert(np.zeros_like(self.data[0]), [pos], other.data[idim == other.hyper][0], axis=axis))
            else:
                tmpData.append(np.insert(np.zeros_like(self.data[0]), [pos], np.zeros_like(other.data[0]), axis=axis))
        return HComplexData(np.array(tmpData), tmpHyper)

    def delete(self, pos, axis):
        """
        Return data with sub-arrays deleted along axis.

        Parameters
        ----------
        pos : int or array of int
            Positions to delete.
        axis : int
            The axis along which to delete the data.

        Returns
        -------
        HComplexData
            A copy of the data with the subarrays deleted.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.delete(self.data, pos, axis), np.copy(self.hyper))

    def concatenate(self, axis):
        """
        Return data with the data concatenated along axis.

        Parameters
        ----------
        axis : int
            The axis along which to concatenate.

        Returns
        -------
        HComplexData
            The concatenated data.
        """
        if axis >= 0:
            axis += 1
        tmpData = np.swapaxes(self.data, 0, 1)
        tmpData = np.concatenate(tmpData, axis)
        return HComplexData(tmpData, np.copy(self.hyper))

    def split(self, sections, axis):
        """
        Return data which is split into multiple sections.

        Parameters
        ----------
        sections : int
            The number of sections into which the data should be split.
            The length of the data along axis should be divisible by the number of sections.
        axis : int
            The axis along which to split the data.

        Returns
        -------
        HComplexData
            The split data.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.swapaxes(np.split(self.data, sections, axis), 0, 1), np.copy(self.hyper))

    def states(self, axis, TPPI=False):
        """
        Converts the data for States or States-TPPI.
        This can only be performed on an axis which is not yet hypercomplex.
        After the operation this axis will be hypercomplex.
        This operation is performed in place.

        Parameters
        ----------
        axis : int
            The axis along which to perform the States transform.
        TPPI : bool, optional
            When True a States-TPPI transform is performed, otherwise a States transform is performed.
            Defaults to False.
        """
        if axis < 0:
            axis = self.ndim() + axis
        if self.isComplex(axis):
            raise HComplexException("Data is already complex in dimension " + str(axis+1))
        if self.shape()[axis] % 2 != 0:
            raise HComplexException("For conversion the data has to have even datapoints in dimension " + str(axis+1))
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
        """
        Converts the data for echo-antiecho.
        This can only be performed on an axis which is not yet hypercomplex.
        After the operation this axis will be hypercomplex.
        This operation is performed in place.

        Parameters
        ----------
        axis : int
            The axis along which to perform the echo-antiecho transform.
        """
        if axis < 0:
            axis = self.ndim() + axis
        if self.isComplex(axis):
            raise HComplexException("Data is already complex in dimension " + str(axis+1))
        if self.shape()[axis] % 2 != 0:
            raise HComplexException("For conversion the data has to have even datapoints in dimension " + str(axis+1))
        slicing1 = (slice(None), ) * (axis+1) + (slice(None, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axis)
        slicing2 = (slice(None), ) * (axis+1) + (slice(1, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axis)
        addHyper = self.hyper + 2**axis
        insertOrder = np.searchsorted(self.hyper, addHyper)
        tmp1 = self.data[slicing1] + self.data[slicing2]
        tmp2 = 1j * (self.data[slicing1] - self.data[slicing2])
        self.data = np.insert(tmp1, insertOrder, tmp2, axis=0)
        self.hyper = np.insert(self.hyper, insertOrder, addHyper)

    def mean(self, axis=-1, **kwargs):
        """
        Compute the arithmetic mean along the specified axis.

        Parameters
        ----------
        axis : int, optional
            The axis along which to calcule the mean values.
            Defaults to the last dimension.
        **kwargs
            Additional keyword arguments are passed to Numpy mean.

        Returns
        -------
        HComplexData
            The mean values.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.mean(self.data, axis=axis, **kwargs), np.copy(self.hyper))

    def sum(self, axis=-1, **kwargs):
        """
        Compute the sum along the specified axis.

        Parameters
        ----------
        axis : int, optional
            The axis along which to calcule the sum.
            Defaults to the last dimension.
        **kwargs
            Additional keyword arguments are passed to Numpy sum.

        Returns
        -------
        HComplexData
            The sum.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.sum(self.data, axis=axis, **kwargs), np.copy(self.hyper))

    def max(self, axis=-1):
        """
        Returns the data maxima along the specified axis.
        The maxima are defined with respect to the real data.
        The hypercomplex imaginary data is ignored in determining the maxima.

        Parameters
        ----------
        axis : int, optional
            The axis along which to calcule the maximum.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            The maxima.
        """
        argVals = np.argmax(self.data[0], axis=axis)
        ind = list(np.indices(argVals.shape))
        ind.insert(axis, argVals)
        return HComplexData(self.data[(slice(None), ) + tuple(ind)], np.copy(self.hyper))

    def min(self, axis=-1):
        """
        Returns the data minima along the specified axis.
        The minima are defined with respect to the real data.
        The hypercomplex imaginary data is ignored in determining the minima.

        Parameters
        ----------
        axis : int, optional
            The axis along which to calcule the minimum.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            The minima.
        """
        argVals = np.argmin(self.data[0], axis=axis)
        ind = list(np.indices(argVals.shape))
        ind.insert(axis, argVals)
        return HComplexData(self.data[(slice(None), ) + tuple(ind)], np.copy(self.hyper))

    def argmax(self, axis=-1):
        """
        Return the positions of the maxima along the specified axis as a hypercomplex object.
        The maxima are defined with respect to the real data.
        The hypercomplex imaginary data is ignored in determining the positions of the maxima.

        Parameters
        ----------
        axis : int, optional
            The axis along which to calcule the postions of the maxima.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            The positions of the maxima.
        """
        return HComplexData(np.argmax(self.data[0], axis=axis))

    def argmin(self, axis=-1):
        """
        Return the positions of the minima along the specified axis as a hypercomplex object.
        The minima are defined with respect to the real data.
        The hypercomplex imaginary data is ignored in determining the positions of the minima.

        Parameters
        ----------
        axis : int, optional
            The axis along which to calcule the postions of the minima.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            The positions of the minima.
        """
        return HComplexData(np.argmin(self.data[0], axis=axis))

    def expand_dims(self, axis=-1):
        """
        Expand the shape of the hypercomplex data.

        Parameters
        ----------
        axis : int, optional
            Position in the expanded axes where the new axis is placed.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the data with an additional dimension.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.expand_dims(self.data, axis), np.copy(self.hyper))

    def append(self, values, axis=-1):
        """
        Returns a copy of the hypercomplex data with data added along a specified axis.

        Parameters
        ----------
        values : HComplexData or array_like
            Data to be appended. It should have the correct shape (the same shape as self.data, excluding axis).
        axis : int, optional
            The axis along which the values are appended.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            The appended data.
        """
        if axis >= 0:
            axis += 1
        if isinstance(values, HComplexData):
            # Fix for unequal hyper
            return HComplexData(np.append(self.data, values.data, axis=axis), np.copy(self.hyper))
        return HComplexData(np.append(self.data, values, axis=axis), np.copy(self.hyper))

    def reshape(self, shape):
        """
        Returns a reshaped version of the hypercomplex data.

        Parameters
        ----------
        shape : tuple of ints
            New shape of the data. Should be compatible with the old shape.

        Returns
        -------
        HComplexData
            The reshaped data.
        """
        newShape = tuple(len(self.data)) + shape
        return HComplexData(self.data.reshape(newShape), np.copy(self.hyper))

    def diff(self, axis=-1):
        """
        Calculate the discrete difference along the given axis.

        Parameters
        ----------
        axis : int, optional
            The axis along which to calcule the difference.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            The difference data.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.diff(self.data, axis=axis), np.copy(self.hyper))

    def cumsum(self, axis=-1):
        """
        Calculate the cumulative sum along the given axis.

        Parameters
        ----------
        axis : int, optional
            The axis along which to calcule the cumulative sum.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            The cumulative sum data.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.cumsum(self.data, axis=axis), np.copy(self.hyper))

    def hilbert(self, axis=-1):
        """
        Performs a Hilbert transform on the data.

        Parameters
        ----------
        axis : int, optional
            The axis over which the Hilbert transform is performed.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the data which has been Hilbert transformed.
        """
        import scipy.signal
        if axis >= 0:
            axis += 1
        # note: scipy.signal.hilbert designed  for use in time domain, 
        # with fft, zeros negative freq, ifft instead of
        # ifft, zeros in time domain and fft required for nmr spectra
        # this results in conjugated spectrum
        tmpData = np.conjugate(scipy.signal.hilbert(np.real(self.data), axis=axis))
        return HComplexData(tmpData, np.copy(self.hyper))

    def regrid(self, newX, oldX, axis=-1):
        """
        Regrid the data from specified x-values to new x-values.
        The interpolation is done using the interp1d function from scipy.interpolate.

        Parameters
        ----------
        newX : array_like
            A 1-D array with the new x-values.
        oldX : array_like
            A 1-D array with the old x-values. It should have the same length as the data along axis.
        axis : int, optional
            The axis along which the interpolation is applied.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the data which has been regrid.
        """
        from scipy import interpolate as intp
        if axis >= 0:
            axis += 1
        tmpData = np.apply_along_axis(lambda data, newX, oldX: intp.interp1d(oldX, data, fill_value=0, bounds_error=False)(newX), axis, self.data, newX, oldX)
        return HComplexData(tmpData, np.copy(self.hyper))

    def resize(self, size, pos, axis=-1):
        """
        Resizes the data along a specified axis.
        When the new size is larger than the old size, zeros are added at a specified position.
        When the new size is smaller, datapoints are removed symmetric around the specified position.

        Parameters
        ----------
        size : int
            The new size of the data along axis.
        pos : int
            The position where zeros should be added or datapoints should be removed.
        axis : int, optional
            The axis along which the data is resized.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the data which has been resized.
        """
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
        """
        Reorders the data along a specified axis.
        The rest of the data is filled with zeros.

        Parameters
        ----------
        pos : array_like
            The positions of the subarrays in the new data.
        newLength : int, optional
            The new length of the data along axis.
            By default one more than the largest value in pos is used.
        axis : int, optional
            The axis along which the data is reordered.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the data which has been resized.
        """
        if axis >= 0:
            axis += 1
        if newLength is None:
            newLength = max(pos) + 1
        if (max(pos) >= newLength) or (min(pos) < 0):
            raise HComplexException("Positions out of bounds in reorder")
        newShape = np.array(self.data.shape)
        newShape[axis] = newLength
        tmpData = np.zeros(newShape, dtype=complex)
        slicing = (slice(None), ) * axis + (pos, )
        tmpData = np.zeros(newShape, dtype=complex)
        tmpData[slicing] = self.data
        return HComplexData(tmpData, np.copy(self.hyper))

    def apply_along_axis(self, func, axis, *args, **kwargs):
        """
        Applies a function to 1-D slices of the data.
        Behaves similar to the Numpy apply_along_axis function.

        Parameters
        ----------
        func : function
            This function should accept 1-D arrays. It is applied to 1-D slices of arr along the specified axis.
        axis : int
            Axis along which the data is sliced.

        Returns
        -------
        HComplexData
            A copy of the data to which the function has been applied.
        """
        if axis >= 0:
            axis += 1
        tmpData = np.apply_along_axis(func, axis, self.data, *args, **kwargs)
        return HComplexData(tmpData, np.copy(self.hyper))

    def roll(self, shift, axis=-1):
        """
        Rolls the data along a given axis.

        Parameters
        ----------
        shift : int
            Number of places to roll the data.
        axis : int, optional
            The axis along which the data is rolled.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the rolled data.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.roll(self.data, shift, axis=axis), np.copy(self.hyper))

    def fft(self, axis=-1):
        """
        Performs a Fast Fourier Transform on the data along a given axis.

        Parameters
        ----------
        axis : int, optional
            The axis over which the Fourier transform is performed.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the data which has been Fourier transformed.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.fft.fft(self.data, axis=axis), np.copy(self.hyper))

    def ifft(self, axis=-1):
        """
        Performs a inverse Fast Fourier Transform on the data along a given axis.

        Parameters
        ----------
        axis : int, optional
            The axis over which the inverse Fourier transform is performed.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the data which has been inverse Fourier transformed.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.fft.ifft(self.data, axis=axis), np.copy(self.hyper))

    def fftshift(self, axis=-1):
        """
        Shift the zero-frequency component to the center of the spectrum.

        Parameters
        ----------
        axis : int, optional
            The axis over which to shift the data.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the shifted data.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.fft.fftshift(self.data, axes=axis), np.copy(self.hyper))

    def ifftshift(self, axis=-1):
        """
        Inverse of fftshift.

        Parameters
        ----------
        axis : int, optional
            The axis over which to shift the data.
            Defaults to the last dimension.

        Returns
        -------
        HComplexData
            A copy of the shifted data.
        """
        if axis >= 0:
            axis += 1
        return HComplexData(np.fft.ifftshift(self.data, axes=axis), np.copy(self.hyper))

    def copy(self):
        """
        Returns a copy of the hypercomplex data.

        Returns
        -------
        HComplexData
            A copy of the data.
        """
        return HComplexData(np.copy(self.data), np.copy(self.hyper))
