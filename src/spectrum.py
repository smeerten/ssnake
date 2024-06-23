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

import copy
import multiprocessing
import itertools
import scipy.optimize
import numpy as np
import nus
import functions as func
import hypercomplex as hc

AUTOPHASETOL = 0.0002 #is ~0.01 degrees


class SpectrumException(Exception):
    pass


#########################################################################
# the generic spectrum class


class Spectrum(object):
    """
    The object that contains all relevant spectral information.
    The functions for processing are methods of this object.
    """

    def __init__(self, data, filePath, freq, sw, spec=None, wholeEcho=None, ref=None, xaxArray=None, customXax=None, history=None, metaData=None, name='', dFilter=None):
        """
        Initializes the Spectrum object.

        Parameters
        ----------
        data : HComplexData
            The NMR data.
        filePath : tuple
            A tuple with the input parameters for the autoLoad function that was used to obtain this data.
        freq : array_like
            The frequency per dimension.
            The list should have the same number of entries as the number of dimensions of data.
        sw : array_like
            The spectral width per dimension.
            The list should have the same number of entries as the number of dimensions of data.
        spec : array_like, optional
            An array with booleans representing whether a dimension is in the frequency domain (True) or in the time domain (False).
            The array should have the same number of entries as the number of dimensions of data.
            By default all dimensions are set to time domain.
        wholeEcho : array_like, optional
            An array with booleans representing whether a dimension is recorded as a whole echo.
            The array should have the same number of entries as the number of dimensions of data.
            By default all dimensions are set to False.
        ref : array_like, optional
            The reference frequency per dimension.
            The array should have the same number of entries as the number of dimensions of data.
            By default the reference frequencies are set to None.
        xaxArray : list of arrays, optional
            The x-axis values per dimension.
            The list should have the same number of entries as the number of dimensions of data.
            By default the x-axes are generated based on the sw values.
        customXax : list of booleans, optional
            Which of the arrays are custom values and have not been automatically generated based on frequency, reference, and spectral width.
            By default all of them are False
        history : list of strings, optional
            The processing history of the data.
            By default the history is set to an empty list.
        metaData : dict, optional
            Contains the metadata.
            The dictionary contains the the keys ('# Scans', 'Acquisition Time [s]', 'Experiment Name','Receiver Gain', 'Recycle Delay [s]', 'Sample', 'Offset [Hz]', 'Time Completed').
            By default all metadata is set to unknown.
        name : str, optional
            A string with the name of the spectrum.
            By default the name is set to an empty string.
        dFilter : float or None, optional
            For a (Bruker) digital filter this value contains the first order phase correction required to correct for the digital filter.
            By default this parameter is set to None.
        """
        self.name = name
        if isinstance(data, hc.HComplexData):
            self.data = data
        else:
            self.data = hc.HComplexData(data)
        self.filePath = filePath
        self.freq = np.array(freq)  # array of center frequency (length is dim, MHz)
        self.sw = np.array(sw, dtype=float)  # array of sweepwidths
        self.dFilter = dFilter # Digital filter first order phase in radian
        self.undoList = []
        self.redoList = []
        self.noUndo = False
        if spec is None:
            self.spec = [0] * self.ndim()
        else:
            self.spec = spec  # int array of length dim where 0 = time domain, 1 = complex spectral
        if wholeEcho is None:
            self.wholeEcho = [False] * self.ndim()
        else:
            self.wholeEcho = wholeEcho  # boolean array of length dim where True indicates a full Echo
        if ref is None:
            self.ref = np.array(self.ndim() * [None])
        else:
            self.ref = np.array(ref, dtype=object)
        if xaxArray is None:
            self.xaxArray = [[] for i in range(self.ndim())]
            self.resetXax()
        else:
            self.xaxArray = xaxArray
            if customXax is not None and len(customXax) == len(self.xaxArray):
                self.customXax = customXax
            else:
                self.customXax = [False] * len(self.xaxArray)
            for i, xax in enumerate(self.xaxArray):
                if xax is None:
                    self.resetXax(i)
        if history is None:
            self.history = []  # list of strings describing all performed operations
        else:
            self.history = history
        if metaData is None:
            self.metaData = {'# Scans': '-', 'Acquisition Time [s]': '-', 'Experiment Name': '-', 'Receiver Gain': '-', 'Recycle Delay [s]': '-',
                             'Sample': '-', 'Offset [Hz]': '-', 'Time Completed': '-'}
        else:
            self.metaData = metaData

    def ndim(self):
        return self.data.ndim()

    def shape(self):
        return self.data.shape()

    def getData(self): # Returns a copy of the data
        return copy.deepcopy(self.data)

    def getHyperData(self, *args):
        return self.data.getHyperData(*args)

    def isComplex(self, *args):
        return self.data.isComplex(*args)

    def rename(self, name):
        """
        Changes the name of the spectrum.

        Parameters
        ----------
        name : str
            The new name.
        """
        self.name = name

    def getHistory(self):
        """
        Returns the history separated by newlines.
        """
        return "\n".join(self.history)

    def addHistory(self, msg):
        """
        Adds a message to the data history.

        Parameters
        ----------
        msg : str
            The message to add to the history list.
        """
        self.history.append(msg)

    def removeFromHistory(self, num=1):
        """
        Gives the history where a given number of messages have been removed from the end of the list.

        Parameters
        ----------
        num : int
            The number of entries to remove.

        Returns
        -------
        list of str
            The history list.
        """
        for i in range(num):
            if self.history:
                val = self.history.pop()
        return val

    def setNoUndo(self, val):
        """
        Sets the "no undo" mode of the data.

        Parameters
        ----------
        val : bool
            When True, the undo and redo lists are cleared and the "no undo" mode is turned on.
        """
        self.noUndo = bool(val)
        if self.noUndo:
            self.undoList = []
            self.redoList = []

    def undo(self):
        """
        Undoes the last operation and puts it in the redo list.

        Returns
        -------
        str
            The undo message.

        Raises
        ------
        SpectrumException
            When the undo list is empty.
        """
        undoFunc = None
        while undoFunc is None and self.undoList:
            undoFunc = self.undoList.pop()
        if undoFunc is None:
            raise SpectrumException("No undo information")
        tmpRedo = self.redoList # Protect the redo list
        undoFunc(self)
        self.redoList = tmpRedo
        self.redoList.append(self.undoList.pop())
        message = self.removeFromHistory(2)
        return "Undo: " + message

    def redo(self):
        """
        Redoes the last undo operation and puts it in the undo list.

        Raises
        ------
        SpectrumException
            When the redo list is empty.
        """
        if not self.redoList:
            raise SpectrumException("No redo information")
        tmpRedo = self.redoList # Protect the redo list
        tmpRedo.pop()(self)
        self.redoList = tmpRedo # Restore the redo list

    def clearUndo(self):
        """
        Clears the undo and redo lists.
        """
        self.undoList = []
        self.redoList = []

    def reload(self):
        """
        Reloads the data based on the filePath of this spectrum.
        """
        import specIO as io
        loadData = io.autoLoad(*self.filePath)
        self.restoreData(loadData, None)

    def checkAxis(self, axis):
        """
        Checks whether a given axis is a valid index for this data.

        Parameters
        ----------
        axis : int
            The index to be tested.

        Returns
        -------
        int
            The index converted to a positive value.

        Raises
        ------
        IndexError
            When the axis is not valid.
        """
        if axis < 0:
            axis += self.ndim()
        if not 0 <= axis < self.ndim():
            raise IndexError("Not a valid axis for Spectrum")
        return axis

    def resetXax(self, axis=None):
        """
        Resets the x-axis based on the spectral width.

        Parameters
        ----------
        axis : int or None, optional
            The axis for which to reset the x-axis.
            When None, all x-axes are reset.
            By default axis is None.
        """
        if axis is not None:
            axis = self.checkAxis(axis)
            val = [axis]
            self.customXax[axis] = False
        else:
            val = range(self.ndim())
            self.customXax = [False] * self.ndim()
        for i in val:
            if self.spec[i] == 0:
                self.xaxArray[i] = np.arange(self.shape()[i]) / (self.sw[i])
            elif self.spec[i] == 1:
                self.xaxArray[i] = np.fft.fftshift(np.fft.fftfreq(self.shape()[i], 1.0 / self.sw[i]))
                if self.ref[i] is not None:
                    self.xaxArray[i] += self.freq[i] - self.ref[i]

    def setXax(self, xax, axis=-1, custom=True):
        """
        Sets the x-axis of a given dimension.

        Parameters
        ----------
        xax : array_like
            The x-axis.
            It should have the same length as the size of the data along dimension axis.
        axis : int, optional
            The dimension for which to set the x-axis.
            By default the last axis is used.
        custom : bool, optional
            Whether the new x-axis is a custom axis.
            True by default

        Raises
        ------
        SpectrumException
            When the length of xax does not match the length of the data
        """
        axis = self.checkAxis(axis)
        if len(xax) != self.shape()[axis]:
            raise SpectrumException("Length of new x-axis does not match length of the data")
        oldXax = self.xaxArray[axis]
        oldCustom = self.customXax[axis]
        self.xaxArray[axis] = xax
        self.customXax[axis] = custom
        self.addHistory("X-axis of dimension " + str(axis + 1) + " was set to " + str(xax).replace('\n', ''))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.setXax(oldXax, axis, oldCustom))

    def insert(self, data, pos, axis=-1):
        """
        Insert data in a given position along a certain dimension.

        Parameters
        ----------
        data : HComplexData or array_like
            The data to insert.
        pos : int
            The index after which to add the data.
        axis : int, optional
            The dimension along which to add the data.
            By default the last dimension is used.

        Raises
        ------
        SpectrumException
            When the inserted spectrum has insufficient dimensions
        """
        if not isinstance(data, hc.HComplexData):
            data = hc.HComplexData(data)
        if self.noUndo:
            returnValue = None
        else:
            if np.all(self.data.hyper == data.hyper) and not self.customXax[axis]: # If both sets have same hyper and it does not have a custom x-axis: easy undo can be used
                returnValue = lambda self: self.delete(range(pos, pos + data.shape()[axis]), axis)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.insert(data, pos, axis))
        axis = self.checkAxis(axis)
        # Check for a change in dimensions
        if len(data.shape()) <= axis:
            raise SpectrumException("Data does not have the correct number of dimensions")
        self.data = self.data.insert(pos, data, axis)
        self.resetXax(axis)
        self.addHistory("Inserted " + str(data.shape()[axis]) + " datapoints in dimension " + str(axis + 1) + " at position " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(returnValue)

    def delete(self, pos, axis=-1):
        """
        Deletes a given range of indices of a certain dimension from the data.

        Parameters
        ----------
        pos : int or array_like
            The indices to remove.
        axis : int, optional
            The dimension along which to delete the data.
            By default the last dimension is used.

        Raises
        ------
        SpectrumException
            When the deletion removes all data from a dimension.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        tmpData = self.data.delete(pos, axis)
        if 0 in tmpData.shape():
            raise SpectrumException('Cannot delete all data')
        self.data = tmpData
        self.xaxArray[axis] = np.delete(self.xaxArray[axis], pos)
        self.customXax[axis] = True
        if isinstance(pos, np.ndarray):
            if pos.ndim == 0:
                pos = int(pos)
        if isinstance(pos, (int, float)):
            length = 1
        else:
            length = len(pos)
        self.addHistory("Removed " + str(length) + " datapoints from dimension " + str(axis + 1) + " at position " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.delete(pos, axis)))

    def __add__(self, other):
        tmpData = copy.deepcopy(self)
        tmpData.add(other)
        return tmpData

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        self.add(other)
        return self

    def add(self, data, axis=None, select=slice(None)):
        """
        Adds given data to the spectrum data.
        The addition follows the Numpy broadcasting rules.

        Parameters
        ----------
        data : Spectrum, HComplexData or ndarray
            The data to add.
        axis : int, optional
            The dimension along which to add the data.
            By default the last dimension is used.
        select : Slice, optional
            An optional selection of the spectrum data on which the addition is performed.
            By default the entire data is used.
        """
        if axis is not None:
            axis = self.checkAxis(axis)
        if isinstance(data, Spectrum):
            data = data.data
        if isinstance(data, np.ndarray):
            if np.prod(data.shape) == 1:
                data = float(data)
        if isinstance(data, np.ndarray) and axis is not None:
            data = data.reshape(data.shape + (1,)*(self.ndim() - axis - 1))
        if not self.noUndo:
            if not isinstance(data, hc.HComplexData):
                returnValue = lambda self: self.subtract(data, axis, select=select)
            elif np.all(self.data.hyper == data.hyper): # If both sets have same hyper: easy subtract can be used for undo
                returnValue = lambda self: self.subtract(data, axis, select=select)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.add(data, axis, select))
        self.data[select] += data
        if isinstance(data, (float, int)):
            self.addHistory("Added " + str(data) + " to data[" + str(select) + "]")
        else:
            self.addHistory("Added to data[" + str(select) + "]")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(returnValue)

    def __sub__(self, other):
        tmpData = copy.deepcopy(self)
        tmpData.subtract(other)
        return tmpData

    def __rsub__(self, other):
        return self.__sub__(other)

    def __isub__(self, other):
        self.subtract(other)
        return self

    def subtract(self, data, axis=None, select=slice(None)):
        """
        Subtract given data from the spectrum data.
        The subtraction follows the Numpy broadcasting rules.

        Parameters
        ----------
        data : Spectrum, HComplexData or ndarray
            The data to subtract.
        axis : int, optional
            The dimension along which to subtract the data.
            By default the last dimension is used.
        select : Slice, optional
            An optional selection of the spectrum data on which the subtraction is performed.
            By default the entire data is used.
        """
        if axis is not None:
            axis = self.checkAxis(axis)
        if isinstance(data, Spectrum):
            data = data.data
        if isinstance(data, np.ndarray):
            if np.prod(data.shape) == 1:
                data = float(data)
        if isinstance(data, np.ndarray) and axis is not None:
            data = data.reshape(data.shape + (1,)*(self.ndim() - axis - 1))
        if not self.noUndo:
            if not isinstance(data, hc.HComplexData):
                returnValue = lambda self: self.add(data, axis, select=select)
            elif np.all(self.data.hyper == data.hyper): #If both sets have same hyper: easy subtract can be used for undo
                returnValue = lambda self: self.add(data, axis, select=select)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.subtract(data, axis, select))
        self.data[select] -= data
        if isinstance(data, (float, int)):
            self.addHistory("Subtracted " + str(data) + " from data[" + str(select) + "]")
        else:
            self.addHistory("Subtracted from data[" + str(select) + "]")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(returnValue)

    def __mul__(self, other):
        tmpData = copy.deepcopy(self)
        tmpData.multiply(other)
        return tmpData

    def __rmul__(self, other):
        return self.__mul__(other)

    def __imul__(self, other):
        self.multiply(other)
        return self

    def multiply(self, data, axis=None, select=slice(None)):
        """
        Multiply given data with the spectrum data.
        The multiplication follows the Numpy broadcasting rules.

        Parameters
        ----------
        data : Spectrum, HComplexData or ndarray
            The data to multiply.
        axis : int, optional
            The dimension along which to multiply the data.
            By default the last dimension is used.
        select : Slice, optional
            An optional selection of the spectrum data on which the multiplication is performed.
            By default the entire data is used.
        """
        if axis is not None:
            axis = self.checkAxis(axis)
        if isinstance(data, Spectrum):
            data = data.data
        if isinstance(data, np.ndarray):
            if np.prod(data.shape) == 1:
                data = float(data)
        if isinstance(data, np.ndarray) and axis is not None:
            data = data.reshape(data.shape + (1,)*(self.ndim() - axis - 1))
        if not self.noUndo:
            if not isinstance(data, hc.HComplexData):
                returnValue = lambda self: self.divide(data, axis, select=select)
            elif np.all(self.data.hyper == data.hyper): # If both sets have same hyper: easy subtract can be used for undo
                returnValue = lambda self: self.divide(data, axis, select=select)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.multiply(data, axis, select))
        self.data[select] *= data
        
        if isinstance(data, (float, int)):
            self.addHistory("Multiplied data[" + str(select) + "] with " + str(data) + " on axis " + str(axis))
        elif isinstance(data, np.ndarray):
            self.addHistory("Multiplied data[" + str(select) + "] with " + str(list(data.flatten())) + " on axis " + str(axis))
        else:
            self.addHistory("Multiplied data[" + str(select) + " on axis " + str(axis))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(returnValue)

    def __div__(self, other):
        tmpData = copy.deepcopy(self)
        tmpData.divide(other)
        return tmpData

    # TODO: implement inverse division

    def __idiv__(self, other):
        self.divide(self, other)
        return self

    def divide(self, data, axis=None, select=slice(None)):
        """
        Divide the spectrum data with the given data.
        The division follows the Numpy broadcasting rules.

        Parameters
        ----------
        data : Spectrum, HComplexData or ndarray
            The data to divide with.
        axis : int, optional
            The dimension along which to divide the data.
            By default the last dimension is used.
        select : Slice, optional
            An optional selection of the spectrum data on which the division is performed.
            By default the entire data is used.
        """
        if axis is not None:
            axis = self.checkAxis(axis)
        if isinstance(data, Spectrum):
            data = data.data
        if isinstance(data, np.ndarray):
            if np.prod(data.shape) == 1:
                data = float(data)
        if isinstance(data, np.ndarray) and axis is not None:
            data = data.reshape(data.shape + (1,)*(self.ndim() - axis - 1))
        if not self.noUndo:
            if not isinstance(data, hc.HComplexData):
                returnValue = lambda self: self.multiply(data, axis, select=select)
            elif np.all(self.data.hyper == data.hyper): #If both sets have same hyper: easy subtract can be used for undo
                returnValue = lambda self: self.multiply(data, axis, select=select)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.divide(data, axis, select))
        self.data[select] /= data
        if isinstance(data, (float, int)):
            self.addHistory("Divided data[" + str(select) + "] with " + str(data))
        else:
            self.addHistory("Divided data[" + str(select) + "]")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(returnValue)

    def normalize(self, mult, scale=1.0, type=0, axis=-1, select=slice(None)):
        """
        Method used by the normalization.

        Parameters
        ----------
        mult : float
            The multiplier to normalize the data.
        scale : float, optional
            The scale to which to set the data.
            By default the scale is 1.
        type : int, optional
            The type of normalization (0=integral, 1=maximum, 2=minimum, only used for the history message).
            By default 0.
        axis : int, optional
            The dimension along which the normalization was performed.
            By default the last dimension is used.
        select : Slice, optional
            An optional selection of the spectrum data on which the normalization was performed.
            By default the entire data is used.

        Raises
        ------
        SpectrumException
            When the multiplication results in an error.
        """
        axis = self.checkAxis(axis)
        try:
            self.data[select] *= mult * scale
        except ValueError as error:
            raise SpectrumException('Normalize: ' + str(error))
        if type == 0:
            self.addHistory("Normalized integral of dimension " + str(axis + 1) + " of data[" + str(select) + "] to " + str(scale))
        elif type == 1:
            self.addHistory("Normalized maximum of dimension " + str(axis + 1) + " of data[" + str(select) + "] to " + str(scale))
        elif type == 2:
            self.addHistory("Normalized minimum of dimension " + str(axis + 1) + " of data[" + str(select) + "] to " + str(scale))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.normalize(1.0 / mult, scale, type, axis, select=select))

    def baselineCorrection(self, baseline, axis=-1, select=slice(None), degree=None, type=None):
        """
        Applies a baseline correction.

        Parameters
        ----------
        baseline : array_like
            The baseline to subtract.
            It follows the Numpy broadcasting rules.
        axis : int, optional
            The dimension along which the baseline correction is performed.
            By default the last dimension is used.
        select : Slice, optional
            An optional selection of the spectrum data on which the baseline correction is performed.
            By default the entire data is used.
        degree : int, optional
            The degree used for the fitting of the baseline, used for history output only.
        type : str, optional
            The type (poly or sin/cos) used for the fitting, used for history output only.
        """
        axis = self.checkAxis(axis)
        baselinetmp = baseline.reshape((self.shape()[axis], ) + (1, ) * (self.ndim() - axis - 1))
        self.data[select] -= baselinetmp
        Message = "Baseline corrected dimension " + str(axis + 1)
        if not isinstance(select, slice):
            Message = Message + " with slice " + str(select)
        elif select != slice(None, None, None):
            Message = Message + " with slice " + str(select)
        if degree and type:
            Message = Message + " ({type}, deg={deg})".format(type=type, deg=str(degree))
        self.addHistory(Message)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.baselineCorrection(-baseline, axis, select=select))

    def concatenate(self, axis=-1):
        """
        Concatenates the data along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension along which to concatenate.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        splitVal = self.shape()[0]
        copyData = None
        if self.data.isComplex(axis) or self.customXax[0] or self.customXax[axis]:
            if not self.noUndo:
                copyData = copy.deepcopy(self)
            if self.data.isComplex(axis):
                self.data = self.data.real(axis)
        invAxis = self.ndim() - axis
        self.data = self.data.concatenate(axis)
        self.data.removeDim(invAxis)
        self.freq = np.delete(self.freq, 0)
        self.sw = np.delete(self.sw, 0)
        self.spec = np.delete(self.spec, 0)
        self.wholeEcho = np.delete(self.wholeEcho, 0)
        self.ref = np.delete(self.ref, 0)
        del self.xaxArray[0]
        del self.customXax[0]
        self.resetXax(axis)
        self.addHistory("Concatenated dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            if copyData is not None:
                self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.concatenate(axis)))
            else:
                self.undoList.append(lambda self: self.split(splitVal, axis))

    def split(self, sections, axis=-1):
        """
        Splits the data into a given number of sections.

        Parameters
        ----------
        sections : int
            The number of sections.
            The length of the data along axis should be divisible by this number.
        axis : int, optional
            The dimension along which to split.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        copyData = None
        if self.customXax[axis]:
            if not self.noUndo:
                copyData = copy.deepcopy(self)
        self.data = self.data.split(sections, axis)
        self.data.insertDim(0)
        self.freq = np.insert(self.freq, 0, self.freq[axis])
        self.sw = np.insert(self.sw, 0, self.sw[axis])
        self.spec = np.insert(self.spec, 0, self.spec[axis])
        self.wholeEcho = np.insert(self.wholeEcho, 0, self.wholeEcho[axis])
        self.ref = np.insert(self.ref, 0, self.ref[axis])
        self.xaxArray.insert(0, [])
        self.customXax.insert(0, False)
        self.resetXax(0)
        self.resetXax(axis + 1)
        self.addHistory("Split dimension " + str(axis + 1) + " into " + str(sections) + " sections")
        self.redoList = []
        if not self.noUndo:
            if copyData is not None:
                self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.split(sections, axis)))
            else:
                self.undoList.append(lambda self: self.concatenate(axis))

    def real(self, axis=-1):
        """
        Sets the data to the real values along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axis = self.checkAxis(axis)
        self.data = self.data.real(axis)
        self.addHistory("Real along dimension " + str(axis+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.real(axis)))

    def imag(self, axis=-1):
        """
        Sets the data to the imaginary values along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axis = self.checkAxis(axis)
        self.data = self.data.imag(axis)
        self.addHistory("Imaginary along dimension " + str(axis+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.imag(axis)))

    def abs(self, axis=-1):
        """
        Sets the data to the absolute values along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axis = self.checkAxis(axis)
        self.data = self.data.abs(axis)
        self.addHistory("Absolute along dimension " + str(axis+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.abs(axis)))

    def conj(self, axis=-1):
        """
        Complex conjugates the data along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        self.data = self.data.conj(axis)
        self.addHistory("Complex conjugate along" + str(axis+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.conj(axis))

    def states(self, axis=-1):
        """
        Converts the data to hypercomplex based on States acquisition along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.states(axis)
        self.resetXax(axis)
        self.addHistory("States conversion on dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.states(axis)))

    def statesTPPI(self, axis=-1):
        """
        Converts the data to hypercomplex based on States-TPPI acquisition along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.states(axis, TPPI=True)
        self.resetXax(axis)
        self.addHistory("States-TPPI conversion on dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.statesTPPI(axis)))

    def echoAntiEcho(self, axis=-1):
        """
        Converts the data to hypercomplex based on echo/antiecho acquisition along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.echoAntiEcho(axis)
        self.resetXax(axis)
        self.addHistory("Echo-antiecho conversion on dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.echoAntiEcho(axis)))

    def subtractAvg(self, pos1=None, pos2=None, axis=-1):
        """
        Subtracts the average values between given indices along a dimension.

        Parameters
        ----------
        pos1 : int, optional
            First index to determine the average.
            0 by default.
        pos2 : int, optional
            Second index to determine the average.
            Length of the data along axis by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.

        Raises
        ------
        SpectrumException
            When pos1 or pos2 are invalid.
        """
        axis = self.checkAxis(axis)
        if pos1 is None:
            pos1 = 0
        if pos2 is None:
            pos2 = self.shape()[axis]
        if not 0 <= pos1 <= self.shape()[axis]:
            raise SpectrumException("Indices not within range")
        if not 0 <= pos2 <= self.shape()[axis]:
            raise SpectrumException("Indices not within range")
        if pos1 == pos2:
            raise SpectrumException("Indices cannot be equal")
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axis + (slice(minPos, maxPos), )
        averages = self.data[slicing].mean(axis=axis, keepdims=True)
        self.data -= averages
        self.addHistory("Subtracted average determined between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.add(averages))

    def matrixManip(self, pos1=None, pos2=None, axis=-1, which=0):
        """
        Performs the matrix manipulation methods such as integrate, sum, average, maximum, minimum, maximum position, and minimum position.

        Parameters
        ----------
        pos1 : int or array_like, optional
            First index/indices of the matrix manipulation.
            0 by default.
        pos2 : int or array_like, optional
            Second index/indices of the matrix manipulation.
            Should have the same length as pos1.
            Length of the data along axis by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        which : int, optional
            The type of matrix manipulation to perform (0=integrate, 1=max, 2=min, 3=maxpos, 4=minpos, 5=sum, 6=average).
            0 by default.

        Raises
        ------
        SpectrumException
            When pos1 and pos2 have unequal lengths or when they contain invalid indices.
        """
        axis = self.checkAxis(axis)
        if pos1 is None:
            pos1 = 0
        if pos2 is None:
            pos2 = self.shape()[axis]
        if isinstance(pos1, int):
            pos1 = np.array([pos1])
            pos2 = np.array([pos2])
        if len(pos1) != len(pos2):
            raise SpectrumException("Length of the two arrays is not equal")
        if len(pos1) == 1:
            keepdims = self.ndim() == 1
        else:
            keepdims = True
        tmpdata = ()
        for i, _ in enumerate(pos1):
            if not 0 <= pos1[i] <= self.shape()[axis]:
                raise SpectrumException("Indices not within range")
            if not 0 <= pos2[i] <= self.shape()[axis]:
                raise SpectrumException("Indices not within range")
            if pos1[i] == pos2[i]:
                raise SpectrumException("Indices cannot be equal")
            minPos = min(pos1[i], pos2[i])
            maxPos = max(pos1[i], pos2[i])
            slicing = (slice(None), ) * axis + (slice(minPos, maxPos), )
            if which == 0:
                if self.spec[axis] == 0:
                    tmp = self.data[slicing].sum(axis=axis) / self.sw[axis]
                else:
                    tmp = self.data[slicing].sum(axis=axis) * self.sw[axis] / (1.0 * self.shape()[axis])
            elif which == 5:
                tmp = self.data[slicing].sum(axis=axis)
            elif which == 1:
                tmp = self.data[slicing].max(axis=axis)
            elif which == 2:
                tmp = self.data[slicing].min(axis=axis)
            elif which == 3:
                tmp = self.data[slicing].argmax(axis=axis)
                maxArgPos = np.array(tmp.data, dtype=int)
                tmpmaxPos = maxArgPos.flatten()
                tmp.data = self.xaxArray[axis][slice(minPos, maxPos)][tmpmaxPos].reshape(maxArgPos.shape)
            elif which == 4:
                tmp = self.data[slicing].argmin(axis=axis)
                minArgPos = np.array(tmp.data, dtype=int)
                tmpminPos = minArgPos.flatten()
                tmp.data = self.xaxArray[axis][slice(minPos, maxPos)][tmpminPos].reshape(minArgPos.shape)
            elif which == 6:
                tmp = self.data[slicing].mean(axis=axis)
            if keepdims:
                tmpdata += (tmp.expand_dims(axis), )
            else:
                tmpdata += (tmp, )
        if not keepdims:
            if self.ndim() == 1:
                self.data = tmpdata[0].reshape((1, ))
                self.resetXax(axis)
            else:
                self.data = tmpdata[0]
                self.freq = np.delete(self.freq, axis)
                self.ref = np.delete(self.ref, axis)
                self.sw = np.delete(self.sw, axis)
                self.spec = np.delete(self.spec, axis)
                self.wholeEcho = np.delete(self.wholeEcho, axis)
                del self.xaxArray[axis]
                del self.customXax[axis]
                self.data.removeDim(axis)
        else:
            self.data = tmpdata[0]
            for extra in tmpdata[1:]:
                self.data = self.data.append(extra, axis=axis)
            self.resetXax(axis)

    def integrate(self, pos1=None, pos2=None, axis=-1):
        """
        Reduce the data to the integrals between given indices.

        Parameters
        ----------
        pos1 : array_like, optional
            First indices to determine the integrals.
            0 by default.
        pos2 : array_like, optional
            Second indices to determine the integrals.
            Should have the same length as pos1.
            Length of the data along axis by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=0)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.integrate(pos1, pos2, axis)))
        self.addHistory("Integrate between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def max(self, pos1=None, pos2=None, axis=-1):
        """
        Reduce the data to the maxima between given indices.

        Parameters
        ----------
        pos1 : array_like, optional
            First indices to determine the maxima.
            0 by default.
        pos2 : array_like, optional
            Second indices to determine the maxima.
            Should have the same length as pos1.
            Length of the data along axis by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=1)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.max(pos1, pos2, axis)))
        self.addHistory("Maximum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def min(self, pos1=None, pos2=None, axis=-1):
        """
        Reduce the data to the minima between given indices.

        Parameters
        ----------
        pos1 : array_like, optional
            First indices to determine the minima.
            0 by default.
        pos2 : array_like, optional
            Second indices to determine the minima.
            Should have the same length as pos1.
            Length of the data along axis by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=2)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.min(pos1, pos2, axis)))
        self.addHistory("Minimum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def argmax(self, pos1=None, pos2=None, axis=-1):
        """
        Reduce the data to the maxima positions between given indices.

        Parameters
        ----------
        pos1 : array_like, optional
            First indices to determine the maxima positions.
            0 by default.
        pos2 : array_like, optional
            Second indices to determine the maxima positions.
            Should have the same length as pos1.
            Length of the data along axis by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=3)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.argmax(pos1, pos2, axis)))
        self.addHistory("Maximum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def argmin(self, pos1=None, pos2=None, axis=-1):
        """
        Reduce the data to the minima positions between given indices.

        Parameters
        ----------
        pos1 : array_like, optional
            First indices to determine the minima positions.
            0 by default.
        pos2 : array_like, optional
            Second indices to determine the minima positions.
            Should have the same length as pos1.
            Length of the data along axis by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=4)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.argmin(pos1, pos2, axis)))
        self.addHistory("Minimum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def sum(self, pos1=None, pos2=None, axis=-1):
        """
        Reduce the data to the sum between given indices.

        Parameters
        ----------
        pos1 : array_like, optional
            First indices to determine the sum.
            0 by default.
        pos2 : array_like, optional
            Second indices to determine the sum.
            Should have the same length as pos1.
            Length of the data along axis by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=5)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.sum(pos1, pos2, axis)))
        self.addHistory("Sum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def average(self, pos1=None, pos2=None, axis=-1):
        """
        Reduce the data to the average between given indices.

        Parameters
        ----------
        pos1 : array_like, optional
            First indices to determine the average.
            0 by default.
        pos2 : array_like, optional
            Second indices to determine the average.
            Should have the same length as pos1.
            Length of the data along axis by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=6)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.average(pos1, pos2, axis)))
        self.addHistory("Average between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def extract(self, pos1=None, pos2=None, axis=-1, step=None):
        """
        Extract the data between two given indices (both indexes included) along a dimension and make it the new data.

        Parameters
        ----------
        pos1 : int, optional
            First index to extract.
            0 by default.
        pos2 : int, optional
            Last index to extract.
            Last index of the data (length-1) along axis by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        step : int, optional
            Step size between indices to extract.
            1 by default
        """
        axis = self.checkAxis(axis)
        if pos1 is None:
            pos1 = 0
        if pos2 is None:
            pos2 = self.shape()[axis] - 1
        if step is None:
            step = 1
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axis + (slice(minPos, maxPos+1, step), )
        if self.spec[axis] == 1:
            oldFxax = self.xaxArray[axis][minPos]
            self.sw[axis] *= (step * np.ceil((maxPos - minPos + 1 )/step)) / (1.0 * self.shape()[axis])
        else:
            self.sw[axis] /= step
        self.data = self.data[slicing]
        if self.spec[axis] == 1:
            newFxax = np.fft.fftshift(np.fft.fftfreq(self.shape()[axis], 1.0 / self.sw[axis]))[0]
            if self.ref[axis] is None:
                self.ref[axis] = self.freq[axis]
            self.freq[axis] = self.ref[axis] - newFxax + oldFxax
        self.resetXax(axis)
        self.addHistory("Extracted part from " + str(minPos) + " to " + str(maxPos) + " with steps of " + str(step) + " of dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.extract(pos1, pos2, axis, step)))

    def fiddle(self, refSpec, lb, axis=-1):
        """
        Performs reference deconvolution using the FIDDLE algorithm.

        Parameters
        ----------
        refSpec : array_like
            The reference spectrum with which to deconvolute the spectrum.
            Should have the same length as the data along axis.
        lb : float
            The linebroadening (in Hz) to be applied during reference deconvolution.
        axis : int, optional
            The dimension.
            By default the last dimension is used.

        Raises
        ------
        SpectrumException
            When the reference FID does not have the same length as the data along axis.
        """
        axis = self.checkAxis(axis)
        axLen = self.shape()[axis]
        if len(refSpec) != axLen:
            raise SpectrumException("Reference FID does not have the correct length")
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        tmpSpec = np.fft.ifftshift(np.real(refSpec))
        pos = np.argmax(tmpSpec)
        refFid = np.fft.ifft(tmpSpec)
        if self.spec[axis] > 0:
            self.__invFourier(axis, tmp=True)
        t = np.arange(axLen) / self.sw[axis]
        idealFid = np.exp(2j * pos * np.pi * np.linspace(0, 1, axLen + 1)[:-1]) * np.exp(-np.pi * lb * t)
        # Make Hilbert transform filter
        h = np.zeros(axLen)
        if axLen % 2 == 0:
            h[0] = h[axLen // 2] = 1
            h[1:axLen // 2] = 2
        else:
            h[0] = 1
            h[1:(axLen + 1) // 2] = 2
        idealFid *= h
        idealFid /= refFid
        idealFid /= np.abs(idealFid[0]) * 2
        self.data *= idealFid
        if self.spec[axis] > 0:
            self.__fourier(axis, tmp=True)
        self.addHistory("FIDDLE over dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.fiddle(refSpec, lb, axis)))

    def diff(self, axis=-1):
        """
        The discrete difference along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.data.diff(axis=axis)
        self.resetXax(axis)
        self.addHistory("Differences over dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.diff(axis)))

    def cumsum(self, axis=-1):
        """
        The cumulative sum along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.data.cumsum(axis=axis)
        self.addHistory("Cumulative sum over dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.cumsum(axis)))

    def flipLR(self, axis=-1):
        """
        Flips the data along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        slicing = (slice(None), ) * axis + (slice(None, None, -1), )
        self.data = self.data[slicing]
        self.addHistory("Flipped dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.flipLR(axis))

    def hilbert(self, axis=-1):
        """
        Performs a Hilbert transform along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        if self.spec[axis] == 0:
            raise SpectrumException(f"Axis {axis+1} must be in frequency domain")
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.icomplexReorder(axis)
        self.data = self.data.hilbert(axis=axis)
        self.data.icomplexReorder(axis)
        self.addHistory("Hilbert transform on dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.hilbert(axis)))

    def autoPhaseAll(self, phaseNum=0, axis=-1):
        """
        Autophases all traces along a given axis individually.

        Parameters
        ----------
        phaseNum : {0, 1}, optional
            Order up to which to perform the autophasing.
            For 0 only zero order phasing is performed, for 1 both zero and first order phasing is performed.
            0 by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        shape = self.data.shape()
        shape = np.delete(shape, axis)
        rangeList = [range(i) for i in shape]
        for i in itertools.product(*rangeList):
            locList = np.insert(i, axis, 0)
            selectList = np.insert(np.array(i, dtype=object), axis, slice(None))
            self.autoPhase(phaseNum, axis, locList, False, selectList)
        if phaseNum == 1:
            self.addHistory("Autophased per trace for 0 + 1 order along axis " + str(axis + 1))
        else:
            self.addHistory("Autophased per trace for 0 order along axis " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.autoPhaseAll(phaseNum, axis)))

    def autoPhase(self, phaseNum=0, axis=-1, locList=None, returnPhases=False, select=slice(None)):
        """
        Autophases a spectrum.

        Parameters
        ----------
        phaseNum : {0, 1}, optional
            Order up to which to perform the autophasing.
            For 0 only zero order phasing is performed, for 1 both zero and first order phasing is performed.
            0 by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        locList : array_like of int, optional
            The indices of the trace to determine the phase values.
            By default the first index of each dimension is used.
        returnPhases : bool, optional
            If True the determined phases are return.
            False by default.
        select : Slice, optional
            An optional selection of the spectrum data on which the phasing is performed.
            By default the entire data is used.

        Raises
        ------
        SpectrumException
            When locList does not have the same length as the number of dimensions or when locList contains invalid indices.
        """
        axis = self.checkAxis(axis)
        if locList is None:
            locList = [0]*self.ndim()
        if len(locList) != self.ndim():
            raise SpectrumException("Data does not have the correct number of dimensions")
        if np.any(locList >= np.array(self.shape())) or np.any(np.array(locList) < 0):
            raise SpectrumException("The location array contains invalid indices")
        locList = np.array(locList, dtype=object)
        locList[axis] = slice(None)
        self.data.icomplexReorder(axis)
        if self.spec[axis] == 0:
            self.__fourier(axis, tmp=True)
        tmp = self.data[locList]
        tmp = tmp.getHyperData(0)   # only optimize on the hyper real data
        x = np.fft.fftshift(np.fft.fftfreq(len(tmp), 1.0 / self.sw[axis])) / self.sw[axis]
        if phaseNum == 0:
            phases = scipy.optimize.minimize(func.ACMEentropy, [0], (tmp, x, False), method='Powell', options={'xtol': AUTOPHASETOL})
            phase0 = phases['x']
            phase1 = 0.0
        elif phaseNum == 1:
            phases = scipy.optimize.minimize(func.ACMEentropy, [0, 0], (tmp, x), method='Powell', options={'xtol': AUTOPHASETOL})
            phase0 = phases['x'][0]
            phase1 = phases['x'][1]
        if self.ref[axis] is None:
            offset = 0
        else:
            offset = self.freq[axis] - self.ref[axis]
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.shape()[axis], 1.0 / self.sw[axis]) + offset) / self.sw[axis] * phase1 * 1j)
        vector = vector.reshape(vector.shape + (1, )*(self.ndim()-axis-1))
        self.data[select] *= np.exp(phase0 * 1j) * vector
        if self.spec[axis] == 0:
            self.__invFourier(axis, tmp=True)
        self.data.icomplexReorder(axis)
        Message = "Autophase: phase0 = " + str(phase0 * 180 / np.pi) + " and phase1 = " + str(phase1 * 180 / np.pi) + " for dimension " + str(axis + 1)
        self.addHistory(Message)
        if returnPhases:
            if not phases['x'].ndim:
                return phases['x'].reshape((1,))
            return phases['x']
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.phase(-phase0, -phase1, 0, axis))

    def __phase(self, phase0, phase1, phase2, offset, axis, select=slice(None)):
        points = np.fft.fftshift(np.fft.fftfreq(self.shape()[axis], 1.0 / self.sw[axis]) + offset) / self.sw[axis]
        vector = np.exp(points * phase1 * 1j + np.power(points, 2) * phase2 * 1j)
        if self.spec[axis] == 0:
            self.__fourier(axis, tmp=True)
        vector = vector.reshape(vector.shape + (1, )*(self.ndim()-axis-1))
        self.data.icomplexReorder(axis)
        self.data[select] *= np.exp(phase0 * 1j) * vector
        self.data.icomplexReorder(axis)
        if self.spec[axis] == 0:
            self.__invFourier(axis, tmp=True)

    def phase(self, phase0=0.0, phase1=0.0, phase2=0.0, axis=-1, select=slice(None), offset=None):
        """
        Phases a spectrum along a given dimension.

        Parameters
        ----------
        phase0 : float, optional
            Zero order phase.
            0.0 by default.
        phase1 : float, optional
            First order phase.
            0.0 by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        select : Slice, optional
            An optional selection of the spectrum data on which the phasing is performed.
            By default the entire data is used.
        offset : float, optional
            The offset frequency for the first order phase correction.
            When set to None, the offset is set to the reference frequency.
            None by default.
        """
        axis = self.checkAxis(axis)
        if offset is None:
            if self.ref[axis] is None:
                offset = 0
            else:
                offset = self.freq[axis] - self.ref[axis]
        self.__phase(phase0, phase1, phase2, offset, axis, select=select)
        Message = "Phasing: phase0 = " + str(phase0 * 180 / np.pi) + " and phase1 = " + str(phase1 * 180 / np.pi) + \
                  " and phase2 = " + str(phase2 * 180 / np.pi) + " for dimension " + str(axis + 1)
        if not isinstance(select, slice):
            Message = Message + " with slice " + str(select)
        elif select != slice(None, None, None):
            Message = Message + " with slice " + str(select)
        self.addHistory(Message)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.phase(-phase0, -phase1, -phase2, axis, select=select))

    def correctDFilter(self, axis=-1):
        """
        Corrects the digital filter via first order phasing along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if self.ref[axis] is None:
            offset = 0
        else:
            offset = self.freq[axis] - self.ref[axis]
        self.__phase(0, self.dFilter, 0, offset, axis)
        Message = "Corrected digital filter"
        self.addHistory(Message)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.phase(0, -self.dFilter, 0, axis, slice(None), offset))

    def apodize(self, lor=None, gauss=None, cos2=[None, None], hamming=None, shift=0.0, shifting=0.0, shiftingAxis=None, axis=-1, select=slice(None), preview=False):
        """
        Apodizes an FID or spectrum along a given dimension.

        Parameters
        ----------
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
        shift : float, optional
            A shift in time of the function.
            A positive value shift the curve to later times.
            By default a shift of 0.0 is used.
        shifting : float, optional
            A shift in time of the function as a function of the x-axis values along shiftingAxis.
            A positive value shift the curve to later times.
            By default a shifting of 0.0 is used.
        shiftingAxis : int, optional
            The dimension for the shifting.
            If this parameter is None, no shifting is applied.
            By default shifting is not used.
        axis : int, optional
            The dimension along which the apodization is applied.
            By default the last dimension is used.
        select : Slice, optional
            An optional selection of the spectrum data on which the apodization is performed.
            By default the entire data is used.
        """
        axis = self.checkAxis(axis)
        if shiftingAxis is None:
            shiftingAxis = 0
            shifting = 0.0
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axLen = self.shape()[axis]
        t = np.arange(0, axLen) / self.sw[axis]
        if shifting != 0.0:
            if self.spec[shiftingAxis]:
                shift1 = shift + shifting * np.arange(self.shape()[shiftingAxis]) / self.sw[shiftingAxis]
            else:
                shift1 = shift + shifting * self.xaxArray[shiftingAxis]
            previewData = np.array([func.apodize(t, s, lor, gauss, cos2, hamming, self.wholeEcho[axis]) for s in shift1])
            if axis < shiftingAxis:
                previewData = np.swapaxes(previewData, 0, 1)
            multShape = np.ones(len(self.shape()), dtype=int)
            multShape[axis] = self.shape()[axis]
            multShape[shiftingAxis] = self.shape()[shiftingAxis]
            previewData = previewData.reshape(multShape)
            if self.spec[axis] > 0:
                self.__invFourier(axis, tmp=True)
            if not isinstance(select, slice):
                previewSelect = np.full(len(previewData.shape), slice(None))
                previewSelect[shiftingAxis] = select[shiftingAxis]
                previewData = previewData[tuple(previewSelect)]
                select = tuple(select)
            self.data[select] *= previewData
            if self.spec[axis] > 0:
                self.__fourier(axis, tmp=True)
        else:
            x = func.apodize(t, shift, lor, gauss, cos2, hamming, self.wholeEcho[axis])
            if preview:
                previewData = [x] * int(np.prod(self.data.shape()) / self.data.shape()[axis])
            if self.spec[axis] > 0:
                self.__invFourier(axis, tmp=True)
            self.data[select] *= x.reshape(x.shape + (1, )*(self.ndim()-axis-1))
            if self.spec[axis] > 0:
                self.__fourier(axis, tmp=True)
        # Create the history message based on the input values.
        Message = "Apodization: "
        if lor is not None:
            Message = Message + "Lorentzian = " + str(lor) + ", "
        if gauss is not None:
            Message = Message + "Gaussian = " + str(gauss) + ", "
        if cos2[0] is not None:
            Message = Message + "Cos2 frequency = " + str(cos2[0]) + ", "
        if cos2[1] is not None:
            Message = Message + "Cos2 phase = " + str(cos2[1]) + ", "
        if hamming is not None:
            Message = Message + "Hamming = " + str(hamming) + ", "
        if shift != 0.0:
            Message = Message + "shift = " + str(shift) + ", "
        if shifting != 0.0:
            Message = Message + "shifting = " + str(shifting) + ", "
        if shiftingAxis != 0:
            Message = Message + "shiftingAxis = " + str(shiftingAxis) + ", "
        if lor is None and gauss is None and cos2 is None and hamming is None:  # If all none, make special message with `zero apodization'
            Message = Message + "zero apodization"
        Message = Message + " for dimension " + str(axis + 1)
        if not isinstance(select, slice):
            Message = Message + " with slice " + str(select)
        elif select != slice(None, None, None):
            Message = Message + " with slice " + str(select)
        self.addHistory(Message)
        if preview:
            return ([t]*len(previewData), previewData)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxis, axis, select=select)))

    def setFreq(self, freq=None, sw=None, axis=-1):
        """
        Sets the frequency or spectral width of a certain dimension to a given value.

        Parameters
        ----------
        freq : float, optional
            The new frequency of axis in Hz
            By default the frequency is unchanged.
        sw : float, optional
            The new spectral width of axis in Hz.
            By default the spectral width is unchanged.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        copyData = None
        if self.customXax[axis] and not self.noUndo:
            copyData = copy.deepcopy(self)
        oldFreq = self.freq[axis]
        oldSw = self.sw[axis]
        if freq is None:
            freq = self.freq[axis]
        if sw is None:
            sw = self.sw[axis]
        self.freq[axis] = float(freq)
        self.sw[axis] = float(sw)
        self.resetXax(axis)
        self.addHistory("Frequency set to " + str(freq * 1e-6) + " MHz and sw set to " + str(sw * 1e-3) + " kHz for dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            if copyData is not None:
                self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.setFreq(freq, sw, axis)))
            else:
                self.undoList.append(lambda self: self.setFreq(oldFreq, oldSw, axis))

    def scaleSw(self, scale, axis=-1):
        """
        Scales the spectral with of a dimension by a given scaling factor.

        Parameters
        ----------
        scale : float
            The scaling factor.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        copyData = None
        if self.customXax[axis] and not self.noUndo:
            copyData = copy.deepcopy(self)
        oldSw = self.sw[axis]
        self.sw[axis] = float(scale) * oldSw
        self.resetXax(axis)
        self.addHistory("Sw scaled by factor " + str(scale) +  " for dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            if copyData is not None:
                self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.scaleSw(scale, axis)))
            else:
                self.undoList.append(lambda self: self.scaleSw(1.0 / scale, axis))

    def setRef(self, ref=None, axis=-1):
        """
        Sets the reference frequency to a given value.

        Parameters
        ----------
        ref : float or None, optional
            The reference frequency in Hz.
            If None, the reference is removed.
            None by default.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        copyData = None
        if self.customXax[axis] and not self.noUndo:
            copyData = copy.deepcopy(self)
        oldRef = self.ref[axis]
        if ref is None:
            self.ref[axis] = None
            self.addHistory("Reference frequency set to 'None' for dimension " + str(axis + 1))
        else:
            self.ref[axis] = float(ref)
            self.addHistory("Reference frequency set to " + str(ref * 1e-6) + " MHz for dimension " + str(axis + 1))
        self.resetXax(axis)
        self.redoList = []
        if not self.noUndo:
            if copyData is not None:
                self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.setRef(ref, axis)))
            else:
                self.undoList.append(lambda self: self.setRef(oldRef, axis))

    def regrid(self, limits, numPoints, axis=-1):
        """
        Regrids the data along a dimension to match a given number of points.

        Parameters
        ----------
        limits : list of float
            The left and right limit of the new x-axis for the regrid.
        numPoints : int
            The number of points the data should have after the regrid.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        newSw = (limits[1] - limits[0]) / (numPoints - 1) * numPoints
        newAxis = np.fft.fftshift(np.fft.fftfreq(numPoints, 1.0 / newSw))
        newAxis = newAxis - (newAxis[0] + newAxis[-1]) / 2 + (limits[0] + limits[-1]) / 2  # Axis with correct min/max
        newFreq = self.freq[axis] + (newAxis[0] + newAxis[-1]) / 2
        if numPoints % 2 == 0:
            newFreq += newSw / numPoints / 2
        self.data = self.data.regrid(newAxis, self.xaxArray[axis], axis)
        self.sw[axis] = newSw
        if self.ref[axis] is None:  # Set new 0 freq to those of the old view, if needed
            self.ref[axis] = self.freq[axis]
        else:
            newFreq += - self.freq[axis] + self.ref[axis]
        self.freq[axis] = newFreq
        self.resetXax(axis)
        self.addHistory("Regrid dimension " + str(axis) + " between " + str(limits[0]) + ' and ' + str(limits[1]) + ' with ' + str(numPoints) + ' points')
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.regrid(limits, numPoints, axis)))

    def setWholeEcho(self, val, axis=-1):
        """
        Sets the Whole Echo attribute for a given dimension.

        Parameters
        ----------
        val : bool
            The new Whole Echo setting.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        self.wholeEcho[axis] = val
        self.addHistory("Whole echo set to " + str(val) + " for dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.setWholeEcho(not val, axis))

    def resize(self, size, pos, axis=-1):
        """
        Resizes the data along a dimension by zerofilling in the time domain.

        Parameters
        ----------
        size : int
            The new size along axis.
        pos : int
            The index after which the zeros are added in the time domain.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.spec[axis]:
            self.__invFourier(axis, tmp=True)
            
        oldsize = self.shape()[axis]
        self.data = self.data.resize(size, pos, axis=axis)
        if self.spec[axis]:
            self.__fourier(axis, tmp=True)
        
        if size <= oldsize and not self.spec[axis]:
            self.xaxArray[axis] = self.xaxArray[axis][:size]
        else:
            self.resetXax(axis)
        self.addHistory("Resized dimension " + str(axis + 1) + " to " + str(size) + " points at position " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.resize(size, pos, axis)))

    def lpsvd(self, nPredict, maxFreq, forward=False, numPoints=None, axis=-1):
        """
        Performs linear prediction using the LPSVD algorithm.

        Parameters
        ----------
        nPredict : int
            The number of datapoints to predict.
        maxFreq : int
            The maximum number of frequencies to take from the SVD.
        forward : bool, optional
            If True, a forward prediction is performed, otherwise a backward prediction is performed.
            False by default.
        numPoints : int, optional
            The number of points to use for SVD.
            For efficiency this number can be reduced.
            By default the entire data is used.
        axis : int, optional
            The dimension.
            By default the last dimension is used.

        Raises
        ------
        SpectrumException
            When the LPSVD algorithm resulted in an error.
        """
        failed = False
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.icomplexReorder(axis)
        if self.spec[axis]:
            self.__invFourier(axis, tmp=True)
        try:
            self.data = self.data.apply_along_axis(func.lpsvd, axis, nPredict, maxFreq, forward, numPoints)
        except Exception:
            failed = True
        if self.spec[axis]:
            self.__fourier(axis, tmp=True)
        self.data.icomplexReorder(axis)
        if failed:
            raise SpectrumException('LPSVD: Could not determine any acceptable values')
        self.resetXax(axis)
        if forward:
            self.addHistory("Forward LPSVD along axis " + str(axis) + " with " + str(nPredict) + " points, max " + str(maxFreq) + " frequencies")
        else:
            self.addHistory("Backward LPSVD along axis " + str(axis) + " with " + str(nPredict) + " points, max " + str(maxFreq) + " frequencies")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.lpsvd(nPredict, maxFreq, forward, numPoints, axis)))

    def setSpec(self, val, axis=-1):
        """
        Sets the spectrum flag for a given dimension.

        Parameters
        ----------
        val : bool
            If True, the dimension is treated as a spectral domain.
            If False, the dimension is treated as a time domain.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        copyData = None
        if self.customXax[axis] and not self.noUndo:
            copyData = copy.deepcopy(self)
        oldVal = self.spec[axis]
        self.spec[axis] = val
        self.resetXax(axis)
        if val:
            self.addHistory("Dimension " + str(axis + 1) + " set to FID")
        else:
            self.addHistory("Dimension " + str(axis + 1) + " set to spectrum")
        self.redoList = []
        if not self.noUndo:
            if copyData is not None:
                self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.setSpec(val, axis)))
            else:
                self.undoList.append(lambda self: self.setSpec(oldVal, axis))

    def swapEcho(self, idx, axis=-1):
        """
        Swaps an echo signal over a given index.

        Parameters
        ----------
        idx : int
            The index over which the data is swapped.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        slicing1 = (slice(None), ) * axis + (slice(None, idx), )
        slicing2 = (slice(None), ) * axis + (slice(idx, None), )
        self.data = self.data[slicing2].append(self.data[slicing1], axis=axis)
        self.wholeEcho[axis] = not self.wholeEcho[axis]
        self.addHistory("Swap echo at position " + str(idx) + " for dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.swapEcho(-idx, axis))

    def shift(self, shift, axis=-1, select=slice(None)):
        """
        Shifts the data a given number of datapoints.

        Parameters
        ----------
        shift : int
            The number of datapoints to shift.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        select : Slice, optional
            An optional selection of the spectrum data on which the shift is performed.
            By default the entire data is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.spec[axis] > 0:
            self.__invFourier(axis, tmp=True)
        mask = np.ones(self.shape()[axis])
        if shift < 0:
            mask[slice(shift, None)] = 0
        else:
            mask[slice(None, shift)] = 0
        self.data[select] = self.data.roll(shift, axis)[select]
        self.data[select] *= mask.reshape(mask.shape + (1,)*(self.ndim()-axis-1))
        if self.spec[axis] > 0:
            self.__fourier(axis, tmp=True)
        Message = "Shifted " + str(shift) + " points in dimension " + str(axis + 1)
        if not isinstance(select, slice):
            Message = Message + " with slice " + str(select)
        elif select != slice(None, None, None):
            Message = Message + " with slice " + str(select)
        self.addHistory(Message)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.shift(shift, axis, select=select)))

    def roll(self, shift, axis=-1, select=slice(None), shift_axis=True):
        """
        Rolls the data by applying a first order phase change in reciprocal space.
        This allows rolling also using non-integer values.

        Parameters
        ----------
        shift : float
            The distance to roll the data along axis. shift is in spectrum point unit but can be float.
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        select : Slice, optional
            An optional selection of the spectrum data on which the shift is performed.
            By default the entire data is used.
        shift_axis : bool, if True, on frequency domain, data and axis are shifted. 
                    If customXax is True, shift xaxArray, else change freq and resetXax
                    If select is not default value, shift_axis is not applied (set to False) 
                    shift_axis has no effect in time domain.
        """
        axis = self.checkAxis(axis)
        if select != slice(None, None, None):
            shift_axis = False
        x_shift = shift * self.sw[axis] / (self.shape()[axis])
        if self.spec[axis] == 0:
            self.__phase(0, -2 * np.pi * shift, 0, 0, axis, select=select)
        else:
            t = np.arange(self.shape()[axis]) / (self.sw[axis])
            freq_step = self.sw[axis] / self.shape()[axis]
            self.__invFourier(axis, tmp=True)
            t = t.reshape(t.shape + (1, )*(self.ndim()-axis-1))
            self.data.icomplexReorder(axis)
            self.data[select] *= np.exp(-1j * t * freq_step * shift * 2 * np.pi)
            self.data.icomplexReorder(axis)
            self.__fourier(axis, tmp=True)
            if self.customXax[axis] == True:
                self.xaxArray[axis] -= x_shift
            else:
                if shift_axis:
                    self.freq[axis] += x_shift
                    self.resetXax()
        Message = "Rolled " + str(shift) + " points in dimension " + str(axis + 1)
        if not isinstance(select, slice):
            Message = Message + " with slice " + str(select)
        elif select != slice(None, None, None):
            Message = Message + " with slice " + str(select)
        if shift_axis:
            Message = Message + f", axis also shifted by {x_shift}"
        else:
            Message = Message + ", axis not shifted"
        self.addHistory(Message)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.roll(-shift, axis, select, shift_axis))

    def align(self, pos1=None, pos2=None, axis=-1):
        """
        Aligns the maxima between given indices along a certain dimension.

        Parameters
        ----------
        pos1 : int
            The first index between which to determine the maximum.
        pos2 : int
            The second index between which to determine the maximum.
        axis : int, optional
            The dimension.
            By default the last dimension is used.

        Raises
        ------
        SpectrumException
            When pos1 or pos2 is invalid.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if pos1 is None:
            pos1 = 0
        if pos2 is None:
            pos2 = self.shape()[axis]
        if not 0 <= pos1 <= self.shape()[axis]:
            raise SpectrumException("Indices not within range")
        if not 0 <= pos2 <= self.shape()[axis]:
            raise SpectrumException("Indices not within range")
        if pos1 == pos2:
            raise SpectrumException("Indices cannot be equal")
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axis + (slice(minPos, maxPos), )
        tmp = self.data[slicing].argmax(axis=axis)
        maxArgPos = -np.array(tmp.data, dtype=int)
        maxArgPos -= maxArgPos.flatten()[0]
        shape = self.data.shape()
        shape = np.delete(shape, axis)
        rangeList = [range(i) for i in shape]
        for i in itertools.product(*rangeList):
            selectList = np.insert(np.array(i, dtype=object), axis, slice(None))
            self.data[selectList] = self.data[selectList].roll(maxArgPos[0][tuple(i)], 0)
        self.addHistory("Maxima aligned between " + str(minPos) + " and " + str(maxPos) + " along axis " + str(axis))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.align(pos1, pos2, axis)))

    def __fourier(self, axis, tmp=False, reorder=None):
        axis = self.checkAxis(axis)
        if reorder is None:
            reorder = [True, True]
        if reorder[0]:
            self.data.icomplexReorder(axis)
        if not self.wholeEcho[axis] and not tmp:
            slicing = (slice(None), ) * axis + (0, )
            self.data[slicing] = self.data[slicing] * 0.5
        self.data = self.data.fft(axis).fftshift(axis)
        if not tmp:
            self.spec[axis] = 1
        if reorder[1]:
            self.data.icomplexReorder(axis)
        self.resetXax(axis)

    def __invFourier(self, axis, tmp=False, reorder=None):
        axis = self.checkAxis(axis)
        if reorder is None:
            reorder = [True, True]
        if reorder[0]:
            self.data.icomplexReorder(axis)
        self.data = self.data.ifftshift(axis).ifft(axis)
        if not self.wholeEcho[axis] and not tmp:
            slicing = (slice(None), ) * axis + (0, )
            self.data[slicing] *= 2.0
        if not tmp:
            self.spec[axis] = 0
        if reorder[1]:
            self.data.icomplexReorder(axis)
        self.resetXax(axis)

    def complexFourier(self, axis=-1):
        """
        Perform a complex Fourier transform along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        if self.spec[axis] == 0:
            self.__fourier(axis)
            self.addHistory("Fourier transform dimension " + str(axis + 1))
        else:
            self.__invFourier(axis)
            self.addHistory("Inverse fourier transform dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.complexFourier(axis))

    def realFourier(self, axis=-1):
        """
        Perform a real Fourier transform along a given dimension.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        """
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axis = self.checkAxis(axis)
        self.data = self.data.real(axis)
        if self.spec[axis] == 0:
            self.__fourier(axis)
            self.addHistory("Real fourier transform dimension " + str(axis + 1))
        else:
            self.__invFourier(axis)
            self.addHistory("Inverse real fourier transform dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.realFourier(axis)))

    def fftshift(self, axis=-1, inv=False):
        """
        Perform an fftshift or inverse fftshift along a given dimension.
        This shifts the zero-frequency component to the center of the spectrum.

        Parameters
        ----------
        axis : int, optional
            The dimension.
            By default the last dimension is used.
        inv : bool, optional
            If True, the inverse fftshift is performed.
            False by default .
        """
        axis = self.checkAxis(axis)
        if inv:
            self.data = self.data.ifftshift(axis=axis)
            self.addHistory("Inverse Fourier shift dimension " + str(axis + 1))
        else:
            self.data = self.data.fftshift(axis=axis)
            self.addHistory("Fourier shift dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.fftshift(axis, not(inv)))

    def shear(self, shear, axis=-1, axis2=-2, toRef=False):
        """
        Shears the data along given axes.

        Parameters
        ----------
        shear : float
            The shearing factor defined in ratio between the x-axes of axis and axis2.
        axis : int, optional
            The dimension along which the shearing is applied.
            By default the last dimension is used.
        axis2 : int, optional
            The shearing axis.
            By default the second last dimension is used.
        toRef : bool, optional
            If True, the shearing is applied relative to the reference frequencies of both dimensions.
            Otherwise the shearing is applied relative to the center of the data.
            False by default.

        Raises
        ------
        SpectrumException
            When axis2 equals axis or when the data does not have more one dimension.
        """
        axis = self.checkAxis(axis)
        axis2 = self.checkAxis(axis2)
        if axis == axis2:
            raise SpectrumException('Both shearing axes cannot be equal')
        if self.ndim() < 2:
            raise SpectrumException("The data does not have enough dimensions for a shearing transformation")
        if self.spec[axis] > 0: #rorder and fft for spec
            self.__invFourier(axis, tmp=True, reorder=[True, False])
        else: #Reorder if FID
            self.data.icomplexReorder(axis)
        shape = self.shape()
        if toRef:
            vec2 = self.xaxArray[axis2]
        else:
            vec2 = np.fft.fftshift(np.fft.fftfreq(shape[axis2], 1 / self.sw[axis2]))
        vec1 = np.linspace(0, shear * 2 * np.pi * shape[axis] / self.sw[axis], shape[axis] + 1)[:-1]
        newShape = [1, ] * self.ndim()
        newShape[axis] = shape[axis]
        newShape[axis2] = shape[axis2]
        if axis > axis2:
            shearMatrix = np.exp(1j * np.outer(vec2, vec1))
        elif axis < axis2:
            shearMatrix = np.exp(1j * np.outer(vec1, vec2))
        self.data *= shearMatrix.reshape(shape)
        if self.spec[axis] > 0:
            self.__fourier(axis, tmp=True, reorder=[False, True])
        else:
            self.data.icomplexReorder(axis)
        self.addHistory("Shearing transform with shearing value " + str(shear) + " over dimensions " + str(axis + 1) + " and " + str(axis2 + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.shear(-shear, axis, axis2, toRef))

    def reorder(self, pos, newLength=None, axis=-1):
        """
        Reorders the data based on a given list of indices.

        Parameters
        ----------
        pos : array_like
            The positions of the subarrays in the new data.
        newLength : int, optional
            The new length of the data along axis.
            By default one more than the largest value in pos is used.
        axis : int, optional
            The axis along which the data is reordered.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.data.reorder(pos, newLength, axis)
        self.resetXax(axis)
        self.addHistory("Reorder dimension " + str(axis + 1) + " to obtain a new length of " + str(newLength) + " with positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.reorder(pos, newLength, axis)))

    def ffm(self, pos, typeVal, axis=-1):
        """
        Uses the fast forward maximum entropy algorithm to reconstruct non-uniform sampled data.

        Parameters
        ----------
        pos : array_like
            A list of indices that are recorded datapoints.
            All other datapoints will be reconstructed.
        typeVal : {0, 1, 2}
            The type of data to be reconstructed.
            0=complex, 1=States or States-TPPI, 2=TPPI.
        axis : int, optional
            The axis along which the data is reconstructed.
            By default the last dimension is used.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        posList = np.delete(range(self.shape()[axis]), pos)  # pos contains the values of fixed points which not to be translated to missing points
        if typeVal == 1:  # type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList / 2), dtype=int)
        if typeVal == 2:  # type is TPPI, for now handle the same as Complex
            pass
        posList = np.unique(posList)
        self.data.icomplexReorder(axis)
        tmpData = self.data.getHyperData(0)
        tmpData = np.rollaxis(tmpData, axis, tmpData.ndim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((int(tmpData.size / tmpShape[-1]), tmpShape[-1]))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(nus.ffm, [(i, posList) for i in tmpData])
        pool.close()
        pool.join()
        tmpData = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axis)
        self.data = hc.HComplexData(tmpData)
        self.__invFourier(axis, tmp=True)  # Transform back to FID
        self.addHistory("Fast Forward Maximum Entropy reconstruction of dimension " + str(axis + 1) + " at positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, None))

    def clean(self, pos, typeVal, axis, gamma, threshold, maxIter):
        """
        Uses the CLEAN algorithm to reconstruct non-uniform sampled data.

        Parameters
        ----------
        pos : array_like
            A list of indices that are recorded datapoints.
            All other datapoints will be reconstructed.
        typeVal : {0, 1, 2}
            The type of data to be reconstructed.
            0=complex, 1=States or States-TPPI, 2=TPPI.
        axis : int
            The axis along which the data is reconstructed.
        gamma : float
            Gamma value of the CLEAN calculation.
        threshold : float
            Stopping limit (0 < x < 1) (stop if residual intensity below this point).
        maxIter : int
            Maximum number of iterations.
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        posList = np.delete(range(self.shape()[axis]), pos)  # pos contains the values of fixed points which not to be translated to missing points
        if typeVal == 1:  # type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList / 2), dtype=int)
        if typeVal == 2:  # type is TPPI, for now handle the same as Complex
            pass
        posList = np.unique(posList)
        self.data.icomplexReorder(axis)
        tmpData = self.data.getHyperData(0)
        tmpData = np.rollaxis(np.fft.fft(tmpData, axis=axis), axis, tmpData.ndim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((int(tmpData.size / tmpShape[-1]), tmpShape[-1]))
        mask = np.ones(tmpShape[-1]) / float(tmpShape[-1])
        mask[posList] = 0.0
        mask = np.fft.fft(mask) # abs or real???
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(nus.clean, [(i, mask, gamma, threshold, maxIter) for i in tmpData])
        pool.close()
        pool.join()
        tmpData = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axis)
        self.data = hc.HComplexData(tmpData)
        self.__invFourier(axis, tmp=True)  # Transform back to FID
        self.addHistory("CLEAN reconstruction (gamma = " + str(gamma) + " , threshold = " + str(threshold) + " , maxIter = " + str(maxIter) + ") " + "of dimension " + str(axis + 1) + " at positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, None))

    def ist(self, pos, typeVal, axis, threshold, maxIter, tracelimit):
        """
        Uses the Iterative Soft Thresholding algorithm to reconstruct non-uniform sampled data.

        Parameters
        ----------
        pos : array_like
            A list of indices that are recorded datapoints.
            All other datapoints will be reconstructed.
        typeVal : {0, 1, 2}
            The type of data to be reconstructed.
            0=complex, 1=States or States-TPPI, 2=TPPI.
        axis : int
            The axis along which the data is reconstructed.
        threshold : float
            threshold. The level (0 < x < 1) at which the data is cut every iteration.
        maxIter : int
            Maximum number of iterations.
        tracelimit : float
            Stopping limit (0 < x < 1) (stop if residual intensity below this point).
        """
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.icomplexReorder(axis)
        tmpData = self.data.getHyperData(0)
        posList = np.delete(range(tmpData.shape[axis]), pos)  # pos contains the values of fixed points which not to be translated to missing points
        if typeVal == 1:  # type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList / 2), dtype=int)
        elif typeVal == 2:  # type is TPPI, for now handle the same as Complex
            pass
        NDmax = np.max(np.max(np.abs(np.real(np.fft.fft(tmpData, axis=axis))))) #Get max of ND matrix
        tmpData = np.rollaxis(tmpData, axis, tmpData.ndim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((int(tmpData.size / tmpShape[-1]), tmpShape[-1]))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(nus.ist, [(i, posList, threshold, maxIter, tracelimit, NDmax) for i in tmpData])
        pool.close()
        pool.join()
        tmpData = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axis)
        self.data = hc.HComplexData(tmpData)
        self.__invFourier(axis, tmp=True)  # Transform back to FID
        self.addHistory("IST reconstruction (threshold = " + str(threshold) + " , maxIter = " + str(maxIter) + " , tracelimit = " + str(tracelimit*100) + ") " + "of dimension " + str(axis + 1) + " at positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, None))

    def getSlice(self, axes, locList, stack=None):
        """
        Generate a new Spectrum object which is a subset of this object.

        Parameters
        ----------
        axes : array_like
            A list with the axes desired in the new object.
        locList : array_like
            The indices to specify the desired data subset.
        stack : array_like, optional
            Slice objects to extract part of an axis.
            Often used for stack or array plots.
            By default the entire data along all axes are used.

        Returns
        -------
        Spectrum
            The subset spectrum.
        """
        locList = np.array(locList, dtype=object)
        if stack is None:
            stack = [slice(None)]*len(axes)
        for i, axis in enumerate(axes):
            axes[i] = self.checkAxis(axis)
        locList[axes] = stack
        sliceData = self.data[locList].copy()
        sliceData.icomplexReorder(axes[-1])
        orderInd = np.argsort(axes)
        sliceData = sliceData.moveaxis(np.arange(sliceData.ndim()), orderInd)
        sliceSpec = copy.deepcopy(Spectrum(hc.HComplexData(sliceData.getHyperData(0)),
                                           self.filePath,
                                           [self.freq[axis] for axis in axes],
                                           [self.sw[axis] for axis in axes],
                                           spec=[self.spec[axis] for axis in axes],
                                           wholeEcho=[self.wholeEcho[axis] for axis in axes],
                                           ref=[self.ref[axis] for axis in axes],
                                           xaxArray=[self.xaxArray[axis][stack[i]] for i, axis in enumerate(axes)],
                                           customXax=[self.customXax[axis] for axis in axes],
                                           history=self.history,
                                           name=self.name))
        sliceSpec.noUndo = True
        return sliceSpec

    def restoreData(self, copyData, returnValue):
        """
        Restore data from an old copy.
        Often used for undo purposes.

        Parameters
        ----------
        copyData : Spectrum
            The old Spectrum object to restore from.
        returnValue
            A return value that should be appended to the undolist.
        """
        if (not self.noUndo) and returnValue is None:
            copyData2 = copy.deepcopy(self)
        self.data = copyData.data
        self.freq = copyData.freq  # array of center frequency (length is dim, MHz)
        self.filePath = copyData.filePath
        self.sw = copyData.sw  # array of sweepwidths
        self.spec = copyData.spec
        self.wholeEcho = copyData.wholeEcho
        self.xaxArray = copyData.xaxArray
        self.customXax = copyData.customXax
        self.ref = copyData.ref
        self.addHistory("Data was restored to a previous state ")
        self.redoList = []
        if (not self.noUndo) and returnValue is None:
            self.undoList.append(lambda self: self.restoreData(copyData2, None))
        else:
            self.undoList.append(returnValue)
