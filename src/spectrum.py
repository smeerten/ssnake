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
import scipy.optimize
import copy
from nus import ffm, clean, ist
import multiprocessing
import reimplement as reim
import functions as func
import hypercomplex as hc
import specIO as io

AUTOPHASETOL = 0.0002 #is ~0.01 degrees


class SpectrumException(Exception):
    pass


#########################################################################
# the generic spectrum class


class Spectrum(object):

    def __init__(self, data, filePath, freq, sw, spec=None, wholeEcho=None, ref=None, xaxArray=None, history=None, metaData = None, name=''):
        self.name = name
        if isinstance(data, hc.HComplexData):
            self.data = data
        else:
            self.data = hc.HComplexData(data)
        self.filePath = filePath
        self.freq = np.array(freq)  # array of center frequency (length is dim, MHz)
        self.sw = np.array(sw,dtype=float)  # array of sweepwidths
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
        if history is None:
            self.history = []  # list of strings describing all performed operations
        else:
            self.history = history
        if metaData is None:
            self.metaData = {'# Scans': '-', 'Acquisition Time [s]': '-', 'Experiment Name': '-','Receiver Gain': '-', 'Recycle Delay [s]': '-',
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
        self.name = name

    def getHistory(self):
        return "\n".join(self.history)

    def addHistory(self, msg):
        self.history.append(msg)

    def removeFromHistory(self, num):
        for i in range(num):
            if len(self.history) > 0:
                val = self.history.pop()
        return val

    def setNoUndo(self,val):
        self.noUndo = bool(val)
        if self.noUndo:
            self.undoList = []
            self.redoList = []

    def undo(self):
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
        if not self.redoList:
            raise SpectrumException("No redo information")
        tmpRedo = self.redoList # Protect the redo list
        tmpRedo.pop()(self)
        self.redoList = tmpRedo # Restore the redo list

    def clearUndo(self):
        self.undoList = []
        self.redoList = []

    def reload(self):
        loadData = io.autoLoad(*self.filePath)
        self.restoreData(loadData, None)
        
    def checkAxis(self, axis):
        if axis < 0:
            axis += self.ndim()
        if not (0 <= axis < self.ndim()):
            raise IndexError("Not a valid axis for Spectrum")
        return axis

    def resetXax(self, axis=None):
        if axis is not None:
            axis = self.checkAxis(axis)
            val = [axis]
        else:
            val = range(self.ndim())
        for i in val:
            if self.spec[i] == 0:
                self.xaxArray[i] = np.arange(self.shape()[i]) / (self.sw[i])
            elif self.spec[i] == 1:
                self.xaxArray[i] = np.fft.fftshift(np.fft.fftfreq(self.shape()[i], 1.0 / self.sw[i]))
                if self.ref[i] is not None:
                    self.xaxArray[i] += self.freq[i] - self.ref[i]

    def setXax(self, xax, axis):
        axis = self.checkAxis(axis)
        if len(xax) != self.shape()[axis]:
            raise SpectrumException("Length of new x-axis does not match length of the data")
        oldXax = self.xaxArray[axis]
        self.xaxArray[axis] = xax
        self.addHistory("X-axis of dimension " + str(axis + 1) + " was set to " + str(xax).replace('\n', ''))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.setXax(oldXax, axis))

    def insert(self, data, pos, axis):
        if not isinstance(data, hc.HComplexData):
            data = hc.HComplexData(data)
        if self.noUndo:
            returnValue = None
        else:
            if self.data.hyper == data.hyper: # If both sets have same hyper: easy undo can be used
                returnValue = lambda self: self.delete(range(pos, pos + data.shape()[axis]), axis)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.insert(data, pos, axis))
        axis = self.checkAxis(axis)
        # Check for a change in dimensions
        self.data = self.data.insert(pos, data, axis)
        self.resetXax(axis)
        self.addHistory("Inserted " + str(data.shape()[axis]) + " datapoints in dimension " + str(axis + 1) + " at position " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(returnValue)

    def delete(self, pos, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        tmpData = self.data.delete(pos, axis)
        if 0 in tmpData.shape():
            raise SpectrumException('Cannot delete all data')
        self.data = tmpData
        self.xaxArray[axis] = np.delete(self.xaxArray[axis], pos)
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
        if isinstance(data, Spectrum):
            data = data.data
        if isinstance(data, np.ndarray) and axis is not None:
            data = data.reshape(data.shape + (1,)*(self.ndim() - axis - 1))
        if not self.noUndo:
            if not isinstance(data, hc.HComplexData):
                returnValue = lambda self: self.subtract(data, axis, select=select)
            elif self.data.hyper == data.hyper: # If both sets have same hyper: easy subtract can be used for undo
                returnValue = lambda self: self.subtract(data, axis, select=select)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.add(data, axis, select))
        self.data[select] += data
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
        if isinstance(data, Spectrum):
            data = data.data
        if isinstance(data, np.ndarray) and axis is not None:
            data = data.reshape(data.shape + (1,)*(self.ndim() - axis - 1))
        if not self.noUndo:
            if not isinstance(data, hc.HComplexData):
                returnValue = lambda self: self.add(data, axis, select=select)
            elif self.data.hyper == data.hyper: #If both sets have same hyper: easy subtract can be used for undo
                returnValue = lambda self: self.add(data, axis, select=select)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.subtract(data, axis, select))
        self.data[select] -= data
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
        if isinstance(data, Spectrum):
            data = data.data
        if isinstance(data, np.ndarray) and axis is not None:
            data = data.reshape(data.shape + (1,)*(self.ndim() - axis - 1))
        if not self.noUndo:
            if not isinstance(data, hc.HComplexData):
                returnValue = lambda self: self.divide(data, axis, select=select)
            elif self.data.hyper == data.hyper: #If both sets have same hyper: easy subtract can be used for undo
                returnValue = lambda self: self.divide(data, axis, select=select)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.multiply(data, axis, select))
        self.data[select] *= data
        self.addHistory("Multiplied data[" + str(select) + "]")
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
        if isinstance(data, Spectrum):
            data = data.data
        if isinstance(data, np.ndarray) and axis is not None:
            data = data.reshape(data.shape + (1,)*(self.ndim() - axis - 1))
        if not self.noUndo:
            if not isinstance(data, hc.HComplexData):
                returnValue = lambda self: self.multiply(data, axis, select=select)
            elif self.data.hyper == data.hyper: #If both sets have same hyper: easy subtract can be used for undo
                returnValue = lambda self: self.multiply(data, axis, select=select)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.divide(data, axis, select))
        self.data[select] /= data
        self.addHistory("Divided by data[" + str(select) + "]")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(returnValue)

    def normalize(self, mult, scale, type, axis, select=slice(None)):
        axis = self.checkAxis(axis)
        try:
            self.data *= mult * scale 
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

    def baselineCorrection(self, baseline, axis, select=slice(None)):
        axis = self.checkAxis(axis)
        baselinetmp = baseline.reshape((self.shape()[axis], ) + (1, ) * (self.ndim() - axis - 1))
        self.data[select] -= baselinetmp
        Message = "Baseline corrected dimension " + str(axis + 1)
        if type(select) is not slice:
            Message = Message + " with slice " + str(select)
        elif select != slice(None, None, None):
            Message = Message + " with slice " + str(select)
        self.addHistory(Message)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.baselineCorrection(-baseline, axis, select=select))

    def concatenate(self, axis):
        axis = self.checkAxis(axis)
        splitVal = self.shape()[axis]
        copyData = None
        if self.data.isComplex(axis):
            if not self.noUndo:
                copyData = copy.deepcopy(self)
            self.data = self.data.real(axis)
        invAxis = self.ndim() - axis
        self.data = self.data.concatenate(axis)
        self.data.removeDim(invAxis)
        self.freq = np.delete(self.freq, axis)
        self.sw = np.delete(self.sw, axis)
        self.spec = np.delete(self.spec, axis)
        self.wholeEcho = np.delete(self.wholeEcho, axis)
        self.ref = np.delete(self.ref, axis)
        del self.xaxArray[axis]
        self.resetXax()
        self.addHistory("Concatenated dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            if copyData is not None:
                self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.concatenate(axis)))
            else:
                self.undoList.append(lambda self: self.split(splitVal, axis))

    def split(self, sections, axis):
        axis = self.checkAxis(axis)
        self.data = self.data.split(sections, axis)
        self.data.insertDim(0)
        self.freq = np.insert(self.freq, 0, self.freq[axis])
        self.sw = np.insert(self.sw, 0, self.sw[axis])
        self.spec = np.insert(self.spec, 0, self.spec[axis])
        self.wholeEcho = np.insert(self.wholeEcho, 0, self.wholeEcho[axis])
        self.ref = np.insert(self.ref, 0, self.ref[axis])
        self.xaxArray.insert(0, [])
        self.resetXax(0)
        self.resetXax(axis + 1)
        self.addHistory("Split dimension " + str(axis + 1) + " into " + str(sections) + " sections")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.concatenate(axis))

    def real(self, axis=-1):
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axis = self.checkAxis(axis)
        self.data = self.data.real(axis)
        self.addHistory("Real along dimension " + str(axis+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.real(axis)))

    def imag(self, axis=-1):
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axis = self.checkAxis(axis)
        self.data = self.data.imag(axis)
        self.addHistory("Imaginary along dimension " + str(axis+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.imag(axis)))

    def abs(self, axis=-1):
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axis = self.checkAxis(axis)
        self.data = self.data.abs(axis)
        self.addHistory("Absolute along dimension " + str(axis+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.abs(axis)))
    
    def conj(self, axis=-1):
        self.data = self.data.conj(axis)
        self.addHistory("Complex conjugate along" + str(axis+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.conj(axis))

    def states(self, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.states(axis)
        self.resetXax(axis)
        self.addHistory("States conversion on dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.states(axis)))

    def statesTPPI(self, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.states(axis, TPPI=True)
        self.resetXax(axis)
        self.addHistory("States-TPPI conversion on dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.statesTPPI(axis)))

    def echoAntiEcho(self, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.echoAntiEcho(axis)
        self.resetXax(axis)
        self.addHistory("Echo-antiecho conversion on dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.echoAntiEcho(axis)))

    def subtractAvg(self, pos1, pos2, axis):
        axis = self.checkAxis(axis)
        if not (0 <= pos1 <= self.shape()[axis]):
            raise SpectrumException("Indices not within range")
        if not (0 <= pos2 <= self.shape()[axis]):
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

    def matrixManip(self, pos1, pos2, axis, which):
        axis = self.checkAxis(axis)
        if isinstance(pos1, int):
            pos1 = np.array([pos1])
            pos2 = np.array([pos2])
        if len(pos1) != len(pos2):
            raise SpectrumException("Length of the two arrays is not equal")
        if len(pos1) == 1:
            if self.ndim() == 1:
                keepdims = True
            else:
                keepdims = False
        else:
            keepdims = True
        tmpdata = ()
        for i in range(len(pos1)):
            if not (0 <= pos1[i] <= self.shape()[axis]):
                raise SpectrumException("Indices not within range")
            if not (0 <= pos2[i] <= self.shape()[axis]):
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
                self.data.removeDim(axis)
        else:
            self.data = tmpdata[0]
            for extra in tmpdata[1:]:
                self.data = self.data.append(extra, axis=axis)
            self.resetXax(axis)

    def integrate(self, pos1, pos2, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=0)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.integrate(pos1, pos2, axis)))
        self.addHistory("Integrate between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def max(self, pos1, pos2, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=1)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.max(pos1, pos2, axis)))
        self.addHistory("Maximum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def min(self, pos1, pos2, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=2)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.min(pos1, pos2, axis)))
        self.addHistory("Minimum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def argmax(self, pos1, pos2, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=3)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.argmax(pos1, pos2, axis)))
        self.addHistory("Maximum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def argmin(self, pos1, pos2, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=4)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.argmin(pos1, pos2, axis)))
        self.addHistory("Minimum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def sum(self, pos1, pos2, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=5)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.sum(pos1, pos2, axis)))
        self.addHistory("Sum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def average(self, pos1, pos2, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axis, which=6)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.average(pos1, pos2, axis)))
        self.addHistory("Average between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axis + 1))

    def extract(self, pos1, pos2, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axis + (slice(minPos, maxPos), )
        if self.spec[axis] == 1:
            oldFxax = self.xaxArray[axis][slice(minPos, maxPos)][0]
            self.sw[axis] *= (maxPos - minPos) / (1.0 * self.shape()[axis])
        self.data = self.data[slicing]
        if self.spec[axis] == 1:
            newFxax = np.fft.fftshift(np.fft.fftfreq(self.shape()[axis], 1.0 / self.sw[axis]))[0]
            if self.ref[axis] is None:
                self.ref[axis] = self.freq[axis]
            self.freq[axis] = self.ref[axis] - newFxax + oldFxax
        self.resetXax(axis)
        self.addHistory("Extracted part between " + str(minPos) + " and " + str(maxPos) + " of dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.extract(pos1, pos2, axis)))

    def fiddle(self, refSpec, lb, axis):
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

    def diff(self, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.data.diff(axis=axis)
        self.resetXax(axis)
        self.addHistory("Differences over dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.diff(axis)))

    def cumsum(self, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.data.cumsum(axis=axis)
        self.addHistory("Cumulative sum over dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.cumsum(axis)))

    def flipLR(self, axis):
        axis = self.checkAxis(axis)
        slicing = (slice(None), ) * axis + (slice(None, None, -1), )
        self.data = self.data[slicing]
        self.addHistory("Flipped dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.flipLR(axis))

    def hilbert(self, axis):
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

    def autoPhase(self, phaseNum, axis, locList, returnPhases=False):
        axis = self.checkAxis(axis)
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
        tmp = tmp.getHyperData(0)
        x = np.fft.fftshift(np.fft.fftfreq(len(tmp), 1.0 / self.sw[axis])) / self.sw[axis]
        # only optimize on the hyper real data
        if phaseNum == 0:
            phases = scipy.optimize.minimize(self.ACMEentropy, [0], (tmp, x, False), method='Powell',options = {'xtol': AUTOPHASETOL})
            phase0 = phases['x']
            phase1 = 0.0
        elif phaseNum == 1:
            phases = scipy.optimize.minimize(self.ACMEentropy, [0, 0], (tmp, x), method='Powell', options = {'xtol': AUTOPHASETOL})
            phase0 = phases['x'][0]
            phase1 = phases['x'][1]
        if self.ref[axis] is None:
            offset = 0
        else:
            offset = self.freq[axis] - self.ref[axis]
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.shape()[axis], 1.0 / self.sw[axis]) + offset) / self.sw[axis] * phase1 * 1j)
        vector = vector.reshape(vector.shape + (1, )*(self.ndim()-axis-1))
        self.data *= np.exp(phase0 * 1j) * vector
        if self.spec[axis] == 0:
            self.__invFourier(axis, tmp=True)
        self.data.icomplexReorder(axis)
        Message = "Autophase: phase0 = " + str(phase0 * 180 / np.pi) + " and phase1 = " + str(phase1 * 180 / np.pi) + " for dimension " + str(axis + 1)
        self.addHistory(Message)
        if returnPhases:
            if phaseNum == 0:
                return [phases['x']]
            else:
                return phases['x']
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.phase(-phase0, -phase1, axis))

    def ACMEentropy(self, phaseIn, data, x, phaseAll=True):
        phase0 = phaseIn[0]
        if phaseAll:
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
        if (np.real(sumas) < 0):
            Pfun = Pfun + sum(as1**2) / 4 / L**2
        return H1 + 1000 * Pfun

    def phase(self, phase0, phase1, axis, select=slice(None)):
        axis = self.checkAxis(axis)
        if self.ref[axis] is None:
            offset = 0
        else:
            offset = self.freq[axis] - self.ref[axis]
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.shape()[axis], 1.0 / self.sw[axis]) + offset) / self.sw[axis] * phase1 * 1j)
        if self.spec[axis] == 0:
            self.__fourier(axis, tmp=True)
        vector = vector.reshape(vector.shape + (1, )*(self.ndim()-axis-1))
        self.data.icomplexReorder(axis)
        self.data[select] *= np.exp(phase0 * 1j) * vector
        self.data.icomplexReorder(axis)
        if self.spec[axis] == 0:
            self.__invFourier(axis, tmp=True)
        Message = "Phasing: phase0 = " + str(phase0 * 180 / np.pi) + " and phase1 = " + str(phase1 * 180 / np.pi) + " for dimension " + str(axis + 1)
        if type(select) is not slice:
            Message = Message + " with slice " + str(select)
        elif select != slice(None, None, None):
            Message = Message + " with slice " + str(select)

        self.addHistory(Message)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.phase(-phase0, -phase1, axis, select=select))

    def apodize(self, lor=None, gauss=None, cos2=[None, None], hamming=None, shift=0.0, shifting=0.0, shiftingAxis=None, axis=-1, select=slice(None), preview=False):
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
            previewData = np.array([func.apodize(t, s, self.sw[axis], axLen, lor, gauss, cos2, hamming, self.wholeEcho[axis]) for s in shift1])
            if axis < shiftingAxis:
                previewData = np.swapaxis(previewData, 0, 1)
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
            x = func.apodize(t, shift, self.sw[axis], axLen, lor, gauss, cos2, hamming, self.wholeEcho[axis])
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
        if type(select) is not slice:
            Message = Message + " with slice " + str(select)
        elif select != slice(None, None, None):
            Message = Message + " with slice " + str(select)
        self.addHistory(Message)
        if preview:
            return ([t]*len(previewData), previewData)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxis, axis, select=select)))

    def setFreq(self, freq, sw, axis):
        axis = self.checkAxis(axis)
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
            self.undoList.append(lambda self: self.setFreq(oldFreq, oldSw, axis))

    def scaleSw(self, scale, axis):
        axis = self.checkAxis(axis)
        oldSw = self.sw[axis]
        self.sw[axis] = float(scale) * oldSw
        self.resetXax(axis)
        self.addHistory("Sw scaled by factor " + str(scale) +  " for dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.scaleSw(1.0 / scale, axis))

    def setRef(self, ref, axis):
        axis = self.checkAxis(axis)
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
            self.undoList.append(lambda self: self.setRef(oldRef, axis))

    def regrid(self, limits, numPoints, axis):
        oldLimits = [self.xaxArray[axis][0], self.xaxArray[axis][-1]]
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

    def setWholeEcho(self, val, axis):
        axis = self.checkAxis(axis)
        self.wholeEcho[axis] = val
        self.addHistory("Whole echo set to " + str(val) + " for dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.setWholeEcho(not val, axis))

    def resize(self, size, pos, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.spec[axis] > 0:
            self.__invFourier(axis, tmp=True)
        self.data = self.data.resize(size, pos, axis=axis)
        if self.spec[axis] > 0:
            self.__fourier(axis, tmp=True)
        self.resetXax(axis)
        self.addHistory("Resized dimension " + str(axis + 1) + " to " + str(size) + " points at position " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.resize(size, pos, axis)))

    def lpsvd(self, nAnalyse, nFreq, nPredict, Direction, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.apply_along_axis(self.lpsvdfunction, axis, nAnalyse, nFreq, nPredict, Direction)
        self.resetXax(axis)
        self.addHistory("LPSVD ")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.lpsvd(nAnalyse, nFreq, nPredict, Direction, axis)))

    def lpsvdfunction(self, data, nAnalyse, nFreq, nPredict, Direction):
        # LPSVD algorithm
        if Direction == 1:  # If backward
            Y = data[0:nAnalyse]
        else:
            Y = data[-nAnalyse:]
        N = len(Y)						# # of complex data points in FID
        L = int(np.floor(N * 3 / 4))						# linear prediction order L = 3/4*N
        A = scipy.linalg.hankel(np.conj(Y[1:N - L + 1]), np.conj(Y[N - L:N]))  # backward prediction data matrix
        h = np.conj(Y[0:N - L])					# backward prediction data vector
        U, S, Vh = np.linalg.svd(A, full_matrices=1)                       # singular value decomposition
        V = np.conj(np.transpose(Vh))
        bias = np.mean(S[nFreq:np.min([N - L - 1, L]) + 1])  # bias compensation
        PolyCoef = np.dot(-V[:, 0:nFreq], np.dot(np.diag(1 / (S[0:nFreq] - bias)), np.dot(np.conj(np.transpose(U[:, 0:nFreq])), h)))  # prediction polynomial coefficients
        s = np.conj(np.log(np.roots(np.append(PolyCoef[::-1], 1))))		# polynomial rooting
        s = s[np.where(s < 0)[0]]
        reconstructed = np.zeros(nPredict, dtype=np.complex128)
        if len(s):  # If there are found frequencies
            Z = np.zeros([N, len(s)], dtype=np.complex)
            for k in range(0, len(s)):
                Z[:, k] = np.exp(s[k])**np.arange(0, N)
            a = np.linalg.lstsq(Z, Y)[0]
            para = np.array([-np.real(s), np.imag(s) / 2 / np.pi, np.abs(a), np.imag(np.log(a / np.abs(a)))])  # WF: reintroduce the scaling factor
            if Direction == 1:  # If backward
                xpredict = np.arange(-nPredict, 0)
                for signal in range(para.shape[1]):
                    reconstructed += para[2, signal] * np.exp(1j * (xpredict * para[1, signal] * 2 * np.pi + para[3, signal])) * np.exp(-xpredict * para[0, signal])
                data = np.concatenate((reconstructed, data))
            else:
                xpredict = np.arange(N, N + nPredict)
                for signal in range(para.shape[1]):
                    reconstructed += para[2, signal] * np.exp(1j * (xpredict * para[1, signal] * 2 * np.pi + para[3, signal])) * np.exp(-xpredict * para[0, signal])
                data = np.concatenate((data, reconstructed))
        else:
            data = np.concatenate((reconstructed, data))
        return data

    def setSpec(self, val, axis):
        axis = self.checkAxis(axis)
        oldVal = self.spec[axis]
        self.spec[axis] = val
        self.resetXax(axis)
        if val:
            self.addHistory("Dimension " + str(axis + 1) + " set to FID")
        else:
            self.addHistory("Dimension " + str(axis + 1) + " set to spectrum")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.setSpec(oldVal, axis))

    def swapEcho(self, idx, axis):
        axis = self.checkAxis(axis)
        slicing1 = (slice(None), ) * axis + (slice(None, idx), )
        slicing2 = (slice(None), ) * axis + (slice(idx, None), )
        self.data = self.data[slicing2].append(self.data[slicing1], axis=axis)
        self.wholeEcho[axis] = not self.wholeEcho[axis]
        self.addHistory("Swap echo at position " + str(idx) + " for dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.swapEcho(-idx, axis))

    def shift(self, shift, axis, select=slice(None), zeros=True):
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
        if zeros:
            self.data[select] *= mask.reshape(mask.shape + (1,)*(self.ndim()-axis-1)) 
        if self.spec[axis] > 0:
            self.__fourier(axis, tmp=True)
        Message = "Shifted " + str(shift) + " points in dimension " + str(axis + 1)

        if type(select) is not slice:
            Message = Message + " with slice " + str(select)
        elif select != slice(None, None, None):
            Message = Message + " with slice " + str(select)
        self.addHistory(Message)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.shift(shift, axis, select=select, zeros=zeros)))
        
    def __fourier(self, axis, tmp=False, reorder=[True,True]):
        axis = self.checkAxis(axis)
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

    def __invFourier(self, axis, tmp=False, reorder=[True,True]):
        axis = self.checkAxis(axis)
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

    def complexFourier(self, axis):
        if self.spec[axis] == 0:
            self.__fourier(axis)
            self.addHistory("Fourier transform dimension " + str(axis + 1))
        else:
            self.__invFourier(axis)
            self.addHistory("Inverse fourier transform dimension " + str(axis + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.complexFourier(axis))

    def realFourier(self, axis):
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

    def fftshift(self, axis, inv=False):
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

    def shear(self, shear, axis, axis2):
        axis = self.checkAxis(axis)
        axis2 = self.checkAxis(axis2)
        if axis == axis2:
            raise SpectrumException('Both shearing axes cannot be equal')
        if self.ndim() < 2:
            raise SpectrumException("The data does not have enough dimensions for a shearing transformation")
        shape = self.shape()
        vec1 = np.linspace(0, shear * 2 * np.pi * shape[axis] / self.sw[axis], shape[axis] + 1)[:-1]
        vec2 = np.fft.fftshift(np.fft.fftfreq(shape[axis2], 1 / self.sw[axis2]))
        newShape = [1, ] * self.ndim()
        newShape[axis] = shape[axis]
        newShape[axis2] = shape[axis2]
        if axis > axis2:
            shearMatrix = np.exp(1j * np.outer(vec2, vec1))
        elif axis < axis2:
            shearMatrix = np.exp(1j * np.outer(vec1, vec2))
        if self.spec[axis] > 0: #rorder and fft for spec
            self.__invFourier(axis, tmp=True, reorder = [True,False])
        else: #Reorder if FID
            self.data.icomplexReorder(axis)
        self.data *= shearMatrix.reshape(shape)
        if self.spec[axis] > 0:
            self.__fourier(axis, tmp=True, reorder=[False,True])
        else:
            self.data.icomplexReorder(axis)
        self.addHistory("Shearing transform with shearing value " + str(shear) + " over dimensions " + str(axis + 1) + " and " + str(axis2 + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.shear(-shear, axis, axis2))

    def reorder(self, pos, newLength, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.data.reorder(pos, newLength, axis)
        self.resetXax(axis)
        self.addHistory("Reorder dimension " + str(axis + 1) + " to obtain a new length of " + str(newLength) + " with positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.reorder(pos, newLength, axis)))

    def ffm(self, pos, typeVal, axis):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        # pos contains the values of fixed points which not to be translated to missing points
        posList = np.delete(range(self.shape()[axis]), pos)
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
        fit = pool.map_async(ffm, [(i, posList) for i in tmpData])
        pool.close()
        pool.join()
        tmpData = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axis)
        self.data = hc.HComplexData(tmpData)
        #Transform back to FID
        self.__invFourier(axis, tmp=True)
        self.addHistory("Fast Forward Maximum Entropy reconstruction of dimension " + str(axis + 1) + " at positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, None))

    def clean(self, pos, typeVal, axis, gamma, threshold, maxIter):
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        # pos contains the values of fixed points which not to be translated to missing points
        posList = np.delete(range(self.shape()[axis]), pos)
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
        fit = pool.map_async(clean, [(i, mask, gamma, threshold, maxIter) for i in tmpData])
        pool.close()
        pool.join()
        tmpData = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axis)
        self.data = hc.HComplexData(tmpData)
        #Transform back to FID
        self.__invFourier(axis, tmp=True)
        self.addHistory("CLEAN reconstruction (gamma = " + str(gamma) + " , threshold = " + str(threshold) + " , maxIter = " + str(maxIter) + ") " + 
        "of dimension " + str(axis + 1) + " at positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, None))
        
    def ist(self,pos, typeVal, axis, threshold, maxIter,tracelimit):
        import scipy.signal
        axis = self.checkAxis(axis)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        # pos contains the values of fixed points which not to be translated to missing points
        self.data.icomplexReorder(axis)
        tmpData = self.data.getHyperData(0)
        posList = np.delete(range(tmpData.shape[axis]), pos)
        if typeVal == 1:  # type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList / 2), dtype=int)
        elif typeVal == 2:  # type is TPPI, for now handle the same as Complex
            pass
        NDmax = np.max(np.max(np.abs(np.real(np.fft.fft(tmpData, axis=axis))))) #Get max of ND matrix
        tmpData = np.rollaxis(tmpData, axis, tmpData.ndim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((int(tmpData.size / tmpShape[-1]), tmpShape[-1]))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(ist, [(i, posList, threshold, maxIter, tracelimit, NDmax) for i in tmpData])
        pool.close()
        pool.join()
        tmpData = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axis)
        self.data = hc.HComplexData(tmpData)
        #Transform back to FID
        self.__invFourier(axis, tmp=True)
        self.addHistory("IST reconstruction (threshold = " + str(threshold) + " , maxIter = " + str(maxIter) + " , tracelimit = " + str(tracelimit*100) + ") " + 
        "of dimension " + str(axis + 1) + " at positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, None))

    def getSlice(self, axes, locList, stack=None):
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
                                           [self.spec[axis] for axis in axes],
                                           [self.wholeEcho[axis] for axis in axes],
                                           [self.ref[axis] for axis in axes],
                                           [self.xaxArray[axis][stack[i]] for i, axis in enumerate(axes)],
                                           self.history,
                                           name=self.name))
        sliceSpec.noUndo = True
        return sliceSpec
                
    def restoreData(self, copyData, returnValue):  # restore data from an old copy for undo purposes
        if (not self.noUndo) and returnValue is None:
            copyData2 = copy.deepcopy(self)
        self.data = copyData.data
        self.freq = copyData.freq  # array of center frequency (length is dim, MHz)
        self.filePath = copyData.filePath
        self.sw = copyData.sw  # array of sweepwidths
        self.spec = copyData.spec
        self.wholeEcho = copyData.wholeEcho
        self.xaxArray = copyData.xaxArray
        self.ref = copyData.ref
        self.addHistory("Data was restored to a previous state ")
        self.redoList = []
        if (not self.noUndo) and returnValue is None:
            self.undoList.append(lambda self: self.restoreData(copyData2, None))
        else:
            self.undoList.append(returnValue)
