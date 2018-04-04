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


#########################################################################
# the generic spectrum class


class Spectrum(object):

    def __init__(self, data, filePath, freq, sw, spec=None, wholeEcho=None, ref=None, xaxArray=None, history=None, msgHandler=None, name=''):
        self.name = name
        if isinstance(data, hc.HComplexData):
            self.data = data
        else:
            self.data = hc.HComplexData(data)
        self.filePath = filePath
        self.freq = np.array(freq)  # array of center frequency (length is dim, MHz)
        self.sw = sw  # array of sweepwidths
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
        self.msgHandler = msgHandler

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

    def dispMsg(self, msg):
        if self.msgHandler is None:
            print(msg)
        else:
            self.msgHandler(msg)

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

    def setNoUndo(val):
        self.noUndo = bool(val)
        if self.noUndo:
            self.undoList = []
            self.redoList = []

    def undo(self):
        undoFunc = None
        while undoFunc is None and self.undoList:
            undoFunc = self.undoList.pop()
        if undoFunc is None:
            self.dispMsg("no undo information")
            return False
        tmpRedo = self.redoList # Protect the redo list
        undoFunc(self)
        self.redoList = tmpRedo
        self.redoList.append(self.undoList.pop())
        message = self.removeFromHistory(2)
        self.dispMsg("Undo: " + message)
        return True

    def redo(self):
        if self.redoList:
            tmpRedo = self.redoList # Protect the redo list
            tmpRedo.pop()(self)
            self.redoList = tmpRedo # Restore the redo list
            return True
        else:
            self.dispMsg("no redo information")
            return False

    def clearUndo(self):
        self.undoList = []
        self.redoList = []

    def reload(self):
        loadData = io.autoLoad(*self.filePath)
        self.restoreData(loadData, None)
        
    def checkAxes(self, axes):
        if axes < 0:
            axes = axes + self.ndim()
        if not (0 <= axes < self.ndim()):
            raise IndexError("Not a valid axis for Spectrum")
        return axes

    def resetXax(self, axes=None):
        if axes is not None:
            axes = self.checkAxes(axes)
            val = [axes]
        else:
            val = range(self.ndim())
        for i in val:
            if self.spec[i] == 0:
                self.xaxArray[i] = np.arange(self.shape()[i]) / (self.sw[i])
            elif self.spec[i] == 1:
                self.xaxArray[i] = np.fft.fftshift(np.fft.fftfreq(self.shape()[i], 1.0 / self.sw[i]))
                if self.ref[i] is not None:
                    self.xaxArray[i] += self.freq[i] - self.ref[i]

    def setXax(self, xax, axes):
        axes = self.checkAxes(axes)
        if len(xax) != self.shape()[axes]:
            self.dispMsg("Length of new x-axis does not match length of the data")
            return
        oldXax = self.xaxArray[axes]
        self.xaxArray[axes] = xax
        self.addHistory("X-axes of dimension " + str(axes + 1) + " was set to " + str(xax).replace('\n', ''))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.setXax(oldXax, axes))

    def insert(self, data, pos, axes):
        if not isinstance(data, hc.HComplexData):
            data = hc.HComplexData(data)
        if self.noUndo:
            returnValue = None
        else:
            if self.data.hyper == data.hyper: # If both sets have same hyper: easy undo can be used
                returnValue = lambda self: self.delete(range(pos, pos + data.shape()[axes]), axes)
            else: # Otherwise: do a deep copy of the class
                copyData = copy.deepcopy(self)
                returnValue = lambda self: self.restoreData(copyData, lambda self: self.insert(data, pos, axes))
        axes = self.checkAxes(axes)
        # Check for a change in dimensions
        self.data = self.data.insert(pos, data, axes)
        self.resetXax(axes)
        self.addHistory("Inserted " + str(data.shape()[axes]) + " datapoints in dimension " + str(axes + 1) + " at position " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(returnValue)

    def delete(self, pos, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        tmpData = self.data.delete(pos, axes)
        if 0 in tmpData.shape():
            self.dispMsg('Cannot delete all data')
            return
        self.data = tmpData
        self.xaxArray[axes] = np.delete(self.xaxArray[axes], pos)
        if isinstance(pos, (int, float)):
            length = 1
        else:
            length = len(pos)
        self.addHistory("Removed " + str(length) + " datapoints from dimension " + str(axes + 1) + " at position " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.delete(pos, axes)))

    def add(self, data, axis=None, select=slice(None)):
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

    def subtract(self, data, axis=None, select=slice(None)):
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

    def multiply(self, data, axis=None, select=slice(None)):
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

    def divide(self, data, axis=None, select=slice(None)):
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

    def normalize(self, mult, scale, type, axes, select=slice(None)):
        axes = self.checkAxes(axes)
        try:
            self.data *= mult * scale 
        except ValueError as error:
            self.dispMsg('Normalize: ' + str(error))
            return None
        if type == 0:
            self.addHistory("Normalized integral of dimension " + str(axes + 1) + " of data[" + str(select) + "] to " + str(scale))
        elif type == 1:
            self.addHistory("Normalized maximum of dimension " + str(axes + 1) + " of data[" + str(select) + "] to " + str(scale))
        elif type == 2:
            self.addHistory("Normalized minimum of dimension " + str(axes + 1) + " of data[" + str(select) + "] to " + str(scale))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.normalize(1.0 / mult, scale, type, axes, select=select))

    def baselineCorrection(self, baseline, axes, select=slice(None)):
        axes = self.checkAxes(axes)
        baselinetmp = baseline.reshape((self.shape()[axes], ) + (1, ) * (self.ndim() - axes - 1))
        self.data[select] -= baselinetmp
        self.addHistory("Baseline corrected dimension " + str(axes + 1) + " of data[" + str(select) + "]")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.baselineCorrection(-baseline, axes, select=select))

    def concatenate(self, axes):
        axes = self.checkAxes(axes)
        splitVal = self.shape()[axes]
        copyData = None
        if self.data.isComplex(axes):
            if not self.noUndo:
                copyData = copy.deepcopy(self)
            self.data = self.data.real(axes)
        invAxis = self.ndim() - axes
        self.data = self.data.concatenate(axes)
        self.data.removeDim(invAxis)
        self.freq = np.delete(self.freq, axes)
        self.sw = np.delete(self.sw, axes)
        self.spec = np.delete(self.spec, axes)
        self.wholeEcho = np.delete(self.wholeEcho, axes)
        self.ref = np.delete(self.ref, axes)
        del self.xaxArray[axes]
        self.resetXax()
        self.addHistory("Concatenated dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            if copyData is not None:
                self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.concatenate(axes)))
            else:
                self.undoList.append(lambda self: self.split(splitVal, axes))

    def split(self, sections, axes):
        axes = self.checkAxes(axes)
        invAxis = self.ndim() - axes
        self.data = self.data.split(sections, axes)
        self.data.insertDim(invAxis)
        self.freq = np.insert(self.freq, 0, self.freq[axes])
        self.sw = np.insert(self.sw, 0, self.sw[axes])
        self.spec = np.insert(self.spec, 0, self.spec[axes])
        self.wholeEcho = np.insert(self.wholeEcho, 0, self.wholeEcho[axes])
        self.ref = np.insert(self.ref, 0, self.ref[axes])
        self.xaxArray.insert(0, [])
        self.resetXax(0)
        self.resetXax(axes + 1)
        self.addHistory("Split dimension " + str(axes + 1) + " into " + str(sections) + " sections")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.concatenate(axes))

    def real(self, axes=-1):
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axes = self.checkAxes(axes)
        self.data = self.data.real(axes)
        self.addHistory("Real along dimension " + str(axes+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.real(axes)))

    def imag(self, axes=-1):
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axes = self.checkAxes(axes)
        self.data = self.data.imag(axes)
        self.addHistory("Imaginary along dimension " + str(axes+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.imag(axes)))

    def abs(self, axes=-1):
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axes = self.checkAxes(axes)
        self.data = self.data.abs(axes)
        self.addHistory("Absolute along dimension " + str(axes+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.abs(axes)))
    
    def conj(self, axis=-1):
        self.data = self.data.conj(axis)
        self.addHistory("Complex conjugate along" + str(axis+1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.conj(axis))

    def states(self, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.shape()[axes] % 2 != 0:
            self.dispMsg("States: data has to be even")
            return
        self.data.states(axes)
        self.resetXax(axes)
        self.addHistory("States conversion on dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.states(axes)))

    def statesTPPI(self, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.shape()[axes] % 2 != 0:
            self.dispMsg("States-TPPI: data has to be even")
            return
        self.data.states(axes, TPPI=True)
        self.resetXax(axes)
        self.addHistory("States-TPPI conversion on dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.statesTPPI(axes)))

    def echoAntiEcho(self, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.shape()[axes] % 2 != 0:
            self.dispMsg("Echo-antiecho: data has to be even")
            return None
        self.data.echoAntiEcho(axes)
        self.resetXax(axes)
        self.addHistory("Echo-antiecho conversion on dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.echoAntiEcho(axes)))

    def subtractAvg(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if not (0 <= pos1 <= self.shape()[axes]):
            self.dispMsg("Indices not within range")
            return
        if not (0 <= pos2 <= self.shape()[axes]):
            self.dispMsg("Indices not within range")
            return
        if pos1 == pos2:
            self.dispMsg("Indices cannot be equal")
            return
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axes + (slice(minPos, maxPos), )
        averages = self.data[slicing].mean(axis=axes, keepdims=True)
        self.data -= averages
        self.addHistory("Subtracted average determined between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.add(averages))

    def matrixManip(self, pos1, pos2, axes, which):
        axes = self.checkAxes(axes)
        if isinstance(pos1, int):
            pos1 = np.array([pos1])
            pos2 = np.array([pos2])
        if len(pos1) != len(pos2):
            self.dispMsg("Length of the two arrays is not equal")
            return
        if len(pos1) == 1:
            keepdims = False
        else:
            keepdims = True
        tmpdata = ()
        for i in range(len(pos1)):
            if not (0 <= pos1[i] <= self.shape()[axes]):
                self.dispMsg("Indices not within range")
                return
            if not (0 <= pos2[i] <= self.shape()[axes]):
                self.dispMsg("Indices not within range")
                return
            if pos1[i] == pos2[i]:
                self.dispMsg("Indices cannot be equal")
                return
            minPos = min(pos1[i], pos2[i])
            maxPos = max(pos1[i], pos2[i])
            slicing = (slice(None), ) * axes + (slice(minPos, maxPos), )
            if which == 0:
                if self.spec[axes] == 0:
                    tmp = self.data[slicing].sum(axis=axes) / self.sw[axes]
                else:
                    tmp = self.data[slicing].sum(axis=axes) * self.sw[axes] / (1.0 * self.shape()[axes])
            elif which == 5:
                tmp = self.data[slicing].sum(axis=axes)
            elif which == 1:
                tmp = self.data[slicing].max(axis=axes)
            elif which == 2:
                tmp = self.data[slicing].min(axis=axes)
            elif which == 3:
                tmp = self.data[slicing].argmax(axis=axes)
                maxArgPos = np.array(tmp.data, dtype=int)
                tmpmaxPos = maxArgPos.flatten()
                tmp.data = self.xaxArray[axes][slice(minPos, maxPos)][tmpmaxPos].reshape(maxArgPos.shape)
            elif which == 4:
                tmp = self.data[slicing].argmin(axis=axes)
                minArgPos = np.array(tmp.data, dtype=int)
                tmpminPos = minArgPos.flatten()
                tmp.data = self.xaxArray[axes][slice(minPos, maxPos)][tmpminPos].reshape(minArgPos.shape)
            elif which == 6:
                tmp = self.data[slicing].mean(axis=axes)
            if keepdims:
                tmpdata += (tmp.expand_dims(axes), )
            else:
                tmpdata += (tmp, )
        if not keepdims:
            if self.ndim() == 1:
                self.data = tmpdata[0].reshape((1, ))
                self.resetXax(axes)
            else:
                self.data = tmpdata[0]
                self.freq = np.delete(self.freq, axes)
                self.ref = np.delete(self.ref, axes)
                self.sw = np.delete(self.sw, axes)
                self.spec = np.delete(self.spec, axes)
                self.wholeEcho = np.delete(self.wholeEcho, axes)
                del self.xaxArray[axes]
                self.data.removeDim(axes)
        else:
            self.data = tmpdata[0]
            for extra in tmpdata[1:]:
                self.data = self.data.append(extra, axis=axes)
            self.resetXax(axes)

    def integrate(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axes, which=0)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.integrate(pos1, pos2, axes)))
        self.addHistory("Integrate between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))

    def max(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axes, which=1)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.max(pos1, pos2, axes)))
        self.addHistory("Maximum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))

    def min(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axes, which=2)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.min(pos1, pos2, axes)))
        self.addHistory("Minimum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))

    def argmax(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axes, which=3)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.argmax(pos1, pos2, axes)))
        self.addHistory("Maximum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))

    def argmin(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axes, which=4)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.argmin(pos1, pos2, axes)))
        self.addHistory("Minimum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))

    def sum(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axes, which=5)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.sum(pos1, pos2, axes)))
        self.addHistory("Sum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))

    def average(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.matrixManip(pos1, pos2, axes, which=6)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.average(pos1, pos2, axes)))
        self.addHistory("Average between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))

    def extract(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axes + (slice(minPos, maxPos), )
        if self.spec[axes] == 1:
            oldFxax = self.xaxArray[axes][slice(minPos, maxPos)][0]
            self.sw[axes] = self.sw[axes] * (maxPos - minPos) / (1.0 * self.shape()[axes])
        self.data = self.data[slicing]
        if self.spec[axes] == 1:
            newFxax = np.fft.fftshift(np.fft.fftfreq(self.shape()[axes], 1.0 / self.sw[axes]))[0]
            if self.ref[axes] is None:
                self.ref[axes] = self.freq[axes]
            self.freq[axes] = self.ref[axes] - newFxax + oldFxax
        self.resetXax(axes)
        self.addHistory("Extracted part between " + str(minPos) + " and " + str(maxPos) + " of dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.extract(axes, pos1, pos2)))

    def fiddle(self, refSpec, lb, axes):
        axes = self.checkAxes(axes)
        axLen = self.shape()[axes]
        if len(refSpec) != axLen:
            self.dispMsg("Reference FID does not have the correct length")
            return
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        tmpSpec = np.fft.ifftshift(np.real(refSpec))
        pos = np.argmax(tmpSpec)
        refFid = np.fft.ifft(tmpSpec)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True)
        t = np.arange(axLen) / self.sw[axes]
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
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True)
        self.addHistory("FIDDLE over dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.fiddle(refSpec, lb, axes)))

    def diff(self, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.data.diff(axis=axes)
        self.resetXax(axes)
        self.addHistory("Differences over dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.diff(axes)))

    def cumsum(self, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.data.cumsum(axis=axes)
        self.addHistory("Cumulative sum over dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.cumsum(axes)))

    def flipLR(self, axes):
        axes = self.checkAxes(axes)
        slicing = (slice(None), ) * axes + (slice(None, None, -1), )
        self.data = self.data[slicing]
        self.addHistory("Flipped dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.flipLR(axes))

    def hilbert(self, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.data.complexReorder(axes)
        self.data = self.data.hilbert(axis=axes)
        self.data = self.data.complexReorder(axes)
        self.addHistory("Hilbert transform on dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.hilbert(axes)))

    def autoPhase(self, phaseNum, axes, locList, returnPhases=False):
        axes = self.checkAxes(axes)
        if len(locList) != self.ndim():
            self.dispMsg("Data does not have the correct number of dimensions")
            return
        if np.any(locList >= np.array(self.shape())) or np.any(np.array(locList) < 0):
            self.dispMsg("The location array contains invalid indices")
            return
        locList = np.array(locList, dtype=object)
        locList[axes] = slice(None)
        self.data = self.data.complexReorder(axes)
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True)
        tmp = self.data[locList]
        tmp = tmp.getHyperData(0)
        x = np.fft.fftshift(np.fft.fftfreq(len(tmp), 1.0 / self.sw[axes])) / self.sw[axes]
        # only optimize on the hyper real data
        if phaseNum == 0:
            phases = scipy.optimize.minimize(self.ACMEentropy, [0], (tmp, x, False), method='Powell',options = {'xtol': AUTOPHASETOL})
            phase0 = phases['x']
            phase1 = 0.0
        elif phaseNum == 1:
            phases = scipy.optimize.minimize(self.ACMEentropy, [0, 0], (tmp, x), method='Powell', options = {'xtol': AUTOPHASETOL})
            phase0 = phases['x'][0]
            phase1 = phases['x'][1]
        if self.ref[axes] is None:
            offset = 0
        else:
            offset = self.freq[axes] - self.ref[axes]
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.shape()[axes], 1.0 / self.sw[axes]) + offset) / self.sw[axes] * phase1 * 1j)
        vector = vector.reshape(vector.shape + (1, )*(self.ndim()-axes-1))
        self.data *= np.exp(phase0 * 1j) * vector
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True, inv=True)
        self.data = self.data.complexReorder(axes)
        Message = "Autophase: phase0 = " + str(phase0 * 180 / np.pi) + " and phase1 = " + str(phase1 * 180 / np.pi) + " for dimension " + str(axes + 1)
        self.addHistory(Message)
        if returnPhases:
            if phaseNum == 0:
                return [phases['x']]
            else:
                return phases['x']
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.phase(-phase0, -phase1, axes))

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

    def phase(self, phase0, phase1, axes, select=slice(None)):
        axes = self.checkAxes(axes)
        if self.ref[axes] is None:
            offset = 0
        else:
            offset = self.freq[axes] - self.ref[axes]
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.shape()[axes], 1.0 / self.sw[axes]) + offset) / self.sw[axes] * phase1 * 1j)
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True)
        vector = vector.reshape(vector.shape + (1, )*(self.ndim()-axes-1))
        self.data = self.data.complexReorder(axes)
        self.data[select] *= np.exp(phase0 * 1j) * vector
        self.data = self.data.complexReorder(axes)
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True, inv=True)
        Message = "Phasing: phase0 = " + str(phase0 * 180 / np.pi) + " and phase1 = " + str(phase1 * 180 / np.pi) + " for dimension " + str(axes + 1)
        if select != slice(None, None, None):
            Message = Message + " of data[" + str(select) + "]"
        self.addHistory(Message)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.phase(-phase0, -phase1, axes, select=select))

    def apodize(self, lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, axes, select=slice(None), preview=False):
        axes = self.checkAxes(axes)
        if shiftingAxes is None:
            shiftingAxes = 0
            shifting = 0
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axLen = self.shape()[axes]
        t = np.arange(0, axLen) / self.sw[axes]
        if shifting != 0.0:
            if self.spec[shiftingAxes]:
                shift1 = shift + shifting * np.arange(self.shape()[shiftingAxes]) / self.sw[shiftingAxes]
            else:
                shift1 = shift + shifting * self.xaxArray[shiftingAxes]
            previewData = np.array([func.apodize(t, s, self.sw[axes], axLen, lor, gauss, cos2, hamming, self.wholeEcho[axes]) for s in shift1])
            if axes < shiftingAxes:
                previewData = np.swapaxis(previewData, 0, 1)
            multShape = np.ones(len(self.shape()), dtype=int)
            multShape[axes] = self.shape()[axes]
            multShape[shiftingAxes] = self.shape()[shiftingAxes]
            previewData = previewData.reshape(multShape)
            if self.spec[axes] > 0:
                self.fourier(axes, tmp=True)
            if not isinstance(select, slice):
                previewSelect = np.full(len(previewData.shape), slice(None))
                previewSelect[shiftingAxes] = select[shiftingAxes]
                previewData = previewData[tuple(previewSelect)]
                select = tuple(select)
            self.data[select] *= previewData
            if self.spec[axes] > 0:
                self.fourier(axes, tmp=True, inv=True)
        else:
            x = func.apodize(t, shift, self.sw[axes], axLen, lor, gauss, cos2, hamming, self.wholeEcho[axes])
            if preview:
                previewData = [x] * int(np.prod(self.data.shape()) / self.data.shape()[axes])
            if self.spec[axes] > 0:
                self.fourier(axes, tmp=True)
            self.data[select] *= x.reshape(x.shape + (1, )*(self.ndim()-axes-1))
            if self.spec[axes] > 0:
                self.fourier(axes, tmp=True, inv=True)
        # Create the history message based on the input values.
        Message = "Apodization: "
        if lor is not None:
            Message = Message + "Lorentzian = " + str(lor) + ", "
        if gauss is not None:
            Message = Message + "Gaussian = " + str(gauss) + ", "
        if cos2 is not None:
            Message = Message + "Cos2 = " + str(cos2) + ", "
        if hamming is not None:
            Message = Message + "Hamming = " + str(hamming) + ", "
        if shift != 0.0:
            Message = Message + "shift = " + str(shift) + ", "
        if shifting != 0:
            Message = Message + "shifting = " + str(shifting) + ", "
        if shiftingAxes != 0:
            Message = Message + "shiftingAxes = " + str(shiftingAxes) + ", "
        if lor is None and gauss is None and cos2 is None and hamming is None:  # If all none, make special message with `zero apodization'
            Message = Message + "zero apodization"
        Message = Message + " for dimension " + str(axes + 1)
        if isinstance(select, slice):
            Message = Message + " of data[" + str(select) + "]"
        self.addHistory(Message)
        if preview:
            return ([t]*len(previewData), previewData)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, axes, select=select)))

    def setFreq(self, freq, sw, axes):
        axes = self.checkAxes(axes)
        oldFreq = self.freq[axes]
        oldSw = self.sw[axes]
        if freq is None:
            freq = self.freq[axes]
        if sw is None:
            sw = self.sw[axes]
        self.freq[axes] = float(freq)
        self.sw[axes] = float(sw)
        self.resetXax(axes)
        self.addHistory("Frequency set to " + str(freq * 1e-6) + " MHz and sw set to " + str(sw * 1e-3) + " kHz for dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.setFreq(oldFreq, oldSw, axes))

    def setRef(self, ref, axes):
        axes = self.checkAxes(axes)
        oldRef = self.ref[axes]
        if ref is None:
            self.ref[axes] = None
            self.addHistory("Reference frequency set to 'None' for dimension " + str(axes + 1))
        else:
            self.ref[axes] = float(ref)
            self.addHistory("Reference frequency set to " + str(ref * 1e-6) + " MHz for dimension " + str(axes + 1))
        self.resetXax(axes)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.setRef(oldRef, axes))

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

    def setWholeEcho(self, val, axes):
        axes = self.checkAxes(axes)
        self.wholeEcho[axes] = val
        self.addHistory("Whole echo set to " + str(val) + " for dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.setWholeEcho(not val, axes))

    def resize(self, size, pos, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True)
        self.data = self.data.resize(size, pos, axis=axes)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True)
        self.resetXax(axes)
        self.addHistory("Resized dimension " + str(axes + 1) + " to " + str(size) + " points at position " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.resize(size, pos, axes)))

    def lpsvd(self, nAnalyse, nFreq, nPredict, Direction, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data.apply_along_axis(self.lpsvdfunction, axes, nAnalyse, nFreq, nPredict, Direction)
        self.resetXax(axes)
        self.addHistory("LPSVD ")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.lpsvd(nAnalyse, nFreq, nPredict, Direction, axes)))

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

    def setSpec(self, val, axes):
        axes = self.checkAxes(axes)
        oldVal = self.spec[axes]
        self.spec[axes] = val
        self.resetXax(axes)
        if val:
            self.addHistory("Dimension " + str(axes + 1) + " set to FID")
        else:
            self.addHistory("Dimension " + str(axes + 1) + " set to spectrum")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.setSpec(oldVal, axes))

    def swapEcho(self, idx, axes):
        axes = self.checkAxes(axes)
        slicing1 = (slice(None), ) * axes + (slice(None, idx), )
        slicing2 = (slice(None), ) * axes + (slice(idx, None), )
        self.data = self.data[slicing2].append(self.data[slicing1], axis=axes)
        self.wholeEcho[axes] = not self.wholeEcho[axes]
        self.addHistory("Swap echo at position " + str(idx) + " for dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.swapEcho(-idx, axes))

    def shift(self, shift, axes, select=slice(None), zeros=True):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True)
        mask = np.ones(self.shape()[axes])
        if shift < 0:
            mask[slice(shift, None)] = 0
        else:
            mask[slice(None, shift)] = 0
        self.data[select] = self.data.roll(shift, axes)[select]
        if zeros:
            self.data[select] *= mask.reshape(mask.shape + (1,)*(self.ndim()-axes-1)) 
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True)
        Message = "Shifted " + str(shift) + " points in dimension " + str(axes + 1)
        if select != slice(None, None, None):
            Message = Message + " of data[" + str(select) + "]"
        self.addHistory(Message)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.shift(shift, axes, select=select, zeros=zeros)))
        
    def fourier(self, axes, tmp=False, inv=False, reorder=[True,True]):
        axes = self.checkAxes(axes)
        if reorder[0]:
            self.data = self.data.complexReorder(axes)
        if np.logical_xor(self.spec[axes], inv) == 0:
            if not self.wholeEcho[axes] and not tmp:
                slicing = (slice(None), ) * axes + (0, )
                self.data[slicing] = self.data[slicing] * 0.5
            self.data = self.data.fft(axes).fftshift(axes)
            if not tmp:
                self.spec[axes] = 1
                self.addHistory("Fourier transform dimension " + str(axes + 1))
        else:
            self.data = self.data.ifftshift(axes).ifft(axes)
            if not self.wholeEcho[axes] and not tmp:
                slicing = (slice(None), ) * axes + (0, )
                self.data[slicing] *= 2.0
            if not tmp:
                self.spec[axes] = 0
                self.addHistory("Inverse Fourier transform dimension " + str(axes + 1))
        if reorder[1]:
            self.data = self.data.complexReorder(axes)
        self.resetXax(axes)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.fourier(axes, tmp=False, inv=False, reorder=[True,True]))

    def realFourier(self, axes):
        axes = self.checkAxes(axes)
        self.data = self.data.real(axes)
        self.fourier(axes)
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.realFourier(axes)))

    def fftshift(self, axes, inv=False):
        axes = self.checkAxes(axes)
        if inv:
            self.data = self.data.ifftshift(axis=axes)
            self.addHistory("Inverse Fourier shift dimension " + str(axes + 1))
        else:
            self.data = self.data.fftshift(axis=axes)
            self.addHistory("Fourier shift dimension " + str(axes + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.fftshift(axes, not(inv)))

    def shear(self, shear, axes, axes2):
        axes = self.checkAxes(axes)
        axes2 = self.checkAxes(axes2)
        if axes == axes2:
            self.dispMsg('Both shearing axes cannot be equal')
            return
        if self.ndim() < 2:
            self.dispMsg("The data does not have enough dimensions for a shearing transformation")
            return
        shape = self.shape()
        vec1 = np.linspace(0, shear * 2 * np.pi * shape[axes] / self.sw[axes], shape[axes] + 1)[:-1]
        vec2 = np.fft.fftshift(np.fft.fftfreq(shape[axes2], 1 / self.sw[axes2]))
        newShape = [1, ] * self.ndim()
        newShape[axes] = shape[axes]
        newShape[axes2] = shape[axes2]
        if axes > axes2:
            shearMatrix = np.exp(1j * np.outer(vec2, vec1))
        elif axes < axes2:
            shearMatrix = np.exp(1j * np.outer(vec1, vec2))
        if self.spec[axes] > 0: #rorder and fft for spec
            self.fourier(axes, tmp=True, reorder = [True,False])
        else: #Reorder if FID
            self.data = self.data.complexReorder(axes)
        self.data *= shearMatrix.reshape(shape)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True, reorder=[False,True])
        else:
            self.data = self.data.complexReorder(axes)
        self.addHistory("Shearing transform with shearing value " + str(shear) + " over dimensions " + str(axes + 1) + " and " + str(axes2 + 1))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.shear(-shear, axes, axes2))

    def reorder(self, pos, newLength, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.data.reorder(pos, newLength, axes)
        self.resetXax(axes)
        self.addHistory("Reorder dimension " + str(axes + 1) + " to obtain a new length of " + str(newLength) + " with positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, lambda self: self.reorder(pos, newLength, axes)))

    def ffm(self, pos, typeVal, axes):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        # pos contains the values of fixed points which not to be translated to missing points
        posList = np.delete(range(self.shape()[axes]), pos)
        if typeVal == 1:  # type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList / 2), dtype=int)
        if typeVal == 2:  # type is TPPI, for now handle the same as Complex
            pass
        posList = np.unique(posList)
        self.data = self.complexReorder(axes)
        tmpData = self.data.getHyperData(0)
        tmpData = np.rollaxis(tmpData, axes, tmpData.ndim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((int(tmpData.size / tmpShape[-1]), tmpShape[-1]))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(ffm, [(i, posList) for i in tmpData])
        pool.close()
        pool.join()
        tmpData = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axes)
        self.data = hc.HComplexData(tmpData)
        #Transform back to FID
        self.fourier(axes, tmp=True, inv=True)
        self.addHistory("Fast Forward Maximum Entropy reconstruction of dimension " + str(axes + 1) + " at positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, None))

    def clean(self, pos, typeVal, axes, gamma, threshold, maxIter):
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        # pos contains the values of fixed points which not to be translated to missing points
        posList = np.delete(range(self.shape()[axes]), pos)
        if typeVal == 1:  # type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList / 2), dtype=int)
        if typeVal == 2:  # type is TPPI, for now handle the same as Complex
            pass
        posList = np.unique(posList)
        self.data = self.data.complexReorder(axes)
        tmpData = self.data.getHyperData(0)
        tmpData = np.rollaxis(np.fft.fft(tmpData, axis=axes), axes, tmpData.ndim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((int(tmpData.size / tmpShape[-1]), tmpShape[-1]))
        mask = np.ones(tmpShape[-1]) / float(tmpShape[-1])
        mask[posList] = 0.0
        mask = np.fft.fft(mask) # abs or real???
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(clean, [(i, mask, gamma, threshold, maxIter) for i in tmpData])
        pool.close()
        pool.join()
        tmpData = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axes)
        self.data = hc.HComplexData(tmpData)
        #Transform back to FID
        self.fourier(axes, tmp=True, inv=True)
        self.addHistory("CLEAN reconstruction (gamma = " + str(gamma) + " , threshold = " + str(threshold) + " , maxIter = " + str(maxIter) + ") " + 
        "of dimension " + str(axes + 1) + " at positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, None))
        
    def ist(self,pos, typeVal, axes, threshold, maxIter,tracelimit):
        import scipy.signal
        axes = self.checkAxes(axes)
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        # pos contains the values of fixed points which not to be translated to missing points
        self.data = self.data.complexReorder(axes)
        tmpData = self.data.getHyperData(0)
        posList = np.delete(range(tmpData.shape[axes]), pos)
        if typeVal == 1:  # type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList / 2), dtype=int)
        elif typeVal == 2:  # type is TPPI, for now handle the same as Complex
            pass
        NDmax = np.max(np.max(np.abs(np.real(np.fft.fft(tmpData, axis=axes))))) #Get max of ND matrix
        tmpData = np.rollaxis(tmpData, axes, tmpData.ndim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((int(tmpData.size / tmpShape[-1]), tmpShape[-1]))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(ist, [(i, posList, threshold, maxIter, tracelimit, NDmax) for i in tmpData])
        pool.close()
        pool.join()
        tmpData = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axes)
        self.data = hc.HComplexData(tmpData)
        #Transform back to FID
        self.fourier(axes, tmp=True, inv=True)
        self.addHistory("IST reconstruction (threshold = " + str(threshold) + " , maxIter = " + str(maxIter) + " , tracelimit = " + str(tracelimit*100) + ") " + 
        "of dimension " + str(axes + 1) + " at positions " + str(pos))
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(lambda self: self.restoreData(copyData, None))

    def getSlice(self, axes, locList, stack=None):
        locList = np.array(locList, dtype=object)
        if stack is None:
            stack = [slice(None)]*len(axes)
        for i, axis in enumerate(axes):
            axes[i] = self.checkAxes(axis)
        locList[axes] = stack
        sliceData = self.data[locList]
        sliceData = sliceData.complexReorder(axes[-1])
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
                                           self.msgHandler,
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
