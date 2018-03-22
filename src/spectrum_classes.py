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
from six import string_types
from spectrumFrame import Plot1DFrame
from safeEval import safeEval
from nus import ffm, clean, ist
import multiprocessing
from matplotlib.pyplot import get_cmap
import matplotlib
import reimplement as reim
import functions as func
import hypercomplex as hc

COLORMAPLIST = ['seismic', 'BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'rainbow', 'jet']
COLORCYCLE = list(matplotlib.rcParams['axes.prop_cycle'])
COLORCONVERTER = matplotlib.colors.ColorConverter()
AUTOPHASETOL = 0.0002 #is ~0.01 degrees


#########################################################################
# the generic spectrum class


class Spectrum(object):

    def __init__(self, name, data, filePath, freq, sw, spec=None, wholeEcho=None, ref=None, xaxArray=None, history=None, msgHandler=None):
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
        if isinstance(select, string_types):
            select = safeEval(select)
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
        if isinstance(select, string_types):
            select = safeEval(select)
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
        if isinstance(select, string_types):
            select = safeEval(select)
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
        if isinstance(select, string_types):
            select = safeEval(select)
        self.data[select] /= data
        self.addHistory("Divided by data[" + str(select) + "]")
        self.redoList = []
        if not self.noUndo:
            self.undoList.append(returnValue)

    def normalize(self, mult, scale, type, axes, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
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
        if isinstance(select, string_types):
            select = safeEval(select)
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
        if len(locList) != self.ndim()-1:
            self.dispMsg("Data does not have the correct number of dimensions")
            return
        if np.any(locList >= np.delete(self.shape(), axes)) or np.any(np.array(locList) < 0):
            self.dispMsg("The location array contains invalid indices")
            return
        self.data = self.data.complexReorder(axes)
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True)
        tmp = self.data[tuple(locList[:axes]) + (slice(None), ) + tuple(locList[axes:])]
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
        if isinstance(select, string_types):
            select = safeEval(select)
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
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if shiftingAxes is None:
            shiftingAxes = 0
            shifting = 0
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axLen = self.shape()[axes]
        t = np.arange(0, axLen) / self.sw[axes]
        if shifting != 0.0:
            previewData = []
            for j in range(self.shape()[shiftingAxes]):
                if self.spec[shiftingAxes]:
                    shift1 = shift + shifting * j / self.sw[shiftingAxes]
                else:
                    shift1 = shift + shifting * self.xaxArray[shiftingAxes][j]
                x = func.apodize(t, shift1, self.sw[axes], axLen, lor, gauss, cos2, hamming, self.wholeEcho[axes])
                if preview:
                    previewData.append(x)
                if self.spec[axes] > 0:
                    self.fourier(axes, tmp=True)
                for i in range(self.shape()[axes]):
                    if axes < shiftingAxes:
                        slicing = (slice(None), ) * axes + (i, )
                    else:
                        slicing = (slice(None), ) * shiftingAxes + (j, )
                    self.data[slicing] *= x[i]
                if self.spec[axes] > 0:
                    self.fourier(axes, tmp=True, inv=True)
        else:
            x = func.apodize(t, shift, self.sw[axes], axLen, lor, gauss, cos2, hamming, self.wholeEcho[axes])
            if preview:
                previewData = [x] * (np.prod(self.data.shape()) / self.data.shape()[axes])
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
        if select != slice(None, None, None):
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
        if isinstance(select, string_types):
            select = safeEval(select)
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

    def getSlice(self, axes, locList):
        axes = self.checkAxes(axes)
        sliceData = self.data[tuple(locList[:axes]) + (slice(None), ) + tuple(locList[axes:])]
        sliceData = sliceData.complexReorder(axes)
        sliceSpec = copy.deepcopy(Spectrum(self.name,
                                           hc.HComplexData(sliceData.getHyperData(0)),
                                           self.filePath,
                                           [self.freq[axes]],
                                           [self.sw[axes]],
                                           [self.spec[axes]],
                                           [self.wholeEcho[axes]],
                                           [self.ref[axes]],
                                           [self.xaxArray[axes]],
                                           self.history,
                                           self.msgHandler))
        sliceSpec.noUndo = True
        return sliceSpec

    def getBlock(self, axes, axes2, locList, stackBegin=None, stackEnd=None, stackStep=None):
        axes = self.checkAxes(axes)
        axes2 = self.checkAxes(axes2)
        stackSlice = reim.floatSlice(stackBegin, stackEnd, stackStep)
        if axes == axes2:
            self.dispMsg("First and second axes are the same")
            return None
        elif axes < axes2:
            sliceData = self.data[tuple(locList[:axes]) + (slice(None), ) + tuple(locList[axes:axes2 - 1]) + (stackSlice, ) + tuple(locList[axes2 - 1:])]
            sliceData = sliceData.complexReorder(axes)
            newData = hc.HComplexData(np.transpose(sliceData.getHyperData(0)))
        elif axes > axes2:
            sliceData = self.data[tuple(locList[:axes2]) + (stackSlice, ) + tuple(locList[axes2:axes - 1]) + (slice(None), ) + tuple(locList[axes - 1:])]
            sliceData = sliceData.complexReorder(axes)
            newData = hc.HComplexData(sliceData.getHyperData(0))            
        sliceSpec = copy.deepcopy(Spectrum(self.name,
                                           newData,
                                           self.filePath,
                                           [self.freq[axes2], self.freq[axes]],
                                           [self.sw[axes2], self.sw[axes]],
                                           [self.spec[axes2], self.spec[axes]],
                                           [self.wholeEcho[axes2], self.wholeEcho[axes]],
                                           [self.ref[axes2], self.ref[axes]],
                                           [self.xaxArray[axes2][stackSlice], self.xaxArray[axes]],
                                           self.history,
                                           self.msgHandler))
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

##################################################################################################
# the class from which the 1d data is displayed, the operations which only edit the content of this class are for previewing


class Current1D(Plot1DFrame):

    X_RESIZE = False
    Y_RESIZE = False

    MARKER = ''
    LINESTYLE = '-'

    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        super(Current1D, self).__init__(root, fig, canvas)
        self.data = data  # the actual spectrum instance
        self.data1D = None  # the data1D
        if duplicateCurrent is None:
            self.axes = len(self.data.shape()) - 1
            self.axes2 = 0
            self.resetLocList()
            self.viewSettings = {"plotType": 0,
                                 "axType": self.root.father.defaultUnits,
                                 "axType2": self.root.father.defaultUnits,
                                 "ppm": self.root.father.defaultPPM,             # display frequency as ppm
                                 "ppm2": self.root.father.defaultPPM,            # display frequency as ppm
                                 "color": self.root.father.defaultColor,
                                 "linewidth": self.root.father.defaultLinewidth,
                                 "grids": self.root.father.defaultGrids,
                                 "colorMap": self.root.father.defaultColorMap,
                                 "contourConst": self.root.father.defaultContourConst,
                                 "contourColors": [self.root.father.defaultPosColor, self.root.father.defaultNegColor],
                                 "diagonalBool": self.root.father.defaultDiagonalBool,
                                 "diagonalMult": self.root.father.defaultDiagonalMult,
                                 "contourSign": 0,
                                 "contourType": 0,
                                 "numLevels": 20,
                                 "minLevels": 0.1,
                                 "maxLevels": 1,
                                 "limitType": 0,
                                 "multiValue": 1.5,
                                 "projTop": 0,
                                 "projRight": 0,
                                 "projLimitsBool": False,
                                 "projLimits": [None, None, None, None],
                                 "projPos": [0, 0],
                                 # CurrentMulti variables
                                 "extraData": [],
                                 "extraLoc": [],
                                 "extraColor": [],
                                 "extraName": [],
                                 "extraAxes": [],
                                 "extraScale": [],
                                 "extraOffset": [],
                                 "extraShift": [],
                                 # CurrentStacked variables
                                 "stackBegin": None,
                                 "stackEnd": None,
                                 "stackStep": None,
                                 "spacing": 0,
            }
            self.upd()  # get the first slice of data
            self.startUp()
        else:
            self.axes = duplicateCurrent.axes
            self.axes2 = duplicateCurrent.axes2
            if isinstance(self, (CurrentStacked, CurrentArrayed, CurrentContour)):
                if self.axes == self.axes2: #If axes are equal, change axes2 to one value below/above
                    if self.axes2 > 0:
                        self.axes2 += -1
                    else:
                        self.axes2 += 1
                if (len(duplicateCurrent.locList) == self.data.ndim() - 2):
                    self.locList = duplicateCurrent.locList
                else:
                    if self.axes < self.axes2:
                        self.locList = np.delete(duplicateCurrent.locList, self.axes2 - 1)
                    else:
                        self.locList = np.delete(duplicateCurrent.locList, self.axes2)
            else:
                if (len(duplicateCurrent.locList) == self.data.ndim() - 1):
                    self.locList = duplicateCurrent.locList
                else:
                    if self.axes < duplicateCurrent.axes2:
                        self.locList = np.insert(duplicateCurrent.locList, duplicateCurrent.axes2 - 1, 0)
                    else:
                        self.locList = np.insert(duplicateCurrent.locList, duplicateCurrent.axes2, 0)
            self.viewSettings = duplicateCurrent.viewSettings
            self.viewSettings.update({"extraData": [],
                                      "extraLoc": [],
                                      "extraColor": [],
                                      "extraName": [],
                                      "extraAxes": [],
                                      "extraScale": [],
                                      "extraOffset": [],
                                      "extraShift": []})
            self.xminlim = duplicateCurrent.xminlim
            self.xmaxlim = duplicateCurrent.xmaxlim
            self.yminlim = duplicateCurrent.yminlim
            self.ymaxlim = duplicateCurrent.ymaxlim
            xReset = self.X_RESIZE or duplicateCurrent.X_RESIZE
            yReset = self.Y_RESIZE or duplicateCurrent.Y_RESIZE
            self.upd()  # get the first slice of data
            self.startUp(xReset, yReset)

    def shape(self):
        return self.data1D.shape()

    def len(self, dim=-1):
        return self.shape()[dim]
    
    def ndim(self):
        return self.data1D.ndim()

    def freq(self, dim=-1):
        return self.data1D.freq[dim]

    def ref(self, dim=-1):
        if self.data1D.ref[dim] is not None:
            return self.data1D.ref[dim]
        else:
            return self.data1D.freq[dim]

    def sw(self, dim=-1):
        return self.data1D.sw[dim]

    def spec(self, dim=-1):
        return self.data1D.spec[dim]

    def xax(self, dim=-1):
        return self.data1D.xaxArray[dim]

    def wholeEcho(self, dim=-1):
        return self.data1D.wholeEcho[dim]

    def getCurrentAxMult(self, axis=-1):
        return self.getAxMult(self.spec(axis), self.viewSettings["axType"], self.viewSettings["ppm"], self.freq(axis), self.ref(axis))
    
    def getDataType(self, data):
        typeList = [np.real, np.imag, np.array, np.abs]
        return typeList[self.viewSettings["plotType"]](data)
    
    def startUp(self, xReset=True, yReset=True):
        self.showFid()  # plot the data
        self.plotReset(xReset, yReset)  # reset the axes limits

    def rename(self, name):
        self.data.rename(name)
        self.canvas.draw()

    def copyCurrent(self, root, fig, canvas, data):
        return Current1D(root, fig, canvas, data, self)

    def upd(self):  # get new data from the data instance
        if self.data.ndim() <= self.axes:
            self.axes = len(self.data.shape()) - 1
        if (len(self.locList) + 1) != self.data.ndim():
            self.resetLocList()
        try:
            self.data1D = self.data.getSlice(self.axes, self.locList)
            if self.data1D is None:
                self.root.rescue()
        except Exception:
            self.resetLocList()
            self.data1D = self.data.getSlice(self.axes, self.locList)
        return True

    def setSlice(self, axes, locList):  # change the slice
        axesSame = True
        if self.axes != axes:
            axesSame = False
            self.axes = axes
        self.locList = locList
        self.upd()
        self.showFid()
        if not axesSame:
            self.plotReset()

    def resetLocList(self):
        self.locList = [0] * (len(self.data.shape()) - 1)

    def getSelect(self):
        tmp = list(self.locList)
        if self.ndim() > 1:
            minVal = min(self.axes, self.axes2)
            maxVal = max(self.axes, self.axes2)
            tmp.insert(minVal, slice(None))
            tmp.insert(maxVal, slice(None))
        else:
            tmp.insert(self.axes, slice(None))
        return tmp

    def setGrids(self, grids):
        self.viewSettings["grids"] = grids

    def setDiagonal(self, diagonalBool=None, diagonalMult=None):
        if diagonalBool is not None:
            self.viewSettings["diagonalBool"] = diagonalBool
        if diagonalMult is not None:
            self.viewSettings["diagonalMult"] = diagonalMult
        self.showFid()

    def isComplex(self, *args):
        return self.data.isComplex(*args)

    def real(self, *args):
        self.root.addMacro(['real', (self.axes - self.data.ndim(), )])
        self.data.real(self.axes)
        self.upd()
        self.showFid()
    
    def imag(self, *args):
        self.root.addMacro(['imag', (self.axes - self.data.ndim(), )])
        self.data.imag(self.axes)
        self.upd()
        self.showFid()
    
    def abs(self, *args):
        self.root.addMacro(['abs', (self.axes - self.data.ndim(), )])
        self.data.abs(self.axes)
        self.upd()
        self.showFid()

    def conj(self, *args):
        self.root.addMacro(['conj', (self.axes - self.data.ndim(), )])
        self.data.conj(self.axes)
        self.upd()
        self.showFid()
    
    def setPhaseInter(self, phase0in, phase1in):  # interactive changing the phase without editing the actual data
        phase0 = float(phase0in)
        phase1 = float(phase1in)
        self.data1D.phase(phase0, phase1, -1)
        self.showFid()
        self.upd()

    def applyPhase(self, phase0, phase1, select=False):  # apply the phase to the actual data
        phase0 = float(phase0)
        phase1 = float(phase1)
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['phase', (phase0, phase1, self.axes - self.data.ndim(), selectSlice)])
        self.data.phase(phase0, phase1, self.axes, selectSlice)
        self.upd()
        self.showFid()

    def fourier(self):  # fourier the actual data and replot
        self.root.addMacro(['fourier', (self.axes - self.data.ndim(), )])
        self.data.fourier(self.axes)
        self.upd()
        if isinstance(self, (CurrentStacked, CurrentArrayed)):
            self.resetSpacing()
        self.showFid()
        self.plotReset()

    def realFourier(self):  # fourier the real data and replot
        self.root.addMacro(['realFourier', (self.axes - self.data.ndim(), )])
        self.data.realFourier(self.axes)
        self.upd()
        if isinstance(self, (CurrentStacked, CurrentArrayed)):
            self.resetSpacing()
        self.showFid()
        self.plotReset()

    def fftshift(self, inv=False):  # fftshift the actual data and replot
        self.root.addMacro(['fftshift', (self.axes - self.data.ndim(), inv)])
        self.data.fftshift(self.axes, inv)
        self.upd()
        self.showFid()

    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=None):  # display the 1D data including the apodization function
        if not type(self) is CurrentContour:
            y = copy.deepcopy(self.data1D.data)
            preview = True
        else:
            preview = False
        if shiftingAxes is None:
            curve = self.data1D.apodize(lor, gauss, cos2, hamming, shift, shifting, None, -1, preview=preview)
        else:
            if (self.ndim() > 1) and (shiftingAxes == self.axes2):
                curve = self.data1D.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, -1, preview=preview)
            else:
                if shiftingAxes == self.axes:
                    self.dispMsg('shiftingAxes cannot be equal to axes')
                    return
                elif shiftingAxes < self.axes:
                    shift += shifting * self.locList[shiftingAxes] / self.data.sw[shiftingAxes]
                else:
                    shift += shifting * self.locList[shiftingAxes - 1] / self.data.sw[shiftingAxes]
                curve = self.data1D.apodize(lor, gauss, cos2, hamming, shift, 0.0, None, -1, preview=preview)
        if not type(self) is CurrentContour:
            if self.spec() == 0:
                tmp = self.getDataType(y.getHyperData(0))
                scale = np.max([np.real(tmp), np.imag(tmp)])
                self.showFid(y, curve[0], scale*np.array(curve[1]), 'g')
            else:
                self.showFid(y)
        else:
            self.showFid()
        self.upd()
    
    def applyApod(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=0, select=False):  # apply the apodization to the actual data
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['apodize', (lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, self.axes - self.data.ndim(), selectSlice)])
        self.data.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, self.axes, selectSlice)
        self.upd()
        self.showFid()

    def setFreq(self, freq, sw):  # set the frequency of the actual data
        self.root.addMacro(['setFreq', (freq, sw, self.axes - self.data.ndim())])
        self.data.setFreq(freq, sw, self.axes)
        self.upd()
        self.showFid()

    def setRef(self, ref):  # set the frequency of the actual data
        oldref = self.ref()
        self.root.addMacro(['setRef', (ref, self.axes - self.data.ndim())])
        self.data.setRef(ref, self.axes)
        if ref is None:
            ref = self.freq()
        val = self.viewSettings["axType"]
        if self.spec() == 1:
            if self.viewSettings["ppm"]:
                self.xminlim = (self.xminlim * oldref * 10**-6 + oldref - ref) / (ref * 10**-6)
                self.xmaxlim = (self.xmaxlim * oldref * 10**-6 + oldref - ref) / (ref * 10**-6)
            else:
                self.xminlim = self.xminlim + (oldref - ref) / 10**(val * 3)  # set new limits, and scale for axis type
                self.xmaxlim = self.xmaxlim + (oldref - ref) / 10**(val * 3)
        self.upd()
        self.showFid()

    def regrid(self, limits, numPoints):
        self.root.addMacro(['regrid', (limits, numPoints, self.axes - self.data.ndim())])
        self.data.regrid(limits, numPoints, self.axes)
        self.upd()
        self.showFid()
        self.plotReset()

    def SN(self, minNoise, maxNoise, minPeak, maxPeak):
        minN = min(minNoise, maxNoise)
        maxN = max(minNoise, maxNoise)
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpData = np.real(self.getDataType(tmpData))
        return (np.max(tmpData[minP:maxP]) / (np.std(tmpData[minN:maxN])))

    def fwhm(self, minPeak, maxPeak, unitType=None):
        from scipy.interpolate import UnivariateSpline
        if unitType is None:
            axType = self.viewSettings["axType"]
            ppm = self.viewSettings["ppm"]
        else:
            axType = unitType
            if unitType == 3:  # ppm
                ppm = 1
            else:
                ppm = 0
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpData = np.real(self.getDataType(tmpData))
        maxPos = np.argmax(tmpData[minP:maxP])
        x = self.xax() * self.getAxMult(self.spec(), axType, ppm, self.freq(), self.ref())
        maxX = x[minP:maxP][maxPos]
        spline = UnivariateSpline(x, tmpData - tmpData[minP:maxP][maxPos] / 2.0, s=0)
        zeroPos = spline.roots()
        left = zeroPos[zeroPos > maxX]
        right = zeroPos[zeroPos < maxX]
        if right.size > 0 and left.size > 0:
            return abs(left[0] - right[-1])
        else:
            return 0.0

    def COM(self, minPeak, maxPeak, unitType=None):  # Centre of Mass
        if unitType is None:
            axType = self.axType
            ppm = self.viewSettings["ppm"]
        else:
            axType = unitType
            if unitType == 3:  # ppm
                ppm = 1
            else:
                ppm = 0
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpAxis = self.xax()
        tmpData = np.real(self.getDataType(tmpData))
        tmpAxis = tmpAxis[minP:maxP] * self.getAxMult(self.spec(), axType, ppm, self.freq(), self.ref())
        tmpData = tmpData[minP:maxP]
        # COM = 1/M *sum(m_i * r_i)
        CentreOM = 1.0 / np.sum(tmpData) * np.sum(tmpData * tmpAxis)
        return CentreOM

    def Integrals(self, minPeak, maxPeak):
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpAxis = self.xax()
        totLen = len(tmpData)
        tmpData = np.real(self.getDataType(tmpData))
        maxim = np.max(np.abs(tmpData))
        tmpAxis = tmpAxis[minP:maxP] 
        tmpData = tmpData[minP:maxP]
        if self.spec() == 0:
            intSum = np.cumsum(tmpData)
            inte = np.sum(tmpData) / self.sw()
        else:
            intSum = np.cumsum(tmpData[-1::-1])[-1::-1]
            inte = np.sum(tmpData) * self.sw() / (1.0 * totLen)
        return inte, tmpAxis, intSum, maxim

    def MaxMin(self, minPeak, maxPeak, type='max'):
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpData = np.real(self.getDataType(tmpData))
        if type == 'max':
            return np.max(tmpData[minP:maxP])
        elif type == 'min':
            return np.min(tmpData[minP:maxP])

    def integralsPreview(self, x, y, maxim):
        xNew = []
        yNew = []
        scale = 0
        for num in range(len(x)):
            if x[num] is not None and y[num] is not None:
                xNew.append(x[num])
                yNew.append(y[num])
                scale = np.max([scale,abs(yNew[-1][0]),abs(yNew[-1][-1])])
        for num in range(len(yNew)):
            yNew[num] = yNew[num] / scale * maxim
        self.showFid(extraX=xNew, extraY=yNew, extraColor='g')

    def resizePreview(self, size, pos):  # set size only on local data
        self.data1D.resize(size, pos, -1)
        self.showFid()
        if not self.spec():
            self.plotReset(True, False)
        self.upd()

    def resize(self, size, pos):  # set size to the actual data
        self.root.addMacro(['resize', (size, pos, self.axes - self.data.ndim())])
        if self.data.noUndo:
            self.data.resize(size, pos, self.axes)
        else:
            self.data.resize(size, pos, self.axes)
        self.upd()
        self.showFid()
        if not self.spec():
            self.plotReset(True, False)

    def lpsvd(self, nAnalyse, nFreq, nPredict, Direction):
        self.root.addMacro(['lpsvd', (nAnalyse, nFreq, nPredict, Direction, self.axes - self.data.ndim())])
        self.data.lpsvd(nAnalyse, nFreq, nPredict, Direction, self.axes)
        self.upd()
        self.showFid()

    def setSpec(self, val):  # change from time to freq domain of the actual data
        self.root.addMacro(['setSpec', (val, self.axes - self.data.ndim())])
        self.data.setSpec(val, self.axes)
        self.upd()
        if isinstance(self, CurrentArrayed):
            self.resetSpacing()
        self.showFid()
        self.plotReset()

    def swapEcho(self, idx):
        self.root.addMacro(['swapEcho', (idx, self.axes - self.data.ndim())])
        self.data.swapEcho(idx, self.axes)
        self.upd()
        self.showFid()

    def swapEchoPreview(self, idx):
        self.data1D.swapEcho(idx, -1)
        self.showFid()
        self.upd()

    def setWholeEcho(self, value):
        valBool = value != 0
        self.root.addMacro(['setWholeEcho', (valBool, self.axes - self.data.ndim())])
        self.data.setWholeEcho(valBool, self.axes)
        self.data1D.setWholeEcho(valBool, self.axes)

    def shift(self, shift, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['shift', (shift, self.axes - self.data.ndim(), selectSlice)])
        self.data.shift(shift, self.axes, selectSlice)
        self.upd()
        self.showFid()

    def shiftPreview(self, shift):
        self.data1D.shift(shift, -1)
        self.showFid()
        self.upd()

    def getdcOffset(self, pos1, pos2):
        minPos = int(min(pos1, pos2))
        maxPos = int(max(pos1, pos2))
        if minPos != maxPos:
            tmpData = self.data1D.data[(len(self.shape()) - 1) * (slice(None), ) + (slice(minPos, maxPos), )]
            return np.mean(tmpData.getHyperData(0))
        else:
            return 0

    def dcOffset(self, offset):
        self.data1D.subtract([offset])
        self.showFid()
        self.upd()

    def baselinePolyFit(self, x, data, bArray, degree):
        import numpy.polynomial.polynomial as poly
        polyCoeff = poly.polyfit(x[bArray], data[bArray], degree)
        return poly.polyval(x, polyCoeff)

    def baselineCorrectionAll(self, degree, removeList, select=False, invert=False):
        tmpAx = np.arange(self.len())
        bArray = np.array([True] * self.len())
        for i in range(int(np.floor(len(removeList) / 2.0))):
            minVal = min(removeList[2 * i], removeList[2 * i + 1])
            maxVal = max(removeList[2 * i], removeList[2 * i + 1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        if invert:
            bArray = np.logical_not(bArray)
        y = np.apply_along_axis(lambda data: self.baselinePolyFit(self.xax(), data, bArray, degree), self.axes, self.data.getHyperData(0))
        y = np.real(self.getDataType(y))
        self.root.addMacro(['subtract', ([y])])
        self.data.subtract([y])

    def baselineCorrection(self, degree, removeList, select=False, invert=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpAx = np.arange(self.len())
        bArray = np.array([True] * self.len())
        for i in range(int(np.floor(len(removeList) / 2.0))):
            minVal = min(removeList[2 * i], removeList[2 * i + 1])
            maxVal = max(removeList[2 * i], removeList[2 * i + 1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        if invert:
            bArray = np.logical_not(bArray)
        y = self.baselinePolyFit(self.xax(), tmpData, bArray, degree)
        y = np.real(self.getDataType(y))
        self.root.addMacro(['baselineCorrection', (y, self.axes - self.data.ndim(), selectSlice)])
        self.data.baselineCorrection(y, self.axes, select=selectSlice)

    def previewBaselineCorrection(self, degree, removeList, invert=False):
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpAx = np.arange(self.len())
        bArray = np.array([True] * self.len())
        for i in range(int(np.floor(len(removeList) / 2.0))):
            minVal = min(removeList[2 * i], removeList[2 * i + 1])
            maxVal = max(removeList[2 * i], removeList[2 * i + 1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        if invert:
            bArray = np.logical_not(bArray)
        check = True
        y = self.baselinePolyFit(self.xax(), tmpData, bArray, degree)
        y = np.real(self.getDataType(y))
        self.data1D.baselineCorrection(y, -1)
        self.resetPreviewRemoveList()
        if check:
            if self.ndim() > 1:
                self.showFid(extraX=[self.xax()], extraY=[y]*self.len(-2), extraColor='g')
            else:
                self.showFid(extraX=[self.xax()], extraY=[y], extraColor='g')
        else:
            self.showFid()
        self.previewRemoveList(removeList, invert)
        self.upd()

    def previewRemoveList(self, removeList, invert=False):
        axMult = self.getAxMult(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"], self.freq(), self.ref())
        self.resetPreviewRemoveList()
        lineColor = 'r'
        if invert:
            lineColor = 'w'
            self.removeListLines.append(self.ax.axvspan(self.xax()[0] * axMult, self.xax()[-1] * axMult, color='r'))
        for i in range(int(np.floor(len(removeList) / 2.0))):
            self.removeListLines.append(self.ax.axvspan(self.xax()[removeList[2 * i]] * axMult, self.xax()[removeList[2 * i + 1]] * axMult, color=lineColor))
        if len(removeList) % 2:
            self.removeListLines.append(self.ax.axvline(self.xax()[removeList[-1]] * axMult, c=lineColor, linestyle='--'))
        self.canvas.draw()

    def resetPreviewRemoveList(self):
        if hasattr(self, 'removeListLines'):
            for i in self.removeListLines:
                i.remove()
        self.removeListLines = []

    def states(self):
        self.root.addMacro(['states', (self.axes - self.data.ndim(), )])
        self.data.states(self.axes)
        self.upd()
        self.showFid()

    def statesTPPI(self):
        self.root.addMacro(['statesTPPI', (self.axes - self.data.ndim(), )])
        self.data.statesTPPI(self.axes)
        self.upd()
        self.showFid()

    def echoAntiEcho(self):
        self.root.addMacro(['echoAntiEcho', (self.axes - self.data.ndim(), )])
        self.data.echoAntiEcho(self.axes)
        self.upd()
        self.showFid()

    def matrixFuncs(self, func, name, pos1, pos2, newSpec=False):
        if newSpec:
            tmpData = copy.deepcopy(self.data)
            func(tmpData, pos1, pos2, self.axes)
            return tmpData
        else:
            self.root.addMacro([name, (pos1, pos2, self.axes - self.data.ndim(), )])
            func(self.data, pos1, pos2, self.axes)
            if self.upd():
                self.showFid()
                self.plotReset()
        
    def integrate(self, pos1, pos2, newSpec=False):
        self.matrixFuncs(lambda obj, *args: obj.integrate(*args), 'integrate', pos1, pos2, newSpec)

    def sum(self, pos1, pos2, newSpec=False):
        self.matrixFuncs(lambda obj, *args: obj.sum(*args), 'sum', pos1, pos2, newSpec)
        
    def max(self, pos1, pos2, newSpec=False):
        self.matrixFuncs(lambda obj, *args: obj.max(*args), 'max', pos1, pos2, newSpec)

    def min(self, pos1, pos2, newSpec=False):
        self.matrixFuncs(lambda obj, *args: obj.min(*args), 'min', pos1, pos2, newSpec)

    def argmax(self, pos1, pos2, newSpec=False):
        self.matrixFuncs(lambda obj, *args: obj.argmax(*args), 'argmax', pos1, pos2, newSpec)

    def argmin(self, pos1, pos2, newSpec=False):
        self.matrixFuncs(lambda obj, *args: obj.argmin(*args), 'argmin', pos1, pos2, newSpec)

    def average(self, pos1, pos2, newSpec=False):
        self.matrixFuncs(lambda obj, *args: obj.average(*args), 'average', pos1, pos2, newSpec)

    def flipLR(self):
        self.data.flipLR(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['flipLR', (self.axes - self.data.ndim(), )])

    def concatenate(self, axes):
        self.data.concatenate(axes)
        self.upd()
        self.showFid()
        self.plotReset()
        self.root.addMacro(['concatenate', (axes - self.data.ndim() - 1, )])

    def split(self, sections):
        self.data.split(sections, self.axes)
        self.upd()
        self.showFid()
        self.plotReset()
        self.root.addMacro(['split', (sections, self.axes - self.data.ndim() + 1)])

    def diff(self):
        self.data.diff(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['diff', (self.axes - self.data.ndim(), )])

    def cumsum(self):
        self.data.cumsum(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['cumsum', (self.axes - self.data.ndim(), )])

    def insert(self, data, pos):
        self.root.addMacro(['insert', (data, pos, self.axes - self.data.ndim())])
        self.data.insert(data, pos, self.axes)
        self.upd()
        self.showFid()
        self.plotReset()

    def delete(self, pos):
        self.data.delete(pos, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['delete', (pos, self.axes - self.data.ndim())])

    def deletePreview(self, pos):
        self.data1D.delete(pos, -1)
        self.showFid()
        self.upd()

    def add(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['add', (data, self.axes - self.data.ndim(), selectSlice)])
        self.data.add(data, self.axes, select=selectSlice)
        self.upd()
        self.showFid()

    def subtract(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['subtract', (data, self.axes - self.data.ndim(), selectSlice)])
        self.data.subtract(data, self.axes, select=selectSlice)
        self.upd()
        self.showFid()

    def multiply(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['multiply', (data, self.axes - self.data.ndim(), selectSlice)])
        self.data.multiply(data, self.axes, select=selectSlice)
        self.upd()
        self.showFid()

    def multiplyPreview(self, data):
        self.data1D.multiply(data, -1)
        self.showFid()
        self.upd()

    def normalize(self, value, scale, type, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['normalize', (value, scale, type, self.axes - self.data.ndim(), selectSlice)])
        self.data.normalize(value, scale, type, self.axes, select=selectSlice)
        self.upd()
        self.showFid()
    
    def subtractAvg(self, pos1, pos2):
        self.root.addMacro(['subtractAvg', (pos1, pos2, self.axes - self.data.ndim())])
        self.data.subtractAvg(pos1, pos2, self.axes)
        self.upd()
        self.showFid()

    def subtractAvgPreview(self, pos1, pos2):
        self.data1D.subtractAvg(pos1, pos2, -1)
        self.showFid()
        self.upd()

    def extract(self, pos1, pos2, newSpec=False):
        if newSpec:
            tmpData = copy.deepcopy(self.data)
            tmpData.extract(pos1, pos2, self.axes)
            return tmpData
        else:
            self.root.addMacro(['extract', (pos1, pos2, self.axes - self.data.ndim())])
            self.data.extract(pos1, pos2, self.axes)
            self.upd()
            self.showFid()
            self.plotReset()

    def fiddle(self, pos1, pos2, lb):
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        refSpec = np.zeros(self.len())
        refSpec[minPos:maxPos] = np.real(tmpData[minPos:maxPos])
        self.root.addMacro(['fiddle', (refSpec.tolist(), lb, self.axes - self.data.ndim())])
        self.data.fiddle(refSpec, lb, self.axes)
        self.upd()
        self.showFid()

    def shearing(self, shear, axes, axes2):
        self.root.addMacro(['shear', (shear, axes - self.data.ndim(), axes2 - self.data.ndim())])
        self.data.shear(shear, axes, axes2)
        self.upd()
        self.showFid()

    def reorder(self, pos, newLength):
        self.root.addMacro(['reorder', (pos, newLength, self.axes - self.data.ndim())])
        self.data.reorder(pos, newLength, self.axes)
        self.upd()
        self.showFid()

    def ffm(self, posList, typeVal):
        self.root.addMacro(['ffm', (posList, typeVal, self.axes - self.data.ndim())])
        self.data.ffm(posList, typeVal, self.axes)
        self.upd()
        self.showFid()

    def clean(self, posList, typeVal, gamma, threshold, maxIter):
        self.root.addMacro(['clean', (posList, typeVal, self.axes - self.data.ndim(), gamma, threshold, maxIter)])
        self.data.clean(posList, typeVal, self.axes, gamma, threshold, maxIter)
        self.upd()
        self.showFid()

    def ist(self, posList, typeVal, threshold, maxIter, tracelimit):
        self.root.addMacro(['ist', (posList, typeVal, self.axes - self.data.ndim(), threshold, maxIter, tracelimit)])
        self.data.ist(posList, typeVal, self.axes, threshold, maxIter, tracelimit)
        self.upd()
        self.showFid()

    def autoPhase(self, phaseNum):
        phases = self.data1D.autoPhase(phaseNum, -1, [0]*(self.ndim()-1), returnPhases=True)
        self.upd()
        return phases

    def directAutoPhase(self, phaseNum):
        tmpLocList = self.locList
        if self.ndim() > 1:
            if hasattr(self, 'stackBegin'):
                val = self.viewSettings["stackBegin"]
            else:
                val = 0
            if self.axes > self.axes2:
                tmpLocList = np.insert(tmpLocList, self.axes2, val)
            else:
                tmpLocList = np.insert(tmpLocList, self.axes2 - 1, val)
        self.root.addMacro(['autoPhase', (phaseNum, self.axes - self.data.ndim(), tmpLocList)])
        self.data.autoPhase(phaseNum, self.axes, tmpLocList)
        self.upd()
        self.showFid()

    def setXaxPreview(self, xax):
        self.data1D.setXax(xax, -1)
        self.showFid()
        self.plotReset()
        self.upd()

    def setXax(self, xax):
        self.root.addMacro(['setXax', (xax, self.axes - self.data.ndim())])
        self.data.setXax(xax, self.axes)
        self.upd()
        self.showFid()
        self.plotReset()

    def setAxType(self, val, update=True):
        oldAxMult = self.getAxMult(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"], self.freq(), self.ref())
        if val == 3:
            self.viewSettings["ppm"] = True
        else:
            self.viewSettings["ppm"] = False
            self.viewSettings["axType"] = val
        newAxMult = self.getAxMult(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"], self.freq(), self.ref())
        self.xminlim = self.xminlim * newAxMult / oldAxMult
        self.xmaxlim = self.xmaxlim * newAxMult / oldAxMult
        if update:
            self.showFid()

    def hilbert(self):
        self.root.addMacro(['hilbert', (self.axes - self.data.ndim(), )])
        self.data.hilbert(self.axes)
        self.upd()
        self.showFid()

    def getColorMap(self):
        return COLORMAPLIST.index(self.viewSettings["colorMap"])

    def setColorMap(self, num):
        self.viewSettings["colorMap"] = COLORMAPLIST[num]

    def setColor(self, color):
        self.viewSettings["color"] = color

    def setLw(self, lw):
        self.viewSettings["linewidth"] = lw

    def setContourColors(self, colors):
        self.viewSettings["contourColors"] = colors

    def setContourConst(self, constant):
        self.viewSettings["contourConst"] = constant

    def getOOM(self):
        absVal = np.max(np.abs(self.data.getHyperData(0)))
        if absVal == 0.0:
            return 1
        else:
            return int(np.floor(np.log10(absVal)))

    def showFid(self, oldData=None, extraX=None, extraY=None, extraColor=None):  # display the 1D data
        self.peakPickReset()
        tmpdata = self.data1D.getHyperData(0)
        self.ax.cla()
        axMult = self.getAxMult(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"], self.freq(), self.ref())
        self.line_xdata = [self.xax() * axMult]
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        if self.len() == 1:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = self.MARKER
            linestyle = self.LINESTYLE
        if oldData is not None:
            tmp = np.real(self.getDataType(oldData.getHyperData(0)))
            self.line_xdata_extra.append(self.xax() * axMult)
            self.line_ydata_extra.append(tmp)
            self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, c='k', alpha=0.2, linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if (extraX is not None):
            for num in range(len(extraX)):
                self.line_xdata_extra.append(extraX[num] * axMult)
                self.line_ydata_extra.append(extraY[num])
                self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker='', linestyle='-', c=extraColor, linewidth=self.viewSettings["linewidth"], picker=True)
        tmpdata = self.getDataType(tmpdata)
        if(self.viewSettings["plotType"] == 2):
            self.line_xdata.append(self.line_xdata[-1])
            self.line_ydata = [np.imag(tmpdata), np.real(tmpdata)]
            self.ax.plot(self.line_xdata[-2], self.line_ydata[-2], marker=marker, linestyle=linestyle, c='r', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
        else:
            self.line_ydata = [np.real(tmpdata)]
        self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
        self.ax.set_xlabel(self.getLabel(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"]))
        if self.logx:
            self.ax.set_xscale('log')
        else:
            self.ax.set_xscale('linear')
            self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.logy:
            self.ax.set_yscale('log')
        else:
            self.ax.set_yscale('linear')
            self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        if self.spec() > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()

    def plotReset(self, xReset=True, yReset=True):  # set the plot limits to min and max values
        showDat = self.data1D.data[0]
        miny = np.min(self.line_ydata)
        maxy = np.max(self.line_ydata)
        for line in self.line_ydata_extra:
            miny = min(miny, min(line))
            maxy = max(maxy, max(line))
        if miny == maxy: # Prevents setting the limits equal
            miny -= 0.01
            maxy += 0.01
        if type(self) is CurrentContour:
            differ = 0
        else:
            differ = 0.05 * (maxy - miny)  # amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny - differ
            self.ymaxlim = maxy + differ
        if self.ref() == 0.0:
            self.viewSettings["ppm"] = False
        minx = np.min(self.line_xdata)
        maxx = np.max(self.line_xdata)
        for line in self.line_xdata_extra:
            minx = min(minx, min(line))
            maxx = max(maxx, max(line))
        if minx == maxx: # Prevents setting the limits equal
            minx -= 0.01
            maxx += 0.01
        if xReset:
            self.xminlim = minx
            self.xmaxlim = maxx
        if self.spec() > 0 and type(self) is not CurrentArrayed:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        if type(self) is CurrentContour: #If contour: reverse y-axis
            if  self.spec(-2) > 0:
                self.ax.set_ylim(self.ymaxlim, self.yminlim)
        self.canvas.draw()


#########################################################################################################
class CurrentScatter(Current1D):

    X_RESIZE = False
    Y_RESIZE = False

    MARKER = 'o'
    LINESTYLE = 'none'


#########################################################################################################
# the class from which the data of multiple spectra is displayed, the operations which only edit the content of this class are for previewing


class CurrentMulti(Current1D):

    X_RESIZE = False
    Y_RESIZE = True

    def setExtraSlice(self, extraNum, axes, locList):  # change the slice
        self.viewSettings["extraAxes"][extraNum] = axes
        self.viewSettings["extraLoc"][extraNum] = locList

    def copyCurrent(self, root, fig, canvas, data):
        return CurrentMulti(root, fig, canvas, data, self)

    def addExtraData(self, data, name):
        self.viewSettings["extraName"].append(name)
        self.viewSettings["extraData"].append(data)
        self.viewSettings["extraLoc"].append([0] * (len(self.viewSettings["extraData"][-1].shape()) - 1))
        self.viewSettings["extraColor"].append(COLORCONVERTER.to_rgb(COLORCYCLE[np.mod(len(self.viewSettings["extraData"]), len(COLORCYCLE))]['color']))  # find a good color system
        self.viewSettings["extraAxes"].append(len(data.shape()) - 1)
        self.viewSettings["extraScale"].append(1.0)
        self.viewSettings["extraOffset"].append(0.0)
        self.viewSettings["extraShift"].append(0.0)
        self.showFid()

    def delExtraData(self, num):
        del self.viewSettings["extraData"][num]
        del self.viewSettings["extraLoc"][num]
        del self.viewSettings["extraColor"][num]
        del self.viewSettings["extraName"][num]
        del self.viewSettings["extraAxes"][num]
        del self.viewSettings["extraScale"][num]
        del self.viewSettings["extraOffset"][num]
        del self.viewSettings["extraShift"][num]
        self.showFid()

    def setExtraColor(self, num, color):
        self.viewSettings["extraColor"][num] = color
        self.showFid()

    def getExtraColor(self, num):
        return tuple(np.array(255 * np.array(self.viewSettings["extraColor"][num]), dtype=int))

    def resetLocList(self):
        super(CurrentMulti, self).resetLocList()
        self.resetExtraLocList()

    def setExtraScale(self, num, scale):
        self.viewSettings["extraScale"][num] = scale
        self.showFid()

    def setExtraOffset(self, num, offset):
        self.viewSettings["extraOffset"][num] = offset
        self.showFid()

    def setExtraShift(self, num, shift):
        self.viewSettings["extraShift"][num] = shift
        self.showFid()

    def resetExtraLocList(self, num=None):
        if num is None:
            for i in range(len(self.viewSettings["extraLoc"])):
                self.viewSettings["extraLoc"][i] = [0] * (len(self.viewSettings["extraData"][i].shape()) - 1)
        else:
            self.viewSettings["extraLoc"][num] = [0] * (len(self.viewSettings["extraData"][num].shape()) - 1)

    def showFid(self, oldData=None, extraX=None, extraY=None, extraColor=None):  # display the 1D data
        self.peakPickReset()
        tmpdata = self.data1D.getHyperData(0)
        self.ax.cla()
        axMult = self.getAxMult(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"], self.freq(), self.ref())
        self.line_xdata = [self.xax() * axMult]
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        for i in range(len(self.viewSettings["extraData"])):
            data = self.viewSettings["extraData"][i]
            try:
                if self.viewSettings["extraData"][i].ndim() <= self.viewSettings["extraAxes"][i]:
                    self.viewSettings["extraAxes"][i] = len(self.viewSettings["extraData"][i].shape()) - 1
                    self.resetExtraLocList(i)
                extraData1D = data.getSlice(self.viewSettings["extraAxes"][i], self.viewSettings["extraLoc"][i])
            except Exception:
                self.resetExtraLocList(i)
                extraData1D = data.getSlice(self.viewSettings["extraAxes"][i], self.viewSettings["extraLoc"][i])
            xax = extraData1D.xaxArray[-1]
            spec = extraData1D.spec[-1]
            freq = extraData1D.freq[-1]
            if extraData1D.ref[-1] is not None:
                ref = extraData1D.ref[-1]
            else:
                ref = extraData1D.freq[-1]
            self.line_xdata_extra.append(self.viewSettings["extraShift"][i] + xax * self.getAxMult(spec, self.viewSettings["axType"], self.viewSettings["ppm"], freq, ref))
            if extraData1D.getHyperData(0).shape[-1] == 1:
                marker = 'o'
                linestyle = 'none'
            else:
                marker = self.MARKER
                linestyle = self.LINESTYLE
            extraData = np.real(self.getDataType(extraData1D.getHyperData(0)))
            self.line_ydata_extra.append(extraData * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i])
            self.ax.plot(self.line_xdata_extra[-1],
                         self.line_ydata_extra[-1],
                         marker=marker, linestyle=linestyle,
                         c=self.viewSettings["extraColor"][i],
                         linewidth=self.viewSettings["linewidth"],
                         label=data.name,
                         picker=True)
        if self.len() == 1:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = ''
            linestyle = '-'
        if oldData is not None:
            tmp = np.real(self.getDataType(oldData.getHyperData(0)))
            self.line_xdata_extra.append(self.xax() * axMult)
            self.line_ydata_extra.append(tmp)
            self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, c='k', alpha=0.2, linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if (extraX is not None):
            for num in range(len(extraX)):
                self.line_xdata_extra.append(extraX[num] * axMult)
                self.line_ydata_extra.append(extraY[num])
                self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], c=extraColor[num], picker=True)
        tmpdata = self.getDataType(tmpdata)
        if(self.viewSettings["plotType"] == 2):
            self.line_xdata.append(self.line_xdata[-1])
            self.line_ydata = [np.imag(tmpdata), np.real(tmpdata)]
            self.ax.plot(self.line_xdata[-2], self.line_ydata[-2], marker=marker, linestyle=linestyle, c='r', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
        else:
            self.line_ydata = [np.real(tmpdata)]
        self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
        self.ax.set_xlabel(self.getLabel(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"]))
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        if self.spec() > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()

#########################################################################################################
# the class from which the stacked data is displayed, the operations which only edit the content of this class are for previewing


class CurrentStacked(Current1D):

    X_RESIZE = False
    Y_RESIZE = True

    def startUp(self, xReset=True, yReset=True):
        self.resetSpacing()
        self.showFid()
        self.plotReset(xReset, yReset)

    def copyCurrent(self, root, fig, canvas, data):
        return CurrentStacked(root, fig, canvas, data, self)

    def upd(self):  # get new data from the data instance
        if self.data.ndim() < 2:
            self.root.rescue()
            return False
        if (len(self.locList) + 2) != self.data.ndim():
            self.resetLocList()
        try:
            self.data1D = self.data.getBlock(self.axes, self.axes2, self.locList, self.viewSettings["stackBegin"], self.viewSettings["stackEnd"], self.viewSettings["stackStep"])
            if updateVar is None:
                self.root.rescue()
                return False
        except Exception:
            self.resetLocList()
            self.data1D = self.data.getBlock(self.axes, self.axes2, self.locList, self.viewSettings["stackBegin"], self.viewSettings["stackEnd"], self.viewSettings["stackStep"])
        return True

    def setBlock(self, axes, axes2, locList, stackBegin=None, stackEnd=None, stackStep=None):  # change the slice
        axesSame = True
        if (self.axes != axes) or (self.axes2 != axes2):
            axesSame = False
        self.axes = axes
        self.axes2 = axes2
        self.viewSettings["stackBegin"] = stackBegin
        self.viewSettings["stackEnd"] = stackEnd
        self.viewSettings["stackStep"] = stackStep
        if stackStep is None:
            if self.data.shape()[self.axes2] > 100:
                self.viewSettings["stackStep"] = 1 + int(self.data.shape()[self.axes2]) / 100
        self.locList = locList
        self.upd()
        self.showFid()
        if not axesSame:
            self.resetSpacing()
            self.plotReset()

    def resetLocList(self):
        self.locList = [0] * (self.data.ndim() - 2)

    def stackSelect(self, stackBegin, stackEnd, stackStep):
        self.viewSettings["stackBegin"] = stackBegin
        self.viewSettings["stackEnd"] = stackEnd
        self.viewSettings["stackStep"] = stackStep
        self.upd()
        self.showFid()
        self.plotReset(self.X_RESIZE, self.Y_RESIZE)

    def setSpacing(self, spacing):
        self.viewSettings["spacing"] = spacing
        self.showFid()
        self.plotReset(self.X_RESIZE, self.Y_RESIZE)

    def resetSpacing(self):
        difference = np.diff(self.data1D.getHyperData(0), axis=0)
        if difference.size == 0:
            self.viewSettings["spacing"] = 0
        else:
            difference = self.getDataType(difference)
            tmpData = self.getDataType(self.data1D.getHyperData(0))
            difference = np.min((np.real(difference), np.imag(difference)))
            amp = np.max((np.real(tmpData), np.imag(tmpData))) - np.min((np.real(tmpData), np.imag(tmpData)))
            self.viewSettings["spacing"] = np.abs(difference) + 0.1 * amp

    def altScroll(self, event):
        self.viewSettings["spacing"] = self.viewSettings["spacing"] * 1.1**event.step
        self.root.sideframe.scrollSpacing(self.viewSettings["spacing"])
        self.showFid()

    def altReset(self):
        self.resetSpacing()
        self.root.sideframe.scrollSpacing(self.viewSettings["spacing"])
        self.showFid()

    def showFid(self, oldData=None, extraX=None, extraY=None, extraColor=None):  # display the 1D data
        self.peakPickReset()
        tmpdata = self.data1D.getHyperData(0)
        self.ax.cla()
        axMult = self.getAxMult(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"], self.freq(), self.ref())
        tmp_line_xdata = self.xax() * axMult
        self.line_xdata = []
        self.line_ydata = []
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        if self.len() == 1:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = self.MARKER
            linestyle = self.LINESTYLE
        if oldData is not None:
            for num in range(len(oldData.getHyperData(0))):
                tmp = np.real(self.getDataType(oldData.getHyperData(0)[num]))
                self.line_xdata_extra.append(tmp_line_xdata)
                self.line_ydata_extra.append(num * self.viewSettings["spacing"] + tmp)
                self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, c='k', alpha=0.2, linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if (extraX is not None):
            tmpx = extraX[0] * axMult
            for num in range(len(extraY)):
                self.line_xdata_extra.append(tmpx)
                self.line_ydata_extra.append(num * self.viewSettings["spacing"] + extraY[num])
                self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], c=extraColor[0], picker=True)
        tmpdata = self.getDataType(tmpdata)
        for num in range(len(tmpdata)):
            if (self.viewSettings["plotType"] == 2):
                self.line_xdata.append(tmp_line_xdata)
                self.line_ydata.append(num * self.viewSettings["spacing"] + np.imag(tmpdata[num]))
                self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c='r', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
            self.line_xdata.append(tmp_line_xdata)
            self.line_ydata.append(num * self.viewSettings["spacing"] + np.real(tmpdata[num]))
            self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
        self.ax.set_xlabel(self.getLabel(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"]))
        if self.spec() > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.canvas.draw()

#########################################################################################################
# the class from which the arrayed data is displayed, the operations which only edit the content of this class are for previewing

class CurrentArrayed(CurrentStacked):

    X_RESIZE = True
    Y_RESIZE = False

    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        if duplicateCurrent is not None:
            if isinstance(duplicateCurrent, CurrentArrayed):
                self.zminlim = duplicateCurrent.zminlim
                self.zmaxlim = duplicateCurrent.zmaxlim
            else:
                # The z-axes limits are in xax units unlike the x-axes and y-axes limits
                axMult = self.getAxMult(duplicateCurrent.spec(),
                                        duplicateCurrent.viewSettings["axType"],
                                        duplicateCurrent.viewSettings["ppm"],
                                        duplicateCurrent.freq(),
                                        duplicateCurrent.ref())
                self.zminlim = (duplicateCurrent.xminlim) / axMult
                self.zmaxlim = (duplicateCurrent.xmaxlim) / axMult
        super(CurrentArrayed, self).__init__(root, fig, canvas, data, duplicateCurrent)

    def startUp(self, xReset=True, yReset=True):
        self.resetSpacing(False)
        self.showFid()
        self.plotReset(xReset, yReset)

    def copyCurrent(self, root, fig, canvas, data):
        return CurrentArrayed(root, fig, canvas, data, self)

    def setAxType2(self, val, update = True):
        oldAxMult = self.getAxMult(self.spec(-2), self.viewSettings["axType2"], self.viewSettings["ppm2"], self.freq(-2), self.ref(-2))
        if val == 3:
            self.viewSettings["ppm2"] = True
        else:
            self.viewSettings["ppm2"] = False
            self.viewSettings["axType2"] = val
        newAxMult = self.getAxMult(self.spec(-2), self.viewSettings["axType2"], self.viewSettings["ppm2"], self.freq(-2), self.ref(-2))
        if update:
            self.showFid()

    def resetSpacing(self, zlims=True):
        if zlims:
            self.zminlim = min(self.xax())
            self.zmaxlim = max(self.xax())
        xaxZlims = (self.xax() > self.zminlim) & (self.xax() < self.zmaxlim)
        self.viewSettings["spacing"] = (self.xax()[xaxZlims][-1] - self.xax()[xaxZlims][0]) * 1.1

    def showFid(self, oldData=None, extraX=None, extraY=None, extraColor=None):  # display the 1D data
        self.peakPickReset()
        tmpdata = self.data1D.getHyperData(0)
        self.ax.cla()
        if self.spec() > 0:
            direc = slice(None, None, -1)
        else:
            direc = slice(None, None, 1)
        xaxZlims = (self.xax() > self.zminlim) & (self.xax() < self.zmaxlim)
        axMult = self.getAxMult(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"], self.freq(), self.ref())
        axMult2 = self.getAxMult(self.spec(-2), self.viewSettings["axType2"], self.viewSettings["ppm2"], self.freq(-2), self.ref(-2))
        self.line_xdata = []
        self.line_ydata = []
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        if self.len() == 1:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = self.MARKER
            linestyle = self.LINESTYLE
        if oldData is not None:
            for num in range(len(oldData.getHyperData(0))):
                tmp = np.real(self.getDataType(oldData.getHyperData(0)[num]))
                self.line_xdata_extra.append((num * self.viewSettings["spacing"] + self.xax()[xaxZlims]) * axMult)
                self.line_ydata_extra.append(tmp[xaxZlims][direc])
                self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, c='k', alpha=0.2, linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if (extraX is not None):
            extraZlims = (extraX[0] > self.zminlim) & (extraX[0] < self.zmaxlim)
            for num in range(len(extraY)):
                self.line_xdata_extra.append((num * self.viewSettings["spacing"] + extraX[0][extraZlims]) * axMult)
                self.line_ydata_extra.append(extraY[num][extraZlims][direc])
                self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], c=extraColor[0], picker=True)
        tmpdata = self.getDataType(tmpdata)
        ticksPos = []
        for num in range(len(tmpdata)):
            if (self.viewSettings["plotType"] == 2):
                self.line_xdata.append((num * self.viewSettings["spacing"] + self.xax()[xaxZlims]) * axMult)
                self.line_ydata.append(np.imag(tmpdata[num][xaxZlims])[direc])
                self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c='r', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
            self.line_xdata.append((num * self.viewSettings["spacing"] + self.xax()[xaxZlims]) * axMult)
            self.line_ydata.append(np.real(tmpdata[num][xaxZlims])[direc])
            self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
            pos = (num * self.viewSettings["spacing"] + 0.5 * (self.xax()[xaxZlims][-1] + self.xax()[xaxZlims][0])) * axMult
            ticksPos.append(pos)
        self.ax.set_xticks(ticksPos)
        self.ax.set_xticklabels([('%#.3g') % x for x in self.xax(-2) * axMult2])
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.canvas.draw()

######################################################################################################


def add_diagonal(axes, mult, *line_args, **line_kwargs):
    identity, = axes.plot([], [], *line_args, **line_kwargs)

    def callback(axes):
        low_x, high_x = axes.get_xlim()
        low_y, high_y = axes.get_ylim()
        if mult == 0:
            low = low_y
            high = high_y
        else:
            low_y = low_y / float(mult)
            high_y = high_y / float(mult)
            low = max(low_x, low_y)
            high = min(high_x, high_y)
        identity.set_data([low, high], [low * mult, high * mult])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes

#########################################################################################################
# the class from which the contour data is displayed, the operations which only edit the content of this class are for previewing


class CurrentContour(CurrentStacked):

    X_RESIZE = False
    Y_RESIZE = True

    def startUp(self, xReset=True, yReset=True):
        self.showFid()
        self.plotReset(xReset, yReset)
    
    def altScroll(self, event):  # Shift scroll scrolls contour limits
        minLevels = self.viewSettings["minLevels"] / 1.1**event.step
        if minLevels > 1:
            minLevels = 1
        if self.viewSettings["maxLevels"] > 1:
            self.viewSettings["maxLevels"] = 1
        self.viewSettings["minLevels"] = minLevels
        self.root.sideframe.minLEntry.setText(format(self.viewSettings["minLevels"] * 100, '.7g'))
        self.root.sideframe.maxLEntry.setText(str(self.viewSettings["maxLevels"] * 100))
        self.plotContour(updateOnly=True)

    def copyCurrent(self, root, fig, canvas, data):
        return CurrentContour(root, fig, canvas, data, self)

    def setBlock(self, axes, axes2, locList):  # change the slice
        axesSame = True
        if (self.axes != axes) or (self.axes2 != axes2):
            axesSame = False
        self.axes = axes
        self.axes2 = axes2
        self.locList = locList
        self.upd()
        self.showFid()
        if not axesSame:
            self.plotReset()

    def setLevels(self, numLevels, maxLevels, minLevels, limitType, contourSign, contourType, multiValue):
        self.viewSettings["numLevels"] = numLevels
        self.viewSettings["maxLevels"] = maxLevels
        self.viewSettings["minLevels"] = minLevels
        self.viewSettings["limitType"] = limitType
        self.viewSettings["contourSign"] = contourSign
        self.viewSettings["contourType"] = contourType
        self.viewSettings["multiValue"] = multiValue
        self.showFid()

    def setProjLimits(self, ProjBool, Limits):
        self.viewSettings["projLimits"] = Limits
        self.viewSettings["projLimitsBool"] = ProjBool

    def setProjPos(self, pos):
        self.viewSettings["projPos"] = pos
    
    def setAxType2(self, val, update = True):
        oldAxMult = self.getAxMult(self.spec(-2), self.viewSettings["axType2"], self.viewSettings["ppm2"], self.freq(-2), self.ref(-2))
        if val == 3:
            self.viewSettings["ppm2"] = True
        else:
            self.viewSettings["ppm2"] = False
            self.viewSettings["axType2"] = val
        newAxMult = self.getAxMult(self.spec(-2), self.viewSettings["axType2"], self.viewSettings["ppm2"], self.freq(-2), self.ref(-2))
        self.yminlim = self.yminlim * newAxMult / oldAxMult
        self.ymaxlim = self.ymaxlim * newAxMult / oldAxMult
        if update:
            self.updateAxes(oldAxMult, newAxMult, 1)

    def setProjType(self, val, direc):
        if direc == 1:
            self.viewSettings["projTop"] = val
        if direc == 2:
            self.viewSettings["projRight"] = val

    def setProjTraces(self, val, direc):
        self.viewSettings["projPos"][direc] = val

    def updateAxes(self,oldAx, newAx, axis):
        scale = newAx / oldAx
        # Scale the path vertices, so no new contours need to be calculated
        cols = self.ax.collections
        for col in cols:
            paths = col.get_paths()
            for path in paths:
                tmp = path.vertices
                tmp[:,axis] = tmp[:,axis] * scale
                path.vertices = tmp
        # Scale the projections
        if axis == 1: # Yaxis
            line = self.y_ax.lines
            line[0].set_ydata(line[0].get_ydata() * scale)
        else:
            line = self.x_ax.lines
            line[0].set_xdata(line[0].get_xdata() * scale)
        # Set the labels
        self.ax.set_xlabel(self.getLabel(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"]))
        self.ax.set_ylabel(self.getLabel(self.spec(-2), self.viewSettings["axType2"], self.viewSettings["ppm2"]))
        # Set the zoom
        if axis == 1:
            ylim = self.ax.get_ylim()
            self.ax.set_ylim(ylim[0] * scale, ylim[1] * scale)
            self.line_ydata = [item*scale for item in self.line_ydata]
        else:
            xlim = self.ax.get_xlim()
            self.ax.set_xlim(xlim[0] * scale, xlim[1] * scale)
            self.line_xdata = [item*scale for item in self.line_xdata]
        self.canvas.draw()

    def showFid(self, oldData=None, extraX=None, extraY=None, extraColor=None, makeContours=True, clearCntr=True):
        # The oldData and extra plots are not displayed in the contourplot for now
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        self.differ = None
        self.peakPickReset()
        tmpdata = self.data1D.getHyperData(0)
        if clearCntr:
            self.ax.cla()
        self.x_ax.cla()
        self.y_ax.cla()
        if self.viewSettings["diagonalBool"]:
            add_diagonal(self.ax, self.viewSettings["diagonalMult"], c='k', ls='--')
        self.line_xdata = [self.xax() * self.getAxMult(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"], self.freq(), self.ref())]
        self.line_ydata = [self.xax(-2) * self.getAxMult(self.spec(-2), self.viewSettings["axType2"], self.viewSettings["ppm2"], self.freq(-2), self.ref(-2))]
        tmpdata = np.real(self.getDataType(tmpdata))
        self.line_zdata = [tmpdata]
        if self.viewSettings["limitType"] == 0:
            self.differ = np.max(np.abs(tmpdata))
        else:
            self.differ = np.max(np.abs(np.ravel(self.data.getHyperData)))
        if makeContours:
            self.plotContour()
        self.showProj()
        self.ax.set_xlabel(self.getLabel(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"]))
        self.ax.set_ylabel(self.getLabel(self.spec(-2), self.viewSettings["axType2"], self.viewSettings["ppm2"]))
        if self.spec():
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        if self.spec(-2):
            self.ax.set_ylim(self.ymaxlim, self.yminlim)
        else:
            self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.x_ax.get_yaxis().get_major_formatter().set_powerlimits((-2, 2))
        self.y_ax.get_xaxis().get_major_formatter().set_powerlimits((-2, 2))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.canvas.draw()

    def plotContour(self, updateOnly=False):  # Plots the contour plot
        X, Y = np.meshgrid(self.line_xdata[-1], self.line_ydata[-1])
        if updateOnly:  # Set some extra stuff if only the contour plot needs updating
            del self.ax.collections[:]  # Clear all plot collections
        if self.viewSettings["contourType"] == 0:  # if linear
            contourLevels = np.linspace(self.viewSettings["minLevels"] * self.differ, self.viewSettings["maxLevels"] * self.differ, self.viewSettings["numLevels"])
        elif self.viewSettings["contourType"] == 1:  # if Multiplier
            contourLevels = [self.viewSettings["minLevels"] * self.differ]
            while contourLevels[-1] < self.viewSettings["maxLevels"] * self.differ and len(contourLevels) < self.viewSettings["numLevels"]:
                contourLevels.append(contourLevels[-1] * self.viewSettings["multiValue"])
            contourLevels = np.array(contourLevels)
        # Trim matrix of unused rows/columns for more efficient contour plotting
        PlotPositive = False
        if self.viewSettings["contourSign"] == 0 or self.viewSettings["contourSign"] == 1:
            if self.line_zdata[-1].shape[0] > 2:  # if size 2 or lower, convolve fails, just take whole data then
                YposMax = np.where(np.convolve(np.max(self.line_zdata[-1], 1) > contourLevels[0], [True, True, True], 'same'))[0]
            else:
                YposMax = np.arange(self.line_zdata[-1].shape[0])
            if YposMax.size > 0:  # if number of positive contours is non-zero
                if self.line_zdata[-1].shape[1] > 2:
                    XposMax = np.where(np.convolve(np.max(self.line_zdata[-1], 0) > contourLevels[0], [True, True, True], 'same'))[0]
                else:
                    XposMax = np.arange(self.line_zdata[-1].shape[1])
                PlotPositive = True
        PlotNegative = False
        if self.viewSettings["contourSign"] == 0 or self.viewSettings["contourSign"] == 2:
            if not self.viewSettings["plotType"] == 3:  # for Absolute plot no negative
                if self.line_zdata[-1].shape[0] > 2:
                    YposMin = np.where(np.convolve(np.min(self.line_zdata[-1], 1) < -contourLevels[0], [True, True, True], 'same'))[0]
                else:
                    YposMin = np.arange(self.line_zdata[-1].shape[0])
                if YposMin.size > 0:  # if number of negative contours is non-zero
                    if self.line_zdata[-1].shape[1] > 2:
                        XposMin = np.where(np.convolve(np.min(self.line_zdata[-1], 0) < -contourLevels[0], [True, True, True], 'same'))[0]
                    else:
                        XposMin = np.arange(self.line_zdata[-1].shape[1])
                    PlotNegative = True
        vmax = max(np.abs(self.viewSettings["minLevels"] * self.differ), np.abs(self.viewSettings["maxLevels"] * self.differ))
        vmin = -vmax
        if self.viewSettings["contourConst"]:
            if PlotPositive:
                self.ax.contour(X[YposMax[:,None],XposMax],Y[YposMax[:,None],XposMax],self.line_zdata[-1][YposMax[:,None],XposMax], colors=self.viewSettings["contourColors"][0], levels=contourLevels, vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], label=self.data.name, linestyles='solid')
            if PlotNegative:
                self.ax.contour(X[YposMin[:,None],XposMin],Y[YposMin[:,None],XposMin],self.line_zdata[-1][YposMin[:,None],XposMin], colors=self.viewSettings["contourColors"][1], levels=-contourLevels[::-1], vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], linestyles='solid')
        else:
            if PlotPositive:
                self.ax.contour(X[YposMax[:,None],XposMax],Y[YposMax[:,None],XposMax],self.line_zdata[-1][YposMax[:,None],XposMax], cmap=get_cmap(self.viewSettings["colorMap"]), levels=contourLevels, vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], label=self.data.name, linestyles='solid')
            if PlotNegative:    
                self.ax.contour(X[YposMin[:,None],XposMin],Y[YposMin[:,None],XposMin],self.line_zdata[-1][YposMin[:,None],XposMin], cmap=get_cmap(self.viewSettings["colorMap"]), levels=-contourLevels[::-1], vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], linestyles='solid')
        if updateOnly:
            self.canvas.draw()

    def showProj(self):
        xLimOld = self.x_ax.get_xlim()
        x = self.line_xdata[-1]  # Get plot data from plot
        yLimOld = self.y_ax.get_ylim()
        y = self.line_ydata[-1]  # Get plot data from plot
        self.x_ax.cla()
        self.y_ax.cla()
        tmpdata = np.real(self.getDataType(self.data1D.getHyperData(0)))
        Limits = self.viewSettings["projLimits"]
        topSlice = (slice(None), slice(None))
        rightSlice = (slice(None), slice(None))
        if self.viewSettings["projLimitsBool"] is True:
            if Limits[0] < Limits[1]:
                topSlice = (slice(Limits[0], Limits[1] + 1), slice(None))
            elif Limits[0] > Limits[1]:
                topSlice = (slice(Limits[1], Limits[0] + 1), slice(None))
            else:
                topSlice = (slice(Limits[0], Limits[0] + 1), slice(None))
            if Limits[2] < Limits[3]:
                rightSlice = (slice(None), slice(Limits[2], Limits[3] + 1))
            elif Limits[2] > Limits[3]:
                rightSlice = (slice(None), slice(Limits[3], Limits[2] + 1))
            else:
                rightSlice = (slice(None), slice(Limits[2], Limits[2] + 1))
        if self.viewSettings["projTop"] == 0:
            xprojdata = np.sum(tmpdata[topSlice], axis=0)
        elif self.viewSettings["projTop"] == 1:
            xprojdata = np.max(tmpdata[topSlice], axis=0)
        elif self.viewSettings["projTop"] == 2:
            xprojdata = np.min(tmpdata[topSlice], axis=0)
        elif self.viewSettings["projTop"] == 4:
            if self.viewSettings["projPos"][0] >= self.data.shape()[self.axes2]:
                self.viewSettings["projPos"][0] = self.data.shape()[self.axes2] - 1
            elif self.viewSettings["projPos"][0] < 0:
                self.viewSettings["projPos"][0] = 0
            xprojdata = tmpdata[self.viewSettings["projPos"][0]]
        if self.viewSettings["projTop"] != 3:
            self.x_ax.plot(x, xprojdata, color=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], picker=True)            
            xmin, xmax = np.min(xprojdata), np.max(xprojdata)
            self.x_ax.set_ylim([xmin - 0.15 * (xmax - xmin), xmax + 0.05 * (xmax - xmin)])  # Set projection limits, and force 15% whitespace below plot
            self.x_ax.set_xlim(xLimOld)
        if self.viewSettings["projRight"] == 0:
            yprojdata = np.sum(tmpdata[rightSlice], axis=1)
        elif self.viewSettings["projRight"] == 1:
            yprojdata = np.max(tmpdata[rightSlice], axis=1)
        elif self.viewSettings["projRight"] == 2:
            yprojdata = np.min(tmpdata[rightSlice], axis=1)
        elif self.viewSettings["projRight"] == 4:
            if self.viewSettings["projPos"][1] >= self.data.shape()[self.axes]:
                self.viewSettings["projPos"][1] = self.data.shape()[self.axes] - 1
            elif self.viewSettings["projPos"][1] < 0:
                self.viewSettings["projPos"][1] = 0
            yprojdata = tmpdata[:, self.viewSettings["projPos"][1]]
        if self.viewSettings["projRight"] != 3:
            self.y_ax.plot(yprojdata, y, color=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], picker=True)
            ymin, ymax = np.min(yprojdata), np.max(yprojdata)
            self.y_ax.set_xlim([ymin - 0.15 * (ymax - ymin), ymax + 0.05 * (ymax - ymin)])  # Set projection limits, and force 15% whitespace below plot
            self.y_ax.set_ylim(yLimOld)
        self.canvas.draw()

    # The peakpicking function needs to be changed for contour plots
    def buttonRelease(self, event):
        if event.button == 1:
            if self.peakPick:
                if self.rect[1] is not None:
                    self.rect[1].remove()
                    self.rect[1] = None
                if self.rect[0] is not None:
                    self.rect[0].remove()
                    self.rect[0] = None
                    self.peakPick = False
                    xdata = self.xax() * self.getAxMult(self.spec(), self.viewSettings["axType"], self.viewSettings["ppm"], self.freq(), self.ref())
                    ydata = self.xax(-2) * self.getAxMult(self.spec(-2), self.viewSettings["axType2"], self.viewSettings["ppm2"], self.freq(-2), self.ref(-2))
                    idx = np.argmin(np.abs(xdata - event.xdata))
                    idy = np.argmin(np.abs(ydata - event.ydata))
                    if self.peakPickFunc is not None:
                        tmpdata = np.real(self.getDataType(self.data1D.getHyperData(0)[idy, idx]))
                        self.peakPickFunc((idx, xdata[idx], tmpdata, idy, ydata[idy]))
                    if not self.peakPick:  # check if peakpicking is still required
                        self.peakPickFunc = None
            else:
                self.leftMouse = False
                if self.rect[0] is not None:
                    self.rect[0].remove()
                if self.rect[1] is not None:
                    self.rect[1].remove()
                if self.rect[2] is not None:
                    self.rect[2].remove()
                if self.rect[3] is not None:
                    self.rect[3].remove()
                self.rect = [None, None, None, None]
                if self.zoomX2 is not None and self.zoomY2 is not None:
                    self.xminlim = min([self.zoomX1, self.zoomX2])
                    self.xmaxlim = max([self.zoomX1, self.zoomX2])
                    self.yminlim = min([self.zoomY1, self.zoomY2])
                    self.ymaxlim = max([self.zoomY1, self.zoomY2])
                    if self.spec() > 0:
                        self.ax.set_xlim(self.xmaxlim, self.xminlim)
                    else:
                        self.ax.set_xlim(self.xminlim, self.xmaxlim)
                    if self.spec(-2) > 0:
                        self.ax.set_ylim(self.ymaxlim, self.yminlim)
                    else:
                        self.ax.set_ylim(self.yminlim, self.ymaxlim)
                self.zoomX1 = None
                self.zoomX2 = None
                self.zoomY1 = None
                self.zoomY2 = None
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw()

