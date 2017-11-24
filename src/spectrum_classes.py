#!/usr/bin/env python

# Copyright 2016 - 2017 Bas van Meerten and Wouter Franssen

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
from mpl_toolkits.mplot3d import proj3d
from six import string_types
from spectrumFrame import Plot1DFrame
from safeEval import safeEval
from nus import ffm, clean, ist
import multiprocessing
from matplotlib.pyplot import get_cmap
import matplotlib
import matplotlib._cntr as cntr
import matplotlib.collections as mcoll
import reimplement as reim
import functions as func

COLORMAPLIST = ['seismic', 'BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'rainbow', 'jet']
COLORCYCLE = list(matplotlib.rcParams['axes.prop_cycle'])
COLORCONVERTER = matplotlib.colors.ColorConverter()

#########################################################################
# the generic data class


class Spectrum(object):

    def __init__(self, name, data, filePath, freq, sw, spec=None, wholeEcho=None, hyper = None, ref=None, xaxArray=None, history=None, msgHandler=None):
        self.name = name
        if isinstance(data, (list)):
            self.data = []
            for item in data:
                self.data.append(np.array(item, dtype=complex))
        else:
            self.data = [np.array(data, dtype=complex)]  # data of dimension dim
        if hyper is None:
            self.hyper = [] #Holds the axes where hypercomplex data exists
        else:
            self.hyper = hyper
        self.filePath = filePath
        self.freq = np.array(freq)  # array of center frequency (length is dim, MHz)
        self.sw = sw  # array of sweepwidths
        self.noUndo = False
        if spec is None:
            self.spec = [0] * self.data.ndim
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
        return self.data[0].ndim

    def shape(self):
        return self.data[0].shape
        
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

    def checkAxes(self, axes):
        if axes < 0:
            axes = axes + self.ndim()
        if not (0 <= axes < self.ndim()):
            self.dispMsg('Not a valid axes')
            return None
        return axes

    def resetXax(self, axes=None):
        if axes is not None:
            axes = self.checkAxes(axes)
            if axes is None:
                return None
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
        if axes is None:
            return None
        if len(xax) != self.shape()[axes]:
            self.dispMsg("Length of new x-axis does not match length of the data")
            return None
        oldXax = self.xaxArray[axes]
        self.xaxArray[axes] = xax
        self.addHistory("X-axes of dimension " + str(axes + 1) + " was set to " + str(xax).replace('\n', ''))
        if self.noUndo:
            return None
        else:
            return lambda self: self.setXax(oldXax, axes)

    def insert(self, data, pos, axes, dataImag=None):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        for index in range(len(self.data)):
            data = np.array(data[index]) 
            if dataImag is not None:
                data += 1j * np.array(dataImag[index])
            oldSize = self.data[index].shape[axes]
            self.data[index] = np.insert(self.data[index], [pos], data, axis=axes)
        self.resetXax(axes)
        self.addHistory("Inserted " + str(self.shape()[axes] - oldSize) + " datapoints in dimension " + str(axes + 1) + " at position " + str(pos))
        if self.noUndo:
            return None
        else:
            return lambda self: self.remove(range(pos, pos + data.shape[axes]), axes)

    def remove(self, pos, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
            returnValue = lambda self: self.restoreData(copyData, lambda self: self.remove(pos, axes))
        tmpdata = []
        for index in range(len(self.data)):
            tmpdata.append(np.delete(self.data[index], pos, axes))
        if (np.array(tmpdata[0].shape) != 0).all():
            self.data = tmpdata
            self.xaxArray[axes] = np.delete(self.xaxArray[axes], pos)
            try:
                length = len(pos)
            except Exception:
                length = 1
            self.addHistory("Removed " + str(length) + " datapoints from dimension " + str(axes + 1) + " at position " + str(pos))
            if self.noUndo:
                return None
            else:
                return lambda self: self.restoreData(copyData, lambda self: self.remove(pos, axes))
        else:
            self.dispMsg('Cannot delete all data')
            return None

    def add(self, data, dataImag= None, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        try:
            for index in range(len(self.data)):
                if dataImag is not None:
                    data[index] = np.array(data[index]) + 1j * np.array(dataImag[index])
                self.data[index][select] = self.data[index][select] + data[index]
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.addHistory("Added to data[" + str(select) + "]")
        if self.noUndo:
            return None
        else:
            return lambda self: self.subtract(data, select=select)

    def subtract(self, data, dataImag=None, select=slice(None), singleHyper = False):
        if isinstance(select, string_types):
            select = safeEval(select)
        try:
            if singleHyper:
                if dataImag is not None:
                    data = np.array(data) + 1j * np.array(dataImag)
                self.data[0][select] = self.data[0][select] - data[0]
            else:
                for index in range(len(self.data)):
                    if dataImag is not None:
                        data[index] = np.array(data[index]) + 1j * np.array(dataImag[index])
                    self.data[index][select] = self.data[index][select] - data[index]
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.addHistory("Subtracted from data[" + str(select) + "]")
        if self.noUndo:
            return None
        else:
            return lambda self: self.add(data, select=select)
    
    def multiplySpec(self, data, dataImag=None, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        try:
            for index in range(len(self.data)):
                if dataImag is not None:
                    data[index] = np.array(data[index]) + 1j * np.array(dataImag[index])
                self.data[index][select] = self.data[index][select] * data[index]
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.addHistory("Multiplied with data[" + str(select) + "]")
        if self.noUndo:
            return None
        else:
            return lambda self: self.divideSpec(data, select=select)

    def divideSpec(self, data, dataImag=None, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        try:
            for index in range(len(self.data)):
                if dataImag is not None:
                    data[index] = np.array(data[index]) + 1j * np.array(dataImag[index])
                self.data[index][select] = self.data[index][select] / data[index]
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.addHistory("Divided by data[" + str(select) + "]")
        if self.noUndo:
            return None
        else:
            return lambda self: self.multiplySpec(data, select=select)

    def multiply(self, mult, axes, multImag=0, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        try:
            mult = np.array(mult) + 1j * np.array(multImag)
            if not self.noUndo:
                copyData = copy.deepcopy(self)
            for index in range(len(self.data)):
                self.data[index][select] = np.apply_along_axis(np.multiply, axes, self.data[index], mult)[select]
        except ValueError as error:
            self.dispMsg('Multiply: ' + str(error))
            return None
        self.addHistory("Multiplied dimension " + str(axes + 1) + " of data[" + str(select) + "]: " + str(mult).replace('\n', ''))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.multiply(mult, axes, select=select))

    def baselineCorrection(self, baseline, axes, baselineImag=0, select=slice(None)):
        hyperView = 0
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        try:
            baseline = np.array(baseline) + 1j * np.array(baselineImag)
            baselinetmp = baseline.reshape((1, ) * axes + (self.shape()[axes], ) + (1, ) * (self.ndim() - axes - 1))
            self.data[hyperView][select] = self.data[hyperView][select] - baselinetmp
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.addHistory("Baseline corrected dimension " + str(axes + 1) + " of data[" + str(select) + "]")
        return lambda self: self.baselineCorrection(-baseline, axes, select=select)

    def concatenate(self, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        splitVal = self.shape()[axes]
        try:
            for index in range(len(self.data)):
                self.data[index] = np.concatenate(self.data[index], axis=axes)
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.data, self.hyper = self.deleteHyper(axes,self.data,self.hyper) #Remove hypercomplex along the axis to be removed
        for i in range(len(self.hyper)):
            if self.hyper[i] > axes:
                self.hyper[i] += -1
        self.freq = np.delete(self.freq, axes)
        self.sw = np.delete(self.sw, axes)
        self.spec = np.delete(self.spec, axes)
        self.wholeEcho = np.delete(self.wholeEcho, axes)
        self.ref = np.delete(self.ref, axes)
        del self.xaxArray[0]
        # self.resetXax(axes)
        self.resetXax()
        self.addHistory("Concatenated dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.split(splitVal, axes)

    def split(self, sections, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        try:
            for index in range(len(self.data)):
                self.data[index] = np.array(np.split(self.data[index], sections, axis=axes))
        except ValueError as error:
            self.dispMsg('Split: ' + str(error))
            return None
        self.freq = np.insert(self.freq, 0, self.freq[axes])
        self.sw = np.insert(self.sw, 0, self.sw[axes])
        self.spec = np.insert(self.spec, 0, self.spec[axes])
        self.wholeEcho = np.insert(self.wholeEcho, 0, self.wholeEcho[axes])
        self.ref = np.insert(self.ref, 0, self.ref[axes])
        self.hyper = [x + 1 for x in self.hyper] #New dim always the new D1. All hyper values must be increased by 1
        self.xaxArray.insert(0, [])
        self.resetXax(0)
        self.resetXax(axes + 1)
        self.addHistory("Split dimension " + str(axes + 1) + " into " + str(sections) + " sections")
        if self.noUndo:
            return None
        else:
            return lambda self: self.concatenate(axes)

    def real(self):
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        for index in range(len(self.data)):
            self.data[index] = np.real(self.data[index])
        self.addHistory("Real")
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.real())

    def imag(self):
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        for index in range(len(self.data)):
            self.data[index] = np.imag(self.data[index])
        self.addHistory("Imaginary")
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.imag())

    def abs(self):
        if not self.noUndo:
            copyData = copy.deepcopy(self)
            returnValue = lambda self: self.restoreData(copyData, lambda self: self.abs())
        for index in range(len(self.data)):
            self.data[index] = np.abs(self.data[index])
        self.addHistory("Absolute")
        return returnValue
    
    def conj(self,axes):
        if self.noUndo:
            return None
        else:
            copyData = copy.deepcopy(self)

        self.data = self.hyperReorder(self.data, axes)
        for index in range(len(self.data)):
            self.data[index] = np.conj(self.data[index])
        self.data = self.hyperReorder(self.data, axes)
        self.addHistory("Complex conjugate dimension " + str(axes))
        return lambda self: self.restoreData(copyData, lambda self: self.abs())

    def states(self, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.shape()[axes] % 2 != 0:
            self.dispMsg("States: data has to be even")
            return None
        tmpdata = self.data
        slicing1 = (slice(None), ) * axes + (slice(None, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        slicing2 = (slice(None), ) * axes + (slice(1, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        self.data = []
        for index in range(len(tmpdata)):
            #tmp = np.real(tmpdata[index][slicing1]) + 1j * np.real(tmpdata[index][slicing2])
            #tmp2 = np.imag(tmpdata[index][slicing1]) + 1j * np.imag(tmpdata[index][slicing2])
            tmp = np.real(tmpdata[index][slicing1]) + 1j * np.imag(tmpdata[index][slicing1])
            tmp2 = np.real(tmpdata[index][slicing2]) + 1j * np.imag(tmpdata[index][slicing2])
            self.data.append(tmp)
            self.data.append(tmp2)
        self.hyper.append(axes)
        self.resetXax(axes)
        self.addHistory("States conversion on dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.states(axes))

    def statesTPPI(self, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.shape()[axes] % 2 != 0:
            self.dispMsg("States-TPPI: data has to be even")
            return None
        tmpdata = self.data
        slicing1 = (slice(None), ) * axes + (slice(None, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        slicing2 = (slice(None), ) * axes + (slice(1, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        self.data = []
        for index in range(len(tmpdata)):
            #tmp = np.real(tmpdata[index][slicing1]) + 1j * np.real(tmpdata[index][slicing2])
            #tmp[slicing2] = -1 * tmp[slicing2]
            #tmp2 = np.imag(tmpdata[index][slicing1]) + 1j * np.imag(tmpdata[index][slicing2])
            #tmp2[slicing2] = -1 * tmp2[slicing2]
            tmp = np.real(tmpdata[index][slicing1]) + 1j * np.imag(tmpdata[index][slicing1])
            tmp[slicing2] = -1 * tmp[slicing2]
            tmp2 = np.real(tmpdata[index][slicing2]) + 1j * np.imag(tmpdata[index][slicing2])
            tmp2[slicing2] = -1 * tmp2[slicing2]
            self.data.append(tmp)
            self.data.append(tmp2)
        self.hyper.append(axes)
        self.resetXax(axes)
        self.addHistory("States-TPPI conversion on dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.statesTPPI(axes))

    def echoAntiEcho(self, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.shape()[axes] % 2 != 0:
            self.dispMsg("Echo-antiecho: data has to be even")
            return None
        tmpdata = self.data
        slicing1 = (slice(None), ) * axes + (slice(None, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        slicing2 = (slice(None), ) * axes + (slice(1, None, 2), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        self.data = []
        for index in range(len(tmpdata)):
            #tmp1 = np.real(tmpdata[index][slicing1] + tmpdata[index][slicing2]) - 1j * np.imag(tmpdata[index][slicing1] - tmpdata[index][slicing2])
            #tmp2 = np.real(tmpdata[index][slicing1] - tmpdata[index][slicing2]) + 1j * np.imag(tmpdata[index][slicing1] + tmpdata[index][slicing2])
            tmp1 = np.real(tmpdata[index][slicing1] + tmpdata[index][slicing2]) + 1j * np.imag(tmpdata[index][slicing1] + tmpdata[index][slicing2])
            tmp2 =  - np.imag(tmpdata[index][slicing1] - tmpdata[index][slicing2]) + 1j * np.real(tmpdata[index][slicing1] - tmpdata[index][slicing2])
            self.data.append(tmp1)
            self.data.append(tmp2)
        self.hyper.append(axes)
        self.resetXax(axes)
        self.addHistory("Echo-antiecho conversion on dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.echoAntiEcho(axes))

    def subtractAvg(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not (0 <= pos1 <= self.shape()[axes]):
            self.dispMsg("Indices not within range")
            return None
        if not (0 <= pos2 <= self.shape()[axes]):
            self.dispMsg("Indices not within range")
            return None
        if pos1 == pos2:
            self.dispMsg("Indices cannot be equal")
            return None
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        for index in range(len(self.data)):
            averages = np.mean(self.data[index][slicing], axis=axes, keepdims=True)
            self.data[index] -= averages
        self.addHistory("Subtracted average determined between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.add(averages)

    def matrixManip(self, pos1, pos2, axes, which):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if isinstance(pos1, int):
            pos1 = np.array([pos1])
            pos2 = np.array([pos2])
        if len(pos1) != len(pos2):
            self.dispMsg("Length of the two arrays is not equal")
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        tmpdata = [() for x in range(2**len(self.hyper))]
        if len(pos1) == 1:
            keepdims = False
        else:
            keepdims = True
        for i in range(len(pos1)):
            if not (0 <= pos1[i] <= self.shape()[axes]):
                self.dispMsg("Indices not within range")
                return None
            if not (0 <= pos2[i] <= self.shape()[axes]):
                self.dispMsg("Indices not within range")
                return None
            if pos1[i] == pos2[i]:
                self.dispMsg("Indices cannot be equal")
                return None
            minPos = min(pos1[i], pos2[i])
            maxPos = max(pos1[i], pos2[i])
            slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), ) * (self.ndim() - 1 - axes)
            for index in range(len(self.data)):
                if which == 0:
                    if self.spec[axes] == 0:
                        tmpdata[index] += (np.sum(self.data[index][slicing], axis=axes, keepdims=keepdims) / self.sw[axes], )
                    else:
                        tmpdata[index] += (np.sum(self.data[index][slicing], axis=axes, keepdims=keepdims) * self.sw[axes] / (1.0 * self.shape()[axes]), )
                elif which == 5:
                    tmpdata[index] += (np.sum(self.data[index][slicing], axis=axes, keepdims=keepdims), )
                elif which == 1:
                    tmpdata[index] += (np.max(self.data[index][slicing], axis=axes, keepdims=keepdims), )
                elif which == 2:
                    tmpdata[index] += (np.min(self.data[index][slicing], axis=axes, keepdims=keepdims), )
                elif which == 3:
                    maxArgPos = np.argmax(np.real(self.data[index][slicing]), axis=axes)
                    tmpmaxPos = maxArgPos.flatten()
                    tmp = self.xaxArray[axes][slice(minPos, maxPos)][tmpmaxPos].reshape(maxArgPos.shape)
                    if keepdims:
                        tmpdata[index] += (np.expand_dims(tmp, axes), )
                    else:
                        tmpdata[index] += (tmp, )
                elif which == 4:
                    minArgPos = np.argmin(np.real(self.data[index][slicing]), axis=axes)
                    tmpminPos = minArgPos.flatten()
                    tmp = self.xaxArray[axes][slice(minPos, maxPos)][tmpminPos].reshape(minArgPos.shape)
                    if keepdims:
                        tmpdata[index] += (np.expand_dims(tmp, axes), )
                    else:
                        tmpdata[index] += (tmp, )
                elif which == 6:
                    tmpdata[index] += (np.mean(self.data[index][slicing], axis=axes, keepdims=keepdims), )
        if which == 0:
            self.addHistory("Integrate between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 5:
            self.addHistory("Sum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 1:
            self.addHistory("Maximum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 2:
            self.addHistory("Minimum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 3:
            self.addHistory("Maximum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 4:
            self.addHistory("Minimum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 6:
            self.addHistory("Average between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        #Remove hyper along this dim if necessary
        if keepdims ==  False and axes in self.hyper:
            tmpdata , self.hyper = self.deleteHyper(axes,tmpdata,self.hyper)
            #Correct hyper for missing dim 
            for i in range(len(self.hyper)):
                if self.hyper[i] > axes:
                    self.hyper[i] += -1
        if len(tmpdata[0]) == 1:
            if self.ndim() == 1:
                self.data = []
                for index in range(len(tmpdata)):
                    self.data.append(np.array([tmpdata[index][0]]))
                self.resetXax(axes)
            else:
                self.data = []
                for index in range(len(tmpdata)):
                    self.data.append(tmpdata[index][0])
                self.freq = np.delete(self.freq, axes)
                self.ref = np.delete(self.ref, axes)
                self.sw = np.delete(self.sw, axes)
                self.spec = np.delete(self.spec, axes)
                self.wholeEcho = np.delete(self.wholeEcho, axes)
                del self.xaxArray[axes]
        else:
            self.data = []
            for index in range(len(tmpdata)):
                self.data.append(np.concatenate(tmpdata[index], axis=axes))
            self.resetXax(axes)
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.matrixManip(pos1, pos2, axes, which))

    def matrixManipNew(self, pos1, pos2, axes, which):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if isinstance(pos1, int):
            pos1 = np.array([pos1])
            pos2 = np.array([pos2])
        if len(pos1) != len(pos2):
            self.dispMsg("Length of the two arrays is not equal")
            return None
        tmpdata = [() for x in range(2**len(self.hyper))]
        if len(pos1) == 1:
            keepdims = False
        else:
            keepdims = True
        for i in range(len(pos1)):
            if not (0 <= pos1[i] <= self.shape()[axes]):
                self.dispMsg("Indices not within range")
                return None
            if not (0 <= pos2[i] <= self.shape()[axes]):
                self.dispMsg("Indices not within range")
                return None
            if pos1[i] == pos2[i]:
                self.dispMsg("Indices cannot be equal")
                return None
            minPos = min(pos1[i], pos2[i])
            maxPos = max(pos1[i], pos2[i])
            slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), ) * (self.ndim() - 1 - axes)
            for index in range(len(self.data)):
                if which == 0:
                    if self.spec[axes] == 0:
                        tmpdata[index] += (np.sum(self.data[index][slicing], axis=axes, keepdims=keepdims) / self.sw[axes], )
                    else:
                        tmpdata[index] += (np.sum(self.data[index][slicing], axis=axes, keepdims=keepdims) * self.sw[axes] / (1.0 * self.shape()[axes]), )
                elif which == 5:
                    tmpdata[index] += (np.sum(self.data[index][slicing], axis=axes, keepdims=keepdims), )
                elif which == 1:
                    tmpdata[index] += (np.max(self.data[index][slicing], axis=axes, keepdims=keepdims), )
                elif which == 2:
                    tmpdata[index] += (np.min(self.data[index][slicing], axis=axes, keepdims=keepdims), )
                elif which == 3:
                    maxArgPos = np.argmax(np.real(self.data[index][slicing]), axis=axes)
                    tmpmaxPos = maxArgPos.flatten()
                    tmp = self.xaxArray[axes][slice(minPos, maxPos)][tmpmaxPos].reshape(maxArgPos.shape)
                    if keepdims:
                        tmpdata[index] += (np.expand_dims(tmp, axes), )
                    else:
                        tmpdata[index] += (tmp, )
                elif which == 4:
                    minArgPos = np.argmin(np.real(self.data[index][slicing]), axis=axes)
                    tmpminPos = minArgPos.flatten()
                    tmp = self.xaxArray[axes][slice(minPos, maxPos)][tmpminPos].reshape(minArgPos.shape)
                    if keepdims:
                        tmpdata += (np.expand_dims(tmp, axes), )
                    else:
                        tmpdata += (tmp, )
                elif which == 6:
                    tmpdata[index] += (np.mean(self.data[index][slicing], axis=axes, keepdims=keepdims), )
        #Remove hyper along this dim if necessary
        tmphyper = copy.deepcopy(self.hyper)
        if keepdims ==  False and axes in self.hyper:
            tmpdata, tmphyper = self.deleteHyper(axes,tmpdata, tmphyper)
            #Correct hyper for missing dim 
            for i in range(len(tmphyper)):
                if tmphyper[i] > axes:
                    tmphyper[i] += -1
        if len(tmpdata[0]) == 1:
            if self.ndim() == 1:
                for index in range(len(tmpdata)):
                    tmpdata[index] = np.array([tmpdata[index][0]])
                newSpec = Spectrum(self.name,
                                   tmpdata,
                                   self.filePath,
                                   copy.deepcopy(self.freq),
                                   copy.deepcopy(self.sw),
                                   copy.deepcopy(self.spec),
                                   copy.deepcopy(self.wholeEcho),
                                   tmphyper,
                                   copy.deepcopy(self.ref),
                                   copy.deepcopy(self.xaxArray),
                                   copy.deepcopy(self.history),
                                   self.msgHandler)
                newSpec.resetXax(axes)
            else:
                tmpXax = copy.deepcopy(self.xaxArray)
                del tmpXax[axes]
                for index in range(len(tmpdata)):
                    tmpdata[index] = tmpdata[index][0]
                newSpec = Spectrum(self.name,
                                   tmpdata,
                                   self.filePath,
                                   copy.deepcopy(np.delete(self.freq, axes)),
                                   copy.deepcopy(np.delete(self.sw, axes)),
                                   copy.deepcopy(np.delete(self.spec, axes)),
                                   copy.deepcopy(np.delete(self.wholeEcho, axes)),
                                   tmphyper,
                                   copy.deepcopy(np.delete(self.ref, axes)),
                                   tmpXax,
                                   copy.deepcopy(self.history),
                                   self.msgHandler)
        else:
            for index in range(len(tmpdata)):
                tmpdata[index] = np.concatenate(tmpdata[index], axis=axes)
            newSpec = Spectrum(self.name,
                               tmpdata,
                               self.filePath,
                               copy.deepcopy(self.freq),
                               copy.deepcopy(self.sw),
                               copy.deepcopy(self.spec),
                               copy.deepcopy(self.wholeEcho),
                               tmphyper,
                               copy.deepcopy(self.ref),
                               copy.deepcopy(self.xaxArray),
                               copy.deepcopy(self.history),
                               self.msgHandler)
            newSpec.resetXax(axes)
        if which == 0:
            newSpec.addHistory("Integrate between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 5:
            newSpec.addHistory("Sum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 1:
            newSpec.addHistory("Maximum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 2:
            newSpec.addHistory("Minimum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 3:
            newSpec.addHistory("Maximum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 4:
            newSpec.addHistory("Minimum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        elif which == 6:
            newSpec.addHistory("Average between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes + 1))
        return newSpec

    def getRegion(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        if self.spec[axes] == 1:
            oldFxax = self.xaxArray[axes][slice(minPos, maxPos)][0]
            self.sw[axes] = self.sw[axes] * (maxPos - minPos) / (1.0 * self.shape()[axes])
        for index in range(len(self.data)):
            self.data[index] = self.data[index][slicing]
        if self.spec[axes] == 1:
            newFxax = np.fft.fftshift(np.fft.fftfreq(self.shape()[axes], 1.0 / self.sw[axes]))[0]
            if self.ref[axes] is None:
                self.ref[axes] = self.freq[axes]
            self.freq[axes] = self.ref[axes] - newFxax + oldFxax
        self.resetXax(axes)
        self.addHistory("Extracted part between " + str(minPos) + " and " + str(maxPos) + " of dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.getRegion(axes, pos1, pos2))

    def getRegionNew(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), ) * (self.data.ndim - 1 - axes)
        tmpsw = copy.deepcopy(self.sw)
        tmpfreq = copy.deepcopy(self.freq)
        tmpref = copy.deepcopy(self.ref)
        if self.spec[axes] == 1:
            oldFxax = self.xaxArray[axes][slice(minPos, maxPos)][0]
            tmpsw[axes] = tmpsw[axes] * (maxPos - minPos) / (1.0 * self.data.shape[axes])
        tmpdata = self.data[slicing]
        if self.spec[axes] == 1:
            newFxax = np.fft.fftshift(np.fft.fftfreq(tmpdata.shape[axes], 1.0 / tmpsw[axes]))[0]
            if tmpref[axes] is None:
                tmpref[axes] = tmpfreq[axes]
            tmpfreq[axes] = tmpref[axes] - newFxax + oldFxax
        newSpec = Spectrum(self.name,
                           tmpdata,
                           self.filePath,
                           tmpfreq,
                           tmpsw,
                           copy.deepcopy(self.spec),
                           copy.deepcopy(self.wholeEcho),
                           tmpref,
                           copy.deepcopy(self.xaxArray),
                           copy.deepcopy(self.history),
                           self.msgHandler)
        newSpec.resetXax(axes)
        newSpec.addHistory("Extracted part between " + str(minPos) + " and " + str(maxPos) + " of dimension " + str(axes + 1))
        return newSpec

    def fiddle(self, refSpec, lb, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        axLen = self.shape()[axes]
        if len(refSpec) != axLen:
            self.dispMsg("Reference FID does not have the correct length")
            return None
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
        for index in range(len(self.data)):
            self.data[index] *= idealFid
        #self.data = np.apply_along_axis(np.multiply, axes, self.data, idealFid)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True)
        self.addHistory("FIDDLE over dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.fiddle(refSpec, lb, axes))

    def diff(self, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        for index in range(len(self.data)):
            self.data[index] = np.diff(self.data[index], axis=axes)
        self.resetXax(axes)
        self.addHistory("Differences over dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.diff(axes))

    def cumsum(self, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        for index in range(len(self.data)):
            self.data[index] = np.cumsum(self.data[index], axis=axes)
        self.addHistory("Cumulative sum over dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.cumsum(axes))

    def flipLR(self, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        slicing = (slice(None), ) * axes + (slice(None, None, -1), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        for index in range(len(self.data)):
            self.data[index] = self.data[index][slicing]
        self.addHistory("Flipped dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.flipLR(axes)

    def hilbert(self, axes):
        import scipy.signal
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        self.data = self.hyperReorder(self.data,axes)
        for index in range(len(self.data)):
            self.data[index] = scipy.signal.hilbert(np.real(self.data[index]), axis=axes)
        self.data = self.hyperReorder(self.data,axes)
        self.addHistory("Hilbert transform on dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.hilbert(axes))

    def autoPhase(self, phaseNum, axes, locList):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if len(locList) != self.ndim()-1:
            self.dispMsg("Data does not have the correct number of dimensions")
            return None
        if np.any(locList >= np.delete(self.shape(), axes)) or np.any(np.array(locList) < 0):
            self.dispMsg("The location array contains invalid indices")
            return None
        tmp = []
        for item in self.data:
            tmp.append(item[tuple(locList[:axes]) + (slice(None), ) + tuple(locList[axes:])])
        if self.spec[axes] == 0:
                tmp = self.fourierLocal(tmp,0, self.axes)
        x = np.fft.fftshift(np.fft.fftfreq(len(tmp[0]), 1.0 / self.sw[axes])) / self.sw[axes]
        tmp = self.hyperReorder(tmp, axes)
        if phaseNum == 0:
            phases = scipy.optimize.minimize(self.ACMEentropy, [0], (tmp, x, False), method='Powell')
            phase0 = phases['x']
            phase1 = 0.0
        elif phaseNum == 1:
            phases = scipy.optimize.minimize(self.ACMEentropy, [0, 0], (tmp, x ), method='Powell')
            phase0 = phases['x'][0]
            phase1 = phases['x'][1]
        tmp = self.hyperReorder(tmp, axes)
        if self.spec == 0:
                tmp = self.fourierLocal(tmp,1, self.axes)
        if self.ref[axes] is None:
            offset = 0
        else:
            offset = self.freq[axes] - self.ref[axes]
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.shape()[axes], 1.0 / self.sw[axes]) + offset) / self.sw[axes] * phase1 * 1j)
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True)
        self.data = self.hyperReorder(self.data, axes)
        for index in range(len(self.data)):
            self.data[index] = self.data[index] * np.exp(phase0 * 1j)
            self.data[index] = np.apply_along_axis(np.multiply, axes, self.data[index], vector)
        self.data = self.hyperReorder(self.data, axes)
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True, inv=True)
        Message = "Autophase: phase0 = " + str(phase0 * 180 / np.pi) + " and phase1 = " + str(phase1 * 180 / np.pi) + " for dimension " + str(axes + 1)
        self.addHistory(Message)
        if self.noUndo:
            return None
        else:
            return lambda self: self.setPhase(-phase0, -phase1, axes)

    def ACMEentropy(self, phaseIn, data, x, phaseAll=True):
        hyperView = 0 #Temp
        phase0 = phaseIn[0]
        if phaseAll:
            phase1 = phaseIn[1]
        else:
            phase1 = 0.0
        L = len(data[hyperView])
        s0 = data[hyperView] * np.exp(1j * (phase0 + phase1 * x))
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

    def fourierLocal(self, fourData, spec, axis, wholeEcho = False):  # fourier the local data for other functions
        ax = len(fourData[0].shape) - 1
        fourData = self.hyperReorder(fourData, axis)
        if spec == 0:
            if not wholeEcho:
                slicing = (slice(None), ) * ax + (0, )
                for index in range(len(fourData)):
                    fourData[index][slicing] = fourData[index][slicing] * 0.5
            for index in range(len(fourData)):
                if ax == 0:
                    fourData[index] = np.fft.fftshift(np.fft.fft(fourData[index]))
                else:
                    fourData[index] = np.fft.fftshift(np.fft.fftn(fourData[index], axes=[axis]), axes=axis)
        else:
            for index in range(len(fourData)):
                if ax == 0:
                    fourData[index] = np.fft.ifft(np.fft.ifftshift(fourData[index]))
                else:
                    fourData[index] = np.fft.ifftn(np.fft.ifftshift(fourData[index], axes=axis), axes=[axis])
            if not wholeEcho:
                slicing = (slice(None), ) * ax + (0, )
                for index in range(len(fourData)):
                    fourData[index][slicing] = fourData[index][slicing] * 2.0
        fourData = self.hyperReorder(fourData, axis)
        return fourData

    def phaseLocal(self, data, sw, offset, phase0, phase1, axis): #Provides a phase function on any data
        tmpdat = self.hyperReorder(data, axis)
        #Data input always as spectrum (calling code should make sure of this)
        if len(data[0].shape) > 1:
            mult = np.repeat([np.exp(np.fft.fftshift(np.fft.fftfreq(len(data[0][0]), 1.0 / sw) + offset) / sw * phase1 * 1j)], len(tmpdat[0]), axis=0)
        else:
            mult = np.exp(np.fft.fftshift(np.fft.fftfreq(len(data[0]), 1.0 / sw) + offset) / sw * phase1 * 1j)
        for index in range(len(tmpdat)):
            tmpdat[index] = tmpdat[index] * mult * np.exp(phase0 * 1j)
        tmpdat = self.hyperReorder(tmpdat, axis)
        return tmpdat

    def setPhase(self, phase0, phase1, axes, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if self.ref[axes] is None:
            offset = 0
        else:
            offset = self.freq[axes] - self.ref[axes]
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.shape()[axes], 1.0 / self.sw[axes]) + offset) / self.sw[axes] * phase1 * 1j)
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True)
        self.data = self.hyperReorder(self.data, axes)
        for index in range(len(self.data)):
            self.data[index][select] = self.data[index][select] * np.exp(phase0 * 1j)
            self.data[index][select] = np.apply_along_axis(np.multiply, axes, self.data[index], vector)[select]
        self.data = self.hyperReorder(self.data, axes)
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True, inv=True)
        Message = "Phasing: phase0 = " + str(phase0 * 180 / np.pi) + " and phase1 = " + str(phase1 * 180 / np.pi) + " for dimension " + str(axes + 1)
        if select != slice(None, None, None):
            Message = Message + " of data[" + str(select) + "]"
        self.addHistory(Message)
        if self.noUndo:
            return None
        else:
            return lambda self: self.setPhase(-phase0, -phase1, axes, select=select)

    def apodize(self, lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, axes, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if shiftingAxes is None:
            shiftingAxes = 0
            shifting = 0
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        axLen = self.shape()[axes]
        t = np.arange(0, axLen) / self.sw[axes]
        if shifting != 0.0:
            for j in range(self.shape()[shiftingAxes]):
                shift1 = shift + shifting * j / self.sw[shiftingAxes]
                x = func.apodize(t, shift1, self.sw[axes], axLen, lor, gauss, cos2, hamming, self.wholeEcho[axes])
                if self.spec[axes] > 0:
                    self.fourier(axes, tmp=True)
                for i in range(self.shape()[axes]):
                    if axes < shiftingAxes:
                        slicing = (slice(None), ) * axes + (i, ) + (slice(None), ) * (shiftingAxes - 1 - axes) + (j, ) + (slice(None), ) * (self.ndim() - 2 - shiftingAxes)
                    else:
                        slicing = (slice(None), ) * shiftingAxes + (j, ) + (slice(None), ) * (axes - 1 - shiftingAxes) + (i, ) + (slice(None), ) * (self.ndim() - 2 - axes)
                    for index in range(len(self.data)): #For all hypercomplex parts
                        self.data[index][slicing] = self.data[index][slicing] * x[i]
                if self.spec[axes] > 0:
                    self.fourier(axes, tmp=True, inv=True)
        else:
            x = func.apodize(t, shift, self.sw[axes], axLen, lor, gauss, cos2, hamming, self.wholeEcho[axes])
            if self.spec[axes] > 0:
                self.fourier(axes, tmp=True)
            for index in range(len(self.data)): #For all hypercomplex parts
                self.data[index][select] = np.apply_along_axis(np.multiply, axes, self.data[index], x)[select]
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
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, axes, select=select))

    def setFreq(self, freq, sw, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
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
        if self.noUndo:
            return None
        else:
            return lambda self: self.setFreq(oldFreq, oldSw, axes)

    def setRef(self, ref, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        oldRef = self.ref[axes]
        if ref is None:
            self.ref[axes] = None
            self.addHistory("Reference frequency set to 'None' for dimension " + str(axes + 1))
        else:
            self.ref[axes] = float(ref)
            self.addHistory("Reference frequency set to " + str(ref * 1e-6) + " MHz for dimension " + str(axes + 1))
        self.resetXax(axes)
        if self.noUndo:
            return None
        else:
            return lambda self: self.setRef(oldRef, axes)

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
        newDat = []
        if len(self.shape()) > 1:
            for index in range(len(self.data)):
                newDat.append(np.apply_along_axis(self.regridFunc, axis, self.data[index],newAxis,self.xaxArray[axis]))
        else:
            for index in range(len(self.data)):
                newDat.append(self.regridFunc(self.data[index], newAxis ,self.xaxArray[axis]))
        self.data = newDat
        self.sw[axis] = newSw
        if self.ref[axis] is None:  # Set new 0 freq to those of the old view, if needed
            self.ref[axis] = self.freq[axis]
        else:
            newFreq += - self.freq[axis] + self.ref[axis]
        self.freq[axis] = newFreq
        self.resetXax(axis)
        self.addHistory("Regrid dimension " + str(axis) + " between " + str(limits[0]) + ' and ' + str(limits[1]) + ' with ' + str(numPoints) + ' points')
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.regrid(limits, numPoints, axis))

    def regridFunc(self, data, newAx, x):
        from scipy import interpolate as intp
        f = intp.interp1d(x, data, fill_value=0, bounds_error=False)
        return f(newAx)

    def setWholeEcho(self, val, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        self.wholeEcho[axes] = val
        self.addHistory("Whole echo set to " + str(val) + " for dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.setWholeEcho(not val, axes)

    def setSize(self, size, pos, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True)
        if size > self.shape()[axes]:
            slicing1 = (slice(None), ) * axes + (slice(None, pos), ) + (slice(None), ) * (self.ndim() - 1 - axes)
            slicing2 = (slice(None), ) * axes + (slice(pos, None), ) + (slice(None), ) * (self.ndim() - 1 - axes)
            for index in range(len(self.data)):
                self.data[index] = np.concatenate((np.pad(self.data[index][slicing1], [(0, 0)] * axes + [(0, size - self.data[index].shape[axes])] + [(0, 0)] * (self.data[index].ndim - axes - 1), 'constant', constant_values=0),
                                        self.data[index][slicing2]), axes)
        else:
            difference = self.shape()[axes] - size
            removeBegin = int(np.floor(difference / 2))
            removeEnd = difference - removeBegin
            if pos < removeBegin:
                slicing = (slice(None), ) * axes + (slice(self.shape()[axes] - size, None), ) + (slice(None), ) * (self.ndim() - 1 - axes)
                for index in range(len(self.data)):
                    self.data[index] = self.data[index][slicing]
            elif self.shape()[axes] - pos < removeEnd:
                slicing = (slice(None), ) * axes + (slice(None, size), ) + (slice(None), ) * (self.ndim() - 1 - axes)
                for index in range(len(self.data)):
                    self.data[index] = self.data[index][slicing]
            else:
                slicing1 = (slice(None), ) * axes + (slice(None, pos - removeBegin), ) + (slice(None), ) * (self.ndim() - 1 - axes)
                slicing2 = (slice(None), ) * axes + (slice(pos + removeEnd, None), ) + (slice(None), ) * (self.ndim() - 1 - axes)
                for index in range(len(self.data)):
                    self.data[index] = np.concatenate((self.data[index][slicing1], self.data[index][slicing2]), axis=axes)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True)
        self.resetXax(axes)
        self.addHistory("Resized dimension " + str(axes + 1) + " to " + str(size) + " points at position " + str(pos))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.setSize(size, pos, axes))

    def setLPSVD(self, nAnalyse, nFreq, nPredict, Direction, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
            returnValue = lambda self: self.restoreData(copyData, lambda self: self.setLPSVD(nAnalyse, nFreq, nPredict, Direction, axes))
        for index in range(len(self.data)):
            self.data[index] = np.apply_along_axis(self.LPSVDfunction, axes, self.data[index], nAnalyse, nFreq, nPredict, Direction)
        self.resetXax(axes)
        self.addHistory("LPSVD ")
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.setLPSVD(nAnalyse, nFreq, nPredict, Direction, axes))

    def LPSVDfunction(self, data, nAnalyse, nFreq, nPredict, Direction):
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

    def changeSpec(self, val, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        oldVal = self.spec[axes]
        self.spec[axes] = val
        self.resetXax(axes)
        if val:
            self.addHistory("Dimension " + str(axes + 1) + " set to FID")
        else:
            self.addHistory("Dimension " + str(axes + 1) + " set to spectrum")
        if self.noUndo:
            return None
        else:
            return lambda self: self.changeSpec(oldVal, axes)

    def swapEcho(self, idx, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        slicing1 = (slice(None), ) * axes + (slice(None, idx), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        slicing2 = (slice(None), ) * axes + (slice(idx, None), ) + (slice(None), ) * (self.ndim() - 1 - axes)
        for index in range(len(self.data)):
            self.data[index] = np.concatenate((self.data[index][slicing2], self.data[index][slicing1]), axes)
        self.wholeEcho[axes] = not self.wholeEcho[axes]
        self.addHistory("Swap echo at position " + str(idx) + " for dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.swapEcho(-idx, axes)

    def shiftData(self, shift, axes, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True)
        mask = np.ones(self.shape()[axes])
        if shift < 0:
            mask[slice(shift, None)] = 0
        else:
            mask[slice(None, shift)] = 0
        for index in range(len(self.data)):
            self.data[index][select] = np.roll(self.data[index], shift, axes)[select]
            self.data[index][select] = np.apply_along_axis(np.multiply, axes, self.data[index], mask)[select]
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True)
        Message = "Shifted " + str(shift) + " points in dimension " + str(axes + 1)
        if select != slice(None, None, None):
            Message = Message + " of data[" + str(select) + "]"
        self.addHistory(Message)
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.shiftData(shift, axes, select=select))

    def hyperReorder(self, data, axis): #A function to reorder the data for a hypercomplex operation
        hyper = [x for x in self.hyper if x == axis]
        if len(hyper) == 0:
            hyper = None
        elif len(hyper) == 1:
            hyper = self.hyper.index(hyper[0])
        else:
            print('error in hyper')
            return
        hyperLen = len(data)
        if hyper == None:
            return data
        else:
            values = np.arange(hyperLen)
            step = 2**(len(self.hyper) - hyper - 1)
            list1 = np.array([],dtype=int)
            list2 = np.array([],dtype=int)
            for index in range(int(hyperLen/step)):
                if index % 2: #if even
                    list2 = np.append(list2,values[0:step])
                else:
                    list1 = np.append(list1,values[0:step])
                values = values[step::]
        for index in range(len(list1)):
            l1 = list1[index]
            l2 = list2[index]
            data[l1], data[l2] = np.real(data[l1]) + 1j*np.real(data[l2]), np.imag(data[l1]) + 1j*np.imag(data[l2])
        return data

    def deleteHyper(self,axis,data,hyper):
        #Deletes hypercomplex data along axis, is any
        #Deletes its entry from self.hyper list.
        if axis in hyper:
            totlen = 2**len(hyper)
            indx = hyper.index(axis)
            step = 2**(len(hyper) - indx - 1)
            boollist = np.array([True,False])
            boollist = np.tile(np.repeat(boollist,totlen / step / 2),step)
            newdat = []
            for i in range(len(boollist)):
                if boollist[i] == True:
                    newdat.append(data[i])
            del hyper[indx]
            data = newdat
        return data, hyper

    def fourier(self, axes, tmp=False, inv=False, reorder = [True,True]):
        axes = self.checkAxes(axes)
        tmpdat = self.data 
        if reorder[0]:
            tmpdat = self.hyperReorder(self.data, axes)
        if axes is None:
            return None
        if np.logical_xor(self.spec[axes], inv) == 0:
            if not self.wholeEcho[axes] and not tmp:
                slicing = (slice(None), ) * axes + (0, ) + (slice(None), ) * (self.ndim() - 1 - axes)
                for index in range(len(tmpdat)):
                    tmpdat[index][slicing] = tmpdat[index][slicing] * 0.5
            for index in range(len(tmpdat)): 
                tmpdat[index] = np.fft.fftshift(np.fft.fftn(tmpdat[index], axes=[axes]), axes=axes)
            if not tmp:
                self.spec[axes] = 1
                self.addHistory("Fourier transform dimension " + str(axes + 1))
        else:
            for index in range(len(tmpdat)): 
                tmpdat[index] = np.fft.ifftn(np.fft.ifftshift(tmpdat[index], axes=axes), axes=[axes])
            if not self.wholeEcho[axes] and not tmp:
                slicing = (slice(None), ) * axes + (0, ) + (slice(None), ) * (self.ndim() - 1 - axes)
                for index in range(len(tmpdat)): 
                    tmpdat[index][slicing] = tmpdat[index][slicing] * 2.0
            if not tmp:
                self.spec[axes] = 0
                self.addHistory("Inverse Fourier transform dimension " + str(axes + 1))
        if reorder[1]:
            self.data = self.hyperReorder(tmpdat, axes)
        self.resetXax(axes)
        if self.noUndo:
            return None
        else:
            return lambda self: self.fourier(axes)

    def realFourier(self, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        if self.spec[axes] == 0:
            if not self.wholeEcho[axes]:
                slicing = (slice(None), ) * axes + (0, ) + (slice(None), ) * (self.ndim() - 1 - axes)
                for index in range(len(self.data)):
                    self.data[index][slicing] = self.data[index][slicing] * 0.5
            for index in range(len(self.data)):
                self.data[index] = np.fft.fftshift(np.fft.fftn(np.real(self.data[index]), axes=[axes]), axes=axes)
            self.spec[axes] = 1
            self.addHistory("Real Fourier transform dimension " + str(axes + 1))
        else:
            for index in range(len(self.data)):
                self.data[index] = np.fft.ifftn(np.fft.ifftshift(np.real(self.data[index]), axes=axes), axes=[axes])
            if not self.wholeEcho[axes]:
                slicing = (slice(None), ) * axes + (0, ) + (slice(None), ) * (self.ndim() - 1 - axes)
                for index in range(len(self.data)):
                    self.data[index][slicing] = self.data[index][slicing] * 2.0
            self.spec[axes] = 0
            self.addHistory("Real inverse Fourier transform dimension " + str(axes + 1))
        self.resetXax(axes)
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.realFourier(axes))

    def fftshift(self, axes, inv=False):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if inv:
            for index in range(len(self.data)):
                self.data[index] = np.fft.ifftshift(self.data[index], axes=[axes])
            self.addHistory("Inverse Fourier shift dimension " + str(axes + 1))
        else:
            for index in range(len(self.data)):
                self.data[index] = np.fft.fftshift(self.data[index], axes=axes)
            self.addHistory("Fourier shift dimension " + str(axes + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.fftshift(axes, not(inv))

    def shear(self, shear, axes, axes2):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        axes2 = self.checkAxes(axes2)
        if axes2 is None:
            return None
        if axes == axes2:
            self.dispMsg('Both shearing axes cannot be equal')
            return None
        if self.ndim() < 2:
            self.dispMsg("The data does not have enough dimensions for a shearing transformation")
            return None
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
            self.data = self.hyperReorder(self.data, axes)
        for index in range(len(self.data)):
            self.data[index] = self.data[index] * shearMatrix.reshape(shape)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True, reorder = [False,True])
        else:
            self.data = self.hyperReorder(self.data, axes)
        self.addHistory("Shearing transform with shearing value " + str(shear) + " over dimensions " + str(axes + 1) + " and " + str(axes2 + 1))
        if self.noUndo:
            return None
        else:
            return lambda self: self.shear(-shear, axes, axes2)

    def reorder(self, pos, newLength, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if newLength is None:
            newLength = max(pos) + 1
        if (max(pos) >= newLength) or (min(pos) < 0):
            self.dispMsg("Reorder: invalid positions")
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        newShape = np.array(self.shape())
        newShape[axes] = newLength
        tmpData = np.zeros(newShape, dtype=complex)
        slicing = (slice(None), ) * axes + (pos, ) + (slice(None), ) * (self.ndim() - 1 - axes)
        for index in range(len(self.data)):
            tmpData = np.zeros(newShape, dtype=complex)
            tmpData[slicing] = self.data[index]
            self.data[index] = tmpData
        self.resetXax(axes)
        self.addHistory("Reorder dimension " + str(axes + 1) + " to obtain a new length of " + str(newLength) + " with positions " + str(pos))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, lambda self: self.reorder(pos, newLength, axes))

    def ffm_1d(self, pos, typeVal, axes):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        # pos contains the values of fixed points which not to be translated to missing points
        posList = np.delete(range(self.shape()[axes]), pos)
        if typeVal == 1:  # type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList / 2), dtype=int)
        if typeVal == 2:  # type is TPPI, for now handle the same as Complex
            pass
        posList = np.unique(posList)
        if axes in self.hyper: #Get good hypercomplex part
            self.data = self.hyperReorder(self.data, axes)[0]
        else:
            self.data = self.data[0]
        tmpData = np.rollaxis(self.data, axes, self.data.ndim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((int(tmpData.size / tmpShape[-1]), tmpShape[-1]))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(ffm, [(i, posList) for i in tmpData])
        pool.close()
        pool.join()
        self.data = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axes)
        #Reconstruct hypercomplex parts
        self.data = self.reconstructHyper(self.data)
        #Transform back to FID
        self.data = self.fourierLocal(self.data, 1, axes)
        self.addHistory("Fast Forward Maximum Entropy reconstruction of dimension " + str(axes + 1) + " at positions " + str(pos))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, None)

    def clean(self, pos, typeVal, axes, gamma, threshold, maxIter):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        # pos contains the values of fixed points which not to be translated to missing points
        posList = np.delete(range(self.shape()[axes]), pos)
        if typeVal == 1:  # type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList / 2), dtype=int)
        if typeVal == 2:  # type is TPPI, for now handle the same as Complex
            pass
        posList = np.unique(posList)
        if axes in self.hyper: #Take correct hypercomplex
            self.data = self.hyperReorder(self.data, axes)[0]
        else:
            self.data = self.data[0]
        tmpData = np.rollaxis(np.fft.fft(self.data, axis=axes), axes, self.data.ndim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((int(tmpData.size / tmpShape[-1]), tmpShape[-1]))
        mask = np.ones(tmpShape[-1]) / float(tmpShape[-1])
        mask[posList] = 0.0
        mask = np.fft.fft(mask)                                                    # abs or real???
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(clean, [(i, mask, gamma, threshold, maxIter) for i in tmpData])
        pool.close()
        pool.join()
        self.data = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axes)
        #Reconstruct hypercomplex parts
        self.data = self.reconstructHyper(self.data)
        #Transform back to FID
        self.data = self.fourierLocal(self.data, 1, axes)
        self.addHistory("CLEAN reconstruction (gamma = " + str(gamma) + " , threshold = " + str(threshold) + " , maxIter = " + str(maxIter) + ") " + 
        "of dimension " + str(axes + 1) + " at positions " + str(pos))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, None)
        
    def ist(self,pos, typeVal, axes, threshold, maxIter,tracelimit)  :
        import scipy.signal
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if not self.noUndo:
            copyData = copy.deepcopy(self)
        # pos contains the values of fixed points which not to be translated to missing points
        if axes in self.hyper:
            self.data = self.hyperReorder(self.data, axes)[0]
        else:
            self.data = self.data[0]
        posList = np.delete(range(self.data.shape[axes]), pos)
        if typeVal == 1:  # type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList / 2), dtype=int)
        elif typeVal == 2:  # type is TPPI, for now handle the same as Complex
            pass
        NDmax = np.max(np.max(np.abs(np.real(np.fft.fft(self.data,axis=axes))))) #Get max of ND matrix
        tmpData = np.rollaxis(self.data, axes, self.data.ndim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((int(tmpData.size / tmpShape[-1]), tmpShape[-1]))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(ist, [(i, posList, threshold, maxIter, tracelimit, NDmax) for i in tmpData])
        pool.close()
        pool.join()
        self.data = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axes)
        #Reconstruct hypercomplex parts
        self.data = self.reconstructHyper(self.data)
        #Transform back to FID
        self.data = self.fourierLocal(self.data, 1, axes)
        self.addHistory("IST reconstruction (threshold = " + str(threshold) + " , maxIter = " + str(maxIter) + " , tracelimit = " + str(tracelimit*100) + ") " + 
        "of dimension " + str(axes + 1) + " at positions " + str(pos))
        if self.noUndo:
            return None
        else:
            return lambda self: self.restoreData(copyData, None)

    def reconstructHyper(self,data):
        #Reconstructs hyper data from R*ndim spectrum
        hyperLen = len(self.hyper)
        totLen = 2**hyperLen
        newData = [ x for x in range(totLen)]
        newData[0] = copy.copy(data)
        if hyperLen != 0: #Construct hyper parts if any
            hilbertBool = np.zeros((totLen,hyperLen)) #Holds which dims need Hilbert transform
            for index in range(len(self.hyper)):
                tmp2 = [0,1]
                step = totLen / (2 ** (index + 1))
                tmp2 = np.tile(np.repeat(tmp2,step),totLen / step / 2)
                hilbertBool[:,index] = tmp2
            for index in range(1,totLen): #For all but the first
                tmp = copy.copy(data) # Get the original data
                for hyper in range(len(hilbertBool[index])):
                    if hilbertBool[index][hyper] == 1.0:
                        tmp = scipy.signal.hilbert(np.real(tmp),axis = self.hyper[hyper])
                        tmp = np.imag(np.conj(tmp))
                newData[index] = tmp
            for index in range(len(newData)):#Do hilbert in the direct dim for all. Axis should not matter
                newData[index] = np.conj(scipy.signal.hilbert(np.real(newData[index]),axis = -1))
        return newData

    def getSlice(self, axes, locList):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        datList = []
        for item in self.data:
            datList.append(item[tuple(locList[:axes]) + (slice(None), ) + tuple(locList[axes:])])
        return copy.deepcopy((datList,
                              self.freq[axes],
                              self.sw[axes],
                              self.spec[axes],
                              self.wholeEcho[axes],
                              self.xaxArray[axes],
                              self.ref[axes]))

    def getBlock(self, axes, axes2, locList, stackBegin=None, stackEnd=None, stackStep=None):
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        axes2 = self.checkAxes(axes2)
        if axes2 is None:
            return None
        stackSlice = reim.floatSlice(stackBegin, stackEnd, stackStep)
        if axes == axes2:
            self.dispMsg("First and second axes are the same")
            return None
        elif axes < axes2:
            datList = []
            for item in self.data:
                datList.append(np.transpose(item[tuple(locList[:axes]) + (slice(None), ) + tuple(locList[axes:axes2 - 1]) + (stackSlice, ) + tuple(locList[axes2 - 1:])]))
            return copy.deepcopy((datList,
                                  self.freq[axes],
                                  self.freq[axes2],
                                  self.sw[axes],
                                  self.sw[axes2],
                                  self.spec[axes],
                                  self.spec[axes2],
                                  self.wholeEcho[axes],
                                  self.wholeEcho[axes2],
                                  self.xaxArray[axes],
                                  self.xaxArray[axes2][stackSlice],
                                  self.ref[axes],
                                  self.ref[axes2]))
        elif axes > axes2:
            datList = []
            for item in self.data:
                datList.append(item[tuple(locList[:axes2]) + (stackSlice, ) + tuple(locList[axes2:axes - 1]) + (slice(None), ) + tuple(locList[axes - 1:])])
            return copy.deepcopy((datList,
                                  self.freq[axes],
                                  self.freq[axes2],
                                  self.sw[axes],
                                  self.sw[axes2],
                                  self.spec[axes],
                                  self.spec[axes2],
                                  self.wholeEcho[axes],
                                  self.wholeEcho[axes2],
                                  self.xaxArray[axes],
                                  self.xaxArray[axes2][stackSlice],
                                  self.ref[axes],
                                  self.ref[axes2]))

    def restoreData(self, copyData, returnValue):  # restore data from an old copy for undo purposes
        if (not self.noUndo) and returnValue is None:
            copyData2 = copy.deepcopy(self)
        self.data = copyData.data
        self.hyper = copyData.hyper
        self.freq = copyData.freq  # array of center frequency (length is dim, MHz)
        self.filePath = copyData.filePath
        self.sw = copyData.sw  # array of sweepwidths
        self.spec = copyData.spec
        self.wholeEcho = copyData.wholeEcho
        self.xaxArray = copyData.xaxArray
        self.ref = copyData.ref
        self.addHistory("Data was restored to a previous state ")
        if self.noUndo:
            return None
        else:
            if returnValue is None:
                return lambda self: self.restoreData(copyData2, None)
            else:
                return returnValue

##################################################################################################
# the class from which the 1d data is displayed, the operations which only edit the content of this class are for previewing


class Current1D(Plot1DFrame):

    X_RESIZE = False
    Y_RESIZE = False

    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        super(Current1D, self).__init__(root, fig, canvas)
        self.xax = None  # x-axis
        self.data = data  # the actual spectrum instance
        self.freq = None  # frequency of the slice
        self.freq2 = None  # frequency of the slice
        self.sw = None  # x-data display
        self.data1D = None  # the data1D
        self.spec = None  # boolean where False=time domain and True=spectral domain
        self.wholeEcho = None
        self.ref = None  # reference frequency
        if duplicateCurrent is None:
            self.ppm = self.root.father.defaultPPM             # display frequency as ppm
            self.ppm2 = self.root.father.defaultPPM             # display frequency as ppm
            self.axes = len(self.data.shape()) - 1
            self.axes2 = 0
            self.resetLocList()
            self.viewSettings = {"plotType": 0,
                                 "axType": self.root.father.defaultUnits,
                                 "axType2": self.root.father.defaultUnits,
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
                                 "multiValue": 1.5,
                                 "projTop": 0,
                                 "projRight": 0,
                                 "projLimitsBool": False,
                                 "projLimits": [None, None, None, None],
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
            self.ppm = duplicateCurrent.ppm
            self.ppm2 = duplicateCurrent.ppm2
            self.axes = duplicateCurrent.axes
            self.axes2 = duplicateCurrent.axes2
            if isinstance(self, (CurrentStacked, CurrentArrayed, CurrentContour)):
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

    def dispMsg(self):
        self.data1D.dispMsg()

    def startUp(self, xReset=True, yReset=True):
        self.plotReset(xReset, yReset)  # reset the axes limits
        self.showFid()  # plot the data

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
            updateVar = self.data.getSlice(self.axes, self.locList)
            if updateVar is None:
                self.root.rescue()
        except Exception:
            self.resetLocList()
            updateVar = self.data.getSlice(self.axes, self.locList)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.sw = updateVar[2]
        self.spec = updateVar[3]
        self.wholeEcho = updateVar[4]
        self.xax = updateVar[5]
        self.ref = updateVar[6]
        if self.ref is None:
            self.ref = self.freq
        self.single = self.data1D[0].shape[-1] == 1
        return True

    def setSlice(self, axes, locList):  # change the slice
        axesSame = True
        if self.axes != axes:
            axesSame = False
            self.axes = axes
        self.locList = locList
        self.upd()
        if not axesSame:
            self.plotReset()
        self.showFid()

    def resetLocList(self):
        self.locList = [0] * (len(self.data.shape()) - 1)

    def getSelect(self):
        tmp = list(self.locList)
        if len(self.data1D[0].shape) > 1:
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
        
    def setPhaseInter(self, phase0in, phase1in):  # interactive changing the phase without editing the actual data
        phase0 = float(phase0in)
        phase1 = float(phase1in)
        self.upd()
        if self.spec == 0:
            tmpdata = self.fourierLocal(self.data1D, 0, self.axes)
        else:
            tmpdata = self.data1D
        if self.ref is None:
            offset = 0
        else:
            offset = +self.freq - self.ref
        tmpdata = self.data.phaseLocal(tmpdata,self.sw,offset,phase0,phase1,self.axes)
        if self.spec == 0:
            tmpdata = self.fourierLocal(tmpdata, 1, self.axes)
        self.data1D = tmpdata
        self.showFid()

    def applyPhase(self, phase0, phase1, select=False):  # apply the phase to the actual data
        phase0 = float(phase0)
        phase1 = float(phase1)
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.setPhase(phase0, phase1, self.axes, selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['phase', (phase0, phase1, self.axes - self.data.ndim(), str(selectSlice))])
        return returnValue

    def fourier(self):  # fourier the actual data and replot
        returnValue = self.data.fourier(self.axes)
        self.upd()
        if isinstance(self, (CurrentStacked, CurrentArrayed)):
            self.resetSpacing()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['fourier', (self.axes - self.data.ndim(), )])
        return returnValue

    def realFourier(self):  # fourier the real data and replot
        returnValue = self.data.realFourier(self.axes)
        self.upd()
        if isinstance(self, (CurrentStacked, CurrentArrayed)):
            self.resetSpacing()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['realFourier', (self.axes - self.data.ndim(), )])
        return returnValue

    def fftshift(self, inv=False):  # fftshift the actual data and replot
        returnValue = self.data.fftshift(self.axes, inv)
        self.upd()
        self.showFid()
        self.root.addMacro(['fftshift', (self.axes - self.data.ndim(), inv)])
        return returnValue

    def fourierLocal(self, fourData, spec, axis):  # fourier the local data for other functions
        #Now links to data structure function
        return self.data.fourierLocal(fourData, spec, axis, self.wholeEcho)

    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=None):  # display the 1D data including the apodization function
        hyperView = 0
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                self.dispMsg('shiftingAxes cannot be equal to axes')
                return
            elif shiftingAxes < self.axes:
                shift += shifting * self.locList[shiftingAxes] / self.data.sw[shiftingAxes]
            else:
                shift += shifting * self.locList[shiftingAxes - 1] / self.data.sw[shiftingAxes]
        length = len(self.data1D[0])
        t = np.arange(0, length) / (self.sw)
        x = func.apodize(t, shift, self.sw, length, lor, gauss, cos2, hamming, self.wholeEcho)
        self.ax.cla()
        y = copy.copy(self.data1D)
        if self.spec == 1:
            y = self.fourierLocal(y, 1, self.axes)
            for index in range(len(y)):
                y[index] = y[index] * x
            y = self.fourierLocal(y, 0, self.axes)
        else:
            for index in range(len(y)):
                y[index] = y[index] * x
        if self.spec == 0:
            if self.viewSettings["plotType"] == 0:
                scale = np.max(np.real(self.data1D[hyperView]))
            elif self.viewSettings["plotType"] == 1:
                scale = np.max(np.imag(self.data1D[hyperView]))
            elif self.viewSettings["plotType"] == 2:
                scale = np.max(np.max(np.real(self.data1D[hyperView])), np.max(np.imag(self.data1D[hyperView])))
            elif self.viewSettings["plotType"] == 3:
                scale = np.max(np.abs(self.data1D[hyperView]))
            self.showFid(y[hyperView], [t], [x * scale], ['g'], old=True)
        else:
            self.showFid(y[hyperView])

    def applyApod(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=0, select=False):  # apply the apodization to the actual data
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, self.axes, selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['apodize', (lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, self.axes - self.data.ndim(), str(selectSlice))])
        return returnValue

    def setFreq(self, freq, sw):  # set the frequency of the actual data
        returnValue = self.data.setFreq(freq, sw, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['freq', (freq, sw, self.axes - self.data.ndim())])
        return returnValue

    def setRef(self, ref):  # set the frequency of the actual data
        oldref = self.ref
        if oldref is None:
            oldref = self.freq
        returnValue = self.data.setRef(ref, self.axes)
        if ref is None:
            ref = self.freq
        val = self.viewSettings["axType"]
        if self.spec == 1:
            if self.ppm:
                self.xminlim = (self.xminlim * oldref * 10**-6 + oldref - ref) / (ref * 10**-6)
                self.xmaxlim = (self.xmaxlim * oldref * 10**-6 + oldref - ref) / (ref * 10**-6)
            else:
                self.xminlim = self.xminlim + (oldref - ref) / 10**(val * 3)  # set new limits, and scale for axis type
                self.xmaxlim = self.xmaxlim + (oldref - ref) / 10**(val * 3)
        self.upd()
        self.showFid()
        self.root.addMacro(['ref', (ref, self.axes - self.data.ndim())])
        return returnValue

    def regrid(self, limits, numPoints):
        ax = self.axes
        returnValue = self.data.regrid(limits, numPoints, ax)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['regrid', (limits, numPoints,ax - self.data.ndim())])
        return returnValue

    def SN(self, minNoise, maxNoise, minPeak, maxPeak):
        hyperView = 0
        minN = min(minNoise, maxNoise)
        maxN = max(minNoise, maxNoise)
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        if len(self.data1D[0].shape) > 1:
            tmpData = self.data1D[hyperView][0]
        else:
            tmpData = self.data1D[hyperView]
        if (self.viewSettings["plotType"] == 0):
            tmpData = np.real(tmpData)
        elif(self.viewSettings["plotType"] == 1):
            tmpData = np.imag(tmpData)
        elif(self.viewSettings["plotType"] == 2):
            tmpData = np.real(tmpData)
        elif(self.viewSettings["plotType"] == 3):
            tmpData = np.abs(tmpData)
        return (np.max(tmpData[minP:maxP]) / (np.std(tmpData[minN:maxN])))

    def fwhm(self, minPeak, maxPeak, unitType=None):
        hyperView = 0
        from scipy.interpolate import UnivariateSpline
        if unitType is None:
            axType = self.viewSettings["axType"]
            ppm = self.ppm
        else:
            axType = unitType
            if unitType == 3:  # ppm
                ppm = 1
            else:
                ppm = 0
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        if len(self.data1D[0].shape) > 1:
            tmpData = self.data1D[hyperView][0]
        else:
            tmpData = self.data1D[hyperView]
        if (self.viewSettings["plotType"] == 0):
            tmpData = np.real(tmpData)
        elif(self.viewSettings["plotType"] == 1):
            tmpData = np.imag(tmpData)
        elif(self.viewSettings["plotType"] == 2):
            tmpData = np.real(tmpData)
        elif(self.viewSettings["plotType"] == 3):
            tmpData = np.abs(tmpData)
        maxPos = np.argmax(tmpData[minP:maxP])
        x = self.xax * self.getAxMult(self.spec, axType, ppm, self.freq, self.ref)
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
        hyperView = 0 
        if unitType is None:
            axType = self.axType
            ppm = self.ppm
        else:
            axType = unitType
            if unitType == 3:  # ppm
                ppm = 1
            else:
                ppm = 0
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        if len(self.data1D[0].shape) > 1:
            tmpData = self.data1D[hyperView][0]
            if len(self.xax.shape) > 1:
                tmpAxis = self.xax[0]
            else:
                tmpAxis = self.xax
        else:
            tmpData = self.data1D[hyperView]
            tmpAxis = self.xax
        if (self.viewSettings["plotType"] == 0):
            tmpData = np.real(tmpData)
        elif(self.viewSettings["plotType"] == 1):
            tmpData = np.imag(tmpData)
        elif(self.viewSettings["plotType"] == 2):
            tmpData = np.real(tmpData)
        elif(self.viewSettings["plotType"] == 3):
            tmpData = np.abs(tmpData)
        tmpAxis = tmpAxis[minP:maxP] * self.getAxMult(self.spec, self.viewSettings["axType"], ppm, self.freq, self.ref)
        tmpData = tmpData[minP:maxP]
        # COM = 1/M *sum(m_i * r_i)
        CentreOM = 1.0 / np.sum(tmpData) * np.sum(tmpData * tmpAxis)
        return CentreOM

    def setSizePreview(self, size, pos):  # set size only on local data
        hyperView = 0
        if len(self.data1D[0].shape) > 1:
            length = len(self.data1D[0][0])
        else:
            length = len(self.data1D[0])
        axes = len(self.data1D[0].shape) - 1
        if self.spec == 1:
            tmpdata = self.fourierLocal(copy.copy(self.data1D), 1, axes)
        else:
            tmpdata = copy.copy(self.data1D)
        if size > length:
            slicing1 = (slice(None), ) * axes + (slice(None, pos), ) + (slice(None), ) * (tmpdata[0].ndim - 1 - axes)
            slicing2 = (slice(None), ) * axes + (slice(pos, None), ) + (slice(None), ) * (tmpdata[0].ndim - 1 - axes)
            for index in range(len(self.data1D)):
                tmpdata[index] = np.concatenate((np.pad(tmpdata[index][slicing1], [(0, 0)] * axes + [(0, size - self.data1D[index].shape[axes])] + [(0, 0)] * (self.data1D[index].ndim - axes - 1),
                    'constant', constant_values=0), tmpdata[index][slicing2]), axes)
        else:
            difference = tmpdata[0].shape[axes] - size
            removeBegin = int(np.floor(difference / 2))
            removeEnd = difference - removeBegin
            if pos < removeBegin:
                slicing = (slice(None), ) * axes + (slice(tmpdata[0].shape[axes] - size, None), ) + (slice(None), ) * (tmpdata[0].ndim - 1 - axes)
                for index in range(len(self.data1D)):
                    tmpdata[index] = tmpdata[index][slicing]
            elif tmpdata[0].shape[axes] - pos < removeEnd:
                slicing = (slice(None), ) * axes + (slice(None, size), ) + (slice(None), ) * (tmpdata[0].ndim - 1 - axes)
                for index in range(len(self.data1D)):
                    tmpdata[index] = tmpdata[index][slicing]
            else:
                slicing1 = (slice(None), ) * axes + (slice(None, pos - removeBegin), ) + (slice(None), ) * (tmpdata[0].ndim - 1 - axes)
                slicing2 = (slice(None), ) * axes + (slice(pos + removeEnd, None), ) + (slice(None), ) * (tmpdata[0].ndim - 1 - axes)
                for index in range(len(self.data1D)):
                    tmpdata[index] = np.concatenate((tmpdata[index][slicing1], tmpdata[index][slicing2]), axis=axes)
        if self.spec == 1:
            self.data1D = self.fourierLocal(tmpdata, 0, axes)
        else:
            self.data1D = tmpdata
        if len(self.data1D[0].shape) > 1:
            length = len(self.data1D[0][0])
        else:
            length = len(self.data1D[0])
        if self.spec == 0:
            self.xax = np.arange(length) / self.sw
        elif self.spec == 1:
            self.xax = np.fft.fftshift(np.fft.fftfreq(length, 1.0 / self.sw)) + self.freq - self.ref
        if not self.spec:
            self.plotReset(True, False)
        self.showFid()
        self.upd()

    def applySize(self, size, pos):  # set size to the actual data
        if self.data.noUndo:
            self.data.setSize(size, pos, self.axes)
            returnValue = None
        else:
            returnValue = self.data.setSize(size, pos, self.axes)
        self.upd()
        if not self.spec:
            self.plotReset(True, False)
        self.showFid()
        self.root.addMacro(['size', (size, pos, self.axes - self.data.ndim())])
        return returnValue

    def applyLPSVD(self, nAnalyse, nFreq, nPredict, Direction):
        returnValue = self.data.setLPSVD(nAnalyse, nFreq, nPredict, Direction, self.axes)
        self.upd()
        self.showFid()
        return returnValue

    def changeSpec(self, val):  # change from time to freq domain of the actual data
        returnValue = self.data.changeSpec(val, self.axes)
        self.upd()
        if isinstance(self, CurrentArrayed):
            self.resetSpacing()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['spec', (val, self.axes - self.data.ndim())])
        return returnValue

    def applySwapEcho(self, idx):
        returnValue = self.data.swapEcho(idx, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['swapecho', (idx, self.axes - self.data.ndim())])
        return returnValue

    def setSwapEchoPreview(self, idx):
        if len(self.data1D[0].shape) > 1:
            for index in range(len(self.data1D)):
                self.data1D[index] = np.concatenate((self.data1D[index][:, idx:], self.data1D[index][:, :idx]), axis=1)
        else:
            for index in range(len(self.data1D)):
                self.data1D[index] = np.concatenate((self.data1D[index][idx:], self.data1D[index][:idx]))
        self.showFid()
        self.upd()

    def setWholeEcho(self, value):
        if value == 0:
            returnValue = self.data.setWholeEcho(False, self.axes)
            self.wholeEcho = False
            self.root.addMacro(['wholeEcho', (False, self.axes - self.data.ndim())])
        else:
            returnValue = self.data.setWholeEcho(True, self.axes)
            self.wholeEcho = True
            self.root.addMacro(['wholeEcho', (True, self.axes - self.data.ndim())])
        return returnValue

    def applyShift(self, shift, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.shiftData(shift, self.axes, selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['shift', (shift, self.axes - self.data.ndim(), str(selectSlice))])
        return returnValue

    def setShiftPreview(self, shift):
        hyperView = 0
        tmpData = self.data1D[hyperView]
        dim = len(self.data1D[0].shape)
        if self.spec > 0:
            tmpData = self.fourierLocal(tmpData, 1, self.axes)
        tmpData = np.roll(tmpData, shift)
        if shift < 0:
            tmpData[(slice(None), ) * (dim - 1) + (slice(shift, None), )] = tmpData[(slice(None), ) * (dim - 1) + (slice(shift, None), )] * 0
        else:
            tmpData[(slice(None), ) * (dim - 1) + (slice(None, shift), )] = tmpData[(slice(None), ) * (dim - 1) + (slice(None, shift), )] * 0
        if self.spec > 0:
            tmpData = self.fourierLocal(tmpData, 0, self.axes)
        self.showFid(tmpData)

    def getdcOffset(self, pos1, pos2):
        hyperView = 0
        minPos = int(min(pos1, pos2))
        maxPos = int(max(pos1, pos2))
        if minPos != maxPos:
            return np.mean(self.data1D[hyperView][(len(self.data1D[hyperView].shape) - 1) * (slice(None), ) + (slice(minPos, maxPos), )])
        else:
            return 0

    def dcOffset(self, offset):
        hyperView = 0
        self.showFid(self.data1D[hyperView] - offset)

    def baselinePolyFit(self, x, data, bArray, degree):
        import numpy.polynomial.polynomial as poly
        polyCoeff = poly.polyfit(x[bArray], data[bArray], degree)
        return poly.polyval(x, polyCoeff)

    def applyBaselineAll(self, degree, removeList, select=False, invert=False):
        hyperView = 0
        tmpAx = np.arange(self.data1D[0].shape[-1])
        bArray = np.array([True] * self.data1D[0].shape[-1])
        for i in range(int(np.floor(len(removeList) / 2.0))):
            minVal = min(removeList[2 * i], removeList[2 * i + 1])
            maxVal = max(removeList[2 * i], removeList[2 * i + 1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        if invert:
            bArray = np.logical_not(bArray)
        try:
            y = np.apply_along_axis(lambda data: self.baselinePolyFit(self.xax, data, bArray, degree), self.axes, self.data.data[hyperView])
            if (self.viewSettings["plotType"] == 0):
                y = np.real(y)
            elif(self.viewSettings["plotType"] == 1):
                y = np.imag(y)
            elif(self.viewSettings["plotType"] == 2):
                y = np.real(y)
            elif(self.viewSettings["plotType"] == 3):
                y = np.abs(y)
            returnValue = self.data.subtract(y,singleHyper = True)
            self.root.addMacro(['subtract', (y.tolist(), None, slice(None), True)])
        except Exception:
            return None
        return returnValue

    def applyBaseline(self, degree, removeList, select=False, invert=False):
        hyperView = 0
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        if len(self.data1D[0].shape) > 1:
            tmpData = self.data1D[hyperView][0]
        else:
            tmpData = self.data1D[hyperView]
        tmpAx = np.arange(self.data1D[0].shape[-1])
        bArray = np.array([True] * self.data1D[0].shape[-1])
        for i in range(int(np.floor(len(removeList) / 2.0))):
            minVal = min(removeList[2 * i], removeList[2 * i + 1])
            maxVal = max(removeList[2 * i], removeList[2 * i + 1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        if invert:
            bArray = np.logical_not(bArray)
        try:
            print('1')
            y = self.baselinePolyFit(self.xax, tmpData, bArray, degree)
            print('2')
            if (self.viewSettings["plotType"] == 0):
                y = np.real(y)
            elif(self.viewSettings["plotType"] == 1):
                y = np.imag(y)
            elif(self.viewSettings["plotType"] == 2):
                y = np.real(y)
            elif(self.viewSettings["plotType"] == 3):
                y = np.abs(y)
            self.root.addMacro(['baselineCorrection', (list(np.real(y)), self.axes - self.data.ndim(), list(np.imag(y)), str(selectSlice))])
            return self.data.baselineCorrection(y, self.axes, select=selectSlice)
        except Exception:
            return None

    def previewBaseline(self, degree, removeList, invert=False):
        hyperView = 0
        if len(self.data1D[0].shape) > 1:
            tmpData = self.data1D[hyperView][0]
        else:
            tmpData = self.data1D[hyperView]
        tmpAx = np.arange(self.data1D[0].shape[-1])
        bArray = np.array([True] * self.data1D[0].shape[-1])
        for i in range(int(np.floor(len(removeList) / 2.0))):
            minVal = min(removeList[2 * i], removeList[2 * i + 1])
            maxVal = max(removeList[2 * i], removeList[2 * i + 1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        if invert:
            bArray = np.logical_not(bArray)
        check = True
        try:
            y = self.baselinePolyFit(self.xax, tmpData, bArray, degree)
            if (self.viewSettings["plotType"] == 0):
                y = np.real(y)
            elif (self.viewSettings["plotType"] == 1):
                y = np.imag(y)
            elif (self.viewSettings["plotType"] == 2):
                y = np.real(y)
            elif (self.viewSettings["plotType"] == 3):
                y = np.abs(y)
        except Exception:
            check = False
        self.resetPreviewRemoveList()
        if check:
            if len(self.data1D[0].shape) > 1:
                self.showFid(self.data1D[hyperView], [self.xax], [y] * self.data1D[0].shape[0], ['g'])
            else:
                self.showFid(self.data1D[hyperView], [self.xax], [y], ['g'])
        else:
            self.showFid()
        self.previewRemoveList(removeList, invert)
        return check

    def previewRemoveList(self, removeList, invert=False):
        axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        self.resetPreviewRemoveList()
        lineColor = 'r'
        if invert:
            lineColor = 'w'
            self.removeListLines.append(self.ax.axvspan(self.xax[0] * axMult, self.xax[-1] * axMult, color='r'))
        for i in range(int(np.floor(len(removeList) / 2.0))):
            self.removeListLines.append(self.ax.axvspan(self.xax[removeList[2 * i]] * axMult, self.xax[removeList[2 * i + 1]] * axMult, color=lineColor))
        if len(removeList) % 2:
            self.removeListLines.append(self.ax.axvline(self.xax[removeList[-1]] * axMult, c=lineColor, linestyle='--'))
        self.canvas.draw()

    def resetPreviewRemoveList(self):
        if hasattr(self, 'removeListLines'):
            for i in self.removeListLines:
                i.remove()
        self.removeListLines = []

    def states(self):
        returnValue = self.data.states(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['states', (self.axes - self.data.ndim(), )])
        return returnValue

    def statesTPPI(self):
        returnValue = self.data.statesTPPI(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['statesTPPI', (self.axes - self.data.ndim(), )])
        return returnValue

    def echoAntiEcho(self):
        returnValue = self.data.echoAntiEcho(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['echoAntiEcho', (self.axes - self.data.ndim(), )])
        return returnValue

    def integrate(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 0)
        else:
            self.root.addMacro(['integrate', (pos1.tolist(), pos2.tolist(), self.axes - self.data.ndim(), )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 0)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def sum(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 5)
        else:
            self.root.addMacro(['sum', (pos1.tolist(), pos2.tolist(), self.axes - self.data.ndim(), )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 5)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def maxMatrix(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 1)
        else:
            self.root.addMacro(['max', (pos1.tolist(), pos2.tolist(), self.axes - self.data.ndim(), )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 1)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def minMatrix(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 2)
        else:
            self.root.addMacro(['min', (pos1.tolist(), pos2.tolist(), self.axes - self.data.ndim(), )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 2)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def argmaxMatrix(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 3)
        else:
            self.root.addMacro(['argmax', (pos1.tolist(), pos2.tolist(), self.axes - self.data.ndim(), )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 3)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def argminMatrix(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 4)
        else:
            self.root.addMacro(['argmin', (pos1.tolist(), pos2.tolist(), self.axes - self.data.ndim(), )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 4)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def average(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 6)
        else:
            self.root.addMacro(['average', (pos1.tolist(), pos2.tolist(), self.axes - self.data.ndim(), )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 6)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def flipLR(self):
        returnValue = self.data.flipLR(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['fliplr', (self.axes - self.data.ndim(), )])
        return returnValue

    def concatenate(self, axes):
        returnValue = self.data.concatenate(axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['concatenate', (axes - self.data.ndim() - 1, )])
        return returnValue

    def split(self, sections):
        returnValue = self.data.split(sections, self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['split', (sections, self.axes - self.data.ndim() + 1)])
        return returnValue

    def diff(self):
        returnValue = self.data.diff(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['diff', (self.axes - self.data.ndim(), )])
        return returnValue

    def cumsum(self):
        returnValue = self.data.cumsum(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['cumsum', (self.axes - self.data.ndim(), )])
        return returnValue

    def insert(self, data, pos):
        returnValue = self.data.insert(data, pos, self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['insert', (np.real(data).tolist(), pos, self.axes - self.data.ndim(), np.imag(data).tolist())])
        return returnValue

    def delete(self, pos):
        returnValue = self.data.remove(pos, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['delete', (pos, self.axes - self.data.ndim())])
        return returnValue

    def deletePreview(self, pos):
        for index in range(len(self.data1D)):
            self.data1D[index] = np.delete(self.data1D[index], pos, axis=len(self.data1D[0].shape) - 1)
        self.xax = np.delete(self.xax, pos)
        if (np.array(self.data1D[0].shape) != 0).all():
            self.showFid()
        self.upd()

    def add(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.add(data, select=selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['add', (np.real(data).tolist(), np.imag(data).tolist(), str(selectSlice))])
        return returnValue

    def subtract(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.subtract(data, select=selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['subtract', (np.real(data).tolist(), np.imag(data).tolist(), str(selectSlice))])
        return returnValue

    def multiplySpec(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.multiplySpec(data, select=selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['multiplySpec', (np.real(data).tolist(), np.imag(data).tolist(), str(selectSlice))])
        return returnValue

    def divideSpec(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.divideSpec(data, select=selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['divideSpec', (np.real(data).tolist(), np.imag(data).tolist(), str(selectSlice))])
        return returnValue

    def multiply(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.multiply(data, self.axes, select=selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['multiply', (np.real(data).tolist(), self.axes - self.data.ndim(), np.imag(data).tolist(), str(selectSlice))])
        return returnValue

    def multiplyPreview(self, data):
        self.upd()
        try:
            for index in range(len(self.data1D)):
                self.data1D[index] = self.data1D[index] * data
            self.showFid()
            return True
        except ValueError as error:
            return error

    def subtractAvg(self, pos1, pos2):
        returnValue = self.data.subtractAvg(pos1, pos2, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['subtractAvg', (pos1, pos2, self.axes - self.data.ndim())])
        return returnValue

    def subtractAvgPreview(self, pos1, pos2):
        self.upd()
        axes = self.data1D[0].ndim - 1
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), ) * (self.data1D[0].ndim - 1 - axes)
        for index in range(len(self.data1D)):
            self.data1D[index] -= np.mean(self.data1D[index][slicing], axis=-1, keepdims=True)
        self.showFid()

    def getRegion(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.getRegionNew(pos1, pos2, self.axes)
        else:
            returnValue = self.data.getRegion(pos1, pos2, self.axes)
            self.upd()
            self.plotReset()
            self.showFid()
            self.root.addMacro(['extract', (pos1, pos2, self.axes - self.data.ndim())])
            return returnValue

    def fiddle(self, pos1, pos2, lb):
        hyperView = 0
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        refSpec = np.zeros(self.data1D[0].shape[-1])
        if len(self.data1D[0].shape) > 1:
            refSpec[minPos:maxPos] = np.real(self.data1D[hyperView][0][minPos:maxPos])
        else:
            refSpec[minPos:maxPos] = np.real(self.data1D[hyperView][minPos:maxPos])
        returnValue = self.data.fiddle(refSpec, lb, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['FIDDLE', (refSpec, lb, self.axes - self.data.ndim())])
        return returnValue

    def shearing(self, shear, axes, axes2):
        returnValue = self.data.shear(shear, axes, axes2)
        self.upd()
        self.showFid()
        self.root.addMacro(['shear', (shear, axes - self.data.ndim(), axes2 - self.data.ndim())])
        return returnValue

    def reorder(self, pos, newLength):
        returnValue = self.data.reorder(pos, newLength, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['reorder', (pos, newLength, self.axes - self.data.ndim())])
        return returnValue

    def ffm(self, posList, typeVal):
        try:
            returnValue = self.data.ffm_1d(posList, typeVal, self.axes)
            self.upd()
            self.showFid()
            self.root.addMacro(['ffm', (posList, typeVal, self.axes - self.data.ndim())])
        except:
            returnValue = None
        return returnValue

    def clean(self, posList, typeVal, gamma, threshold, maxIter):
        try:
            returnValue = self.data.clean(posList, typeVal, self.axes, gamma, threshold, maxIter)
            self.upd()
            self.showFid()
            self.root.addMacro(['clean', (posList, typeVal, self.axes - self.data.ndim(), gamma, threshold, maxIter)])
        except:
            returnValue = None
        return returnValue

    def ist(self, posList, typeVal, threshold, maxIter, tracelimit):
        try:
            returnValue = self.data.ist(posList, typeVal, self.axes, threshold, maxIter, tracelimit)
            self.upd()
            self.showFid()
            self.root.addMacro(['ist', (posList, typeVal, self.axes - self.data.ndim(), threshold, maxIter, tracelimit)])
        except Exception:
            returnValue = None
        return returnValue

    def autoPhase(self, phaseNum):
        self.upd()
        hyperView = 0
        if len(self.data1D[0].shape) > 1:
            tmp = []
            for item in self.data1D:
                tmp.append(item[0])
        else:
            tmp = self.data1D
        if self.spec == 0:
            tmp = self.fourierLocal(tmp,0, self.axes)
        x = np.fft.fftshift(np.fft.fftfreq(len(tmp[0]), 1.0 / self.sw)) / self.sw
        tmp = self.data.hyperReorder(tmp, self.axes)
        if phaseNum == 0:
            phases = scipy.optimize.minimize(self.data.ACMEentropy, [0], (tmp, x, False), method='Powell')
            phases = [phases['x']]
        elif phaseNum == 1:
            phases = scipy.optimize.minimize(self.data.ACMEentropy, [0, 0], (tmp, x ), method='Powell')
            phases = phases['x']
        tmp = self.data.hyperReorder(tmp, self.axes)
        if self.spec == 0:
                tmp = self.fourierLocal(tmp,1, self.axes)
        return phases

    def directAutoPhase(self, phaseNum):
        tmpLocList = self.locList
        if len(self.data1D[0].shape) > 1:
            if hasattr(self, 'stackBegin'):
                val = self.viewSettings["stackBegin"]
            else:
                val = 0
            if self.axes > self.axes2:
                tmpLocList = np.insert(tmpLocList, self.axes2, val)
            else:
                tmpLocList = np.insert(tmpLocList, self.axes2 - 1, val)
        returnValue = self.data.autoPhase(phaseNum, self.axes, tmpLocList)
        self.upd()
        self.showFid()
        self.root.addMacro(['autoPhase', (phaseNum, self.axes - self.data.ndim(), tmpLocList)])
        return returnValue

    def setXaxPreview(self, xax):
        self.xax = xax
        self.plotReset()
        self.showFid()
        self.upd()

    def setXax(self, xax):
        returnVal = self.data.setXax(xax, self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['setxax', (xax, self.axes - self.data.ndim())])
        return returnVal

    def setAxType(self, val):
        oldAxMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        if val == 3:
            self.ppm = True
        else:
            self.ppm = False
            self.viewSettings["axType"] = val
        newAxMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        self.xminlim = self.xminlim * newAxMult / oldAxMult
        self.xmaxlim = self.xmaxlim * newAxMult / oldAxMult
        self.showFid()

    def hilbert(self):
        returnValue = self.data.hilbert(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['hilbert', (self.axes - self.data.ndim(), )])
        return returnValue

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
        absVal = np.max(np.abs(self.data.data))
        if absVal == 0.0:
            return 1
        else:
            return int(np.floor(np.log10(absVal)))

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False, output=None):  # display the 1D data
        hyperView = 0
        self.peakPickReset()
        if tmpdata is None:
            tmpdata = self.data1D[0]
        self.ax.cla()
        axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        self.line_xdata = self.xax * axMult
        if self.single:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = ''
            linestyle = '-'
        if old:
            if (self.viewSettings["plotType"] == 0):
                oldData = np.real(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 1):
                oldData = np.imag(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 2):
                oldData = np.real(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 3):
                oldData = np.abs(self.data1D[hyperView])
            self.ax.plot(self.line_xdata, oldData, marker=marker, linestyle=linestyle, c='k', alpha=0.2, linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if (extraX is not None):
            for num in range(len(extraX)):
                self.ax.plot(extraX[num] * axMult, extraY[num], marker=marker, linestyle=linestyle, c=extraColor[num], linewidth=self.viewSettings["linewidth"], picker=True)
        if (self.viewSettings["plotType"] == 0):
            self.line_ydata = np.real(tmpdata)
        elif(self.viewSettings["plotType"] == 1):
            self.line_ydata = np.imag(tmpdata)
        elif(self.viewSettings["plotType"] == 2):
            self.line_ydata = np.real(tmpdata)
            self.ax.plot(self.line_xdata, np.imag(tmpdata), marker=marker, linestyle=linestyle, c='r', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
        elif(self.viewSettings["plotType"] == 3):
            self.line_ydata = np.abs(tmpdata)
        self.ax.plot(self.line_xdata, self.line_ydata, marker=marker, linestyle=linestyle, c=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
        self.ax.set_xlabel(self.getLabel(self.spec, self.viewSettings["axType"], self.ppm))
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        if self.spec > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        if output is not None:
            self.canvas.print_figure(output)
        self.canvas.draw()

    def plotReset(self, xReset=True, yReset=True):  # set the plot limits to min and max values
        showDat = self.data1D[0]
        if self.viewSettings["plotType"] == 0:
            miny = min(np.real(showDat))
            maxy = max(np.real(showDat))
        elif self.viewSettings["plotType"] == 1:
            miny = min(np.imag(showDat))
            maxy = max(np.imag(showDat))
        elif self.viewSettings["plotType"] == 2:
            miny = min(min(np.real(showDat)), min(np.imag(showDat)))
            maxy = max(max(np.real(showDat)), max(np.imag(showDat)))
        elif self.viewSettings["plotType"] == 3:
            miny = min(np.abs(showDat))
            maxy = max(np.abs(showDat))
        else:
            miny = -1
            maxy = 1
        differ = 0.05 * (maxy - miny)  # amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny - differ
            self.ymaxlim = maxy + differ
        if self.ref == 0.0:
            self.ppm = False
        axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        if xReset:
            self.xminlim = min(self.xax * axMult)
            self.xmaxlim = max(self.xax * axMult)
        if self.spec > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)


#########################################################################################################
class CurrentScatter(Current1D):

    X_RESIZE = False
    Y_RESIZE = False

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False, output=None):  # display the 1D data
        hyperView = 0
        self.peakPickReset()
        if tmpdata is None:
            tmpdata = self.data1D[hyperView]
        self.ax.cla()
        axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        self.line_xdata = self.xax * axMult
        if old:
            if (self.viewSettings["plotType"] == 0):
                oldData = np.real(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 1):
                oldData = np.imag(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 2):
                oldData = np.real(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 3):
                oldData = np.abs(self.data1D[hyperView])
            self.ax.plot(self.line_xdata, oldData, marker='o', linestyle='none', c='k', alpha=0.2, label=self.data.name + '_old', picker=True)
        if (extraX is not None):
            for num in range(len(extraX)):
                self.ax.plot(extraX[num] * axMult, extraY[num], marker='o', linestyle='none', c=extraColor[num])
        if (self.viewSettings["plotType"] == 0):
            self.line_ydata = np.real(tmpdata)
        elif(self.viewSettings["plotType"] == 1):
            self.line_ydata = np.imag(tmpdata)
        elif(self.viewSettings["plotType"] == 2):
            self.line_ydata = np.real(tmpdata)
            self.ax.plot(self.line_xdata, np.imag(tmpdata), marker='o', linestyle='none', c='r', label=self.data.name + '_imag', picker=True)
        elif(self.viewSettings["plotType"] == 3):
            self.line_ydata = np.abs(tmpdata)
        self.ax.plot(self.line_xdata, self.line_ydata, marker='o', linestyle='none', c=self.viewSettings["color"], label=self.data.name, picker=True)
        self.ax.set_xlabel(self.getLabel(self.spec, self.viewSettings["axType"], self.ppm))
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        if self.spec > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        if output is not None:
            self.canvas.print_figure(output)
        self.canvas.draw()

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
        self.locList = [0] * (len(self.data.shape()) - 1)
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

    def plotReset(self, xReset=True, yReset=True):  # set the plot limits to min and max values
        hyperView = 0
        if self.viewSettings["plotType"] == 0:
            miny = np.min(np.real(self.data1D[hyperView]))
            maxy = np.max(np.real(self.data1D[hyperView]))
        elif self.viewSettings["plotType"] == 1:
            miny = np.min(np.imag(self.data1D[hyperView]))
            maxy = np.max(np.imag(self.data1D[hyperView]))
        elif self.viewSettings["plotType"] == 2:
            miny = np.min((np.min(np.real(self.data1D[hyperView])), np.min(np.imag(self.data1D[hyperView]))))
            maxy = np.max((np.max(np.real(self.data1D[hyperView])), np.max(np.imag(self.data1D[hyperView]))))
        elif self.viewSettings["plotType"] == 3:
            miny = np.min(np.abs(self.data1D[hyperView]))
            maxy = np.max(np.abs(self.data1D[hyperView]))
        else:
            miny = -1
            maxy = 1
        if self.ref == 0.0:
            self.ppm = False
        axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        minx = min(self.xax * axMult)
        maxx = max(self.xax * axMult)
        for i in range(len(self.viewSettings["extraData"])):
            data = self.viewSettings["extraData"][i]
            try:
                if data.ndim() <= self.viewSettings["extraAxes"][i]:
                    self.viewSettings["extraAxes"][i] = len(data.shape()) - 1
                    self.resetExtraLocList(i)
                updateVar = data.getSlice(self.viewSettings["extraAxes"][i], self.viewSettings["extraLoc"][i])
            except Exception:
                self.resetExtraLocList(i)
                updateVar = data.getSlice(self.viewSettings["extraAxes"][i], self.viewSettings["extraLoc"][i])
            data1D = updateVar[0]
            spec = updateVar[3]
            xax = updateVar[5]
            ref = updateVar[6]
            axMult = self.getAxMult(spec, self.viewSettings["axType"], self.ppm, data.freq[self.viewSettings["extraAxes"][i]], ref)
            maxx = max(max(xax * axMult), maxx)
            minx = min(min(xax * axMult), minx)
            if self.viewSettings["plotType"] == 0:
                miny = min(np.min(np.real(data1D[hyperView]) * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i]), miny)
                maxy = max(np.max(np.real(data1D[hyperView]) * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i]), maxy)
            elif self.viewSettings["plotType"] == 1:
                miny = min(np.min(np.imag(data1D[hyperView]) * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i]), miny)
                maxy = max(np.max(np.imag(data1D[hyperView]) * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i]), maxy)
            elif self.viewSettings["plotType"] == 2:
                miny = min(np.min((np.min(np.real(data1D[hyperView]) * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i]), np.min(np.imag(data1D[hyperView]) * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i]))), miny)
                maxy = max(np.max((np.max(np.real(data1D[hyperView]) * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i]), np.max(np.imag(data1D[hyperView]) * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i]))), maxy)
            elif self.viewSettings["plotType"] == 3:
                miny = min(np.min(np.abs(data1D[hyperView]) * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i]), miny)
                maxy = max(np.max(np.abs(data1D[hyperView]) * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i]), maxy)
        differ = 0.05 * (maxy - miny)  # amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny - differ
            self.ymaxlim = maxy + differ
        if xReset:
            self.xminlim = minx
            self.xmaxlim = maxx
        if self.spec > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False, output=None):  # display the 1D data
        hyperView = 0
        self.peakPickReset()
        self.ax.cla()
        for i in range(len(self.viewSettings["extraData"])):
            data = self.viewSettings["extraData"][i]
            try:
                if self.viewSettings["extraData"][i].ndim() <= self.viewSettings["extraAxes"][i]:
                    self.viewSettings["extraAxes"][i] = len(self.viewSettings["extraData"][i].shape()) - 1
                    self.resetExtraLocList(i)
                updateVar = data.getSlice(self.viewSettings["extraAxes"][i], self.viewSettings["extraLoc"][i])
            except Exception:
                self.resetExtraLocList(i)
                updateVar = data.getSlice(self.viewSettings["extraAxes"][i], self.viewSettings["extraLoc"][i])
            data1D = updateVar[0]
            spec = updateVar[3]
            xax = updateVar[5]
            ref = updateVar[6]
            line_xdata = xax * self.getAxMult(spec, self.viewSettings["axType"], self.ppm, data.freq[self.viewSettings["extraAxes"][i]], ref)
            if len(data1D[0]) == 1:
                marker = 'o'
                linestyle = 'none'
            else:
                marker = ''
                linestyle = '-'
            if (self.viewSettings["plotType"] == 0):
                extraData = np.real(data1D[hyperView])
            elif(self.viewSettings["plotType"] == 1):
                extraData = np.imag(data1D[hyperView])
            elif(self.viewSettings["plotType"] == 2):
                extraData = np.real(data1D[hyperView])
            elif(self.viewSettings["plotType"] == 3):
                extraData = np.abs(data1D[hyperView])
            self.ax.plot(line_xdata + self.viewSettings["extraShift"][i],
                         extraData * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i],
                         marker=marker, linestyle=linestyle,
                         c=self.viewSettings["extraColor"][i],
                         linewidth=self.viewSettings["linewidth"],
                         label=data.name,
                         picker=True)
        if tmpdata is None:
            tmpdata = self.data1D[hyperView]
        axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        self.line_xdata = self.xax * axMult
        if self.single:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = ''
            linestyle = '-'
        if old:
            if (self.viewSettings["plotType"] == 0):
                oldData = np.real(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 1):
                oldData = np.imag(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 2):
                oldData = np.real(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 3):
                oldData = np.abs(self.data1D[hyperView])
            self.ax.plot(self.line_xdata, oldData, marker=marker, linestyle=linestyle, c='k', alpha=0.2, linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if (extraX is not None):
            for num in range(len(extraX)):
                self.ax.plot(extraX[num] * axMult, extraY[num], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], c=extraColor[num], picker=True)
        if (self.viewSettings["plotType"] == 0):
            self.line_ydata = np.real(tmpdata)
        elif(self.viewSettings["plotType"] == 1):
            self.line_ydata = np.imag(tmpdata)
        elif(self.viewSettings["plotType"] == 2):
            self.line_ydata = np.real(tmpdata)
            self.ax.plot(self.line_xdata, np.imag(tmpdata), marker=marker, linestyle=linestyle, c='r', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
        elif(self.viewSettings["plotType"] == 3):
            self.line_ydata = np.abs(tmpdata)
        self.ax.plot(self.line_xdata, self.line_ydata, marker=marker, linestyle=linestyle, c=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
        self.ax.set_xlabel(self.getLabel(self.spec, self.viewSettings["axType"], self.ppm))
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        if self.spec > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        if output is not None:
            self.canvas.print_figure(output)
        self.canvas.draw()

#########################################################################################################
# the class from which the stacked data is displayed, the operations which only edit the content of this class are for previewing


class CurrentStacked(Current1D):

    X_RESIZE = False
    Y_RESIZE = True

    def startUp(self, xReset=True, yReset=True):
        self.resetSpacing()
        self.plotReset(xReset, yReset)
        self.showFid()

    def copyCurrent(self, root, fig, canvas, data):
        return CurrentStacked(root, fig, canvas, data, self)

    def upd(self):  # get new data from the data instance
        if self.data.ndim() < 2:
            self.root.rescue()
            return False
        if (len(self.locList) + 2) != self.data.ndim():
            self.resetLocList()
        try:
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList, self.viewSettings["stackBegin"], self.viewSettings["stackEnd"], self.viewSettings["stackStep"])
            if updateVar is None:
                self.root.rescue()
                return False
        except Exception:
            self.resetLocList()
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList, self.viewSettings["stackBegin"], self.viewSettings["stackEnd"], self.viewSettings["stackStep"])
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.freq2 = updateVar[2]
        self.sw = updateVar[3]
        self.sw2 = updateVar[4]
        self.spec = updateVar[5]
        self.spec2 = updateVar[6]
        self.wholeEcho = updateVar[7]
        self.wholeEcho2 = updateVar[8]
        self.xax = updateVar[9]
        self.xax2 = updateVar[10]
        self.ref = updateVar[11]
        self.ref2 = updateVar[12]
        if self.ref is None:
            self.ref = self.freq
        if self.ref2 is None:
            self.ref2 = self.freq2
        self.single = self.data1D[0].shape[-1] == 1
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
        if not axesSame:
            self.resetSpacing()
            self.plotReset()
        self.showFid()

    def resetLocList(self):
        self.locList = [0] * (len(self.data.shape()) - 2)

    def stackSelect(self, stackBegin, stackEnd, stackStep):
        self.viewSettings["stackBegin"] = stackBegin
        self.viewSettings["stackEnd"] = stackEnd
        self.viewSettings["stackStep"] = stackStep
        self.upd()
        self.plotReset(False, True)
        self.showFid()

    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=None):  # display the 1D data including the apodization function
        hyperView = 0
        t = np.arange(0, len(self.data1D[0][0])) / (self.sw)
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                self.dispMsg('shiftingAxes cannot be equal to axes')
            elif shiftingAxes == self.axes2:
                ar = np.arange(self.data.shape()[self.axes2])[reim.floatSlice(self.viewSettings["stackBegin"], self.viewSettings["stackEnd"], self.viewSettings["stackStep"])]
                x = np.ones((len(ar), len(self.data1D[0][0])))
                for i in range(len(ar)):
                    shift1 = shift + shifting * ar[i] / self.data.sw[shiftingAxes]
                    x2 = func.apodize(t, shift1, self.sw, len(self.data1D[0][0]), lor, gauss, cos2, hamming, self.wholeEcho)
                    x[i] = x2
            else:
                if (shiftingAxes < self.axes) and (shiftingAxes < self.axes2):
                    shift += shifting * self.locList[shiftingAxes] * self.data.shape()[shiftingAxes] / self.data.sw[shiftingAxes]
                elif (shiftingAxes > self.axes) and (shiftingAxes > self.axes2):
                    shift += shifting * self.locList[shiftingAxes - 2] * self.data.shape()[shiftingAxes] / self.data.sw[shiftingAxes]
                else:
                    shift += shifting * self.locList[shiftingAxes - 1] * self.data.shape()[shiftingAxes] / self.data.sw[shiftingAxes]
                x = func.apodize(t, shift, self.sw, len(self.data1D[0][0]), lor, gauss, cos2, hamming, self.wholeEcho)
                x = np.repeat([x], len(self.data1D[0]), axis=0)
        else:
            x = func.apodize(t, shift, self.sw, len(self.data1D[0][0]), lor, gauss, cos2, hamming, self.wholeEcho)
            x = np.repeat([x], len(self.data1D[0]), axis=0)
        y = copy.copy(self.data1D)
        self.ax.cla()
        if self.spec == 1:
            y = self.fourierLocal(y, 1,self.axes)
            for index in range(len(y)):
                y[index] = y[index] * x
            y = self.fourierLocal(y, 0,self.axes)
        else:
            for index in range(len(y)):
                y[index] = y[index] * x
        if self.spec == 0:
            if self.viewSettings["plotType"] == 0:
                scale = np.max(np.real(self.data1D[hyperView]))
            elif self.viewSettings["plotType"] == 1:
                scale = np.max(np.imag(self.data1D[hyperView]))
            elif self.viewSettings["plotType"] == 2:
                scale = np.max((np.max(np.real(self.data1D[hyperView])), np.max(np.imag(self.data1D[hyperView]))))
            elif self.viewSettings["plotType"] == 3:
                scale = np.max(np.abs(self.data1D[hyperView]))
            self.showFid(y[hyperView], [t], x * scale, ['g'], old=True)
        else:
            self.showFid(y[hyperView])

    def setSpacing(self, spacing):
        self.viewSettings["spacing"] = spacing
        self.plotReset(False, True)
        self.showFid()

    def resetSpacing(self):
        hyperView = 0
        difference = np.diff(self.data1D[hyperView], axis=0)
        if difference.size == 0:
            self.viewSettings["spacing"] = 0
        else:
            if self.viewSettings["plotType"] == 0:
                difference = np.min(np.real(difference))
                amp = np.max(np.real(self.data1D[hyperView])) - np.min(np.real(self.data1D[hyperView]))
            elif self.viewSettings["plotType"] == 1:
                difference = np.min(np.imag(difference))
                amp = np.max(np.imag(self.data1D[hyperView])) - np.min(np.imag(self.data1D[hyperView]))
            elif self.viewSettings["plotType"] == 2:
                difference = np.min((np.real(difference), np.imag(difference)))
                amp = np.max((np.real(self.data1D[hyperView]), np.imag(self.data1D[hyperView]))) - np.min((np.real(self.data1D[hyperView]), np.imag(self.data1D[hyperView])))
            elif self.viewSettings["plotType"] == 3:
                difference = np.min(np.abs(difference))
                amp = np.max(np.abs(self.data1D[hyperView])) - np.min(np.abs(self.data1D[hyperView]))
            self.viewSettings["spacing"] = np.abs(difference) + 0.1 * amp

    def altScroll(self, event):
        self.viewSettings["spacing"] = self.viewSettings["spacing"] * 1.1**event.step
        self.root.sideframe.scrollSpacing(self.viewSettings["spacing"])
        self.showFid()

    def altReset(self):
        self.resetSpacing()
        self.root.sideframe.scrollSpacing(self.viewSettings["spacing"])
        self.showFid()

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False):  # display the 1D data
        hyperView = 0
        self.peakPickReset()
        if tmpdata is None:
            tmpdata = self.data1D[hyperView]
        self.ax.cla()
        axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        self.line_xdata = self.xax * axMult
        if self.single:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = ''
            linestyle = '-'
        if old:
            for num in range(len(self.data1D[0])):
                if (self.viewSettings["plotType"] == 0):
                    tmpData = np.real(self.data1D[hyperView][num])
                elif(self.viewSettings["plotType"] == 1):
                    tmpData = np.imag(self.data1D[hyperView][num])
                elif(self.viewSettings["plotType"] == 2):
                    tmpData = np.real(self.data1D[hyperView][num])
                elif(self.viewSettings["plotType"] == 3):
                    tmpData = np.abs(self.data1D[hyperView][num])
                self.ax.plot(self.line_xdata, num * self.viewSettings["spacing"] + tmpData, marker=marker, linestyle=linestyle, c='k', alpha=0.2, linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if (extraX is not None):
            for num in range(len(extraY)):
                self.ax.plot(extraX[0] * axMult, num * self.viewSettings["spacing"] + extraY[num], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], c=extraColor[0], picker=True)
        if (self.viewSettings["plotType"] == 0):
            tmpdata = np.real(tmpdata)
        elif (self.viewSettings["plotType"] == 1):
            tmpdata = np.imag(tmpdata)
        elif(self.viewSettings["plotType"] == 3):
            tmpdata = np.abs(tmpdata)
        self.line_ydata = np.real(tmpdata[0])
        for num in range(len(tmpdata)):
            if (self.viewSettings["plotType"] == 2):
                self.ax.plot(self.line_xdata, num * self.viewSettings["spacing"] + np.imag(tmpdata[num]), marker=marker, linestyle=linestyle, c='r', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
            self.ax.plot(self.line_xdata, num * self.viewSettings["spacing"] + np.real(tmpdata[num]), marker=marker, linestyle=linestyle, c=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
        self.ax.set_xlabel(self.getLabel(self.spec, self.viewSettings["axType"], self.ppm))
        if self.spec > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.canvas.draw()

    def plotReset(self, xReset=True, yReset=True):  # set the plot limits to min and max values
        hyperView = 0
        self.ax = self.fig.gca()
        incr = np.repeat(np.arange(len(self.data1D[0])).reshape((len(self.data1D[0]), 1)), len(self.data1D[0][0]), axis=1) * self.viewSettings["spacing"]
        if self.viewSettings["plotType"] == 0:
            miny = np.min(np.real(self.data1D[hyperView]) + incr)
            maxy = np.max(np.real(self.data1D[hyperView]) + incr)
        elif self.viewSettings["plotType"] == 1:
            miny = np.min(np.imag(self.data1D[hyperView]) + incr)
            maxy = np.max(np.imag(self.data1D[hyperView]) + incr)
        elif self.viewSettings["plotType"] == 2:
            miny = np.min((np.min(np.real(self.data1D[hyperView]) + incr), np.min(np.imag(self.data1D[hyperView]) + incr)))
            maxy = np.max((np.max(np.real(self.data1D[hyperView]) + incr), np.max(np.imag(self.data1D[hyperView]) + incr)))
        elif self.viewSettings["plotType"] == 3:
            miny = np.min(np.abs(self.data1D[hyperView]) + incr)
            maxy = np.max(np.abs(self.data1D[hyperView]) + incr)
        else:
            miny = -1
            maxy = 1
        differ = 0.05 * (maxy - miny)  # amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny - differ
            self.ymaxlim = maxy + differ
        if self.ref == 0.0:
            self.ppm = False
        axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        if xReset:
            self.xminlim = min(self.xax * axMult)
            self.xmaxlim = max(self.xax * axMult)
        if self.spec > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)

#########################################################################################################
# the class from which the arrayed data is displayed, the operations which only edit the content of this class are for previewing

class CurrentArrayed(Current1D):

    X_RESIZE = True
    Y_RESIZE = False

    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        if duplicateCurrent is not None:
            if isinstance(duplicateCurrent, CurrentArrayed):
                self.zminlim = duplicateCurrent.zminlim
                self.zmaxlim = duplicateCurrent.zmaxlim
            else:
                # The z-axes limits are in xax units unlike the x-axes and y-axes limits
                axMult = self.getAxMult(duplicateCurrent.spec,
                                        duplicateCurrent.viewSettings["axType"],
                                        duplicateCurrent.ppm,
                                        duplicateCurrent.freq,
                                        duplicateCurrent.ref)
                self.zminlim = (duplicateCurrent.xminlim) / axMult
                self.zmaxlim = (duplicateCurrent.xmaxlim) / axMult
        super(CurrentArrayed, self).__init__(root, fig, canvas, data, duplicateCurrent)

    def startUp(self, xReset=True, yReset=True):
        self.resetSpacing(False)
        self.plotReset(xReset, yReset)
        self.showFid()

    def copyCurrent(self, root, fig, canvas, data):
        return CurrentArrayed(root, fig, canvas, data, self)

    def upd(self):  # get new data from the data instance
        if self.data.ndim() < 2:
            self.root.rescue()
            return False
        if (len(self.locList) + 2) != self.data.ndim():
            self.resetLocList()
        try:
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList, self.viewSettings["stackBegin"], self.viewSettings["stackEnd"], self.viewSettings["stackStep"])
            if updateVar is None:
                self.root.rescue()
                return False
        except Exception:
            self.resetLocList()
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList, self.viewSettings["stackBegin"], self.viewSettings["stackEnd"], self.viewSettings["stackStep"])
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.freq2 = updateVar[2]
        self.sw = updateVar[3]
        self.sw2 = updateVar[4]
        self.spec = updateVar[5]
        self.spec2 = updateVar[6]
        self.wholeEcho = updateVar[7]
        self.wholeEcho2 = updateVar[8]
        self.xax = updateVar[9]
        self.xax2 = updateVar[10]
        self.ref = updateVar[11]
        self.ref2 = updateVar[12]
        if self.ref is None:
            self.ref = self.freq
        if self.ref2 is None:
            self.ref2 = self.freq2
        self.single = self.data1D[0].shape[-1] == 1
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
        if self.data.shape()[self.axes2] > 100:
            self.viewSettings["stackStep"] = 1 + int(self.data.shape()[self.axes2]) / 100
        self.locList = locList
        self.upd()
        if not axesSame:
            self.resetSpacing()
            self.plotReset()
        self.showFid()

    def resetLocList(self):
        self.locList = [0] * (len(self.data.shape()) - 2)

    def setAxType2(self, val):
        oldAxMult = self.getAxMult(self.spec2, self.viewSettings["axType2"], self.ppm2, self.freq2, self.ref2)
        if val == 3:
            self.ppm2 = True
        else:
            self.ppm2 = False
            self.viewSettings["axType2"] = val
        newAxMult = self.getAxMult(self.spec2, self.viewSettings["axType2"], self.ppm2, self.freq2, self.ref2)
        self.showFid()

    def stackSelect(self, stackBegin, stackEnd, stackStep):
        self.viewSettings["stackBegin"] = stackBegin
        self.viewSettings["stackEnd"] = stackEnd
        self.viewSettings["stackStep"] = stackStep
        self.upd()
        self.plotReset(True, False)
        self.showFid()

    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=None):  # display the 1D data including the apodization function
        hyperView = 0
        t = np.arange(0, len(self.data1D[0][0])) / (self.sw)
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                self.dispMsg('shiftingAxes cannot be equal to axes')
            elif shiftingAxes == self.axes2:
                ar = np.arange(self.data.shape()[self.axes2])[slice(self.viewSettings["stackBegin"], self.viewSettings["stackEnd"], self.viewSettings["stackStep"])]
                x = np.ones((len(ar), len(self.data1D[0][0])))
                for i in range(len(ar)):
                    shift1 = shift + shifting * ar[i] / self.data.sw[shiftingAxes]
                    x2 = func.apodize(t, shift1, self.sw, len(self.data1D[0][0]), lor, gauss, cos2, hamming, self.wholeEcho)
                    x[i] = x2
            else:
                if (shiftingAxes < self.axes) and (shiftingAxes < self.axes2):
                    shift += shifting * self.locList[shiftingAxes] * self.data.shape()[shiftingAxes] / self.data.sw[shiftingAxes]
                elif (shiftingAxes > self.axes) and (shiftingAxes > self.axes2):
                    shift += shifting * self.locList[shiftingAxes - 2] * self.data.shape()[shiftingAxes] / self.data.sw[shiftingAxes]
                else:
                    shift += shifting * self.locList[shiftingAxes - 1] * self.data.shape()[shiftingAxes] / self.data.sw[shiftingAxes]
                x = func.apodize(t, shift, self.sw, len(self.data1D[0][0]), lor, gauss, cos2, hamming, self.wholeEcho)
                x = np.repeat([x], len(self.data1D[0]), axis=0)
        else:
            x = func.apodize(t, shift, self.sw, len(self.data1D[0][0]), lor, gauss, cos2, hamming, self.wholeEcho)
            x = np.repeat([x], len(self.data1D[0]), axis=0)
        y = copy.copy(self.data1D)
        self.ax.cla()
        if self.spec == 1:
            y = self.fourierLocal(y, 1,self.axes)
            for index in range(len(y)):
                y[index] = y[index] * x
            y = self.fourierLocal(y, 0, self.axes)
        else:
            for index in range(len(y)):
                y[index] = y[index] * x
        if self.spec == 0:
            if self.viewSettings["plotType"] == 0:
                tmpData = np.max(np.real(self.data1D[hyperView]))
            elif self.viewSettings["plotType"] == 1:
                tmpData = np.max(np.imag(self.data1D[hyperView]))
            elif self.viewSettings["plotType"] == 2:
                tmpData = np.max((np.max(np.real(self.data1D[hyperView])), np.max(np.imag(self.data1D[hyperView]))))
            elif self.viewSettings["plotType"] == 3:
                tmpData = np.max(np.abs(self.data1D[hyperView]))
            self.showFid(y[hyperView], [t], x * tmpData, ['g'], old=True)
        else:
            self.showFid(y[hyperView])
   
    def setSpacing(self, spacing):
        self.viewSettings["spacing"] = spacing
        self.plotReset(True, False)
        self.showFid()

    def resetSpacing(self, zlims=True):
        if zlims:
            self.zminlim = min(self.xax)
            self.zmaxlim = max(self.xax)
        xaxZlims = (self.xax > self.zminlim) & (self.xax < self.zmaxlim)
        self.viewSettings["spacing"] = (self.xax[xaxZlims][-1] - self.xax[xaxZlims][0]) * 1.1

    def altScroll(self, event):
        self.viewSettings["spacing"] = self.viewSettings["spacing"] * 1.1**event.step
        self.root.sideframe.scrollSpacing(self.viewSettings["spacing"])
        self.showFid()

    def altReset(self):
        self.resetSpacing()
        self.root.sideframe.scrollSpacing(self.viewSettings["spacing"])
        self.showFid()

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False):  # display the 1D data
        hyperView = 0
        self.peakPickReset()
        if tmpdata is None:
            tmpdata = self.data1D[hyperView]
        self.ax.cla()
        if self.spec > 0:
            direc = slice(None, None, -1)
        else:
            direc = slice(None, None, 1)
        xaxZlims = (self.xax > self.zminlim) & (self.xax < self.zmaxlim)
        axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        axMult2 = self.getAxMult(self.spec2, self.viewSettings["axType2"], self.ppm2, self.freq2, self.ref2)
        if self.single:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = ''
            linestyle = '-'
        if old:
            if (self.viewSettings["plotType"] == 0):
                oldData = np.real(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 1):
                oldData = np.imag(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 2):
                oldData = np.real(self.data1D[hyperView])
            elif(self.viewSettings["plotType"] == 3):
                oldData = np.abs(self.data1D[hyperView])
            for num in range(len(self.data1D[hyperView])):
                self.ax.plot((num * self.viewSettings["spacing"] + self.xax[xaxZlims]) * axMult, oldData[num][xaxZlims][direc], marker=marker, linestyle=linestyle, c='k', alpha=0.2, linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if (extraX is not None):
            extraZlims = (extraX[0] > self.zminlim) & (extraX[0] < self.zmaxlim)
            for num in range(len(extraY)):
                self.ax.plot((num * self.viewSettings["spacing"] + extraX[0][extraZlims]) * axMult, extraY[num][extraZlims][direc], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], c=extraColor[0], picker=True)
        if (self.viewSettings["plotType"] == 0):
            tmpdata = np.real(tmpdata)
        elif(self.viewSettings["plotType"] == 1):
            tmpdata = np.imag(tmpdata)
        elif(self.viewSettings["plotType"] == 3):
            tmpdata = np.abs(tmpdata)
        self.line_xdata = []
        self.line_ydata = []
        ticksPos = []
        for num in range(len(tmpdata)):
            if (self.viewSettings["plotType"] == 2):
                self.ax.plot((num * self.viewSettings["spacing"] + self.xax[xaxZlims]) * axMult, np.imag(tmpdata[num][xaxZlims])[direc], marker=marker, linestyle=linestyle, c='r', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
            self.line_xdata = np.append(self.line_xdata, (num * self.viewSettings["spacing"] + self.xax[xaxZlims]) * axMult)
            self.line_ydata = np.append(self.line_ydata, np.real(tmpdata[num][xaxZlims])[direc])
            self.ax.plot((num * self.viewSettings["spacing"] + self.xax[xaxZlims]) * axMult, np.real(tmpdata[num][xaxZlims])[direc], marker=marker, linestyle=linestyle, c=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
            pos = (num * self.viewSettings["spacing"] + 0.5 * (self.xax[xaxZlims][-1] + self.xax[xaxZlims][0])) * axMult
            ticksPos.append(pos)
        self.ax.set_xticks(ticksPos)
        self.ax.set_xticklabels([('%#.3g') % x for x in self.xax2 * axMult2])
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.canvas.draw()

    def plotReset(self, xReset=True, yReset=True):  # set the plot limits to min and max values
        hyperView = 0
        if self.viewSettings["plotType"] == 0:
            miny = np.min(np.real(self.data1D[hyperView]))
            maxy = np.max(np.real(self.data1D[hyperView]))
        elif self.viewSettings["plotType"] == 1:
            miny = np.min(np.imag(self.data1D[hyperView]))
            maxy = np.max(np.imag(self.data1D[hyperView]))
        elif self.viewSettings["plotType"] == 2:
            miny = np.min((np.min(np.real(self.data1D[hyperView])), np.min(np.imag(self.data1D[hyperView]))))
            maxy = np.max((np.max(np.real(self.data1D[hyperView])), np.max(np.imag(self.data1D[hyperView]))))
        elif self.viewSettings["plotType"] == 3:
            miny = np.min(np.abs(self.data1D[hyperView]))
            maxy = np.max(np.abs(self.data1D[hyperView]))
        else:
            miny = -1
            maxy = 1
        differ = 0.05 * (maxy - miny)  # amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny - differ
            self.ymaxlim = maxy + differ
        if xReset:
            xaxZlims = (self.xax > self.zminlim) & (self.xax < self.zmaxlim)
            if self.ref == 0.0:
                self.ppm = False
            axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
            self.xminlim = min(self.xax[xaxZlims] * axMult)
            self.xmaxlim = (max(self.xax[xaxZlims]) + (len(self.data1D[0]) - 1) * self.viewSettings["spacing"]) * axMult
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
  
    def stackSelect(self, stackBegin, stackEnd, stackStep):
        self.viewSettings["stackBegin"] = stackBegin
        self.viewSettings["stackEnd"] = stackEnd
        self.viewSettings["stackStep"] = stackStep
        self.upd()
        self.plotReset(True, False)
        self.showFid()

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


class CurrentContour(Current1D):

    X_RESIZE = False
    Y_RESIZE = True

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

    def upd(self):  # get new data from the data instance
        if self.data.ndim() < 2:
            self.root.rescue()
            return False
        if (len(self.locList) + 2) != self.data.ndim():
            self.resetLocList()
        try:
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList)
            if updateVar is None:
                self.root.rescue()
                return False
        except Exception:
            self.resetLocList()
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.freq2 = updateVar[2]
        self.sw = updateVar[3]
        self.sw2 = updateVar[4]
        self.spec = updateVar[5]
        self.spec2 = updateVar[6]
        self.wholeEcho = updateVar[7]
        self.wholeEcho2 = updateVar[8]
        self.xax = updateVar[9]
        self.xax2 = updateVar[10]
        self.ref = updateVar[11]
        self.ref2 = updateVar[12]
        if self.ref is None:
            self.ref = self.freq
        if self.ref2 is None:
            self.ref2 = self.freq2
        self.single = self.data1D[0].shape[-1] == 1
        return True

    def setBlock(self, axes, axes2, locList):  # change the slice
        axesSame = True
        if (self.axes != axes) or (self.axes2 != axes2):
            axesSame = False
        self.axes = axes
        self.axes2 = axes2
        self.locList = locList
        self.upd()
        if not axesSame:
            self.plotReset()
        self.showFid()

    def setLevels(self, numLevels, maxLevels, minLevels, contourSign, contourType, multiValue):
        self.viewSettings["numLevels"] = numLevels
        self.viewSettings["maxLevels"] = maxLevels
        self.viewSettings["minLevels"] = minLevels
        self.viewSettings["contourSign"] = contourSign
        self.viewSettings["contourType"] = contourType
        self.viewSettings["multiValue"] = multiValue
        self.showFid()

    def setProjLimits(self, ProjBool, Limits):
        self.viewSettings["projLimits"] = Limits
        self.viewSettings["projLimitsBool"] = ProjBool

    def resetLocList(self):
        self.locList = [0] * (len(self.data.shape()) - 2)

    def setAxType2(self, val):
        oldAxMult = self.getAxMult(self.spec2, self.viewSettings["axType2"], self.ppm2, self.freq2, self.ref2)
        if val == 3:
            self.ppm2 = True
        else:
            self.ppm2 = False
            self.viewSettings["axType2"] = val
        newAxMult = self.getAxMult(self.spec2, self.viewSettings["axType2"], self.ppm2, self.freq2, self.ref2)
        self.yminlim = self.yminlim * newAxMult / oldAxMult
        self.ymaxlim = self.ymaxlim * newAxMult / oldAxMult
        self.showFid()

    def setProjType(self, val, direc):
        if direc == 1:
            self.viewSettings["projTop"] = val
        if direc == 2:
            self.viewSettings["projRight"] = val

    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=None):  # display the 1D data including the apodization function
        hyperView = 0
        t = np.arange(0, len(self.data1D[0][0])) / (self.sw)
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                self.dispMsg('shiftingAxes cannot be equal to axes')
            elif shiftingAxes == self.axes2:
                ar = np.arange(self.data.shape()[self.axes2])
                x = np.ones((len(ar), len(self.data1D[0][0])))
                for i in range(len(ar)):
                    shift1 = shift + shifting * ar[i] / self.data.sw[shiftingAxes]
                    x2 = func.apodize(t, shift1, self.sw, len(self.data1D[0][0]), lor, gauss, cos2, hamming, self.wholeEcho)
                    x[i] = x2
            else:
                if (shiftingAxes < self.axes) and (shiftingAxes < self.axes2):
                    shift += shifting * self.locList[shiftingAxes] * self.data.shape()[shiftingAxes] / self.data.sw[shiftingAxes]
                elif (shiftingAxes > self.axes) and (shiftingAxes > self.axes2):
                    shift += shifting * self.locList[shiftingAxes - 2] * self.data.shape()[shiftingAxes] / self.data.sw[shiftingAxes]
                else:
                    shift += shifting * self.locList[shiftingAxes - 1] * self.data.shape()[shiftingAxes] / self.data.sw[shiftingAxes]
                x = func.apodize(t, shift, self.sw, len(self.data1D[0][0]), lor, gauss, cos2, hamming, self.wholeEcho)
                x = np.repeat([x], len(self.data1D), axis=0)
        else:
            x = func.apodize(t, shift, self.sw, len(self.data1D[0][0]), lor, gauss, cos2, hamming, self.wholeEcho)
            x = np.repeat([x], len(self.data1D[0]), axis=0)
        y = copy.copy(self.data1D)
        self.ax.cla()
        if self.spec == 1:
            y = self.fourierLocal(y, 1, self.axes)
            for index in range(len(y)):
                y[index] = y[index] * x
            y = self.fourierLocal(y, 0, self.axes)
        else:
            for index in range(len(y)):
                y[index] = y[index] * x
        self.showFid(y[hyperView])

    def showFid(self, tmpdata=None):  # display the 1D data
        self.differ = None
        self.peakPickReset()
        hyperView = 0
        if tmpdata is None:
            self.tmpdata = self.data1D[hyperView]
        else:
            self.tmpdata = tmpdata
        self.ax.cla()
        self.x_ax.cla()
        self.y_ax.cla()
        if self.viewSettings["diagonalBool"]:
            add_diagonal(self.ax, self.viewSettings["diagonalMult"], c='k', ls='--')
        self.x = self.xax * self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        self.line_xdata = self.x
        self.y = self.xax2 * self.getAxMult(self.spec2, self.viewSettings["axType2"], self.ppm2, self.freq2, self.ref2)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        if (self.viewSettings["plotType"] == 0):
            self.tmpdata = np.real(self.tmpdata)
        elif(self.viewSettings["plotType"] == 1):
            self.tmpdata = np.imag(self.tmpdata)
        elif(self.viewSettings["plotType"] == 2):
            self.tmpdata = np.real(self.tmpdata)
        elif(self.viewSettings["plotType"] == 3):
            self.tmpdata = np.abs(self.tmpdata)
        self.differ = np.max(np.abs(self.tmpdata))
        self.plotContour(X=self.X, Y=self.Y)
        self.showProj()
        self.ax.set_xlabel(self.getLabel(self.spec, self.viewSettings["axType"], self.ppm))
        self.ax.set_ylabel(self.getLabel(self.spec2, self.viewSettings["axType2"], self.ppm2))
        if self.spec:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        if self.spec2:
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

    def plotContour(self, X=False, Y=False, updateOnly=False):  # Plots the contour plot
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
            if self.tmpdata.shape[0] > 2:  # if size 2 or lower, convolve fails, just take whole data then
                YposMax = np.where(np.convolve(np.max(self.tmpdata, 1) > contourLevels[0], [True, True, True], 'same'))[0]
            else:
                YposMax = np.arange(self.tmpdata.shape[0])
            if YposMax.size > 0:  # if number of positive contours is non-zero
                if self.tmpdata.shape[1] > 2:
                    XposMax = np.where(np.convolve(np.max(self.tmpdata, 0) > contourLevels[0], [True, True, True], 'same'))[0]
                else:
                    XposMax = np.arange(self.tmpdata.shape[1])
                PlotPositive = True
        PlotNegative = False
        if self.viewSettings["contourSign"] == 0 or self.viewSettings["contourSign"] == 2:
            if not self.viewSettings["plotType"] == 3:  # for Absolute plot no negative
                if self.tmpdata.shape[0] > 2:
                    YposMin = np.where(np.convolve(np.min(self.tmpdata, 1) < -contourLevels[0], [True, True, True], 'same'))[0]
                else:
                    YposMin = np.arange(self.tmpdata.shape[0])
                if YposMin.size > 0:  # if number of negative contours is non-zero
                    if self.tmpdata.shape[1] > 2:
                        XposMin = np.where(np.convolve(np.min(self.tmpdata, 0) < -contourLevels[0], [True, True, True], 'same'))[0]
                    else:
                        XposMin = np.arange(self.tmpdata.shape[1])
                    PlotNegative = True
        def contourTrace(level, color):
            level = c.trace(level)
            segs = level[:len(level) // 2]
            col = mcoll.LineCollection(segs)
            col.set_label(self.data.name)
            col.set_linewidth(self.viewSettings["linewidth"])
            col.set_linestyle('solid')
            col.set_color(color)
            return col
        if self.viewSettings["contourConst"]:
            collections = []
            if PlotPositive:
                c = cntr.Cntr(self.X[YposMax[:, None], XposMax], self.Y[YposMax[:, None], XposMax], self.tmpdata[YposMax[:, None], XposMax])
                for level in contourLevels:
                    collections.append(contourTrace(level, self.viewSettings["contourColors"][0]))
            if PlotNegative:
                c = cntr.Cntr(self.X[YposMin[:, None], XposMin], self.Y[YposMin[:, None], XposMin], self.tmpdata[YposMin[:, None], XposMin])
                for level in -contourLevels[::-1]:
                    collections.append(contourTrace(level, self.viewSettings["contourColors"][1]))
            for col in collections:  # plot all
                self.ax.add_collection(col)
        else:
            vmax = max(np.abs(self.viewSettings["minLevels"] * self.differ), np.abs(self.viewSettings["maxLevels"] * self.differ))
            vmin = -vmax
            colorMap = get_cmap(self.viewSettings["colorMap"])
            collections = []
            if PlotPositive:
                c = cntr.Cntr(self.X[YposMax[:, None], XposMax], self.Y[YposMax[:, None], XposMax], self.tmpdata[YposMax[:, None], XposMax])
                for level in contourLevels:
                    clevel = colorMap((level - vmin) / (vmax - vmin))
                    collections.append(contourTrace(level, clevel))
            if PlotNegative:
                c = cntr.Cntr(self.X[YposMin[:, None], XposMin], self.Y[YposMin[:, None], XposMin], self.tmpdata[YposMin[:, None], XposMin])
                for level in -contourLevels[::-1]:
                    clevel = colorMap((level - vmin) / (vmax - vmin))
                    collections.append(contourTrace(level, clevel))
            for col in collections:  # plot all
                self.ax.add_collection(col)
        if updateOnly:
            self.canvas.draw()

    def showProj(self):
        xLimOld = self.x_ax.get_xlim()
        x = self.x  # Get plot data from plot
        yLimOld = self.y_ax.get_ylim()
        y = self.y  # Get plot data from plot
        self.x_ax.cla()
        self.y_ax.cla()
        tmpdata = self.data1D[0]
        if (self.viewSettings["plotType"] == 0):
            tmpdata = np.real(tmpdata)
        elif(self.viewSettings["plotType"] == 1):
            tmpdata = np.imag(tmpdata)
        elif(self.viewSettings["plotType"] == 2):
            tmpdata = np.real(tmpdata)
        elif(self.viewSettings["plotType"] == 3):
            tmpdata = np.abs(tmpdata)
        Limits = self.viewSettings["projLimits"]
        topSlice = np.s_[:, :]
        rightSlice = np.s_[:, :]
        if self.viewSettings["projLimitsBool"] is True:
            if Limits[0] < Limits[1]:
                topSlice = np.s_[Limits[0]:Limits[1] + 1, :]
            elif Limits[0] > Limits[1]:
                topSlice = np.s_[Limits[1]:Limits[0] + 1, :]
            else:
                topSlice = np.s_[Limits[0]:Limits[0] + 1, :]
            if Limits[2] < Limits[3]:
                rightSlice = np.s_[:, Limits[2]:Limits[3] + 1]
            elif Limits[2] > Limits[3]:
                rightSlice = np.s_[:, Limits[3]:Limits[2] + 1]
            else:
                rightSlice = np.s_[:, Limits[2]:Limits[2] + 1]
        if self.viewSettings["projTop"] == 0:
            xprojdata = np.sum(tmpdata[topSlice], axis=0)
        elif self.viewSettings["projTop"] == 1:
            xprojdata = np.max(tmpdata[topSlice], axis=0)
        elif self.viewSettings["projTop"] == 2:
            xprojdata = np.min(tmpdata[topSlice], axis=0)
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
        if self.viewSettings["projRight"] != 3:
            self.y_ax.plot(yprojdata, y, color=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], picker=True)
            ymin, ymax = np.min(yprojdata), np.max(yprojdata)
            self.y_ax.set_xlim([ymin - 0.15 * (ymax - ymin), ymax + 0.05 * (ymax - ymin)])  # Set projection limits, and force 15% whitespace below plot
            self.y_ax.set_ylim(yLimOld)
        self.canvas.draw()

    def plotReset(self, xReset=True, yReset=True):  # set the plot limits to min and max values
        if self.ref == 0.0:
            self.ppm = False
        axMult = self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
        if xReset:
            self.xminlim = min(self.xax * axMult)
            self.xmaxlim = max(self.xax * axMult)
        if self.ref2 == 0.0:
            self.ppm2 = False
        axMult2 = self.getAxMult(self.spec2, self.viewSettings["axType2"], self.ppm2, self.freq2, self.ref2)
        if yReset:
            self.yminlim = min(self.xax2 * axMult2)
            self.ymaxlim = max(self.xax2 * axMult2)
        if self.spec:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        if self.spec2:
            self.ax.set_ylim(self.ymaxlim, self.yminlim)
        else:
            self.ax.set_ylim(self.yminlim, self.ymaxlim)

    # The peakpicking function needs to be changed for contour plots
    def buttonRelease(self, event):
        hyperView = 0
        if event.button == 1:
            if self.peakPick:
                if self.rect[1] is not None:
                    self.rect[1].remove()
                    self.rect[1] = None
                if self.rect[0] is not None:
                    self.rect[0].remove()
                    self.rect[0] = None
                    self.peakPick = False
                    xdata = self.xax * self.getAxMult(self.spec, self.viewSettings["axType"], self.ppm, self.freq, self.ref)
                    ydata = self.xax2 * self.getAxMult(self.spec2, self.viewSettings["axType2"], self.ppm2, self.freq2, self.ref2)
                    idx = np.argmin(np.abs(xdata - event.xdata))
                    idy = np.argmin(np.abs(ydata - event.ydata))
                    if self.peakPickFunc is not None:
                        if (self.viewSettings["plotType"] == 0):
                            tmpdata = np.real(self.data1D[hyperView][idy, idx])
                        elif(self.viewSettings["plotType"] == 1):
                            tmpdata = np.imag(self.data1D[hyperView][idy, idx])
                        elif(self.viewSettings["plotType"] == 2):
                            tmpdata = np.real(self.data1D[hyperView][idy, idx])
                        elif(self.viewSettings["plotType"] == 3):
                            tmpdata = np.abs(self.data1D[hyperView][idy, idx])
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
                    if self.spec > 0:
                        self.ax.set_xlim(self.xmaxlim, self.xminlim)
                    else:
                        self.ax.set_xlim(self.xminlim, self.xmaxlim)
                    if self.spec2 > 0:
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
