#!/usr/bin/env python

# Copyright 2016 Bas van Meerten and Wouter Franssen

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
import scipy.optimize
import copy
from mpl_toolkits.mplot3d import proj3d
from six import string_types
from spectrumFrame import Plot1DFrame
from safeEval import safeEval
from nus import *
import multiprocessing
from matplotlib import cm
from matplotlib.pyplot import get_cmap
import math

COLORMAPLIST = ['seismic', 'BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'rainbow', 'jet']

#########################################################################
#the generic data class
class Spectrum:
    def __init__(self, name, data, loadFunc , freq , sw , spec=None, wholeEcho=None, ref=None, xaxArray=None, history=None, msgHandler=None):
        self.name = name
        self.dim = len(data.shape)                    #number of dimensions
        self.data = np.array(data, dtype=complex)      #data of dimension dim
        self.loadFunc = loadFunc
        self.freq = np.array(freq)                              #array of center frequency (length is dim, MHz)
        self.sw = sw                                  #array of sweepwidths
        if spec is None:
            self.spec = [0] * self.dim
        else:
            self.spec = spec                              #int array of length dim where 0 = time domain, 1 = complex spectral
        if wholeEcho is None:
            self.wholeEcho = [False] * self.dim
        else:
            self.wholeEcho = wholeEcho                    #boolean array of length dim where True indicates a full Echo
        if ref is None:
            self.ref = self.dim * [None]
        else:
            self.ref = ref
        if xaxArray is None:
            self.xaxArray = [[] for i in range(self.dim)]
            self.resetXax()
        else:
            self.xaxArray = xaxArray
        if history is None:
            self.history = []                            #list of strings describing all performed operations
        else:
            self.history = history
        self.msgHandler = msgHandler

    def dispMsg(self, msg):
        if self.msgHandler == None:
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
                self.history.pop()
                   
    def checkAxes(self, axes):
        if axes < 0:
            axes = axes + self.dim
        if not (0 <= axes < self.dim):
            self.dispMsg('Not a valid axes')
            return None
        return axes
            
    def reload(self, mainProgram):
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.reload(mainProgram))
        self.restoreData(self.loadFunc(mainProgram), None)
        return returnValue
            
    def resetXax(self, axes=None):  
        if axes is not None:
            axes = self.checkAxes(axes)
            if axes == None:
                return None
            val = [axes]
        else:
            val = range(self.dim)
        for i in val:
            if self.spec[i] == 0:
                self.xaxArray[i] = np.arange(self.data.shape[i]) / (self.sw[i])
            elif self.spec[i] == 1:
                self.xaxArray[i] = np.fft.fftshift(np.fft.fftfreq(self.data.shape[i], 1.0/self.sw[i]))
                if self.ref[i] is not None:
                    self.xaxArray[i] += self.freq[i] - self.ref[i]
                    
    def setXax(self, xax, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        if len(xax) != self.data.shape[axes]:
            self.dispMsg("Length of new x-axis does not match length of the data")
            return None
        oldXax = self.xaxArray[axes]
        self.xaxArray[axes] = xax
        self.addHistory("X-axes of dimension "+str(axes+1)+" was set to "+str(xax).replace('\n', ''))
        return lambda self: self.setXax(oldXax, axes)
                
    def insert(self, data, pos, axes, dataImag=0):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        data = np.array(data) + 1j*np.array(dataImag)
        oldSize = self.data.shape[axes]
        self.data = np.insert(self.data, [pos], data, axis=axes)
        self.resetXax(axes)
        self.addHistory("Inserted " + str(self.data.shape[axes]-oldSize) + " datapoints in dimension "+str(axes+1)+" at position "+str(pos))
        return lambda self: self.remove(range(pos, pos+data.shape[axes]), axes)

    def remove(self, pos, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.remove(pos, axes))
        tmpdata = np.delete(self.data, pos, axes)
        if (np.array(tmpdata.shape) != 0).all():
            self.data = tmpdata
            self.xaxArray[axes] = np.delete(self.xaxArray[axes], pos)
            self.addHistory("Removed "+ str(len(pos)) + " datapoints from dimension "+str(axes+1)+" at position "+str(pos))
            return returnValue
        else:
            self.dispMsg('Cannot delete all data')
            return None

    def add(self, data, dataImag=0, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        try:
            data = np.array(data) + 1j*np.array(dataImag)
            self.data[select] = self.data[select] + data
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.addHistory("Added to data["+str(select)+"]: "+str(data).replace('\n', ''))
        return lambda self: self.subtract(data, select=select)
        
    def subtract(self, data, dataImag=0, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        try:
            data = np.array(data) + 1j*np.array(dataImag)
            self.data[select] = self.data[select] - data
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.addHistory("Subtracted from data["+str(select)+"]: "+str(data).replace('\n', ''))
        return lambda self: self.add(data, select=select)

    def multiply(self, mult, axes, multImag=0, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        try:
            mult = np.array(mult) + 1j*np.array(multImag)
            copyData = copy.deepcopy(self)
            returnValue = lambda self: self.restoreData(copyData, lambda self: self.multiply(mult, axes, select=select))
            self.data[select] = np.apply_along_axis(np.multiply, axes, self.data, mult)[select]
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.addHistory("Multiplied dimension "+str(axes+1)+" of data["+str(select)+"]: "+str(mult).replace('\n', ''))
        return returnValue

    def baselineCorrection(self, baseline, axes, baselineImag = 0, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        try:
            baseline = np.array(baseline) + 1j*np.array(baselineImag)
            baselinetmp = baseline.reshape((1, )*axes+(self.data.shape[axes], )+(1, )*(self.dim-axes-1))
            self.data[select] = self.data[select] - baselinetmp
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.addHistory("Baseline corrected dimension "+str(axes+1)+" of data["+str(select)+"]")
        return lambda self: self.baselineCorrection(-baseline, axes, select=select) 
    
    def concatenate(self, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        splitVal = self.data.shape[axes]
        try:
            self.data = np.concatenate(self.data, axis=axes)
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.dim = len(self.data.shape)
        self.freq = np.delete(self.freq, axes)
        self.sw = np.delete(self.sw, axes)
        self.spec = np.delete(self.spec, axes)
        self.wholeEcho =  np.delete(self.wholeEcho, axes)
        self.ref = np.delete(self.ref, axes)
        del self.xaxArray[0]
        #self.resetXax(axes)
        self.resetXax()
        self.addHistory("Concatenated dimension "+str(axes+1))
        return lambda self: self.split(splitVal, axes)
    
    def split(self, sections, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        try:
            self.data = np.array(np.split(self.data, sections, axis=axes))
        except ValueError as error:
            self.dispMsg(str(error))
            return None
        self.dim = len(self.data.shape)
        self.freq = np.insert(self.freq, 0, self.freq[axes])
        self.sw = np.insert(self.sw, 0 , self.sw[axes])
        self.spec = np.insert(self.spec, 0 , self.spec[axes])
        self.wholeEcho =  np.insert(self.wholeEcho, 0 , self.wholeEcho[axes])
        self.ref = np.insert(self.ref, 0 , self.ref[axes])
        self.xaxArray.insert(0, [])
        self.resetXax(0)
        self.resetXax(axes+1)
        self.addHistory("Split dimension " + str(axes+1) + " into " + str(sections) + " sections")
        return lambda self: self.concatenate(axes)
    
    def real(self):
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.real())
        self.data = np.real(self.data)
        self.addHistory("Real")
        return returnValue

    def imag(self):
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.imag())
        self.data = np.imag(self.data)
        self.addHistory("Imaginary")
        return returnValue

    def abs(self):
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.abs())
        self.data = np.abs(self.data)
        self.addHistory("Absolute")
        return returnValue
    
    def states(self, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.states(axes))
        if self.data.shape[axes]%2 != 0:
            self.dispMsg("data has to be even for States")
            return None
        tmpdata = np.real(self.data)
        slicing1 = (slice(None), ) * axes + (slice(None, None, 2), ) + (slice(None), )*(self.dim-1-axes)
        slicing2 = (slice(None), ) * axes + (slice(1, None, 2), ) + (slice(None), )*(self.dim-1-axes)
        tmpdata = tmpdata[slicing1]+1j*tmpdata[slicing2]
        self.data = tmpdata
        self.resetXax(axes)
        self.addHistory("States conversion on dimension " + str(axes+1))
        return returnValue
    
    def statesTPPI(self, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.statesTPPI(axes))
        if self.data.shape[axes]%2 != 0:
            self.dispMsg("data has to be even for States-TPPI")
            return None
        tmpdata = np.real(self.data)
        slicing1 = (slice(None), ) * axes + (slice(None, None, 2), ) + (slice(None), )*(self.dim-1-axes)
        slicing2 = (slice(None), ) * axes + (slice(1, None, 2), ) + (slice(None), )*(self.dim-1-axes)
        tmpdata = tmpdata[slicing1]+1j*tmpdata[slicing2]
        tmpdata[slicing2] = -1*tmpdata[slicing2]
        self.data = tmpdata
        self.resetXax(axes)
        self.addHistory("States-TPPI conversion on dimension " + str(axes+1))
        return returnValue

    def echoAntiEcho(self, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.echoAntiEcho(axes))
        if self.data.shape[axes]%2 != 0:
            self.dispMsg("data has to be even for echo-antiecho")
            return None
        slicing1 = (slice(None), ) * axes + (slice(None, None, 2), ) + (slice(None), )*(self.dim-1-axes)
        slicing2 = (slice(None), ) * axes + (slice(1, None, 2), ) + (slice(None), )*(self.dim-1-axes)
        tmpdata = np.real(self.data[slicing1]+self.data[slicing2])-1j*np.real(self.data[slicing1]-self.data[slicing2])
        self.data = tmpdata
        self.resetXax(axes)
        self.addHistory("Echo-antiecho conversion on dimension " + str(axes+1))
        return returnValue

    def subtractAvg(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        if not (0 <= pos1 <= self.data.shape[axes]):
            self.dispMsg("Indices not within range")
            return None
        if not (0 <= pos2 <= self.data.shape[axes]):
            self.dispMsg("Indices not within range")
            return None
        if pos1 == pos2:
            self.dispMsg("Indices cannot be equal")
            return None
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), )*(self.dim-1-axes)
        averages = np.mean(self.data[slicing], axis=axes, keepdims=True)
        self.data -= averages
        self.addHistory("Subtracted average determined between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        return lambda self: self.add(averages) 
    
    def matrixManip(self, pos1, pos2, axes, which):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        if isinstance(pos1, int):
            pos1 = np.array([pos1])
            pos2 = np.array([pos2])
        if len(pos1) != len(pos2):
            self.dispMsg("Length of the two arrays is not equal")
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.matrixManip(pos1, pos2, axes, which))
        tmpdata = ()
        if len(pos1) == 1:
            keepdims = False
        else:
            keepdims = True
        for i in range(len(pos1)):
            if not (0 <= pos1[i] <= self.data.shape[axes]):
                self.dispMsg("Indices not within range")
                return None
            if not (0 <= pos2[i] <= self.data.shape[axes]):
                self.dispMsg("Indices not within range")
                return None
            if pos1[i] == pos2[i]:
                self.dispMsg("Indices cannot be equal")
                return None
            minPos = min(pos1[i], pos2[i])
            maxPos = max(pos1[i], pos2[i])
            slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), )*(self.dim-1-axes)
            if which == 0:
                if self.spec[axes] == 0:
                    tmpdata += (np.sum(self.data[slicing], axis=axes, keepdims=keepdims)/self.sw[axes], )
                else:
                    tmpdata += (np.sum(self.data[slicing], axis=axes, keepdims=keepdims)*self.sw[axes]/(1.0*self.data.shape[axes]), )
            elif which == 5:
                tmpdata += (np.sum(self.data[slicing], axis=axes, keepdims=keepdims), )
            elif which == 1:
                tmpdata += (np.amax(self.data[slicing], axis=axes, keepdims=keepdims), )
            elif which == 2:
                tmpdata += (np.amin(self.data[slicing], axis=axes, keepdims=keepdims), )
            elif which == 3:
                maxArgPos = np.argmax(np.real(self.data[slicing]), axis=axes)
                tmpmaxPos = maxArgPos.flatten()
                tmp = self.xaxArray[axes][slice(minPos, maxPos)][tmpmaxPos].reshape(maxArgPos.shape)
                if keepdims:
                    tmpdata += (np.expand_dims(tmp, axes), )
                else:
                    tmpdata += (tmp, )
            elif which == 4:
                minArgPos = np.argmin(np.real(self.data[slicing]), axis=axes)
                tmpminPos = minArgPos.flatten()
                tmp = self.xaxArray[axes][slice(minPos, maxPos)][tmpminPos].reshape(minArgPos.shape)
                if keepdims:
                    tmpdata += (np.expand_dims(tmp, axes), )
                else:
                    tmpdata += (tmp, )
            elif which == 6:
                tmpdata += (np.mean(self.data[slicing], axis=axes, keepdims=keepdims), )
        if which == 0:
            self.addHistory("Integrate between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 5:
            self.addHistory("Sum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 1:
            self.addHistory("Maximum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 2:        
            self.addHistory("Minimum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 3:
            self.addHistory("Maximum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 4:
            self.addHistory("Minimum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 6:
            self.addHistory("Average between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        if len(tmpdata) == 1:
            if self.dim == 1:
                self.data = np.array([tmpdata[0]])
                self.resetXax(axes)
            else:
                self.data = tmpdata[0]
                self.dim = self.dim - 1
                self.freq = np.delete(self.freq, axes)
                self.ref = np.delete(self.ref, axes)
                self.sw = np.delete(self.sw, axes)
                self.spec = np.delete(self.spec, axes)
                self.wholeEcho = np.delete(self.wholeEcho, axes)
                del self.xaxArray[axes]
        else:
            self.data = np.concatenate(tmpdata, axis=axes)
            self.resetXax(axes)
        return returnValue

    def matrixManipNew(self, pos1, pos2, axes, which):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        if isinstance(pos1, int):
            pos1 = np.array([pos1])
            pos2 = np.array([pos2])
        if len(pos1) != len(pos2):
            self.dispMsg("Length of the two arrays is not equal")
            return None
        tmpdata = ()
        if len(pos1) == 1:
            keepdims = False
        else:
            keepdims = True
        for i in range(len(pos1)):
            if not (0 <= pos1[i] <= self.data.shape[axes]):
                self.dispMsg("Indices not within range")
                return None
            if not (0 <= pos2[i] <= self.data.shape[axes]):
                self.dispMsg("Indices not within range")
                return None
            if pos1[i] == pos2[i]:
                self.dispMsg("Indices cannot be equal")
                return None
            minPos = min(pos1[i], pos2[i])
            maxPos = max(pos1[i], pos2[i])
            slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), )*(self.dim-1-axes)
            if which == 0:
                if self.spec[axes] == 0:
                    tmpdata += (np.sum(self.data[slicing], axis=axes, keepdims=keepdims)/self.sw[axes], )
                else:
                    tmpdata += (np.sum(self.data[slicing], axis=axes, keepdims=keepdims)*self.sw[axes]/(1.0*self.data.shape[axes]), )
            elif which == 5:
                tmpdata += (np.sum(self.data[slicing], axis=axes, keepdims=keepdims), )
            elif which == 1:
                tmpdata += (np.amax(self.data[slicing], axis=axes, keepdims=keepdims), )
            elif which == 2:
                tmpdata += (np.amin(self.data[slicing], axis=axes, keepdims=keepdims), )
            elif which == 3:
                maxArgPos = np.argmax(np.real(self.data[slicing]), axis=axes)
                tmpmaxPos = maxArgPos.flatten()
                tmp = self.xaxArray[axes][slice(minPos, maxPos)][tmpmaxPos].reshape(maxArgPos.shape)
                if keepdims:
                    tmpdata += (np.expand_dims(tmp, axes), )
                else:
                    tmpdata += (tmp, )
            elif which == 4:
                minArgPos = np.argmin(np.real(self.data[slicing]), axis=axes)
                tmpminPos = minArgPos.flatten()
                tmp = self.xaxArray[axes][slice(minPos, maxPos)][tmpminPos].reshape(minArgPos.shape)
                if keepdims:
                    tmpdata += (np.expand_dims(tmp, axes), )
                else:
                    tmpdata += (tmp, )
            elif which == 6:
                tmpdata += (np.mean(self.data[slicing], axis=axes, keepdims=keepdims), )
        if len(tmpdata) == 1:
            if self.dim == 1:
                tmpdata = np.array([tmpdata[0]])
                newSpec = Spectrum(self.name, tmpdata, self.loadFunc, copy.deepcopy(self.freq), copy.deepcopy(self.sw), copy.deepcopy(self.spec), copy.deepcopy(self.wholeEcho), copy.deepcopy(self.ref), copy.deepcopy(self.xaxArray), copy.deepcopy(self.history), self.msgHandler)
                newSpec.resetXax(axes)
            else:
                tmpXax = copy.deepcopy(self.xaxArray)
                del tmpXax[axes]
                newSpec = Spectrum(self.name, tmpdata[0], self.loadFunc, copy.deepcopy(np.delete(self.freq, axes)), copy.deepcopy(np.delete(self.sw, axes)), copy.deepcopy(np.delete(self.spec, axes)), copy.deepcopy(np.delete(self.wholeEcho, axes)), copy.deepcopy(np.delete(self.ref, axes)), tmpXax, copy.deepcopy(self.history), self.msgHandler)
        else:
            tmpdata = np.concatenate(tmpdata, axis=axes)
            newSpec = Spectrum(self.name, tmpdata, self.loadFunc, copy.deepcopy(self.freq), copy.deepcopy(self.sw), copy.deepcopy(self.spec), copy.deepcopy(self.wholeEcho), copy.deepcopy(self.ref), copy.deepcopy(self.xaxArray), copy.deepcopy(self.history), self.msgHandler)
            newSpec.resetXax(axes)
        if which == 0:
            newSpec.addHistory("Integrate between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 5:
            newSpec.addHistory("Sum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 1:
            newSpec.addHistory("Maximum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 2:        
            newSpec.addHistory("Minimum between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 3:
            newSpec.addHistory("Maximum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))
        elif which == 4:
            newSpec.addHistory("Minimum position between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1)) 
        elif which == 6:
            newSpec.addHistory("Average between " + str(pos1) + " and " + str(pos2) + " of dimension " + str(axes+1))          
        return newSpec
    
    def getRegion(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.getRegion(axes, pos1, pos2))
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), )*(self.dim-1-axes)
        if self.spec[axes] == 1:
            oldFxax = self.xaxArray[axes][slice(minPos, maxPos)][0]
            self.sw[axes] = self.sw[axes]*(maxPos-minPos)/(1.0*self.data.shape[axes])
        self.data = self.data[slicing]
        if self.spec[axes] == 1:
            newFxax = np.fft.fftshift(np.fft.fftfreq(self.data.shape[axes], 1.0/self.sw[axes]))[0]
            if self.ref[axes] is None:
                self.ref[axes] = self.freq[axes]
            self.freq[axes] = self.freq[axes] - newFxax+oldFxax
        self.resetXax(axes)
        self.addHistory("Extracted part between " + str(minPos) + " and " + str(maxPos) + " of dimension " + str(axes+1))
        return returnValue
    
    def getRegionNew(self, pos1, pos2, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), )*(self.dim-1-axes)
        tmpsw = copy.deepcopy(self.sw)
        tmpfreq = copy.deepcopy(self.freq)
        tmpref = copy.deepcopy(self.ref)
        if self.spec[axes] == 1:
            oldFxax = self.xaxArray[axes][slice(minPos, maxPos)][0]
            tmpsw[axes] = tmpsw[axes]*(maxPos-minPos)/(1.0*self.data.shape[axes])
        tmpdata = self.data[slicing]
        if self.spec[axes] == 1:
            newFxax = np.fft.fftshift(np.fft.fftfreq(tmpdata.shape[axes], 1.0/tmpsw[axes]))[0]
            if tmpref[axes] is None:
                tmpref[axes] = tmpfreq[axes]
            tmpfreq[axes] = tmpfreq[axes] - newFxax+oldFxax
        newSpec = Spectrum(self.name, tmpdata, self.loadFunc, tmpfreq, tmpsw, copy.deepcopy(self.spec), copy.deepcopy(self.wholeEcho), tmpref, copy.deepcopy(self.xaxArray), copy.deepcopy(self.history), self.msgHandler)
        newSpec.resetXax(axes)
        newSpec.addHistory("Extracted part between " + str(minPos) + " and " + str(maxPos) + " of dimension " + str(axes+1))
        return newSpec

    def fiddle(self, refSpec, lb,  axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        axLen = self.data.shape[axes]
        if len(refSpec) != axLen:
            self.dispMsg("Reference FID does not have the correct length")
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.fiddle(refSpec, lb, axes))
        tmpSpec = np.fft.ifftshift(np.real(refSpec))
        pos = np.argmax(tmpSpec)
        refFid = np.fft.ifft(tmpSpec)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True)
        t = np.arange(axLen)/self.sw[axes]
        idealFid = np.exp(2j*pos*np.pi*np.linspace(0, 1, axLen+1)[:-1])*np.exp(-np.pi*lb*t)
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
        idealFid /= np.abs(idealFid[0])*2
        self.data *= idealFid
        #self.data = np.apply_along_axis(np.multiply, axes, self.data, idealFid)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True)
        self.addHistory("FIDDLE over dimension " + str(axes+1))
        return returnValue
    
    def diff(self, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.diff(axes))
        self.data = np.diff(self.data, axis=axes)
        self.resetXax(axes)
        self.addHistory("Differences over dimension " + str(axes+1))
        return returnValue

    def cumsum(self, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.cumsum(axes))
        self.data = np.cumsum(self.data, axis=axes)
        self.addHistory("Cumulative sum over dimension " + str(axes+1))
        return returnValue
    
    def flipLR(self, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        slicing = (slice(None), ) * axes + (slice(None, None, -1), ) + (slice(None), )*(self.dim-1-axes)
        self.data = self.data[slicing]
        self.addHistory("Flipped dimension " + str(axes+1))
        return lambda self: self.flipLR(axes)
    
    def hilbert(self, axes):
        import scipy.signal
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.hilbert(axes))
        self.data = scipy.signal.hilbert(np.real(self.data), axis=axes)
        self.addHistory("Hilbert transform on dimension " + str(axes+1))
        return returnValue

    def setPhase(self, phase0, phase1, axes, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if axes is None:
            return None
        if self.ref[axes] is None:
            offset = 0
        else:
            offset = self.freq[axes]-self.ref[axes]
        vector = np.exp(np.fft.fftshift(np.fft.fftfreq(self.data.shape[axes], 1.0/self.sw[axes])+offset)/self.sw[axes]*phase1*1j)
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True)
        self.data[select] = self.data[select]*np.exp(phase0*1j)
        self.data[select] = np.apply_along_axis(np.multiply, axes, self.data, vector)[select]
        if self.spec[axes] == 0:
            self.fourier(axes, tmp=True, inv=True)
            
        Message = "Phasing: phase0 = " + str(phase0*180/np.pi) + " and phase1 = " + str(phase1*180/np.pi) + " for dimension "+ str(axes+1)
        if select != slice(None,None,None):
            Message = Message + " of data["+str(select)+"]"
        self.addHistory(Message)
        return lambda self: self.setPhase(-phase0, -phase1, axes, select=select)

    def apodize(self, lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, axes, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        if shiftingAxes == None:
            shiftingAxes = 0
            shifting = 0
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, axes, select=select))
        axLen = self.data.shape[axes]
        t = np.arange(0, axLen)/self.sw[axes]
        if shifting != 0.0:
            for j in range(self.data.shape[shiftingAxes]):
                shift1 = shift + shifting*j*self.data.shape[shiftingAxes]/self.sw[shiftingAxes]
                t2 = t - shift1
                x = np.ones(axLen)
                if lor is not None:
                    x = x*np.exp(-np.pi*lor*abs(t2))
                if gauss is not None:
                    x = x*np.exp(-((np.pi*gauss*t2)**2)/(4*np.log(2)))
                if cos2 is not None:
                    x = x*(np.cos(cos2*(-0.5*shift1*np.pi*self.sw[axes]/axLen+np.linspace(0, 0.5*np.pi, axLen)))**2)
                if hamming is not None:
                    alpha = 0.53836 # constant for hamming window
                    x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw[axes]/axLen+np.linspace(0, np.pi, axLen))))
                if self.wholeEcho[axes]:
                    x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
                if self.spec[axes] > 0:
                    self.fourier(axes, tmp=True)
                for i in range(self.data.shape[axes]):
                    if axes < shiftingAxes:
                        slicing = (slice(None), ) * axes + (i, ) + (slice(None), )*(shiftingAxes-1-axes) + (j, ) + (slice(None), )*(self.dim-2-shiftingAxes)
                    else:
                        slicing = (slice(None), ) * shiftingAxes + (j, ) + (slice(None), )*(axes-1-shiftingAxes) + (i, ) + (slice(None), )*(self.dim-2-axes)
                    self.data[slicing] = self.data[slicing]*x[i]
                if self.spec[axes] > 0:
                    self.fourier(axes, tmp=True, inv=True)
        else:
            t2 = t - shift
            x = np.ones(axLen)
            if lor is not None:
                x = x*np.exp(-np.pi*lor*abs(t2))
            if gauss is not None:
                x = x*np.exp(-((np.pi*gauss*t2)**2)/(4*np.log(2)))
            if cos2 is not None:
                x = x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw[axes]/axLen+np.linspace(0, 0.5*np.pi, axLen)))**2)
            if hamming is not None:
                alpha = 0.53836 # constant for hamming window
                x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw[axes]/axLen+np.linspace(0, np.pi, axLen))))
            if self.wholeEcho[axes]:
                x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
            if self.spec[axes] > 0:
                self.fourier(axes, tmp=True)
            self.data[select] = np.apply_along_axis(np.multiply, axes, self.data, x)[select]
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
        if lor is None and gauss is None and cos2 is None and hamming is None: #If all none, make special message with `zero apodization'
            Message = Message + "zero apodization"
        Message = Message + " for dimension " + str(axes+1)
        if select != slice(None,None,None):
            Message = Message + " of data["+str(select)+"]"
        self.addHistory(Message)
        return returnValue

    def setFreq(self, freq, sw, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        oldFreq = self.freq[axes]
        oldSw = self.sw[axes]
        self.freq[axes] = float(freq)
        self.sw[axes] = float(sw)
        self.resetXax(axes)
        self.addHistory("Frequency set to "+str(freq*1e-6)+ " MHz and sw set to "+ str(sw*1e-3) + " kHz for dimension " + str(axes+1))
        return lambda self: self.setFreq(oldFreq, oldSw, axes)
    
    def setRef(self, ref, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        oldRef = self.ref[axes]
        if ref is None:
           self.ref[axes] = None
           self.addHistory("Reference frequency set to 'None' for dimension " + str(axes+1))
        else:
            self.ref[axes] = float(ref)
            self.addHistory("Reference frequency set to "+str(ref*1e-6)+ " MHz for dimension " + str(axes+1))
        self.resetXax(axes)
        
        return lambda self: self.setRef(oldRef, axes)

    def setWholeEcho(self, val, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        self.wholeEcho[axes] = val
        self.addHistory("Whole echo set to "+str(val)+ " for dimension " + str(axes+1))
        return lambda self: self.setWholeEcho(not val, axes)
    
    def setSize(self, size, pos, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.setSize(size, pos, axes))
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True)
        if size > self.data.shape[axes]:
            slicing1  = (slice(None), ) * axes + (slice(None, pos), ) + (slice(None), )*(self.dim-1-axes)
            slicing2  = (slice(None), ) * axes + (slice(pos, None), ) + (slice(None), )*(self.dim-1-axes)
            self.data = np.concatenate((np.pad(self.data[slicing1], [(0, 0)]*axes+[(0, size-self.data.shape[axes])]+[(0, 0)]*(self.dim-axes-1), 'constant', constant_values=0), self.data[slicing2]), axes)
        else:
            difference = self.data.shape[axes] - size
            removeBegin = int(np.floor(difference/2))
            removeEnd = difference - removeBegin
            if pos < removeBegin:
                slicing = (slice(None), ) * axes + (slice(self.data.shape[axes]-size, None), ) + (slice(None), )*(self.dim-1-axes)
                self.data = self.data[slicing]
            elif self.data.shape[axes]-pos < removeEnd:
                slicing = (slice(None), ) * axes + (slice(None, size), ) + (slice(None), )*(self.dim-1-axes)
                self.data = self.data[slicing]
            else:
                slicing1  = (slice(None), ) * axes + (slice(None, pos-removeBegin), ) + (slice(None), )*(self.dim-1-axes)
                slicing2  = (slice(None), ) * axes + (slice(pos+removeEnd, None), ) + (slice(None), )*(self.dim-1-axes)
                self.data = np.concatenate((self.data[slicing1], self.data[slicing2]), axis=axes)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True)
        self.dim = len(self.data.shape)
        self.resetXax(axes)
        self.addHistory("Resized dimension " + str(axes+1) + " to " + str(size) + " points at position " +str(pos))
        return returnValue

    def setLPSVD(self,nAnalyse,nFreq,nPredict, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.setLPSVD(nAnalyse,nFreq,nPredict, axes))

        self.data = np.apply_along_axis(self.LPSVDfunction,axes,self.data,nAnalyse,nFreq,nPredict)

        self.resetXax(axes)
        self.addHistory("LPSVD ")
        return returnValue
        
    def LPSVDfunction(self,data,nAnalyse,nFreq,nPredict):
         #LPSVD algorithm
        Y=data[0:nAnalyse]
        N=len(Y)						# # of complex data points in FID
        L=math.floor(N*3/4)						# linear prediction order L = 3/4*N
        A=scipy.linalg.hankel(np.conj(Y[1:N-L+1]),np.conj(Y[N-L:N]))	# backward prediction data matrix
        h=np.conj(Y[0:N-L])					# backward prediction data vector
        U,S,Vh = np.linalg.svd(A,full_matrices=1)                       # singular value decomposition
        V=np.conj(np.transpose(Vh))
        bias=np.mean(S[nFreq:np.min([N-L-1,L])+1])	# bias compensation
        
        PolyCoef=np.dot(-V[:,0:nFreq],np.dot(np.diag(1/(S[0:nFreq]-bias)),np.dot(np.conj(np.transpose(U[:,0:nFreq])),h)))	# prediction polynomial coefficients
        s=np.conj(np.log(np.roots(np.append(PolyCoef[::-1],1))))		# polynomial rooting
        s = s[np.where(s<0)[0]]
        
        Z=np.zeros([N,len(s)],dtype=np.complex)
        for k in range(0,len(s)):
            Z[:,k]=np.exp(s[k])**np.arange(0,N)
        
        
        a=np.linalg.lstsq(Z, Y)[0]
        
        para=np.array([-np.real(s), np.imag(s)/2/np.pi, np.abs(a), np.imag(np.log(a/np.abs(a)))]) #WF: reintroduce the scaling factor
        
        xpredict = np.arange(-nPredict,0)
        reconstructed = np.zeros(nPredict,dtype=np.complex128)
        for signal in range( para.shape[1]):
            reconstructed += para[2,signal]*np.exp(1j*(xpredict*para[1,signal]*2*np.pi+para[3,signal]))*np.exp(-xpredict*para[0,signal])
            
        
        
        data = np.concatenate((reconstructed,data))
        return data
        
        
        
    def changeSpec(self, val, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        oldVal = self.spec[axes]
        self.spec[axes] = val
        self.resetXax(axes)
        if val:
            self.addHistory("Dimension " + str(axes+1) + " set to FID")
        else:
            self.addHistory("Dimension " + str(axes+1) + " set to spectrum")
        return lambda self: self.changeSpec(oldVal, axes)

    def swapEcho(self, idx, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        slicing1 = (slice(None), ) * axes + (slice(None, idx), ) + (slice(None), )*(self.dim-1-axes)
        slicing2 = (slice(None), ) * axes + (slice(idx, None), ) + (slice(None), )*(self.dim-1-axes)
        self.data = np.concatenate((self.data[slicing2], self.data[slicing1]), axes)
        self.wholeEcho[axes] = not self.wholeEcho[axes]
        self.addHistory("Swap echo at position "+str(idx)+" for dimension " + str(axes+1))
        return lambda self: self.swapEcho(-idx, axes)
            
    def shiftData(self, shift, axes, select=slice(None)):
        if isinstance(select, string_types):
            select = safeEval(select)
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.shiftData(shift, axes, select=select))
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True)
        self.data[select] = np.roll(self.data, shift, axes)[select]
        mask = np.ones(self.data.shape[axes])
        if shift < 0:
            mask[slice(shift, None)] = 0
        else:
            mask[slice(None, shift)] = 0
        self.data[select] = np.apply_along_axis(np.multiply, axes, self.data, mask)[select]
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True)
            
        Message = "Shifted "+str(shift)+" points in dimension " + str(axes+1)
        if select != slice(None,None,None):
            Message = Message + " of data["+str(select)+"]"
        self.addHistory(Message)
        return returnValue
#    
#    def LPSVD(self, axes):
#        axes = self.checkAxes(axes)
#        if axes == None:
#            return None
#        M = 1 #Number of frequencies
#        y = self.data[axes][0:100]
#        N = y.shape[0]						# # of complex data points in FID
#        L = math.floor(N*3/4)						# linear prediction order L = 3/4*N
#        A = scipy.linalg.hankel(np.conj(y[1:N-L+1]), np.conj(y[N-L:N]))	# backward prediction data matrix
#        h = np.conj(y[0:N-L])					# backward prediction data vector
#        U, S, V = np.linalg.svd(A)					# singular value decomposition
#        S = np.diag(S)
#        bias = np.mean(S[M:np.min([N-L-1, L])])	# bias compensation
#        PolyCoef = np.dot(-V[:, 0:M], np.dot(np.diag(1/(S[0:M]-bias)), np.dot(np.transpose(U[:, 0:M]), h)))	# prediction polynomial coefficients
#        s = np.conj(np.log(np.roots(np.append(PolyCoef[::-1], 1))))		# polynomial rooting
#        s = s[np.where(np.real(s)<0)[0]]
#        a = 1

    def fourier(self, axes, tmp=False, inv=False):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        if np.logical_xor(self.spec[axes], inv) == 0:
            if not self.wholeEcho[axes] and not tmp:
                slicing = (slice(None), ) * axes + (0, ) + (slice(None), )*(self.dim-1-axes)
                self.data[slicing] = self.data[slicing]*0.5
            self.data = np.fft.fftshift(np.fft.fftn(self.data, axes=[axes]), axes=axes)
            if not tmp:
                self.spec[axes] = 1
                self.addHistory("Fourier transform dimension " + str(axes+1))
        else:
            self.data = np.fft.ifftn(np.fft.ifftshift(self.data, axes=axes), axes=[axes])
            if not self.wholeEcho[axes] and not tmp:
                slicing = (slice(None), ) * axes + (0, ) + (slice(None), )*(self.dim-1-axes)
                self.data[slicing] = self.data[slicing]*2.0
            if not tmp:
                self.spec[axes] = 0
                self.addHistory("Inverse Fourier transform dimension " + str(axes+1))
        self.resetXax(axes)
        return lambda self: self.fourier(axes)

    def realFourier(self, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.realFourier(axes))
        if self.spec[axes] == 0:
            if not self.wholeEcho[axes]:
                slicing = (slice(None), ) * axes + (0, ) + (slice(None), )*(self.dim-1-axes)
                self.data[slicing] = self.data[slicing]*0.5
            self.data = np.fft.fftshift(np.fft.fftn(np.real(self.data), axes=[axes]), axes=axes)
            self.spec[axes] = 1
            self.addHistory("Real Fourier transform dimension " + str(axes+1))
        else:
            self.data = np.fft.ifftn(np.fft.ifftshift(np.real(self.data), axes=axes), axes=[axes])
            if not self.wholeEcho[axes]:
                slicing = (slice(None), ) * axes + (0, ) + (slice(None), )*(self.dim-1-axes)
                self.data[slicing] = self.data[slicing]*2.0
            self.spec[axes] = 0
            self.addHistory("Real inverse Fourier transform dimension " + str(axes+1))
        self.resetXax(axes)
        return returnValue

    def fftshift(self, axes, inv=False):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        if inv:
            self.data = np.fft.ifftshift(self.data, axes=[axes])
            self.addHistory("Inverse Fourier shift dimension " + str(axes+1))
        else:
            self.data = np.fft.fftshift(self.data, axes=axes)
            self.addHistory("Fourier shift dimension " + str(axes+1))
        return lambda self: self.fftshift(axes, not(inv))

    def shear(self, shear, axes, axes2):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        axes2 = self.checkAxes(axes2)
        if axes2 == None:
            return None
        if axes == axes2:
            self.dispMsg('Both shearing axes cannot be equal')
            return None
        if self.dim < 2:
            self.dispMsg("The data does not have enough dimensions for a shearing transformation")
            return None
        shape = self.data.shape
        vec1 = np.linspace(0, shear*2*np.pi*shape[axes]/self.sw[axes], shape[axes]+1)[:-1]
        vec2 = np.fft.fftshift(np.fft.fftfreq(shape[axes2], 1/self.sw[axes2]))
        newShape = [1, ]*self.dim
        newShape[axes] = shape[axes]
        newShape[axes2] = shape[axes2]
        if axes > axes2:
            shearMatrix = np.exp(1j*np.outer(vec2, vec1))
        elif axes < axes2:
            shearMatrix = np.exp(1j*np.outer(vec1, vec2))
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True)
        self.data = self.data*shearMatrix.reshape(shape)
        if self.spec[axes] > 0:
            self.fourier(axes, tmp=True, inv=True)
        self.addHistory("Shearing transform with shearing value "+str(shear)+" over dimensions " + str(axes+1) + " and " + str(axes2+1))
        return lambda self: self.shear(-shear, axes, axes2)

    def reorder(self, pos, newLength, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        if newLength is None:
            newLength = max(pos)+1
        if (max(pos) >= newLength) or (min(pos)< 0):
            self.dispMsg("Invalid positions")
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, lambda self: self.reorder(pos, newLength, axes))
        newShape = np.array(self.data.shape)
        newShape[axes] = newLength
        tmpData = np.zeros(newShape, dtype=complex)
        slicing = (slice(None), ) * axes + (pos, ) + (slice(None), )*(self.dim-1-axes)
        tmpData[slicing] = self.data
        self.data = tmpData
        self.resetXax(axes)
        self.addHistory("Reorder dimension " + str(axes+1) + " to obtain a new length of "+str(newLength)+" with positions "+str(pos))
        return returnValue
    
    def ffm_1d(self, pos, typeVal, axes):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, None)
        #pos contains the values of fixed points which not to be translated to missing points
        posList = np.delete(range(self.data.shape[axes]), pos)
        if typeVal == 1: #type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList/2), dtype=int)
        if typeVal == 2: #type is TPPI, for now handle the same as Complex
            pass
        posList = np.unique(posList)
        tmpData = np.rollaxis(self.data, axes, self.dim)
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((tmpData.size/tmpShape[-1], tmpShape[-1]))
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(ffm, [(i, posList) for i in tmpData])
        pool.close()
        pool.join()
        self.data = np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axes)
        self.addHistory("Fast Forward Maximum Entropy reconstruction of dimension " + str(axes+1) + " at positions "+str(pos))
        return returnValue
    
    def clean(self, pos, typeVal, axes, gamma, threshold, maxIter):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        copyData = copy.deepcopy(self)
        returnValue = lambda self: self.restoreData(copyData, None)
        #pos contains the values of fixed points which not to be translated to missing points
        posList = np.delete(range(self.data.shape[axes]), pos)
        if typeVal == 1: #type is States or States-TPPI, the positions need to be divided by 2
            posList = np.array(np.floor(posList/2), dtype=int)
        if typeVal == 2: #type is TPPI, for now handle the same as Complex
            pass
        posList = np.unique(posList)
        tmpData = np.rollaxis(np.fft.fft(self.data, axis=axes), axes, self.dim) 
        tmpShape = tmpData.shape
        tmpData = tmpData.reshape((tmpData.size/tmpShape[-1], tmpShape[-1]))
        mask = np.ones(tmpShape[-1])/float(tmpShape[-1])
        #mask = np.exp(-np.pi*lb*np.arange(tmpShape[-1]) / (self.sw[axes]))/float(tmpShape[-1])
        mask[posList] = 0.0
        mask = np.fft.fft(mask)                                                    # abs or real???
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        fit = pool.map_async(clean, [(i, mask, gamma, threshold, maxIter) for i in tmpData])
        pool.close()
        pool.join()
        self.data = np.fft.ifft(np.rollaxis(np.array(fit.get()).reshape(tmpShape), -1, axes), axis=axes)
        self.addHistory("CLEAN reconstruction of dimension " + str(axes+1) + " at positions "+str(pos))
        return returnValue
    
    def getSlice(self, axes, locList):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        return copy.deepcopy((self.data[tuple(locList[:axes])+(slice(None), )+tuple(locList[axes:])], self.freq[axes], self.sw[axes], self.spec[axes], self.wholeEcho[axes], self.xaxArray[axes], self.ref[axes]))

    def getBlock(self, axes, axes2, locList, stackBegin=None, stackEnd=None, stackStep=None):
        axes = self.checkAxes(axes)
        if axes == None:
            return None
        axes2 = self.checkAxes(axes2)
        if axes2 == None:
            return None
        stackSlice = slice(stackBegin, stackEnd, stackStep)
        if axes == axes2:
            self.dispMsg("First and second axes are the same")
            return None
        elif axes < axes2:
            return copy.deepcopy((np.transpose(self.data[tuple(locList[:axes])+(slice(None), )+tuple(locList[axes:axes2-1])+(stackSlice, )+tuple(locList[axes2-1:])]), self.freq[axes], self.freq[axes2], self.sw[axes], self.sw[axes2], self.spec[axes], self.spec[axes2], self.wholeEcho[axes], self.wholeEcho[axes2], self.xaxArray[axes], self.xaxArray[axes2][stackSlice], self.ref[axes], self.ref[axes2]))
        elif axes > axes2:
            return copy.deepcopy((self.data[tuple(locList[:axes2])+(stackSlice, )+tuple(locList[axes2:axes-1])+(slice(None), )+tuple(locList[axes-1:])], self.freq[axes], self.freq[axes2], self.sw[axes], self.sw[axes2], self.spec[axes], self.spec[axes2], self.wholeEcho[axes], self.wholeEcho[axes2], self.xaxArray[axes], self.xaxArray[axes2][stackSlice], self.ref[axes], self.ref[axes2]))
        
    def restoreData(self, copyData, returnValue): # restore data from an old copy for undo purposes
        if returnValue is None:
            copyData2 = copy.deepcopy(self)
            returnValue = lambda self: self.restoreData(copyData2, None)
        self.data = copyData.data
        self.dim = len(self.data.shape)                    #number of dimensions
        self.freq = copyData.freq                              #array of center frequency (length is dim, MHz)
        self.loadFunc = copyData.loadFunc
        self.sw = copyData.sw                                  #array of sweepwidths
        self.spec = copyData.spec
        self.wholeEcho = copyData.wholeEcho
        self.xaxArray = copyData.xaxArray
        self.ref = copyData.ref
        self.addHistory("Data was restored to a previous state ")
        return returnValue

##################################################################################################
#the class from which the 1d data is displayed, the operations which only edit the content of this class are for previewing
class Current1D(Plot1DFrame):

    X_RESIZE = False
    Y_RESIZE = False
    
    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        Plot1DFrame.__init__(self, root, fig, canvas)
        self.xax = None               #x-axis
        self.data = data              #the actual spectrum instance
        self.freq = None              #frequency of the slice 
        self.sw = None                #x-data display
        self.data1D = None            #the data1D
        self.spec = None              #boolean where False=time domain and True=spectral domain
        self.wholeEcho = None
        self.ref = None               #reference frequency
        if duplicateCurrent is None:
            self.ppm = False             # display frequency as ppm
            self.axes = len(self.data.data.shape)-1
            self.resetLocList() 
            self.plotType = 0
            self.axType = 1
            self.color = self.root.father.defaultColor                  # color of the main line
            self.linewidth = self.root.father.defaultLinewidth  
            self.grids = self.root.father.defaultGrids                  # display x and y grid
            self.colorMap = self.root.father.defaultColorMap            # colormap for contour like plots
            self.contourConst = self.root.father.defaultContourConst    # bool contour levels have constant color
            self.contourColors = [self.root.father.defaultPosColor, self.root.father.defaultNegColor]  # The colors of the constant color contours
            self.upd()                                                  # get the first slice of data
            self.fig.suptitle(self.data.name)
            self.startUp()
        else:
            self.ppm = duplicateCurrent.ppm
            self.axes = duplicateCurrent.axes
            if isinstance(self, (CurrentStacked, CurrentArrayed, CurrentContour, CurrentSkewed)):
                if (len(duplicateCurrent.locList) == self.data.dim-2):
                    self.locList = duplicateCurrent.locList
                else:
                    if self.axes < self.axes2:
                        self.locList = np.delete(duplicateCurrent.locList, self.axes2-1)
                    else:
                        self.locList = np.delete(duplicateCurrent.locList, self.axes2)
            else:
                if (len(duplicateCurrent.locList) == self.data.dim-1):
                    self.locList = duplicateCurrent.locList
                else:
                    if self.axes < duplicateCurrent.axes2:
                        self.locList = np.insert(duplicateCurrent.locList, duplicateCurrent.axes2-1, 0)
                    else:
                        self.locList = np.insert(duplicateCurrent.locList, duplicateCurrent.axes2, 0)
            self.color = duplicateCurrent.color
            self.linewidth = duplicateCurrent.linewidth
            self.colorMap = duplicateCurrent.colorMap
            self.plotType = duplicateCurrent.plotType
            self.axType = duplicateCurrent.axType
            self.grids = duplicateCurrent.grids
            self.contourConst = duplicateCurrent.contourConst
            self.contourColors = duplicateCurrent.contourColors
            self.xminlim = duplicateCurrent.xminlim
            self.xmaxlim = duplicateCurrent.xmaxlim
            self.yminlim = duplicateCurrent.yminlim
            self.ymaxlim = duplicateCurrent.ymaxlim
            xReset = self.X_RESIZE or duplicateCurrent.X_RESIZE
            yReset = self.Y_RESIZE or duplicateCurrent.Y_RESIZE
            self.upd()   #get the first slice of data
            self.fig.suptitle(self.data.name)
            self.startUp(xReset, yReset)
            
    def dispMsg():
        self.data1D.dispMsg()
            
    def startUp(self, xReset=True, yReset=True):
        self.plotReset(xReset, yReset) #reset the axes limits
        self.showFid() #plot the data

    def rename(self, name):
        self.data.rename(name)
        self.fig.suptitle(name)
        self.canvas.draw()
        
    def copyCurrent(self, root, fig, canvas, data):
        return Current1D(root, fig, canvas, data, self)
        
    def upd(self): #get new data from the data instance
        if self.data.dim <= self.axes:
            self.axes = len(self.data.data.shape)-1
        if (len(self.locList)+1) != self.data.dim:
            self.resetLocList()
        try:
            updateVar = self.data.getSlice(self.axes, self.locList)
        except:
            self.resetLocList()
            updateVar = self.data.getSlice(self.axes, self.locList)
        self.data1D = updateVar[0]
        self.freq = updateVar[1]
        self.sw = updateVar[2]
        self.spec = updateVar[3]
        self.wholeEcho = updateVar[4]
        self.xax = updateVar[5]
        self.ref = updateVar[6]
        if self.ref == None:
            self.ref = self.freq
        self.single = self.data1D.shape[-1] == 1
        return True
        
    def setSlice(self, axes, locList): #change the slice 
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
        self.locList = [0]*(len(self.data.data.shape)-1)

    def getSelect(self):
        tmp = list(self.locList)
        if len(self.data1D.shape) > 1:
            minVal = min(self.axes, self.axes2)
            maxVal = max(self.axes, self.axes2)
            tmp.insert(minVal, slice(None))
            tmp.insert(maxVal, slice(None))
        else:
            tmp.insert(self.axes, slice(None))
        return tmp

    def setGrids(self, grids):
        self.grids = grids
    
    def setPhaseInter(self, phase0in, phase1in): #interactive changing the phase without editing the actual data
        phase0 = float(phase0in)
        phase1 = float(phase1in)
        self.upd()
        if self.spec == 0:
            tmpdata = self.fourierLocal(self.data1D, 0)
        else:
            tmpdata = self.data1D
        tmpdata = tmpdata*np.exp(phase0*1j)
        if self.ref is None:
            offset = 0
        else:
            offset = +self.freq-self.ref
        if len(self.data1D.shape) > 1:
            mult = np.repeat([np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata[0]), 1.0/self.sw)+offset)/self.sw*phase1*1j)], len(tmpdata), axis=0)
        else:
            mult = np.exp(np.fft.fftshift(np.fft.fftfreq(len(tmpdata), 1.0/self.sw)+offset)/self.sw*phase1*1j)
        tmpdata = tmpdata*mult
        if self.spec == 0:
            tmpdata = self.fourierLocal(tmpdata, 1)
        self.data1D = tmpdata
        self.showFid()

    def applyPhase(self, phase0, phase1, select=False):# apply the phase to the actual data
        phase0 = float(phase0)
        phase1 = float(phase1)
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.setPhase(phase0, phase1, self.axes, selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['phase', (phase0, phase1, self.axes-self.data.dim, str(selectSlice))])
        return returnValue

    def fourier(self): #fourier the actual data and replot
        returnValue = self.data.fourier(self.axes)
        self.upd()
        if isinstance(self, (CurrentStacked, CurrentArrayed)):
            self.resetSpacing()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['fourier', (self.axes-self.data.dim, )])
        return returnValue

    def realFourier(self): #fourier the real data and replot
        returnValue = self.data.realFourier(self.axes)
        self.upd()
        if isinstance(self, (CurrentStacked, CurrentArrayed)):
            self.resetSpacing()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['realFourier', (self.axes-self.data.dim, )])
        return returnValue

    def fftshift(self, inv=False): #fftshift the actual data and replot
        returnValue = self.data.fftshift(self.axes, inv)
        self.upd()
        self.showFid()
        self.root.addMacro(['fftshift', (self.axes-self.data.dim, inv)])
        return returnValue
    
    def fourierLocal(self, fourData, spec): #fourier the local data for other functions
        ax = len(fourData.shape)-1
        a = fourData
        if spec == 0:
            if not self.wholeEcho:
                slicing = (slice(None), ) * ax + (0, )
                fourData[slicing] = fourData[slicing]*0.5
            fourData = np.fft.fftshift(np.fft.fftn(fourData, axes=[ax]), axes=ax)
        else:
            fourData = np.fft.ifftn(np.fft.ifftshift(fourData, axes=ax), axes=[ax])
            if not self.wholeEcho:
                slicing = (slice(None), ) * ax + (0, ) 
                fourData[slicing] = fourData[slicing]*2.0
        return fourData

    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None , shift=0.0, shifting=0.0, shiftingAxes=None): #display the 1D data including the apodization function
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                self.dispMsg('shiftingAxes cannot be equal to axes')
                return
            elif shiftingAxes < self.axes:
                shift += shifting*self.locList[shiftingAxes]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
            else:
                shift += shifting*self.locList[shiftingAxes-1]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
        length = len(self.data1D)
        t = np.arange(0, length)/(self.sw)
        t2 = t-shift
        x = np.ones(length)
        if lor is not None:
            x = x*np.exp(-np.pi*lor*abs(t2))
        if gauss is not None:
            x = x*np.exp(-((np.pi*gauss*t2)**2)/(4*np.log(2)))
        if cos2 is not None:
            x = x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/length+np.linspace(0, 0.5*np.pi, len(self.data1D))))**2)
        if hamming is not None:
            alpha = 0.53836 # constant for hamming window
            x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/length+np.linspace(0, np.pi, length))))
        if self.wholeEcho:
            x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
        self.ax.cla()
        y = self.data1D
        if self.spec == 1:
            y = np.fft.ifft(np.fft.ifftshift(y))
            y = y*x
            y = np.fft.fftshift(np.fft.fft(y))
        else:
            y = y*x
        if self.spec == 0:
            if self.plotType == 0:
                self.showFid(y, [t], [x*max(np.real(self.data1D))], ['g'], old=True)
            elif self.plotType == 1:
                self.showFid(y, [t], [x*max(np.imag(self.data1D))], ['g'], old=True)
            elif self.plotType == 2:
                self.showFid(y, [t], [x*max(max(np.real(self.data1D)), max(np.imag(self.data1D)))], ['g'], old=True)
            elif self.plotType == 3:
                self.showFid(y, [t], [x*max(np.abs(self.data1D))], ['g'], old=True)
        else:
            self.showFid(y)

    def applyApod(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=0, select=False): #apply the apodization to the actual data
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, self.axes, selectSlice)
        self.upd() 
        self.showFid()
        self.root.addMacro(['apodize', (lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, self.axes-self.data.dim, str(selectSlice))])
        return returnValue

    def setFreq(self, freq, sw): #set the frequency of the actual data
        returnValue = self.data.setFreq(freq, sw, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['freq', (freq, sw, self.axes-self.data.dim)])
        return returnValue

    def setRef(self, ref): #set the frequency of the actual data
        returnValue = self.data.setRef(ref, self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['ref', (ref, self.axes-self.data.dim)])
        return returnValue

    def SN(self, minNoise, maxNoise, minPeak, maxPeak):
        minN = min(minNoise, maxNoise)
        maxN = max(minNoise, maxNoise)
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        if len(self.data1D.shape) > 1:
            tmpData = self.data1D[0]
        else:
            tmpData = self.data1D
        if (self.plotType == 0):
            tmpData = np.real(tmpData)
        elif(self.plotType == 1):
            tmpData = np.imag(tmpData)
        elif(self.plotType == 2):
            tmpData = np.real(tmpData)
        elif(self.plotType == 3):
            tmpData = np.abs(tmpData)
        return (np.amax(tmpData[minP:maxP])/(np.std(tmpData[minN:maxN])))
    
    def fwhm(self, minPeak, maxPeak):
        from scipy.interpolate import UnivariateSpline
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        if len(self.data1D.shape) > 1:
            tmpData = self.data1D[0]
        else:
            tmpData = self.data1D
        if (self.plotType == 0):
            tmpData = np.real(tmpData)
        elif(self.plotType == 1):
            tmpData = np.imag(tmpData)
        elif(self.plotType == 2):
            tmpData = np.real(tmpData)
        elif(self.plotType == 3):
            tmpData = np.abs(tmpData)
        maxPos = np.argmax(tmpData[minP:maxP])
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        x = self.xax*axMult
        maxX = x[minP:maxP][maxPos]
        spline = UnivariateSpline(x, tmpData-tmpData[minP:maxP][maxPos]/2.0, s=0)
        zeroPos = spline.roots()
        left = zeroPos[zeroPos>maxX]
        right = zeroPos[zeroPos<maxX]
        if right.size>0 and left.size>0:
            return abs(left[0]-right[-1])
        else:
            return 0.0
            
    def COM(self, minPeak, maxPeak): #Centre of Mass
        from scipy.interpolate import UnivariateSpline
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        if len(self.data1D.shape) > 1:
            tmpData = self.data1D[0]
            tmpAxis = self.xax[0]
        else:
            tmpData = self.data1D
            tmpAxis = self.xax
        if (self.plotType == 0):
            tmpData = np.real(tmpData)
        elif(self.plotType == 1):
            tmpData = np.imag(tmpData)
        elif(self.plotType == 2):
            tmpData = np.real(tmpData)
        elif(self.plotType == 3):
            tmpData = np.abs(tmpData)
        
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        tmpAxis = tmpAxis[minP:maxP]*axMult
        
        tmpData = tmpData[minP:maxP]
        #COM = 1/M *sum(m_i * r_i)
        CentreOM = 1.0/np.sum(tmpData) * np.sum(tmpData*tmpAxis)
        return CentreOM
        
    
    def setSizePreview(self, size, pos): #set size only on local data
        if len(self.data1D.shape) > 1:
            length = len(self.data1D[0])
        else:
            length = len(self.data1D)
        if self.spec == 1:
            tmpdata = self.fourierLocal(self.data1D, 1)
        else:
            tmpdata = self.data1D
        axes = len(self.data1D.shape)-1
        if size > length:
            slicing1  = (slice(None), ) * axes + (slice(None, pos), ) + (slice(None), )*(tmpdata.ndim-1-axes)
            slicing2  = (slice(None), ) * axes + (slice(pos, None), ) + (slice(None), )*(tmpdata.ndim-1-axes)
            tmpdata = np.concatenate((np.pad(tmpdata[slicing1], [(0, 0)]*axes+[(0, size-tmpdata.shape[axes])]+[(0, 0)]*(tmpdata.ndim-axes-1), 'constant', constant_values=0), tmpdata[slicing2]), axes)
        else:
            difference = tmpdata.shape[axes] - size
            removeBegin = int(np.floor(difference/2))
            removeEnd = difference - removeBegin
            if pos < removeBegin:
                slicing = (slice(None), ) * axes + (slice(tmpdata.shape[axes]-size, None), ) + (slice(None), )*(tmpdata.ndim-1-axes)
                tmpdata = tmpdata[slicing]
            elif tmpdata.shape[axes]-pos < removeEnd:
                slicing = (slice(None), ) * axes + (slice(None, size), ) + (slice(None), )*(tmpdata.ndim-1-axes)
                tmpdata = tmpdata[slicing]
            else:
                slicing1  = (slice(None), ) * axes + (slice(None, pos-removeBegin), ) + (slice(None), )*(tmpdata.ndim-1-axes)
                slicing2  = (slice(None), ) * axes + (slice(pos+removeEnd, None), ) + (slice(None), )*(tmpdata.ndim-1-axes)
                tmpdata = np.concatenate((tmpdata[slicing1], tmpdata[slicing2]), axis=axes)
        if self.spec == 1:
            self.data1D = self.fourierLocal(tmpdata, 0)
        else:
            self.data1D = tmpdata
        if len(self.data1D.shape) > 1:
            length = len(self.data1D[0])
        else:
            length = len(self.data1D)
        if self.spec == 0:
            self.xax = np.arange(length)/self.sw
        elif self.spec == 1:
            self.xax = np.fft.fftshift(np.fft.fftfreq(length, 1.0/self.sw))+self.freq-self.ref
        if not self.spec:
            self.plotReset(True, False)
        self.showFid()
        self.upd()

    def applySize(self, size, pos): #set size to the actual data
        returnValue = self.data.setSize(size, pos, self.axes)
        self.upd()
        if not self.spec:
            self.plotReset(True, False)
        self.showFid()
        self.root.addMacro(['size', (size, pos, self.axes-self.data.dim)])
        return returnValue
        
    def applyLPSVD(self, nAnalyse, nFreq, nPredict):
        returnValue = self.data.setLPSVD(nAnalyse, nFreq, nPredict, self.axes)
        self.upd()
        self.showFid()
        return returnValue
        
    def changeSpec(self, val): #change from time to freq domain of the actual data
        returnValue = self.data.changeSpec(val, self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['spec', (val, self.axes-self.data.dim)])
        return returnValue
    
    def applySwapEcho(self, idx):
        returnValue = self.data.swapEcho(idx, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['swapecho', (idx, self.axes-self.data.dim)])
        return returnValue

    def setSwapEchoPreview(self, idx):
        if len(self.data1D.shape) > 1:
            self.data1D = np.concatenate((self.data1D[:, idx:], self.data1D[:, :idx]), axis=1)
        else:
            self.data1D = np.concatenate((self.data1D[idx:], self.data1D[:idx]))
        self.showFid()
        self.upd()

    def setWholeEcho(self, value):
        if value == 0:
            returnValue = self.data.setWholeEcho(False, self.axes)
            self.wholeEcho = False
            self.root.addMacro(['wholeEcho', (False, self.axes-self.data.dim)])
        else:
            returnValue = self.data.setWholeEcho(True, self.axes)
            self.wholeEcho = True
            self.root.addMacro(['wholeEcho', (True, self.axes-self.data.dim)])
        return returnValue

    def applyShift(self, shift, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.shiftData(shift, self.axes, selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['shift', (shift, self.axes-self.data.dim, str(selectSlice))])
        return returnValue

    def setShiftPreview(self, shift):
        tmpData = self.data1D
        dim = len(self.data1D.shape)
        if self.spec > 0:
            tmpData = self.fourierLocal(tmpData, 1)
        tmpData = np.roll(tmpData, shift)
        if shift<0:
            tmpData[(slice(None), )*(dim-1)+(slice(shift, None), )] = tmpData[(slice(None), )*(dim-1)+(slice(shift, None), )]*0
        else:
            tmpData[(slice(None), )*(dim-1)+(slice(None, shift), )] = tmpData[(slice(None), )*(dim-1)+(slice(None, shift), )]*0
        if self.spec > 0:
           tmpData = self.fourierLocal(tmpData, 0)
        self.showFid(tmpData)
        
    def getdcOffset(self, pos1, pos2):
        minPos = int(min(pos1, pos2))
        maxPos = int(max(pos1, pos2))
        if minPos != maxPos:
            return np.mean(self.data1D[(len(self.data1D.shape)-1)*(slice(None), )+(slice(minPos, maxPos), )])
        else:
            return 0
            
    def dcOffset(self, offset):
        self.showFid(self.data1D-offset)

    def applyBaseline(self, degree, removeList, select=False):
        import numpy.polynomial.polynomial as poly
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        if len(self.data1D.shape) > 1:
            tmpData = self.data1D[0]
        else:
            tmpData = self.data1D
        tmpAx = np.arange(self.data1D.shape[-1])
        bArray = np.array([True]*self.data1D.shape[-1])
        for i in range(int(np.floor(len(removeList)/2.0))):
            minVal = min(removeList[2*i], removeList[2*i+1])
            maxVal = max(removeList[2*i], removeList[2*i+1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        polyCoeff = poly.polyfit(self.xax[bArray], tmpData[bArray], degree)
        y = poly.polyval(self.xax, polyCoeff)
        self.root.addMacro(['baselineCorrection', (list(np.real(y)), self.axes-self.data.dim, list(np.imag(y)), str(selectSlice))])
        return self.data.baselineCorrection(y, self.axes, select=selectSlice)
    
    def previewBaseline(self, degree, removeList):
        import numpy.polynomial.polynomial as poly
        if len(self.data1D.shape) > 1:
            tmpData = self.data1D[0]
        else:
            tmpData = self.data1D
        tmpAx = np.arange(self.data1D.shape[-1])
        bArray = np.array([True]*self.data1D.shape[-1])
        for i in range(int(np.floor(len(removeList)/2.0))):
            minVal = min(removeList[2*i], removeList[2*i+1])
            maxVal = max(removeList[2*i], removeList[2*i+1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        polyCoeff = poly.polyfit(self.xax[bArray], tmpData[bArray], degree)
        y = poly.polyval(self.xax, polyCoeff)
        if (self.plotType == 0):
            y = np.real(y)
        elif (self.plotType == 1):
            y = np.imag(y)
        elif (self.plotType == 2):
            y = np.real(y)
        elif (self.plotType == 3):
            y = np.abs(y)
        self.resetPreviewRemoveList()
        if len(self.data1D.shape) > 1:
            self.showFid(self.data1D, [self.xax], [y]*self.data1D.shape[0], ['g'])
        else:
            self.showFid(self.data1D, [self.xax], [y], ['g'])
        self.previewRemoveList(removeList)
    
    def previewRemoveList(self, removeList):
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.resetPreviewRemoveList()
        self.removeListLines = []
        for i in range(int(np.floor(len(removeList)/2.0))):
            self.removeListLines.append(self.ax.axvspan(self.xax[removeList[2*i]]*axMult, self.xax[removeList[2*i+1]]*axMult, color='r'))
        if len(removeList)%2:
            self.removeListLines.append(self.ax.axvline(self.xax[removeList[-1]]*axMult, c='r', linestyle='--'))
        self.canvas.draw()

    def resetPreviewRemoveList(self):
        if hasattr(self, 'removeListLines'):
            for i in self.removeListLines:
                i.remove()
            del self.removeListLines
    
    def states(self):
        returnValue = self.data.states(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['states', (self.axes-self.data.dim, )])
        return returnValue
    
    def statesTPPI(self):
        returnValue = self.data.statesTPPI(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['statesTPPI', (self.axes-self.data.dim, )])
        return returnValue

    def echoAntiEcho(self):
        returnValue = self.data.echoAntiEcho(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['echoAntiEcho', (self.axes-self.data.dim, )])
        return returnValue
    
    def integrate(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 0)
        else:
            self.root.addMacro(['integrate', (list(pos1), list(pos2), self.axes-self.data.dim, )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 0)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def sum(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 5)
        else:
            self.root.addMacro(['sum', (list(pos1), list(pos2), self.axes-self.data.dim, )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 5)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def maxMatrix(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 1)
        else:
            self.root.addMacro(['max', (list(pos1), list(pos2), self.axes-self.data.dim, )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 1)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue
    
    def minMatrix(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 2)
        else:
            self.root.addMacro(['min', (list(pos1), list(pos2), self.axes-self.data.dim, )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 2)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue
    
    def argmaxMatrix(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 3)
        else:
            self.root.addMacro(['argmax', (list(pos1), list(pos2), self.axes-self.data.dim, )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 3)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def argminMatrix(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 4)
        else:
            self.root.addMacro(['argmin', (list(pos1), list(pos2), self.axes-self.data.dim, )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 4)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue
    
    def average(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.matrixManipNew(pos1, pos2, self.axes, 6)
        else:
            self.root.addMacro(['average', (list(pos1), list(pos2), self.axes-self.data.dim, )])
            returnValue = self.data.matrixManip(pos1, pos2, self.axes, 6)
            if self.upd():
                self.plotReset()
                self.showFid()
            return returnValue

    def flipLR(self):
        returnValue = self.data.flipLR(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['fliplr', (self.axes-self.data.dim, )])
        return returnValue

    def concatenate(self, axes):
        returnValue = self.data.concatenate(axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['concatenate', (axes-self.data.dim+1, )])
        return returnValue

    def split(self, sections):
        returnValue = self.data.split(sections, self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['split', (sections, self.axes-self.data.dim-1)])
        return returnValue

    def diff(self):
        returnValue = self.data.diff(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['diff', (self.axes-self.data.dim)])
        return returnValue

    def cumsum(self):
        returnValue = self.data.cumsum(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['cumsum', (self.axes-self.data.dim)])
        return returnValue
    
    def insert(self, data, pos):
        returnValue = self.data.insert(data, pos, self.axes)
        self.upd()
        self.plotReset()
        self.showFid()
        self.root.addMacro(['insert', (np.real(data).tolist(), pos, self.axes-self.data.dim, np.imag(data).tolist())])
        return returnValue
    
    def delete(self, pos):
        returnValue = self.data.remove(pos, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['delete', (pos, self.axes-self.data.dim)])
        return returnValue

    def deletePreview(self, pos):
        self.data1D = np.delete(self.data1D, pos, axis=len(self.data1D.shape)-1)
        self.xax = np.delete(self.xax, pos)
        if (np.array(self.data1D.shape) != 0).all():
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
    
    def multiply(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        returnValue = self.data.multiply(data, self.axes, select=selectSlice)
        self.upd()
        self.showFid()
        self.root.addMacro(['multiply', (np.real(data).tolist(), self.axes-self.data.dim, np.imag(data).tolist(), str(selectSlice))])
        return returnValue
    
    def multiplyPreview(self, data):
        self.upd()
        self.data1D = self.data1D*data
        self.showFid()

    def subtractAvg(self, pos1, pos2):
        returnValue = self.data.subtractAvg(pos1, pos2, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['subtractAvg', (pos1, pos2, self.axes-self.data.dim)])
        return returnValue

    def subtractAvgPreview(self, pos1, pos2):
        self.upd()
        axes = self.data1D.ndim - 1
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        slicing = (slice(None), ) * axes + (slice(minPos, maxPos), ) + (slice(None), )*(self.data1D.ndim-1-axes)
        self.data1D -= np.mean(self.data1D[slicing], axis=-1, keepdims=True)
        self.showFid()
        
    def getRegion(self, pos1, pos2, newSpec=False):
        if newSpec:
            return self.data.getRegionNew(pos1, pos2, self.axes)
        else:
            returnValue = self.data.getRegion(pos1, pos2, self.axes)
            self.upd()
            self.plotReset()
            self.showFid()
            self.root.addMacro(['extract', (pos1, pos2, self.axes-self.data.dim)])
            return returnValue

    def fiddle(self, pos1, pos2, lb):
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        refSpec = np.zeros(self.data1D.shape[-1])
        if len(self.data1D.shape) > 1:
            refSpec[minPos:maxPos] = np.real(self.data1D[0][minPos:maxPos])
        else:
            refSpec[minPos:maxPos] = np.real(self.data1D[minPos:maxPos])
        returnValue = self.data.fiddle(refSpec,  lb,  self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['FIDDLE', (refSpec,  lb,  self.axes-self.data.dim)])
        return returnValue
        
    def shearing(self, shear, axes, axes2):
        returnValue = self.data.shear(shear, axes, axes2)
        self.upd()
        self.showFid()
        self.root.addMacro(['shear', (shear, axes-self.data.dim, axes2-self.data.dim)])
        return returnValue

    def reorder(self, pos, newLength):
        returnValue = self.data.reorder(pos, newLength, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['reorder', (pos, newLength, self.axes-self.data.dim)])
        return returnValue

    def ffm(self, posList, typeVal):
        returnValue = self.data.ffm_1d(posList, typeVal, self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['ffm', (posList, typeVal, self.axes-self.data.dim)])
        return returnValue
    
    def clean(self, posList, typeVal, gamma, threshold, maxIter):
        returnValue = self.data.clean(posList, typeVal, self.axes, gamma, threshold, maxIter)
        self.upd()
        self.showFid()
        self.root.addMacro(['clean', (posList, typeVal, self.axes-self.data.dim, gamma, threshold, maxIter)])
        return returnValue
    
    def ACMEentropy(self, phaseIn, phaseAll=True):
        if len(self.data1D.shape) > 1:
            tmp = self.data1D[0]
        else:
            tmp = self.data1D
        phase0 = phaseIn[0]
        if phaseAll:
            phase1 = phaseIn[1]
        else:
            phase1 = 0.0
        L = len(tmp)
        if self.spec>0:
            x = np.fft.fftshift(np.fft.fftfreq(L, 1.0/self.sw))/self.sw
            s0 = tmp*np.exp(1j*(phase0+phase1*x))
        else:
            s0 = np.fft.fftshift(np.fft.fft(tmp))*np.exp(1j*(phase0+phase1*x))
        s2 = np.real(s0)
        ds1 = np.abs((s2[3:L]-s2[1:L-2])/2.0)
        p1 = ds1/sum(ds1)
        p1[np.where(p1 == 0)] = 1
        h1  = -p1*np.log(p1)
        H1  = sum(h1)
        Pfun = 0.0
        as1 = s2 - np.abs(s2)
        sumas   = sum(as1)
        if (np.real(sumas) < 0): 
            Pfun = Pfun + sum(as1**2)/4/L**2
        return H1+1000*Pfun 

    def autoPhase(self, phaseNum):
        self.upd()
        if phaseNum == 0:
            phases = scipy.optimize.minimize(self.ACMEentropy, [0], (False, ), method='Powell')
            phases = [phases['x']]
        elif phaseNum == 1:
            phases = scipy.optimize.minimize(self.ACMEentropy, [0, 0], method='Powell')
            phases = phases['x']
        return phases

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
        self.root.addMacro(['setxax', (xax, self.axes-self.data.dim)])
        return returnVal

    def setAxType(self, val):
        if self.spec == 1:
            if self.ppm:
                oldAxMult = 1e6/self.ref
            else:
                oldAxMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            oldAxMult = 1000.0**self.axType
        if self.spec == 1:
            if val == 3:
                newAxMult = 1e6/self.ref
            else:
                newAxMult = 1.0/(1000.0**val)
        elif self.spec == 0:
            newAxMult = 1000.0**val 
        if val == 3:
            self.ppm = True
        else:
            self.ppm = False
            self.axType = val
        self.xminlim = self.xminlim * newAxMult / oldAxMult
        self.xmaxlim = self.xmaxlim * newAxMult / oldAxMult
        self.showFid()

    def hilbert(self):
        returnValue = self.data.hilbert(self.axes)
        self.upd()
        self.showFid()
        self.root.addMacro(['hilbert', (self.axes-self.data.dim, )])
        return returnValue
    
    def getDisplayedData(self):
        if len(self.data1D.shape) > 1:
            tmp = self.data1D[0]
        else:
            tmp = self.data1D
        if self.plotType == 0:
            return np.real(tmp)
        elif self.plotType == 1:
            return np.imag(tmp)
        elif self.plotType == 2:
            return np.real(tmp)
        elif self.plotType == 3:
            return np.abs(tmp)      

    def getColorMap(self):
        return COLORMAPLIST.index(self.colorMap)
        
    def setColorMap(self, num):
        self.colorMap = COLORMAPLIST[num]

    def setColor(self, color):
        self.color = color

    def setLw(self, lw):
        self.linewidth = lw

    def setContourColors(self, colors):
        self.contourColors = colors

    def setContourConst(self, constant):
        self.contourConst = constant
        
    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False, output=None): #display the 1D data
        self.peakPickReset()
        if tmpdata is None:
            tmpdata = self.data1D
        self.ax.cla()
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.line_xdata = self.xax*axMult
        if old:
            if (self.plotType == 0):
                if self.single:
                    self.ax.scatter(self.line_xdata, np.real(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
                else:
                    self.ax.plot(self.line_xdata, np.real(self.data1D), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 1):
                if self.single:
                    self.ax.scatter(self.line_xdata, np.imag(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
                else:
                    self.ax.plot(self.line_xdata, np.imag(self.data1D), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 2):
                if self.single:
                    self.ax.scatter(self.line_xdata, np.real(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
                else:
                    self.ax.plot(self.line_xdata, np.real(self.data1D), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 3):
                if self.single:
                    self.ax.scatter(self.line_xdata, np.abs(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
                else:
                    self.ax.plot(self.line_xdata, np.abs(self.data1D), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
        if (extraX is not None):
            for num in range(len(extraX)):
                if self.single:
                    self.ax.scatter(extraX[num]*axMult, extraY[num], c=extraColor[num])
                else:
                    self.ax.plot(extraX[num]*axMult, extraY[num], c=extraColor[num], linewidth=self.linewidth)
        if (self.plotType == 0):
            self.line_ydata = np.real(tmpdata)
            if self.single:
                self.ax.scatter(self.line_xdata, np.real(tmpdata), c=self.color, label=self.data.name)
            else:
                self.ax.plot(self.line_xdata, np.real(tmpdata), c=self.color, linewidth=self.linewidth, label=self.data.name)
        elif(self.plotType == 1):
            self.line_ydata = np.imag(tmpdata)
            if self.single:
                self.ax.scatter(self.line_xdata, np.imag(tmpdata), c=self.color, label=self.data.name)
            else:
                self.ax.plot(self.line_xdata, np.imag(tmpdata), c=self.color, linewidth=self.linewidth, label=self.data.name)
        elif(self.plotType == 2):
            self.line_ydata = np.real(tmpdata)
            if self.single:
                self.ax.scatter(self.line_xdata, np.imag(tmpdata), c='r', label=self.data.name+'_imag')
                self.ax.scatter(self.line_xdata, np.real(tmpdata), c=self.color, label=self.data.name)
            else:
                self.ax.plot(self.line_xdata, np.imag(tmpdata), c='r', linewidth=self.linewidth, label=self.data.name+'_imag')
                self.ax.plot(self.line_xdata, np.real(tmpdata), c=self.color, linewidth=self.linewidth, label=self.data.name)
        elif(self.plotType == 3):
            self.line_ydata = np.abs(tmpdata)
            if self.single:
                self.ax.scatter(self.line_xdata, np.abs(tmpdata), c=self.color, label=self.data.name)
            else:
                self.ax.plot(self.line_xdata, np.abs(tmpdata), c=self.color, linewidth=self.linewidth, label=self.data.name)
        if self.spec == 0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(u'Time [\u03BCs]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec ==1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.grids[0])
        self.ax.yaxis.grid(self.grids[1])
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        if output is not None:
            self.canvas.print_figure(output)
        self.canvas.draw()

    def plotReset(self, xReset=True, yReset=True): #set the plot limits to min and max values
        if self.plotType == 0:
            miny = min(np.real(self.data1D))
            maxy = max(np.real(self.data1D))
        elif self.plotType == 1:
            miny = min(np.imag(self.data1D))
            maxy = max(np.imag(self.data1D))
        elif self.plotType == 2:
            miny = min(min(np.real(self.data1D)), min(np.imag(self.data1D)))
            maxy = max(max(np.real(self.data1D)), max(np.imag(self.data1D)))
        elif self.plotType == 3:
            miny = min(np.abs(self.data1D))
            maxy = max(np.abs(self.data1D))
        else:
            miny = -1
            maxy = 1
        differ = 0.05*(maxy-miny) #amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny-differ
            self.ymaxlim = maxy+differ
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if xReset:
            self.xminlim = min(self.xax*axMult)
            self.xmaxlim = max(self.xax*axMult)
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)


#########################################################################################################
class CurrentScatter(Current1D):

    X_RESIZE = False
    Y_RESIZE = False
    
    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        Current1D.__init__(self, root, fig, canvas, data, duplicateCurrent)

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False, output=None): #display the 1D data
        self.peakPickReset()
        if tmpdata is None:
            tmpdata = self.data1D
        self.ax.cla()
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.line_xdata = self.xax*axMult
        if old:
            if (self.plotType == 0):
                self.ax.scatter(self.line_xdata, np.real(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
            elif(self.plotType == 1):
                self.ax.scatter(self.line_xdata, np.imag(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
            elif(self.plotType == 2):
                self.ax.scatter(self.line_xdata, np.real(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
            elif(self.plotType == 3):
                self.ax.scatter(self.line_xdata, np.abs(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
        if (extraX is not None):
            for num in range(len(extraX)):
                self.ax.scatter(extraX[num]*axMult, extraY[num], c=extraColor[num])
        if (self.plotType == 0):
            self.line_ydata = np.real(tmpdata)
            self.ax.scatter(self.line_xdata, np.real(tmpdata), c=self.color, label=self.data.name)
        elif(self.plotType == 1):
            self.line_ydata = np.imag(tmpdata)
            self.ax.scatter(self.line_xdata, np.imag(tmpdata), c=self.color, label=self.data.name)
        elif(self.plotType == 2):
            self.line_ydata = np.real(tmpdata)
            self.ax.scatter(self.line_xdata, np.imag(tmpdata), c='r', label=self.data.name+'_imag')
            self.ax.scatter(self.line_xdata, np.real(tmpdata), c=self.color, label=self.data.name)
        elif(self.plotType == 3):
            self.line_ydata = np.abs(tmpdata)
            self.ax.scatter(self.line_xdata, np.abs(tmpdata), c=self.color, label=self.data.name)
        if self.spec == 0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(u'Time [\u03BCs]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.grids[0])
        self.ax.yaxis.grid(self.grids[1])
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        if output is not None:
            self.canvas.print_figure(output)
        self.canvas.draw()

#########################################################################################################
#the class from which the data of multiple spectra is displayed, the operations which only edit the content of this class are for previewing
class CurrentMulti(Current1D):

    X_RESIZE = False
    Y_RESIZE = True
    
    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        if hasattr(duplicateCurrent, 'extraData'):
            self.extraData = duplicateCurrent.extraData
        else:
            self.extraData = []
        if hasattr(duplicateCurrent, 'extraLoc'):
            self.extraLoc = duplicateCurrent.extraLoc
        else:
            self.extraLoc = []
        if hasattr(duplicateCurrent, 'extraColor'):
            self.extraColor = duplicateCurrent.extraColor
        else:
            self.extraColor = []
        if hasattr(duplicateCurrent, 'extraName'):
            self.extraName = duplicateCurrent.extraName
        else:
            self.extraName = []
        if hasattr(duplicateCurrent, 'extraAxes'):
            self.extraAxes = duplicateCurrent.extraAxes
        else:
            self.extraAxes = [] 
        Current1D.__init__(self, root, fig, canvas, data, duplicateCurrent)

    def setExtraSlice(self, extraNum, axes, locList): #change the slice
        self.extraAxes[extraNum] = axes
        self.extraLoc[extraNum] = locList
        #self.showFid()
        
    def copyCurrent(self, root, fig, canvas, data):
        return CurrentMulti(root, fig, canvas, data, self)
    
    def addExtraData(self, data, name):
        self.extraName.append(name)
        self.extraData.append(data)
        self.extraLoc.append([0]*(len(self.extraData[-1].data.shape)-1))
        self.extraColor.append((0, 0, 0, 1))                                        #find a good color system
        self.extraAxes.append(len(data.data.shape)-1)
        self.showFid()

    def delExtraData(self, num):
        del self.extraData[num]
        del self.extraLoc[num]
        del self.extraColor[num]
        del self.extraName[num]
        del self.extraAxes[num]
        self.showFid()

    def setExtraColor(self, num, color):
        self.extraColor[num] = color
        self.showFid()

    def getExtraColor(self, num):
        return tuple(np.array(255*np.array(self.extraColor[num]), dtype=int))
        
    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-1)
        self.resetExtraLocList()
        
    def resetExtraLocList(self, num=None):
        if num is None:
            for i in range(len(self.extraLoc)):
                self.extraLoc[i] = [0]*(len(self.extraData[i].data.shape)-1)
        else:
            self.extraLoc[num] = [0]*(len(self.extraData[num].data.shape)-1)
                
    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False, output=None): #display the 1D data
        self.peakPickReset()
        self.ax.cla()
        for i in range(len(self.extraData)):
            data = self.extraData[i]
            try:
                if self.extraData[i].dim <= self.extraAxes[i]:
                    self.extraAxes[i] = len(self.extraData[i].data.shape)-1
                    self.resetExtraLocList(i)
                updateVar = data.getSlice(self.extraAxes[i], self.extraLoc[i])
            except:
                self.resetExtraLocList(i)
                updateVar = data.getSlice(self.extraAxes[i], self.extraLoc[i])
            data1D = updateVar[0]
            spec = updateVar[3]
            xax = updateVar[5]
            ref = updateVar[6]
            if ref is None:
                ref = data.freq[self.extraAxes[i]]
            if spec == 1:
                if self.ppm:
                    axMult = 1e6/ref
                else:
                    axMult = 1.0/(1000.0**self.axType)
            elif spec == 0:
                axMult = 1000.0**self.axType
            line_xdata = xax*axMult
            if (self.plotType == 0):
                if len(data1D) == 1:
                    self.ax.scatter(line_xdata, np.real(data1D), c=self.extraColor[i], label=data.name)
                else:
                    self.ax.plot(line_xdata, np.real(data1D), c=self.extraColor[i], linewidth=self.linewidth, label=data.name)
            elif(self.plotType == 1):
                if len(data1D) == 1:
                    self.ax.scatter(line_xdata, np.imag(data1D), c=self.extraColor[i], label=data.name)
                else:
                    self.ax.plot(line_xdata, np.imag(data1D), c=self.extraColor[i], linewidth=self.linewidth, label=data.name)
            elif(self.plotType == 2):
                if len(data1D) == 1:
                    self.ax.scatter(line_xdata, np.real(data1D), c=self.extraColor[i], label=data.name)
                else:
                    self.ax.plot(line_xdata, np.real(data1D), c=self.extraColor[i], linewidth=self.linewidth, label=data.name)
            elif(self.plotType == 3):
                if len(data1D) == 1:
                    self.ax.scatter(line_xdata, np.abs(data1D), c=self.extraColor[i], label=data.name)
                else:
                    self.ax.plot(line_xdata, np.abs(data1D), c=self.extraColor[i], linewidth=self.linewidth, label=data.name)
        if tmpdata is None:
            tmpdata = self.data1D
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.line_xdata = self.xax*axMult
        if old:
            if (self.plotType == 0):
                if self.single:
                    self.ax.scatter(self.line_xdata, np.real(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
                else:
                    self.ax.plot(self.line_xdata, np.real(self.data1D), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 1):
                if self.single:
                    self.ax.scatter(self.line_xdata, np.imag(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
                else:
                    self.ax.plot(self.line_xdata, np.imag(self.data1D), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 2):
                if self.single:
                    self.ax.scatter(self.line_xdata, np.real(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
                else:
                    self.ax.plot(self.line_xdata, np.real(self.data1D), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 3):
                if self.single:
                    self.ax.scatter(self.line_xdata, np.abs(self.data1D), c='k', alpha=0.2, label=self.data.name+'_old')
                else:
                    self.ax.plot(self.line_xdata, np.abs(self.data1D), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
        if (extraX is not None):
            for num in range(len(extraX)):
                if self.single:
                    self.ax.scatter(extraX[num]*axMult, extraY[num], c=extraColor[num])
                else:
                    self.ax.plot(extraX[num]*axMult, extraY[num], linewidth=self.linewidth, c=extraColor[num])
        if (self.plotType == 0):
            self.line_ydata = np.real(tmpdata)
            if self.single:
                self.ax.scatter(self.line_xdata, np.real(tmpdata), c=self.color, label=self.data.name)
            else:
                self.ax.plot(self.line_xdata, np.real(tmpdata), c=self.color, linewidth=self.linewidth, label=self.data.name)
        elif(self.plotType == 1):
            self.line_ydata = np.imag(tmpdata)
            if self.single:
                self.ax.scatter(self.line_xdata, np.imag(tmpdata), c=self.color, label=self.data.name)
            else:
                self.ax.plot(self.line_xdata, np.imag(tmpdata), c=self.color, linewidth=self.linewidth, label=self.data.name)
        elif(self.plotType == 2):
            self.line_ydata = np.real(tmpdata)
            if self.single:
                self.ax.scatter(self.line_xdata, np.imag(tmpdata), c='r', label=self.data.name+'_imag')
                self.ax.scatter(self.line_xdata, np.real(tmpdata), c=self.color, label=self.data.name)
            else:
                self.ax.plot(self.line_xdata, np.imag(tmpdata), c='r', linewidth=self.linewidth, label=self.data.name+'_imag')
                self.ax.plot(self.line_xdata, np.real(tmpdata), c=self.color, linewidth=self.linewidth, label=self.data.name)
        elif(self.plotType == 3):
            self.line_ydata = np.abs(tmpdata)
            if self.single:
                self.ax.scatter(self.line_xdata, np.abs(tmpdata), c=self.color, label=self.data.name)
            else:
                self.ax.plot(self.line_xdata, np.abs(tmpdata), c=self.color, linewidth=self.linewidth, label=self.data.name)
        if self.spec == 0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(u'Time [\u03BCs]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.grids[0])
        self.ax.yaxis.grid(self.grids[1])
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        if output is not None:
            self.canvas.print_figure(output)
        self.canvas.draw()
        
#########################################################################################################
#the class from which the stacked data is displayed, the operations which only edit the content of this class are for previewing
class CurrentStacked(Current1D):

    X_RESIZE = False
    Y_RESIZE = True
    
    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        self.data = data
        if hasattr(duplicateCurrent, 'axes2'):
            self.axes2 = duplicateCurrent.axes2
        else:
            self.axes2 = len(self.data.data.shape)-2
            if hasattr(duplicateCurrent, 'axes'):
                if self.axes2 == duplicateCurrent.axes:
                    self.axes2 = (self.axes2-1) % self.data.dim
        if hasattr(duplicateCurrent, 'stackBegin'):
            self.stackBegin = duplicateCurrent.stackBegin
        else:
            self.stackBegin = None
        if hasattr(duplicateCurrent, 'stackEnd'):
            self.stackEnd = duplicateCurrent.stackEnd
        else:
            self.stackEnd = None
        if hasattr(duplicateCurrent, 'stackStep'):
            self.stackStep = duplicateCurrent.stackStep
        else:
            self.stackStep = None
            if self.data.data.shape[self.axes2] > 100:
                self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
        self.spacing = 0
        Current1D.__init__(self, root, fig, canvas, data, duplicateCurrent)
        #self.startUp()

    def startUp(self, xReset=True, yReset=True):
        self.resetSpacing()
        self.plotReset(xReset, yReset)
        self.showFid()
        
    def copyCurrent(self, root, fig, canvas, data):
        return CurrentStacked(root, fig, canvas, data, self)
        
    def upd(self): #get new data from the data instance
        if self.data.dim < 2:
            self.root.rescue()
            return False
        if (len(self.locList)+2) != self.data.dim:
            self.resetLocList()
        try:
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList, self.stackBegin, self.stackEnd, self.stackStep)
        except:
            self.resetLocList()
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList, self.stackBegin, self.stackEnd, self.stackStep)
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
        self.single = self.data1D.shape[-1] == 1
        return True

    def setBlock(self, axes, axes2, locList, stackBegin=None, stackEnd=None, stackStep=None): #change the slice
        axesSame = True
        if (self.axes != axes) or (self.axes2 != axes2) :
            axesSame = False
        self.axes = axes
        self.axes2 = axes2
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        if stackStep == None:
            if self.data.data.shape[self.axes2] > 100:
                self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
        self.locList = locList
        self.upd()
        if not axesSame:
            self.resetSpacing()
            self.plotReset()
        self.showFid()

    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-2)
        
    def stackSelect(self, stackBegin, stackEnd, stackStep):
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        self.upd()
        self.plotReset(False, True)
        self.showFid()

    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=None): #display the 1D data including the apodization function
        t = np.arange(0, len(self.data1D[0]))/(self.sw)
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                self.dispMsg('shiftingAxes cannot be equal to axes')
            elif shiftingAxes == self.axes2:
                ar = np.arange(self.data.data.shape[self.axes2])[slice(self.stackBegin, self.stackEnd, self.stackStep)]
                x = np.ones((len(ar), len(self.data1D[0])))
                for i in range(len(ar)):
                    shift1 = shift + shifting*ar[i]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                    t2 = t - shift1
                    x2 = np.ones(len(self.data1D[0]))
                    if lor is not None:
                        x2 = x2*np.exp(-np.pi*lor*abs(t2))
                    if gauss is not None:
                        x2 = x2*np.exp(-((np.pi*gauss*t2)**2)/(4*np.log(2)))
                    if cos2 is not None:
                        x2 = x2*(np.cos(cos2*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
                    if hamming is not None:
                        alpha = 0.53836 # constant for hamming window
                        x2 = x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, np.pi, len(self.data1D[0])))))
                    if self.wholeEcho:
                        x2[-1:-(len(x2)/2+1):-1] = x2[:len(x2)/2]
                    x[i] = x2
            else:
                if (shiftingAxes < self.axes) and (shiftingAxes < self.axes2):
                    shift += shifting*self.locList[shiftingAxes]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                elif (shiftingAxes > self.axes) and (shiftingAxes > self.axes2):
                    shift += shifting*self.locList[shiftingAxes-2]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                else:
                    shift += shifting*self.locList[shiftingAxes-1]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                t2 = t - shift
                x = np.ones(len(self.data1D[0]))
                if lor is not None:
                    x = x*np.exp(-lor*abs(t2))
                if gauss is not None:
                    x = x*np.exp(-(gauss*t2)**2)
                if cos2 is not None:
                    x = x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
                if hamming is not None:
                    alpha = 0.53836 # constant for hamming window
                    x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, np.pi, len(self.data1D[0])))))
                if self.wholeEcho:
                    x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
                x = np.repeat([x], len(self.data1D), axis=0)
        else:
            t2 = t - shift
            x = np.ones(len(self.data1D[0]))
            if lor is not None:
                x = x*np.exp(-lor*abs(t2))
            if gauss is not None:
                x = x*np.exp(-(gauss*t2)**2)
            if cos2 is not None:
                x = x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
            if hamming is not None:
                alpha = 0.53836 # constant for hamming window
                x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, np.pi, len(self.data1D[0])))))
            if self.wholeEcho:
                x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
            x = np.repeat([x], len(self.data1D), axis=0)
        y = self.data1D
        self.ax.cla()
        if self.spec == 1:
            y = np.fft.ifftn(np.fft.ifftshift(y, axes=1), axes=[1])
            y = y*x
            y = np.fft.fftshift(np.fft.fftn(y, axes=[1]), axes=1)
        else:
            y = y*x
        if self.spec == 0:
            if self.plotType == 0:
                self.showFid(y, [t], x*np.amax(np.real(self.data1D)), ['g'], old=True)
            elif self.plotType == 1:
                self.showFid(y, [t], x*np.amax(np.imag(self.data1D)), ['g'], old=True)
            elif self.plotType == 2:
                self.showFid(y, [t], x*np.amax(np.amax(np.real(self.data1D)), np.amax(np.imag(self.data1D))), ['g'], old=True)
            elif self.plotType == 3:
                self.showFid(y, [t], x*np.amax(np.abs(self.data1D)), ['g'], old=True)
        else:
            self.showFid(y)

    def setSpacing(self, spacing):
        self.spacing = spacing
        self.plotReset(False, True)
        self.showFid()

    def resetSpacing(self):
        difference = np.diff(self.data1D, axis=0)
        if difference.size == 0:
            self.spacing = 0
        else:
            if self.plotType == 0:
                difference = np.amin(np.real(difference))
                amp = np.amax(np.real(self.data1D))-np.amin(np.real(self.data1D))
            elif self.plotType == 1:
                difference = np.amin(np.imag(difference))
                amp = np.amax(np.imag(self.data1D))-np.amin(np.imag(self.data1D))
            elif self.plotType == 2:
                difference = np.amin((np.real(difference), np.imag(difference)))
                amp = np.amax((np.real(self.data1D), np.imag(self.data1D)))-np.amin((np.real(self.data1D), np.imag(self.data1D)))
            elif self.plotType == 3:
                difference = np.amin(np.abs(difference))
                amp = np.amax(np.abs(self.data1D))-np.amin(np.abs(self.data1D))
            self.spacing = np.abs(difference) + 0.1*amp   
            
    def altScroll(self, event):
        self.spacing = self.spacing*1.1**event.step
        self.root.sideframe.scrollSpacing(self.spacing)
        self.showFid()

    def altReset(self):
        self.resetSpacing()
        self.root.sideframe.scrollSpacing(self.spacing)
        self.showFid()
        
    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False): #display the 1D data
        self.peakPickReset()
        if tmpdata is None:
            tmpdata = self.data1D
        self.ax.cla()
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.line_xdata = self.xax*axMult
        if old:
            if (self.plotType == 0):
                for num in range(len(self.data1D)):
                    if self.single:
                        self.ax.scatter(self.line_xdata, num*self.spacing+np.real(self.data1D[num]), c='k', alpha=0.2, label=self.data.name+'_old')
                    else:
                        self.ax.plot(self.line_xdata, num*self.spacing+np.real(self.data1D[num]), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 1):
                for num in range(len(self.data1D)):
                    if self.single:
                        self.ax.scatter(self.line_xdata, num*self.spacing+np.imag(self.data1D[num]), c='k', alpha=0.2, label=self.data.name+'_old')
                    else:
                        self.ax.plot(self.line_xdata, num*self.spacing+np.imag(self.data1D[num]), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 2):
                for num in range(len(self.data1D)):
                    if self.single:
                        self.ax.scatter(self.line_xdata, num*self.spacing+np.real(self.data1D[num]), c='k', alpha=0.2, label=self.data.name+'_old')
                    else:
                        self.ax.plot(self.line_xdata, num*self.spacing+np.real(self.data1D[num]), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 3):
                for num in range(len(self.data1D)):
                    if self.single:
                        self.ax.scatter(self.line_xdata, num*self.spacing+np.abs(self.data1D[num]), c='k', alpha=0.2, label=self.data.name+'_old')
                    else:
                        self.ax.plot(self.line_xdata, num*self.spacing+np.abs(self.data1D[num]), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
        if (extraX is not None):
            for num in range(len(extraY)):
                if self.single:
                    self.ax.scatter(extraX[0]*axMult, num*self.spacing+extraY[num], c=extraColor[0])
                else:
                    self.ax.plot(extraX[0]*axMult, num*self.spacing+extraY[num], linewidth=self.linewidth, c=extraColor[0])
        if (self.plotType == 0):
            tmpdata = np.real(tmpdata)
        elif (self.plotType == 1):
            tmpdata = np.imag(tmpdata)
        elif(self.plotType == 3):
            tmpdata = np.abs(tmpdata)
        self.line_ydata = np.real(tmpdata[0])
        if self.single:
            for num in range(len(tmpdata)):
                if (self.plotType == 2):
                    self.ax.scatter(self.line_xdata, num*self.spacing+np.imag(tmpdata[num]), c='r', label=self.data.name+'_imag')
                self.ax.scatter(self.line_xdata, num*self.spacing+np.real(tmpdata[num]), c=self.color, label=self.data.name)
        else:
            for num in range(len(tmpdata)):
                if (self.plotType == 2):
                    self.ax.plot(self.line_xdata, num*self.spacing+np.imag(tmpdata[num]), c='r', linewidth=self.linewidth, label=self.data.name+'_imag')
                self.ax.plot(self.line_xdata, num*self.spacing+np.real(tmpdata[num]), c=self.color, linewidth=self.linewidth, label=self.data.name)  
        if self.spec == 0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(u'Time [\u03BCs]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.grids[0])
        self.ax.yaxis.grid(self.grids[1])
        self.canvas.draw()

    def plotReset(self, xReset=True, yReset=True): #set the plot limits to min and max values
        self.ax = self.fig.gca()
        incr = np.repeat(np.arange(len(self.data1D)).reshape((len(self.data1D), 1)), len(self.data1D[0]), axis=1)*self.spacing
        if self.plotType == 0:
            miny = np.amin(np.real(self.data1D)+incr)
            maxy = np.amax(np.real(self.data1D)+incr)
        elif self.plotType == 1:
            miny = np.amin(np.imag(self.data1D)+incr)
            maxy = np.amax(np.imag(self.data1D)+incr)
        elif self.plotType == 2:
            miny = np.amin((np.amin(np.real(self.data1D)+incr), np.amin(np.imag(self.data1D)+incr)))
            maxy = np.amax((np.amax(np.real(self.data1D)+incr), np.amax(np.imag(self.data1D)+incr)))
        elif self.plotType == 3:
            miny = np.amin(np.abs(self.data1D)+incr)
            maxy = np.amax(np.abs(self.data1D)+incr)
        else:
            miny = -1
            maxy = 1
        differ = 0.05*(maxy-miny) #amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny-differ
            self.ymaxlim = maxy+differ
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if xReset:
            self.xminlim = min(self.xax*axMult)
            self.xmaxlim = max(self.xax*axMult)
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)

#########################################################################################################
#the class from which the arrayed data is displayed, the operations which only edit the content of this class are for previewing
class CurrentArrayed(Current1D):

    X_RESIZE = True
    Y_RESIZE = False
    
    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        self.data = data
        if hasattr(duplicateCurrent, 'axes2'):
            self.axes2 = duplicateCurrent.axes2
        else:
            self.axes2 = len(self.data.data.shape)-2
            if hasattr(duplicateCurrent, 'axes'):
                if self.axes2 == duplicateCurrent.axes:
                    self.axes2 = (self.axes2-1) % self.data.dim
        if hasattr(duplicateCurrent, 'stackBegin'):
            self.stackBegin = duplicateCurrent.stackBegin
        else:
            self.stackBegin = None
        if hasattr(duplicateCurrent, 'stackEnd'):
            self.stackEnd = duplicateCurrent.stackEnd
        else:
            self.stackEnd = None
        if hasattr(duplicateCurrent, 'stackStep'):
            self.stackStep = duplicateCurrent.stackStep
        else:
            self.stackStep = None
            if self.data.data.shape[self.axes2] > 100:
                self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
        self.spacing = 0
        if duplicateCurrent is not None:
            if isinstance(duplicateCurrent, CurrentArrayed):
                self.zminlim = duplicateCurrent.zminlim
                self.zmaxlim = duplicateCurrent.zmaxlim
            else:
                #The z-axes limits are in xax units unlike the x-axes and y-axes limits
                if duplicateCurrent.spec == 1:
                    if duplicateCurrent.ppm:
                        axMult = 1e6/duplicateCurrent.ref
                    else:
                        axMult = 1.0/(1000.0**duplicateCurrent.axType)
                elif duplicateCurrent.spec == 0:
                    axMult = 1000.0**duplicateCurrent.axType
                self.zminlim = (duplicateCurrent.xminlim)/axMult
                self.zmaxlim = (duplicateCurrent.xmaxlim)/axMult
        Current1D.__init__(self, root, fig, canvas, data, duplicateCurrent)

    def startUp(self, xReset=True, yReset=True):
        self.resetSpacing(False)
        self.plotReset(xReset, yReset)
        self.showFid()

    def copyCurrent(self, root, fig, canvas, data):
        return CurrentArrayed(root, fig, canvas, data, self)
        
    def upd(self): #get new data from the data instance
        if self.data.dim < 2:
            self.root.rescue()
            return False
        if (len(self.locList)+2) != self.data.dim:
            self.resetLocList()
        try:
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList, self.stackBegin, self.stackEnd, self.stackStep)
        except:
            self.resetLocList()
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList, self.stackBegin, self.stackEnd, self.stackStep)
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
        self.single = self.data1D.shape[-1] == 1
        return True
 
    def setBlock(self, axes, axes2, locList, stackBegin=None, stackEnd=None, stackStep=None): #change the slice
        axesSame = True
        if (self.axes != axes) or (self.axes2 != axes2) :
            axesSame = False
        self.axes = axes
        self.axes2 = axes2
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        if self.data.data.shape[self.axes2] > 100:
            self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
        self.locList = locList
        self.upd()
        if not axesSame:
            self.resetSpacing()
            self.plotReset()
        self.showFid()

    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-2)
        
    def stackSelect(self, stackBegin, stackEnd, stackStep):
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        self.upd()
        self.plotReset(True, False)
        self.showFid()

    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=None): #display the 1D data including the apodization function
        t = np.arange(0, len(self.data1D[0]))/(self.sw)
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                self.dispMsg('shiftingAxes cannot be equal to axes')
            elif shiftingAxes == self.axes2:
                ar = np.arange(self.data.data.shape[self.axes2])[slice(self.stackBegin, self.stackEnd, self.stackStep)]
                x = np.ones((len(ar), len(self.data1D[0])))
                for i in range(len(ar)):
                    shift1 = shift + shifting*ar[i]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                    t2 = t - shift1
                    x2 = np.ones(len(self.data1D[0]))
                    if lor is not None:
                        x2 = x2*np.exp(-np.pi*lor*abs(t2))
                    if gauss is not None:
                        x2 = x2*np.exp(-((np.pi*gauss*t2)**2)/(4*np.log(2)))
                    if cos2 is not None:
                        x2 = x2*(np.cos(cos2*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
                    if hamming is not None:
                        alpha = 0.53836 # constant for hamming window
                        x2 = x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, np.pi, len(self.data1D[0])))))
                    if self.wholeEcho:
                        x2[-1:-(len(x2)/2+1):-1] = x2[:len(x2)/2]
                    x[i] = x2
            else:
                if (shiftingAxes < self.axes) and (shiftingAxes < self.axes2):
                    shift += shifting*self.locList[shiftingAxes]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                elif (shiftingAxes > self.axes) and (shiftingAxes > self.axes2):
                    shift += shifting*self.locList[shiftingAxes-2]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                else:
                    shift += shifting*self.locList[shiftingAxes-1]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                t2 = t - shift
                x = np.ones(len(self.data1D[0]))
                if lor is not None:
                    x = x*np.exp(-lor*abs(t2))
                if gauss is not None:
                    x = x*np.exp(-(gauss*t2)**2)
                if cos2 is not None:
                    x = x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
                if hamming is not None:
                    alpha = 0.53836 # constant for hamming window
                    x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D)+np.linspace(0, np.pi, len(self.data1D)))))
                if self.wholeEcho:
                    x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
                x = np.repeat([x], len(self.data1D), axis=0)
        else:
            t2 = t - shift
            x = np.ones(len(self.data1D[0]))
            if lor is not None:
                x = x*np.exp(-lor*abs(t2))
            if gauss is not None:
                x = x*np.exp(-(gauss*t2)**2)
            if cos2 is not None:
                x = x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
            if hamming is not None:
                alpha = 0.53836 # constant for hamming window
                x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D)+np.linspace(0, np.pi, len(self.data1D)))))
            if self.wholeEcho:
                x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
            x = np.repeat([x], len(self.data1D), axis=0)
        y = self.data1D
        self.ax.cla()
        if self.spec == 1:
            y = np.fft.ifftn(np.fft.ifftshift(y, axes=1), axes=[1])
            y = y*x
            y = np.fft.fftshift(np.fft.fftn(y, axes=[1]), axes=1)
        else:
            y = y*x
        if self.spec == 0:
            if self.plotType == 0:
                self.showFid(y, [t], x*np.amax(np.real(self.data1D)), ['g'], old=True)
            elif self.plotType == 1:
                self.showFid(y, [t], x*np.amax(np.imag(self.data1D)), ['g'], old=True)
            elif self.plotType == 2:
                self.showFid(y, [t], x*np.amax(np.amax(np.real(self.data1D)), np.amax(np.imag(self.data1D))), ['g'], old=True)
            elif self.plotType == 3:
                self.showFid(y, [t], x*np.amax(np.abs(self.data1D)), ['g'], old=True)
        else:
            self.showFid(y)

    def setSpacing(self, spacing):
        self.spacing = spacing
        self.plotReset(True, False)
        self.showFid()

    def resetSpacing(self, zlims=True):
        if zlims:
            self.zminlim = min(self.xax)
            self.zmaxlim = max(self.xax)
        xaxZlims = (self.xax > self.zminlim) & (self.xax < self.zmaxlim)
        self.spacing = (self.xax[xaxZlims][-1]-self.xax[xaxZlims][0])*1.1

    def altScroll(self, event):
        self.spacing = self.spacing*1.1**event.step
        self.root.sideframe.scrollSpacing(self.spacing)
        self.showFid()

    def altReset(self):
        self.resetSpacing()
        self.root.sideframe.scrollSpacing(self.spacing)
        self.showFid()
        
    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False): #display the 1D data
        self.peakPickReset()
        if tmpdata is None:
            tmpdata = self.data1D
        self.ax.cla()
        if self.spec > 0:
            direc = slice(None, None, -1)
        else:
            direc = slice(None, None, 1) 
        xaxZlims = (self.xax > self.zminlim) & (self.xax < self.zmaxlim)
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if old:
            if (self.plotType == 0):
                oldData = np.real(self.data1D)
            elif(self.plotType == 1):
                oldData = np.imag(self.data1D)
            elif(self.plotType == 2):
                oldData = np.real(self.data1D)
            elif(self.plotType == 3):
                oldData = np.abs(self.data1D)
            for num in range(len(self.data1D)):
                if self.single:
                    self.ax.scatter((num*self.spacing+self.xax[xaxZlims])*axMult, oldData[num][xaxZlims][direc], c='k', alpha=0.2, label=self.data.name+'_old')
                else:
                    self.ax.plot((num*self.spacing+self.xax[xaxZlims])*axMult, oldData[num][xaxZlims][direc], c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
        if (extraX is not None):
            extraZlims = (extraX[0] > self.zminlim) & (extraX[0] < self.zmaxlim)
            for num in range(len(extraY)):
                if self.single:
                    self.ax.scatter((num*self.spacing+extraX[0][extraZlims])*axMult, extraY[num][extraZlims][direc], c=extraColor[0])
                else:
                    self.ax.plot((num*self.spacing+extraX[0][extraZlims])*axMult, extraY[num][extraZlims][direc], linewidth=self.linewidth, c=extraColor[0])
        if (self.plotType == 0):
            tmpdata = np.real(tmpdata)
        elif(self.plotType == 1):
            tmpdata = np.imag(tmpdata)
        elif(self.plotType == 3):
            tmpdata = np.abs(tmpdata)
        self.line_xdata = []
        self.line_ydata = []
        if self.single:
            for num in range(len(tmpdata)):
                if (self.plotType == 2):
                    self.ax.scatter((num*self.spacing+self.xax[xaxZlims])*axMult, np.imag(tmpdata[num][xaxZlims])[direc], c='r', label=self.data.name+'_imag')
                self.line_xdata = np.append(self.line_xdata, (num*self.spacing+self.xax[xaxZlims])*axMult)
                self.line_ydata = np.append(self.line_ydata, np.real(tmpdata[num][xaxZlims])[direc])
                self.ax.scatter((num*self.spacing+self.xax[xaxZlims])*axMult, np.real(tmpdata[num][xaxZlims])[direc], c=self.color, label=self.data.name)
        else:
            for num in range(len(tmpdata)):
                if (self.plotType == 2):
                    self.ax.plot((num*self.spacing+self.xax[xaxZlims])*axMult, np.imag(tmpdata[num])[direc], c='r', linewidth=self.linewidth, label=self.data.name+'_imag')
                self.line_xdata = np.append(self.line_xdata, (num*self.spacing+self.xax[xaxZlims])*axMult)
                self.line_ydata = np.append(self.line_ydata, np.real(tmpdata[num][xaxZlims])[direc])
                self.ax.plot((num*self.spacing+self.xax[xaxZlims])*axMult, np.real(tmpdata[num][xaxZlims])[direc], c=self.color, linewidth=self.linewidth, label=self.data.name)
        if self.spec == 0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(u'Time [\u03BCs]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.grids[0])
        self.ax.yaxis.grid(self.grids[1])
        self.canvas.draw()

    def plotReset(self, xReset=True, yReset=True): #set the plot limits to min and max values
        if self.plotType == 0:
            miny = np.amin(np.real(self.data1D))
            maxy = np.amax(np.real(self.data1D))
        elif self.plotType == 1:
            miny = np.amin(np.imag(self.data1D))
            maxy = np.amax(np.imag(self.data1D))
        elif self.plotType == 2:
            miny = np.amin((np.amin(np.real(self.data1D)), np.amin(np.imag(self.data1D))))
            maxy = np.amax((np.amax(np.real(self.data1D)), np.amax(np.imag(self.data1D))))
        elif self.plotType == 3:
            miny = np.amin(np.abs(self.data1D))
            maxy = np.amax(np.abs(self.data1D))
        else:
            miny = -1
            maxy = 1
        differ = 0.05*(maxy-miny) #amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny-differ
            self.ymaxlim = maxy+differ
        if xReset:
            xaxZlims = (self.xax > self.zminlim) & (self.xax < self.zmaxlim)
            if self.spec == 1:
                if self.ppm:
                    axMult = 1e6/self.ref
                else:
                    axMult = 1.0/(1000.0**self.axType)
            elif self.spec == 0:
                axMult = 1000.0**self.axType
            self.xminlim = min(self.xax[xaxZlims]*axMult)
            self.xmaxlim = (max(self.xax[xaxZlims])+(len(self.data1D)-1)*self.spacing)*axMult
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)

#########################################################################################################
#the class from which the contour data is displayed, the operations which only edit the content of this class are for previewing
class CurrentContour(Current1D):

    X_RESIZE = False
    Y_RESIZE = True
    
    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        self.data = data
        if hasattr(duplicateCurrent, 'axes2'):
            self.axes2 = duplicateCurrent.axes2
        else:
            self.axes2 = len(self.data.data.shape)-2
            if hasattr(duplicateCurrent, 'axes'):
                if self.axes2 == duplicateCurrent.axes:
                    self.axes2 = (self.axes2-1) % self.data.dim
        if hasattr(duplicateCurrent, 'axType2'):
            self.axType2 = duplicateCurrent.axType2
        else:
            self.axType2 = 1
        if hasattr(duplicateCurrent, 'ppm2'):
            self.ppm2 = duplicateCurrent.ppm2
        else:
            self.ppm2 = False
        if hasattr(duplicateCurrent, 'numLevels'):
            self.numLevels = duplicateCurrent.numLevels
        else:
            self.numLevels = 20
        if hasattr(duplicateCurrent, 'minLevels'):
            self.minLevels = duplicateCurrent.minLevels
        else:
            self.minLevels = 0.1
        if hasattr(duplicateCurrent, 'maxLevels'):
            self.maxLevels = duplicateCurrent.maxLevels
        else:
            self.maxLevels = 1.0
        if hasattr(duplicateCurrent, 'projType1'):
            self.projType1 = duplicateCurrent.projType1
        else:
            self.projType1 = 0
        if hasattr(duplicateCurrent, 'projType2'):
            self.projType2 = duplicateCurrent.projType2
        else:
            self.projType2 = 0
        Current1D.__init__(self, root, fig, canvas, data, duplicateCurrent)
        
    def copyCurrent(self, root, fig, canvas, data):
        return CurrentContour(root, fig, canvas, data, self)
        
    def upd(self): #get new data from the data instance
        if self.data.dim < 2:
            self.root.rescue()
            return False
        if (len(self.locList)+2) != self.data.dim:
            self.resetLocList()
        try:
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList)
        except:
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
        self.single = self.data1D.shape[-1] == 1
        return True

    def setBlock(self, axes, axes2, locList): #change the slice
        axesSame = True
        if (self.axes != axes) or (self.axes2 != axes2) :
            axesSame = False
        self.axes = axes
        self.axes2 = axes2
        self.locList = locList
        self.upd()
        if not axesSame:
            self.plotReset()
        self.showFid()

    def setLevels(self, numLevels, maxLevels, minLevels):
        self.numLevels = numLevels
        self.maxLevels = maxLevels
        self.minLevels = minLevels
        self.showFid()
        
    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-2)
        
    def setAxType2(self, val):
        if self.spec2 == 1:
            if self.ppm2:
                oldAxMult = 1e6/self.ref2
            else:
                oldAxMult = 1.0/(1000.0**self.axType2)
        elif self.spec2 == 0:
            oldAxMult = 1000.0**self.axType2
        if self.spec2 == 1:
            if val == 3:
                newAxMult = 1e6/self.ref2
            else:
                newAxMult = 1.0/(1000.0**val)
        elif self.spec2 == 0:
            newAxMult = 1000.0**val 
        if val == 3:
            self.ppm2 = True
        else:
            self.ppm2 = False
            self.axType2 = val
        self.yminlim = self.yminlim * newAxMult / oldAxMult
        self.ymaxlim = self.ymaxlim * newAxMult / oldAxMult
        self.showFid()

    def setProjType(self, val, direc):
        if direc == 1:
            self.projType1 = val
        if direc == 2:
            self.projType2 = val
        
    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=None): #display the 1D data including the apodization function
        t = np.arange(0, len(self.data1D[0]))/(self.sw)
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                self.dispMsg('shiftingAxes cannot be equal to axes')
            elif shiftingAxes == self.axes2:
                ar = np.arange(self.data.data.shape[self.axes2])
                x = np.ones((len(ar), len(self.data1D[0])))
                for i in range(len(ar)):
                    shift1 = shift + shifting*ar[i]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                    t2 = t - shift1
                    x2 = np.ones(len(self.data1D[0]))
                    if lor is not None:
                        x2 = x2*np.exp(-np.pi*lor*abs(t2))
                    if gauss is not None:
                        x2 = x2*np.exp(-((np.pi*gauss*t2)**2)/(4*np.log(2)))
                    if cos2 is not None:
                        x2 = x2*(np.cos(cos2*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
                    if hamming is not None:
                        alpha = 0.53836 # constant for hamming window
                        x2 = x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, np.pi, len(self.data1D[0])))))
                    if self.wholeEcho:
                        x2[-1:-(len(x2)/2+1):-1] = x2[:len(x2)/2]
                    x[i] = x2
            else:
                if (shiftingAxes < self.axes) and (shiftingAxes < self.axes2):
                    shift += shifting*self.locList[shiftingAxes]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                elif (shiftingAxes > self.axes) and (shiftingAxes > self.axes2):
                    shift += shifting*self.locList[shiftingAxes-2]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                else:
                    shift += shifting*self.locList[shiftingAxes-1]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                t2 = t - shift
                x = np.ones(len(self.data1D[0]))
                if lor is not None:
                    x = x*np.exp(-lor*abs(t2))
                if gauss is not None:
                    x = x*np.exp(-(gauss*t2)**2)
                if cos2 is not None:
                    x = x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
                if hamming is not None:
                    alpha = 0.53836 # constant for hamming window
                    x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, np.pi, len(self.data1D[0])))))
                if self.wholeEcho:
                    x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
                x = np.repeat([x], len(self.data1D), axis=0)
        else:
            t2 = t - shift
            x = np.ones(len(self.data1D[0]))
            if lor is not None:
                x = x*np.exp(-lor*abs(t2))
            if gauss is not None:
                x = x*np.exp(-(gauss*t2)**2)
            if cos2 is not None:
                x = x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
            if hamming is not None:
                alpha = 0.53836 # constant for hamming window
                x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, np.pi, len(self.data1D[0])))))
            if self.wholeEcho:
                x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
            x = np.repeat([x], len(self.data1D), axis=0)
        y = self.data1D
        self.ax.cla()
        if self.spec == 1:
            y = np.fft.ifftn(np.fft.ifftshift(y, axes=1), axes=[1])
            y = y*x
            y = np.fft.fftshift(np.fft.fftn(y, axes=[1]), axes=1)
        else:
            y = y*x
        self.showFid(y)
        
    def showFid(self, tmpdata=None): #display the 1D data
        self.peakPickReset()
        if tmpdata is None:
            tmpdata = self.data1D
        self.ax.cla()
        self.x_ax.cla()
        self.y_ax.cla()
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if self.spec2 == 1:
            if self.ppm2:
                axMult2 = 1e6/self.ref2
            else:
                axMult2 = 1.0/(1000.0**self.axType2)
        elif self.spec2 == 0:
            axMult2 = 1000.0**self.axType2
        x = self.xax*axMult
        self.line_xdata = x
        y = self.xax2*axMult2
        X, Y = np.meshgrid(x, y)
        if (self.plotType == 0):
            tmpdata = np.real(tmpdata)
        elif(self.plotType == 1):
            tmpdata = np.imag(tmpdata)
        elif(self.plotType == 2):
            tmpdata = np.real(tmpdata)
        elif(self.plotType == 3):
            tmpdata = np.abs(tmpdata)   
        differ = np.amax(np.abs(tmpdata))
        contourLevels = np.linspace(self.minLevels*differ, self.maxLevels*differ, self.numLevels)
        vmax = max(np.abs(self.minLevels*differ), np.abs(self.maxLevels*differ))
        vmin = -vmax
        if self.contourConst:
            self.ax.contour(X, Y, tmpdata, colors=self.contourColors[0], levels=contourLevels, vmax=vmax, vmin=vmin, linewidths=self.linewidth, label=self.data.name, linestyles='solid')
            self.ax.contour(X, Y, tmpdata, colors=self.contourColors[1], levels=-contourLevels[::-1], vmax=vmax, vmin=vmin, linewidths=self.linewidth, linestyles='solid')
        else:
            self.ax.contour(X, Y, tmpdata, cmap=get_cmap(self.colorMap), levels=contourLevels, vmax=vmax, vmin=vmin, linewidths=self.linewidth, label=self.data.name, linestyles='solid')
            self.ax.contour(X, Y, tmpdata, cmap=get_cmap(self.colorMap), levels=-contourLevels[::-1], vmax=vmax, vmin=vmin, linewidths=self.linewidth, linestyles='solid')
        self.line_ydata = tmpdata[0]
        if self.projType1 == 0:
            self.x_ax.plot(x, np.sum(tmpdata, axis=0), color=self.color, linewidth=self.linewidth)
        elif self.projType1 == 1:
            self.x_ax.plot(x, np.max(tmpdata, axis=0), color=self.color, linewidth=self.linewidth)
        elif self.projType1 == 2:
            self.x_ax.plot(x, np.min(tmpdata, axis=0), color=self.color, linewidth=self.linewidth)
        if self.projType2 == 0:
            self.y_ax.plot(np.sum(tmpdata, axis=1), y, color=self.color, linewidth=self.linewidth)
        elif self.projType2 == 1:
            self.y_ax.plot(np.max(tmpdata, axis=1), y, color=self.color, linewidth=self.linewidth)
        elif self.projType2 == 2:
            self.y_ax.plot(np.min(tmpdata, axis=1), y, color=self.color, linewidth=self.linewidth, )
        if self.spec == 0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(u'Time [\u03BCs]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        if self.spec2 == 0:
            if self.axType2 == 0:
                self.ax.set_ylabel('Time [s]')
            elif self.axType2 == 1:
                self.ax.set_ylabel('Time [ms]')
            elif self.axType2 == 2:
                self.ax.set_ylabel(u'Time [\u03BCs]')
            else:
                self.ax.set_ylabel('User defined')
        elif self.spec2 == 1:
            if self.ppm2:
                self.ax.set_ylabel('Frequency [ppm]')
            else:
                if self.axType2 == 0:
                    self.ax.set_ylabel('Frequency [Hz]')
                elif self.axType2 == 1:
                    self.ax.set_ylabel('Frequency [kHz]')
                elif self.axType2 == 2:
                    self.ax.set_ylabel('Frequency [MHz]')
                else:
                    self.ax.set_ylabel('User defined')
        else:
            self.ax.set_ylabel('')
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
        self.ax.xaxis.grid(self.grids[0])
        self.ax.yaxis.grid(self.grids[1])
        self.canvas.draw()

    def plotReset(self, xReset=True, yReset=True): #set the plot limits to min and max values
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if xReset:
            self.xminlim = min(self.xax*axMult)
            self.xmaxlim = max(self.xax*axMult)
        if self.spec2 == 1:
            if self.ppm2:
                axMult2 = 1e6/self.ref2
            else:
                axMult2 = 1.0/(1000.0**self.axType2)
        elif self.spec2 == 0:
            axMult2 = 1000.0**self.axType2
        if yReset:
            self.yminlim = min(self.xax2*axMult2)
            self.ymaxlim = max(self.xax2*axMult2)
        if self.spec:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        if self.spec2:
            self.ax.set_ylim(self.ymaxlim, self.yminlim)
        else:
            self.ax.set_ylim(self.yminlim, self.ymaxlim)

    #The peakpicking function needs to be changed for contour plots
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
                    if self.spec == 1:
                        if self.ppm:
                            axMult = 1e6/self.ref
                        else:
                            axMult = 1.0/(1000.0**self.axType)
                    elif self.spec == 0:
                        axMult = 1000.0**self.axType
                    xdata = self.xax*axMult
                    if self.spec2 == 1:
                        if self.ppm2:
                            axMult2 = 1e6/self.ref2
                        else:
                            axMult2 = 1.0/(1000.0**self.axType2)
                    elif self.spec2 == 0:
                        axMult2 = 1000.0**self.axType2
                    ydata = self.xax2*axMult2
                    idx = np.argmin(np.abs(xdata-event.xdata))
                    idy = np.argmin(np.abs(ydata-event.ydata))
                    if self.peakPickFunc is not None:
                        if (self.plotType == 0):
                            tmpdata = np.real(self.data1D[idy, idx])
                        elif(self.plotType == 1):
                            tmpdata = np.imag(self.data1D[idy, idx])
                        elif(self.plotType == 2):
                            tmpdata = np.real(self.data1D[idy, idx])
                        elif(self.plotType == 3):
                            tmpdata = np.abs(self.data1D[idy, idx])
                        self.peakPickFunc((idx, xdata[idx], tmpdata, ydata[idy]))
                    if not self.peakPick: #check if peakpicking is still required
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


#########################################################################################################
#The skewed plot class
class CurrentSkewed(Current1D):

    X_RESIZE = False
    Y_RESIZE = True
    
    def __init__(self, root, fig, canvas, data, duplicateCurrent= None):
        self.data = data
        if hasattr(duplicateCurrent, 'axes2'):
            self.axes2 = duplicateCurrent.axes2
        else:
            self.axes2 = len(self.data.data.shape)-2
            if hasattr(duplicateCurrent, 'axes'):
                if self.axes2 == duplicateCurrent.axes:
                    self.axes2 = (self.axes2-1) % self.data.dim
        if hasattr(duplicateCurrent, 'stackBegin'):
            self.stackBegin = duplicateCurrent.stackBegin
        else:
            self.stackBegin = None
        if hasattr(duplicateCurrent, 'stackEnd'):
            self.stackEnd = duplicateCurrent.stackEnd
        else:
            self.stackEnd = None
        if hasattr(duplicateCurrent, 'stackStep'):
            self.stackStep = duplicateCurrent.stackStep
        else:
            self.stackStep = None
            if self.data.data.shape[self.axes2] > 100:
                self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
        if hasattr(duplicateCurrent, 'axType2'):
            self.axType2 = duplicateCurrent.axType2
        else:
            self.axType2 = 1
        if hasattr(duplicateCurrent, 'ppm2'):
            self.ppm2 = duplicateCurrent.ppm2
        else:
            self.ppm2 = False
        Current1D.__init__(self, root, fig, canvas, data, duplicateCurrent)

    def startUp(self, xReset=True, yReset=True):
        self.altReset()
        self.plotReset(xReset, yReset)
        self.setSkewed(-0.2, 70)
        
    def copyCurrent(self, root, fig, canvas, data):
        return CurrentSkewed(root, fig, canvas, data, self)
    
    def upd(self): #get new data from the data instance
        if self.data.dim < 2:
            self.root.rescue()
            return False
        if (len(self.locList)+2) != self.data.dim:
            self.resetLocList()
        try:
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList, self.stackBegin, self.stackEnd, self.stackStep)
        except:
            self.resetLocList()
            updateVar = self.data.getBlock(self.axes, self.axes2, self.locList, self.stackBegin, self.stackEnd, self.stackStep)
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
        self.single = self.data1D.shape[-1] == 1
        return True
 
    def setBlock(self, axes, axes2, locList, stackBegin=None, stackEnd=None, stackStep=None): #change the slice
        axesSame = True
        if (self.axes != axes) or (self.axes2 != axes2) :
            axesSame = False
        self.axes = axes
        self.axes2 = axes2
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        if self.data.data.shape[self.axes2] > 100:
            self.stackStep = 1+int(self.data.data.shape[self.axes2])/100
        self.locList = locList
        self.upd()
        if not axesSame:
            self.plotReset()
        self.showFid()

    def setSkewed(self, skewed, elevation):
        self.skewed = skewed
        self.elevation = elevation
        proj3d.persp_transformation = lambda zfront, zback : np.array([[1, 0, 0, 0], 
                                                                       [skewed , 1.0, 0, 0], 
                                                                       [0, 0, zfront, 0], 
                                                                       [0, 0, -0.00001, zback]])
        self.ax.view_init(elev=self.elevation, azim=180*np.arctan(skewed/np.sin(np.pi*elevation/180))/np.pi-90)
        self.showFid()
        
    def resetLocList(self):
        self.locList = [0]*(len(self.data.data.shape)-2)
        
    def stackSelect(self, stackBegin, stackEnd, stackStep):
        self.stackBegin = stackBegin
        self.stackEnd = stackEnd
        self.stackStep = stackStep
        self.upd()
        self.plotReset(False, True)
        self.showFid()

    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxes=None): #display the 1D data including the apodization function
        t = np.arange(0, len(self.data1D[0]))/(self.sw)
        if shiftingAxes is not None:
            if shiftingAxes == self.axes:
                self.dispMsg('shiftingAxes cannot be equal to axes')
            elif shiftingAxes == self.axes2:
                ar = np.arange(self.data.data.shape[self.axes2])[slice(self.stackBegin, self.stackEnd, self.stackStep)]
                x = np.ones((len(ar), len(self.data1D[0])))
                for i in range(len(ar)):
                    shift1 = shift + shifting*ar[i]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                    t2 = t - shift1
                    x2 = np.ones(len(self.data1D[0]))
                    if lor is not None:
                        x2 = x2*np.exp(-np.pi*lor*abs(t2))
                    if gauss is not None:
                        x2 = x2*np.exp(-((np.pi*gauss*t2)**2)/(4*np.log(2)))
                    if cos2 is not None:
                        x2 = x2*(np.cos(cos2*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
                    if hamming is not None:
                        alpha = 0.53836 # constant for hamming window
                        x2 = x2*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift1*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, np.pi, len(self.data1D[0])))))
                    if self.wholeEcho:
                        x2[-1:-(len(x2)/2+1):-1] = x2[:len(x2)/2]
                    x[i] = x2
            else:
                if (shiftingAxes < self.axes) and (shiftingAxes < self.axes2):
                    shift += shifting*self.locList[shiftingAxes]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                elif (shiftingAxes > self.axes) and (shiftingAxes > self.axes2):
                    shift += shifting*self.locList[shiftingAxes-2]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                else:
                    shift += shifting*self.locList[shiftingAxes-1]*self.data.data.shape[shiftingAxes]/self.data.sw[shiftingAxes]
                t2 = t - shift
                x = np.ones(len(self.data1D[0]))
                if lor is not None:
                    x = x*np.exp(-lor*abs(t2))
                if gauss is not None:
                    x = x*np.exp(-(gauss*t2)**2)
                if cos2 is not None:
                    x = x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
                if hamming is not None:
                    alpha = 0.53836 # constant for hamming window
                    x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, np.pi, len(self.data1D[0])))))
                if self.wholeEcho:
                    x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
                x = np.repeat([x], len(self.data1D), axis=0)
        else:
            t2 = t - shift
            x = np.ones(len(self.data1D[0]))
            if lor is not None:
                x = x*np.exp(-lor*abs(t2))
            if gauss is not None:
                x = x*np.exp(-(gauss*t2)**2)
            if cos2 is not None:
                x = x*(np.cos(cos2*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, 0.5*np.pi, len(self.data1D[0]))))**2)
            if hamming is not None:
                alpha = 0.53836 # constant for hamming window
                x = x*(alpha+(1-alpha)*np.cos(hamming*(-0.5*shift*np.pi*self.sw/len(self.data1D[0])+np.linspace(0, np.pi, len(self.data1D[0])))))
            if self.wholeEcho:
                x[-1:-(len(x)/2+1):-1] = x[:len(x)/2]
            x = np.repeat([x], len(self.data1D), axis=0)
        y = self.data1D
        self.ax.cla()
        if self.spec == 1:
            y = np.fft.ifftn(np.fft.ifftshift(y, axes=1), axes=[1])
            y = y*x
            y = np.fft.fftshift(np.fft.fftn(y, axes=[1]), axes=1)
        else:
            y = y*x
        if self.spec == 0:
            if self.plotType == 0:
                self.showFid(y, [t], x*np.amax(np.real(self.data1D)), ['g'], old=True)
            elif self.plotType == 1:
                self.showFid(y, [t], x*np.amax(np.imag(self.data1D)), ['g'], old=True)
            elif self.plotType == 2:
                self.showFid(y, [t], x*np.amax(np.amax(np.real(self.data1D)), np.amax(np.imag(self.data1D))), ['g'], old=True)
            elif self.plotType == 3:
                self.showFid(y, [t], x*np.amax(np.abs(self.data1D)), ['g'], old=True)
        else:
            self.showFid(y)

    def showFid(self, tmpdata=None, extraX=None, extraY=None, extraColor=None, old=False): #display the 1D data
        self.peakPickReset()
        if tmpdata is None:
            tmpdata = self.data1D
        self.ax.cla()
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if self.spec2 == 1:
            if self.ppm2:
                axMult2 = 1e6/self.ref2
            else:
                axMult2 = 1.0/(1000.0**self.axType2)
        elif self.spec2 == 0:
            axMult2 = 1000.0**self.axType2
        x = self.xax*axMult
        self.line_xdata = x
        y = self.xax2*axMult2
        if old:
            if (self.plotType == 0):
                for num in range(len(self.data1D)):
                    self.ax.plot(x, y[num]*np.ones(len(x)), np.real(self.data1D[num]), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 1):
                for num in range(len(self.data1D)):
                    self.ax.plot(x, y[num]*np.ones(len(x)), np.imag(self.data1D[num]), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 2):
                for num in range(len(self.data1D)):
                    self.ax.plot(x, y[num]*np.ones(len(x)), np.real(self.data1D[num]), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
            elif(self.plotType == 3):
                for num in range(len(self.data1D)):
                    self.ax.plot(x, y[num]*np.ones(len(x)), np.abs(self.data1D[num]), c='k', alpha=0.2, linewidth=self.linewidth, label=self.data.name+'_old')
        if (extraX is not None):
            for num in range(len(extraY)):
                self.ax.plot(extraX[0]*axMult, y[num]*np.ones(len(extraX[0])), extraY[num], c= extraColor[0], linewidth=self.linewidth)
        if (self.plotType == 0):
            self.line_ydata = np.real(tmpdata[0])
            for num in range(len(tmpdata)):
                self.ax.plot(x, y[num]*np.ones(len(x)), np.real(tmpdata[num]), c=self.color, linewidth=self.linewidth, label=self.data.name)
        elif(self.plotType == 1):
            self.line_ydata = np.imag(tmpdata[0])
            for num in range(len(tmpdata)):
                self.ax.plot(x, y[num]*np.ones(len(x)), np.imag(tmpdata[num]), c=self.color, linewidth=self.linewidth, label=self.data.name)
        elif(self.plotType == 2):
            self.line_ydata = np.real(tmpdata[0])
            for num in range(len(tmpdata)):
                self.ax.plot(x, y[num]*np.ones(len(x)), np.imag(tmpdata[num]), c='r', linewidth=self.linewidth, label=self.data.name+'_imag')
                self.ax.plot(x, y[num]*np.ones(len(x)), np.real(tmpdata[num]), c=self.color, linewidth=self.linewidth, label=self.data.name)
        elif(self.plotType == 3):
            self.line_ydata = np.abs(tmpdata[0])
            for num in range(len(tmpdata)):
                self.ax.plot(x, y[num]*np.ones(len(x)), np.abs(tmpdata[num]), c=self.color, linewidth=self.linewidth, label=self.data.name)
        if self.spec == 0:
            if self.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.axType == 2:
                self.ax.set_xlabel(u'Time [\u03BCs]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        if self.spec2 == 0:
            if self.axType2 == 0:
                self.ax.set_ylabel('Time [s]')
            elif self.axType2 == 1:
                self.ax.set_ylabel('Time [ms]')
            elif self.axType2 == 2:
                self.ax.set_ylabel(u'Time [\u03BCs]')
            else:
                self.ax.set_ylabel('User defined')
        elif self.spec2 == 1:
            if self.ppm2:
                self.ax.set_ylabel('Frequency [ppm]')
            else:
                if self.axType2 == 0:
                    self.ax.set_ylabel('Frequency [Hz]')
                elif self.axType2 == 1:
                    self.ax.set_ylabel('Frequency [kHz]')
                elif self.axType2 == 2:
                    self.ax.set_ylabel('Frequency [MHz]')
                else:
                    self.ax.set_ylabel('User defined')
        else:
            self.ax.set_ylabel('')
        if self.spec:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        if self.spec2:
            self.ax.set_ylim(self.ymaxlim, self.yminlim)
        else:
            self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.set_zlim(self.zminlim, self.zmaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.w_zaxis.line.set_lw(0.)
        self.ax.set_zticks([])
        self.ax.grid(False)
        self.ax.xaxis.pane.set_edgecolor('white')
        self.ax.yaxis.pane.set_edgecolor('white')
        self.ax.zaxis.pane.set_edgecolor('white')
        self.ax.xaxis.pane.fill = False
        self.ax.yaxis.pane.fill = False
        self.ax.zaxis.pane.fill = False
        self.canvas.draw()

    def plotReset(self, xReset=True, yReset=True): #set the plot limits to min and max values
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if xReset:
            self.xminlim = min(self.xax*axMult)
            self.xmaxlim = max(self.xax*axMult)
        if self.spec2 == 1:
            if self.ppm2:
                axMult2 = 1e6/self.ref2
            else:
                axMult2 = 1.0/(1000.0**self.axType2)
        elif self.spec2 == 0:
            axMult2 = 1000.0**self.axType2
        if yReset:
            self.yminlim = min(self.xax2*axMult2)
            self.ymaxlim = max(self.xax2*axMult2)
        if self.spec:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        if self.spec2:
            self.ax.set_ylim(self.ymaxlim, self.yminlim)
        else:
            self.ax.set_ylim(self.yminlim, self.ymaxlim)

    def altScroll(self, event):
        middle = (self.zmaxlim+self.zminlim)/2.0
        width = self.zmaxlim-self.zminlim
        width = width*0.9**event.step
        self.zmaxlim = middle+width/2.0
        self.zminlim = middle-width/2.0
        self.ax.set_zlim(self.zminlim, self.zmaxlim)
        self.canvas.draw()
        
    def altReset(self):
        if self.plotType == 0:
            minz = np.amin(np.real(self.data1D))
            maxz = np.amax(np.real(self.data1D))
        elif self.plotType == 1:
            minz = np.amin(np.imag(self.data1D))
            maxz = np.amax(np.imag(self.data1D))
        elif self.plotType == 2:
            minz = np.amin(np.amin(np.real(self.data1D)), np.amin(np.imag(self.data1D)))
            maxz = np.amax(np.amax(np.real(self.data1D)), np.amax(np.imag(self.data1D)))
        elif self.plotType == 3:
            minz = np.amin(np.abs(self.data1D))
            maxz = np.amax(np.abs(self.data1D))
        else:
            minz = -1
            maxz = 1
        differ = 0.05*(maxz-minz) #amount to add to show all datapoints (10%)
        self.zminlim = minz-differ
        self.zmaxlim = maxz+differ
        self.ax.set_zlim(self.zminlim, self.zmaxlim)
        self.canvas.draw()
