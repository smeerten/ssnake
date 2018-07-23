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
import copy
from matplotlib.pyplot import get_cmap
import matplotlib
import matplotlib.ticker as ticker
import spectrum as sc
from spectrumFrame import PlotFrame
import reimplement as reim

COLORMAPLIST = ['seismic', 'BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'rainbow', 'jet']
COLORCYCLE = list(matplotlib.rcParams['axes.prop_cycle'])
COLORCONVERTER = matplotlib.colors.ColorConverter()

MINXNUMTICKS = 12
MINYNUMTICKS = 8


##################################################################################################
# the class from which the 1d data is displayed, the operations which only edit the content of this class are for previewing


class Current1D(PlotFrame):

    X_RESIZE = False
    Y_RESIZE = False
    NDIM_PLOT = 1  # Number of dimensions in the plot
    MARKER = ''
    LINESTYLE = '-'

    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        super(Current1D, self).__init__(root, fig, canvas)
        self.data = data  # the actual spectrum instance
        self.data1D = None  # the data1D
        if duplicateCurrent is None:
            self.axes = np.array([len(self.data.shape()) - 1], dtype=int)
            self.resetLocList()
            self.viewSettings = {"plotType": 0,
                                 "showTitle": self.root.father.defaultShowTitle,
                                 "axType": np.array([self.root.father.defaultUnits] * self.NDIM_PLOT, dtype=int),
                                 "ppm": np.array([self.root.father.defaultPPM] * self.NDIM_PLOT, dtype=bool),   # display frequency as ppm
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
            if self.viewSettings['showTitle']:
                self.fig.suptitle(self.data.name)
            else:
                self.fig.suptitle('')
            self.startUp()
        else:
            self.axes = self.fixAxes(copy.deepcopy(duplicateCurrent.axes))
            self.locList = copy.deepcopy(duplicateCurrent.locList)
            self.viewSettings = copy.deepcopy(duplicateCurrent.viewSettings)
            self.viewSettings.update({"extraData": [],
                                      "extraLoc": [],
                                      "extraColor": [],
                                      "extraName": [],
                                      "extraAxes": [],
                                      "extraScale": [],
                                      "extraOffset": [],
                                      "extraShift": []})
            if len(self.viewSettings['axType']) != self.NDIM_PLOT:
                diff = self.NDIM_PLOT - len(self.viewSettings['axType'])
                self.viewSettings['axType'] = np.append(np.array([self.root.father.defaultUnits] * diff, dtype=int), self.viewSettings['axType'][-self.NDIM_PLOT:])
                self.viewSettings['ppm'] = np.append(np.array([self.root.father.defaultPPM] * diff, dtype=bool), self.viewSettings['ppm'][-self.NDIM_PLOT:])
            if type(self) not in (CurrentStacked, CurrentArrayed) or type(self) not in (CurrentStacked, CurrentArrayed):
                self.viewSettings.update({"stackBegin": None, "stackEnd": None, "stackStep": None})
            self.xminlim = duplicateCurrent.xminlim
            self.xmaxlim = duplicateCurrent.xmaxlim
            self.yminlim = duplicateCurrent.yminlim
            self.ymaxlim = duplicateCurrent.ymaxlim
            xReset = self.X_RESIZE or duplicateCurrent.X_RESIZE
            yReset = self.Y_RESIZE or duplicateCurrent.Y_RESIZE
            self.upd()  # get the first slice of data
            if self.viewSettings['showTitle']:
                self.fig.suptitle(self.data.name)
            else:
                self.fig.suptitle('')
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
        return self.getAxMult(self.spec(axis), self.getAxType(), self.getppm(), self.freq(axis), self.ref(axis))

    def getAxType(self, num=-1):
        return self.viewSettings["axType"][num]

    def getppm(self, num=-1):
        return self.viewSettings["ppm"][num]

    def getDataType(self, data):
        typeList = [np.real, np.imag, np.array, np.abs]
        return typeList[self.viewSettings["plotType"]](data)
    
    def startUp(self, xReset=True, yReset=True):
        self.showFid()  # plot the data
        self.plotReset(xReset, yReset)  # reset the axes limits

    def fixAxes(self, axes):
        if len(axes) != self.NDIM_PLOT:
            fullAxes = np.arange(self.data.ndim())
            fullAxes = np.delete(fullAxes, axes)
            diff = self.NDIM_PLOT - len(axes)
            return np.append(fullAxes[-diff:], axes[-self.NDIM_PLOT:])
        else:
            return axes
        
    def rename(self, name):
        self.data.rename(name)
        if self.viewSettings['showTitle']:
            self.fig.suptitle(self.data.name)
        else:
            self.fig.suptitle('')
        self.canvas.draw()

    def copyCurrent(self, root, fig, canvas, data):
        return Current1D(root, fig, canvas, data, self)

    def upd(self):  # get new data from the data instance
        if self.data.ndim() <= self.axes[-1]:
            self.axes = np.array([len(self.data.shape()) - 1])
        if len(self.locList) != self.data.ndim():
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
        if not np.array_equal(self.axes, axes):
            axesSame = False
            self.axes = axes
        self.locList = locList
        self.upd()
        self.showFid()
        if not axesSame:
            self.plotReset()

    def resetLocList(self):
        self.locList = [0] * self.data.ndim()

    def getSelect(self):
        tmp = np.array(self.locList, dtype=object)
        tmp[self.axes] = slice(None)
        return tmp

    def setGrids(self, grids):
        self.viewSettings["grids"] = grids

    def setDiagonal(self, diagonalBool=None, diagonalMult=None):
        if diagonalBool is not None:
            self.viewSettings["diagonalBool"] = diagonalBool
        if diagonalMult is not None:
            self.viewSettings["diagonalMult"] = diagonalMult
        self.showFid()

    def reload(self, *args):
        self.data.reload()
        self.upd()
        self.showFid()
        self.plotReset()
        self.root.addMacro(['reload'])
        
    def isComplex(self, *args):
        return self.data.isComplex(*args)

    def setNoUndo(self,val):
        self.root.addMacro(['setNoUndo', (val,)])
        self.data.setNoUndo(val)

    def real(self, *args):
        self.root.addMacro(['real', (self.axes[-1] - self.data.ndim(), )])
        self.data.real(self.axes[-1])
        self.upd()
        self.showFid()
    
    def imag(self, *args):
        self.root.addMacro(['imag', (self.axes[-1] - self.data.ndim(), )])
        self.data.imag(self.axes[-1])
        self.upd()
        self.showFid()
    
    def abs(self, *args):
        self.root.addMacro(['abs', (self.axes[-1] - self.data.ndim(), )])
        self.data.abs(self.axes[-1])
        self.upd()
        self.showFid()

    def conj(self, *args):
        self.root.addMacro(['conj', (self.axes[-1] - self.data.ndim(), )])
        self.data.conj(self.axes[-1])
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
        self.root.addMacro(['phase', (phase0, phase1, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.phase(phase0, phase1, self.axes[-1], selectSlice)
        self.upd()
        self.showFid()

    def complexFourier(self):  # fourier the actual data and replot
        self.root.addMacro(['complexFourier', (self.axes[-1] - self.data.ndim(), )])
        self.data.complexFourier(self.axes[-1])
        self.upd()
        if isinstance(self, (CurrentStacked, CurrentArrayed)):
            self.resetSpacing()
        self.showFid()
        self.plotReset()

    def realFourier(self):  # fourier the real data and replot
        self.root.addMacro(['realFourier', (self.axes[-1] - self.data.ndim(), )])
        self.data.realFourier(self.axes[-1])
        self.upd()
        if isinstance(self, (CurrentStacked, CurrentArrayed)):
            self.resetSpacing()
        self.showFid()
        self.plotReset()

    def fftshift(self, inv=False):  # fftshift the actual data and replot
        self.root.addMacro(['fftshift', (self.axes[-1] - self.data.ndim(), inv)])
        self.data.fftshift(self.axes[-1], inv)
        self.upd()
        self.showFid()

    def apodPreview(self, lor=None, gauss=None, cos2= [None, None], hamming=None, shift=0.0, shifting=0.0, shiftingAxis=None):  # display the 1D data including the apodization function
        y = self.data1D.data.copy()
        preview = True
        if shiftingAxis is None:
            curve = self.data1D.apodize(lor, gauss, cos2, hamming, shift, shifting, None, -1, preview=preview)
        else:
            if (self.ndim() > 1) and (shiftingAxis == self.axes[-2]):
                curve = self.data1D.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxis, -1, preview=preview)
            else:
                if shiftingAxis == self.axes[-1]:
                    raise sc.SpectrumException('shiftingAxis cannot be equal to axis')
                shift += shifting * self.locList[shiftingAxis] / self.data.sw[shiftingAxis]
                curve = self.data1D.apodize(lor, gauss, cos2, hamming, shift, 0.0, None, -1, preview=preview)
        if self.spec() == 0 and not isinstance(self, CurrentContour):
            tmp = self.getDataType(y.getHyperData(0))
            scale = np.max([np.real(tmp), np.imag(tmp)])
            self.showFid(y, curve[0], scale*np.array(curve[1]), extraColor='g')
        else:
            self.showFid(y)
        self.upd()
    
    def applyApod(self, lor=None, gauss=None, cos2= [None,None], hamming=None, shift=0.0, shifting=0.0, shiftingAxis=0, select=False):  # apply the apodization to the actual data
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['apodize', (lor, gauss, cos2, hamming, shift, shifting, shiftingAxis, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxis, self.axes[-1], selectSlice)
        self.upd()
        self.showFid()

    def setFreq(self, freq, sw):  # set the frequency of the actual data
        self.root.addMacro(['setFreq', (freq, sw, self.axes[-1] - self.data.ndim())])
        self.data.setFreq(freq, sw, self.axes[-1])
        self.upd()
        self.showFid()

    def scaleSw(self,scale):
        self.root.addMacro(['scaleSw', (scale, self.axes[-1] - self.data.ndim())])
        self.data.scaleSw(scale, self.axes[-1])
        self.upd()
        self.showFid()

    def setRef(self, ref):  # set the frequency of the actual data
        oldref = self.ref()
        self.root.addMacro(['setRef', (ref, self.axes[-1] - self.data.ndim())])
        self.data.setRef(ref, self.axes[-1])
        if ref is None:
            ref = self.freq()
        val = self.getAxType()
        if self.spec() == 1:
            if self.getppm():
                self.xminlim = (self.xminlim * oldref * 10**-6 + oldref - ref) / (ref * 10**-6)
                self.xmaxlim = (self.xmaxlim * oldref * 10**-6 + oldref - ref) / (ref * 10**-6)
            else:
                self.xminlim = self.xminlim + (oldref - ref) / 10**(val * 3)  # set new limits, and scale for axis type
                self.xmaxlim = self.xmaxlim + (oldref - ref) / 10**(val * 3)
        self.upd()
        self.showFid()

    def regrid(self, limits, numPoints):
        self.root.addMacro(['regrid', (limits, numPoints, self.axes[-1] - self.data.ndim())])
        self.data.regrid(limits, numPoints, self.axes[-1])
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

    def fwhm(self, minPeak, maxPeak, level, unitType=None):
        from scipy.interpolate import UnivariateSpline
        if unitType is None:
            axType = self.getAxType()
            ppm = self.getppm()
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
        spline = UnivariateSpline(x, tmpData - level * tmpData[minP:maxP][maxPos], s=0)
        zeroPos = spline.roots()
        left = zeroPos[zeroPos > maxX]
        right = zeroPos[zeroPos < maxX]
        if right.size > 0 and left.size > 0:
            return abs(left[0] - right[-1])
        else:
            return 0.0

    def COM(self, minPeak, maxPeak, unitType=None):  # Centre of Mass
        dim = len(minPeak)
        ppm = []
        axType = []
        if unitType is None:
            for i in range(dim):
                axType.append(self.getAxType(-(i + 1)))
                ppm.append(self.getppm(-(i+1)))
        else:
            axType = unitType
            for i in range(dim):
                if unitType[i] == 3:
                    ppm.append(1)
                else:
                    ppm.append(0)
        for i in range(dim):
            if minPeak[i] > maxPeak[i]:
                minPeak[i], maxPeak[i] = maxPeak[i], minPeak[i]
        tmpData = self.data1D.getHyperData(0)
        tmpData = np.real(self.getDataType(tmpData))
        #Slice data
        slc = ()
        for i in range(dim):
            slc = slc + (slice(minPeak[-(i+1)],maxPeak[-(i+1)],None),)
        tmpData = tmpData[slc]

        COM = []
        axes = list(range(dim))
        for i in range(dim):
            sumAxes = list(axes)
            del sumAxes[-(i+1)]
            tmpAxis = self.xax(-(i+1))
            tmpAxis = tmpAxis[minPeak[i]:maxPeak[i]] * self.getAxMult(self.spec(-(i+1)), axType[i], ppm[i], self.freq(-(i+1)), self.ref(-(i+1)))
            tmpDataNew = np.sum(tmpData,tuple(sumAxes))
            # COM = 1/M *sum(m_i * r_i)
            COM.append(1.0 / np.sum(tmpDataNew) * np.sum(tmpDataNew * tmpAxis))
        return COM

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
        self.root.addMacro(['resize', (size, pos, self.axes[-1] - self.data.ndim())])
        self.data.resize(size, pos, self.axes[-1])
        self.upd()
        self.showFid()
        if not self.spec():
            self.plotReset(True, False)

    def lpsvd(self, nAnalyse, nFreq, nPredict, Direction):
        self.root.addMacro(['lpsvd', (nAnalyse, nFreq, nPredict, Direction, self.axes[-1] - self.data.ndim())])
        self.data.lpsvd(nAnalyse, nFreq, nPredict, Direction, self.axes[-1])
        self.upd()
        self.showFid()

    def setSpec(self, val):  # change from time to freq domain of the actual data
        self.root.addMacro(['setSpec', (val, self.axes[-1] - self.data.ndim())])
        self.data.setSpec(val, self.axes[-1])
        self.upd()
        if isinstance(self, CurrentArrayed):
            self.resetSpacing()
        self.showFid()
        self.plotReset()

    def swapEcho(self, idx):
        self.root.addMacro(['swapEcho', (idx, self.axes[-1] - self.data.ndim())])
        self.data.swapEcho(idx, self.axes[-1])
        self.upd()
        self.showFid()

    def swapEchoPreview(self, idx):
        self.data1D.swapEcho(idx, -1)
        self.showFid()
        self.upd()

    def setWholeEcho(self, value):
        valBool = value != 0
        self.root.addMacro(['setWholeEcho', (valBool, self.axes[-1] - self.data.ndim())])
        self.data.setWholeEcho(valBool, self.axes[-1])
        self.upd()

    def shift(self, shift, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['shift', (shift, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.shift(shift, self.axes[-1], selectSlice)
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
        y = np.apply_along_axis(lambda data: self.baselinePolyFit(self.xax(), data, bArray, degree), self.axes[-1], self.data.getHyperData(0))
        y = np.real(self.getDataType(y))
        self.root.addMacro(['subtract', (y)])
        self.data.subtract(y)

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
        self.root.addMacro(['baselineCorrection', (y, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.baselineCorrection(y, self.axes[-1], select=selectSlice)

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
        y = self.baselinePolyFit(self.xax(), tmpData, bArray, degree)
        y = np.real(self.getDataType(y))
        self.resetPreviewRemoveList()
        if self.NDIM_PLOT > 1:
            if isinstance(self, CurrentContour):
                self.showFid()
            else:
                self.showFid(extraX=[self.xax()], extraY=[y]*self.len(-2), extraColor='g')
        else:
            self.showFid(extraX=[self.xax()], extraY=[y], extraColor='g')
        self.previewRemoveList(removeList, invert)
        self.upd()

    def previewRemoveList(self, removeList, invert=False):
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
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
        self.root.addMacro(['states', (self.axes[-1] - self.data.ndim(), )])
        self.data.states(self.axes[-1])
        self.upd()
        self.showFid()

    def statesTPPI(self):
        self.root.addMacro(['statesTPPI', (self.axes[-1] - self.data.ndim(), )])
        self.data.statesTPPI(self.axes[-1])
        self.upd()
        self.showFid()

    def echoAntiEcho(self):
        self.root.addMacro(['echoAntiEcho', (self.axes[-1] - self.data.ndim(), )])
        self.data.echoAntiEcho(self.axes[-1])
        self.upd()
        self.showFid()

    def matrixFuncs(self, func, name, pos1, pos2, newSpec=False):
        if newSpec:
            tmpData = copy.deepcopy(self.data)
            func(tmpData, pos1, pos2, self.axes[-1])
            return tmpData
        else:
            self.root.addMacro([name, (pos1, pos2, self.axes[-1] - self.data.ndim(), )])
            func(self.data, pos1, pos2, self.axes[-1])
            if self.upd():
                self.showFid()
                self.plotReset()
        
    def integrate(self, pos1, pos2, newSpec=False):
        return self.matrixFuncs(lambda obj, *args: obj.integrate(*args), 'integrate', pos1, pos2, newSpec)

    def sum(self, pos1, pos2, newSpec=False):
        return self.matrixFuncs(lambda obj, *args: obj.sum(*args), 'sum', pos1, pos2, newSpec)
        
    def max(self, pos1, pos2, newSpec=False):
        return self.matrixFuncs(lambda obj, *args: obj.max(*args), 'max', pos1, pos2, newSpec)

    def min(self, pos1, pos2, newSpec=False):
        return self.matrixFuncs(lambda obj, *args: obj.min(*args), 'min', pos1, pos2, newSpec)

    def argmax(self, pos1, pos2, newSpec=False):
        return self.matrixFuncs(lambda obj, *args: obj.argmax(*args), 'argmax', pos1, pos2, newSpec)

    def argmin(self, pos1, pos2, newSpec=False):
        return self.matrixFuncs(lambda obj, *args: obj.argmin(*args), 'argmin', pos1, pos2, newSpec)

    def average(self, pos1, pos2, newSpec=False):
        return self.matrixFuncs(lambda obj, *args: obj.average(*args), 'average', pos1, pos2, newSpec)

    def flipLR(self):
        self.data.flipLR(self.axes[-1])
        self.upd()
        self.showFid()
        self.root.addMacro(['flipLR', (self.axes[-1] - self.data.ndim(), )])

    def concatenate(self, axes):
        self.data.concatenate(axes)
        self.upd()
        self.showFid()
        self.plotReset()
        self.root.addMacro(['concatenate', (axes - self.data.ndim() - 1, )])

    def split(self, sections):
        self.data.split(sections, self.axes[-1])
        self.upd()
        self.showFid()
        self.plotReset()
        self.root.addMacro(['split', (sections, self.axes[-1] - self.data.ndim() + 1)])

    def diff(self):
        self.data.diff(self.axes[-1])
        self.upd()
        self.showFid()
        self.root.addMacro(['diff', (self.axes[-1] - self.data.ndim(), )])

    def cumsum(self):
        self.data.cumsum(self.axes[-1])
        self.upd()
        self.showFid()
        self.root.addMacro(['cumsum', (self.axes[-1] - self.data.ndim(), )])

    def insert(self, data, pos):
        self.root.addMacro(['insert', (data, pos, self.axes[-1] - self.data.ndim())])
        self.data.insert(data, pos, self.axes[-1])
        self.upd()
        self.showFid()
        self.plotReset()

    def delete(self, pos):
        self.data.delete(pos, self.axes[-1])
        self.upd()
        self.showFid()
        self.root.addMacro(['delete', (pos, self.axes[-1] - self.data.ndim())])

    def deletePreview(self, pos):
        self.data1D.delete(pos, -1)
        self.showFid()
        self.upd()

    def add(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['add', (data, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.add(data, self.axes[-1], select=selectSlice)
        self.upd()
        self.showFid()

    def subtract(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['subtract', (data, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.subtract(data, self.axes[-1], select=selectSlice)
        self.upd()
        self.showFid()

    def multiply(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['multiply', (data, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.multiply(data, self.axes[-1], select=selectSlice)
        self.upd()
        self.showFid()

    def multiplyPreview(self, data):
        self.data1D.multiply(data, -1)
        self.showFid()
        self.upd()
        
    def divide(self, data, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['divide', (data, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.divide(data, self.axes[-1], select=selectSlice)
        self.upd()
        self.showFid()

    def normalize(self, value, scale, type, select=False):
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['normalize', (value, scale, type, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.normalize(value, scale, type, self.axes[-1], select=selectSlice)
        self.upd()
        self.showFid()
    
    def subtractAvg(self, pos1, pos2):
        self.root.addMacro(['subtractAvg', (pos1, pos2, self.axes[-1] - self.data.ndim())])
        self.data.subtractAvg(pos1, pos2, self.axes[-1])
        self.upd()
        self.showFid()

    def subtractAvgPreview(self, pos1, pos2):
        self.data1D.subtractAvg(pos1, pos2, -1)
        self.showFid()
        self.upd()

    def extract(self, pos1, pos2, newSpec=False):
        if newSpec:
            tmpData = copy.deepcopy(self.data)
            tmpData.extract(pos1, pos2, self.axes[-1])
            return tmpData
        else:
            self.root.addMacro(['extract', (pos1, pos2, self.axes[-1] - self.data.ndim())])
            self.data.extract(pos1, pos2, self.axes[-1])
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
        self.root.addMacro(['fiddle', (refSpec.tolist(), lb, self.axes[-1] - self.data.ndim())])
        self.data.fiddle(refSpec, lb, self.axes[-1])
        self.upd()
        self.showFid()

    def shearing(self, shear, axes, axes2):
        self.root.addMacro(['shear', (shear, axes - self.data.ndim(), axes2 - self.data.ndim())])
        self.data.shear(shear, axes, axes2)
        self.upd()
        self.showFid()

    def reorder(self, pos, newLength):
        self.root.addMacro(['reorder', (pos, newLength, self.axes[-1] - self.data.ndim())])
        self.data.reorder(pos, newLength, self.axes[-1])
        self.upd()
        self.showFid()

    def ffm(self, posList, typeVal):
        self.root.addMacro(['ffm', (posList, typeVal, self.axes[-1] - self.data.ndim())])
        self.data.ffm(posList, typeVal, self.axes[-1])
        self.upd()
        self.showFid()

    def clean(self, posList, typeVal, gamma, threshold, maxIter):
        self.root.addMacro(['clean', (posList, typeVal, self.axes[-1] - self.data.ndim(), gamma, threshold, maxIter)])
        self.data.clean(posList, typeVal, self.axes[-1], gamma, threshold, maxIter)
        self.upd()
        self.showFid()

    def ist(self, posList, typeVal, threshold, maxIter, tracelimit):
        self.root.addMacro(['ist', (posList, typeVal, self.axes[-1] - self.data.ndim(), threshold, maxIter, tracelimit)])
        self.data.ist(posList, typeVal, self.axes[-1], threshold, maxIter, tracelimit)
        self.upd()
        self.showFid()

    def autoPhase(self, phaseNum):
        phases = self.data1D.autoPhase(phaseNum, -1, [0]*self.ndim(), returnPhases=True)
        self.upd()
        return phases

    def directAutoPhase(self, phaseNum):
        tmpLocList = copy.copy(self.locList)
        if self.ndim() > 1:
            tmpLocList[self.axes[:-1]] = self.viewSettings["stackBegin"]
        self.root.addMacro(['autoPhase', (phaseNum, self.axes[-1] - self.data.ndim(), tmpLocList)])
        self.data.autoPhase(phaseNum, self.axes[-1], tmpLocList)
        self.upd()
        self.showFid()

    def setXaxPreview(self, xax):
        self.data1D.setXax(xax, -1)
        self.showFid()
        self.plotReset()
        self.upd()

    def setXax(self, xax):
        self.root.addMacro(['setXax', (xax, self.axes[-1] - self.data.ndim())])
        self.data.setXax(xax, self.axes[-1])
        self.upd()
        self.showFid()
        self.plotReset()

    def setAxType(self, val, update=True, num=-1):
        oldAxMult = self.getAxMult(self.spec(num), self.getAxType(num), self.getppm(num), self.freq(num), self.ref(num))
        if val == 3:
            self.viewSettings["ppm"][num] = True
        else:
            self.viewSettings["ppm"][num] = False
            self.viewSettings["axType"][num] = val
        newAxMult = self.getAxMult(self.spec(num), self.getAxType(num), self.getppm(num), self.freq(num), self.ref(num))
        if num == -1 or num == (self.ndim()-1):
            self.xminlim = self.xminlim * newAxMult / oldAxMult
            self.xmaxlim = self.xmaxlim * newAxMult / oldAxMult
        elif num == -2 or num == (self.ndim()-2):
            self.yminlim = self.yminlim * newAxMult / oldAxMult
            self.ymaxlim = self.ymaxlim * newAxMult / oldAxMult
        if update:
            self.showFid()

    def hilbert(self):
        self.root.addMacro(['hilbert', (self.axes[-1] - self.data.ndim(), )])
        self.data.hilbert(self.axes[-1])
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
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
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
        if extraX is not None:
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
        self.ax.set_xlabel(self.getLabel(self.spec(), self.axes[-1], self.getAxType(), self.getppm()))
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
        self.setTicks()
        self.canvas.draw()

    def setTicks(self,Xset = True,Yset = True):
        if  matplotlib.__version__[0] == '2':
            if Xset:
                self.ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins='auto', steps=[1,2,2.5,5,10], min_n_ticks=MINXNUMTICKS))
            if Yset:
                self.ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins='auto', steps=[1,2,2.5,5,10], min_n_ticks=MINYNUMTICKS))

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
        if isinstance(self, CurrentContour):
            differ = 0
        else:
            differ = 0.05 * (maxy - miny)  # amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny - differ
            self.ymaxlim = maxy + differ
        if self.ref() == 0.0:
            self.viewSettings["ppm"][-1] = False
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
        if isinstance(self, CurrentContour): #If contour: reverse y-axis
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
        self.viewSettings["extraLoc"].append([0] * (len(self.viewSettings["extraData"][-1].shape()) ))
        self.viewSettings["extraColor"].append(COLORCONVERTER.to_rgb(COLORCYCLE[np.mod(len(self.viewSettings["extraData"]), len(COLORCYCLE))]['color']))  # find a good color system
        self.viewSettings["extraAxes"].append([len(data.shape()) - 1])
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
                self.viewSettings["extraLoc"][i] = [0] * (len(self.viewSettings["extraData"][i].shape()))
        else:
            self.viewSettings["extraLoc"][num] = [0] * (len(self.viewSettings["extraData"][num].shape()))

    def showFid(self, oldData=None, extraX=None, extraY=None, extraColor=None):  # display the 1D data
        self.peakPickReset()
        tmpdata = self.data1D.getHyperData(0)
        self.ax.cla()
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        self.line_xdata = [self.xax() * axMult]
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        for i in range(len(self.viewSettings["extraData"])):
            data = self.viewSettings["extraData"][i]
            try:
                if self.viewSettings["extraData"][i].ndim() <= self.viewSettings["extraAxes"][i][-1]:
                    self.viewSettings["extraAxes"][i][-1] = len(self.viewSettings["extraData"][i].shape()) - 1
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
            self.line_xdata_extra.append(self.viewSettings["extraShift"][i] + xax * self.getAxMult(spec, self.getAxType(), self.getppm(), freq, ref))
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
        self.ax.set_xlabel(self.getLabel(self.spec(), self.axes[-1], self.getAxType(), self.getppm()))
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.setTicks()
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
    ZERO_SCROLL_ALLOWED = False
    NDIM_PLOT = 2

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
        if len(self.locList) != self.data.ndim():
            self.resetLocList()
        stack = [reim.floatSlice(self.viewSettings["stackBegin"], self.viewSettings["stackEnd"], self.viewSettings["stackStep"]), slice(None)]
        try:
            self.data1D = self.data.getSlice(self.axes, self.locList, stack)
            if self.data1D is None:
                self.root.rescue()
        except Exception:
            self.resetLocList()
            self.data1D = self.data.getSlice(self.axes, self.locList, stack)
        return True

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
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
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
        self.ax.set_xlabel(self.getLabel(self.spec(), self.axes[-1], self.getAxType(), self.getppm()))
        if self.spec() > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.setTicks()
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.canvas.draw()

#########################################################################################################
# the class from which the arrayed data is displayed, the operations which only edit the content of this class are for previewing

class CurrentArrayed(CurrentStacked):

    X_RESIZE = True
    Y_RESIZE = False
    INVERT_X = False
    ZERO_SCROLL_ALLOWED = True

    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        if duplicateCurrent is not None:
            if isinstance(duplicateCurrent, CurrentArrayed):
                self.zminlim = duplicateCurrent.zminlim
                self.zmaxlim = duplicateCurrent.zmaxlim
            else:
                # The z-axes limits are in xax units unlike the x-axes and y-axes limits
                axMult = self.getAxMult(duplicateCurrent.spec(),
                                        duplicateCurrent.getAxType(),
                                        duplicateCurrent.getppm(),
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
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        axMult2 = self.getAxMult(self.spec(-2), self.getAxType(-2), self.getppm(-2), self.freq(-2), self.ref(-2))
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
        self.setTicks(Xset = False)
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
    GRID_PLOT = True
    INVERT_Y = True
    ZERO_SCROLL_ALLOWED = False

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
        self.showFid()
        #self.plotContour(updateOnly=True)

    def copyCurrent(self, root, fig, canvas, data):
        return CurrentContour(root, fig, canvas, data, self)

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
        self.ax.set_xlabel(self.getLabel(self.spec(), self.axes[-1], self.getAxType(), self.getppm()))
        self.ax.set_ylabel(self.getLabel(self.spec(-2), self.axes[-2], self.getAxType(-2), self.getppm(-2)))
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

    def showFid(self, oldData=None, extraX=None, extraY=None, extraZ=None, extraColor=None):
        # The oldData and extra plots are not displayed in the contourplot for now
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        self.line_zdata_extra = []
        self.differ = None
        self.peakPickReset()
        tmpdata = np.real(self.getDataType(self.data1D.getHyperData(0)))
        if self.viewSettings["limitType"] == 0:
            self.differ = np.max(np.abs(tmpdata))
        else:
            self.differ = np.max(np.abs(np.ravel(self.data.getHyperData)))
        self.ax.cla()
        self.clearProj()
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        axMult2 = self.getAxMult(self.spec(-2), self.getAxType(-2), self.getppm(-2), self.freq(-2), self.ref(-2))
        if self.viewSettings["diagonalBool"]:
            add_diagonal(self.ax, self.viewSettings["diagonalMult"], c='k', ls='--')
        if oldData is not None:
            tmp = np.real(self.getDataType(oldData.getHyperData(0)))
            self.line_xdata_extra.append(self.xax() * axMult)
            self.line_ydata_extra.append(self.xax(-2) * axMult2)
            self.line_zdata_extra.append(tmp)
            self.plotContour(self.line_xdata_extra[-1], self.line_ydata_extra[-1], self.line_zdata_extra[-1], color=['k','k'], alpha=0.2)
            self.showProj(self.line_xdata_extra[-1], self.line_ydata_extra[-1], self.line_zdata_extra[-1], 'k')
        if extraX is not None:
            for num in range(len(extraX)):
                self.line_xdata_extra.append(extraX[num] * axMult)
                self.line_ydata_extra.append(extraY[num] * axMult2)
                self.line_zdata_extra.append(extraZ[num])
                self.plotContour(self.line_xdata_extra[-1], self.line_ydata_extra[-1], self.line_zdata_extra[-1], color=[extraColor[num], extraColor[num]])
                self.showProj(self.line_xdata_extra[-1], self.line_ydata_extra[-1], self.line_zdata_extra[-1], extraColor[num])
        self.line_xdata = [self.xax() * axMult]
        self.line_ydata = [self.xax(-2) * axMult2]
        self.line_zdata = [tmpdata]
        if self.viewSettings["limitType"] == 0:
            self.differ = np.max(np.abs(tmpdata))
        else:
            self.differ = np.max(np.abs(np.ravel(self.data.getHyperData)))
        self.plotContour(self.line_xdata[-1], self.line_ydata[-1], self.line_zdata[-1])
        self.showProj(self.line_xdata[-1], self.line_ydata[-1], self.line_zdata[-1])
        self.ax.set_xlabel(self.getLabel(self.spec(), self.axes[-1], self.getAxType(), self.getppm()))
        self.ax.set_ylabel(self.getLabel(self.spec(-2), self.axes[-2], self.getAxType(-2), self.getppm(-2)))
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
        self.setTicks()
        self.canvas.draw()

    def plotContour(self, line_xdata, line_ydata, line_zdata, color=None, alpha=1, updateOnly=False):  # Plots the contour plot
        if color is None and self.viewSettings["contourConst"]:
            color = self.viewSettings["contourColors"]
        X, Y = np.meshgrid(line_xdata, line_ydata)
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
            if line_zdata.shape[0] > 2:  # if size 2 or lower, convolve fails, just take whole data then
                YposMax = np.where(np.convolve(np.max(line_zdata, 1) > contourLevels[0], [True, True, True], 'same'))[0]
            else:
                YposMax = np.arange(line_zdata.shape[0])
            if YposMax.size > 0:  # if number of positive contours is non-zero
                if line_zdata.shape[1] > 2:
                    XposMax = np.where(np.convolve(np.max(line_zdata, 0) > contourLevels[0], [True, True, True], 'same'))[0]
                else:
                    XposMax = np.arange(line_zdata.shape[1])
                PlotPositive = True
        PlotNegative = False
        if self.viewSettings["contourSign"] == 0 or self.viewSettings["contourSign"] == 2:
            if not self.viewSettings["plotType"] == 3:  # for Absolute plot no negative
                if line_zdata.shape[0] > 2:
                    YposMin = np.where(np.convolve(np.min(line_zdata, 1) < -contourLevels[0], [True, True, True], 'same'))[0]
                else:
                    YposMin = np.arange(line_zdata.shape[0])
                if YposMin.size > 0:  # if number of negative contours is non-zero
                    if line_zdata.shape[1] > 2:
                        XposMin = np.where(np.convolve(np.min(line_zdata, 0) < -contourLevels[0], [True, True, True], 'same'))[0]
                    else:
                        XposMin = np.arange(line_zdata.shape[1])
                    PlotNegative = True
        vmax = max(np.abs(self.viewSettings["minLevels"] * self.differ), np.abs(self.viewSettings["maxLevels"] * self.differ))
        vmin = -vmax
        if color is not None:
            if PlotPositive:
                self.ax.contour(X[YposMax[:,None],XposMax],Y[YposMax[:,None],XposMax],line_zdata[YposMax[:,None],XposMax], colors=color[0], alpha=alpha, levels=contourLevels, vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], linestyles='solid')
            if PlotNegative:
                self.ax.contour(X[YposMin[:,None],XposMin],Y[YposMin[:,None],XposMin],line_zdata[YposMin[:,None],XposMin], colors=color[1], alpha=alpha, levels=-contourLevels[::-1], vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], linestyles='solid')
        else:
            if PlotPositive:
                self.ax.contour(X[YposMax[:,None],XposMax],Y[YposMax[:,None],XposMax],line_zdata[YposMax[:,None],XposMax], cmap=get_cmap(self.viewSettings["colorMap"]), levels=contourLevels, vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], linestyles='solid')
            if PlotNegative:    
                self.ax.contour(X[YposMin[:,None],XposMin],Y[YposMin[:,None],XposMin],line_zdata[YposMin[:,None],XposMin], cmap=get_cmap(self.viewSettings["colorMap"]), levels=-contourLevels[::-1], vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], linestyles='solid')
        self.setTicks()
        if updateOnly:
            self.canvas.draw()

    def clearProj(self):
        self.x_ax.cla()
        self.y_ax.cla()

    def showProj(self, line_xdata=None, line_ydata=None, line_zdata=None, color=None):
        xLimOld = self.x_ax.get_xlim()
        if line_xdata is None:
            x = self.line_xdata[-1]
            y = self.line_ydata[-1]
            tmpdata = np.real(self.getDataType(self.data1D.getHyperData(0)))
        else:
            x = line_xdata
            y = line_ydata
            tmpdata = line_zdata
        yLimOld = self.y_ax.get_ylim()
        if color is None:
            color = self.viewSettings["color"]
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
            if self.viewSettings["projPos"][0] >= self.data.shape()[self.axes[-2]]:
                self.viewSettings["projPos"][0] = self.data.shape()[self.axes[-2]] - 1
            elif self.viewSettings["projPos"][0] < 0:
                self.viewSettings["projPos"][0] = 0
            xprojdata = tmpdata[self.viewSettings["projPos"][0]]
        if self.viewSettings["projTop"] != 3:
            self.x_ax.plot(x, xprojdata, color=color, linewidth=self.viewSettings["linewidth"], picker=True)            
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
            if self.viewSettings["projPos"][1] >= self.data.shape()[self.axes[-1]]:
                self.viewSettings["projPos"][1] = self.data.shape()[self.axes[-1]] - 1
            elif self.viewSettings["projPos"][1] < 0:
                self.viewSettings["projPos"][1] = 0
            yprojdata = tmpdata[:, self.viewSettings["projPos"][1]]
        if self.viewSettings["projRight"] != 3:
            self.y_ax.plot(yprojdata, y, color=color, linewidth=self.viewSettings["linewidth"], picker=True)
            ymin, ymax = np.min(yprojdata), np.max(yprojdata)
            self.y_ax.set_xlim([ymin - 0.15 * (ymax - ymin), ymax + 0.05 * (ymax - ymin)])  # Set projection limits, and force 15% whitespace below plot
            self.y_ax.set_ylim(yLimOld)
        self.setTicks()
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
                    xdata = self.xax() * self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
                    ydata = self.xax(-2) * self.getAxMult(self.spec(-2), self.getAxType(-2), self.getppm(-2), self.freq(-2), self.ref(-2))
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

