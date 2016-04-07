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
from PyQt4 import QtGui, QtCore
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import scipy.optimize
import os.path
import copy
from safeEval import *
from spectrumFrame import Plot1DFrame
from zcw import *
from widgetClasses import *

pi = np.pi

##############################################################################
class FittingWindow(QtGui.QWidget):
    #Inherited by the fitting windows
    def __init__(self, mainProgram, oldMainWindow):
        QtGui.QWidget.__init__(self, mainProgram)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        grid = QtGui.QGridLayout(self)
        grid.addWidget(self.canvas, 0, 0)
        self.setup()
        grid.addWidget(self.paramframe, 1, 0)
        grid.setColumnStretch(0, 1)
        grid.setRowStretch(0, 1)
        self.grid = grid
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)

    def setup(self):
        pass
        
    def createNewData(self, data, axes, store=False):
        masterData = self.get_masterData()
        if store:
            self.mainProgram.dataFromFit(data, [masterData.freq[axes], masterData.freq[axes]], [masterData.sw[axes], masterData.sw[axes]] , [False, masterData.spec[axes]], [False, masterData.wholeEcho[axes]], [None, masterData.ref[axes]], [np.arange(len(data)), masterData.xaxArray[axes]], axes)
        else:
            self.mainProgram.dataFromFit(data, copy.deepcopy(masterData.freq), copy.deepcopy(masterData.sw) , copy.deepcopy(masterData.spec), copy.deepcopy(masterData.wholeEcho), copy.deepcopy(masterData.ref), copy.deepcopy(masterData.xaxArray), axes)
        
    def rename(self, name):
        self.fig.suptitle(name)
        self.canvas.draw()
        self.oldMainWindow.rename(name)
        
    def buttonPress(self, event):
        self.current.buttonPress(event)

    def buttonRelease(self, event):
        self.current.buttonRelease(event)

    def pan(self, event):
        self.current.pan(event)

    def scroll(self, event):
        self.current.scroll(event)
        
    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()
        
    def kill(self):
        for i in reversed(range(self.grid.count())): 
            self.grid.itemAt(i).widget().deleteLater()
        self.grid.deleteLater()
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        del self.paramframe
        del self.fig
        del self.canvas
        self.deleteLater()
        
    def cancel(self):
        for i in reversed(range(self.grid.count())): 
            self.grid.itemAt(i).widget().deleteLater()
        self.grid.deleteLater()
        del self.current
        del self.paramframe
        del self.canvas
        del self.fig
        self.mainProgram.closeFitWindow(self.oldMainWindow)
        self.deleteLater()

##############################################################################
class IntegralsWindow(FittingWindow): 
    def __init__(self, mainProgram, oldMainWindow):
        FittingWindow.__init__(self, mainProgram, oldMainWindow)

    def setup(self):
        self.current = IntegralsFrame(self, self.fig, self.canvas, self.oldMainWindow.current)
        self.paramframe = IntegralsParamFrame(self.current, self)
        
#################################################################################   
class IntegralsFrame(Plot1DFrame): 
    def __init__(self, rootwindow, fig, canvas, current):
        Plot1DFrame.__init__(self, rootwindow, fig, canvas)
        self.data1D = current.getDisplayedData()
        self.current = current
        self.spec = self.current.spec
        self.xax = self.current.xax
        self.plotType = 0
        self.rootwindow = rootwindow
        self.pickNum = 0
        self.pickWidth = False
        self.plotReset()
        self.showPlot()
        
    def plotReset(self): 
        a = self.fig.gca()
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
        differ = 0.05*(maxy-miny) 
        self.yminlim = miny-differ
        self.ymaxlim = maxy+differ
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xminlim = min(self.xax*axMult)
        self.xmaxlim = max(self.xax*axMult)
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim, self.xminlim)
        else:
            a.set_xlim(self.xminlim, self.xmaxlim)
        a.set_ylim(self.yminlim, self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        a = self.fig.gca()
        a.cla()
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.line_xdata = self.xax*axMult
        self.line_ydata = self.data1D
        a.plot(self.xax*axMult, self.data1D)
        if tmpAx is not None:
            a.plot(tmpAx*axMult, tmpdata)
        for i in range(len(tmpAx2)):
            a.plot(tmpAx2[i]*axMult, tmpdata2[i])
        if self.spec == 0:
            if self.current.axType == 0:
                a.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                a.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                a.set_xlabel(r'Time [$\mu$s]')
            else:
                a.set_xlabel('User defined')
        elif self.spec == 1:
            if self.current.ppm:
                a.set_xlabel('Frequency [ppm]')
            else:
                if self.current.axType == 0:
                    a.set_xlabel('Frequency [Hz]')
                elif self.current.axType == 1:
                    a.set_xlabel('Frequency [kHz]')
                elif self.current.axType == 2:
                    a.set_xlabel('Frequency [MHz]')
                else:
                    a.set_xlabel('User defined')
        else:
            a.set_xlabel('')
        a.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        a.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim, self.xminlim)
        else:
            a.set_xlim(self.xminlim, self.xmaxlim)
        a.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()

    def togglePick(self, var):
        self.peakPickReset()
        if var == 1:
            self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos) 
            self.peakPick = True
        else:
            self.peakPickFunc = None
            self.peakPick = False
        
    def pickDeconv(self, pos):
        self.rootwindow.paramframe.addValue(pos[1])
        self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos) 
        self.peakPick = True 

#################################################################################
class IntegralsParamFrame(QtGui.QWidget): 
    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        self.frame1 = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1, 0, 0)
        grid.addLayout(self.frame2, 0, 1)
        grid.addLayout(self.frame3, 0, 2)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton, 0, 0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton, 1, 0)
        resetButton = QtGui.QPushButton("Reset")
        resetButton.clicked.connect(self.reset)
        self.frame1.addWidget(resetButton, 0, 1)
        self.pickTick = QtGui.QCheckBox("Pick")
        self.pickTick.stateChanged.connect(self.togglePick)
        self.frame1.addWidget(self.pickTick, 1, 1)
        self.absIntTick = QtGui.QCheckBox("Relative integrals")
        self.absIntTick.setChecked(True)
        self.absIntTick.stateChanged.connect(self.fit)
        self.frame1.addWidget(self.absIntTick, 2, 0, 1, 2)
        self.frame1.setColumnStretch(10, 1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.frame3.addWidget(QLabel("Min bounds:"), 1, 0)
        self.frame3.addWidget(QLabel("Max bounds:"), 1, 1)
        self.frame3.addWidget(QLabel("Integral:"), 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        if self.parent.spec == 1:
            if self.parent.current.ppm:
                axMult = 1e6/self.parent.current.ref
            else:
                axMult = 1.0/(1000.0**self.parent.current.axType)
        elif self.parent.spec == 0:
            axMult = 1000.0**self.parent.current.axType
        self.xax = self.parent.current.xax * axMult
        self.integralIter = 0
        self.refVal = None
        self.minValues = np.array([min(self.xax)]) #dummy variables
        self.maxValues = np.array([max(self.xax)]) #dummy variables
        self.intValues = np.array([1.0]) #dummy variables
        self.minEntries = []
        self.maxEntries = []
        self.intEntries = []
        self.deleteButtons = []
        self.integralIter = 0
        self.entryCount = 1
        self.first = True
        if self.parent.current.plotType == 0:
            self.maxy = max(np.real(self.parent.current.data1D))
            self.diffy = self.maxy - min(np.real(self.parent.current.data1D))
        elif self.parent.current.plotType == 1:
            self.maxy = max(np.imag(self.parent.current.data1D))
            self.diffy = self.maxy - min(np.imag(self.parent.current.data1D))
        elif self.parent.current.plotType == 2:
            self.maxy = max(np.real(self.parent.current.data1D))
            self.diffy = self.maxy - min(np.real(self.parent.current.data1D))
        elif self.parent.current.plotType == 3:
            self.maxy = max(np.abs(self.parent.current.data1D))
            self.diffy = self.maxy - min(np.abs(self.parent.current.data1D))
        self.minEntries.append(QtGui.QLineEdit())
        self.minEntries[0].setAlignment(QtCore.Qt.AlignHCenter)
        self.minEntries[0].editingFinished.connect(lambda self=self: self.setVal(self.minEntries[0], True))
        self.frame3.addWidget(self.minEntries[0], 2, 0)
        self.maxEntries.append(QtGui.QLineEdit())
        self.maxEntries[0].setAlignment(QtCore.Qt.AlignHCenter)
        self.maxEntries[0].editingFinished.connect(lambda self=self: self.setVal(self.maxEntries[0], False))
        self.frame3.addWidget(self.maxEntries[0], 2, 1)
        self.intEntries.append(QtGui.QLineEdit())
        self.intEntries[0].setAlignment(QtCore.Qt.AlignHCenter)
        self.intEntries[0].editingFinished.connect(lambda self=self: self.setRef(self.intEntries[0]))
        self.frame3.addWidget(self.intEntries[0], 2, 2)
        self.deleteButtons.append(QtGui.QPushButton("X"))
        self.deleteButtons[0].clicked.connect(lambda extra, self=self: self.deleteEntry(self.deleteButtons[0]))
        self.frame3.addWidget(self.deleteButtons[0], 2, 3)
        self.minEntries.append(QtGui.QLineEdit())
        self.reset()
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        
    def reset(self):
        for i in range(self.integralIter+1):
            self.deleteEntry(num=0, reset=True)
        self.refVal = None
        self.pickTick.setChecked(True)
        self.togglePick()
        self.parent.showPlot()

    def addValue(self, value):
        if self.first:
            self.minValues[self.integralIter] = value
            self.minEntries[self.integralIter].setText("%#.3g" % value)
            self.first = False
        else:
            tmp = self.minValues[self.integralIter]
            self.minValues[self.integralIter] = min(value, tmp)
            self.maxValues[self.integralIter] = max(value, tmp)
            self.minValues = np.append(self.minValues, min(self.xax))
            self.maxValues = np.append(self.maxValues, max(self.xax))
            self.intValues = np.append(self.intValues, 1)
            self.minEntries[self.integralIter].setText("%#.3g" % min(value, tmp))
            self.maxEntries[self.integralIter].setText("%#.3g" % max(value, tmp))
            self.integralIter += 1
            self.minEntries.append(QtGui.QLineEdit())
            self.minEntries[self.integralIter].setAlignment(QtCore.Qt.AlignHCenter)
            self.minEntries[self.integralIter].editingFinished.connect(lambda self=self, tmp=self.minEntries[self.integralIter]: self.setVal(tmp, True))
            self.frame3.addWidget(self.minEntries[self.integralIter], 2+self.entryCount, 0)
            self.maxEntries.append(QtGui.QLineEdit())
            self.maxEntries[self.integralIter].setAlignment(QtCore.Qt.AlignHCenter)
            self.maxEntries[self.integralIter].editingFinished.connect(lambda self=self, tmp=self.maxEntries[self.integralIter]: self.setVal(tmp, False))
            self.frame3.addWidget(self.maxEntries[self.integralIter], 2+self.entryCount, 1)
            self.intEntries.append(QtGui.QLineEdit())
            self.intEntries[self.integralIter].setAlignment(QtCore.Qt.AlignHCenter)
            self.intEntries[self.integralIter].editingFinished.connect(lambda self=self, tmp=self.intEntries[self.integralIter]: self.setRef(tmp))
            self.frame3.addWidget(self.intEntries[self.integralIter], 2+self.entryCount, 2)
            self.deleteButtons.append(QtGui.QPushButton("X"))
            self.deleteButtons[self.integralIter].clicked.connect(lambda extra, self=self, tmp=self.deleteButtons[self.integralIter]: self.deleteEntry(tmp))
            self.frame3.addWidget(self.deleteButtons[self.integralIter], 2+self.entryCount, 3)
            self.entryCount += 1
            self.first = True
            self.fit()

    def deleteEntry(self, button=None, num=None, reset=False):
        if num is None:
            num = self.deleteButtons.index(button)
        if num == self.integralIter:
            self.minValues[num] = min(self.xax)
            self.maxValues[num] = max(self.xax)
            self.intValues[num] = 1.0
            self.minEntries[num].setText("")
            self.maxEntries[num].setText("")
            self.first = True
            return
        self.frame3.removeWidget(self.maxEntries[num])
        self.frame3.removeWidget(self.minEntries[num])
        self.frame3.removeWidget(self.deleteButtons[num])
        self.frame3.removeWidget(self.intEntries[num])
        self.maxEntries[num].deleteLater()
        self.minEntries[num].deleteLater()
        self.deleteButtons[num].deleteLater()
        self.intEntries[num].deleteLater()
        self.maxEntries.pop(num)
        self.minEntries.pop(num)
        self.deleteButtons.pop(num)
        self.intEntries.pop(num)
        self.minValues = np.delete(self.minValues, num)
        self.maxValues = np.delete(self.maxValues, num)
        self.intValues = np.delete(self.intValues, num)
        self.integralIter -= 1
        if not reset:
            self.fit()
            
    def togglePick(self):
        self.parent.togglePick(self.pickTick.isChecked())

    def setVal(self, entry, isMin=False):
        inp = safeEval(entry.text())
        if inp is None:
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        if inp < min(self.xax):
            inp = min(self.xax)
        if inp > max(self.xax):
            inp = max(self.xax) 
        if isMin:
            num = self.minEntries.index(entry)
            self.minValues[num] = min(inp, self.maxValues[num])
            self.maxValues[num] = max(inp, self.maxValues[num])
        else:
            num = self.maxEntries.index(entry)
            self.maxValues[num] = max(inp, self.minValues[num])
            self.minValues[num] = min(inp, self.minValues[num])
        if num == self.integralIter:
            self.integralIter += 1
            self.minValues = np.append(self.minValues, min(self.xax))
            self.maxValues = np.append(self.maxValues, max(self.xax))
            self.intValues = np.append(self.intValues, 1)
            self.minEntries.append(QtGui.QLineEdit())
            self.minEntries[self.integralIter].setAlignment(QtCore.Qt.AlignHCenter)
            self.minEntries[self.integralIter].editingFinished.connect(lambda self=self, tmp=self.minEntries[self.integralIter]: self.setVal(tmp, True))
            self.frame3.addWidget(self.minEntries[self.integralIter], 2+self.entryCount, 0)
            self.maxEntries.append(QtGui.QLineEdit())
            self.maxEntries[self.integralIter].setAlignment(QtCore.Qt.AlignHCenter)
            self.maxEntries[self.integralIter].editingFinished.connect(lambda self=self, tmp=self.maxEntries[self.integralIter]: self.setVal(tmp, False))
            self.frame3.addWidget(self.maxEntries[self.integralIter], 2+self.entryCount, 1)
            self.intEntries.append(QtGui.QLineEdit())
            self.intEntries[self.integralIter].setAlignment(QtCore.Qt.AlignHCenter)
            self.intEntries[self.integralIter].editingFinished.connect(lambda self=self, tmp=self.intEntries[self.integralIter]: self.setRef(tmp))
            self.frame3.addWidget(self.intEntries[self.integralIter], 2+self.entryCount, 2)
            self.deleteButtons.append(QtGui.QPushButton("X"))
            self.deleteButtons[self.integralIter].clicked.connect(lambda extra, self=self, tmp=self.deleteButtons[self.integralIter]: self.deleteEntry(tmp))
            self.frame3.addWidget(self.deleteButtons[self.integralIter], 2+self.entryCount, 3)
            self.entryCount += 1
            self.first = True
        self.minEntries[num].setText("%#.3g" %self.minValues[num])
        self.maxEntries[num].setText("%#.3g" %self.maxValues[num])
        self.fit()
            
    def setRef(self, entry):
        num = self.intEntries.index(entry)
        inp = safeEval(entry.text())
        if inp is None:
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        self.refVal =  self.intValues[num] / float(inp)
        self.fit()
    
    def displayInt(self):
        if self.refVal is None:
            self.refVal = self.intValues[0]
        if self.absIntTick.isChecked():
            tmpInts = self.intValues/float(self.refVal)
        else:
            tmpInts = self.intValues
        for i in range(self.integralIter):
            self.intEntries[i].setText("%#.3g" % tmpInts[i])
    
    def fit(self, *args):
        x = []
        y = []
        if self.integralIter == 0:
            self.parent.showPlot()
            return
        for i in range(self.integralIter):
            tmpx = self.parent.current.xax[(self.minValues[i]<self.xax)&(self.maxValues[i]>self.xax)]
            if self.parent.current.plotType == 0:
                tmpy = np.real(self.parent.current.data1D[(self.minValues[i]<self.xax)&(self.maxValues[i]>self.xax)])
            elif self.parent.current.plotType == 1:
                tmpy = np.imag(self.parent.current.data1D[(self.minValues[i]<self.xax)&(self.maxValues[i]>self.xax)])
            elif self.parent.current.plotType == 2:
                tmpy = np.real(self.parent.current.data1D[(self.minValues[i]<self.xax)&(self.maxValues[i]>self.xax)])
            elif self.parent.current.plotType == 3:
                tmpy = np.abs(self.parent.current.data1D[(self.minValues[i]<self.xax)&(self.maxValues[i]>self.xax)])
            self.intValues[i] = np.sum(tmpy)*self.parent.current.sw/float(self.parent.current.data1D.shape[-1])
            if self.parent.spec == 1:
                x = np.append(x, tmpx[::-1])
                y = np.append(y, np.cumsum(tmpy[::-1]))
            else:
                x = np.append(x, tmpx)
                y = np.append(y, np.cumsum(tmpy))
            x = np.append(x, float('nan'))
            y = np.append(y, float('nan'))
        self.displayInt()
        y = y / (max(y)-min(y)) * self.diffy
        self.parent.showPlot(x, y)
        
##############################################################################
class RelaxWindow(FittingWindow): 
    def __init__(self, mainProgram, oldMainWindow):
        FittingWindow.__init__(self, mainProgram, oldMainWindow)

    def setup(self):
        self.current = RelaxFrame(self, self.fig, self.canvas, self.oldMainWindow.current)
        self.paramframe = RelaxParamFrame(self.current, self)
        
#################################################################################   
class RelaxFrame(Plot1DFrame): 
    def __init__(self, rootwindow, fig, canvas, current):
        self.ref = current.ref
        self.axType = current.axType
        self.freq = current.freq
        self.xax = current.xax
        self.data1D = current.getDisplayedData()
        self.plotType = 0
        self.logx = 0
        self.logy = 0
        Plot1DFrame.__init__(self, rootwindow, fig, canvas)
        self.current = current
        self.rootwindow = rootwindow
        self.plotReset()
        self.showPlot()

    def plotReset(self): 
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
        differ = 0.05*(maxy-miny) 
        self.yminlim = miny-differ
        self.ymaxlim = maxy+differ
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.xminlim = min(self.xax*axMult)
        self.xmaxlim = max(self.xax*axMult)
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)

    def showPlot(self, tmpAx=None, tmpdata=None): 
        self.ax.cla()
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if tmpAx is not None:
            self.ax.plot(tmpAx*axMult, tmpdata)
        self.ax.scatter(self.xax*axMult, self.data1D)
        if self.logx == 0:
            self.ax.set_xscale('linear')
        else:
           self.ax.set_xscale('log')
        if self.logy == 0:
            self.ax.set_yscale('linear')
        else:
           self.ax.set_yscale('log') 
        if self.spec == 0:
            if self.current.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.current.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.current.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.current.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.current.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        if self.logx == 0:
            self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.logy == 0:
            self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.canvas.draw()
        
    def scroll(self, event):
        if self.rightMouse:
            if self.logx == 0:
                middle = (self.xmaxlim+self.xminlim)/2.0
                width = self.xmaxlim-self.xminlim
                width = width*0.9**event.step
                self.xmaxlim = middle+width/2.0
                self.xminlim = middle-width/2.0
                if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                    self.ax.set_xlim(self.xmaxlim, self.xminlim)
                else:
                    self.ax.set_xlim(self.xminlim, self.xmaxlim)
            else:
                middle = (np.log(self.xmaxlim)+np.log(self.xminlim))/2.0
                width = np.log(self.xmaxlim)-np.log(self.xminlim)
                width = width*0.9**event.step
                self.xmaxlim = np.exp(middle+width/2.0)
                self.xminlim = np.exp(middle-width/2.0)
                if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                    self.ax.set_xlim(self.xmaxlim, self.xminlim)
                else:
                    self.ax.set_xlim(self.xminlim, self.xmaxlim)
        else:
            if self.logy == 0:
                middle = (self.ymaxlim+self.yminlim)/2.0
                width = self.ymaxlim-self.yminlim
                width = width*0.9**event.step
                self.ymaxlim = middle+width/2.0
                self.yminlim = middle-width/2.0
                if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
                    self.ax.set_ylim(self.ymaxlim, self.yminlim)
                else:
                    self.ax.set_ylim(self.yminlim, self.ymaxlim)
            else:
                middle = (np.log(self.ymaxlim)+np.log(self.yminlim))/2.0
                width = np.log(self.ymaxlim)-np.log(self.yminlim)
                width = width*0.9**event.step
                self.ymaxlim = np.exp(middle+width/2.0)
                self.yminlim = np.exp(middle-width/2.0)
                if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
                    self.ax.set_ylim(self.ymaxlim, self.yminlim)
                else:
                    self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()
        
    def buttonRelease(self, event):
        if event.button == 1:
            if self.peakPick:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                    self.rect[0] = None
                    self.peakPick = False
                    idx = np.argmin(np.abs(self.line_xdata-event.xdata))
                    if self.peakPickFunc is not None:
                        self.peakPickFunc((idx, self.line_xdata[idx], self.line_ydata[idx]))
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
                    if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                        self.ax.set_xlim(self.xmaxlim, self.xminlim)
                    else:
                        self.ax.set_xlim(self.xminlim, self.xmaxlim)
                    if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
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

    def pan(self, event):
        if self.rightMouse and self.panX is not None and self.panY is not None:
            if self.logx == 0 and self.logy == 0:
                inv = self.ax.transData.inverted()
                point = inv.transform((event.x, event.y))
                x = point[0]
                y = point[1]
            else:
                x = event.xdata
                y = event.ydata
                if x is None or y is None:
                    return
            if self.logx == 0:
                diffx = x-self.panX
                self.xmaxlim = self.xmaxlim-diffx
                self.xminlim = self.xminlim-diffx
            else:
                diffx = np.log(x)-np.log(self.panX)
                self.xmaxlim = np.exp(np.log(self.xmaxlim)-diffx)
                self.xminlim = np.exp(np.log(self.xminlim)-diffx)
            if self.logy == 0:
                diffy = y-self.panY
                self.ymaxlim = self.ymaxlim-diffy
                self.yminlim = self.yminlim-diffy
            else:
                diffy = np.log(y)-np.log(self.panY)
                self.ymaxlim = np.exp(np.log(self.ymaxlim)-diffy)
                self.yminlim = np.exp(np.log(self.yminlim)-diffy)
            if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                self.ax.set_xlim(self.xmaxlim, self.xminlim)
            else:
                self.ax.set_xlim(self.xminlim, self.xmaxlim)
            if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
                self.ax.set_ylim(self.ymaxlim, self.yminlim)
            else:
                self.ax.set_ylim(self.yminlim, self.ymaxlim)
            self.canvas.draw()
        elif self.peakPick:
            if self.rect[0] is not None:
                self.rect[0].remove()
                self.rect[0] = None
            if event.xdata is not None:
                self.rect[0] = self.ax.axvline(event.xdata, c='k', linestyle='--')
            self.canvas.draw()
        elif self.leftMouse and (self.zoomX1 is not None) and (self.zoomY1 is not None):
            if self.logx == 0 and self.logy == 0:
                inv = self.ax.transData.inverted()
                point = inv.transform((event.x, event.y))
                self.zoomX2 = point[0]
                self.zoomY2 = point[1]
            else:
                self.zoomX2 = event.xdata
                self.zoomY2 = event.ydata
                if self.zoomX2 is None or self.zoomY2 is None:
                    return
            if self.rect[0] is not None:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                if self.rect[1] is not None:
                    self.rect[1].remove()
                if self.rect[2] is not None:
                    self.rect[2].remove()
                if self.rect[3] is not None:
                    self.rect[3].remove()
                self.rect = [None, None, None, None]
            self.rect[0], = self.ax.plot([self.zoomX1, self.zoomX2], [self.zoomY2, self.zoomY2], 'k', clip_on= False)
            self.rect[1], = self.ax.plot([self.zoomX1, self.zoomX2], [self.zoomY1, self.zoomY1], 'k', clip_on=False)
            self.rect[2], = self.ax.plot([self.zoomX1, self.zoomX1], [self.zoomY1, self.zoomY2], 'k', clip_on=False)
            self.rect[3], = self.ax.plot([self.zoomX2, self.zoomX2], [self.zoomY1, self.zoomY2], 'k', clip_on=False)
            self.canvas.draw()
            
    def setLog(self, logx, logy):
        self.logx = logx
        self.logy = logy
        if self.logx == 0:
            self.ax.set_xscale('linear')
        else:
           self.ax.set_xscale('log')
        if self.logy == 0:
            self.ax.set_yscale('linear')
        else:
           self.ax.set_yscale('log') 
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()

#################################################################################
class RelaxParamFrame(QtGui.QWidget): 
    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1, 0, 0)
        grid.addLayout(self.optframe, 0, 1)
        grid.addLayout(self.frame2, 0, 2)
        grid.addLayout(self.frame3, 0, 3)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton, 0, 0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton, 1, 0)
        fitAllButton = QtGui.QPushButton("Fit all")
        fitAllButton.clicked.connect(self.fitAll)
        self.frame1.addWidget(fitAllButton, 2, 0)
        copyResultButton = QtGui.QPushButton("Copy result")
        copyResultButton.clicked.connect(lambda: self.sim(True))
        self.frame1.addWidget(copyResultButton, 3, 0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton, 4, 0)
        self.frame1.setColumnStretch(10, 1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel("Amplitude:"), 0, 0, 1, 2)
        self.ampTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.ampTick, 1, 0)
        self.ampEntry = QtGui.QLineEdit()
        self.ampEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ampEntry.setText("%#.3g" % np.amax(self.parent.data1D))
        self.frame2.addWidget(self.ampEntry, 1, 1)
        self.frame2.addWidget(QLabel("Constant:"), 2, 0, 1, 2)
        self.constTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.constTick, 3, 0)
        self.constEntry = QtGui.QLineEdit()
        self.constEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.constEntry.setText("1.0")
        self.frame2.addWidget(self.constEntry, 3, 1)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1', '2', '3', '4'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        self.frame3.addWidget(QLabel("Coefficient:"), 1, 0, 1, 2)
        self.frame3.addWidget(QLabel("T [s]:"), 1, 2, 1, 2)
        self.frame3.setColumnStretch(10, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.xlog = QtGui.QCheckBox('x-log')
        self.xlog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.xlog, 0, 0, QtCore.Qt.AlignTop)
        self.ylog = QtGui.QCheckBox('y-log')
        self.ylog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.ylog, 1, 0, QtCore.Qt.AlignTop)
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        self.coeffTicks = []
        self.coeffEntries = []
        self.t1Ticks = []
        self.t1Entries = []
        for i in range(4):
            self.coeffTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.coeffTicks[i], i+2, 0)
            self.coeffEntries.append(QtGui.QLineEdit())
            self.coeffEntries[i].setText("-1.0")
            self.coeffEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.coeffEntries[i], i+2, 1)
            self.t1Ticks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.t1Ticks[i], i+2, 2)
            self.t1Entries.append(QtGui.QLineEdit())
            self.t1Entries[i].setText("1.0")
            self.t1Entries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.t1Entries[i], i+2, 3)
            if i > 0:
                self.coeffTicks[i].hide()
                self.coeffEntries[i].hide()
                self.t1Ticks[i].hide()
                self.t1Entries[i].hide()

    def setLog(self, *args):
        self.parent.setLog(self.xlog.isChecked(), self.ylog.isChecked())
                
    def changeNum(self, *args):
        val = self.numExp.currentIndex()+1
        for i in range(4):
            if i < val:
                self.coeffTicks[i].show()
                self.coeffEntries[i].show()
                self.t1Ticks[i].show()
                self.t1Entries[i].show()
            else:
                self.coeffTicks[i].hide()
                self.coeffEntries[i].hide()
                self.t1Ticks[i].hide()
                self.t1Entries[i].hide()

    def fitFunc(self, x, *param):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        testFunc = np.zeros(len(x))
        if struc[0]:
            amplitude = param[0]
            param = np.delete(param, [0])
        else:
            amplitude = argu[0]
            argu = np.delete(argu, [0])
        if struc[1]:
            constant = param[0]
            param = np.delete(param, [0])
        else:
            constant = argu[0]
            argu = np.delete(argu, [0])
        for i in range(1, numExp+1):
            if struc[2*i]:
                coeff = param[0]
                param = np.delete(param, [0])
            else:
                coeff = argu[0]
                argu = np.delete(argu, [0])
            if struc[2*i+1]:
                T1 = param[0]
                param = np.delete(param, [0])
            else:
                T1 = argu[0]
                argu = np.delete(argu, [0])
            testFunc += coeff*np.exp(-1.0*x/T1) 
        return amplitude*(constant+testFunc)

    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.ampEntry.text())
        if inp is None:
            return False
        self.ampEntry.setText('%#.3g' % inp)
        inp = safeEval(self.constEntry.text())
        if inp is None:
            return False
        self.constEntry.setText('%#.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.coeffEntries[i].text())
            if inp is None:
                return False
            self.coeffEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.t1Entries[i].text())
            if inp is None:
                return False
            self.t1Entries[i].setText('%#.3g' % inp)
        return True
    
    def fit(self, *args):
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outCoeff = np.zeros(numExp)
        outT1 = np.zeros(numExp)
        if not self.ampTick.isChecked():
            guess.append(float(self.ampEntry.text()))
            struc.append(True)
        else:
            outAmp = float(self.ampEntry.text())
            argu.append(outAmp)
            struc.append(False)
        if not self.constTick.isChecked():
            guess.append(float(self.constEntry.text()))
            struc.append(True)
        else:
            outConst = float(self.constEntry.text())
            argu.append(outConst)
            struc.append(False)
        for i in range(numExp):
            if not self.coeffTicks[i].isChecked():
                guess.append(float(self.coeffEntries[i].text()))
                struc.append(True)
            else:
                outCoeff[i] = float(self.coeffEntries[i].text())
                argu.append(outCoeff[i])
                struc.append(False)
            if not self.t1Ticks[i].isChecked():
                guess.append(float(self.t1Entries[i].text()))
                struc.append(True)
            else:
                outT1[i] = safeEval(self.t1Entries[i].text())
                argu.append(outT1[i])
                struc.append(False)
        self.args = (numExp, struc, argu)
        fitVal = scipy.optimize.curve_fit(self.fitFunc, self.parent.xax, self.parent.data1D, guess)
        counter = 0
        if struc[0]:
            self.ampEntry.setText('%#.3g' % fitVal[0][counter])
            outAmp = fitVal[0][counter]
            counter += 1
        if struc[1]:
            self.constEntry.setText('%#.3g' % fitVal[0][counter])
            outConst = fitVal[0][counter]
            counter += 1
        for i in range(1, numExp+1):
            if struc[2*i]:
                self.coeffEntries[i-1].setText('%#.3g' % fitVal[0][counter])
                outCoeff[i-1] = fitVal[0][counter]
                counter += 1
            if struc[2*i+1]:
                self.t1Entries[i-1].setText('%#.3g' % fitVal[0][counter])
                outT1[i-1] = fitVal[0][counter]
                counter += 1
        self.disp(outAmp, outConst, outCoeff, outT1)

    def fitAll(self, *args):
        FitAllSelectionWindow(self, ["Amplitude", "Constant", "Coefficient", "T"])
        
    def fitAllFunc(self, outputs):
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outCoeff = np.zeros(numExp)
        outT1 = np.zeros(numExp)
        if not self.ampTick.isChecked():
            guess.append(float(self.ampEntry.text()))
            struc.append(True)
        else:
            outAmp = float(self.ampEntry.text())
            argu.append(outAmp)
            struc.append(False)
        if not self.constTick.isChecked():
            guess.append(float(self.constEntry.text()))
            struc.append(True)
        else:
            outConst = float(self.constEntry.text())
            argu.append(outConst)
            struc.append(False)
        for i in range(numExp):
            if not self.coeffTicks[i].isChecked():
                guess.append(float(self.coeffEntries[i].text()))
                struc.append(True)
            else:
                outCoeff[i] = float(self.coeffEntries[i].text())
                argu.append(outCoeff[i])
                struc.append(False)
            if not self.t1Ticks[i].isChecked():
                guess.append(float(self.t1Entries[i].text()))
                struc.append(True)
            else:
                outT1[i] = safeEval(self.t1Entries[i].text())
                argu.append(outT1[i])
                struc.append(False)
        self.args = (numExp, struc, argu)
        fullData = self.parent.current.data.data
        axes = self.parent.current.axes
        dataShape = fullData.shape
        dataShape2 = np.delete(dataShape, axes)
        rolledData = np.rollaxis(fullData, axes)
        intOutputs = np.array(outputs, dtype=int)
        numOutputs = np.sum(intOutputs[:2]) + numExp*np.sum(intOutputs[2:])
        outputData = np.zeros((np.product(dataShape2), numOutputs), dtype=complex)
        counter2 = 0
        for j in rolledData.reshape(dataShape[axes], np.product(dataShape2)).T:
            try:
                fitVal = scipy.optimize.curve_fit(self.fitFunc, self.parent.xax, np.real(j), guess)
            except:
                fitVal = [[0]*10]
            counter = 0
            if struc[0]:
                outAmp = fitVal[0][counter]
                counter += 1
            if struc[1]:
                outConst = fitVal[0][counter]
                counter += 1
            for i in range(1, numExp+1):
                if struc[2*i]:
                    outCoeff[i-1] = fitVal[0][counter]
                    counter += 1
                if struc[2*i+1]:
                    outT1[i-1] = fitVal[0][counter]
                    counter += 1
            outputArray = []
            if outputs[0]:
                outputArray = np.concatenate((outputArray, [outAmp]))
            if outputs[1]:
                outputArray = np.concatenate((outputArray, [outConst]))
            if outputs[2]:
                outputArray = np.concatenate((outputArray, outCoeff))
            if outputs[3]:
                outputArray = np.concatenate((outputArray, outT1))
            outputData[counter2] = outputArray
            counter2 += 1
        newShape = np.concatenate((np.array(dataShape2), [numOutputs]))
        self.rootwindow.createNewData(np.rollaxis(outputData.reshape(newShape), -1, axes), axes)
        
    def sim(self, store=False):
        numExp = self.numExp.currentIndex()+1
        outAmp = safeEval(self.ampEntry.text())
        outConst = safeEval(self.constEntry.text())
        if outAmp is None or outConst is None:
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        outCoeff = []
        outT1 = []
        for i in range(numExp):
            outCoeff.append(safeEval(self.coeffEntries[i].text()))
            outT1.append(safeEval(self.t1Entries[i].text()))
            if outCoeff[i] is None or outT1[i] is None:
                self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
                return
        if store:
            self.copyResult(outAmp, outConst, outCoeff, outT1)
        else:
            self.disp(outAmp, outConst, outCoeff, outT1)

    def copyResult(self, outAmp, outConst, outCoeff, outT1):
        outCurve = np.zeros((len(outCoeff)+2, len(self.parent.xax)))
        tmp = np.zeros(len(self.parent.xax))
        for i in range(len(outCoeff)):
            outCurve[i] = outAmp * (outConst + outCoeff[i]*np.exp(-self.parent.xax/outT1[i]))
            tmp += outCurve[i] 
        outCurve[len(outCoeff)] = tmp - (len(outCoeff)-1)*outAmp * outConst
        outCurve[len(outCoeff)+1] = self.parent.data1D
        self.rootwindow.createNewData(np.array(outCurve), self.parent.current.axes, True)
        
    def disp(self, outAmp, outConst, outCoeff, outT1):
        numCurve = 100 #number of points in output curve
        outCurve = np.zeros(numCurve)
        if self.xlog.isChecked():
            x = np.logspace(np.log(min(self.parent.xax)), np.log(max(self.parent.xax)), numCurve)
        else:
            x = np.linspace(min(self.parent.xax), max(self.parent.xax), numCurve)
        for i in range(len(outCoeff)):
            outCurve += outCoeff[i]*np.exp(-x/outT1[i])
        self.parent.showPlot(x, outAmp*(outConst+outCurve))

##############################################################################
class DiffusionWindow(FittingWindow): 
    def __init__(self, mainProgram, oldMainWindow):
        FittingWindow.__init__(self, mainProgram, oldMainWindow)

    def setup(self):
        self.current = DiffusionFrame(self, self.fig, self.canvas, self.oldMainWindow.current)
        self.paramframe = DiffusionParamFrame(self.current, self)
    
#################################################################################   
class DiffusionFrame(Plot1DFrame): 
    def __init__(self, rootwindow, fig, canvas, current):
        self.ref = current.ref
        self.axType = current.axType
        self.freq = current.freq
        self.xax = current.xax
        self.data1D = current.getDisplayedData()
        self.plotType = 0
        self.logx = 0
        self.logy = 0
        Plot1DFrame.__init__(self, rootwindow, fig, canvas)
        self.current = current
        self.rootwindow = rootwindow
        self.plotReset()
        self.showPlot()

    def plotReset(self): 
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
        differ = 0.05*(maxy-miny) 
        self.yminlim = miny-differ
        self.ymaxlim = maxy+differ
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.xminlim = min(self.xax*axMult)
        self.xmaxlim = max(self.xax*axMult)
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)

    def showPlot(self, tmpAx=None, tmpdata=None): 
        self.ax.cla()
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if tmpAx is not None:
            self.ax.plot(tmpAx*axMult, tmpdata)
        self.ax.scatter(self.xax*axMult, self.data1D)
        if self.logx == 0:
            self.ax.set_xscale('linear')
        else:
           self.ax.set_xscale('log')
        if self.logy == 0:
            self.ax.set_yscale('linear')
        else:
           self.ax.set_yscale('log') 
        if self.spec == 0:
            if self.current.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.current.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.current.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.current.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.current.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        if self.logx == 0:
            self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.logy == 0:
            self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.canvas.draw()
        
    def scroll(self, event):
        if self.rightMouse:
            if self.logx == 0:
                middle = (self.xmaxlim+self.xminlim)/2.0
                width = self.xmaxlim-self.xminlim
                width = width*0.9**event.step
                self.xmaxlim = middle+width/2.0
                self.xminlim = middle-width/2.0
                if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                    self.ax.set_xlim(self.xmaxlim, self.xminlim)
                else:
                    self.ax.set_xlim(self.xminlim, self.xmaxlim)
            else:
                middle = (np.log(self.xmaxlim)+np.log(self.xminlim))/2.0
                width = np.log(self.xmaxlim)-np.log(self.xminlim)
                width = width*0.9**event.step
                self.xmaxlim = np.exp(middle+width/2.0)
                self.xminlim = np.exp(middle-width/2.0)
                if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                    self.ax.set_xlim(self.xmaxlim, self.xminlim)
                else:
                    self.ax.set_xlim(self.xminlim, self.xmaxlim)
        else:
            if self.logy == 0:
                middle = (self.ymaxlim+self.yminlim)/2.0
                width = self.ymaxlim-self.yminlim
                width = width*0.9**event.step
                self.ymaxlim = middle+width/2.0
                self.yminlim = middle-width/2.0
                if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
                    self.ax.set_ylim(self.ymaxlim, self.yminlim)
                else:
                    self.ax.set_ylim(self.yminlim, self.ymaxlim)
            else:
                middle = (np.log(self.ymaxlim)+np.log(self.yminlim))/2.0
                width = np.log(self.ymaxlim)-np.log(self.yminlim)
                width = width*0.9**event.step
                self.ymaxlim = np.exp(middle+width/2.0)
                self.yminlim = np.exp(middle-width/2.0)
                if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
                    self.ax.set_ylim(self.ymaxlim, self.yminlim)
                else:
                    self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()
        
    def buttonRelease(self, event):
        if event.button == 1:
            if self.peakPick:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                    self.rect[0] = None
                    self.peakPick = False
                    idx = np.argmin(np.abs(self.line_xdata-event.xdata))
                    if self.peakPickFunc is not None:
                        self.peakPickFunc((idx, self.line_xdata[idx], self.line_ydata[idx]))
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
                    if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                        self.ax.set_xlim(self.xmaxlim, self.xminlim)
                    else:
                        self.ax.set_xlim(self.xminlim, self.xmaxlim)
                    if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
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

    def pan(self, event):
        if self.rightMouse and self.panX is not None and self.panY is not None:
            if self.logx == 0 and self.logy == 0:
                inv = self.ax.transData.inverted()
                point = inv.transform((event.x, event.y))
                x = point[0]
                y = point[1]
            else:
                x = event.xdata
                y = event.ydata
                if x is None or y is None:
                    return
            if self.logx == 0:
                diffx = x-self.panX
                self.xmaxlim = self.xmaxlim-diffx
                self.xminlim = self.xminlim-diffx
            else:
                diffx = np.log(x)-np.log(self.panX)
                self.xmaxlim = np.exp(np.log(self.xmaxlim)-diffx)
                self.xminlim = np.exp(np.log(self.xminlim)-diffx)
            if self.logy == 0:
                diffy = y-self.panY
                self.ymaxlim = self.ymaxlim-diffy
                self.yminlim = self.yminlim-diffy
            else:
                diffy = np.log(y)-np.log(self.panY)
                self.ymaxlim = np.exp(np.log(self.ymaxlim)-diffy)
                self.yminlim = np.exp(np.log(self.yminlim)-diffy)
            if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                self.ax.set_xlim(self.xmaxlim, self.xminlim)
            else:
                self.ax.set_xlim(self.xminlim, self.xmaxlim)
            if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
                self.ax.set_ylim(self.ymaxlim, self.yminlim)
            else:
                self.ax.set_ylim(self.yminlim, self.ymaxlim)
            self.canvas.draw()
        elif self.peakPick:
            if self.rect[0] is not None:
                self.rect[0].remove()
                self.rect[0] = None
            if event.xdata is not None:
                self.rect[0] = self.ax.axvline(event.xdata, c='k', linestyle='--')
            self.canvas.draw()
        elif self.leftMouse and (self.zoomX1 is not None) and (self.zoomY1 is not None):
            if self.logx == 0 and self.logy == 0:
                inv = self.ax.transData.inverted()
                point = inv.transform((event.x, event.y))
                self.zoomX2 = point[0]
                self.zoomY2 = point[1]
            else:
                self.zoomX2 = event.xdata
                self.zoomY2 = event.ydata
                if self.zoomX2 is None or self.zoomY2 is None:
                    return
            if self.rect[0] is not None:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                if self.rect[1] is not None:
                    self.rect[1].remove()
                if self.rect[2] is not None:
                    self.rect[2].remove()
                if self.rect[3] is not None:
                    self.rect[3].remove()
                self.rect = [None, None, None, None]
            self.rect[0], = self.ax.plot([self.zoomX1, self.zoomX2], [self.zoomY2, self.zoomY2], 'k', clip_on=False)
            self.rect[1], = self.ax.plot([self.zoomX1, self.zoomX2], [self.zoomY1, self.zoomY1], 'k', clip_on=False)
            self.rect[2], = self.ax.plot([self.zoomX1, self.zoomX1], [self.zoomY1, self.zoomY2], 'k', clip_on=False)
            self.rect[3], = self.ax.plot([self.zoomX2, self.zoomX2], [self.zoomY1, self.zoomY2], 'k', clip_on=False)
            self.canvas.draw()
            
    def setLog(self, logx, logy):
        self.logx = logx
        self.logy = logy
        if self.logx == 0:
            self.ax.set_xscale('linear')
        else:
           self.ax.set_xscale('log')
        if self.logy == 0:
            self.ax.set_yscale('linear')
        else:
           self.ax.set_yscale('log') 
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()

#################################################################################
class DiffusionParamFrame(QtGui.QWidget): 
    def __init__(self, parent, rootwindow): 
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        self.frame4 = QtGui.QGridLayout()
        grid.addLayout(self.frame1, 0, 0)
        grid.addLayout(self.optframe, 0, 1)
        grid.addLayout(self.frame2, 0, 2)
        grid.addLayout(self.frame3, 0, 3)
        grid.addLayout(self.frame4, 0, 4)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton, 0, 0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton, 1, 0)
        fitAllButton = QtGui.QPushButton("Fit all")
        fitAllButton.clicked.connect(self.fitAll)
        self.frame1.addWidget(fitAllButton, 2, 0)
        copyResultButton = QtGui.QPushButton("Copy result")
        copyResultButton.clicked.connect(lambda: self.sim(True))
        self.frame1.addWidget(copyResultButton, 3, 0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton, 4, 0)
        self.frame1.setColumnStretch(10, 1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel(u"\u03b3 [MHz/T]:"), 0, 0)
        self.gammaEntry = QtGui.QLineEdit()
        self.gammaEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.gammaEntry.setText("42.576")
        self.frame2.addWidget(self.gammaEntry, 1, 0)
        self.frame2.addWidget(QLabel(u"\u03b4 [s]:"), 2, 0)
        self.deltaEntry = QtGui.QLineEdit()
        self.deltaEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.deltaEntry.setText("1.0")
        self.frame2.addWidget(self.deltaEntry, 3, 0)
        self.frame2.addWidget(QLabel(u"\u0394 [s]:"), 4, 0)
        self.triangleEntry = QtGui.QLineEdit()
        self.triangleEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.triangleEntry.setText("1.0")
        self.frame2.addWidget(self.triangleEntry, 5, 0)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.frame3.addWidget(QLabel("Amplitude:"), 0, 0, 1, 2)
        self.ampTick = QtGui.QCheckBox('')
        self.frame3.addWidget(self.ampTick, 1, 0)
        self.ampEntry = QtGui.QLineEdit()
        self.ampEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ampEntry.setText("%#.3g" % np.amax(self.parent.data1D))
        self.frame3.addWidget(self.ampEntry, 1, 1)
        self.frame3.addWidget(QLabel("Constant:"), 2, 0, 1, 2)
        self.constTick = QtGui.QCheckBox('')
        self.constTick.setChecked(True)
        self.frame3.addWidget(self.constTick, 3, 0)
        self.constEntry = QtGui.QLineEdit()
        self.constEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.constEntry.setText("0.0")
        self.frame3.addWidget(self.constEntry, 3, 1)
        self.frame3.setColumnStretch(10, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1', '2', '3', '4'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame4.addWidget(self.numExp, 0, 0, 1, 2)
        self.frame4.addWidget(QLabel("Coefficient:"), 1, 0, 1, 2)
        self.frame4.addWidget(QLabel("D [m^2/s]:"), 1, 2, 1, 2)
        self.frame4.setColumnStretch(20, 1)
        self.frame4.setAlignment(QtCore.Qt.AlignTop)
        self.xlog = QtGui.QCheckBox('x-log')
        self.xlog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.xlog, 0, 0)
        self.ylog = QtGui.QCheckBox('y-log')
        self.ylog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.ylog, 1, 0)
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        self.coeffEntries = []
        self.coeffTicks = []
        self.dEntries = []
        self.dTicks = []
        for i in range(4):
            self.coeffTicks.append(QtGui.QCheckBox(''))
            self.frame4.addWidget(self.coeffTicks[i], i+2, 0)
            self.coeffEntries.append(QtGui.QLineEdit())
            self.coeffEntries[i].setText("1.0")
            self.coeffEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame4.addWidget(self.coeffEntries[i], i+2, 1)
            self.dTicks.append(QtGui.QCheckBox(''))
            self.frame4.addWidget(self.dTicks[i], i+2, 2)
            self.dEntries.append(QtGui.QLineEdit())
            self.dEntries[i].setText("1.0e-9")
            self.dEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame4.addWidget(self.dEntries[i], i+2, 3)
            if i > 0:
                self.coeffTicks[i].hide()
                self.coeffEntries[i].hide()
                self.dTicks[i].hide()
                self.dEntries[i].hide()
                
    def setLog(self, *args):
        self.parent.setLog(self.xlog.isChecked(), self.ylog.isChecked())
                
    def changeNum(self, *args):
        val = self.numExp.currentIndex()+1
        for i in range(4):
            if i < val:
                self.coeffTicks[i].show()
                self.coeffEntries[i].show()
                self.dTicks[i].show()
                self.dEntries[i].show()
            else:
                self.coeffTicks[i].hide()
                self.coeffEntries[i].hide()
                self.dTicks[i].hide()
                self.dEntries[i].hide()

    def fitFunc(self, x, *param):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        gamma = self.args[3]
        delta = self.args[4]
        triangle = self.args[5]
        testFunc = np.zeros(len(self.parent.data1D))
        if struc[0]:
            amplitude = param[0]
            param = np.delete(param, [0])
        else:
            amplitude = argu[0]
            argu = np.delete(argu, [0])
        if struc[1]:
            constant = param[0]
            param = np.delete(param, [0])
        else:
            constant = argu[0]
            argu = np.delete(argu, [0])
        for i in range(1, numExp+1):
            if struc[2*i]:
                coeff = param[0]
                param = np.delete(param, [0])
            else:
                coeff = argu[0]
                argu = np.delete(argu, [0])
            if struc[2*i+1]:
                D = param[0]
                param = np.delete(param, [0])
            else:
                D = argu[0]
                argu = np.delete(argu, [0])
            testFunc += coeff*np.exp(-(gamma*delta*x)**2*D*(triangle-delta/3.0)) 
        return amplitude*(constant+testFunc)
    
    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.gammaEntry.text())
        if inp is None:
            return False
        self.gammaEntry.setText('%#.3g' % inp)
        inp = safeEval(self.deltaEntry.text())
        if inp is None:
            return False
        self.deltaEntry.setText('%#.3g' % inp)
        inp = safeEval(self.triangleEntry.text())
        if inp is None:
            return False
        self.triangleEntry.setText('%#.3g' % inp)
        inp = safeEval(self.ampEntry.text())
        if inp is None:
            return False
        self.ampEntry.setText('%#.3g' % inp)
        inp = safeEval(self.constEntry.text())
        if inp is None:
            return False
        self.constEntry.setText('%#.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.coeffEntries[i].text())
            if inp is None:
                return False
            self.coeffEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.dEntries[i].text())
            if inp is None:
                return False
            self.dEntries[i].setText('%#.3g' % inp)
        return True
    
    def fit(self, *args):
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outCoeff = np.zeros(numExp)
        outD = np.zeros(numExp)
        gamma = float(self.gammaEntry.text())
        delta = float(self.deltaEntry.text())
        triangle = float(self.triangleEntry.text())
        if not self.ampTick.isChecked():
            guess.append(float(self.ampEntry.text()))
            struc.append(True)
        else:
            outAmp = float(self.ampEntry.text())
            argu.append(outAmp)
            struc.append(False)
        if not self.constTick.isChecked():
            guess.append(float(self.constEntry.text()))
            struc.append(True)
        else:
            outConst = float(self.constEntry.text())
            argu.append(outConst)
            struc.append(False)
        for i in range(numExp):
            if not self.coeffTicks[i].isChecked():
                guess.append(float(self.coeffEntries[i].text()))
                struc.append(True)
            else:
                outCoeff[i] = float(self.coeffEntries[i].text())
                argu.append(outCoeff[i])
                struc.append(False)
            if not self.dTicks[i].isChecked():
                guess.append(float(self.dEntries[i].text()))
                struc.append(True)
            else:
                outD[i] = float(self.dEntries[i].text())
                argu.append(outD[i])
                struc.append(False)
        self.args = (numExp, struc, argu, gamma, delta, triangle)
        fitVal = scipy.optimize.curve_fit(self.fitFunc, self.parent.xax, self.parent.data1D, guess)
        counter = 0
        if struc[0]:
            self.ampEntry.setText('%#.3g' % fitVal[0][counter])
            outAmp = fitVal[0][counter]
            counter += 1
        if struc[1]:
            self.constEntry.setText('%#.3g' % fitVal[0][counter])
            outConst = fitVal[0][counter]
            counter += 1
        for i in range(1, numExp+1):
            if struc[2*i]:
                self.coeffEntries[i-1].setText('%#.3g' % fitVal[0][counter])
                outCoeff[i-1] = fitVal[0][counter]
                counter += 1
            if struc[2*i+1]:
                self.dEntries[i-1].setText('%#.3g' % fitVal[0][counter])
                outD[i-1] = fitVal[0][counter]
                counter += 1
        self.disp(outAmp, outConst, outCoeff, outD, gamma, delta, triangle)

    def fitAll(self, *args):
        FitAllSelectionWindow(self, ["Amplitude", "Constant", "Coefficient", "D"])
        
    def fitAllFunc(self, outputs):
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outCoeff = np.zeros(numExp)
        outD = np.zeros(numExp)
        gamma = float(self.gammaEntry.text())
        delta = float(self.deltaEntry.text())
        triangle = float(self.triangleEntry.text())
        if not self.ampTick.isChecked():
            guess.append(float(self.ampEntry.text()))
            struc.append(True)
        else:
            outAmp = float(self.ampEntry.text())
            argu.append(outAmp)
            struc.append(False)
        if not self.constTick.isChecked():
            guess.append(float(self.constEntry.text()))
            struc.append(True)
        else:
            outConst = float(self.constEntry.text())
            argu.append(outConst)
            struc.append(False)
        for i in range(numExp):
            if not self.coeffTicks[i].isChecked():
                guess.append(float(self.coeffEntries[i].text()))
                struc.append(True)
            else:
                outCoeff[i] = float(self.coeffEntries[i].text())
                argu.append(outCoeff[i])
                struc.append(False)
            if not self.dTicks[i].isChecked():
                guess.append(float(self.dEntries[i].text()))
                struc.append(True)
            else:
                outD[i] = safeEval(self.dEntries[i].text())
                argu.append(outD[i])
                struc.append(False)
        self.args = (numExp, struc, argu, gamma, delta, triangle)
        fullData = self.parent.current.data.data
        axes = self.parent.current.axes
        dataShape = fullData.shape
        dataShape2 = np.delete(dataShape, axes)
        rolledData = np.rollaxis(fullData, axes)
        intOutputs = np.array(outputs, dtype=int)
        numOutputs = np.sum(intOutputs[:2]) + numExp*np.sum(intOutputs[2:])
        outputData = np.zeros((np.product(dataShape2), numOutputs), dtype=complex)
        counter2 = 0
        for j in rolledData.reshape(dataShape[axes], np.product(dataShape2)).T:
            try:
                fitVal = scipy.optimize.curve_fit(self.fitFunc, self.parent.xax, np.real(j), guess)
            except:
                fitVal = [[0]*10]
            counter = 0
            if struc[0]:
                outAmp = fitVal[0][counter]
                counter += 1
            if struc[1]:
                outConst = fitVal[0][counter]
                counter += 1
            for i in range(1, numExp+1):
                if struc[2*i]:
                    outCoeff[i-1] = fitVal[0][counter]
                    counter += 1
                if struc[2*i+1]:
                    outD[i-1] = fitVal[0][counter]
                    counter += 1
            outputArray = []
            if outputs[0]:
                outputArray = np.concatenate((outputArray, [outAmp]))
            if outputs[1]:
                outputArray = np.concatenate((outputArray, [outConst]))
            if outputs[2]:
                outputArray = np.concatenate((outputArray, outCoeff))
            if outputs[3]:
                outputArray = np.concatenate((outputArray, outD))
            outputData[counter2] = outputArray
            counter2 += 1
        newShape = np.concatenate((np.array(dataShape2), [numOutputs]))
        self.rootwindow.createNewData(np.rollaxis(outputData.reshape(newShape), -1, axes), axes)
        
    def sim(self, store=False):
        numExp = self.numExp.currentIndex()+1
        outAmp = safeEval(self.ampEntry.text())
        outConst = safeEval(self.constEntry.text())
        gamma = safeEval(self.gammaEntry.text())
        delta = safeEval(self.deltaEntry.text())
        triangle = safeEval(self.triangleEntry.text())
        if not np.isfinite([outAmp, outConst, gamma, delta, triangle]).all():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        outCoeff = []
        outD = []
        for i in range(numExp):
            outCoeff.append(safeEval(self.coeffEntries[i].text()))
            outD.append(safeEval(self.dEntries[i].text()))
            if outCoeff[i] is None or outD[i] is None:
                self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
                return
        if store:
            self.copyResult(outAmp, outConst, outCoeff, outD, gamma, delta, triangle)
        else: 
            self.disp(outAmp, outConst, outCoeff, outD, gamma, delta, triangle)

    def copyResult(self, outAmp, outConst, outCoeff, outD, gamma, delta, triangle):
        outCurve = np.zeros((len(outCoeff)+2, len(self.parent.xax)))
        tmp = np.zeros(len(self.parent.xax))
        for i in range(len(outCoeff)):
            outCurve[i] = outCoeff[i]*np.exp(-(gamma*delta*self.parent.xax)**2*outD[i]*(triangle-delta/3.0))
            tmp += outCurve[i]
        outCurve[len(outCoeff)] = tmp
        outCurve[len(outCoeff)+1] = self.parent.data1D
        self.rootwindow.createNewData(np.array(outCurve), self.parent.current.axes, True)
        
    def disp(self, outAmp, outConst, outCoeff, outD, gamma, delta, triangle):
        numCurve = 100 
        outCurve = np.zeros(numCurve)
        if self.xlog.isChecked():
            x = np.logspace(np.log(min(self.parent.xax)), np.log(max(self.parent.xax)), numCurve)
        else:
            x = np.linspace(min(self.parent.xax), max(self.parent.xax), numCurve)
        for i in range(len(outCoeff)):
            outCurve += outCoeff[i]*np.exp(-(gamma*delta*x)**2*outD[i]*(triangle-delta/3.0))
        self.parent.showPlot(x, outAmp*(outConst+outCurve))
        
##############################################################################
class PeakDeconvWindow(FittingWindow): 
    def __init__(self, mainProgram, oldMainWindow):
        FittingWindow.__init__(self, mainProgram, oldMainWindow)

    def setup(self):
        self.current = PeakDeconvFrame(self, self.fig, self.canvas, self.oldMainWindow.current)
        self.paramframe = PeakDeconvParamFrame(self.current, self)

#################################################################################   
class PeakDeconvFrame(Plot1DFrame): 
    def __init__(self, rootwindow, fig, canvas, current):
        Plot1DFrame.__init__(self, rootwindow, fig, canvas)
        self.data1D = current.getDisplayedData()
        self.current = current
        self.spec = self.current.spec
        self.xax = self.current.xax
        self.plotType = 0
        self.rootwindow = rootwindow
        self.pickNum = 0
        self.pickWidth = False
        self.plotReset()
        self.showPlot()

    def plotReset(self): 
        a = self.fig.gca()
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
        differ = 0.05*(maxy-miny) 
        self.yminlim = miny-differ
        self.ymaxlim = maxy+differ
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xminlim = min(self.xax*axMult)
        self.xmaxlim = max(self.xax*axMult)
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim, self.xminlim)
        else:
            a.set_xlim(self.xminlim, self.xmaxlim)
        a.set_ylim(self.yminlim, self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        a = self.fig.gca()
        a.cla()
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.line_xdata = self.xax*axMult
        self.line_ydata = self.data1D
        a.plot(self.xax*axMult, self.data1D)
        if tmpAx is not None:
            a.plot(tmpAx*axMult, tmpdata)
        for i in range(len(tmpAx2)):
            a.plot(tmpAx2[i]*axMult, tmpdata2[i])
        if self.spec == 0:
            if self.current.axType == 0:
                a.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                a.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                a.set_xlabel(r'Time [$\mu$s]')
            else:
                a.set_xlabel('User defined')
        elif self.spec == 1:
            if self.current.ppm:
                a.set_xlabel('Frequency [ppm]')
            else:
                if self.current.axType == 0:
                    a.set_xlabel('Frequency [Hz]')
                elif self.current.axType == 1:
                    a.set_xlabel('Frequency [kHz]')
                elif self.current.axType == 2:
                    a.set_xlabel('Frequency [MHz]')
                else:
                    a.set_xlabel('User defined')
        else:
            a.set_xlabel('')
        a.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        a.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim, self.xminlim)
        else:
            a.set_xlim(self.xminlim, self.xmaxlim)
        a.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()

    def togglePick(self, var):
        self.peakPickReset()
        if var == 1:
            self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos) 
            self.peakPick = True
        else:
            self.peakPickFunc = None
            self.peakPick = False
        
    def pickDeconv(self, pos):
        if self.pickWidth:
            if self.current.spec == 1:
                if self.current.ppm:
                    axMult = 1e6/self.current.ref
                else:
                    axMult = 1.0/(1000.0**self.current.axType)
            elif self.current.spec == 0:
                axMult = 1000.0**self.current.axType
            width = (2*abs(float(self.rootwindow.paramframe.posEntries[self.pickNum].text())-pos[1]))/axMult
            self.rootwindow.paramframe.ampEntries[self.pickNum].setText("%#.3g" %(float(self.rootwindow.paramframe.ampEntries[self.pickNum].text())*width))
            self.rootwindow.paramframe.lorEntries[self.pickNum].setText("%#.3g" % abs(width))
            self.pickNum += 1
            self.pickWidth = False
        else:
            self.rootwindow.paramframe.posEntries[self.pickNum].setText("%#.3g" %pos[1])
            left = pos[0] - 10 
            if left < 0:
                left = 0
            right = pos[0] + 10
            if right >= len(self.data1D) :
                right = len(self.data1D)-1
            self.rootwindow.paramframe.ampEntries[self.pickNum].setText("%#.3g" %(pos[2]*np.pi*0.5))
            if self.pickNum < 10:
                self.rootwindow.paramframe.numExp.setCurrentIndex(self.pickNum)
                self.rootwindow.paramframe.changeNum()
            self.pickWidth = True
        if self.pickNum < 10:
            self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos) 
            self.peakPick = True 

#################################################################################
class PeakDeconvParamFrame(QtGui.QWidget): 
    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        if self.parent.current.spec == 1:
            self.axAdd = self.parent.current.freq-self.parent.current.ref
            if self.parent.current.ppm:
                self.axMult = 1e6/self.parent.current.ref
            else:
                self.axMult = 1.0/(1000.0**self.parent.current.axType)
        elif self.parent.current.spec == 0:
            self.axMult = 1000.0**self.parent.current.axType
            self.axAdd = 0
        self.frame1 = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1, 0, 0)
        grid.addLayout(self.frame2, 0, 1)
        grid.addLayout(self.frame3, 0, 2)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton, 0, 0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton, 1, 0)
        fitAllButton = QtGui.QPushButton("Fit all")
        fitAllButton.clicked.connect(self.fitAll)
        self.frame1.addWidget(fitAllButton, 2, 0)
        copyResultButton = QtGui.QPushButton("Copy result")
        copyResultButton.clicked.connect(lambda: self.sim(True))
        self.frame1.addWidget(copyResultButton, 3, 0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton, 4, 0)
        resetButton = QtGui.QPushButton("Reset")
        resetButton.clicked.connect(self.reset)
        self.frame1.addWidget(resetButton, 0, 1)
        self.pickTick = QtGui.QCheckBox("Pick")
        self.pickTick.stateChanged.connect(self.togglePick)
        self.frame1.addWidget(self.pickTick, 1, 1)
        self.frame1.setColumnStretch(10, 1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel("Bgrnd:"), 0, 0, 1, 2)
        self.bgrndTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.bgrndTick, 1, 0)
        self.bgrndEntry = QtGui.QLineEdit()
        self.bgrndEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.bgrndEntry.setText("0.0")
        self.frame2.addWidget(self.bgrndEntry, 1, 1)
        self.frame2.addWidget(QLabel("Slope:"), 2, 0, 1, 2)
        self.slopeTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.slopeTick, 3, 0)
        self.slopeEntry = QtGui.QLineEdit()
        self.slopeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.slopeEntry.setText("0.0")
        self.frame2.addWidget(self.slopeEntry, 3, 1)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        self.frame3.addWidget(QLabel("Position:"), 1, 0, 1, 2)
        self.frame3.addWidget(QLabel("Integral:"), 1, 2, 1, 2)
        self.frame3.addWidget(QLabel("Lorentz [Hz]:"), 1, 4, 1, 2)
        self.frame3.addWidget(QLabel("Gauss [Hz]:"), 1, 6, 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.posTicks = []
        self.posEntries = []
        self.ampTicks = []
        self.ampEntries = []
        self.lorTicks = []
        self.lorEntries = []
        self.gaussTicks = []
        self.gaussEntries = []
        for i in range(10):
            self.posTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.posTicks[i], i+2, 0)
            self.posEntries.append(QtGui.QLineEdit())
            self.posEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.posEntries[i], i+2, 1)
            self.ampTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.ampTicks[i], i+2, 2)
            self.ampEntries.append(QtGui.QLineEdit())
            self.ampEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.ampEntries[i], i+2, 3)
            self.lorTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.lorTicks[i], i+2, 4)
            self.lorEntries.append(QtGui.QLineEdit())
            self.lorEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.lorEntries[i], i+2, 5)
            self.gaussTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.gaussTicks[i], i+2, 6)
            self.gaussEntries.append(QtGui.QLineEdit())
            self.gaussEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.gaussEntries[i], i+2, 7)
        self.reset()
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)

    def reset(self):
        self.parent.pickNum = 0
        self.bgrndEntry.setText("0.0")
        self.bgrndTick.setChecked(True)
        self.slopeEntry.setText("0.0")
        self.slopeTick.setChecked(True)
        self.numExp.setCurrentIndex(0)
        self.pickTick.setChecked(True)
        for i in range(10):
            self.posEntries[i].setText("0.0")
            self.posTicks[i].setChecked(False)
            self.ampEntries[i].setText("1.0")
            self.ampTicks[i].setChecked(False)
            self.lorEntries[i].setText("1.0")
            self.lorTicks[i].setChecked(False)
            self.gaussEntries[i].setText("0.0")
            self.gaussTicks[i].setChecked(True)
            if i > 0:
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()
        self.togglePick()
        self.parent.pickWidth = False
        self.parent.showPlot()

    def changeNum(self, *args):
        val = self.numExp.currentIndex()+1
        for i in range(10):
            if i < val:
                self.posTicks[i].show()
                self.posEntries[i].show()
                self.ampTicks[i].show()
                self.ampEntries[i].show()
                self.lorTicks[i].show()
                self.lorEntries[i].show()
                self.gaussTicks[i].show()
                self.gaussEntries[i].show()
            else:
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()

    def togglePick(self):
        self.parent.togglePick(self.pickTick.isChecked())
                
    def fitFunc(self, x, *param):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        testFunc = np.zeros(len(self.parent.data1D))
        if struc[0]:
            bgrnd = param[0]
            param = np.delete(param, [0])
        else:
            bgrnd = argu[0]
            argu = np.delete(argu, [0])
        if struc[1]:
            slope = param[0]
            param = np.delete(param, [0])
        else:
            slope = argu[0]
            argu = np.delete(argu, [0])
        for i in range(numExp):
            if struc[4*i+2]:
                pos = param[0]
                param = np.delete(param, [0])
            else:
                pos = argu[0]
                argu = np.delete(argu, [0])
            if struc[4*i+3]:
                amp = param[0]
                param = np.delete(param, [0])
            else:
                amp = argu[0]
                argu = np.delete(argu, [0])
            if struc[4*i+4]:
                width = abs(param[0])
                param = np.delete(param, [0])
            else:
                width = argu[0]
                argu = np.delete(argu, [0])
            if struc[4*i+5]:
                gauss = abs(param[0])
                param = np.delete(param, [0])
            else:
                gauss = argu[0]
                argu = np.delete(argu, [0])
            t = np.arange(len(x))/self.parent.current.sw
            timeSignal = np.exp(1j*2*np.pi*t*(pos/self.axMult-self.axAdd))*np.exp(-np.pi*width*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))*2/self.parent.current.sw
            timeSignal[0] = timeSignal[0]*0.5
            testFunc += amp*np.real(np.fft.fftshift(np.fft.fft(timeSignal)))
        testFunc += bgrnd+slope*x
        return testFunc

    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.bgrndEntry.text())
        if inp is None:
            return False
        self.bgrndEntry.setText('%#.3g' % inp)
        inp = safeEval(self.slopeEntry.text())
        if inp is None:
            return False
        self.slopeEntry.setText('%#.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.posEntries[i].text())
            if inp is None:
                return False
            self.posEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.ampEntries[i].text())
            if inp is None:
                return False
            self.ampEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.lorEntries[i].text())
            if inp is None:
                return False
            self.lorEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.gaussEntries[i].text())
            if inp is None:
                return False
            self.gaussEntries[i].setText('%#.3g' % inp)
        return True
    
    def fit(self, *args):
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = float(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            if not self.posTicks[i].isChecked():
                guess.append(float(self.posEntries[i].text()))
                struc.append(True)
            else:
                outPos[i] = float(self.posEntries[i].text())
                argu.append(outPos[i])
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].isChecked():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp, struc, argu)
        fitVal = scipy.optimize.curve_fit(self.fitFunc, self.parent.xax, self.parent.data1D, p0=guess)
        counter = 0
        if struc[0]:
            self.bgrndEntry.setText('%#.3g' % fitVal[0][counter])
            outBgrnd = fitVal[0][counter]
            counter += 1
        if struc[1]:
            self.slopeEntry.setText('%#.3g' % fitVal[0][counter])
            outSlope = fitVal[0][counter]
            counter += 1
        for i in range(numExp):
            if struc[4*i+2]:
                self.posEntries[i].setText('%#.3g' % fitVal[0][counter])
                outPos[i] = fitVal[0][counter]
                counter += 1
            if struc[4*i+3]:
                self.ampEntries[i].setText('%#.3g' % fitVal[0][counter])
                outAmp[i] = fitVal[0][counter]
                counter += 1
            if struc[4*i+4]:
                self.lorEntries[i].setText('%#.3g' % abs(fitVal[0][counter]))
                outWidth[i] = abs(fitVal[0][counter])
                counter += 1
            if struc[4*i+5]:
                self.gaussEntries[i].setText('%#.3g' % abs(fitVal[0][counter]))
                outGauss[i] = abs(fitVal[0][counter])
                counter += 1
        self.disp(outBgrnd, outSlope, outAmp, outPos, outWidth, outGauss)

    def fitAll(self, *args):
        FitAllSelectionWindow(self, ["Background", "Slope", "Position", "Integral", "Lorentz", "Gauss"])
        
    def fitAllFunc(self, outputs):
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = float(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            if not self.posTicks[i].isChecked():
                guess.append(float(self.posEntries[i].text()))
                struc.append(True)
            else:
                outPos[i] = float(self.posEntries[i].text())
                argu.append(outPos[i])
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].isChecked():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp, struc, argu)
        fullData = self.parent.current.data.data
        axes = self.parent.current.axes
        dataShape = fullData.shape
        dataShape2 = np.delete(dataShape, axes)
        rolledData = np.rollaxis(fullData, axes)
        intOutputs = np.array(outputs, dtype=int)
        numOutputs = np.sum(intOutputs[:2]) + numExp*np.sum(intOutputs[2:])
        outputData = np.zeros((np.product(dataShape2), numOutputs), dtype=complex)
        counter2 = 0
        for j in rolledData.reshape(dataShape[axes], np.product(dataShape2)).T:
            try:
                fitVal = scipy.optimize.curve_fit(self.fitFunc, self.parent.xax, np.real(j), guess)
            except:
                fitVal = [[0]*10]
        
            counter = 0
            if struc[0]:
                outBgrnd = fitVal[0][counter]
                counter += 1
            if struc[1]:
                outSlope = fitVal[0][counter]
                counter += 1
            for i in range(numExp):
                if struc[4*i+2]:
                    outPos[i] = fitVal[0][counter]
                    counter += 1
                if struc[4*i+3]:
                    outAmp[i] = fitVal[0][counter]
                    counter += 1
                if struc[4*i+4]:
                    outWidth[i] = abs(fitVal[0][counter])
                    counter += 1
                if struc[4*i+5]:
                    outGauss[i] = abs(fitVal[0][counter])
                    counter += 1
            outputArray = []
            if outputs[0]:
                outputArray = np.concatenate((outputArray, [outBgrnd]))
            if outputs[1]:
                outputArray = np.concatenate((outputArray, [outSlope]))
            if outputs[2]:
                outputArray = np.concatenate((outputArray, outPos))
            if outputs[3]:
                outputArray = np.concatenate((outputArray, outAmp))
            if outputs[4]:
                outputArray = np.concatenate((outputArray, outWidth))
            if outputs[5]:
                outputArray = np.concatenate((outputArray, outGauss))
            outputData[counter2] = outputArray
            counter2 += 1
        newShape = np.concatenate((np.array(dataShape2), [numOutputs]))
        self.rootwindow.createNewData(np.rollaxis(outputData.reshape(newShape), -1, axes), axes)

    def sim(self, store=False):
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        outBgrnd = safeEval(self.bgrndEntry.text())
        outSlope = safeEval(self.slopeEntry.text())
        if outBgrnd is None or outSlope is None:
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        for i in range(numExp):
            outPos[i] = safeEval(self.posEntries[i].text())
            outAmp[i] = safeEval(self.ampEntries[i].text())
            outWidth[i] = abs(safeEval(self.lorEntries[i].text()))
            outGauss[i] = abs(safeEval(self.gaussEntries[i].text()))
            if not np.isfinite([outPos[i], outAmp[i], outWidth[i], outGauss[i]]).all():
                self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
                return
            self.lorEntries[i].setText('%#.3g' % outWidth[i])
            self.gaussEntries[i].setText('%#.3g' % outGauss[i])
        self.disp(outBgrnd, outSlope, outAmp, outPos, outWidth, outGauss, store)
        
    def disp(self, outBgrnd, outSlope, outAmp, outPos, outWidth, outGauss, store=False):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x = []
        t = np.arange(len(tmpx))/self.parent.current.sw
        for i in range(len(outAmp)):
            x.append(tmpx)
            timeSignal = np.exp(1j*2*np.pi*t*(outPos[i]/self.axMult-self.axAdd))*np.exp(-np.pi*outWidth[i]*t)*np.exp(-((np.pi*outGauss[i]*t)**2)/(4*np.log(2)))*2/self.parent.current.sw
            timeSignal[0] = timeSignal[0]*0.5
            y = outAmp[i]*np.real(np.fft.fftshift(np.fft.fft(timeSignal)))
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        if store:
            outCurvePart.append(outCurve)
            outCurvePart.append(self.parent.data1D)
            self.rootwindow.createNewData(np.array(outCurvePart), self.parent.current.axes, True)
        else:
            self.parent.showPlot(tmpx, outCurve, x, outCurvePart)

##############################################################################
class TensorDeconvWindow(FittingWindow): 
    def __init__(self, mainProgram, oldMainWindow):
        FittingWindow.__init__(self, mainProgram, oldMainWindow)

    def setup(self):
        self.current = TensorDeconvFrame(self, self.fig, self.canvas, self.oldMainWindow.current)
        self.paramframe = TensorDeconvParamFrame(self.current, self)
        
#################################################################################   
class TensorDeconvFrame(Plot1DFrame): 
    def __init__(self, rootwindow, fig, canvas, current):
        Plot1DFrame.__init__(self, rootwindow, fig, canvas)
        self.data1D = current.getDisplayedData()
        self.current = current
        self.spec = self.current.spec
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xax = self.current.xax*axMult
        self.plotType = 0
        self.rootwindow = rootwindow
        self.pickNum = 0
        self.pickNum2 = 0
        self.plotReset()
        self.showPlot()

    def plotReset(self): #set the plot limits to min and max values
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
        self.yminlim = miny-differ
        self.ymaxlim = maxy+differ
        self.xminlim = min(self.xax)
        self.xmaxlim = max(self.xax)
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        self.ax.cla()
        self.line_xdata = self.xax
        self.line_ydata = self.data1D
        self.ax.plot(self.xax, self.data1D)
        self.ax.plot(tmpAx, tmpdata)
        for i in range(len(tmpAx2)):
            self.ax.plot(tmpAx2[i], tmpdata2[i])
        if self.spec == 0:
            if self.current.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.current.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.current.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.current.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.current.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()

    def togglePick(self, var):
        self.peakPickReset()
        if var == 1:
            self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos) 
            self.peakPick = True
        else:
            self.peakPickFunc = None
            self.peakPick = False
        
    def pickDeconv(self, pos):
        if self.pickNum2 == 0:
            if self.pickNum < 10:
                self.rootwindow.paramframe.numExp.setCurrentIndex(self.pickNum)
                self.rootwindow.paramframe.changeNum()
            self.rootwindow.paramframe.t11Entries[self.pickNum].setText("%.3g" %self.current.xax[pos[0]])
            self.pickNum2 = 1
        elif self.pickNum2 == 1:
            self.rootwindow.paramframe.t22Entries[self.pickNum].setText("%.3g" %self.current.xax[pos[0]])
            self.pickNum2 = 2
        elif self.pickNum2 == 2:
            self.rootwindow.paramframe.t33Entries[self.pickNum].setText("%.3g" %self.current.xax[pos[0]])
            self.pickNum2 = 0
            self.pickNum += 1
        if self.pickNum < 10:
            self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos) 
            self.peakPick = True 

#################################################################################
class TensorDeconvParamFrame(QtGui.QWidget): 
    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        self.cheng = 15
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        if self.parent.current.spec == 1:
            self.axAdd = self.parent.current.freq-self.parent.current.ref
        elif self.parent.current.spec == 0:
            self.axAdd = 0
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1, 0, 0)
        grid.addLayout(self.optframe, 0, 1)
        grid.addLayout(self.frame2, 0, 2)
        grid.addLayout(self.frame3, 0, 3)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton, 0, 0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton, 1, 0)
        fitAllButton = QtGui.QPushButton("Fit all")
        fitAllButton.clicked.connect(self.fitAll)
        self.frame1.addWidget(fitAllButton, 2, 0)
        copyResultButton = QtGui.QPushButton("Copy result")
        copyResultButton.clicked.connect(lambda: self.sim(True))
        self.frame1.addWidget(copyResultButton, 3, 0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton, 4, 0)
        resetButton = QtGui.QPushButton("Reset")
        resetButton.clicked.connect(self.reset)
        self.frame1.addWidget(resetButton, 0, 1)
        self.pickTick = QtGui.QCheckBox("Pick")
        self.pickTick.stateChanged.connect(self.togglePick)
        self.frame1.addWidget(self.pickTick, 1, 1)
        self.frame1.setColumnStretch(10, 1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.optframe.addWidget(QLabel("Cheng:"), 0, 0)
        self.chengEntry = QtGui.QLineEdit()
        self.chengEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.chengEntry.setText(str(self.cheng))
        self.optframe.addWidget(self.chengEntry, 1, 0)
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel("Bgrnd:"), 0, 0, 1, 2)
        self.bgrndTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.bgrndTick, 1, 0)
        self.bgrndEntry = QtGui.QLineEdit()
        self.bgrndEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.bgrndEntry.setText("0.0")
        self.frame2.addWidget(self.bgrndEntry, 1, 1)
        self.frame2.addWidget(QLabel("Slope:"), 2, 0, 1, 2)
        self.slopeTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.slopeTick, 3, 0)
        self.slopeEntry = QtGui.QLineEdit()
        self.slopeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.slopeEntry.setText("0.0")
        self.frame2.addWidget(self.slopeEntry, 3, 1)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        self.frame3.addWidget(QLabel("T11:"), 1, 0, 1, 2)
        self.frame3.addWidget(QLabel("T22:"), 1, 2, 1, 2)
        self.frame3.addWidget(QLabel("T33:"), 1, 4, 1, 2)
        self.frame3.addWidget(QLabel("Integral:"), 1, 6, 1, 2)
        self.frame3.addWidget(QLabel("Lorentz [Hz]:"), 1, 8, 1, 2)
        self.frame3.addWidget(QLabel("Gauss [Hz]:"), 1, 10, 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.t11Entries = []
        self.t11Ticks = []
        self.t22Entries = []
        self.t22Ticks = []
        self.t33Entries = []
        self.t33Ticks = []
        self.ampEntries = []
        self.ampTicks = []
        self.lorEntries = []
        self.lorTicks = []
        self.gaussEntries = []
        self.gaussTicks = []
        for i in range(10):
            self.t11Ticks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.t11Ticks[i], i+2, 0)
            self.t11Entries.append(QtGui.QLineEdit())
            self.t11Entries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.t11Entries[i], i+2, 1)
            self.t22Ticks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.t22Ticks[i], i+2, 2)
            self.t22Entries.append(QtGui.QLineEdit())
            self.t22Entries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.t22Entries[i], i+2, 3)
            self.t33Ticks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.t33Ticks[i], i+2, 4)
            self.t33Entries.append(QtGui.QLineEdit())
            self.t33Entries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.t33Entries[i], i+2, 5)
            self.ampTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.ampTicks[i], i+2, 6)
            self.ampEntries.append(QtGui.QLineEdit())
            self.ampEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.ampEntries[i], i+2, 7)
            self.lorTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.lorTicks[i], i+2, 8)
            self.lorEntries.append(QtGui.QLineEdit())
            self.lorEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.lorEntries[i], i+2, 9)
            self.gaussTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.gaussTicks[i], i+2, 10)
            self.gaussEntries.append(QtGui.QLineEdit())
            self.gaussEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.gaussEntries[i], i+2, 11)
        self.reset()
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        
    def reset(self):
        self.parent.pickNum = 0
        self.parent.pickNum2 = 0
        self.bgrndEntry.setText("0.0")
        self.bgrndTick.setChecked(True)
        self.slopeEntry.setText("0.0")
        self.slopeTick.setChecked(True)
        self.cheng = 15
        self.chengEntry.setText(str(self.cheng))
        self.pickTick.setChecked(True)
        for i in range(10):
            self.t11Entries[i].setText("0.0")
            self.t11Ticks[i].setChecked(False)
            self.t22Entries[i].setText("0.0")
            self.t22Ticks[i].setChecked(False)
            self.t33Entries[i].setText("0.0")
            self.t33Ticks[i].setChecked(False)
            self.ampEntries[i].setText("1.0")
            self.ampTicks[i].setChecked(False)
            self.lorEntries[i].setText("10.0")
            self.lorTicks[i].setChecked(True)
            self.gaussEntries[i].setText("0.0")
            self.gaussTicks[i].setChecked(True)
            if i > 0:
                self.t11Ticks[i].hide()
                self.t11Entries[i].hide()
                self.t22Ticks[i].hide()
                self.t22Entries[i].hide()
                self.t33Ticks[i].hide()
                self.t33Entries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()
        self.togglePick()
        self.parent.showPlot()

    def togglePick(self):
        self.parent.togglePick(self.pickTick.isChecked())
                
    def setCheng(self, *args):
        inp = safeEval(self.chengEntry.text())
        if inp is None:
            self.cheng = 15
        else:
            self.cheng = int(inp)
        self.chengEntry.setText(str(self.cheng))
                
    def changeNum(self, *args):
        val = self.numExp.currentIndex()+1
        for i in range(10):
            if i < val:
                self.t11Ticks[i].show()
                self.t11Entries[i].show()
                self.t22Ticks[i].show()
                self.t22Entries[i].show()
                self.t33Ticks[i].show()
                self.t33Entries[i].show()
                self.ampTicks[i].show()
                self.ampEntries[i].show()
                self.lorTicks[i].show()
                self.lorEntries[i].show()
                self.gaussTicks[i].show()
                self.gaussEntries[i].show()
            else:
                self.t11Ticks[i].hide()
                self.t11Entries[i].hide()
                self.t22Ticks[i].hide()
                self.t22Entries[i].hide()
                self.t33Ticks[i].hide()
                self.t33Entries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()

    def tensorFunc(self, x, t11, t22, t33, lor, gauss):
        t11 = t11*self.multt11
        t22 = t22*self.multt22
        t33 = t33*self.multt33
        v = t11+t22+t33-self.axAdd
        length = len(x)
        t = np.arange(length)/self.parent.current.sw
        final = np.zeros(length)
        mult = v/(self.parent.current.sw)*length
        x1 = np.array(np.round(mult)+np.floor(length/2.0), dtype=int)
        weight = self.weight[np.logical_and(x1>=0, x1<length)]
        x1 = x1[np.logical_and(x1>=0, x1<length)]
        final = np.bincount(x1, weight, length)
        apod = np.exp(-np.pi*lor*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
        apod[-1:-(len(apod)/2+1):-1] = apod[:len(apod)/2]
        I = np.real(np.fft.fft(np.fft.ifft(final)*apod))
        I = I/self.parent.current.sw*len(I)
        return I
                
    def fitFunc(self, param, x, y):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        testFunc = np.zeros(len(self.parent.data1D))
        if struc[0]:
            bgrnd = param[0]
            param = np.delete(param, [0])
        else:
            bgrnd = argu[0]
            argu = np.delete(argu, [0])
        if struc[1]:
            slope = param[0]
            param = np.delete(param, [0])
        else:
            slope = argu[0]
            argu = np.delete(argu, [0])
        for i in range(numExp):
            if struc[6*i+2]:
                t11 = param[0]
                param = np.delete(param, [0])
            else:
                t11 = argu[0]
                argu = np.delete(argu, [0])
            if struc[6*i+3]:
                t22 = param[0]
                param = np.delete(param, [0])
            else:
                t22 = argu[0]
                argu = np.delete(argu, [0])
            if struc[6*i+4]:
                t33 = param[0]
                param = np.delete(param, [0])
            else:
                t33 = argu[0]
                argu = np.delete(argu, [0])
            if struc[6*i+5]:
                amp = param[0]
                param = np.delete(param, [0])
            else:
                amp = argu[0]
                argu = np.delete(argu, [0])
            if struc[6*i+6]:
                width = abs(param[0])
                param = np.delete(param, [0])
            else:
                width = argu[0]
                argu = np.delete(argu, [0])
            if struc[6*i+7]:
                gauss = abs(param[0])
                param = np.delete(param, [0])
            else:
                gauss = argu[0]
                argu = np.delete(argu, [0])
            testFunc += amp*self.tensorFunc(x, t11, t22, t33, width, gauss)
        testFunc += bgrnd+slope*x
        return np.sum((np.real(testFunc)-y)**2)

    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.bgrndEntry.text())
        if inp is None:
            return False
        self.bgrndEntry.setText('%#.3g' % inp)
        inp = safeEval(self.slopeEntry.text())
        if inp is None:
            return False
        self.slopeEntry.setText('%#.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.t11Entries[i].text())
            if inp is None:
                return False
            self.t11Entries[i].setText('%#.3g' % inp)
            inp = safeEval(self.t22Entries[i].text())
            if inp is None:
                return False
            self.t22Entries[i].setText('%#.3g' % inp)
            inp = safeEval(self.t33Entries[i].text())
            if inp is None:
                return False
            self.t33Entries[i].setText('%#.3g' % inp)
            inp = safeEval(self.ampEntries[i].text())
            if inp is None:
                return False
            self.ampEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.lorEntries[i].text())
            if inp is None:
                return False
            self.lorEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.gaussEntries[i].text())
            if inp is None:
                return False
            self.gaussEntries[i].setText('%#.3g' % inp)
        return True
        
    def fit(self, *args):
        self.setCheng()
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outt11 = np.zeros(numExp)
        outt22 = np.zeros(numExp)
        outt33 = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = float(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            if not self.t11Ticks[i].isChecked():
                guess.append(float(self.t11Entries[i].text()))
                struc.append(True)
            else:
                outt11[i] = float(self.t11Entries[i].text())
                argu.append(outt11[i])
                struc.append(False)
            if not self.t22Ticks[i].isChecked():
                guess.append(float(self.t22Entries[i].text()))
                struc.append(True)
            else:
                outt22[i] = float(self.t22Entries[i].text())
                argu.append(outt22[i])
                struc.append(False)
            if not self.t33Ticks[i].isChecked():
                guess.append(float(self.t33Entries[i].text()))
                struc.append(True)
            else:
                outt33[i] = float(self.t33Entries[i].text())
                argu.append(outt33[i])
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].isChecked():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp, struc, argu)
        phi, theta, self.weight = zcw_angles(self.cheng, symm=2)
        self.multt11 = np.sin(theta)**2*np.cos(phi)**2
        self.multt22 = np.sin(theta)**2*np.sin(phi)**2
        self.multt33 = np.cos(theta)**2
        fitVal = scipy.optimize.fmin(self.fitFunc, guess, args=(self.parent.xax, np.real(self.parent.data1D)))
        counter = 0
        if struc[0]:
            self.bgrndEntry.setText('%.3g' % fitVal[counter])
            outBgrnd = fitVal[counter]
            counter += 1
        if struc[1]:
            self.slopeEntry.setText('%.3g' % fitVal[counter])
            outSlope = fitVal[counter]
            counter += 1
        for i in range(numExp):
            if struc[6*i+2]:
                self.t11Entries[i].setText('%.3g' % fitVal[counter])
                outt11[i] = fitVal[counter]
                counter += 1
            if struc[6*i+3]:
                self.t22Entries[i].setText('%.3g' % fitVal[counter])
                outt22[i] = fitVal[counter]
                counter += 1
            if struc[6*i+4]:
                self.t33Entries[i].setText('%.3g' % fitVal[counter])
                outt33[i] = fitVal[counter]
                counter += 1
            if struc[6*i+5]:
                self.ampEntries[i].setText('%.3g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[6*i+6]:
                self.lorEntries[i].setText('%.3g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
            if struc[6*i+7]:
                self.gaussEntries[i].setText('%.3g' % abs(fitVal[counter]))
                outGauss[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd, outSlope, outt11, outt22, outt33, outAmp, outWidth, outGauss)

    def fitAll(self, *args):
        FitAllSelectionWindow(self, ["Background", "Slope", "T11", "T22", "T33", "Integral", "Lorentz", "Gauss"])

    def fitAllFunc(self, outputs):
        self.setCheng()
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outt11 = np.zeros(numExp)
        outt22 = np.zeros(numExp)
        outt33 = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = float(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            if not self.t11Ticks[i].isChecked():
                guess.append(float(self.t11Entries[i].text()))
                struc.append(True)
            else:
                outt11[i] = float(self.t11Entries[i].text())
                argu.append(outt11[i])
                struc.append(False)
            if not self.t22Ticks[i].isChecked():
                guess.append(float(self.t22Entries[i].text()))
                struc.append(True)
            else:
                outt22[i] = float(self.t22Entries[i].text())
                argu.append(outt22[i])
                struc.append(False)
            if not self.t33Ticks[i].isChecked():
                guess.append(float(self.t33Entries[i].text()))
                struc.append(True)
            else:
                outt33[i] = float(self.t33Entries[i].text())
                argu.append(outt33[i])
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].isChecked():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp, struc, argu)
        phi, theta, self.weight = zcw_angles(self.cheng, symm=2)
        self.multt11 = np.sin(theta)**2*np.cos(phi)**2
        self.multt22 = np.sin(theta)**2*np.sin(phi)**2
        self.multt33 = np.cos(theta)**2
        fullData = self.parent.current.data.data
        axes = self.parent.current.axes
        dataShape = fullData.shape
        dataShape2 = np.delete(dataShape, axes)
        rolledData = np.rollaxis(fullData, axes)
        intOutputs = np.array(outputs, dtype=int)
        numOutputs = np.sum(intOutputs[:2]) + numExp*np.sum(intOutputs[2:])
        outputData = np.zeros((np.product(dataShape2), numOutputs), dtype=complex)
        counter2 = 0
        for j in rolledData.reshape(dataShape[axes], np.product(dataShape2)).T:
            try:
                fitVal = scipy.optimize.fmin(self.fitFunc, guess, args=(self.parent.xax, np.real(j)))
            except:
                fitVal = [[0]*10]
            counter = 0
            if struc[0]:
                outBgrnd = fitVal[counter]
                counter += 1
            if struc[1]:
                outSlope = fitVal[counter]
                counter += 1
            for i in range(numExp):
                if struc[6*i+2]:
                    outt11[i] = fitVal[counter]
                    counter += 1
                if struc[6*i+3]:
                    outt22[i] = fitVal[counter]
                    counter += 1
                if struc[6*i+4]:
                    outt33[i] = fitVal[counter]
                    counter += 1
                if struc[6*i+5]:
                    outAmp[i] = fitVal[counter]
                    counter += 1
                if struc[6*i+6]:
                    outWidth[i] = abs(fitVal[counter])
                    counter += 1
                if struc[6*i+7]:
                    outGauss[i] = abs(fitVal[counter])
                    counter += 1
            outputArray = []
            if outputs[0]:
                outputArray = np.concatenate((outputArray, [outBgrnd]))
            if outputs[1]:
                outputArray = np.concatenate((outputArray, [outSlope]))
            if outputs[2]:
                outputArray = np.concatenate((outputArray, outt11))
            if outputs[3]:
                outputArray = np.concatenate((outputArray, outt22))
            if outputs[4]:
                outputArray = np.concatenate((outputArray, outt33))
            if outputs[5]:
                outputArray = np.concatenate((outputArray, outAmp)) 
            if outputs[6]:
                outputArray = np.concatenate((outputArray, outWidth))
            if outputs[7]:
                outputArray = np.concatenate((outputArray, outGauss))
            outputData[counter2] = outputArray
            counter2 += 1
        newShape = np.concatenate((np.array(dataShape2), [numOutputs]))
        self.rootwindow.createNewData(np.rollaxis(outputData.reshape(newShape), -1, axes), axes)
        
    def sim(self, store=False):
        self.setCheng()
        numExp = self.numExp.currentIndex()+1
        bgrnd = safeEval(self.bgrndEntry.text())
        slope = safeEval(self.slopeEntry.text())
        if bgrnd is None or slope is None:
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        t11 = np.zeros(numExp)
        t22 = np.zeros(numExp)
        t33 = np.zeros(numExp)
        amp = np.zeros(numExp)
        width = np.zeros(numExp)
        gauss = np.zeros(numExp)
        for i in range(numExp):
            t11[i] = safeEval(self.t11Entries[i].text())
            t22[i] = safeEval(self.t22Entries[i].text())
            t33[i] = safeEval(self.t33Entries[i].text())
            amp[i] = safeEval(self.ampEntries[i].text())
            width[i] = safeEval(self.lorEntries[i].text())
            gauss[i] = safeEval(self.gaussEntries[i].text())
            if not np.isfinite([t11[i], t22[i], t33[i], amp[i], width[i], gauss[i]]).all():
                self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
                return
        phi, theta, self.weight = zcw_angles(self.cheng, symm=2)
        self.multt11 = np.sin(theta)**2*np.cos(phi)**2
        self.multt22 = np.sin(theta)**2*np.sin(phi)**2
        self.multt33 = np.cos(theta)**2
        self.disp(bgrnd, slope, t11, t22, t33, amp, width, gauss, store)
        
    def disp(self, outBgrnd, outSlope, outt11, outt22, outt33, outAmp, outWidth, outGauss, store=False):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x = []
        for i in range(len(outt11)):
            x.append(tmpx)
            y =  outAmp[i]*self.tensorFunc(tmpx, outt11[i], outt22[i], outt33[i], outWidth[i], outGauss[i])
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        if store:
            outCurvePart.append(outCurve)
            outCurvePart.append(self.parent.data1D)
            self.rootwindow.createNewData(np.array(outCurvePart), self.parent.current.axes, True)
        else:
            self.parent.showPlot(tmpx, outCurve, x, outCurvePart)
        
##############################################################################
class HerzfeldBergerWindow(FittingWindow): 
    def __init__(self, mainProgram, oldMainWindow):
        FittingWindow.__init__(self, mainProgram, oldMainWindow)

    def setup(self):
        self.current = HerzfeldBergerFrame(self, self.fig, self.canvas, self.oldMainWindow.current)
        self.paramframe = HerzfeldBergerParamFrame(self.current, self)
        
#################################################################################   
class HerzfeldBergerFrame(Plot1DFrame): 
    def __init__(self, rootwindow, fig, canvas, current):
        Plot1DFrame.__init__(self, rootwindow, fig, canvas)
        self.data1D = current.getDisplayedData()
        self.current = current
        self.spec = self.current.spec
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xax = self.current.xax*axMult
        self.plotType = 0
        self.rootwindow = rootwindow
        self.plotReset()
        self.showPlot()

    def plotReset(self): #set the plot limits to min and max values
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
        self.yminlim = miny-differ
        self.ymaxlim = maxy+differ
        self.xminlim = min(self.xax)
        self.xmaxlim = max(self.xax)
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        self.ax.cla()
        self.line_xdata = self.xax
        self.line_ydata = self.data1D
        self.ax.plot(self.xax, self.data1D)
        self.ax.plot(tmpAx, tmpdata)
        for i in range(len(tmpAx2)):
            self.ax.plot(tmpAx2[i], tmpdata2[i])
        if self.spec == 0:
            if self.current.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.current.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.current.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.current.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.current.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()

    def togglePick(self, var):
        self.peakPickReset()
        if var == 1:
            self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos) 
            self.peakPick = True
        else:
            self.peakPickFunc = None
            self.peakPick = False
        
    def pickDeconv(self, pos):
        self.rootwindow.paramframe.addValue(pos[0])
        self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos) 
        self.peakPick = True 

#################################################################################
class HerzfeldBergerParamFrame(QtGui.QWidget): 
    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        self.cheng = 15
        self.NSTEPS = 30
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        if self.parent.current.spec == 1:
            self.axAdd = self.parent.current.freq-self.parent.current.ref
        elif self.parent.current.spec == 0:
            self.axAdd = 0
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        self.frame4 = QtGui.QGridLayout()
        grid.addLayout(self.frame1, 0, 0, 2, 1)
        grid.addLayout(self.optframe, 0, 1, 2, 1)
        grid.addLayout(self.frame2, 0, 2)
        grid.addLayout(self.frame3, 0, 3)
        grid.addLayout(self.frame4, 1, 3)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton, 0, 0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton, 1, 0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton, 3, 0)
        resetButton = QtGui.QPushButton("Reset")
        resetButton.clicked.connect(self.reset)
        self.frame1.addWidget(resetButton, 0, 1)
        self.pickTick = QtGui.QCheckBox("Pick")
        self.pickTick.stateChanged.connect(self.togglePick)
        self.frame1.addWidget(self.pickTick, 1, 1)
        self.frame1.setColumnStretch(10, 1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.optframe.addWidget(QLabel("Cheng:"), 0, 0)
        self.chengEntry = QtGui.QLineEdit()
        self.chengEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.chengEntry.setText(str(self.cheng))
        self.optframe.addWidget(self.chengEntry, 1, 0)
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.optframe.addWidget(QLabel("Spinning speed [kHz]:"), 2, 0)
        self.spinEntry = QtGui.QLineEdit()
        self.spinEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.spinEntry.setText("30.0")
        self.optframe.addWidget(self.spinEntry, 3, 0)
        #self.frame2.setColumnStretch(10, 1)
        #self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.frame3.addWidget(QLabel(u"\u03B4 [ppm]:"), 1, 0, 1, 2)
        self.frame3.addWidget(QLabel(u"\u03B7:"), 1, 2, 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.deltaEntry = QtGui.QLineEdit()
        self.deltaEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.deltaEntry.setText("10.0")
        self.frame3.addWidget(self.deltaEntry, 2, 1)
        self.deltaTick = QtGui.QCheckBox('')
        self.frame3.addWidget(self.deltaTick, 2, 0)
        self.etaEntry = QtGui.QLineEdit()
        self.etaEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.etaEntry.setText("0.00")
        self.frame3.addWidget(self.etaEntry, 2, 3)
        self.etaTick = QtGui.QCheckBox('')
        self.frame3.addWidget(self.etaTick, 2, 2)
        self.frame4.addWidget(QLabel("Sideband:"), 0, 0)
        self.frame4.addWidget(QLabel("Integral:"), 1, 0)
        self.frame4.addWidget(QLabel("Result:"), 2, 0)
        self.frame4.setColumnStretch(100, 1)
        self.sidebandList = []
        self.sidebandEntries = []
        self.integralList = []
        self.integralEntries = []
        self.resultList = []
        self.resultLabels = []
        self.tmpPos = None
        self.reset()
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)

    def reset(self):
        self.sidebandList = []
        self.integralList = []
        self.resultList = []
        for i in range(len(self.sidebandEntries)):
            self.frame4.removeWidget(self.sidebandEntries[i])
            self.frame4.removeWidget(self.integralEntries[i])
            self.frame4.removeWidget(self.resultLabels[i])
            self.sidebandEntries[i].deleteLater()
            self.integralEntries[i].deleteLater()
            self.resultLabels[i].deleteLater()
        self.sidebandEntries = []
        self.integralEntries = []
        self.resultLabels = []
        self.pickTick.setChecked(True)
        self.togglePick()
        self.parent.showPlot()

    def addValue(self, value):
        if self.tmpPos is None:
            self.tmpPos = value
        else:
            n = len(self.integralList)
            self.sidebandList.append((-1)**n*n//2)
            self.integralList.append(np.sum(self.parent.data1D[min(self.tmpPos, value):max(self.tmpPos, value)])*self.parent.current.sw/float(self.parent.data1D.shape[-1]))
            self.tmpPos = None
            self.sidebandEntries.append(QtGui.QSpinBox(self, ))
            self.sidebandEntries[-1].setMinimum(-99)
            self.sidebandEntries[-1].setValue(self.sidebandList[-1])
            self.frame4.addWidget(self.sidebandEntries[-1], 0, len(self.sidebandEntries))
            self.integralEntries.append(QtGui.QLineEdit())
            self.integralEntries[-1].setText('%#.5g' % self.integralList[-1])
            self.frame4.addWidget(self.integralEntries[-1], 1, len(self.sidebandEntries))
            self.resultLabels.append(QtGui.QLabel())
            self.frame4.addWidget(self.resultLabels[-1], 2, len(self.sidebandEntries))

    def togglePick(self):
        self.parent.togglePick(self.pickTick.isChecked())
                
    def setCheng(self, *args):
        inp = safeEval(self.chengEntry.text())
        if inp is None:
            self.cheng = 15
        else:
            self.cheng = int(inp)
        self.chengEntry.setText(str(self.cheng))

    def hbFunc(self, omega0, delta, eta):
        omegars =  omega0*delta*(self.C1  + self.C2 +  eta*(self.C1eta + self.C2eta + self.S1 + self.S2 ))
        #QTrs = np.exp(-1j*np.cumsum(omegars, axis=1)*self.tresolution)
        nsteps = self.C1.shape[1]
        QTrs = np.concatenate([np.ones([self.C1.shape[0],1]),np.exp(-1j*np.cumsum(omegars, axis=1)*self.tresolution)[:,:-1]],1)
        for j in range(1,nsteps):
            QTrs[:,j] = np.exp(-1j*np.sum(omegars[:,0:j]*self.tresolution,1))
        rhoT0sr = np.conj(QTrs)
        #calculate the gamma-averaged FID over 1 rotor period for all crystallites
        favrs = np.zeros(self.NSTEPS, dtype=complex)
        for j in range(self.NSTEPS):
            favrs[j] += np.sum(self.weight * np.sum(rhoT0sr * np.roll(QTrs, -j, axis=1), 1) / self.NSTEPS**2)
        #calculate the sideband intensities by doing an FT and pick the ones that are needed further
        sidebands = np.real(np.fft.fft(favrs))
        return sidebands
                
    def fitFunc(self, param, x, y):
        struc = self.args[0]
        argu = self.args[1]
        omega0 = self.args[2]
        if struc[0]:
            delta = param[0]
            param = np.delete(param, [0])
        else:
            delta = argu[0]
            argu = np.delete(argu, [0])
        if struc[1]:
            eta = param[0]
            param = np.delete(param, [0])
        else:
            eta = argu[0]
            argu = np.delete(argu, [0])
        testFunc = self.hbFunc(omega0, delta, eta)
        testFunc = testFunc[x]/np.sum(testFunc[x])*np.sum(self.integralList)
        return np.sum((testFunc-y)**2)

    def disp(self, outDelta, outEta):
        testFunc = self.hbFunc(self.parent.current.freq*np.pi*2, outDelta, outEta)
        results = testFunc[self.sidebandList]
        results /= np.sum(results)
        for i in range(len(self.resultLabels)):
            self.resultLabels[i].setText('%#.5g' % (results[i]*np.sum(self.integralList)))

    def checkInputs(self):
        inp = safeEval(self.deltaEntry.text())
        if inp is None:
            return False
        self.deltaEntry.setText('%#.3g' % inp)
        inp = safeEval(self.etaEntry.text())
        if inp is None:
            return False
        self.etaEntry.setText('%#.3g' % inp)
        if np.unique(self.sidebandList).size != len(self.sidebandList):
            self.rootwindow.mainProgram.dispMsg("Multiple sidebands have the same index")
            return False
        for i in range(len(self.sidebandList)):
            inp = safeEval(self.integralEntries[i].text())
            self.integralList[i] = inp
            if inp is None:
                return False
            self.integralEntries[i].setText('%#.5g' % inp)
        return True
        
    def fit(self, *args):
        if len(self.integralList) < 2:
            self.rootwindow.mainProgram.dispMsg("Not enough integrals selected")
            return
        self.setCheng()
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        omegar = float(self.spinEntry.text())*1e3*np.pi*2
        if not self.deltaTick.isChecked():
            guess.append(float(self.deltaEntry.text())*1e-6)
            struc.append(True)
        else:
            outDelta = float(self.deltaEntry.text())*1e-6
            argu.append(outDelta)
            struc.append(False)
        if not self.etaTick.isChecked():
            guess.append(float(self.etaEntry.text()))
            struc.append(True)
        else:
            outEta = float(self.etaEntry.text())
            argu.append(outEta)
            struc.append(False)
        self.args = (struc, argu, self.parent.current.freq*np.pi*2)
        theta, phi, self.weight = zcw_angles(self.cheng, symm=2) #Theta and phi switched as algorithm actually uses alpha and beta as input.
        sinPhi = np.sin(phi)
        cosPhi = np.cos(phi)
        sin2Theta = np.sin(2*theta)
        cos2Theta = np.cos(2*theta)
        self.tresolution = 2*np.pi/omegar/self.NSTEPS
        self.t = np.linspace(0, self.tresolution*(self.NSTEPS-1), self.NSTEPS)
        cosOmegarT = np.cos(omegar*self.t)
        cos2OmegarT = np.cos(2*omegar*self.t)
        self.S1 = np.array([np.sqrt(2)/3 *  sinPhi * sin2Theta]).transpose()* np.sin(omegar*self.t)
        self.S2 =  np.array([cosPhi * sin2Theta/3]).transpose()* np.sin(2*omegar*self.t)
        self.C1 = np.array([np.sqrt(2)/3 *  sinPhi * cosPhi * 3 ]).transpose()* cosOmegarT
        self.C1eta = np.transpose([cos2Theta / 3.0]) * self.C1
        self.C2 = np.array([-1.0 / 3 *  3/2*sinPhi**2 ]).transpose()* cos2OmegarT
        self.C2eta = np.array([1.0 / 3 /2*(1+cosPhi**2)*cos2Theta ]).transpose()* cos2OmegarT
        fitVal = scipy.optimize.fmin(self.fitFunc, guess, args=(self.sidebandList, np.real(self.integralList)))
        counter = 0
        if struc[0]:
            self.deltaEntry.setText('%.3g' % (fitVal[counter]*1e6))
            outDelta = fitVal[counter]
            counter += 1
        if struc[1]:
            self.etaEntry.setText('%.3g' % fitVal[counter])
            outEta = fitVal[counter]
            counter += 1
        self.disp(outDelta, outEta)
 
    def sim(self):
        if len(self.integralList) < 2:
            self.rootwindow.mainProgram.dispMsg("Not enough integrals selected")
            return
        self.setCheng()
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        omegar = float(self.spinEntry.text())*1e3*np.pi*2
        outDelta = float(self.deltaEntry.text())*1e-6
        outEta = float(self.etaEntry.text())
        theta, phi, self.weight = zcw_angles(self.cheng, symm=2) #Theta and phi switched as algorithm actually uses alpha and beta as input.
        sinPhi = np.sin(phi)
        cosPhi = np.cos(phi)
        sin2Theta = np.sin(2*theta)
        cos2Theta = np.cos(2*theta)
        self.tresolution = 2*np.pi/omegar/self.NSTEPS
        self.t = np.linspace(0, self.tresolution*(self.NSTEPS-1), self.NSTEPS)
        cosOmegarT = np.cos(omegar*self.t)
        cos2OmegarT = np.cos(2*omegar*self.t)
        self.S1 = np.array([np.sqrt(2)/3 *  sinPhi * sin2Theta]).transpose()* np.sin(omegar*self.t)
        self.S2 =  np.array([cosPhi * sin2Theta/3]).transpose()* np.sin(2*omegar*self.t)
        self.C1 = np.array([np.sqrt(2)/3 *  sinPhi * cosPhi * 3 ]).transpose()* cosOmegarT
        self.C1eta = np.transpose([cos2Theta / 3.0]) * self.C1
        self.C2 = np.array([-1.0 / 3 *  3/2*sinPhi**2 ]).transpose()* cos2OmegarT
        self.C2eta = np.array([1.0 / 3 /2*(1+cosPhi**2)*cos2Theta ]).transpose()* cos2OmegarT
        self.disp(outDelta, outEta)
       
##############################################################################
class Quad1DeconvWindow(FittingWindow): 
    def __init__(self, mainProgram, oldMainWindow):
        FittingWindow.__init__(self, mainProgram, oldMainWindow)

    def setup(self):
        self.current = Quad1DeconvFrame(self, self.fig, self.canvas, self.oldMainWindow.current)
        self.paramframe = Quad1DeconvParamFrame(self.current, self)

#################################################################################   
class Quad1DeconvFrame(Plot1DFrame):
    def __init__(self, rootwindow, fig, canvas, current):
        Plot1DFrame.__init__(self, rootwindow, fig, canvas)
        self.data1D = current.getDisplayedData()
        self.current = current
        self.spec = self.current.spec
        self.xax = self.current.xax
        self.plotType = 0
        self.rootwindow = rootwindow
        self.pickNum = 0
        self.pickNum2 = 0
        self.plotReset()
        self.showPlot()

    def plotReset(self): #set the plot limits to min and max values
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
        self.yminlim = miny-differ
        self.ymaxlim = maxy+differ
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xminlim = min(self.xax*axMult)
        self.xmaxlim = max(self.xax*axMult)
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        self.ax.cla()
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.line_xdata = self.xax*axMult
        self.line_ydata = self.data1D
        self.ax.plot(self.xax*axMult, self.data1D)
        if tmpAx is not None:
            self.ax.plot(tmpAx*axMult, tmpdata)
        for i in range(len(tmpAx2)):
            self.ax.plot(tmpAx2[i]*axMult, tmpdata2[i])
        if self.spec == 0:
            if self.current.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec == 1:
            if self.current.ppm:
                self.ax.set_xlabel('Frequency [ppm]')
            else:
                if self.current.axType == 0:
                    self.ax.set_xlabel('Frequency [Hz]')
                elif self.current.axType == 1:
                    self.ax.set_xlabel('Frequency [kHz]')
                elif self.current.axType == 2:
                    self.ax.set_xlabel('Frequency [MHz]')
                else:
                    self.ax.set_xlabel('User defined')
        else:
            self.ax.set_xlabel('')
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()

#################################################################################
class Quad1DeconvParamFrame(QtGui.QWidget):
    
    Ioptions = ['1', '3/2', '2', '5/2', '3', '7/2', '4', '9/2']
    
    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        self.cheng = 15
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        if self.parent.current.spec == 1:
            self.axAdd = self.parent.current.freq-self.parent.current.ref
        elif self.parent.current.spec == 0:
            self.axAdd = 0
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1, 0, 0)
        grid.addLayout(self.optframe, 0, 1)
        grid.addLayout(self.frame2, 0, 2)
        grid.addLayout(self.frame3, 0, 3)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton, 0, 0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton, 1, 0)
        fitAllButton = QtGui.QPushButton("Fit all")
        fitAllButton.clicked.connect(self.fitAll)
        self.frame1.addWidget(fitAllButton, 2, 0)
        copyResultButton = QtGui.QPushButton("Copy result")
        copyResultButton.clicked.connect(lambda: self.sim(True))
        self.frame1.addWidget(copyResultButton, 3, 0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton, 4, 0)
        self.frame1.setColumnStretch(10, 1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.optframe.addWidget(QLabel("Cheng:"), 0, 0)
        self.chengEntry = QtGui.QLineEdit()
        self.chengEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.chengEntry.setText(str(self.cheng))
        self.optframe.addWidget(self.chengEntry, 1, 0)
        self.optframe.addWidget(QLabel("I:"), 0, 1)
        self.IEntry = QtGui.QComboBox()
        self.IEntry.addItems(self.Ioptions)
        self.IEntry.setCurrentIndex(1)
        self.optframe.addWidget(self.IEntry, 1, 1)
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel("Bgrnd:"), 0, 0, 1, 2)
        self.bgrndTick = QtGui.QCheckBox('')
        self.bgrndTick.setChecked(True)
        self.frame2.addWidget(self.bgrndTick, 1, 0)
        self.bgrndEntry = QtGui.QLineEdit()
        self.bgrndEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.bgrndEntry.setText("0.0")
        self.frame2.addWidget(self.bgrndEntry, 1, 1)
        self.frame2.addWidget(QLabel("Slope:"), 2, 0, 1, 2)
        self.slopeTick = QtGui.QCheckBox('')
        self.slopeTick.setChecked(True)
        self.frame2.addWidget(self.slopeTick, 3, 0)
        self.slopeEntry = QtGui.QLineEdit()
        self.slopeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.slopeEntry.setText("0.0")
        self.frame2.addWidget(self.slopeEntry, 3, 1)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        self.frame3.addWidget(QLabel("Pos:"), 1, 0, 1, 2)
        self.frame3.addWidget(QLabel("Cq [MHz]:"), 1, 2, 1, 2)
        self.frame3.addWidget(QLabel(u"\u03b7:"), 1, 4, 1, 2)
        self.frame3.addWidget(QLabel("integral:"), 1, 6, 1, 2)
        self.frame3.addWidget(QLabel("Lorentz [Hz]:"), 1, 8, 1, 2)
        self.frame3.addWidget(QLabel("Gauss [Hz]:"), 1, 10, 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.posEntries = []
        self.posTicks = []
        self.cqEntries = []
        self.cqTicks = []
        self.etaEntries = []
        self.etaTicks = []
        self.ampEntries = []
        self.ampTicks = []
        self.lorEntries = []
        self.lorTicks = []
        self.gaussEntries = []
        self.gaussTicks = []
        for i in range(10):
            self.posTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.posTicks[i], i+2, 0)
            self.posEntries.append(QtGui.QLineEdit())
            self.posEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.posEntries[i].setText("0.0")
            self.frame3.addWidget(self.posEntries[i], i+2, 1)
            self.cqTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.cqTicks[i], i+2, 2)
            self.cqEntries.append(QtGui.QLineEdit())
            self.cqEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.cqEntries[i].setText("0.0")
            self.frame3.addWidget(self.cqEntries[i], i+2, 3)
            self.etaTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.etaTicks[i], i+2, 4)
            self.etaEntries.append(QtGui.QLineEdit())
            self.etaEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.etaEntries[i].setText("0.0")
            self.frame3.addWidget(self.etaEntries[i], i+2, 5)
            self.ampTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.ampTicks[i], i+2, 6)
            self.ampEntries.append(QtGui.QLineEdit())
            self.ampEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.ampEntries[i].setText("1.0")
            self.frame3.addWidget(self.ampEntries[i], i+2, 7)
            self.lorTicks.append(QtGui.QCheckBox(''))
            self.lorTicks[i].setChecked(True)
            self.frame3.addWidget(self.lorTicks[i], i+2, 8)
            self.lorEntries.append(QtGui.QLineEdit())
            self.lorEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.lorEntries[i].setText("10.0")
            self.frame3.addWidget(self.lorEntries[i], i+2, 9)
            self.gaussTicks.append(QtGui.QCheckBox(''))
            self.gaussTicks[i].setChecked(True)
            self.frame3.addWidget(self.gaussTicks[i], i+2, 10)
            self.gaussEntries.append(QtGui.QLineEdit())
            self.gaussEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.gaussEntries[i].setText("0.0")
            self.frame3.addWidget(self.gaussEntries[i], i+2, 11)
            if i > 0:
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.cqTicks[i].hide()
                self.cqEntries[i].hide()
                self.etaTicks[i].hide()
                self.etaEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        
    def checkI(self, I):
        return I*0.5+1
                
    def setCheng(self, *args):
        inp = safeEval(self.chengEntry.text())
        if inp is None:
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        self.cheng = int(inp)
        self.chengEntry.setText(str(self.cheng))
                
    def changeNum(self, *args):
        val = self.numExp.currentIndex()+1
        for i in range(10):
            if i < val:
                self.posTicks[i].show()
                self.posEntries[i].show()
                self.cqTicks[i].show()
                self.cqEntries[i].show()
                self.etaTicks[i].show()
                self.etaEntries[i].show()
                self.ampTicks[i].show()
                self.ampEntries[i].show()
                self.lorTicks[i].show()
                self.lorEntries[i].show()
                self.gaussTicks[i].show()
                self.gaussEntries[i].show()
            else:
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.cqTicks[i].hide()
                self.cqEntries[i].hide()
                self.etaTicks[i].hide()
                self.etaEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()

    def tensorFunc(self, x, I, pos, cq, eta, width, gauss):
        m = np.arange(-I, I)
        v = []
        weights = []
        pos = pos - self.axAdd
        for i in m:
            tmp = (cq/(4*I*(2*I-1))*(I*(I+1)-3*(i+1)**2))-(cq/(4*I*(2*I-1))*(I*(I+1)-3*(i)**2))
            v = np.append(v, tmp*(self.angleStuff1-eta*self.angleStuff2)+pos)
            weights = np.append(weights, self.weight)
        length = len(x)
        t = np.arange(length)/self.parent.current.sw
        final = np.zeros(length)
        mult = v/(self.parent.current.sw)*length
        x1 = np.array(np.round(mult)+np.floor(length/2), dtype=int)
        weights = weights[np.logical_and(x1>=0, x1<length)]
        x1 = x1[np.logical_and(x1>=0, x1<length)]
        final = np.bincount(x1, weights, length)
        apod = np.exp(-np.pi*width*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
        apod[-1:-(len(apod)/2+1):-1] = apod[:len(apod)/2]
        inten = np.real(np.fft.fft(np.fft.ifft(final)*apod))
        inten = inten/self.parent.current.sw*len(inten)
        return inten
                
    def fitFunc(self, param, x, y):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        I = self.args[3]
        testFunc = np.zeros(len(self.parent.data1D))
        if struc[0]:
            bgrnd = param[0]
            param = np.delete(param, [0])
        else:
            bgrnd = argu[0]
            argu = np.delete(argu, [0])
        if struc[1]:
            slope = param[0]
            param = np.delete(param, [0])
        else:
            slope = argu[0]
            argu = np.delete(argu, [0])
        for i in range(numExp):
            if struc[6*i+2]:
                pos = param[0]
                param = np.delete(param, [0])
            else:
                pos = argu[0]
                argu = np.delete(argu, [0])
            if struc[6*i+3]:
                cq = param[0]
                param = np.delete(param, [0])
            else:
                cq = argu[0]
                argu = np.delete(argu, [0])
            if struc[6*i+4]:
                eta = param[0]
                param = np.delete(param, [0])
            else:
                eta = argu[0]
                argu = np.delete(argu, [0])
            if struc[6*i+5]:
                amp = param[0]
                param = np.delete(param, [0])
            else:
                amp = argu[0]
                argu = np.delete(argu, [0])
            if struc[6*i+6]:
                width = abs(param[0])
                param = np.delete(param, [0])
            else:
                width = argu[0]
                argu = np.delete(argu, [0])
            if struc[6*i+7]:
                gauss = abs(param[0])
                param = np.delete(param, [0])
            else:
                gauss = argu[0]
                argu = np.delete(argu, [0])
            testFunc += amp*self.tensorFunc(x, I, pos, cq, eta, width, gauss)
        testFunc += bgrnd+slope*x
        return np.sum((np.real(testFunc)-y)**2)

    def setAngleStuff(self):
        phi, theta, self.weight = zcw_angles(self.cheng, symm=2)
        self.angleStuff1 = 0.5*(3*np.cos(theta)**2-1)
        self.angleStuff2 = 0.5*np.cos(2*phi)*(np.sin(theta)**2)
        
    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.bgrndEntry.text())
        if inp is None:
            return False
        self.bgrndEntry.setText('%#.3g' % inp)
        inp = safeEval(self.slopeEntry.text())
        if inp is None:
            return False
        self.slopeEntry.setText('%#.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.posEntries[i].text())
            if inp is None:
                return False
            self.posEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.cqEntries[i].text())
            if inp is None:
                return False
            self.cqEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.etaEntries[i].text())
            if inp is None:
                return False
            self.etaEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.ampEntries[i].text())
            if inp is None:
                return False
            self.ampEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.lorEntries[i].text())
            if inp is None:
                return False
            self.lorEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.gaussEntries[i].text())
            if inp is None:
                return False
            self.gaussEntries[i].setText('%#.3g' % inp)
        return True
    
    def fit(self, *args):
        self.setCheng()
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outCq = np.zeros(numExp)
        outEta = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        I = self.checkI(self.IEntry.currentIndex())
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = float(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            if not self.posTicks[i].isChecked():
                guess.append(float(self.posEntries[i].text()))
                struc.append(True)
            else:
                outPos[i] = float(self.posEntries[i].text())
                argu.append(outPos[i])
                struc.append(False)
            if not self.cqTicks[i].isChecked():
                guess.append(float(self.cqEntries[i].text())*1e6)
                struc.append(True)
            else:
                outCq[i] = float(self.cqEntries[i].text())*1e6
                argu.append(outCq[i])
                struc.append(False)
            if not self.etaTicks[i].isChecked():
                guess.append(float(self.etaEntries[i].text()))
                struc.append(True)
            else:
                outEta[i] = float(self.etaEntries[i].text())
                argu.append(outEta[i])
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].isChecked():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp, struc, argu, I)
        self.setAngleStuff()
        fitVal = scipy.optimize.fmin(self.fitFunc, guess, args=(self.parent.xax, np.real(self.parent.data1D)))
        counter = 0
        if struc[0]:
            self.bgrndEntry.setText('%.3g' % fitVal[counter])
            outBgrnd = fitVal[counter]
            counter += 1
        if struc[1]:
            self.slopeEntry.setText('%.3g' % fitVal[counter])
            outSlope = fitVal[counter]
            counter += 1
        for i in range(numExp):
            if struc[6*i+2]:
                self.posEntries[i].setText('%.3g' % fitVal[counter])
                outPos[i] = fitVal[counter]
                counter += 1
            if struc[6*i+3]:
                self.cqEntries[i].setText('%.3g' % (fitVal[counter]*1e-6))
                outCq[i] = fitVal[counter]
                counter += 1
            if struc[6*i+4]:
                self.etaEntries[i].setText('%.3g' % abs(fitVal[counter]))
                outEta[i] = fitVal[counter]
                counter += 1
            if struc[6*i+5]:
                self.ampEntries[i].setText('%.3g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[6*i+6]:
                self.lorEntries[i].setText('%.3g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
            if struc[6*i+7]:
                self.gaussEntries[i].setText('%.3g' % abs(fitVal[counter]))
                outGauss[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd, outSlope, I, outPos, outCq, outEta, outAmp, outWidth, outGauss)

    def fitAll(self, *args):
        FitAllSelectionWindow(self, ["Background", "Slope", "Position", "Cq", u"\u03b7", "Integral", "Lorentz", "Gauss"])

    def fitAllFunc(self, outputs):
        self.setCheng()
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outCq = np.zeros(numExp)
        outEta = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        I = self.checkI(self.IEntry.currentIndex())
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = float(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            if not self.posTicks[i].isChecked():
                guess.append(float(self.posEntries[i].text()))
                struc.append(True)
            else:
                outPos[i] = float(self.posEntries[i].text())
                argu.append(outPos[i])
                struc.append(False)
            if not self.cqTicks[i].isChecked():
                guess.append(float(self.cqEntries[i].text())*1e6)
                struc.append(True)
            else:
                outCq[i] = float(self.cqEntries[i].text())*1e6
                argu.append(outCq[i])
                struc.append(False)
            if not self.etaTicks[i].isChecked():
                guess.append(float(self.etaEntries[i].text()))
                struc.append(True)
            else:
                outEta[i] = float(self.etaEntries[i].text())
                argu.append(outEta[i])
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].isChecked():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp, struc, argu, I)
        self.setAngleStuff()
        fullData = self.parent.current.data.data
        axes = self.parent.current.axes
        dataShape = fullData.shape
        dataShape2 = np.delete(dataShape, axes)
        rolledData = np.rollaxis(fullData, axes)
        intOutputs = np.array(outputs, dtype=int)
        numOutputs = np.sum(intOutputs[:2]) + numExp*np.sum(intOutputs[2:])
        outputData = np.zeros((np.product(dataShape2), numOutputs), dtype=complex)
        counter2 = 0
        for j in rolledData.reshape(dataShape[axes], np.product(dataShape2)).T:
            try:
                fitVal = scipy.optimize.fmin(self.fitFunc, guess, args=(self.parent.xax, np.real(j)))
            except:
                fitVal = [[0]*10]
            counter = 0
            if struc[0]:
                outBgrnd = fitVal[counter]
                counter += 1
            if struc[1]:
                outSlope = fitVal[counter]
                counter += 1
            for i in range(numExp):
                if struc[6*i+2]:
                    outPos[i] = fitVal[counter]
                    counter += 1
                if struc[6*i+3]:
                    outCq[i] = fitVal[counter]
                    counter += 1
                if struc[6*i+4]:
                    outEta[i] = fitVal[counter]
                    counter += 1
                if struc[6*i+5]:
                    outAmp[i] = fitVal[counter]
                    counter += 1
                if struc[6*i+6]:
                    outWidth[i] = abs(fitVal[counter])
                    counter += 1
                if struc[6*i+7]:
                    outGauss[i] = abs(fitVal[counter])
                    counter += 1
            outputArray = []
            if outputs[0]:
                outputArray = np.concatenate((outputArray, [outBgrnd]))
            if outputs[1]:
                outputArray = np.concatenate((outputArray, [outSlope]))
            if outputs[2]:
                outputArray = np.concatenate((outputArray, outPos))
            if outputs[3]:
                outputArray = np.concatenate((outputArray, outCq))
            if outputs[4]:
                outputArray = np.concatenate((outputArray, outEta))
            if outputs[5]:
                outputArray = np.concatenate((outputArray, outAmp)) 
            if outputs[6]:
                outputArray = np.concatenate((outputArray, outWidth))
            if outputs[7]:
                outputArray = np.concatenate((outputArray, outGauss))
            outputData[counter2] = outputArray
            counter2 += 1
        newShape = np.concatenate((np.array(dataShape2), [numOutputs]))
        self.rootwindow.createNewData(np.rollaxis(outputData.reshape(newShape), -1, axes), axes)
        
    def sim(self, store=False):
        self.setCheng()
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        numExp = self.numExp.currentIndex()+1
        bgrnd = float(self.bgrndEntry.text())
        slope = float(self.slopeEntry.text())
        pos = np.zeros(numExp)
        cq = np.zeros(numExp)
        eta = np.zeros(numExp)
        amp = np.zeros(numExp)
        width = np.zeros(numExp)
        gauss = np.zeros(numExp)
        I = self.checkI(self.IEntry.currentIndex())
        for i in range(numExp):
            pos[i] = float(self.posEntries[i].text())
            cq[i] = float(self.cqEntries[i].text())*1e6
            eta[i] = float(self.etaEntries[i].text())
            amp[i] = float(self.ampEntries[i].text())
            width[i] = float(self.lorEntries[i].text())
            gauss[i] = float(self.gaussEntries[i].text())
        self.setAngleStuff()
        self.disp(bgrnd, slope, I, pos, cq, eta, amp, width, gauss, store)

    def disp(self, outBgrnd, outSlope, outI, outPos, outCq, outEta, outAmp, outWidth, outGauss, store=False):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x = []
        for i in range(len(outPos)):
            x.append(tmpx)
            y =  outAmp[i]*self.tensorFunc(tmpx, outI, outPos[i], outCq[i], outEta[i], outWidth[i], outGauss[i])
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        if store:
            outCurvePart.append(outCurve)
            outCurvePart.append(self.parent.data1D)
            self.rootwindow.createNewData(np.array(outCurvePart), self.parent.current.axes, True)
        else:
            self.parent.showPlot(tmpx, outCurve, x, outCurvePart)
        
##############################################################################
class Quad2DeconvWindow(FittingWindow): 
    def __init__(self, mainProgram, oldMainWindow, mas=False):
        self.mas = mas
        FittingWindow.__init__(self, mainProgram, oldMainWindow)

    def setup(self):
        self.current = Quad1DeconvFrame(self, self.fig, self.canvas, self.oldMainWindow.current)
        if self.mas:
            self.paramframe = Quad2MASDeconvParamFrame(self.current, self)
        else:
            self.paramframe = Quad2StaticDeconvParamFrame(self.current, self)
         
#################################################################################
class Quad2StaticDeconvParamFrame(Quad1DeconvParamFrame): 

    Ioptions = ['3/2', '5/2', '7/2', '9/2']
    
    def __init__(self, parent, rootwindow):
        Quad1DeconvParamFrame.__init__(self, parent, rootwindow)
        
    def checkI(self, I):
        return I*1.0+1.5
    
    def setAngleStuff(self):
        phi, theta, self.weight = zcw_angles(self.cheng, symm=2)
        self.angleStuff1 = -27/8.0*np.cos(theta)**4+15/4.0*np.cos(theta)**2-3/8.0
        self.angleStuff2 = (-9/4.0*np.cos(theta)**4+2*np.cos(theta)**2+1/4.0)*np.cos(2*phi)
        self.angleStuff3 = -1/2.0*np.cos(theta)**2+1/3.0+(-3/8.0*np.cos(theta)**4+3/4.0*np.cos(theta)**2-3/8.0)*np.cos(2*phi)**2

    def tensorFunc(self, x, I, pos, cq, eta, width, gauss):
        pos = pos - self.axAdd
        # v = -cq**2/(6.0*self.parent.current.freq)*(I*(I+1)-3/4.0)*(self.angleStuff1+self.angleStuff2*eta+self.angleStuff3*eta**2)+pos
        v = -1/(6*self.parent.current.freq)*(3*cq/(2*I*(2*I-1)))**2*(I*(I+1)-3.0/4)*(self.angleStuff1+self.angleStuff2*eta+self.angleStuff3*eta**2)+pos
        length = len(x)
        t = np.arange(length)/self.parent.current.sw
        final = np.zeros(length)
        mult = v/(self.parent.current.sw)*length
        x1 = np.array(np.round(mult)+np.floor(length/2), dtype=int)
        weights = self.weight[np.logical_and(x1>=0, x1<length)]
        x1 = x1[np.logical_and(x1>=0, x1<length)]
        final = np.bincount(x1, weights, length)
        apod = np.exp(-np.pi*width*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
        apod[-1:-(len(apod)/2+1):-1] = apod[:len(apod)/2]
        inten = np.real(np.fft.fft(np.fft.ifft(final)*apod))
        inten = inten/self.parent.current.sw/len(inten)
        return inten
    
#################################################################################
class Quad2MASDeconvParamFrame(Quad2StaticDeconvParamFrame): 
    Ioptions = ['3/2', '5/2', '7/2', '9/2']
    
    def __init__(self, parent, rootwindow):
        Quad2StaticDeconvParamFrame.__init__(self, parent, rootwindow)
        
    def setAngleStuff(self):
        phi, theta, self.weight = zcw_angles(self.cheng, symm=2)
        self.angleStuff1 = 21/16.0*np.cos(theta)**4-9/8.0*np.cos(theta)**2+5/16.0
        self.angleStuff2 = (-7/8.0*np.cos(theta)**4+np.cos(theta)**2-1/8.0)*np.cos(2*phi)
        self.angleStuff3 = 1/12.0*np.cos(theta)**2+(+7/48.0*np.cos(theta)**4-7/24.0*np.cos(theta)**2+7/48.0)*np.cos(2*phi)**2

##############################################################################
class Quad2CzjzekWindow(FittingWindow): 
    def __init__(self, mainProgram, oldMainWindow, mas=False):
        self.mas = mas
        FittingWindow.__init__(self, mainProgram, oldMainWindow)

    def setup(self):
        self.current = Quad1DeconvFrame(self, self.fig, self.canvas, self.oldMainWindow.current)
        if self.mas:
            self.paramframe = Quad2MASCzjzekParamFrame(self.current, self)
        else:
            self.paramframe = Quad2StaticCzjzekParamFrame(self.current, self)

#################################################################################
class Quad2StaticCzjzekParamFrame(QtGui.QWidget): 

    Ioptions = ['3/2', '5/2', '7/2', '9/2']

    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        self.cheng = 15
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        if self.parent.current.spec == 1:
            self.axAdd = self.parent.current.freq-self.parent.current.ref
        elif self.parent.current.spec == 0:
            self.axAdd = 0
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1, 0, 0)
        grid.addLayout(self.optframe, 0, 1)
        grid.addLayout(self.frame2, 0, 2)
        grid.addLayout(self.frame3, 0, 3)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton, 0, 0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton, 1, 0)
        fitAllButton = QtGui.QPushButton("Fit all")
        fitAllButton.clicked.connect(self.fitAll)
        self.frame1.addWidget(fitAllButton, 2, 0)
        copyResultButton = QtGui.QPushButton("Copy result")
        copyResultButton.clicked.connect(lambda: self.sim(True))
        self.frame1.addWidget(copyResultButton, 3, 0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton, 4, 0)
        self.frame1.setColumnStretch(10, 1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.optframe.addWidget(QLabel("Cheng:"), 0, 0)
        self.chengEntry = QtGui.QLineEdit()
        self.chengEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.chengEntry.setText(str(self.cheng))
        self.optframe.addWidget(self.chengEntry, 1, 0)
        self.optframe.addWidget(QLabel("I:"), 0, 1)
        self.IEntry = QtGui.QComboBox()
        self.IEntry.addItems(self.Ioptions)
        self.IEntry.setCurrentIndex(1)
        self.optframe.addWidget(self.IEntry, 1, 1)
        self.optframe.addWidget(QLabel("Cq grid size:"), 2, 0)
        self.cqGridEntry = QtGui.QLineEdit()
        self.cqGridEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.cqGridEntry.setText("50")
        self.cqGridEntry.returnPressed.connect(self.setGrid)
        self.optframe.addWidget(self.cqGridEntry, 3, 0)
        self.optframe.addWidget(QLabel(u"\u03b7 grid size:"), 4, 0)
        self.etaGridEntry = QtGui.QLineEdit()
        self.etaGridEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.etaGridEntry.setText("10")
        self.etaGridEntry.returnPressed.connect(self.setGrid)
        self.optframe.addWidget(self.etaGridEntry, 5, 0)
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel("Bgrnd:"), 0, 0, 1, 2)
        self.bgrndTick = QtGui.QCheckBox('')
        self.bgrndTick.setChecked(True)
        self.frame2.addWidget(self.bgrndTick, 1, 0)
        self.bgrndEntry = QtGui.QLineEdit()
        self.bgrndEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.bgrndEntry.setText("0.0")
        self.frame2.addWidget(self.bgrndEntry, 1, 1)
        self.frame2.addWidget(QLabel("Slope:"), 2, 0, 1, 2)
        self.slopeTick = QtGui.QCheckBox('')
        self.slopeTick.setChecked(True)
        self.frame2.addWidget(self.slopeTick, 3, 0)
        self.slopeEntry = QtGui.QLineEdit()
        self.slopeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.slopeEntry.setText("0.0")
        self.frame2.addWidget(self.slopeEntry, 3, 1)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0)
        self.frame3.addWidget(QLabel("d:"), 1, 0)
        self.frame3.addWidget(QLabel("Pos:"), 1, 1, 1, 2)
        self.frame3.addWidget(QLabel(u"\u03c3 [MHz]:"), 1, 3, 1, 2)
        self.frame3.addWidget(QLabel("Integral:"), 1, 5, 1, 2)
        self.frame3.addWidget(QLabel("Lorentz [Hz]:"), 1, 7, 1, 2)
        self.frame3.addWidget(QLabel("Gauss [Hz]:"), 1, 9, 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.posEntries = []
        self.posTicks = []
        self.sigmaEntries = []
        self.sigmaTicks = []
        self.dEntries = []
        self.ampEntries = []
        self.ampTicks = []
        self.lorEntries = []
        self.lorTicks = []
        self.gaussEntries = []
        self.gaussTicks = []
        for i in range(10):
            self.dEntries.append(QtGui.QLineEdit())
            self.dEntries[i].setText("5")
            self.frame3.addWidget(self.dEntries[i], i+2, 0)
            self.posTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.posTicks[i], i+2, 1)
            self.posEntries.append(QtGui.QLineEdit())
            self.posEntries[i].setText("0.0")
            self.frame3.addWidget(self.posEntries[i], i+2, 2)
            self.sigmaTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.sigmaTicks[i], i+2, 3)
            self.sigmaEntries.append(QtGui.QLineEdit())
            self.sigmaEntries[i].setText("1.0")
            self.frame3.addWidget(self.sigmaEntries[i], i+2, 4)
            self.ampTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.ampTicks[i], i+2, 5)
            self.ampEntries.append(QtGui.QLineEdit())
            self.ampEntries[i].setText("1.0")
            self.frame3.addWidget(self.ampEntries[i], i+2, 6)
            self.lorTicks.append(QtGui.QCheckBox(''))
            self.lorTicks[i].setChecked(True)
            self.frame3.addWidget(self.lorTicks[i], i+2, 7)
            self.lorEntries.append(QtGui.QLineEdit())
            self.lorEntries[i].setText("10.0")
            self.frame3.addWidget(self.lorEntries[i], i+2, 8)
            self.gaussTicks.append(QtGui.QCheckBox(''))
            self.gaussTicks[i].setChecked(True)
            self.frame3.addWidget(self.gaussTicks[i], i+2, 9)
            self.gaussEntries.append(QtGui.QLineEdit())
            self.gaussEntries[i].setText("0.0")
            self.frame3.addWidget(self.gaussEntries[i], i+2, 10)

            if i > 0:
                self.dEntries[i].hide()
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.sigmaTicks[i].hide()
                self.sigmaEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)

    def checkI(self, I):
        return I*1.0+1.5
                
    def setCheng(self, *args):
        inp = safeEval(self.chengEntry.text())
        if inp is None:
            self.cheng = 15
        else:
            self.cheng = int(inp)
        self.chengEntry.setText(str(self.cheng))

    def setGrid(self, *args):
        inp = safeEval(self.cqGridEntry.text())
        if inp is None:
            return False
        self.cqGridEntry.setText(str(int(inp)))
        inp = safeEval(self.etaGridEntry.text())
        if inp is None:
            return False
        self.etaGridEntry.setText(str(int(inp)))
        return True
        
    def changeNum(self, *args):
        val = self.numExp.currentIndex()+1
        for i in range(10):
            if i < val:
                self.dEntries[i].show()
                self.posTicks[i].show()
                self.posEntries[i].show()
                self.sigmaTicks[i].show()
                self.sigmaEntries[i].show()
                self.ampTicks[i].show()
                self.ampEntries[i].show()
                self.lorTicks[i].show()
                self.lorEntries[i].show()
                self.gaussTicks[i].show()
                self.gaussEntries[i].show()
            else:
                self.dEntries[i].hide()
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.sigmaTicks[i].hide()
                self.sigmaEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()

    def bincounting(self, x1, weight, length):
        weights = weight[np.logical_and(x1>=0, x1<length)]
        x1 = x1[np.logical_and(x1>=0, x1<length)]
        return np.fft.ifft(np.bincount(x1, weights, length))
    
    def genLib(self, length, I, maxCq, numCq, numEta):
        self.cq, self.eta = np.meshgrid(np.linspace(0, maxCq, numCq), np.linspace(0, 1, numEta))
        cq = self.cq[..., None]
        eta = self.eta[..., None]
       # v = -cq**2/(6.0*self.parent.current.freq)*(I*(I+1)-3/4.0)*(self.angleStuff1+self.angleStuff2*eta+self.angleStuff3*eta**2)
        v = -1/(6.0*self.parent.current.freq)*(3*cq/(2*I*(2*I-1)))**2*(I*(I+1)-3.0/4)*(self.angleStuff1+self.angleStuff2*eta+self.angleStuff3*eta**2)
        mult = v/(self.parent.current.sw)*length
        x1 = np.array(np.round(mult)+np.floor(length/2), dtype=int)
        self.lib = np.apply_along_axis(self.bincounting, 2, x1, self.weight, length)

    def tensorFunc(self, sigma, d, pos, width, gauss):
        import scipy.ndimage
        pos = pos - self.axAdd
        cq = self.cq
        eta = self.eta
        czjzek = cq**(d-1)*eta/(np.sqrt(2*np.pi)*sigma**d)*(1-eta**2/9.0)*np.exp(-cq**2/(2.0*sigma**2)*(1+eta**2/3.0))
        czjzek = czjzek/np.sum(czjzek)
        fid = np.sum(self.lib*czjzek[..., None], axis=(0, 1))
        t = np.arange(len(fid))/self.parent.current.sw
        apod = np.exp(-np.pi*width*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
        apod[-1:-(len(apod)/2+1):-1] = apod[:len(apod)/2]
        spectrum = scipy.ndimage.interpolation.shift(np.real(np.fft.fft(fid*apod)), len(fid)*pos/self.parent.current.sw)
        spectrum = spectrum/self.parent.current.sw*len(spectrum)
        return spectrum
    
    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.bgrndEntry.text())
        if inp is None:
            return False
        self.bgrndEntry.setText('%#.3g' % inp)
        inp = safeEval(self.slopeEntry.text())
        if inp is None:
            return False
        self.slopeEntry.setText('%#.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.dEntries[i].text())
            if inp is None:
                return False
            if inp < 1:
                inp = 1
            elif inp > 5:
                inp = 5
            self.dEntries[i].setText(str(int(inp)))
            inp = safeEval(self.posEntries[i].text())
            if inp is None:
                return False
            self.posEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.sigmaEntries[i].text())
            if inp is None:
                return False
            self.sigmaEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.ampEntries[i].text())
            if inp is None:
                return False
            self.ampEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.lorEntries[i].text())
            if inp is None:
                return False
            self.lorEntries[i].setText('%#.3g' % inp)
            inp = safeEval(self.gaussEntries[i].text())
            if inp is None:
                return False
            self.gaussEntries[i].setText('%#.3g' % inp)
        return True
    
    def fitFunc(self, param, x, y):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        testFunc = np.zeros(len(self.parent.data1D))
        if struc[0]:
            bgrnd = param[0]
            param = np.delete(param, [0])
        else:
            bgrnd = argu[0]
            argu = np.delete(argu, [0])
        if struc[1]:
            slope = param[0]
            param = np.delete(param, [0])
        else:
            slope = argu[0]
            argu = np.delete(argu, [0])
        for i in range(numExp):
            d = argu[0]
            argu = np.delete(argu, [0])
            if struc[5*i+2]:
                pos = param[0]
                param = np.delete(param, [0])
            else:
                pos = argu[0]
                argu = np.delete(argu, [0])
            if struc[5*i+3]:
                sigma = param[0]
                param = np.delete(param, [0])
            else:
                sigma = argu[0]
                argu = np.delete(argu, [0])
            if struc[5*i+4]:
                amp = param[0]
                param = np.delete(param, [0])
            else:
                amp = argu[0]
                argu = np.delete(argu, [0])
            if struc[5*i+5]:
                width = abs(param[0])
                param = np.delete(param, [0])
            else:
                width = argu[0]
                argu = np.delete(argu, [0])
            if struc[5*i+6]:
                gauss = abs(param[0])
                param = np.delete(param, [0])
            else:
                gauss = argu[0]
                argu = np.delete(argu, [0])
            testFunc += amp*self.tensorFunc(sigma, d, pos, width, gauss)
        testFunc += bgrnd+slope*x
        return np.sum((np.real(testFunc)-y)**2)

    def setAngleStuff(self):
        phi, theta, self.weight = zcw_angles(self.cheng, symm=2)
        self.angleStuff1 = -27/8.0*np.cos(theta)**4+15/4.0*np.cos(theta)**2-3/8.0
        self.angleStuff2 = (-9/4.0*np.cos(theta)**4+2*np.cos(theta)**2+1/4.0)*np.cos(2*phi)
        self.angleStuff3 = -1/2.0*np.cos(theta)**2+1/3.0+(-3/8.0*np.cos(theta)**4+3/4.0*np.cos(theta)**2-3/8.0)*np.cos(2*phi)**2
        
    def fit(self, *args):
        self.setCheng()
        if not self.setGrid():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        maxCq = 0.0
        I = self.checkI(self.IEntry.currentIndex())
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outSigma = np.zeros(numExp)
        outD = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = safeEval(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            inp = int(self.dEntries[i].text())
            if inp < 1:
                inp = 1
            elif inp > 5:
                inp = 5
            argu.append(inp)
            outD[i] = inp
            self.dEntries[i].setText('%.3g' % inp)
            if not self.posTicks[i].isChecked():
                guess.append(float(self.posEntries[i].text()))
                struc.append(True)
            else:
                outPos[i] = float(self.posEntries[i].text())
                argu.append(outPos[i])
                struc.append(False)
            if not self.sigmaTicks[i].isChecked():
                inp = float(self.sigmaEntries[i].text())
                maxCq = max(maxCq, inp*1e6)
                guess.append(inp*1e6)
                struc.append(True)
            else:
                inp = float(self.sigmaEntries[i].text())
                maxCq = max(maxCq, inp*1e6)
                argu.append(inp*1e6)
                outSigma[i] = inp
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].isChecked():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp, struc, argu)
        self.setAngleStuff()
        numCq = int(self.cqGridEntry.text())
        numEta = int(self.etaGridEntry.text())
        self.genLib(len(self.parent.xax), I, maxCq*4.0, numCq, numEta)
        fitVal = scipy.optimize.fmin(self.fitFunc, guess, args=(self.parent.xax, np.real(self.parent.data1D)))
        counter = 0
        if struc[0]:
            self.bgrndEntry.setText('%.3g' % fitVal[counter])
            outBgrnd = fitVal[counter]
            counter += 1
        if struc[1]:
            self.slopeEntry.setText('%.3g' % fitVal[counter])
            outSlope = fitVal[counter]
            counter += 1
        for i in range(numExp):
            if struc[5*i+2]:
                self.posEntries[i].setText('%.3g' % fitVal[counter])
                outPos[i] = fitVal[counter]
                counter += 1
            if struc[5*i+3]:
                self.sigmaEntries[i].setText('%.3g' % (fitVal[counter]*1e-6))
                outSigma[i] = fitVal[counter]
                counter += 1
            if struc[5*i+4]:
                self.ampEntries[i].setText('%.3g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[5*i+5]:
                self.lorEntries[i].setText('%.3g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
            if struc[5*i+6]:
                self.gaussEntries[i].setText('%.3g' % abs(fitVal[counter]))
                outGauss[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd, outSlope, outPos, outSigma, outD, outAmp, outWidth, outGauss)

    def fitAll(self, *args):
        FitAllSelectionWindow(self, ["Background", "Slope", "Position", u"\u03c3", "Integral", "Lorentz", "Gauss"])

    def fitAllFunc(self, outputs):
        self.setCheng()
        if not self.setGrid():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        maxCq = 0.0
        I = self.checkI(self.IEntry.currentIndex())
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outSigma = np.zeros(numExp)
        outD = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = safeEval(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            inp = int(self.dEntries[i].text())
            if inp < 1:
                inp = 1
            elif inp > 5:
                inp = 5
            argu.append(inp)
            outD[i] = inp
            self.dEntries[i].setText('%.3g' % inp)
            if not self.posTicks[i].isChecked():
                guess.append(float(self.posEntries[i].text()))
                struc.append(True)
            else:
                outPos[i] = float(self.posEntries[i].text())
                argu.append(outPos[i])
                struc.append(False)
            if not self.sigmaTicks[i].isChecked():
                inp = float(self.sigmaEntries[i].text())
                maxCq = max(maxCq, inp*1e6)
                guess.append(inp*1e6)
                struc.append(True)
            else:
                inp = float(self.sigmaEntries[i].text())
                maxCq = max(maxCq, inp*1e6)
                argu.append(inp*1e6)
                outSigma[i] = inp
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].text():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp, struc, argu)
        self.setAngleStuff()
        numCq = int(self.cqGridEntry.text())
        numEta = int(self.etaGridEntry.text())
        self.genLib(len(self.parent.xax), I, maxCq*4.0, numCq, numEta)
        fullData = self.parent.current.data.data
        axes = self.parent.current.axes
        dataShape = fullData.shape
        dataShape2 = np.delete(dataShape, axes)
        rolledData = np.rollaxis(fullData, axes)
        intOutputs = np.array(outputs, dtype=int)
        numOutputs = np.sum(intOutputs[:2]) + numExp*np.sum(intOutputs[2:])
        outputData = np.zeros((np.product(dataShape2), numOutputs), dtype=complex)
        counter2 = 0
        for j in rolledData.reshape(dataShape[axes], np.product(dataShape2)).T:
            try:
                fitVal = scipy.optimize.fmin(self.fitFunc, guess, args=(self.parent.xax, np.real(j)))
            except:
                fitVal = [[0]*10]
            counter = 0
            if struc[0]:
                outBgrnd = fitVal[counter]
                counter += 1
            if struc[1]:
                outSlope = fitVal[counter]
                counter += 1
            for i in range(numExp):
                if struc[5*i+2]:
                    outPos[i] = fitVal[counter]
                    counter += 1
                if struc[5*i+3]:
                    outSigma[i] = fitVal[counter]
                    counter += 1
                if struc[5*i+4]:
                    outAmp[i] = fitVal[counter]
                    counter += 1
                if struc[5*i+5]:
                    outWidth[i] = abs(fitVal[counter])
                    counter += 1
                if struc[5*i+6]:
                    outGauss[i] = abs(fitVal[counter])
                    counter += 1
            outputArray = []
            if outputs[0]:
                outputArray = np.concatenate((outputArray, [outBgrnd]))
            if outputs[1]:
                outputArray = np.concatenate((outputArray, [outSlope]))
            if outputs[2]:
                outputArray = np.concatenate((outputArray, outPos))
            if outputs[3]:
                outputArray = np.concatenate((outputArray, outSigma))
            if outputs[4]:
                outputArray = np.concatenate((outputArray, outAmp)) 
            if outputs[5]:
                outputArray = np.concatenate((outputArray, outWidth))
            if outputs[6]:
                outputArray = np.concatenate((outputArray, outGauss))
            outputData[counter2] = outputArray
            counter2 += 1
        newShape = np.concatenate((np.array(dataShape2), [numOutputs]))
        self.rootwindow.createNewData(np.rollaxis(outputData.reshape(newShape), -1, axes), axes)
        
    def sim(self, store=False):
        self.setCheng()
        if not self.setGrid():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        if not self.checkInputs():
            self.rootwindow.mainProgram.dispMsg("One of the inputs is not valid")
            return
        numExp = self.numExp.currentIndex()+1
        bgrnd = float(self.bgrndEntry.text())
        slope = float(self.slopeEntry.text())
        pos = np.zeros(numExp)
        sigma = np.zeros(numExp)
        d = np.zeros(numExp)
        amp = np.zeros(numExp)
        width = np.zeros(numExp)
        gauss = np.zeros(numExp)
        I = self.checkI(self.IEntry.currentIndex())
        for i in range(numExp):
            pos[i] = safeEval(self.posEntries[i].text())
            sigma[i] = safeEval(self.sigmaEntries[i].text())*1e6
            d[i] = safeEval(self.dEntries[i].text())
            amp[i] = safeEval(self.ampEntries[i].text())
            width[i] = safeEval(self.lorEntries[i].text())
            gauss[i] = safeEval(self.gaussEntries[i].text())
        self.setAngleStuff()
        numCq = int(self.cqGridEntry.text())
        numEta = int(self.etaGridEntry.text())
        self.genLib(len(self.parent.xax), I, max(sigma)*4.0, numCq, numEta)
        self.disp(bgrnd, slope, pos, sigma, d, amp, width, gauss, store)


    def disp(self, outBgrnd, outSlope, outPos, outSigma, outD, outAmp, outWidth, outGauss, store=False):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x = []
        for i in range(len(outPos)):
            x.append(tmpx)
            y =  outAmp[i]*self.tensorFunc(outSigma[i], outD[i], outPos[i], outWidth[i], outGauss[i])
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        if store:
            outCurvePart.append(outCurve)
            outCurvePart.append(self.parent.data1D)
            self.rootwindow.createNewData(np.array(outCurvePart), self.parent.current.axes, True)
        else:
            self.parent.showPlot(tmpx, outCurve, x, outCurvePart)
        self.parent.showPlot(tmpx, outCurve, x, outCurvePart)
        
#################################################################################
class Quad2MASCzjzekParamFrame(Quad2StaticCzjzekParamFrame): 

    Ioptions = ['3/2', '5/2', '7/2', '9/2']
    
    def __init__(self, parent, rootwindow):
        Quad2StaticCzjzekParamFrame.__init__(self, parent, rootwindow)
        
    def setAngleStuff(self):
        phi, theta, self.weight = zcw_angles(self.cheng, symm=2)
        self.angleStuff1 = 21/16.0*np.cos(theta)**4-9/8.0*np.cos(theta)**2+5/16.0
        self.angleStuff2 = (-7/8.0*np.cos(theta)**4+np.cos(theta)**2-1/8.0)*np.cos(2*phi)
        self.angleStuff3 = 1/12.0*np.cos(theta)**2+(+7/48.0*np.cos(theta)**4-7/24.0*np.cos(theta)**2+7/48.0)*np.cos(2*phi)**2
         
######################################################################

class FitAllSelectionWindow(QtGui.QWidget): 
    def __init__(self, parent, fitNames):
        QtGui.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Select output")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        
        self.ticks = []
        for i in range(len(fitNames)):
            self.ticks.append(QtGui.QCheckBox(fitNames[i]))
            grid.addWidget(self.ticks[i], i, 0)
        
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.fit)
        layout.addWidget(okButton, 1, 1)
        self.show()
        self.setFixedSize(self.size())
        self.setGeometry(self.frameSize().width()-self.geometry().width(), self.frameSize().height()-self.geometry().height(), 0, 0)
        
    def fit(self):
        returnVals = []
        for i in self.ticks:
            returnVals.append(i.isChecked())
        self.deleteLater()
        self.father.fitAllFunc(np.array(returnVals, dtype=bool))
