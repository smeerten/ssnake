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
try:
    from PyQt4 import QtGui, QtCore
    from PyQt4 import QtGui as QtWidgets
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
except ImportError:
    from PyQt5 import QtGui, QtCore, QtWidgets
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import scipy.optimize
import multiprocessing
import copy
import re
import time
import tempfile
import os
import subprocess
import shutil
from safeEval import safeEval
from views import Current1D, CurrentContour
import widgetClasses as wc
import functions as func
import simFunctions as simFunc
import specIO as io
from ssNake import SideFrame
import extendedCzjzek as eCzjz

pi = np.pi
stopDict = {}  # Global dictionary with stopping commands for fits


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def checkLinkTuple(inp):
    if len(inp) == 2:
        inp += (1, 0, 0)
    elif len(inp) == 3:
        inp += (0, 0)
    return inp

#############################################################################################


class TabFittingWindow(QtWidgets.QWidget):

    PRECIS = 4
    MINMETHOD = 'Powell'
    NUMFEVAL = 150

    def __init__(self, father, oldMainWindow):
        super(TabFittingWindow, self).__init__(father)
        self.father = father
        self.oldMainWindow = oldMainWindow
        self.subFitWindows = []
        self.tabs = QtWidgets.QTabWidget(self)
        self.tabs.setTabPosition(2)
        self.mainFitWindow = FittingWindow(father, oldMainWindow, self)
        self.current = self.mainFitWindow.current
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.closeTab)
        self.tabs.addTab(self.mainFitWindow, 'Spectrum')
        grid3 = QtWidgets.QGridLayout(self)
        grid3.addWidget(self.tabs, 0, 0)
        grid3.setColumnStretch(0, 1)
        grid3.setRowStretch(0, 1)

    def getTabNames(self):
        return [self.tabs.tabText(i) for i in range(self.tabs.count())]

    def getCurrentTabName(self):
        return self.tabs.tabText(self.tabs.currentIndex())

    def addSpectrum(self):
        text = QtWidgets.QInputDialog.getItem(self, "Select spectrum to add", "Spectrum name:", self.father.workspaceNames, 0, False)
        if text[1]:
            self.subFitWindows.append(FittingWindow(self.father, self.father.workspaces[self.father.workspaceNames.index(text[0])], self, False))
            self.tabs.addTab(self.subFitWindows[-1], str(text[0]))
            self.tabs.setCurrentIndex(len(self.subFitWindows))

    def removeSpectrum(self, spec):
        num = self.subFitWindows.index(spec)
        self.tabs.removeTab(num + 1)
        del self.subFitWindows[num]

    def closeTab(self, num):
        if num > 0:
            self.tabs.removeTab(num)
            del self.subFitWindows[num - 1]
        else:
            self.mainFitWindow.paramframe.closeWindow()

    def fit(self):
        value = self.mainFitWindow.paramframe.getFitParams()
        if value is None:
            return
        else:
            xax, data1D, guess, args, out = value
        xax = [xax]
        out = [out]
        nameList = ['Spectrum']
        selectList = [slice(0, len(guess))]
        for i in range(len(self.subFitWindows)):
            xax_tmp, data1D_tmp, guess_tmp, args_tmp, out_tmp = self.subFitWindows[i].paramframe.getFitParams()
            out.append(out_tmp)
            xax.append(xax_tmp)
            nameList.append('bla')
            selectList.append(slice(len(guess), len(guess) + len(guess_tmp)))
            data1D = np.append(data1D, data1D_tmp)
            guess += guess_tmp
            new_args = ()
            for n in range(len(args)):
                new_args += (args[n] + args_tmp[n],)
            args = new_args  # tuples are immutable
        new_args = (nameList, selectList) + args
        allFitVal = self.mainFitWindow.paramframe.fit(xax, data1D, guess, new_args)
        if allFitVal is None:
            return
        allFitVal = allFitVal['x']
        fitVal = []
        for length in selectList:
            if allFitVal.ndim is 0:
                fitVal.append(np.array([allFitVal]))
            else:
                fitVal.append(allFitVal[length])
        args_out = []
        for n in range(len(args)):
            args_out.append([args[n][0]])
        self.mainFitWindow.paramframe.setResults(fitVal[0], args_out, out[0])
        for i in (range(len(self.subFitWindows))):
            args_out = []
            for n in range(len(args)):
                args_out.append([args[n][i + 1]])
            self.subFitWindows[i].paramframe.setResults(fitVal[i + 1], args_out, out[i + 1])

    def getNum(self, paramfitwindow):
        fitwindow = paramfitwindow.rootwindow
        if fitwindow is self.mainFitWindow:
            return 0
        else:
            return self.subFitWindows.index(fitwindow) + 1

    def getParams(self):
        params = [self.mainFitWindow.paramframe.getSimParams()]
        for window in self.subFitWindows:
            tmp_params = window.paramframe.getSimParams()
            params = np.append(params, [tmp_params], axis=0)
        if params[0] is None:
            return None
        return params
            
    def disp(self, *args, **kwargs):
        params = self.getParams()
        if params is None:
            return
        self.mainFitWindow.paramframe.disp(params, 0, *args, **kwargs)
        for i in range(len(self.subFitWindows)):
            self.subFitWindows[i].paramframe.disp(params, i + 1, *args, **kwargs)

    def get_masterData(self):
        return self.oldMainWindow.get_masterData()

    def get_current(self):
        return self.oldMainWindow.get_current()

    def kill(self):
        self.mainFitWindow.kill()

##############################################################################


class FitCopySettingsWindow(QtWidgets.QWidget):

    def __init__(self, parent, returnFunction, single=False):
        super(FitCopySettingsWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.returnFunction = returnFunction
        self.setWindowTitle("Settings")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        self.allTraces = QtWidgets.QCheckBox("Export all traces")
        if not single:
            grid.addWidget(self.allTraces, 0, 0)
        self.original = QtWidgets.QCheckBox("Include original")
        self.original.setChecked(True)
        grid.addWidget(self.original, 1, 0)
        self.subFits = QtWidgets.QCheckBox("Include subfits")
        self.subFits.setChecked(True)
        grid.addWidget(self.subFits, 2, 0)
        self.difference = QtWidgets.QCheckBox("Include difference")
        grid.addWidget(self.difference, 3, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok", self)
        okButton.clicked.connect(self.applyAndClose)
        okButton.setFocus()
        layout.addWidget(okButton, 2, 1)
        grid.setRowStretch(100, 1)
        self.show()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def closeEvent(self, *args):
        self.deleteLater()

    def applyAndClose(self, *args):
        self.deleteLater()
        self.returnFunction([self.allTraces.isChecked(), self.original.isChecked(), self.subFits.isChecked(), self.difference.isChecked()])

##################################################################################################


class ParamCopySettingsWindow(QtWidgets.QWidget):

    def __init__(self, parent, paramNames, returnFunction, single=False):
        super(ParamCopySettingsWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.returnFunction = returnFunction
        self.setWindowTitle("Parameters to export")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        self.allTraces = QtWidgets.QCheckBox("Export all traces")
        if not single:
            grid.addWidget(self.allTraces, 0, 0)
        self.exportList = []
        for i in range(len(paramNames)):
            self.exportList.append(QtWidgets.QCheckBox(paramNames[i]))
            grid.addWidget(self.exportList[-1], i + 1, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok", self)
        okButton.clicked.connect(self.applyAndClose)
        okButton.setFocus()
        layout.addWidget(okButton, 2, 1)
        grid.setRowStretch(100, 1)
        self.show()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def closeEvent(self, *args):
        self.deleteLater()

    def applyAndClose(self, *args):
        self.deleteLater()
        answers = []
        for checkbox in self.exportList:
            answers.append(checkbox.isChecked())
        self.returnFunction(self.allTraces.isChecked(), answers)

################################################################################


class FittingWindow(QtWidgets.QWidget):
    # Inherited by the fitting windows

    def __init__(self, father, oldMainWindow, tabWindow, isMain=True):
        super(FittingWindow, self).__init__(father)
        self.isMain = isMain
        self.father = father
        self.oldMainWindow = oldMainWindow
        self.tabWindow = tabWindow
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        grid = QtWidgets.QGridLayout(self)
        grid.addWidget(self.canvas, 0, 0)
        self.current = self.tabWindow.CURRENTWINDOW(self, self.fig, self.canvas, self.oldMainWindow.get_current())
        self.paramframe = self.tabWindow.PARAMFRAME(self.current, self, isMain=self.isMain)
        grid.addWidget(self.paramframe, 1, 0, 1, 2)
        grid.setColumnStretch(0, 1)
        grid.setRowStretch(0, 1)
        self.grid = grid
        self.fittingSideFrame = FittingSideFrame(self)
        self.grid.addWidget(self.fittingSideFrame, 0, 1)
        self.canvas.mpl_connect('button_press_event', self.buttonPress)
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)

    def updAllFrames(self, *args):
        pass
        
    def fit(self):
        self.tabWindow.fit()
        self.paramframe.togglePick()

    def sim(self, *args, **kwargs):
        self.tabWindow.disp(**kwargs)
        self.paramframe.togglePick()

    def getCurrentTabName(self, *args, **kwargs):
        return self.tabWindow.getCurrentTabName(*args, **kwargs)

    def getTabNames(self, *args, **kwargs):
        return self.tabWindow.getTabNames(*args, **kwargs)

    def getParams(self, *args, **kwargs):
        return self.tabWindow.getParams(*args, **kwargs)

    def getNum(self, *args, **kwargs):
        return self.tabWindow.getNum(*args, **kwargs)

    def addSpectrum(self):
        self.tabWindow.addSpectrum()

    def removeSpectrum(self):
        self.tabWindow.removeSpectrum(self)

    def createNewData(self, data, axes, params=False, fitAll=False):
        masterData = self.get_masterData()
        if fitAll:
            if params:
                self.father.dataFromFit(data,
                                        masterData.filePath,
                                        np.append([masterData.freq[axes], masterData.freq[axes]], np.delete(masterData.freq, axes)),
                                        np.append([masterData.sw[axes], masterData.sw[axes]], np.delete(masterData.sw, axes)),
                                        np.append([False, False], np.delete(masterData.spec, axes)),
                                        np.append([False, False], np.delete(masterData.wholeEcho, axes)),
                                        np.append([None, None], np.delete(masterData.ref, axes)),
                                        None,
                                        None)
            else:
                self.father.dataFromFit(data,
                                        masterData.filePath,
                                        np.append(masterData.freq[axes], masterData.freq),
                                        np.append(masterData.sw[axes], masterData.sw),
                                        np.append(False, masterData.spec),
                                        np.append(False, masterData.wholeEcho),
                                        np.append(None, masterData.ref),
                                        None,
                                        None)
        else:
            if params:
                self.father.dataFromFit(data,
                                        masterData.filePath,
                                        [masterData.freq[axes], masterData.freq[axes]],
                                        [masterData.sw[axes], masterData.sw[axes]],
                                        [False, False],
                                        [False, False],
                                        [None, None],
                                        [np.arange(data.shape[0]), np.arange(data.shape[1])],
                                        0)
            else:
                self.father.dataFromFit(data,
                                        masterData.filePath,
                                        np.append(masterData.freq[axes], self.current.data1D.freq),
                                        np.append(masterData.sw[axes], self.current.data1D.sw),
                                        np.append(False, self.current.data1D.spec),
                                        np.append(False, self.current.data1D.wholeEcho),
                                        np.append(None, self.current.data1D.ref),
                                        [np.arange(len(data))] + self.current.data1D.xaxArray,
                                        0)

    def rename(self, name):
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
        self.father.closeFitWindow(self.oldMainWindow)
        self.deleteLater()

##############################################################################


class FittingSideFrame(SideFrame):

    FITTING = True


#################################################################################


class FitPlotFrame(Current1D):

    MARKER = ''
    LINESTYLE = '-'
    FITNUM = 10  # Standard number of fits

    def __init__(self, rootwindow, fig, canvas, current):
        self.data = current.data
        tmp = np.array(current.data.shape(), dtype=int)
        tmp = np.delete(tmp, self.fixAxes(current.axes))
        self.fitDataList = np.full(tmp, None, dtype=object)
        self.fitPickNumList = np.zeros(tmp, dtype=int)
        self.rootwindow = rootwindow
        super(FitPlotFrame, self).__init__(rootwindow, fig, canvas, current.data, current)

    def getRedLocList(self):
        return tuple(np.delete(self.locList, self.axes))
        
    def setSlice(self, axes, locList):
        self.rootwindow.paramframe.checkInputs()
        self.pickWidth = False
        super(FitPlotFrame, self).setSlice(axes, locList)
        self.rootwindow.paramframe.checkFitParamList(self.getRedLocList())
        self.rootwindow.paramframe.dispParams()
        self.rootwindow.paramframe.togglePick()

    def getData1D(self):
        return np.real(self.getDataType(self.data1D.getHyperData(0)))

    def showFid(self):
        extraX = []
        extraY = []
        self.locList = np.array(self.locList, dtype=int)
        if self.fitDataList[self.getRedLocList()] is not None:
            tmp = self.fitDataList[self.getRedLocList()]
            extraX.append(tmp[0])
            extraY.append(tmp[1])
            for i in range(len(tmp[2])):
                extraX.append(tmp[2][i])
                extraY.append(tmp[3][i])
        super(FitPlotFrame, self).showFid(extraX=extraX, extraY=extraY)

#################################################################################


class AbstractParamFrame(QtWidgets.QWidget):

    TICKS = True  # Fitting parameters can be fixed by checkboxes

    def __init__(self, parent, rootwindow, isMain=True):
        super(AbstractParamFrame, self).__init__(rootwindow)
        self.parent = parent
        self.FITNUM = self.parent.FITNUM
        self.rootwindow = rootwindow
        self.isMain = isMain # display fitting buttons
        tmp = np.array(self.parent.data.shape(), dtype=int)
        tmp = np.delete(tmp, self.parent.axes)
        self.fitParamList = np.zeros(tmp, dtype=object)
        self.fitNumList = np.zeros(tmp, dtype=int)
        grid = QtWidgets.QGridLayout(self)
        self.setLayout(grid)
        self.axMult = self.parent.getAxMult(self.parent.spec(),
                                            self.parent.getAxType(),
                                            self.parent.getppm(),
                                            self.parent.freq(),
                                            self.parent.ref())
        if self.parent.spec() == 1:
            if self.parent.viewSettings["ppm"][-1]:
                self.axUnit = 'ppm'
            else:
                axUnits = ['Hz', 'kHz', 'MHz']
                self.axUnit = axUnits[self.parent.getAxType()]
        elif self.parent.spec() == 0:
            axUnits = ['s', 'ms', u"\u03bcs"]
            self.axUnit = axUnits[self.parent.getAxType()]
        self.frame1 = QtWidgets.QGridLayout()
        self.optframe = QtWidgets.QGridLayout()
        self.frame2 = QtWidgets.QGridLayout()
        self.frame3 = QtWidgets.QGridLayout()
        self.frame4 = QtWidgets.QGridLayout()
        grid.addLayout(self.frame1, 0, 0)
        grid.addLayout(self.optframe, 0, 1)
        grid.addLayout(self.frame2, 0, 2)
        grid.addLayout(self.frame3, 0, 3)
        grid.addLayout(self.frame4, 0, 4)
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        simButton = QtWidgets.QPushButton("Sim")
        simButton.clicked.connect(self.rootwindow.sim)
        self.frame1.addWidget(simButton, 0, 0)
        fitButton = QtWidgets.QPushButton("Fit")
        fitButton.clicked.connect(self.rootwindow.fit)
        self.frame1.addWidget(fitButton, 1, 0)
        self.stopButton = QtWidgets.QPushButton("Stop")
        self.stopButton.clicked.connect(self.stopMP)
        self.stopButton.setStyleSheet('background-color: green') 
        self.frame1.addWidget(self.stopButton, 1, 0)
        self.stopButton.hide()
        self.process1 = None
        self.queue = None
        fitAllButton = QtWidgets.QPushButton("Fit all")
        fitAllButton.clicked.connect(self.fitAll)
        self.frame1.addWidget(fitAllButton, 2, 0)
        self.stopAllButton = QtWidgets.QPushButton("Stop all")
        self.stopAllButton.clicked.connect(self.stopAll)
        self.stopAllButton.setStyleSheet('background-color: green') 
        self.frame1.addWidget(self.stopAllButton, 2, 0)
        self.stopAllButton.hide()
        copyParamsButton = QtWidgets.QPushButton("Copy parameters")
        copyParamsButton.clicked.connect(self.copyParams)
        self.frame1.addWidget(copyParamsButton, 3, 0)
        copyResultButton = QtWidgets.QPushButton("Result to workspace")
        copyResultButton.clicked.connect(self.resultToWorkspaceWindow)
        self.frame1.addWidget(copyResultButton, 4, 0)
        copyParamButton = QtWidgets.QPushButton("Param. to workspace")
        copyParamButton.clicked.connect(self.paramToWorkspaceWindow)
        self.frame1.addWidget(copyParamButton, 5, 0)
        addSpecButton = QtWidgets.QPushButton("Add spectrum")
        addSpecButton.clicked.connect(self.rootwindow.addSpectrum)
        self.frame1.addWidget(addSpecButton, 6, 0)
        if self.isMain:
            cancelButton = QtWidgets.QPushButton("&Cancel")
            cancelButton.clicked.connect(self.closeWindow)
        else:
            cancelButton = QtWidgets.QPushButton("&Delete")
            cancelButton.clicked.connect(self.rootwindow.removeSpectrum)
        prefButton = QtWidgets.QPushButton("Preferences")
        prefButton.clicked.connect(self.createPrefWindow)
        self.frame1.addWidget(prefButton, 0, 1)
        self.frame1.addWidget(cancelButton, 7, 0)
        self.frame1.setColumnStretch(10, 1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.checkFitParamList(self.getRedLocList())

    def getRedLocList(self):
        return self.parent.getRedLocList()

    def togglePick(self):
        # Dummy function for fitting routines which require peak picking
        pass
    
    def checkFitParamList(self, locList):
        locList = tuple(locList)
        if not self.fitParamList[locList]:
            self.fitParamList[locList] = self.defaultValues(0)

    def closeWindow(self, *args):
        self.stopMP()
        self.rootwindow.cancel()

    def copyParams(self):
        self.checkInputs()
        locList = self.getRedLocList()
        self.checkFitParamList(locList)
        for elem in np.nditer(self.fitParamList, flags=["refs_ok"], op_flags=['readwrite']):
            elem[...] = copy.deepcopy(self.fitParamList[locList])
        for elem in np.nditer(self.fitNumList, flags=["refs_ok"], op_flags=['readwrite']):
            elem[...] = self.fitNumList[locList]

    def dispParams(self):
        locList = self.getRedLocList()
        val = self.fitNumList[locList] + 1
        for name in self.SINGLENAMES:
            if isinstance(self.fitParamList[locList][name][0], (float, int)):
                self.entries[name][0].setText(('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % self.fitParamList[locList][name][0])
            if self.TICKS:
                self.ticks[name][0].setChecked(self.fitParamList[locList][name][1])
        self.setNumExp()
        for i in range(self.FITNUM):
            for name in self.MULTINAMES:
                if isinstance(self.fitParamList[locList][name][i][0], (float, int)):
                    self.entries[name][i].setText(('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % self.fitParamList[locList][name][i][0])
                if self.TICKS:
                    self.ticks[name][i].setChecked(self.fitParamList[locList][name][i][1])
                if i < val:
                    if self.TICKS:
                        self.ticks[name][i].show()
                    self.entries[name][i].show()
                else:
                    if self.TICKS:
                        self.ticks[name][i].hide()
                    self.entries[name][i].hide()

    def setNumExp(self):
        locList = self.getRedLocList()
        self.numExp.setCurrentIndex(self.fitNumList[locList])

    def changeNum(self, *args):
        val = self.numExp.currentIndex() + 1
        locList = self.getRedLocList()
        self.fitNumList[locList] = self.numExp.currentIndex()
        for i in range(self.FITNUM):
            for name in self.MULTINAMES:
                if i < val:
                    if self.TICKS:
                        self.ticks[name][i].show()
                    self.entries[name][i].show()
                else:
                    if self.TICKS:
                        self.ticks[name][i].hide()
                    self.entries[name][i].hide()

    def checkInputs(self):
        locList = self.getRedLocList()
        numExp = self.getNumExp()
        for name in self.SINGLENAMES:
            if self.TICKS:
                self.fitParamList[locList][name][1] = self.ticks[name][0].isChecked()
            inp = safeEval(self.entries[name][0].text())
            if inp is None:
                return False
            elif isinstance(inp, float):
                self.entries[name][0].setText(('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % inp)
            else:
                self.entries[name][0].setText(str(inp))
            self.fitParamList[locList][name][0] = inp
        for i in range(numExp):
            for name in self.MULTINAMES:
                if self.TICKS:
                    self.fitParamList[locList][name][i][1] = self.ticks[name][i].isChecked()
                inp = safeEval(self.entries[name][i].text())
                if inp is None:
                    return False
                elif isinstance(inp, float):
                    self.entries[name][i].setText(('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % inp)
                else:
                    self.entries[name][i].setText(str(inp))
                self.fitParamList[locList][name][i][0] = inp
        return True

    def getNumExp(self):
        return self.numExp.currentIndex() + 1

    def getExtraParams(self, out):
        return (out, [])

    def getFitParams(self):
        if not self.checkInputs():
            self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
            return
        struc = {}
        for name in (self.SINGLENAMES + self.MULTINAMES):
            struc[name] = []
        guess = []
        argu = []
        numExp = self.getNumExp()
        locList = self.getRedLocList()
        out = {}
        for name in self.SINGLENAMES:
            out[name] = [0.0]
        for name in self.MULTINAMES:
            out[name] = np.zeros(numExp)
        for name in self.SINGLENAMES:
            if isfloat(self.entries[name][0].text()):
                if not self.fitParamList[locList][name][1]:
                    guess.append(float(self.entries[name][0].text()))
                    struc[name].append((1, len(guess) - 1))
                else:
                    out[name][0] = float(self.entries[name][0].text())
                    argu.append(out[name][0])
                    struc[name].append((0, len(argu) - 1))
            else:
                struc[name].append((2, checkLinkTuple(safeEval(self.entries[name][0].text()))))
        for i in range(numExp):
            for name in self.MULTINAMES:
                if isfloat(self.entries[name][i].text()):
                    if not self.fitParamList[locList][name][i][1]:
                        guess.append(float(self.entries[name][i].text()))
                        struc[name].append((1, len(guess) - 1))
                    else:
                        out[name][i] = float(self.entries[name][i].text())
                        argu.append(out[name][i])
                        struc[name].append((0, len(argu) - 1))
                else:
                    struc[name].append((2, checkLinkTuple(safeEval(self.entries[name][i].text()))))
        out, extraArgu = self.getExtraParams(out)
        argu.append(extraArgu)
        args = ([numExp], [struc], [argu], [self.parent.data1D.sw], [self.axMult])
        return (self.parent.data1D.xaxArray, self.parent.getData1D(), guess, args, out)

    def fit(self, xax, data1D, guess, args):
        self.queue = multiprocessing.Queue()
        self.process1 = multiprocessing.Process(target=self.FITFUNC, args=(xax, data1D, guess, args, self.queue, self.rootwindow.tabWindow.MINMETHOD, self.rootwindow.tabWindow.NUMFEVAL))
        self.process1.start()
        self.running = True
        self.stopButton.show()
        while self.running:
            if not self.queue.empty():
                self.running = False
            QtWidgets.qApp.processEvents()
            time.sleep(0.1)
        if self.queue is None:
            return
        fitVal = self.queue.get(timeout=2)
        self.stopMP()
        if fitVal is None:
            self.rootwindow.father.dispMsg('Optimal parameters not found')
            return
        return fitVal

    def setResults(self, fitVal, args, out):
        locList = self.getRedLocList()
        numExp = args[0][0]
        struc = args[1][0]
        for name in self.SINGLENAMES:
            if struc[name][0][0] == 1:
                self.fitParamList[locList][name][0] = fitVal[struc[name][0][1]]
        for i in range(numExp):
            for name in self.MULTINAMES:
                if struc[name][i][0] == 1:
                    self.fitParamList[locList][name][i][0] = fitVal[struc[name][i][1]]
        self.dispParams()
        self.rootwindow.sim()

    def stopMP(self, *args):
        if self.queue is not None:
            self.process1.terminate()
            self.queue.close()
            self.queue.join_thread()
            self.process1.join()
        self.queue = None
        self.process1 = None
        self.running = False
        self.stopButton.hide()

    def stopAll(self, *args):
        self.runningAll = False
        self.stopMP()
        self.stopAllButton.hide()

    def fitAll(self, *args):
        self.runningAll = True
        self.stopAllButton.show()
        tmp = np.array(self.parent.data.shape())
        tmp[self.parent.axes] = 1
        tmp2 = ()
        for i in tmp:
            tmp2 += (np.arange(i),)
        grid = np.array([i.flatten() for i in np.meshgrid(*tmp2)]).T
        for i in grid:
            QtWidgets.qApp.processEvents()
            if self.runningAll is False:
                break
            self.parent.setSlice(self.parent.axes, i)
            self.rootwindow.fit()
            self.rootwindow.fittingSideFrame.upd()
        self.stopAllButton.hide()

    def getSimParams(self):
        if not self.checkInputs():
            self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
            return
        numExp = self.getNumExp()
        out = {}
        for name in self.SINGLENAMES:
            out[name] = [0.0]
        for name in self.MULTINAMES:
            out[name] = [0.0] * numExp
        for name in self.SINGLENAMES:
            inp = safeEval(self.entries[name][0].text())
            out[name][0] = inp
        for i in range(numExp):
            for name in self.MULTINAMES:
                inp = safeEval(self.entries[name][i].text())
                out[name][i] = inp
        out, tmp = self.getExtraParams(out)
        return out

    def paramToWorkspaceWindow(self):
        paramNameList = self.SINGLENAMES + self.MULTINAMES
        if self.parent.data.ndim() == 1:
            single = True
        else:
            single = False
        ParamCopySettingsWindow(self, paramNameList, lambda allTraces, settings, self=self: self.paramToWorkspace(allTraces, settings), single)

    def paramToWorkspace(self, allTraces, settings):
        if not self.checkInputs():
            self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
            return
        paramNameList = np.array(self.SINGLENAMES + self.MULTINAMES, dtype=object)
        locList = self.getRedLocList()
        if not np.any(settings):
            return
        names = paramNameList[settings]
        params = self.rootwindow.getParams()
        if allTraces:
            num = self.rootwindow.getNum(self)
            maxNum = np.max(self.fitNumList)+1
            tmp = np.array(self.parent.data.shape(), dtype=int)
            tmp = np.delete(tmp, self.parent.axes)
            data = np.zeros((sum(settings), maxNum) + tuple(tmp))
            tmp2 = ()
            for i in tmp:
                tmp2 += (np.arange(i),)
            grid = np.array([i.flatten() for i in np.meshgrid(*tmp2)]).T
            for i in grid:
                self.checkFitParamList(tuple(i))
                for j in range(len(names)):
                    if names[j] in self.SINGLENAMES:
                        inp = self.fitParamList[tuple(i)][names[j]][0]
                        if isinstance(inp, tuple):
                            inp = checkLinkTuple(inp)
                            if inp[4] is num:
                                inp = inp[2] * self.fitParamList[tuple(i)][inp[0]][inp[1]][0] + inp[3]
                            else:
                                inp = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                        data[(j,) + (slice(None),) + tuple(i)].fill(inp)
                    else:
                        tmpInp = self.fitParamList[tuple(i)][names[j]].T[0][:(self.fitNumList[tuple(i)] + 1)]
                        for n in range(len(tmpInp)):
                            inp = tmpInp[n]
                            if isinstance(inp, tuple):
                                inp = checkLinkTuple(inp)
                                if inp[4] is num:
                                    inp = inp[2] * self.fitParamList[tuple(i)][inp[0]][inp[1]][0] + inp[3]
                                else:
                                    inp = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                            data[(j,) + (slice(None),) + tuple(i)][n] = inp
            self.rootwindow.createNewData(data, self.parent.axes[-1], True, True)
        else:
            data = np.zeros((sum(settings), self.fitNumList[locList] + 1))
            for i in range(len(names)):
                if names[i] in self.SINGLENAMES:
                    inp = self.fitParamList[locList][names[i]][0]
                    if isinstance(inp, tuple):
                        inp = checkLinkTuple(inp)
                        inp = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                    data[i].fill(inp)
                else:
                    tmpInp = self.fitParamList[locList][names[i]].T[0][:(self.fitNumList[locList] + 1)]
                    for j in range(len(tmpInp)):
                        inp = tmpInp[j]
                        if isinstance(inp, tuple):
                            inp = checkLinkTuple(inp)
                            inp = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                        data[i][j] = inp
            self.rootwindow.createNewData(data, self.parent.axes[-1], True)

    def resultToWorkspaceWindow(self):
        if self.parent.data.ndim() == 1:
            single = True
        else:
            single = False
        FitCopySettingsWindow(self, lambda settings, self=self: self.resultToWorkspace(settings), single)

    def resultToWorkspace(self, settings):
        if settings is None:
            return
        if settings[0]:
            oldLocList = self.parent.locList
            maxNum = np.max(self.fitNumList) + 1
            extraLength = 1
            if settings[1]:
                extraLength += 1
            if settings[2]:
                extraLength += maxNum
            if settings[3]:
                extraLength += 1
            data = np.zeros((extraLength,) + self.parent.data.shape())
            tmp = np.array(self.parent.data.shape(), dtype=int)
            tmp[self.parent.axes] = 1
            tmp2 = ()
            for i in tmp:
                tmp2 += (np.arange(i),)
            grid = np.array([i.flatten() for i in np.meshgrid(*tmp2)]).T
            for i in grid:
                self.parent.setSlice(self.parent.axes, i)
                j = np.delete(i, self.parent.axes)
                data[(slice(None),) + tuple(j[:self.parent.axes[-1]]) + (slice(None),) + tuple(j[self.parent.axes[-1]:])] = self.prepareResultToWorkspace(settings, maxNum)
            self.parent.setSlice(self.parent.axes, oldLocList)
            self.rootwindow.createNewData(data, self.parent.axes[-1], False, True)
        else:
            data = self.prepareResultToWorkspace(settings)
            self.rootwindow.createNewData(data, self.parent.axes[-1], False)

    def prepareResultToWorkspace(self, settings, minLength=1):
        self.calculateResultsToWorkspace(True)
        locList = self.getRedLocList()
        fitData = self.parent.fitDataList[locList]
        if fitData is None:
            fitData = [np.zeros(len(self.parent.getData1D())), np.zeros(len(self.parent.getData1D())), np.zeros(len(self.parent.getData1D())), np.array([np.zeros(len(self.parent.getData1D()))] * minLength)]
        outCurvePart = []
        if settings[1]:
            outCurvePart.append(self.parent.getData1D())
        if settings[2]:
            for i in fitData[3]:
                outCurvePart.append(i)
            if len(fitData[3]) < minLength:
                for i in range(minLength - len(fitData[3])):
                    outCurvePart.append(np.zeros(len(self.parent.getData1D())))
        if settings[3]:
            outCurvePart.append(self.parent.getData1D() - fitData[1])
        outCurvePart.append(fitData[1])
        self.calculateResultsToWorkspace(False)
        return np.array(outCurvePart)

    def calculateResultsToWorkspace(self, *args):
        # Some fitting methods need to recalculate the curves before exporting
        pass

    def createPrefWindow(self, *args):
        PrefWindow(self.rootwindow.tabWindow)

##############################################################################


class PrefWindow(QtWidgets.QWidget):

    METHODLIST = ['Powell', 'Nelder-Mead']

    def __init__(self, parent):
        super(PrefWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Preferences")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Min. method:"), 0, 0)
        self.minmethodBox = QtWidgets.QComboBox(self)
        self.minmethodBox.addItems(self.METHODLIST)
        self.minmethodBox.setCurrentIndex(self.METHODLIST.index(self.father.MINMETHOD))
        grid.addWidget(self.minmethodBox, 0, 1)
        grid.addWidget(wc.QLabel("Significant digits:"), 1, 0)
        self.precisBox = QtWidgets.QSpinBox(self)
        self.precisBox.setValue(self.father.PRECIS)
        grid.addWidget(self.precisBox, 1, 1)
        grid.addWidget(wc.QLabel("# evaluations:"), 2, 0)
        self.numFevalBox = QtWidgets.QSpinBox(self)
        self.numFevalBox.setMaximum(100000)
        self.numFevalBox.setMinimum(1)
        self.numFevalBox.setValue(self.father.NUMFEVAL)
        grid.addWidget(self.numFevalBox, 2, 1)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 4, 0)
        okButton = QtWidgets.QPushButton("&Ok", self)
        okButton.clicked.connect(self.applyAndClose)
        okButton.setFocus()
        layout.addWidget(okButton, 4, 1)
        grid.setRowStretch(100, 1)
        self.show()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def closeEvent(self, *args):
        self.deleteLater()

    def applyAndClose(self, *args):
        self.father.PRECIS = self.precisBox.value()
        self.father.MINMETHOD = self.METHODLIST[self.minmethodBox.currentIndex()]
        self.father.NUMFEVAL = self.numFevalBox.value()
        self.closeEvent()

##############################################################################


class RelaxWindow(TabFittingWindow):

    def __init__(self, father, oldMainWindow):
        self.CURRENTWINDOW = RelaxFrame
        self.PARAMFRAME = RelaxParamFrame
        super(RelaxWindow, self).__init__(father, oldMainWindow)

#################################################################################


class RelaxFrame(FitPlotFrame):

    MARKER = 'o'
    LINESTYLE = 'none'
    FITNUM = 4  # Maximum number of fits

#################################################################################


class RelaxParamFrame(AbstractParamFrame):

    def __init__(self, parent, rootwindow, isMain=True):
        self.SINGLENAMES = ['amp', 'const']
        self.MULTINAMES = ['coeff', 't']
        self.PARAMTEXT = {'amp': 'Amplitude', 'const': 'Constant', 'coeff': 'Coefficient', 't': 'Relaxation time'}
        self.FITFUNC = relaxationmpFit
        super(RelaxParamFrame, self).__init__(parent, rootwindow, isMain)
        locList = self.getRedLocList()
        self.ticks = {'amp': [], 'const': [], 'coeff': [], 't': []}
        self.entries = {'amp': [], 'const': [], 'coeff': [], 't': []}
        self.frame2.addWidget(wc.QLabel("Amplitude:"), 0, 0, 1, 2)
        self.ticks['amp'].append(QtWidgets.QCheckBox(''))
        self.frame2.addWidget(self.ticks['amp'][-1], 1, 0)
        self.entries['amp'].append(wc.FitQLineEdit(self, 'amp', ('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % self.fitParamList[locList]['amp'][0]))
        self.frame2.addWidget(self.entries['amp'][-1], 1, 1)
        self.frame2.addWidget(wc.QLabel("Constant:"), 2, 0, 1, 2)
        self.ticks['const'].append(QtWidgets.QCheckBox(''))
        self.frame2.addWidget(self.ticks['const'][-1], 3, 0)
        self.entries['const'].append(wc.FitQLineEdit(self, 'const', ('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % self.fitParamList[locList]['const'][0]))
        self.frame2.addWidget(self.entries['const'][-1], 3, 1)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtWidgets.QComboBox()
        self.numExp.addItems(['1', '2', '3', '4'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        self.frame3.addWidget(wc.QLabel("Coefficient:"), 1, 0, 1, 2)
        self.frame3.addWidget(wc.QLabel("T [s]:"), 1, 2, 1, 2)
        self.frame3.setColumnStretch(10, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.xlog = QtWidgets.QCheckBox('x-log')
        self.xlog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.xlog, 0, 0, QtCore.Qt.AlignTop)
        self.ylog = QtWidgets.QCheckBox('y-log')
        self.ylog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.ylog, 1, 0, QtCore.Qt.AlignTop)
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        for i in range(self.FITNUM):
            for j in range(len(self.MULTINAMES)):
                self.ticks[self.MULTINAMES[j]].append(QtWidgets.QCheckBox(''))
                self.ticks[self.MULTINAMES[j]][i].setChecked(self.fitParamList[locList][self.MULTINAMES[j]][i][1])
                self.frame3.addWidget(self.ticks[self.MULTINAMES[j]][i], i + 2, 2 * j)
                self.entries[self.MULTINAMES[j]].append(wc.FitQLineEdit(self, self.MULTINAMES[j], ('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % self.fitParamList[locList][self.MULTINAMES[j]][i][0]))
                self.frame3.addWidget(self.entries[self.MULTINAMES[j]][i], i + 2, 2 * j + 1)
                if i > 0:
                    self.ticks[self.MULTINAMES[j]][i].hide()
                    self.entries[self.MULTINAMES[j]][i].hide()

    def defaultValues(self, inp):
        if not inp:
            return {'amp': [np.max(self.parent.getData1D()), False],
                    'const': [1.0, False],
                    'coeff': np.repeat([np.array([-1.0, False], dtype=object)], self.FITNUM, axis=0),
                    't': np.repeat([np.array([1.0, False], dtype=object)], self.FITNUM, axis=0)}
        else:
            return inp

    def setLog(self, *args):
        self.parent.setLog(self.xlog.isChecked(), self.ylog.isChecked())
        self.rootwindow.sim()

    def disp(self, params, num, display=True, prepExport=False):
        out = params[num]
        for name in ['amp', 'const']:
            inp = out[name][0]
            if isinstance(inp, tuple):
                inp = checkLinkTuple(inp)
                out[name][0] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
        numExp = len(out['coeff'])
        for i in range(numExp):
            for name in ['coeff', 't']:
                inp = out[name][i]
                if isinstance(inp, tuple):
                    inp = checkLinkTuple(inp)
                    out[name][i] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                if not np.isfinite(out[name][i]):
                    self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
                    return
        tmpx = self.parent.xax()
        if prepExport:
            x = tmpx
            outCurve = out['const'][0] * np.ones(len(x))
        else:
            numCurve = 100  # number of points in output curve
            outCurve = out['const'][0] * np.ones(numCurve)
            if self.xlog.isChecked():
                x = np.logspace(np.log(min(tmpx)), np.log(max(tmpx)), numCurve)
            else:
                x = np.linspace(min(tmpx), max(tmpx), numCurve)
        for i in range(len(out['coeff'])):
            outCurve += out['coeff'][i] * np.exp(-x / out['t'][i])
        locList = self.getRedLocList()
        self.parent.fitDataList[locList] = [x, out['amp'][0] * outCurve, [], []]
        if display:
            self.parent.showFid()

    def calculateResultsToWorkspace(self, prepExport):
        self.rootwindow.sim(display=False, prepExport=prepExport)

#############################################################################


def relaxationmpFit(xax, data1D, guess, args, queue, minmethod, numfeval):
    try:
        fitVal = scipy.optimize.minimize(lambda *param: np.sum((data1D - relaxationfitFunc(param, xax, args))**2), guess, method=minmethod, options = {'maxfev': numfeval})
    except Exception:
        fitVal = None
    queue.put(fitVal)


def relaxationfitFunc(params, allX, args):
    params = params[0]
    specName = args[0]
    specSlices = args[1]
    allParam = []
    for length in specSlices:
        allParam.append(params[length])
    allStruc = args[3]
    allArgu = args[4]
    fullTestFunc = []
    for n in range(len(allX)):
        x = allX[n][-1]
        testFunc = np.zeros(len(x))
        param = allParam[n]
        numExp = args[2][n]
        struc = args[3][n]
        argu = args[4][n]
        sw = args[5][n][-1]
        axMult = args[6][n]
        parameters = {}
        for name in ['amp', 'const']:
            if struc[name][0][0] == 1:
                parameters[name] = param[struc[name][0][1]]
            elif struc[name][0][0] == 0:
                parameters[name] = argu[struc[name][0][1]]
            else:
                altStruc = struc[name][0][1]
                if struc[altStruc[0]][altStruc[1]][0] == 1:
                    parameters[name] = altStruc[2] * allParam[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                elif struc[altStruc[0]][altStruc[1]][0] == 0:
                    parameters[name] = altStruc[2] * allArgu[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
        for i in range(numExp):
            for name in ['coeff', 't']:
                if struc[name][i][0] == 1:
                    parameters[name] = param[struc[name][i][1]]
                elif struc[name][i][0] == 0:
                    parameters[name] = argu[struc[name][i][1]]
                else:
                    altStruc = struc[name][i][1]
                    strucTarget = allStruc[altStruc[4]]
                    if strucTarget[altStruc[0]][altStruc[1]][0] == 1:
                        parameters[name] = altStruc[2] * allParam[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                    elif strucTarget[altStruc[0]][altStruc[1]][0] == 0:
                        parameters[name] = altStruc[2] * allArgu[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
            testFunc += parameters['coeff'] * np.exp(-1.0 * x / parameters['t'])
        fullTestFunc = np.append(fullTestFunc, parameters['amp'] * (parameters['const'] + testFunc))
    return fullTestFunc

##############################################################################


class DiffusionWindow(TabFittingWindow):

    def __init__(self, father, oldMainWindow):
        self.CURRENTWINDOW = RelaxFrame
        self.PARAMFRAME = DiffusionParamFrame
        super(DiffusionWindow, self).__init__(father, oldMainWindow)

#################################################################################


class DiffusionParamFrame(AbstractParamFrame):

    def __init__(self, parent, rootwindow, isMain=True):
        self.SINGLENAMES = ['amp', 'const']
        self.MULTINAMES = ['coeff', 'd']
        self.PARAMTEXT = {'amp': 'Amplitude', 'const': 'Constant', 'coeff': 'Coefficient', 'd': 'Diffusion constant'}
        self.FITFUNC = diffusionmpFit
        super(DiffusionParamFrame, self).__init__(parent, rootwindow, isMain)
        locList = self.getRedLocList()
        self.ticks = {'amp': [], 'const': [], 'coeff': [], 'd': []}
        self.entries = {'amp': [], 'const': [], 'coeff': [], 'd': []}
        self.frame2.addWidget(wc.QLabel(u"\u03b3 [MHz/T]:"), 0, 0)
        self.gammaEntry = wc.QLineEdit("42.576")
        self.frame2.addWidget(self.gammaEntry, 1, 0)
        self.frame2.addWidget(wc.QLabel(u"\u03b4 [s]:"), 2, 0)
        self.deltaEntry = wc.QLineEdit("1.0")
        self.frame2.addWidget(self.deltaEntry, 3, 0)
        self.frame2.addWidget(wc.QLabel(u"\u0394 [s]:"), 4, 0)
        self.triangleEntry = wc.QLineEdit("1.0")
        self.frame2.addWidget(self.triangleEntry, 5, 0)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.frame3.addWidget(wc.QLabel("Amplitude:"), 0, 0, 1, 2)
        self.ticks['amp'].append(QtWidgets.QCheckBox(''))
        self.frame3.addWidget(self.ticks['amp'][-1], 1, 0)
        self.entries['amp'].append(wc.FitQLineEdit(self, 'amp', ('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % np.max(self.parent.getData1D())))
        self.frame3.addWidget(self.entries['amp'][-1], 1, 1)
        self.frame3.addWidget(wc.QLabel("Constant:"), 2, 0, 1, 2)
        self.ticks['const'].append(QtWidgets.QCheckBox(''))
        self.ticks['const'][-1].setChecked(True)
        self.frame3.addWidget(self.ticks['const'][-1], 3, 0)
        self.entries['const'].append(wc.FitQLineEdit(self, 'const', "0.0"))
        self.frame3.addWidget(self.entries['const'][-1], 3, 1)
        self.frame3.setColumnStretch(10, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtWidgets.QComboBox()
        self.numExp.addItems(['1', '2', '3', '4'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame4.addWidget(self.numExp, 0, 0, 1, 2)
        self.frame4.addWidget(wc.QLabel("Coefficient:"), 1, 0, 1, 2)
        self.frame4.addWidget(wc.QLabel("D [m^2/s]:"), 1, 2, 1, 2)
        self.frame4.setColumnStretch(20, 1)
        self.frame4.setAlignment(QtCore.Qt.AlignTop)
        self.xlog = QtWidgets.QCheckBox('x-log')
        self.xlog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.xlog, 0, 0)
        self.ylog = QtWidgets.QCheckBox('y-log')
        self.ylog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.ylog, 1, 0)
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.coeffEntries = []
        self.coeffTicks = []
        self.dEntries = []
        self.dTicks = []
        for i in range(self.FITNUM):
            for j in range(len(self.MULTINAMES)):
                self.ticks[self.MULTINAMES[j]].append(QtWidgets.QCheckBox(''))
                self.ticks[self.MULTINAMES[j]][i].setChecked(self.fitParamList[locList][self.MULTINAMES[j]][i][1])
                self.frame4.addWidget(self.ticks[self.MULTINAMES[j]][i], i + 2, 2 * j)
                self.entries[self.MULTINAMES[j]].append(wc.FitQLineEdit(self, self.MULTINAMES[j], ('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % self.fitParamList[locList][self.MULTINAMES[j]][i][0]))
                self.frame4.addWidget(self.entries[self.MULTINAMES[j]][i], i + 2, 2 * j + 1)
                if i > 0:
                    self.ticks[self.MULTINAMES[j]][i].hide()
                    self.entries[self.MULTINAMES[j]][i].hide()

    def defaultValues(self, inp):
        if not inp:
            return {'amp': [np.max(self.parent.getData1D()), False],
                    'const': [0.0, False],
                    'coeff': np.repeat([np.array([1.0, False], dtype=object)], self.FITNUM, axis=0),
                    'd': np.repeat([np.array([1.0e-9, False], dtype=object)], self.FITNUM, axis=0)}
        else:
            return inp

    def setLog(self, *args):
        self.parent.setLog(self.xlog.isChecked(), self.ylog.isChecked())
        self.rootwindow.sim()

    def getExtraParams(self, out):
        out['gamma'] = [safeEval(self.gammaEntry.text())]
        out['delta'] = [safeEval(self.deltaEntry.text())]
        out['triangle'] = [safeEval(self.triangleEntry.text())]
        return (out, [out['gamma'][-1], out['delta'][-1], out['triangle'][-1]])

    def disp(self, params, num, display=True, prepExport=False):
        out = params[num]
        for name in ['amp', 'const', 'gamma', 'delta', 'triangle']:
            inp = out[name][0]
            if isinstance(inp, tuple):
                inp = checkLinkTuple(inp)
                out[name][0] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
        numExp = len(out['coeff'])
        for i in range(numExp):
            for name in ['coeff', 'd']:
                inp = out[name][i]
                if isinstance(inp, tuple):
                    inp = checkLinkTuple(inp)
                    out[name][i] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                if not np.isfinite(out[name][i]):
                    self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
                    return
        tmpx = self.parent.xax()
        if prepExport:
            x = tmpx
            outCurve = out['const'][0] * np.ones(len(x))
        else:
            numCurve = 100  # number of points in output curve
            outCurve = out['const'][0] * np.ones(numCurve)
            if self.xlog.isChecked():
                x = np.logspace(np.log(min(tmpx)), np.log(max(tmpx)), numCurve)
            else:
                x = np.linspace(min(tmpx), max(tmpx), numCurve)
        for i in range(len(out['coeff'])):
            outCurve += out['coeff'][i] * np.exp(-(out['gamma'][0] * out['delta'][0] * x)**2 * out['d'][i] * (out['triangle'][0] - out['delta'][0] / 3.0))
        locList = self.getRedLocList()
        self.parent.fitDataList[locList] = [x, out['amp'][0] * outCurve, [], []]
        if display:
            self.parent.showFid()

    def calculateResultsToWorkspace(self, prepExport):
        self.rootwindow.sim(display=False, prepExport=prepExport)

##############################################################################


def diffusionmpFit(xax, data1D, guess, args, queue, minmethod, numfeval):
    try:
        fitVal = scipy.optimize.minimize(lambda *param: np.sum((data1D - diffusionfitFunc(param, xax, args))**2), guess, method=minmethod, options = {'maxfev': numfeval})
    except Exception:
        fitVal = None
    queue.put(fitVal)


def diffusionfitFunc(params, allX, args):
    params = params[0]
    specName = args[0]
    specSlices = args[1]
    allParam = []
    for length in specSlices:
        allParam.append(params[length])
    allStruc = args[3]
    allArgu = args[4]
    fullTestFunc = []
    for n in range(len(allX)):
        x = allX[n][-1]
        testFunc = np.zeros(len(x))
        param = allParam[n]
        numExp = args[2][n]
        struc = args[3][n]
        argu = args[4][n]
        sw = args[5][n][-1]
        axMult = args[6][n]
        parameters = {}
        parameters['gamma'] = argu[-1][0]
        parameters['delta'] = argu[-1][1]
        parameters['triangle'] = argu[-1][2]
        for name in ['amp', 'const']:
            if struc[name][0][0] == 1:
                parameters[name] = param[struc[name][0][1]]
            elif struc[name][0][0] == 0:
                parameters[name] = argu[struc[name][0][1]]
            else:
                altStruc = struc[name][0][1]
                if struc[altStruc[0]][altStruc[1]][0] == 1:
                    parameters[name] = altStruc[2] * allParam[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                elif struc[altStruc[0]][altStruc[1]][0] == 0:
                    parameters[name] = altStruc[2] * allArgu[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
        for i in range(numExp):
            for name in ['coeff', 'd']:
                if struc[name][i][0] == 1:
                    parameters[name] = param[struc[name][i][1]]
                elif struc[name][i][0] == 0:
                    parameters[name] = argu[struc[name][i][1]]
                else:
                    altStruc = struc[name][i][1]
                    strucTarget = allStruc[altStruc[4]]
                    if strucTarget[altStruc[0]][altStruc[1]][0] == 1:
                        parameters[name] = altStruc[2] * allParam[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                    elif strucTarget[altStruc[0]][altStruc[1]][0] == 0:
                        parameters[name] = altStruc[2] * allArgu[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
            testFunc += parameters['coeff'] * np.exp(-(parameters['gamma'] * parameters['delta'] * x)**2 * parameters['d'] * (parameters['triangle'] - parameters['delta'] / 3.0))
        fullTestFunc = np.append(fullTestFunc, parameters['amp'] * (parameters['const'] + testFunc))
    return fullTestFunc

##############################################################################


class PeakDeconvWindow(TabFittingWindow):
    
    def __init__(self, father, oldMainWindow):
        self.CURRENTWINDOW = PeakDeconvFrame
        self.PARAMFRAME = PeakDeconvParamFrame
        super(PeakDeconvWindow, self).__init__(father, oldMainWindow)

#################################################################################


class PeakDeconvFrame(FitPlotFrame):
    
    FITNUM = 10

    def togglePick(self, var):
        self.peakPickReset()
        if var == 1 and self.fitPickNumList[self.getRedLocList()] < self.FITNUM:
            self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos)
            self.peakPick = True
        else:
            self.peakPickFunc = None
            self.peakPick = False

    def pickDeconv(self, pos):
        locList = self.getRedLocList()
        pickNum = self.fitPickNumList[locList]
        if self.pickWidth:
            axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
            width = (2 * abs(float(self.rootwindow.paramframe.entries['pos'][pickNum].text()) - pos[1])) / self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
            self.rootwindow.paramframe.entries['amp'][pickNum].setText(('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % (float(self.rootwindow.paramframe.entries['amp'][pickNum].text()) * width))
            self.rootwindow.paramframe.entries['lor'][pickNum].setText(('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % abs(width))
            self.fitPickNumList[locList] += 1
            self.pickWidth = False
            self.rootwindow.sim()
        else:
            self.rootwindow.paramframe.entries['pos'][pickNum].setText(('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % pos[1])
            left = pos[0] - self.FITNUM
            if left < 0:
                left = 0
            right = pos[0] + self.FITNUM
            if right >= self.len():
                right = self.len() - 1
            self.rootwindow.paramframe.entries['amp'][pickNum].setText(('%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g') % (pos[2] * np.pi * 0.5))
            if pickNum < self.FITNUM:
                self.rootwindow.paramframe.numExp.setCurrentIndex(pickNum)
                self.rootwindow.paramframe.changeNum()
            self.pickWidth = True
        if pickNum < self.FITNUM:
            self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos)
            self.peakPick = True

#################################################################################


class PeakDeconvParamFrame(AbstractParamFrame):

    def __init__(self, parent, rootwindow, isMain=True):
        self.SINGLENAMES = ['bgrnd']
        self.MULTINAMES = ['pos', 'amp', 'lor', 'gauss']
        self.PARAMTEXT = {'bgrnd': 'Background', 'pos': 'Position', 'amp': 'Integral', 'lor': 'Lorentz', 'gauss': 'Gauss'}
        self.FITFUNC = peakDeconvmpFit
        # Get full integral
        self.fullInt = np.sum(parent.getData1D()) * parent.sw() / float(len(parent.getData1D()))
        super(PeakDeconvParamFrame, self).__init__(parent, rootwindow, isMain)
        resetButton = QtWidgets.QPushButton("Reset")
        resetButton.clicked.connect(self.reset)
        self.frame1.addWidget(resetButton, 1, 1)
        self.pickTick = QtWidgets.QCheckBox("Pick")
        self.pickTick.stateChanged.connect(self.togglePick)
        self.frame1.addWidget(self.pickTick, 2, 1)
        self.ticks = {'bgrnd': [], 'pos': [], 'amp': [], 'lor': [], 'gauss': []}
        self.entries = {'bgrnd': [], 'pos': [], 'amp': [], 'lor': [], 'gauss': [], 'method': []}
        self.frame2.addWidget(wc.QLabel("Bgrnd:"), 0, 0, 1, 2)
        self.ticks['bgrnd'].append(QtWidgets.QCheckBox(''))
        self.frame2.addWidget(self.ticks['bgrnd'][0], 1, 0)
        self.entries['bgrnd'].append(wc.FitQLineEdit(self, 'bgrnd', "0.0"))
        self.frame2.addWidget(self.entries['bgrnd'][0], 1, 1)
        self.frame2.addWidget(wc.QLabel("Method:"), 4, 0, 1, 2)
        self.entries['method'].append(QtWidgets.QComboBox())
        self.entries['method'][0].addItems(['Exact', 'Approx'])
        self.frame2.addWidget(self.entries['method'][0], 5, 1)
        self.frame2.setColumnStretch(self.FITNUM, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtWidgets.QComboBox()
        self.numExp.addItems([str(x + 1) for x in range(self.FITNUM)])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        self.frame3.addWidget(wc.QLabel("Position [" + self.axUnit + "]:"), 1, 0, 1, 2)
        self.frame3.addWidget(wc.QLabel("Integral:"), 1, 2, 1, 2)
        self.frame3.addWidget(wc.QLabel("Lorentz [Hz]:"), 1, 4, 1, 2)
        self.frame3.addWidget(wc.QLabel("Gauss [Hz]:"), 1, 6, 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        for i in range(self.FITNUM):
            for j in range(len(self.MULTINAMES)):
                self.ticks[self.MULTINAMES[j]].append(QtWidgets.QCheckBox(''))
                self.frame3.addWidget(self.ticks[self.MULTINAMES[j]][i], i + 2, 2 * j)
                self.entries[self.MULTINAMES[j]].append(wc.FitQLineEdit(self, self.MULTINAMES[j]))
                self.frame3.addWidget(self.entries[self.MULTINAMES[j]][i], i + 2, 2 * j + 1)
        self.reset()

    def defaultValues(self, inp):
        if not inp:
            return {'bgrnd': [0.0, True],
                    'pos': np.repeat([np.array([0.0, False], dtype=object)], self.FITNUM, axis=0),
                    'amp': np.repeat([np.array([self.fullInt, False], dtype=object)], self.FITNUM, axis=0),
                    'lor': np.repeat([np.array([1.0, False], dtype=object)], self.FITNUM, axis=0),
                    'gauss': np.repeat([np.array([0.0, True], dtype=object)], self.FITNUM, axis=0)}
        else:
            return inp

    def reset(self):
        locList = self.getRedLocList()
        self.fitNumList[locList] = 0
        for name in ['bgrnd']:
            self.fitParamList[locList][name] = [0.0, True]
        self.pickTick.setChecked(True)
        defaults = {'pos': [0.0, False], 'amp': [self.fullInt, False], 'lor': [1.0, False], 'gauss': [0.0, True]}
        for i in range(self.FITNUM):
            for name in defaults.keys():
                self.fitParamList[locList][name][i] = defaults[name]
        self.togglePick()
        self.parent.pickWidth = False
        self.parent.fitPickNumList[locList] = 0
        self.dispParams()

    def togglePick(self):
        self.parent.togglePick(self.pickTick.isChecked())

    def getExtraParams(self, out):
        out['method'] = [self.entries['method'][0].currentIndex()]
        return (out, [out['method'][-1]])

    def disp(self, params, num):
        out = params[num]
        for name in self.SINGLENAMES:
            inp = out[name][0]
            if isinstance(inp, tuple):
                inp = checkLinkTuple(inp)
                out[name][0] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
        numExp = len(out['pos'])
        for i in range(numExp):
            for name in self.MULTINAMES:
                inp = out[name][i]
                if isinstance(inp, tuple):
                    inp = checkLinkTuple(inp)
                    out[name][i] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                if not np.isfinite(out[name][i]):
                    self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
                    return
        tmpx = self.parent.xax()
        bgrnd = out['bgrnd'][0]
        outCurve = bgrnd * np.ones(len(tmpx))
        outCurvePart = []
        x = []
        for i in range(len(out['amp'])):
            x.append(tmpx)
            pos = out['pos'][i] / self.axMult
            y = simFunc.voigtLine(tmpx, pos, out['lor'][i], out['gauss'][i], out['amp'][i], out['method'][0])
            outCurvePart.append(bgrnd + y)
            outCurve += y
        locList = self.getRedLocList()
        self.parent.fitDataList[locList] = [tmpx, outCurve, x, outCurvePart]
        self.parent.showFid()

##############################################################################


def peakDeconvmpFit(xax, data1D, guess, args, queue, minmethod, numfeval):
    try:
        fitVal = scipy.optimize.minimize(lambda *param: np.sum((data1D - peakDeconvfitFunc(param, xax, args))**2), guess, method=minmethod, options = {'maxfev': numfeval})
    except Exception:
        fitVal = None
    queue.put(fitVal)


def peakDeconvfitFunc(params, allX, args):
    params = params[0]
    specName = args[0]
    specSlices = args[1]
    allParam = []
    for length in specSlices:
        allParam.append(params[length])
    allStruc = args[3]
    allArgu = args[4]
    fullTestFunc = []
    for n in range(len(allX)):
        x = allX[n][-1]
        testFunc = np.zeros(len(x))
        param = allParam[n]
        numExp = args[2][n]
        struc = args[3][n]
        argu = args[4][n]
        sw = args[5][n][-1]
        axMult = args[6][n]
        parameters = {}
        parameters['method'] = argu[-1][0]
        for name in ['bgrnd']:
            if struc[name][0][0] == 1:
                parameters[name] = param[struc[name][0][1]]
            elif struc[name][0][0] == 0:
                parameters[name] = argu[struc[name][0][1]]
            else:
                altStruc = struc[name][0][1]
                if struc[altStruc[0]][altStruc[1]][0] == 1:
                    parameters[name] = altStruc[2] * allParam[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                elif struc[altStruc[0]][altStruc[1]][0] == 0:
                    parameters[name] = altStruc[2] * allArgu[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
        for i in range(numExp):
            for name in ['pos', 'amp', 'lor', 'gauss']:
                if struc[name][i][0] == 1:
                    parameters[name] = param[struc[name][i][1]]
                elif struc[name][i][0] == 0:
                    parameters[name] = argu[struc[name][i][1]]
                else:
                    altStruc = struc[name][i][1]
                    strucTarget = allStruc[altStruc[4]]
                    if strucTarget[altStruc[0]][altStruc[1]][0] == 1:
                        parameters[name] = altStruc[2] * allParam[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                    elif strucTarget[altStruc[0]][altStruc[1]][0] == 0:
                        parameters[name] = altStruc[2] * allArgu[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
            pos = parameters['pos'] / axMult
            testFunc += simFunc.voigtLine(x, pos, parameters['lor'], parameters['gauss'], parameters['amp'], parameters['method'])
        testFunc += parameters['bgrnd']
        fullTestFunc = np.append(fullTestFunc, testFunc)
    return fullTestFunc

##############################################################################


class TensorDeconvWindow(TabFittingWindow):

    def __init__(self, father, oldMainWindow):
        self.CURRENTWINDOW = TensorDeconvFrame
        self.PARAMFRAME = TensorDeconvParamFrame
        super(TensorDeconvWindow, self).__init__(father, oldMainWindow)

#####################################################################################


class TensorDeconvFrame(FitPlotFrame):

    FITNUM = 10  # Maximum number of fits
    
    def __init__(self, rootwindow, fig, canvas, current):
        self.pickNum = 0
        self.pickNum2 = 0
        super(TensorDeconvFrame, self).__init__(rootwindow, fig, canvas, current)

    def togglePick(self, var):
        self.peakPickReset()
        if var == 1:
            self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos)
            self.peakPick = True
        else:
            self.peakPickFunc = None
            self.peakPick = False

    def pickDeconv(self, pos):
        printStr = "%#." + str(self.rootwindow.tabWindow.PRECIS) + "g"
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        if self.pickNum2 == 0:
            if self.pickNum < self.FITNUM:
                self.rootwindow.paramframe.numExp.setCurrentIndex(self.pickNum)
                self.rootwindow.paramframe.changeNum()
            self.rootwindow.paramframe.entries['t11'][self.pickNum].setText(printStr % (self.xax()[pos[0]] * axMult))
            self.pickNum2 = 1
        elif self.pickNum2 == 1:
            self.rootwindow.paramframe.entries['t22'][self.pickNum].setText(printStr % (self.xax()[pos[0]] * axMult))
            self.pickNum2 = 2
        elif self.pickNum2 == 2:
            self.rootwindow.paramframe.entries['t33'][self.pickNum].setText(printStr % (self.xax()[pos[0]] * axMult))
            self.pickNum2 = 0
            self.pickNum += 1
        if self.pickNum < self.FITNUM:
            self.peakPickFunc = lambda pos, self=self: self.pickDeconv(pos)
            self.peakPick = True

#################################################################################


class TensorDeconvParamFrame(AbstractParamFrame):

    def __init__(self, parent, rootwindow, isMain=True):
        self.SINGLENAMES = ['bgrnd', 'spinspeed']
        self.MULTINAMES = ['t11', 't22', 't33', 'amp', 'lor', 'gauss']
        self.PARAMTEXT = {'bgrnd': 'Background', 'spinspeed': 'Spinning Speed', 't11': 'T11', 't22': 'T22', 't33': 'T33', 'amp': 'Integral', 'lor': 'Lorentz', 'gauss': 'Gauss'}
        self.FITFUNC = tensorDeconvmpFit

        # Get full integral
        self.fullInt = np.sum(parent.getData1D()) * parent.sw() / float(len(parent.getData1D()))

        self.cheng = 15
        super(TensorDeconvParamFrame, self).__init__(parent, rootwindow, isMain)
        resetButton = QtWidgets.QPushButton("Reset")
        resetButton.clicked.connect(self.reset)
        self.frame1.addWidget(resetButton, 1, 1)
        self.pickTick = QtWidgets.QCheckBox("Pick")
        self.pickTick.stateChanged.connect(self.togglePick)
        self.frame1.addWidget(self.pickTick, 2, 1)
        self.ticks = {'bgrnd': [], 'spinspeed': [], 't11': [], 't22': [], 't33': [], 'amp': [], 'lor': [], 'gauss': []}
        self.entries = {'bgrnd': [], 'spinspeed': [], 't11': [], 't22': [], 't33': [], 'amp': [], 'lor': [], 'gauss': [], 'shiftdef': [], 'cheng': [], 'mas': [], 'numssb': []}

        self.optframe.addWidget(wc.QLabel("Cheng:"), 0, 0)
        self.entries['cheng'].append(QtWidgets.QSpinBox())
        self.entries['cheng'][-1].setAlignment(QtCore.Qt.AlignHCenter)
        self.entries['cheng'][-1].setValue(self.cheng)
        self.optframe.addWidget(self.entries['cheng'][-1], 1, 0)
        self.entries['mas'].append(QtWidgets.QCheckBox('Spinning'))
        self.entries['mas'][-1].stateChanged.connect(self.MASChange)
        self.optframe.addWidget(self.entries['mas'][-1], 2, 0)
        self.sidebandLabel = wc.QLabel("# sidebands:")
        self.optframe.addWidget(self.sidebandLabel, 3, 0)
        self.sidebandLabel.setEnabled(False)
        self.entries['numssb'].append(QtWidgets.QSpinBox())
        self.entries['numssb'][-1].setAlignment(QtCore.Qt.AlignHCenter)
        self.entries['numssb'][-1].setValue(32)
        self.entries['numssb'][-1].setMaximum(100000)
        self.entries['numssb'][-1].setEnabled(False)
        self.optframe.addWidget(self.entries['numssb'][-1], 4, 0)
        self.shiftDefType = 0  # variable to remember the selected tensor type
        self.optframe.addWidget(wc.QLabel("Definition:"), 5, 0)
        self.entries['shiftdef'].append(QtWidgets.QComboBox())
        self.entries['shiftdef'][-1].addItems([u'\u03b411 - \u03b422 - \u03b433',
                                               u'\u03b4xx - \u03b4yy - \u03b4zz',
                                               u'\u03b4iso - \u03b4aniso - \u03b7',
                                               u'\u03b4iso - \u03a9 - \u03b7'])
        self.entries['shiftdef'][-1].currentIndexChanged.connect(self.changeShiftDef)
        self.optframe.addWidget(self.entries['shiftdef'][-1], 6, 0)
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.spinLabel = wc.QLabel("Spin. speed [kHz]:")
        self.spinLabel.setEnabled(False)
        self.frame2.addWidget(self.spinLabel, 0, 0, 1, 2)
        self.ticks['spinspeed'].append(QtWidgets.QCheckBox(''))
        self.ticks['spinspeed'][-1].setEnabled(False)
        self.frame2.addWidget(self.ticks['spinspeed'][-1], 1, 0)
        self.entries['spinspeed'].append(wc.FitQLineEdit(self, 'spinspeed', "10.0"))
        self.frame2.addWidget(self.entries['spinspeed'][-1], 1, 1)
        self.entries['spinspeed'][-1].setEnabled(False)
        self.frame2.addWidget(wc.QLabel("Bgrnd:"), 2, 0, 1, 2)
        self.ticks['bgrnd'].append(QtWidgets.QCheckBox(''))
        self.frame2.addWidget(self.ticks['bgrnd'][-1], 3, 0)
        self.entries['bgrnd'].append(wc.FitQLineEdit(self, 'bgrnd', "0.0"))
        self.frame2.addWidget(self.entries['bgrnd'][-1], 3, 1)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtWidgets.QComboBox()
        self.numExp.addItems([str(x + 1) for x in range(self.FITNUM)])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        if self.parent.viewSettings["ppm"]:
            axUnit = 'ppm'
        else:
            axUnit = ['Hz', 'kHz', 'MHz'][self.parent.getAxType()]
        # Labels
        self.label11 = wc.QLabel(u'\u03b4' + '<sub>11</sub> [' + axUnit + '] :')
        self.label22 = wc.QLabel(u'\u03b4' + '<sub>22</sub> [' + axUnit + '] :')
        self.label33 = wc.QLabel(u'\u03b4' + '<sub>33</sub> [' + axUnit + '] :')
        self.frame3.addWidget(self.label11, 1, 0, 1, 2)
        self.frame3.addWidget(self.label22, 1, 2, 1, 2)
        self.frame3.addWidget(self.label33, 1, 4, 1, 2)
        self.labelxx = wc.QLabel(u'\u03b4' + '<sub>xx</sub> [' + axUnit + '] :')
        self.labelyy = wc.QLabel(u'\u03b4' + '<sub>yy</sub> [' + axUnit + '] :')
        self.labelzz = wc.QLabel(u'\u03b4' + '<sub>zz</sub> [' + axUnit + '] :')
        self.labelxx.hide()
        self.labelyy.hide()
        self.labelzz.hide()
        self.frame3.addWidget(self.labelxx, 1, 0, 1, 2)
        self.frame3.addWidget(self.labelyy, 1, 2, 1, 2)
        self.frame3.addWidget(self.labelzz, 1, 4, 1, 2)
        self.labeliso = wc.QLabel(u'\u03b4' + '<sub>iso</sub> [' + axUnit + '] :')
        self.labelaniso = wc.QLabel(u'\u03b4' + '<sub>aniso</sub> [' + axUnit + '] :')
        self.labeleta = wc.QLabel(u'\u03b7:')
        self.labeliso.hide()
        self.labelaniso.hide()
        self.labeleta.hide()
        self.frame3.addWidget(self.labeliso, 1, 0, 1, 2)
        self.frame3.addWidget(self.labelaniso, 1, 2, 1, 2)
        self.frame3.addWidget(self.labeleta, 1, 4, 1, 2)
        self.labeliso2 = wc.QLabel(u'\u03b4' + '<sub>iso</sub> [' + axUnit + '] :')
        self.labelspan = wc.QLabel(u'\u03a9 [' + axUnit + '] :')
        self.labelskew = wc.QLabel(u'\u03ba:')
        self.labeliso2.hide()
        self.labelspan.hide()
        self.labelskew.hide()
        self.frame3.addWidget(self.labeliso2, 1, 0, 1, 2)
        self.frame3.addWidget(self.labelspan, 1, 2, 1, 2)
        self.frame3.addWidget(self.labelskew, 1, 4, 1, 2)
        self.frame3.addWidget(wc.QLabel("Integral:"), 1, 6, 1, 2)
        self.frame3.addWidget(wc.QLabel("Lorentz [Hz]:"), 1, 8, 1, 2)
        self.frame3.addWidget(wc.QLabel("Gauss [Hz]:"), 1, 10, 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        for i in range(self.FITNUM):
            for j in range(len(self.MULTINAMES)):
                self.ticks[self.MULTINAMES[j]].append(QtWidgets.QCheckBox(''))
                self.frame3.addWidget(self.ticks[self.MULTINAMES[j]][i], i + 2, 2 * j)
                self.entries[self.MULTINAMES[j]].append(wc.FitQLineEdit(self, self.MULTINAMES[j]))
                self.frame3.addWidget(self.entries[self.MULTINAMES[j]][i], i + 2, 2 * j + 1)
        self.reset()

    def MASChange(self, state):
        if state:  # When turned on
            self.entries['spinspeed'][-1].setEnabled(True)
            self.ticks['spinspeed'][-1].setEnabled(True)
            self.spinLabel.setEnabled(True)
            self.entries['numssb'][-1].setEnabled(True)
            self.sidebandLabel.setEnabled(True)
        else:
            self.entries['spinspeed'][-1].setEnabled(False)
            self.ticks['spinspeed'][-1].setEnabled(False)
            self.spinLabel.setEnabled(False)
            self.entries['numssb'][-1].setEnabled(False)
            self.sidebandLabel.setEnabled(False)

    def defaultValues(self, inp):
        if not inp:
            return {'bgrnd': [0.0, True],
                    'spinspeed': [10.0, True],
                    't11': np.repeat([np.array([0.0, False], dtype=object)], self.FITNUM, axis=0),
                    't22': np.repeat([np.array([0.0, False], dtype=object)], self.FITNUM, axis=0),
                    't33': np.repeat([np.array([0.0, False], dtype=object)], self.FITNUM, axis=0),
                    'amp': np.repeat([np.array([self.fullInt, False], dtype=object)], self.FITNUM, axis=0),
                    'lor': np.repeat([np.array([1.0, False], dtype=object)], self.FITNUM, axis=0),
                    'gauss': np.repeat([np.array([0.0, True], dtype=object)], self.FITNUM, axis=0)}
        else:
            return inp

    def reset(self):
        locList = self.getRedLocList()
        self.parent.pickNum = 0
        self.parent.pickNum2 = 0
        self.cheng = 15
        self.entries['cheng'][-1].setValue(self.cheng)
        self.fitNumList[locList] = 0
        self.fitParamList[locList] = self.defaultValues(0)
        self.pickTick.setChecked(True)
        self.togglePick()
        self.dispParams()

    def togglePick(self):
        self.parent.togglePick(self.pickTick.isChecked())

    def changeShiftDef(self):
        NewType = self.entries['shiftdef'][-1].currentIndex()
        OldType = self.shiftDefType
        if NewType == 0:
            self.label11.show()
            self.label22.show()
            self.label33.show()
            self.labelxx.hide()
            self.labelyy.hide()
            self.labelzz.hide()
            self.labeliso.hide()
            self.labelaniso.hide()
            self.labeleta.hide()
            self.labeliso2.hide()
            self.labelspan.hide()
            self.labelskew.hide()
            self.pickTick.setChecked(True)
            self.pickTick.show()
        elif NewType == 1:
            self.label11.hide()
            self.label22.hide()
            self.label33.hide()
            self.labelxx.show()
            self.labelyy.show()
            self.labelzz.show()
            self.labeliso.hide()
            self.labelaniso.hide()
            self.labeleta.hide()
            self.labeliso2.hide()
            self.labelspan.hide()
            self.labelskew.hide()
            self.pickTick.setChecked(True)
            self.pickTick.show()
        elif NewType == 2:
            self.label11.hide()
            self.label22.hide()
            self.label33.hide()
            self.labelxx.hide()
            self.labelyy.hide()
            self.labelzz.hide()
            self.labeliso.show()
            self.labelaniso.show()
            self.labeleta.show()
            self.labeliso2.hide()
            self.labelspan.hide()
            self.labelskew.hide()
            self.pickTick.setChecked(False)
            self.pickTick.hide()
        elif NewType == 3:
            self.label11.hide()
            self.label22.hide()
            self.label33.hide()
            self.labelxx.hide()
            self.labelyy.hide()
            self.labelzz.hide()
            self.labeliso.hide()
            self.labelaniso.hide()
            self.labeleta.hide()
            self.labeliso2.show()
            self.labelspan.show()
            self.labelskew.show()
            self.pickTick.setChecked(False)
            self.pickTick.hide()
        val = self.numExp.currentIndex() + 1
        tensorList = []
        for i in range(10):  # Convert input
            if i < val:
                T11 = safeEval(self.entries['t11'][i].text())
                T22 = safeEval(self.entries['t22'][i].text())
                T33 = safeEval(self.entries['t33'][i].text())
                startTensor = [T11, T22, T33]
                if None in startTensor:
                    self.entries['shiftdef'][-1].setCurrentIndex(OldType)  # error, reset to old view
                    self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
                    return
                Tensors = func.shiftConversion(startTensor, OldType)
                for element in range(3):  # Check for `ND' s
                    if isinstance(Tensors[NewType][element], str):
                        Tensors[NewType][element] = 0
                tensorList.append(Tensors)
        printStr = '%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g'
        for i in range(10):  # Print output if not stopped before
            if i < val:
                self.entries['t11'][i].setText(printStr % tensorList[i][NewType][0])
                self.entries['t22'][i].setText(printStr % tensorList[i][NewType][1])
                self.entries['t33'][i].setText(printStr % tensorList[i][NewType][2])
        self.shiftDefType = NewType

    def getExtraParams(self, out):
        out['mas'] = [self.entries['mas'][-1].isChecked()]
        cheng = self.entries['cheng'][-1].value()
        out['cheng'] = [cheng]
        out['shiftdef'] = [self.entries['shiftdef'][0].currentIndex()]
        if out['mas'][0]:
            out['numssb'] = [self.entries['numssb'][0].value()]
            return (out, [out['mas'][-1], out['cheng'][-1], out['shiftdef'][-1], out['numssb'][-1]])
        else:
            weight, angleStuff = simFunc.csaAngleStuff(cheng)
            out['weight'] = [weight]
            out['multt'] =  angleStuff
            return (out, [out['mas'][-1], out['multt'][-1], out['weight'][-1], out['shiftdef'][-1]])

    def disp(self, params, num):
        out = params[num]
        for name in self.SINGLENAMES:
            inp = out[name][0]
            if isinstance(inp, tuple):
                inp = checkLinkTuple(inp)
                out[name][0] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
        numExp = len(out[self.MULTINAMES[0]])
        for i in range(numExp):
            for name in self.MULTINAMES:
                inp = out[name][i]
                if isinstance(inp, tuple):
                    inp = checkLinkTuple(inp)
                    out[name][i] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                if not np.isfinite(out[name][i]):
                    self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
                    return
        tmpx = self.parent.xax()
        bgrnd = out['bgrnd'][0]
        outCurve = bgrnd * np.ones(len(tmpx))
        outCurvePart = []
        x = []
        for i in range(len(out['amp'])):
            x.append(tmpx)
            t11 = out['t11'][i] / self.axMult
            t22 = out['t22'][i] / self.axMult
            if out['shiftdef'][-1] in [0, 1]:
                t33 = out['t33'][i] / self.axMult
            else:
                t33 = out['t33'][i]
            tensor = func.shiftConversion([t11, t22, t33], out['shiftdef'][-1])
            if out['mas'][0]:
                y = out['amp'][i] * simFunc.tensorMASDeconvtensorFunc(tmpx, tensor[2][0], tensor[2][1], tensor[2][2], out['lor'][i], out['gauss'][i], self.parent.sw(), out['spinspeed'][0], out['cheng'][0], out['numssb'][-1])
            else:
                y = out['amp'][i] * simFunc.tensorDeconvtensorFunc(tmpx, tensor[0][0], tensor[0][1], tensor[0][2], out['lor'][i], out['gauss'][i], out['multt'][0], self.parent.sw(), out['weight'][0])
            y = np.real(np.fft.fftn(y))
            outCurvePart.append(bgrnd + y)
            outCurve += y
        locList = self.getRedLocList()
        self.parent.fitDataList[locList] = [tmpx, outCurve, x, outCurvePart]
        self.parent.showFid()
        self.changeShiftDef()  # Reformat output to correct display

##############################################################################


def tensorDeconvmpFit(xax, data1D, guess, args, queue, minmethod, numfeval):
    try:
        fitVal = scipy.optimize.minimize(lambda *param: np.sum((data1D - tensorDeconvfitFunc(param, xax, args))**2), guess, method=minmethod, options = {'maxfev': numfeval})
    except Exception:
        fitVal = None
    queue.put(fitVal)


def tensorDeconvfitFunc(params, allX, args):
    params = params[0]
    specName = args[0]
    specSlices = args[1]
    allParam = []
    for length in specSlices:
        allParam.append(params[length])
    allStruc = args[3]
    allArgu = args[4]
    fullTestFunc = []
    for n in range(len(allX)):
        x = allX[n][-1]
        testFunc = np.zeros(len(x))
        param = allParam[n]
        numExp = args[2][n]
        struc = args[3][n]
        argu = args[4][n]
        sw = args[5][n][-1]
        axMult = args[6][n]
        parameters = {}
        mas = argu[-1][0]
        if mas:
            parameters['cheng'] = argu[-1][1]
            parameters['shiftdef'] = argu[-1][2]
            parameters['numssb'] = argu[-1][3]
        else:
            parameters['multt'] = argu[-1][1]
            parameters['weight'] = argu[-1][2]
            parameters['shiftdef'] = argu[-1][3]
        for name in ['spinspeed', 'bgrnd']:
            if struc[name][0][0] == 1:
                parameters[name] = param[struc[name][0][1]]
            elif struc[name][0][0] == 0:
                parameters[name] = argu[struc[name][0][1]]
            else:
                altStruc = struc[name][0][1]
                if struc[altStruc[0]][altStruc[1]][0] == 1:
                    parameters[name] = altStruc[2] * allParam[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                elif struc[altStruc[0]][altStruc[1]][0] == 0:
                    parameters[name] = altStruc[2] * allArgu[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
        for i in range(numExp):
            for name in ['t11', 't22', 't33', 'amp', 'lor', 'gauss']:
                if struc[name][i][0] == 1:
                    parameters[name] = param[struc[name][i][1]]
                elif struc[name][i][0] == 0:
                    parameters[name] = argu[struc[name][i][1]]
                else:
                    altStruc = struc[name][i][1]
                    strucTarget = allStruc[altStruc[4]]
                    if strucTarget[altStruc[0]][altStruc[1]][0] == 1:
                        parameters[name] = altStruc[2] * allParam[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                    elif strucTarget[altStruc[0]][altStruc[1]][0] == 0:
                        parameters[name] = altStruc[2] * allArgu[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
            t11 = parameters['t11'] / axMult
            t22 = parameters['t22'] / axMult
            if parameters['shiftdef'] in [0, 1]:
                t33 = parameters['t33'] / axMult
            else:
                t33 = parameters['t33']
            tensor = func.shiftConversion([t11, t22, t33], parameters['shiftdef'])
            if mas:
                testFunc += parameters['amp'] * simFunc.tensorMASDeconvtensorFunc(x, tensor[2][0], tensor[2][1], tensor[2][2], parameters['lor'], parameters['gauss'], sw, parameters['spinspeed'], parameters['cheng'], parameters['numssb'])
            else:
                testFunc += parameters['amp'] * simFunc.tensorDeconvtensorFunc(x, tensor[0][0], tensor[0][1], tensor[0][2], parameters['lor'], parameters['gauss'], parameters['multt'], sw, parameters['weight'])
        testFunc = np.real(np.fft.fftn(testFunc))
        testFunc += parameters['bgrnd']
        fullTestFunc = np.append(fullTestFunc, testFunc)
    return fullTestFunc


##############################################################################


class Quad1DeconvWindow(TabFittingWindow):

    def __init__(self, father, oldMainWindow):
        self.CURRENTWINDOW = Quad1DeconvFrame
        self.PARAMFRAME = Quad1DeconvParamFrame
        super(Quad1DeconvWindow, self).__init__(father, oldMainWindow)

#################################################################################


class Quad1DeconvFrame(FitPlotFrame):

    FITNUM = 10  # Maximum number of fits

#################################################################################


class Quad1DeconvParamFrame(AbstractParamFrame):

    Ioptions = ['1', '3/2', '2', '5/2', '3', '7/2', '4', '9/2']
    Ivalues = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]

    def __init__(self, parent, rootwindow, isMain=True):
        self.SINGLENAMES = ['bgrnd', 'spinspeed']
        self.MULTINAMES = ['pos', 'cq', 'eta', 'amp', 'lor', 'gauss']
        self.PARAMTEXT = {'bgrnd': 'Background', 'spinspeed': 'Spinning Speed', 'pos': 'Position', 'cq': 'Cq', 'eta': 'eta', 'amp': 'Integral', 'lor': 'Lorentz', 'gauss': 'Gauss'}
        self.FITFUNC = quad1mpFit
        self.setAngleStuff = simFunc.quad1DeconvsetAngleStuff
        self.tensorFunc = simFunc.quad1DeconvtensorFunc
        self.cheng = 15
        # Get full integral
        self.fullInt = np.sum(parent.getData1D()) * parent.sw() / float(len(parent.getData1D()))
        super(Quad1DeconvParamFrame, self).__init__(parent, rootwindow, isMain)
        self.ticks = {'bgrnd': [], 'spinspeed': [], 'pos': [], 'cq': [], 'eta': [], 'amp': [], 'lor': [], 'gauss': []}
        self.entries = {'bgrnd': [], 'spinspeed': [], 'pos': [], 'cq': [], 'eta': [], 'amp': [], 'lor': [], 'gauss': [], 'shiftdef': [], 'cheng': [], 'I': [], 'mas': [], 'numssb': []}
        self.optframe.addWidget(wc.QLabel("Cheng:"), 0, 0)
        self.entries['cheng'].append(QtWidgets.QSpinBox())
        self.entries['cheng'][-1].setAlignment(QtCore.Qt.AlignHCenter)
        self.entries['cheng'][-1].setValue(self.cheng)
        self.optframe.addWidget(self.entries['cheng'][-1], 1, 0)
        self.optframe.addWidget(wc.QLabel("I:"), 2, 0)
        self.entries['I'].append(QtWidgets.QComboBox())
        self.entries['I'][-1].addItems(self.Ioptions)
        self.entries['I'][-1].setCurrentIndex(1)
        self.optframe.addWidget(self.entries['I'][-1], 3, 0)
        self.entries['mas'].append(QtWidgets.QCheckBox('Spinning'))
        self.entries['mas'][-1].stateChanged.connect(self.MASChange)
        self.optframe.addWidget(self.entries['mas'][-1], 4, 0)
        self.sidebandLabel = wc.QLabel("# sidebands:")
        self.sidebandLabel.setEnabled(False)
        self.optframe.addWidget(self.sidebandLabel, 5, 0)
        self.entries['numssb'].append(QtWidgets.QSpinBox())
        self.entries['numssb'][-1].setAlignment(QtCore.Qt.AlignHCenter)
        self.entries['numssb'][-1].setValue(32)
        self.entries['numssb'][-1].setMaximum(100000)
        self.entries['numssb'][-1].setEnabled(False)
        self.optframe.addWidget(self.entries['numssb'][-1], 6, 0)
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.spinLabel = wc.QLabel("Spin. speed [kHz]:")
        self.frame2.addWidget(self.spinLabel, 0, 0, 1, 2)
        self.spinLabel.setEnabled(False)
        self.ticks['spinspeed'].append(QtWidgets.QCheckBox(''))
        self.frame2.addWidget(self.ticks['spinspeed'][-1], 1, 0)
        self.ticks['spinspeed'][-1].setEnabled(False)
        self.entries['spinspeed'].append(wc.QLineEdit("10.0"))
        self.frame2.addWidget(self.entries['spinspeed'][-1], 1, 1)
        self.entries['spinspeed'][-1].setEnabled(False)
        self.frame2.addWidget(wc.QLabel("Bgrnd:"), 2, 0, 1, 2)
        self.ticks['bgrnd'].append(QtWidgets.QCheckBox(''))
        self.frame2.addWidget(self.ticks['bgrnd'][-1], 3, 0)
        self.entries['bgrnd'].append(wc.QLineEdit("0.0"))
        self.frame2.addWidget(self.entries['bgrnd'][-1], 3, 1)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtWidgets.QComboBox()
        self.numExp.addItems([str(x + 1) for x in range(self.FITNUM)])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        if self.parent.viewSettings["ppm"]:
            axUnit = 'ppm'
        else:
            axUnit = ['Hz', 'kHz', 'MHz'][self.parent.getAxType()]
        # Labels
        self.labelpos = wc.QLabel(u'Position [' + axUnit + ']:')
        self.labelcq = wc.QLabel(u'C<sub>Q</sub> [MHz]:')
        self.labeleta = wc.QLabel(u'\u03B7:')
        self.frame3.addWidget(self.labelpos, 1, 0, 1, 2)
        self.frame3.addWidget(self.labelcq, 1, 2, 1, 2)
        self.frame3.addWidget(self.labeleta, 1, 4, 1, 2)
        self.frame3.addWidget(wc.QLabel("Integral:"), 1, 6, 1, 2)
        self.frame3.addWidget(wc.QLabel("Lorentz [Hz]:"), 1, 8, 1, 2)
        self.frame3.addWidget(wc.QLabel("Gauss [Hz]:"), 1, 10, 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        for i in range(self.FITNUM):
            for j in range(len(self.MULTINAMES)):
                self.ticks[self.MULTINAMES[j]].append(QtWidgets.QCheckBox(''))
                self.frame3.addWidget(self.ticks[self.MULTINAMES[j]][i], i + 2, 2 * j)
                self.entries[self.MULTINAMES[j]].append(wc.FitQLineEdit(self, self.MULTINAMES[j]))
                self.frame3.addWidget(self.entries[self.MULTINAMES[j]][i], i + 2, 2 * j + 1)
        self.dispParams()

    def MASChange(self, state):
        if state:  # When turned on
            self.entries['spinspeed'][-1].setEnabled(True)
            self.ticks['spinspeed'][-1].setEnabled(True)
            self.spinLabel.setEnabled(True)
            self.entries['numssb'][-1].setEnabled(True)
            self.sidebandLabel.setEnabled(True)
        else:
            self.entries['spinspeed'][-1].setEnabled(False)
            self.ticks['spinspeed'][-1].setEnabled(False)
            self.spinLabel.setEnabled(False)
            self.entries['numssb'][-1].setEnabled(False)
            self.sidebandLabel.setEnabled(False)

    def defaultValues(self, inp):
        if not inp:
            return {'bgrnd': [0.0, True],
                    'spinspeed': [10.0, True],
                    'pos': np.repeat([np.array([0.0, False], dtype=object)], self.FITNUM, axis=0),
                    'cq': np.repeat([np.array([1.0, False], dtype=object)], self.FITNUM, axis=0),
                    'eta': np.repeat([np.array([0.0, False], dtype=object)], self.FITNUM, axis=0),
                    'amp': np.repeat([np.array([self.fullInt, False], dtype=object)], self.FITNUM, axis=0),
                    'lor': np.repeat([np.array([1.0, False], dtype=object)], self.FITNUM, axis=0),
                    'gauss': np.repeat([np.array([0.0, True], dtype=object)], self.FITNUM, axis=0)}
        else:
            return inp

    def checkI(self, I):
        return I * 0.5 + 1

    def getExtraParams(self, out):
        out['mas'] = [self.entries['mas'][-1].isChecked()]
        out['I'] = [self.checkI(self.entries['I'][-1].currentIndex())]
        cheng = safeEval(self.entries['cheng'][-1].text())
        out['cheng'] = [cheng]
        if out['mas'][0]:
            out['numssb'] = [self.entries['numssb'][0].value()]
            return (out, [out['mas'][-1], out['I'][-1], out['cheng'][-1], out['numssb'][-1]])
        else:
            weight, angleStuff = self.setAngleStuff(cheng)
            out['weight'] = [weight]
            out['anglestuff'] = [angleStuff]
            out['tensorfunc'] = [self.tensorFunc]
            out['freq'] = [self.parent.freq()]
            return (out, [out['mas'][-1], out['I'][-1], out['weight'][-1], out['anglestuff'][-1], out['tensorfunc'][-1], out['freq'][-1]])


    def checkParam(self):
         val = self.numExp.currentIndex() + 1
         printStr = '%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g'
         for i in range(10):  # Print output if not stopped before
            if i < val:
                try:
                    cq = float(safeEval(self.entries['cq'][i].text()))
                    eta = float(safeEval(self.entries['eta'][i].text()))
                    Res = func.quadConversion([cq,eta], 1, 0)[0]
                    self.entries['cq'][i].setText(printStr % Res[0])
                    self.entries['eta'][i].setText(printStr % Res[1])
                except Exception:
                    return

    def disp(self, params, num):
        out = params[num]
        for name in self.SINGLENAMES:
            inp = out[name][0]
            if isinstance(inp, tuple):
                inp = checkLinkTuple(inp)
                out[name][0] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
        numExp = len(out[self.MULTINAMES[0]])
        for i in range(numExp):
            for name in self.MULTINAMES:
                inp = out[name][i]
                if isinstance(inp, tuple):
                    inp = checkLinkTuple(inp)
                    out[name][i] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                if not np.isfinite(out[name][i]):
                    self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
                    return
        tmpx = self.parent.xax()
        bgrnd = out['bgrnd'][0]
        outCurve = bgrnd * np.ones(len(tmpx))
        outCurvePart = []
        x = []
        for i in range(len(out['amp'])):
            x.append(tmpx)
            if out['mas'][0]:
                y = out['amp'][i] * simFunc.quad1MASFunc(tmpx, out['pos'][i]/self.axMult, out['cq'][i], out['eta'][i], out['lor'][i], out['gauss'][i], self.parent.sw(), out['spinspeed'][0], out['cheng'][0], out['I'][0], out['numssb'][0])
            else:
                y = out['amp'][i] * self.tensorFunc(tmpx, out['I'][0], out['pos'][i]/self.axMult, out['cq'][i], out['eta'][i], out['lor'][i], out['gauss'][i], out['anglestuff'][0], self.parent.freq(), self.parent.sw(), out['weight'][0])
            y = np.real(np.fft.fftn(y))
            outCurvePart.append(bgrnd + y)
            outCurve += y
        locList = self.getRedLocList()
        self.parent.fitDataList[locList] = [tmpx, outCurve, x, outCurvePart]
        self.parent.showFid()
        self.checkParam()


##############################################################################


def quad1mpFit(xax, data1D, guess, args, queue, minmethod, numfeval):
    try:
        fitVal = scipy.optimize.minimize(lambda *param: np.sum((data1D - quad1fitFunc(param, xax, args))**2), guess, method=minmethod, options = {'maxfev': numfeval})
    except Exception:
        fitVal = None
    queue.put(fitVal)


def quad1fitFunc(params, allX, args):
    params = params[0]
    specName = args[0]
    specSlices = args[1]
    allParam = []
    for length in specSlices:
        allParam.append(params[length])
    allStruc = args[3]
    allArgu = args[4]
    fullTestFunc = []
    for n in range(len(allX)):
        x = allX[n][-1]
        testFunc = np.zeros(len(x))
        param = allParam[n]
        numExp = args[2][n]
        struc = args[3][n]
        argu = args[4][n]
        sw = args[5][n][-1]
        axMult = args[6][n]
        parameters = {}
        mas = argu[-1][0]
        parameters['I'] = argu[-1][1]
        if mas:
            parameters['cheng'] = argu[-1][2]
            parameters['numssb'] = argu[-1][3]
        else:
            parameters['weight'] = argu[-1][2]
            parameters['anglestuff'] = argu[-1][3]
            parameters['tensorfunc'] = argu[-1][4]
            parameters['freq'] = argu[-1][5]
        for name in ['spinspeed', 'bgrnd']:
            if struc[name][0][0] == 1:
                parameters[name] = param[struc[name][0][1]]
            elif struc[name][0][0] == 0:
                parameters[name] = argu[struc[name][0][1]]
            else:
                altStruc = struc[name][0][1]
                if struc[altStruc[0]][altStruc[1]][0] == 1:
                    parameters[name] = altStruc[2] * allParam[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                elif struc[altStruc[0]][altStruc[1]][0] == 0:
                    parameters[name] = altStruc[2] * allArgu[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
        for i in range(numExp):
            for name in ['pos', 'cq', 'eta', 'amp', 'lor', 'gauss']:
                if struc[name][i][0] == 1:
                    parameters[name] = param[struc[name][i][1]]
                elif struc[name][i][0] == 0:
                    parameters[name] = argu[struc[name][i][1]]
                else:
                    altStruc = struc[name][i][1]
                    strucTarget = allStruc[altStruc[4]]
                    if strucTarget[altStruc[0]][altStruc[1]][0] == 1:
                        parameters[name] = altStruc[2] * allParam[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                    elif strucTarget[altStruc[0]][altStruc[1]][0] == 0:
                        parameters[name] = altStruc[2] * allArgu[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
            if mas:
                testFunc += parameters['amp'] * simFunc.quad1MASFunc(x, parameters['pos']/axMult, parameters['cq'], parameters['eta'], parameters['lor'], parameters['gauss'], sw, parameters['spinspeed'], parameters['cheng'], parameters['I'], parameters['numssb'])
            else:
                testFunc += parameters['amp'] * parameters['tensorfunc'](x, parameters['I'], parameters['pos']/axMult, parameters['cq'], parameters['eta'], parameters['lor'], parameters['gauss'], parameters['anglestuff'], parameters['freq'], sw, parameters['weight'])
            testFunc = np.real(np.fft.fftn(testFunc))
            testFunc += parameters['bgrnd']
        fullTestFunc = np.append(fullTestFunc, testFunc)
    return fullTestFunc


##############################################################################


class Quad2DeconvWindow(TabFittingWindow):

    def __init__(self, father, oldMainWindow):
        self.CURRENTWINDOW = Quad1DeconvFrame
        self.PARAMFRAME = Quad2DeconvParamFrame
        super(Quad2DeconvWindow, self).__init__(father, oldMainWindow)

#################################################################################


class Quad2DeconvParamFrame(Quad1DeconvParamFrame):

    Ioptions = ['3/2', '5/2', '7/2', '9/2']

    def __init__(self, parent, rootwindow, isMain=True):
        super(Quad2DeconvParamFrame, self).__init__(parent, rootwindow, isMain)
        self.setAngleStuff = simFunc.quad2StaticsetAngleStuff
        self.tensorFunc = simFunc.quad2tensorFunc
        self.entries['I'][-1].setCurrentIndex(0)
        self.spinLabel.hide()
        self.ticks['spinspeed'][-1].hide()
        self.entries['spinspeed'][-1].hide()
        self.sidebandLabel.hide()
        self.entries['numssb'][-1].hide()

    def getExtraParams(self, out):
        out['mas'] = [0]  # Second order quadrupole MAS is calculated as if it is static
        if self.entries['mas'][-1].isChecked():
            self.setAngleStuff = simFunc.quad2MASsetAngleStuff
        else:
            self.setAngleStuff = simFunc.quad2StaticsetAngleStuff
        out['I'] = [self.checkI(self.entries['I'][-1].currentIndex())]
        cheng = safeEval(self.entries['cheng'][-1].text())
        out['cheng'] = [cheng]
        weight, angleStuff = self.setAngleStuff(cheng)
        out['weight'] = [weight]
        out['anglestuff'] = [angleStuff]
        out['tensorfunc'] = [self.tensorFunc]
        out['freq'] = [self.parent.freq()]
        return (out, [out['mas'][-1], out['I'][-1], out['weight'][-1], out['anglestuff'][-1], out['tensorfunc'][-1], out['freq'][-1]])

    def checkI(self, I):
        return I * 1.0 + 1.5

##############################################################################


class CzjzekPrefWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        super(CzjzekPrefWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Grid Settings")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel(u"\u03c9<sub>Q</sub> grid size:"), 0, 0)
        self.wqsteps = QtWidgets.QSpinBox()
        self.wqsteps.setMinimum(2)
        self.wqsteps.setMaximum(1000)
        self.wqsteps.setAlignment(QtCore.Qt.AlignHCenter) 
        self.wqsteps.setValue(self.father.wqsteps)
        grid.addWidget(self.wqsteps, 1, 0)
        grid.addWidget(wc.QLabel(u"\u03b7 grid size:"), 2, 0)
        self.etasteps = QtWidgets.QSpinBox()
        self.etasteps.setMinimum(2)
        self.etasteps.setMaximum(1000)
        self.etasteps.setAlignment(QtCore.Qt.AlignHCenter) 
        self.etasteps.setValue(self.father.etasteps)
        grid.addWidget(self.etasteps, 3, 0)
        grid.addWidget(wc.QLabel(u"\u03BD<sub>Q</sub><sup>min</sup> [MHz]:"), 4, 0)
        self.wqmin = wc.QLineEdit(str(self.father.wqmin), self.checkWq)
        grid.addWidget(self.wqmin, 5, 0)
        grid.addWidget(wc.QLabel(u"\u03BD<sub>Q</sub><sup>max</sup> [MHz]:"), 6, 0)
        self.wqmax = wc.QLineEdit(str(self.father.wqmax), self.checkWq)
        grid.addWidget(self.wqmax, 7, 0)
        grid.addWidget(wc.QLabel(u"\u03B7<sup>min</sup>:"), 8, 0)
        self.etamin = wc.QLineEdit(str(self.father.etamin), self.checkEta)
        grid.addWidget(self.etamin, 9, 0)
        grid.addWidget(wc.QLabel(u"\u03B7<sup>max</sup>:"), 10, 0)
        self.etamax = wc.QLineEdit(str(self.father.etamax), self.checkEta)
        grid.addWidget(self.etamax, 11, 0)
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        grid.addWidget(self.canvas, 0, 1, 13,4)
        self.ax = self.fig.add_subplot(111)
   
        self.simButton = QtWidgets.QPushButton("Sim", parent=self)
        self.simButton.clicked.connect(self.plotDist)
        grid.addWidget(self.simButton, 13, 3)


        grid.addWidget(wc.QLabel("Site:"), 13, 1)
        self.site = QtWidgets.QSpinBox()
        self.site.setMinimum(1)
        self.site.setMaximum(self.father.numExp.currentIndex() + 1)
        self.site.setAlignment(QtCore.Qt.AlignHCenter) 
        grid.addWidget(self.site, 13, 2)

        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 4, 0)
        okButton = QtWidgets.QPushButton("&Ok", self)
        okButton.clicked.connect(self.applyAndClose)
        okButton.setFocus()
        layout.addWidget(okButton, 4, 1)
        grid.setRowStretch(12, 1)
        grid.setColumnStretch(4, 1)
        #self.plotDist()
        self.show()
        self.resize(800, 600)

    def plotDist(self):
        self.ax.cla()
        wqsteps = self.wqsteps.value()
        etasteps = self.etasteps.value()
        wqmax = safeEval(self.wqmax.text(),type='FI')
        wqmin = safeEval(self.wqmin.text(),type='FI')
        etamax = safeEval(self.etamax.text(),type='FI')
        etamin = safeEval(self.etamin.text(),type='FI')
        wq, eta = np.meshgrid(np.linspace(wqmin, wqmax, wqsteps), np.linspace(etamin, etamax, etasteps))
        #self.MULTINAMES = ['d', 'pos', 'sigma','wq0', 'eta0', 'amp', 'lor', 'gauss']
        method = self.father.entries['method'][0].currentIndex()
        site = self.site.value() - 1
        d = safeEval(self.father.entries[self.father.MULTINAMES[0]][site].text(),type='FI')
        sigma = safeEval(self.father.entries[self.father.MULTINAMES[2]][site].text(),type='FI')
        wq0 = safeEval(self.father.entries[self.father.MULTINAMES[3]][site].text(),type='FI')
        eta0 = safeEval(self.father.entries[self.father.MULTINAMES[4]][site].text(),type='FI')
        if method == 0 or (eta0 == 0.0 and wq0 == 0.0):
            czjzek = simFunc.czjzekIntensities(sigma, d, wq.flatten(), eta.flatten(), wq0, eta0)
        elif method == 1:
            czjzek = eCzjz.getInts(sigma, d, eta0, wq0, wq.flatten(), eta.flatten())
        czjzek = czjzek.reshape(etasteps,wqsteps)

        self.ax.contour(wq.transpose(),eta.transpose(),czjzek.transpose(),10)
        self.ax.set_xlabel(u"\u03BD$_Q$ [MHz]")
        self.ax.set_ylabel(u"\u03B7")
        self.canvas.draw()



    def checkWq(self):
        inp = safeEval(self.wqmax.text(),type='FI')
        if inp is None:
            return False
        self.wqmax.setText(str(float(inp)))
        inp = safeEval(self.wqmin.text(),type='FI')
        if inp is None:
            return False
        self.wqmin.setText(str(float(inp)))
        return True

    def checkEta(self):
        inp = safeEval(self.etamax.text(),type='FI')
        if inp is None:
            return False
        if inp < 0.0 or inp > 1.0:
            return False
        self.etamax.setText(str(abs(float(inp))))
        inp = safeEval(self.etamin.text(),type='FI')
        if inp is None:
            return False
        if inp < 0.0 or inp > 1.0:
            return False
        self.etamin.setText(str(abs(float(inp))))


    def closeEvent(self, *args):
        self.deleteLater()

    def applyAndClose(self, *args):
        self.father.wqsteps = self.wqsteps.value()
        self.father.etasteps = self.etasteps.value()
        inp = safeEval(self.wqmax.text(),type='FI')
        if inp is None:
            self.father.rootwindow.father.dispMsg(u"\u03BD_Q_max value not valid.")
            return
        self.father.wqmax = abs(safeEval(self.wqmax.text()))
        inp = abs(safeEval(self.wqmin.text(),type='FI'))
        if inp is None:
            self.father.rootwindow.father.dispMsg(u"\03BD_Q_min value not valid.")
            return
        self.father.wqmin = abs(safeEval(self.wqmin.text()))
        #eta
        inp = safeEval(self.etamax.text(),type='FI')
        if inp is None:
            self.father.rootwindow.father.dispMsg(u"\u03B7_max value not valid.")
            return
        if inp < 0.0 or inp > 1.0:
            self.father.rootwindow.father.dispMsg(u"\u03B7_max value not valid.")
            return 
        self.father.etamax = abs(float(inp))
        inp = safeEval(self.etamin.text(),type='FI')
        if inp is None:
            self.father.rootwindow.father.dispMsg(u"\u03B7_min value not valid.")
            return
        if inp < 0.0 or inp > 1.0:
            self.father.rootwindow.father.dispMsg(u"\u03B7_min value not valid.")
            return
        self.father.etamin = abs(float(inp))
        self.closeEvent()

##############################################################################

class Quad2CzjzekWindow(TabFittingWindow):

    def __init__(self, father, oldMainWindow):
        self.CURRENTWINDOW = Quad1DeconvFrame
        self.PARAMFRAME = Quad2CzjzekParamFrame
        super(Quad2CzjzekWindow, self).__init__(father, oldMainWindow)

#################################################################################


class Quad2CzjzekParamFrame(AbstractParamFrame):

    Ioptions = ['3/2', '5/2', '7/2', '9/2']

    def __init__(self, parent, rootwindow, isMain=True):
        self.SINGLENAMES = ['bgrnd']
        self.MULTINAMES = ['d', 'pos', 'sigma','wq0', 'eta0', 'amp', 'lor', 'gauss']
        self.PARAMTEXT = {'bgrnd': 'Background', 'd': 'd parameter', 'pos': 'Position', 'sigma': 'Sigma', 'wq0': 'Wq0', 'eta0': 'Eta0', 'amp': 'Integral', 'lor': 'Lorentz', 'gauss': 'Gauss'}
        self.FITFUNC = quad2CzjzekmpFit
        self.wqsteps = 50
        self.etasteps = 10
        self.wqmax = 4
        self.wqmin = 0
        self.etamax = 1
        self.etamin = 0
        # Get full integral
        self.fullInt = np.sum(parent.getData1D()) * parent.sw() / float(len(parent.getData1D()))
        super(Quad2CzjzekParamFrame, self).__init__(parent, rootwindow, isMain)
        self.cheng = 15
        self.ticks = {'bgrnd': [], 'pos': [], 'd': [], 'sigma': [], 'wq0': [], 'eta0': [], 'amp': [], 'lor': [], 'gauss': []}
        self.entries = {'bgrnd': [],'method': [], 'pos': [], 'd': [], 'sigma': [], 'wq0': [], 'eta0': [],'amp': [], 'lor': [], 'gauss': [], 'method': [], 'cheng': [], 'I': [], 'wqgrid': [], 'etagrid': [], 'wqmax': [], 'mas': []}
        czjzekPrefButton = QtWidgets.QPushButton("Grid Settings")
        czjzekPrefButton.clicked.connect(self.createCzjzekPrefWindow)
        self.frame1.addWidget(czjzekPrefButton, 1, 1)

        self.optframe.addWidget(wc.QLabel("Cheng:"), 0, 0)
        self.entries['cheng'].append(QtWidgets.QSpinBox())
        self.entries['cheng'][-1].setAlignment(QtCore.Qt.AlignHCenter)
        self.entries['cheng'][-1].setValue(self.cheng)
        self.optframe.addWidget(self.entries['cheng'][-1], 1, 0)
        self.optframe.addWidget(wc.QLabel("I:"), 2, 0)
        self.entries['I'].append(QtWidgets.QComboBox())
        self.entries['I'][-1].addItems(self.Ioptions)
        self.entries['I'][-1].setCurrentIndex(0)
        self.optframe.addWidget(self.entries['I'][-1], 3, 0)
        self.entries['mas'].append(QtWidgets.QCheckBox('Spinning'))
        self.optframe.addWidget(self.entries['mas'][-1], 4, 0)
        loadLibButton = QtWidgets.QPushButton("Load Library")
        loadLibButton.clicked.connect(self.loadLib)
        self.optframe.addWidget(loadLibButton, 11, 0)
        self.extLibCheck = QtWidgets.QCheckBox("Ext. Library")
        self.extLibCheck.setEnabled(False)
        self.optframe.addWidget(self.extLibCheck, 12, 0)
        self.optframe.setColumnStretch(21, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(wc.QLabel("Bgrnd:"), 0, 0, 1, 2)
        self.ticks['bgrnd'].append(QtWidgets.QCheckBox(''))
        self.ticks['bgrnd'][-1].setChecked(True)
        self.frame2.addWidget(self.ticks['bgrnd'][-1], 1, 0)
        self.entries['bgrnd'].append(wc.FitQLineEdit(self, 'bgrnd', "0.0"))
        self.entries['method'].append(QtWidgets.QComboBox())
        self.frame2.addWidget(wc.QLabel("Type:"), 4, 1)
        self.entries['method'][0].addItems(['Normal', 'Extended'])
        self.entries['method'][0].currentIndexChanged.connect(self.changeType)
        self.frame2.addWidget(self.entries['method'][0], 5, 1)
        self.frame2.addWidget(self.entries['bgrnd'][-1], 1, 1)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtWidgets.QComboBox()
        self.numExp.addItems([str(x + 1) for x in range(self.FITNUM)])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        self.frame3.addWidget(wc.QLabel("d:"), 1, 0, 1, 2)
        if self.parent.viewSettings["ppm"]:
            axUnit = 'ppm'
        else:
            axUnit = ['Hz', 'kHz', 'MHz'][self.parent.getAxType()]
        self.frame3.addWidget(wc.QLabel("Pos [" + axUnit + "]:"), 1, 2, 1, 2)
        self.frame3.addWidget(wc.QLabel(u"\u03c3 [MHz]:"), 1, 4, 1, 2)
        self.frame3.addWidget(wc.QLabel(u"\u03BD<sub>Q</sub>0 [MHz]:"), 1, 6, 1, 2)
        self.frame3.addWidget(wc.QLabel(u"\u03B70:"), 1, 8, 1, 2)
        self.frame3.addWidget(wc.QLabel("Integral:"), 1, 10, 1, 2)
        self.frame3.addWidget(wc.QLabel("Lorentz [Hz]:"), 1, 12, 1, 2)
        self.frame3.addWidget(wc.QLabel("Gauss [Hz]:"), 1, 14, 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        for i in range(self.FITNUM):
            for j in range(len(self.MULTINAMES)):
                self.ticks[self.MULTINAMES[j]].append(QtWidgets.QCheckBox(''))
                self.frame3.addWidget(self.ticks[self.MULTINAMES[j]][i], i + 2, 2 * j)
                self.entries[self.MULTINAMES[j]].append(wc.FitQLineEdit(self, self.MULTINAMES[j]))
                self.frame3.addWidget(self.entries[self.MULTINAMES[j]][i], i + 2, 2 * j + 1)
        self.changeType(0)
        self.dispParams()

    def changeType(self,index):
        if index == 0:
            for i in range(self.FITNUM): 
                self.entries[self.MULTINAMES[3]][i].setEnabled(False)
                self.entries[self.MULTINAMES[4]][i].setEnabled(False)
                self.ticks[self.MULTINAMES[3]][i].setChecked(True)
                self.ticks[self.MULTINAMES[4]][i].setChecked(True)
                self.ticks[self.MULTINAMES[3]][i].setEnabled(False)
                self.ticks[self.MULTINAMES[4]][i].setEnabled(False)
        elif index == 1:
            for i in range(self.FITNUM): 
                self.entries[self.MULTINAMES[3]][i].setEnabled(True)
                self.entries[self.MULTINAMES[4]][i].setEnabled(True)
                self.ticks[self.MULTINAMES[3]][i].setEnabled(True)
                self.ticks[self.MULTINAMES[4]][i].setEnabled(True)

    def createCzjzekPrefWindow(self, *args):
        CzjzekPrefWindow(self)

    def defaultValues(self, inp):
        if not inp:
            return {'bgrnd': [0.0, True],
                    'pos': np.repeat([np.array([0.0, False], dtype=object)], self.FITNUM, axis=0),
                    'd': np.repeat([np.array([5.0, True], dtype=object)], self.FITNUM, axis=0),
                    'sigma': np.repeat([np.array([1.0, False], dtype=object)], self.FITNUM, axis=0),
                    'wq0': np.repeat([np.array([0.0, True], dtype=object)], self.FITNUM, axis=0),
                    'eta0': np.repeat([np.array([0.0, True], dtype=object)], self.FITNUM, axis=0),
                    'amp': np.repeat([np.array([self.fullInt, False], dtype=object)], self.FITNUM, axis=0),
                    'lor': np.repeat([np.array([10.0, False], dtype=object)], self.FITNUM, axis=0),
                    'gauss': np.repeat([np.array([0.0, True], dtype=object)], self.FITNUM, axis=0)}
        else:
            return inp

    def checkI(self, I):
        return I * 1.0 + 1.5

    def setGrid(self, *args):
        inp = safeEval(self.entries['wqmax'][-1].text())
        if inp is None:
            return False
        self.entries['wqmax'][-1].setText(str(float(inp)))
        return True

    def bincounting(self, x1, weight, length):
        weights = weight[np.logical_and(x1 >= 0, x1 < length)]
        x1 = x1[np.logical_and(x1 >= 0, x1 < length)]
        fid = np.fft.ifft(np.bincount(x1, weights, length))
        #fid[0] *= 1
        return fid

    def genLib(self, length, I, minWq, maxWq,minEta,maxEta, numWq, numEta, angleStuff, freq, sw, weight):
        wq_return, eta_return = np.meshgrid(np.linspace(minWq, maxWq, numWq), np.linspace(minEta, maxEta, numEta))
        wq = wq_return[..., None]
        eta = eta_return[..., None]
        v = -1 / (6.0 * freq) * wq**2 * (I * (I + 1) - 3.0 / 4) * (angleStuff[0] + angleStuff[1] * eta + angleStuff[2] * eta**2)
        mult = v / sw * length
        x1 = np.array(np.round(mult) + np.floor(length / 2), dtype=int)
        lib = np.apply_along_axis(self.bincounting, 2, x1, weight, length)
        lib = lib.reshape((np.prod(lib.shape[:-1]), lib.shape[-1]))
        wq_return = wq_return.flatten()
        eta_return = eta_return.flatten()        
        return lib, wq_return, eta_return

    def loadLib(self):
        dirName = self.rootwindow.father.loadFitLibDir()
        nameList = os.listdir(dirName)
        cq = []
        eta = []
        data = []
        for name in nameList:
            matchName = re.search("-(\d+\.\d+)-(\d+\.\d+)\.(fid|spe)$", name)
            if matchName:
                eta.append(float(matchName.group(1)))
                cq.append(float(matchName.group(2)))
                fullName = os.path.join(dirName, name)
                libData = io.autoLoad(fullName)
                if libData.ndim() is not 1:
                    self.rootwindow.father.dispMsg("A spectrum in the library is not a 1D spectrum.")
                    continue
                if not libData.spec[0]:
                    libData.fourier(0)
                libData.regrid([self.parent.xax()[0], self.parent.xax()[-1]], len(self.parent.xax()), 0)
                libData.fftshift(0)
                libData.fourier(0)
                data.append(np.real(libData.data[0]))
        cq = np.array(cq) * 1e6
        eta = np.array(eta)
        data = np.array(data)
        numWq = len(np.unique(cq))
        numEta = len(np.unique(eta))
        if len(cq) != numWq * numEta:
            self.rootwindow.father.dispMsg("Library to be loaded is not of a rectangular grid in Cq and eta.")
            return
        sortIndex = np.lexsort((cq, eta))
        self.cqLib = cq[sortIndex]#.reshape((numEta, numWq))
        self.etaLib = eta[sortIndex]#.reshape((numEta, numWq))
        self.lib = data[sortIndex]#.reshape((numEta, numWq, len(data[0])))
        self.extLibCheck.setEnabled(True)
        self.extLibCheck.setChecked(True)

    def getExtraParams(self, out):
        mas = self.entries['mas'][-1].isChecked()
        wqMax = self.wqmax
        wqMin = self.wqmin
        etaMax = self.etamax
        etaMin = self.etamin
        I = self.checkI(self.entries['I'][-1].currentIndex())
        numWq = self.wqsteps
        numEta = self.etasteps
        if mas:
            weight, angleStuff = simFunc.quad2MASsetAngleStuff(self.entries['cheng'][-1].value())
        else:
            weight, angleStuff = simFunc.quad2StaticsetAngleStuff(self.entries['cheng'][-1].value())
        if self.extLibCheck.isChecked():
            lib = self.lib
            wq = self.cqLib
            eta = self.etaLib
        else:
            lib, wq, eta = self.genLib(len(self.parent.xax()), I, wqMin * 1e6, wqMax * 1e6, etaMin, etaMax, numWq, numEta, angleStuff, self.parent.freq(), self.parent.sw(), weight)
        out['I'] = [I]
        out['lib'] = [lib]
        out['wq'] = [wq]
        out['eta'] = [eta]
        out['freq'] = [self.parent.freq()]
        out['method'] = [self.entries['method'][0].currentIndex()]
        return (out, [I, lib, wq, eta, self.parent.freq(),self.entries['method'][0].currentIndex()])

    def disp(self, params, num):
        out = params[num]
        for name in self.SINGLENAMES:
            inp = out[name][0]
            if isinstance(inp, tuple):
                inp = checkLinkTuple(inp)
                out[name][0] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
        numExp = len(out[self.MULTINAMES[0]])
        for i in range(numExp):
            for name in self.MULTINAMES:
                inp = out[name][i]
                if isinstance(inp, tuple):
                    inp = checkLinkTuple(inp)
                    out[name][i] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                if not np.isfinite(out[name][i]):
                    self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
                    return
        tmpx = self.parent.xax()
        bgrnd = out['bgrnd'][0]
        method = out['method'][0]
        outCurve = bgrnd * np.ones(len(tmpx))
        outCurvePart = []
        x = []
        for i in range(len(out['pos'])):
            x.append(tmpx)
            if method == 0:
                y = out['amp'][i] * simFunc.quad2CzjzektensorFunc(tmpx, out['sigma'][i], out['d'][i], out['pos'][i]/self.axMult, out['lor'][i], out['gauss'][i], out['wq'][0], out['eta'][0], out['lib'][0], self.parent.freq(), self.parent.sw())
            if method == 1:
                y = out['amp'][i] * simFunc.quad2CzjzektensorFunc(tmpx, out['sigma'][i], out['d'][i], out['pos'][i]/self.axMult, out['lor'][i], out['gauss'][i], out['wq'][0], out['eta'][0], out['lib'][0], self.parent.freq(), self.parent.sw(), out['wq0'][i] * 1e6, out['eta0'][i])
            outCurvePart.append(bgrnd + y)
            outCurve += y
        locList = self.getRedLocList()
        self.parent.fitDataList[locList] = [tmpx, outCurve, x, outCurvePart]
        self.parent.showFid()

#################################################################################


def quad2CzjzekmpFit(xax, data1D, guess, args, queue, minmethod, numfeval):
    #try:
    fitVal = scipy.optimize.minimize(lambda *param: np.sum((data1D - quad2CzjzekfitFunc(param, xax, args))**2), guess, method=minmethod, options = {'maxfev': numfeval})
    #except Exception:
    #    fitVal = None
    queue.put(fitVal)


def quad2CzjzekfitFunc(params, allX, args):
    params = params[0]
    specName = args[0]
    specSlices = args[1]
    allParam = []
    for length in specSlices:
        allParam.append(params[length])
    allStruc = args[3]
    allArgu = args[4]
    fullTestFunc = []
    for n in range(len(allX)):
        x = allX[n][-1]
        testFunc = np.zeros(len(x))
        param = allParam[n]
        numExp = args[2][n]
        struc = args[3][n]
        argu = args[4][n]
        sw = args[5][n][-1]
        axMult = args[6][n]
        parameters = {}
        parameters['I'] = argu[-1][0]
        parameters['lib'] = argu[-1][1]
        parameters['wq'] = argu[-1][2]#.flatten()
        parameters['eta'] = argu[-1][3]#.flatten()
        parameters['freq'] = argu[-1][4]
        parameters['method'] = argu[-1][5]
        for name in ['bgrnd']:
            if struc[name][0][0] == 1:
                parameters[name] = param[struc[name][0][1]]
            elif struc[name][0][0] == 0:
                parameters[name] = argu[struc[name][0][1]]
            else:
                altStruc = struc[name][0][1]
                if struc[altStruc[0]][altStruc[1]][0] == 1:
                    parameters[name] = altStruc[2] * allParam[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                elif struc[altStruc[0]][altStruc[1]][0] == 0:
                    parameters[name] = altStruc[2] * allArgu[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
        for i in range(numExp):
            for name in ['d', 'pos', 'sigma', 'amp', 'lor', 'gauss', 'wq0', 'eta0']:
                if struc[name][i][0] == 1:
                    parameters[name] = param[struc[name][i][1]]
                elif struc[name][i][0] == 0:
                    parameters[name] = argu[struc[name][i][1]]
                else:
                    altStruc = struc[name][i][1]
                    strucTarget = allStruc[altStruc[4]]
                    if strucTarget[altStruc[0]][altStruc[1]][0] == 1:
                        parameters[name] = altStruc[2] * allParam[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                    elif strucTarget[altStruc[0]][altStruc[1]][0] == 0:
                        parameters[name] = altStruc[2] * allArgu[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
            if parameters['method'] == 0:
                testFunc += parameters['amp'] * simFunc.quad2CzjzektensorFunc(x, parameters['sigma'], parameters['d'], parameters['pos']/axMult, parameters['lor'], parameters['gauss'], parameters['wq'], parameters['eta'], parameters['lib'], parameters['freq'], sw)
            elif parameters['method'] == 1:
                testFunc += parameters['amp'] * simFunc.quad2CzjzektensorFunc(x, parameters['sigma'], parameters['d'], parameters['pos']/axMult, parameters['lor'], parameters['gauss'], parameters['wq'], parameters['eta'], parameters['lib'], parameters['freq'], sw, parameters['wq0'] * 1e6, parameters['eta0'])
        testFunc += parameters['bgrnd']
        fullTestFunc = np.append(fullTestFunc, testFunc)
    return fullTestFunc

#################################################################################

class SIMPSONDeconvWindow(TabFittingWindow):

    def __init__(self, father, oldMainWindow):
        self.CURRENTWINDOW = SIMPSONDeconvFrame
        self.PARAMFRAME = SIMPSONDeconvParamFrame
        super(SIMPSONDeconvWindow, self).__init__(father, oldMainWindow)

#################################################################################


class SIMPSONDeconvFrame(FitPlotFrame):

    FITNUM = 10  # Maximum number of fits

#################################################################################


class SIMPSONDeconvParamFrame(AbstractParamFrame):

    def __init__(self, parent, rootwindow, isMain=True):
        self.SINGLENAMES = []
        self.MULTINAMES = []
        self.PARAMTEXT = {}
        self.numExp = QtWidgets.QComboBox()
        self.script = None
        self.txtOutput = ["", ""]
        self.FITFUNC = SIMPSONDeconvmpFit
        # Get full integral
        super(SIMPSONDeconvParamFrame, self).__init__(parent, rootwindow, isMain)
        resetButton = QtWidgets.QPushButton("Reset")
        resetButton.clicked.connect(self.reset)
        self.frame1.addWidget(resetButton, 1, 1)
        loadButton = QtWidgets.QPushButton("Load Script")
        loadButton.clicked.connect(self.loadScript)
        self.frame1.addWidget(loadButton, 2, 1)
        outputButton = QtWidgets.QPushButton("Output")
        outputButton.clicked.connect(self.txtOutputWindow)
        self.frame1.addWidget(outputButton, 3, 1)
        self.commandLine = wc.QLineEdit("simpson")
        self.frame1.addWidget(self.commandLine, 4, 1)
        self.numExp = QtWidgets.QComboBox()
        self.numExp.addItems([str(x + 1) for x in range(self.FITNUM)])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        self.labels = {}
        self.ticks = {}
        self.entries = {}
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.reset()

    def txtOutputWindow(self):
        TxtOutputWindow(self.rootwindow, self.txtOutput[0], self.txtOutput[1])

    def defaultValues(self, inp):
        if not inp:
            val = {}
            for name in self.MULTINAMES:
                if name is "amp":
                    num = 1.0
                else:
                    num = 0.0
                val[name] = np.repeat([np.array([num, False], dtype=object)], self.FITNUM, axis=0)
            return val
        else:
            return inp

    def reset(self):
        locList = self.getRedLocList()
        self.fitParamList[locList] = self.defaultValues(0)
        self.dispParams()

    def loadScript(self):
        fileName = self.rootwindow.father.loadSIMPSONScript()
        with open(fileName, "r") as myfile:
            inFile = myfile.read()
        matches = np.unique(re.findall("(@\w+)", inFile))
        self.script = inFile
        self.MULTINAMES = [e[1:] for e in matches]
        self.MULTINAMES.insert(0, "amp")
        self.MULTINAMES.extend(["lor", "gauss"])
        for n in self.PARAMTEXT.keys():
            self.labels[n][0].deleteLater()
            for i in range(self.FITNUM):
                self.ticks[n][i].deleteLater()
                self.entries[n][i].deleteLater()
        self.PARAMTEXT = {}
        self.labels = {}
        self.ticks = {}
        self.entries = {}
        for i in range(len(self.MULTINAMES)):
            name = self.MULTINAMES[i]
            self.PARAMTEXT[name] = name
            self.labels[name] = [wc.QLabel(name)]
            self.frame3.addWidget(self.labels[name][0], 1, 2*i, 1, 2)
            self.ticks[name] = []
            self.entries[name] = []
            for j in range(self.FITNUM):
                self.ticks[name].append(QtWidgets.QCheckBox(''))
                self.frame3.addWidget(self.ticks[name][j], 2+j, 2*i)
                self.entries[name].append(wc.FitQLineEdit(self, name))
                self.frame3.addWidget(self.entries[name][j], 2+j, 2*i+1)
        self.reset()

    def getExtraParams(self, out):
        out["nameList"] = [self.MULTINAMES]
        out["command"] = [self.commandLine.text()]
        out["script"] = [self.script]
        out["txtOutput"] = [self.txtOutput]
        out["spec"] = [self.parent.spec()]
        return (out, [out["nameList"][-1], out["command"][-1], out["script"][-1], out["txtOutput"][-1], out["spec"][-1]])

    def disp(self, params, num):
        out = params[num]
        if len(self.MULTINAMES) > 0: 
            numExp = len(out[self.MULTINAMES[0]])
            for i in range(numExp):
                for name in self.MULTINAMES:
                    inp = out[name][i]
                    if isinstance(inp, tuple):
                        inp = checkLinkTuple(inp)
                        out[name][i] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                    if not np.isfinite(out[name][i]):
                        self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
                        return
            tmpx = self.parent.xax()
            bgrnd = np.zeros(len(tmpx))
            outCurve = bgrnd
            outCurvePart = []
            x = []
            for i in range(numExp):
                x.append(tmpx)
                inputPar = {}
                for name in self.MULTINAMES:
                    inputPar[name] = out[name][i]
                y = SIMPSONRunScript(out["command"][0], out["script"][0], inputPar, tmpx, out["txtOutput"][0], out["spec"][0])
                if y is None:
                    self.rootwindow.father.dispMsg("Fitting: The script didn't output anything", 'red')
                    return
                outCurvePart.append(bgrnd + y)
                outCurve += y
            locList = self.getRedLocList()

            self.parent.fitDataList[locList] = [tmpx, outCurve, x, outCurvePart]
            self.parent.showFid()

##############################################################################


def SIMPSONDeconvmpFit(xax, data1D, guess, args, queue, minmethod, numfeval):
    try:
        fitVal = scipy.optimize.minimize(lambda *param: np.sum((data1D - SIMPSONDeconvfitFunc(param, xax, args))**2), guess, method=minmethod, options = {'maxfev': numfeval})
    except Exception:
        fitVal = None
    queue.put(fitVal)


def SIMPSONDeconvfitFunc(params, allX, args):
    params = params[0]
    specName = args[0]
    specSlices = args[1]
    allParam = []
    for length in specSlices:
        allParam.append(params[length])
    allStruc = args[3]
    allArgu = args[4]
    fullTestFunc = []
    for n in range(len(allX)):
        x = allX[n][-1]
        testFunc = np.zeros(len(x))
        param = allParam[n]
        numExp = args[2][n]
        struc = args[3][n]
        argu = args[4][n]
        sw = args[5][n][-1]
        axMult = args[6][n]
        nameList = argu[-1][0]
        command = argu[-1][1]
        script = argu[-1][2]
        txtOutput = argu[-1][3]
        spec = argu[-1][4]
        parameters = {}
        for i in range(numExp):
            for name in nameList:
                if struc[name][i][0] == 1:
                    parameters[name] = param[struc[name][i][1]]
                elif struc[name][i][0] == 0:
                    parameters[name] = argu[struc[name][i][1]]
                else:
                    altStruc = struc[name][i][1]
                    strucTarget = allStruc[altStruc[4]]
                    if strucTarget[altStruc[0]][altStruc[1]][0] == 1:
                        parameters[name] = altStruc[2] * allParam[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                    elif strucTarget[altStruc[0]][altStruc[1]][0] == 0:
                        parameters[name] = altStruc[2] * allArgu[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
            testFunc += SIMPSONRunScript(command, script, parameters, x, txtOutput, spec)
        fullTestFunc = np.append(fullTestFunc, testFunc)
    return fullTestFunc


def SIMPSONRunScript(command, script, parameters, xax, output=None, spec=True):
    if script is None:
        return None
    for elem in parameters.keys():
        script = script.replace('@' + elem, str(parameters[elem]))
    directory_name = tempfile.mkdtemp()
    inputFileName = "simpsonScript.in"
    fullPath = os.path.join(directory_name, inputFileName)
    with open(fullPath, "w") as text_file:
        text_file.write(script)
    process = subprocess.Popen(command + ' ' + fullPath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=directory_name)
    if output:
        output[0], output[1] = process.communicate()
    else:
        process.wait()
    fileList = os.listdir(directory_name)
    fileList.remove(inputFileName)
    if not fileList:
        shutil.rmtree(directory_name, ignore_errors=True)
        return None
    outputFileName = fileList[0]
    masterData = io.autoLoad(os.path.join(directory_name, outputFileName))
    masterData.noUndo = True
    masterData.apodize(parameters["lor"], parameters["gauss"], 0, 0, 0, 0, 0, 0)
    if masterData.spec[0] != spec:
        masterData.fourier(0)
    masterData.regrid([xax[0], xax[-1]], len(xax), 0)
    shutil.rmtree(directory_name, ignore_errors=True)
    return parameters["amp"] * np.real(masterData.getHyperData(0))

class TxtOutputWindow(wc.ToolWindows):

    NAME = "Script Output"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent, txt, errTxt):
        super(TxtOutputWindow, self).__init__(parent)
        self.cancelButton.hide()
        self.grid.addWidget(QtWidgets.QLabel("Output:"), 1, 0)
        self.valEntry = QtWidgets.QTextEdit()
        self.valEntry.setReadOnly(True)
        self.valEntry.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.valEntry.setText(txt.decode())
        self.grid.addWidget(self.valEntry, 2, 0)
        self.grid.addWidget(QtWidgets.QLabel("Errors:"), 3, 0)
        self.errEntry = QtWidgets.QTextEdit()
        self.errEntry.setReadOnly(True)
        self.errEntry.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.errEntry.setText(errTxt.decode())
        self.grid.addWidget(self.errEntry, 4, 0)
        self.resize(550, 700)


#################################################################################


class FunctionFitWindow(TabFittingWindow):

    def __init__(self, father, oldMainWindow):
        self.CURRENTWINDOW = FunctionFitFrame
        self.PARAMFRAME = FunctionFitParamFrame
        super(FunctionFitWindow, self).__init__(father, oldMainWindow)

#################################################################################


class FunctionFitFrame(FitPlotFrame):

    FITNUM = 10  # Maximum number of fits

#################################################################################


class FunctionFitParamFrame(AbstractParamFrame):

    def __init__(self, parent, rootwindow, isMain=True):
        self.SINGLENAMES = []
        self.MULTINAMES = []
        self.PARAMTEXT = {}
        self.numExp = QtWidgets.QComboBox()
        self.function = ""
        self.FITFUNC = functionmpFit
        super(FunctionFitParamFrame, self).__init__(parent, rootwindow, isMain)
        resetButton = QtWidgets.QPushButton("Reset")
        resetButton.clicked.connect(self.reset)
        self.frame1.addWidget(resetButton, 1, 1)
        functionButton = QtWidgets.QPushButton("Input Function")
        functionButton.clicked.connect(self.functionInput)
        self.frame1.addWidget(functionButton, 2, 1)
        self.numExp = QtWidgets.QComboBox()
        self.numExp.addItems([str(x + 1) for x in range(self.FITNUM)])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        self.labels = {}
        self.ticks = {}
        self.entries = {}
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.reset()

    def defaultValues(self, inp):
        if not inp:
            val = {}
            for name in self.MULTINAMES:
                val[name] = np.repeat([np.array([0.0, False], dtype=object)], self.FITNUM, axis=0)
            return val
        else:
            return inp

    def reset(self):
        locList = self.getRedLocList()
        self.fitParamList[locList] = self.defaultValues(0)
        self.dispParams()

    def functionInput(self):
        FunctionInputWindow(self, self.function)

    def functionInputSetup(self):
        matches = np.unique(re.findall("(@\w+)", self.function))
        self.MULTINAMES = [e[1:] for e in matches]
        for n in self.PARAMTEXT.keys():
            self.labels[n][0].deleteLater()
            for i in range(self.FITNUM):
                self.ticks[n][i].deleteLater()
                self.entries[n][i].deleteLater()
        self.PARAMTEXT = {}
        self.labels = {}
        self.ticks = {}
        self.entries = {}
        for i in range(len(self.MULTINAMES)):
            name = self.MULTINAMES[i]
            self.PARAMTEXT[name] = name
            self.labels[name] = [wc.QLabel(name)]
            self.frame3.addWidget(self.labels[name][0], 1, 2*i, 1, 2)
            self.ticks[name] = []
            self.entries[name] = []
            for j in range(self.FITNUM):
                self.ticks[name].append(QtWidgets.QCheckBox(''))
                self.frame3.addWidget(self.ticks[name][j], 2+j, 2*i)
                self.entries[name].append(wc.FitQLineEdit(self, name))
                self.frame3.addWidget(self.entries[name][j], 2+j, 2*i+1)
        self.reset()

    def getExtraParams(self, out):
        out["nameList"] = [self.MULTINAMES]
        out["function"] = [self.function]
        return (out, [out["nameList"][-1], out["function"][-1]])

    def disp(self, params, num):
        out = params[num]
        if len(self.MULTINAMES) > 0: #If there are elements
            numExp = len(out[self.MULTINAMES[0]])
            for i in range(numExp):
                for name in self.MULTINAMES:
                    inp = out[name][i]
                    if isinstance(inp, tuple):
                        inp = checkLinkTuple(inp)
                        out[name][i] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                    if not np.isfinite(out[name][i]):
                        self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
                        return
            tmpx = self.parent.xax()
            bgrnd = np.zeros(len(tmpx))
            outCurve = bgrnd
            outCurvePart = []
            x = []
            for i in range(numExp):
                x.append(tmpx)
                inputPar = {}
                for name in self.MULTINAMES:
                    inputPar[name] = out[name][i]
                y = functionRun(out["function"][0], inputPar, tmpx)
                if y is None:
                    self.rootwindow.father.dispMsg("Fitting: The script didn't output anything", 'red')
                    return
                outCurvePart.append(bgrnd + y)
                outCurve += y
            locList = self.getRedLocList()
            self.parent.fitDataList[locList] = [tmpx, outCurve, x, outCurvePart]
            self.parent.showFid()

##############################################################################


def functionmpFit(xax, data1D, guess, args, queue, minmethod, numfeval):
    try:
        fitVal = scipy.optimize.minimize(lambda *param: np.sum((data1D - functionfitFunc(param, xax, args))**2), guess, method=minmethod, options = {'maxfev': numfeval})
    except Exception:
        fitVal = None
    queue.put(fitVal)


def functionfitFunc(params, allX, args):
    params = params[0]
    specName = args[0]
    specSlices = args[1]
    allParam = []
    for length in specSlices:
        allParam.append(params[length])
    allStruc = args[3]
    allArgu = args[4]
    fullTestFunc = []
    for n in range(len(allX)):
        x = allX[n][-1]
        testFunc = np.zeros(len(x))
        param = allParam[n]
        numExp = args[2][n]
        struc = args[3][n]
        argu = args[4][n]
        sw = args[5][n][-1]
        axMult = args[6][n]
        nameList = argu[-1][0]
        function = argu[-1][1]
        parameters = {}
        for i in range(numExp):
            for name in nameList:
                if struc[name][i][0] == 1:
                    parameters[name] = param[struc[name][i][1]]
                elif struc[name][i][0] == 0:
                    parameters[name] = argu[struc[name][i][1]]
                else:
                    altStruc = struc[name][i][1]
                    strucTarget = allStruc[altStruc[4]]
                    if strucTarget[altStruc[0]][altStruc[1]][0] == 1:
                        parameters[name] = altStruc[2] * allParam[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                    elif strucTarget[altStruc[0]][altStruc[1]][0] == 0:
                        parameters[name] = altStruc[2] * allArgu[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
            testFunc += functionRun(function, parameters, x)
        fullTestFunc = np.append(fullTestFunc, testFunc)
    return fullTestFunc


def functionRun(function, parameters, xax):
    for elem in parameters.keys():
        function = function.replace('@' + elem, str(parameters[elem]))
    return safeEval(function, length=len(xax), x=xax)


class FunctionInputWindow(wc.ToolWindows):

    NAME = "Fitting Function"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent, txt):
        super(FunctionInputWindow, self).__init__(parent)
        self.grid.addWidget(QtWidgets.QLabel("Function:"), 1, 0)
        self.valEntry = QtWidgets.QTextEdit()
        self.valEntry.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.valEntry.setText(txt)
        self.grid.addWidget(self.valEntry, 2, 0)
        self.resize(550, 400)

    def applyFunc(self):
        self.father.function = self.valEntry.toPlainText()
        self.father.functionInputSetup()

    def closeEvent(self, *args):
        self.deleteLater()

#################################################################################


class FitContourFrame(CurrentContour):

    MARKER = ''
    LINESTYLE = '-'
    FITNUM = 10  # Standard number of fits

    def __init__(self, rootwindow, fig, canvas, current):
        self.data = current.data
        tmp = np.array(current.data.shape(), dtype=int)
        tmp = np.delete(tmp, self.fixAxes(current.axes))
        self.fitDataList = np.full(tmp, None, dtype=object)
        self.fitPickNumList = np.zeros(tmp, dtype=int)
        self.rootwindow = rootwindow
        super(FitContourFrame, self).__init__(rootwindow, fig, canvas, current.data, current)

    def getRedLocList(self):
        return tuple(np.delete(self.locList, self.axes))
        
    def setSlice(self, axes, locList):
        self.rootwindow.paramframe.checkInputs()
        self.pickWidth = False
        super(FitContourFrame, self).setSlice(axes, locList)
        self.rootwindow.paramframe.checkFitParamList(self.getRedLocList())
        self.rootwindow.paramframe.dispParams()
        self.rootwindow.paramframe.togglePick()

    def getData1D(self):
        return np.real(self.getDataType(self.data1D.getHyperData(0)))

    def showFid(self):
        extraX = []
        extraY = []
        extraZ = []
        self.locList = np.array(self.locList, dtype=int)
        if self.fitDataList[self.getRedLocList()] is not None:
            tmp = self.fitDataList[self.getRedLocList()]
            extraX.append(tmp[0][1])
            extraY.append(tmp[0][0])
            extraZ.append(tmp[1])
            for i in range(len(tmp[2])):
                extraX.append(tmp[2][i][1])
                extraY.append(tmp[2][i][0])
                extraZ.append(tmp[3][i])
        super(FitContourFrame, self).showFid(extraX=extraX, extraY=extraY, extraZ=extraZ, extraColor=['C'+str(x%10) for x in range(len(extraX))])

        
##############################################################################


class MqmasDeconvWindow(TabFittingWindow):

    def __init__(self, father, oldMainWindow):
        self.CURRENTWINDOW = MqmasDeconvFrame
        self.PARAMFRAME = MqmasDeconvParamFrame
        super(MqmasDeconvWindow, self).__init__(father, oldMainWindow)

#################################################################################


class MqmasDeconvFrame(FitContourFrame):

    FITNUM = 10  # Maximum number of fits

#################################################################################


class MqmasDeconvParamFrame(AbstractParamFrame):

    Ioptions = ['3/2', '5/2', '7/2', '9/2']
    Ivalues = [1.5, 2.5, 3.5, 4.5]
    MQvalues = [3, 5, 7, 9]

    def __init__(self, parent, rootwindow, isMain=True):
        self.SINGLENAMES = ['bgrnd']
        self.MULTINAMES = ['pos', 'cq', 'eta', 'amp', 'lor2', 'gauss2', 'lor1', 'gauss1']
        self.PARAMTEXT = {'bgrnd': 'Background', 'pos': 'Position', 'cq': 'Cq', 'eta': 'eta', 'amp': 'Integral', 'lor2': 'Lorentz 2', 'gauss2': 'Gauss 2', 'lor1': 'Lorentz 1', 'gauss1': 'Gauss 1'}
        self.FITFUNC = mqmasmpFit
        self.setAngleStuff = simFunc.mqmasAngleStuff
        self.tensorFunc = simFunc.mqmasFunc
        self.cheng = 15
        # Get full integral
        self.fullInt = np.sum(parent.getData1D()) * parent.sw() / float(parent.getData1D().shape[-1]) * parent.sw(-2) / float(parent.getData1D().shape[-2])
        super(MqmasDeconvParamFrame, self).__init__(parent, rootwindow, isMain)
        self.ticks = {'bgrnd': [], 'pos': [], 'cq': [], 'eta': [], 'amp': [], 'lor2': [], 'gauss2': [], 'lor1': [], 'gauss1': []}
        self.entries = {'bgrnd': [], 'pos': [], 'cq': [], 'eta': [], 'amp': [], 'lor2': [], 'gauss2': [], 'lor1': [], 'gauss1': [], 'cheng': [], 'I': [], 'MQ': [], 'shear': [], 'scale': []}
        self.optframe.addWidget(wc.QLabel("Cheng:"), 0, 0)
        self.entries['cheng'].append(QtWidgets.QSpinBox())
        self.entries['cheng'][-1].setAlignment(QtCore.Qt.AlignHCenter)
        self.entries['cheng'][-1].setValue(self.cheng)
        self.optframe.addWidget(self.entries['cheng'][-1], 1, 0)
        self.optframe.addWidget(wc.QLabel("I:"), 2, 0)
        self.entries['I'].append(QtWidgets.QComboBox())
        self.entries['I'][-1].addItems(self.Ioptions)
        self.entries['I'][-1].setCurrentIndex(1)
        self.optframe.addWidget(self.entries['I'][-1], 3, 0)
        self.optframe.addWidget(wc.QLabel("MQ:"), 4, 0)
        self.entries['MQ'].append(QtWidgets.QComboBox())
        self.entries['MQ'][-1].addItems([str(i) for i in self.MQvalues])
        self.entries['MQ'][-1].setCurrentIndex(0)
        self.optframe.addWidget(self.entries['MQ'][-1], 5, 0)

        self.optframe.addWidget(wc.QLabel("Shear:"), 6, 0)
        self.entries['shear'].append(wc.QLineEdit("0.0"))
        self.optframe.addWidget(self.entries['shear'][-1], 7, 0)

        self.optframe.addWidget(wc.QLabel("Scale sw:"), 8, 0)
        self.entries['scale'].append(wc.QLineEdit("1.0"))
        self.optframe.addWidget(self.entries['scale'][-1], 9, 0)
        
        self.optframe.setColumnStretch(10, 1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(wc.QLabel("Bgrnd:"), 2, 0, 1, 2)
        self.ticks['bgrnd'].append(QtWidgets.QCheckBox(''))
        self.frame2.addWidget(self.ticks['bgrnd'][-1], 3, 0)
        self.entries['bgrnd'].append(wc.QLineEdit("0.0"))
        self.frame2.addWidget(self.entries['bgrnd'][-1], 3, 1)
        self.frame2.setColumnStretch(10, 1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtWidgets.QComboBox()
        self.numExp.addItems([str(x + 1) for x in range(self.FITNUM)])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp, 0, 0, 1, 2)
        if self.parent.viewSettings["ppm"][-1]:
            axUnit = 'ppm'
        else:
            axUnit = ['Hz', 'kHz', 'MHz'][self.parent.getAxType()]
        # Labels
        self.labelpos = wc.QLabel(u'Position [' + axUnit + ']:')
        self.labelcq = wc.QLabel(u'C<sub>Q</sub> [MHz]:')
        self.labeleta = wc.QLabel(u'\u03B7:')
        self.frame3.addWidget(self.labelpos, 1, 0, 1, 2)
        self.frame3.addWidget(self.labelcq, 1, 2, 1, 2)
        self.frame3.addWidget(self.labeleta, 1, 4, 1, 2)
        self.frame3.addWidget(wc.QLabel("Integral:"), 1, 6, 1, 2)
        self.frame3.addWidget(wc.QLabel("Lorentz 2 [Hz]:"), 1, 8, 1, 2)
        self.frame3.addWidget(wc.QLabel("Gauss 2 [Hz]:"), 1, 10, 1, 2)
        self.frame3.addWidget(wc.QLabel("Lorentz 1 [Hz]:"), 1, 12, 1, 2)
        self.frame3.addWidget(wc.QLabel("Gauss 1 [Hz]:"), 1, 14, 1, 2)
        self.frame3.setColumnStretch(20, 1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        for i in range(self.FITNUM):
            for j in range(len(self.MULTINAMES)):
                self.ticks[self.MULTINAMES[j]].append(QtWidgets.QCheckBox(''))
                self.frame3.addWidget(self.ticks[self.MULTINAMES[j]][i], i + 2, 2 * j)
                self.entries[self.MULTINAMES[j]].append(wc.FitQLineEdit(self, self.MULTINAMES[j]))
                self.frame3.addWidget(self.entries[self.MULTINAMES[j]][i], i + 2, 2 * j + 1)
        self.dispParams()

    def defaultValues(self, inp):
        if not inp:
            return {'bgrnd': [0.0, True],
                    'pos': np.repeat([np.array([0.0, False], dtype=object)], self.FITNUM, axis=0),
                    'cq': np.repeat([np.array([1.0, False], dtype=object)], self.FITNUM, axis=0),
                    'eta': np.repeat([np.array([0.0, False], dtype=object)], self.FITNUM, axis=0),
                    'amp': np.repeat([np.array([self.fullInt, False], dtype=object)], self.FITNUM, axis=0),
                    'lor2': np.repeat([np.array([1.0, False], dtype=object)], self.FITNUM, axis=0),
                    'gauss2': np.repeat([np.array([0.0, True], dtype=object)], self.FITNUM, axis=0),
                    'lor1': np.repeat([np.array([1.0, False], dtype=object)], self.FITNUM, axis=0),
                    'gauss1': np.repeat([np.array([0.0, True], dtype=object)], self.FITNUM, axis=0)}
        else:
            return inp

    def checkI(self, I):
        return I + 3/2.0

    def getExtraParams(self, out):
        if self.entries['MQ'][-1].currentIndex() > self.checkI(self.entries['I'][-1].currentIndex()):
            raise RuntimeError("MQ cannot be larger than I")
        out['I'] = [self.Ivalues[self.entries['I'][-1].currentIndex()]]
        out['MQ'] = [self.MQvalues[self.entries['MQ'][-1].currentIndex()]]
        out['shear'] = [safeEval(self.entries['shear'][-1].text())]
        out['scale'] = [safeEval(self.entries['scale'][-1].text())]
        cheng = safeEval(self.entries['cheng'][-1].text())
        out['cheng'] = [cheng]
        weight, angleStuff = self.setAngleStuff(cheng)
        out['weight'] = [weight]
        out['anglestuff'] = [angleStuff]
        out['tensorfunc'] = [self.tensorFunc]
        out['freq'] = [[self.parent.freq(-2), self.parent.freq()]]
        return (out, [out['I'][-1], out['MQ'][-1], out['shear'][-1], out['scale'][-1], out['weight'][-1], out['anglestuff'][-1], out['tensorfunc'][-1], out['freq'][-1]])

    def checkParam(self):
         val = self.numExp.currentIndex() + 1
         printStr = '%#.' + str(self.rootwindow.tabWindow.PRECIS) + 'g'
         for i in range(10):  # Print output if not stopped before
            if i < val:
                try:
                    cq = float(safeEval(self.entries['cq'][i].text()))
                    eta = float(safeEval(self.entries['eta'][i].text()))
                    Res = func.quadConversion([cq,eta], 1, 0)[0]
                    self.entries['cq'][i].setText(printStr % Res[0])
                    self.entries['eta'][i].setText(printStr % Res[1])
                except Exception:
                    return

    def disp(self, params, num):
        out = params[num]
        for name in self.SINGLENAMES:
            inp = out[name][0]
            if isinstance(inp, tuple):
                inp = checkLinkTuple(inp)
                out[name][0] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
        numExp = len(out[self.MULTINAMES[0]])
        for i in range(numExp):
            for name in self.MULTINAMES:
                inp = out[name][i]
                if isinstance(inp, tuple):
                    inp = checkLinkTuple(inp)
                    out[name][i] = inp[2] * params[inp[4]][inp[0]][inp[1]] + inp[3]
                if not np.isfinite(out[name][i]):
                    self.rootwindow.father.dispMsg("Fitting: One of the inputs is not valid")
                    return
        tmpx = [self.parent.xax(-2), self.parent.xax()]
        bgrnd = out['bgrnd'][0]
        outCurve = bgrnd
        outCurvePart = []
        x = []
        for i in range(len(out['amp'])):
            x.append(tmpx)
            y = out['amp'][i] * self.tensorFunc(tmpx, out['I'][0], out['MQ'][0], out['shear'][0], out['scale'][0], out['pos'][i]/self.axMult, out['cq'][i], out['eta'][i], [out['lor1'][i], out['lor2'][i]], [out['gauss1'][i], out['gauss2'][i]], out['anglestuff'][0], [self.parent.freq(-2), self.parent.freq()], [self.parent.sw(-2), self.parent.sw()], out['weight'][0])
            y = np.real(np.fft.fftshift(np.fft.fft(y, axis=0), axes=0))
            outCurvePart.append(bgrnd + y)
            outCurve += y
        locList = self.getRedLocList()
        self.parent.fitDataList[locList] = [tmpx, outCurve, x, outCurvePart]
        self.parent.showFid()
        self.checkParam()

##############################################################################


def mqmasmpFit(xax, data1D, guess, args, queue, minmethod, numfeval):
    try:
        fitVal = scipy.optimize.minimize(lambda *param: np.sum((data1D - mqmasfitFunc(param, xax, args))**2), guess, method=minmethod, options = {'maxfev': numfeval})
    except Exception:
        fitVal = None
    queue.put(fitVal)


def mqmasfitFunc(params, allX, args):
    params = params[0]
    specName = args[0]
    specSlices = args[1]
    allParam = []
    for length in specSlices:
        allParam.append(params[length])
    allStruc = args[3]
    allArgu = args[4]
    fullTestFunc = []
    for n in range(len(allX)):
        x = allX[n]
        testFunc = np.zeros([len(i) for i in x], dtype=complex)
        param = allParam[n]
        numExp = args[2][n]
        struc = args[3][n]
        argu = args[4][n]
        sw = args[5][n]
        axMult = args[6][n]
        parameters = {}
        parameters['I'] = argu[-1][0]
        parameters['MQ'] = argu[-1][1]
        parameters['shear'] = argu[-1][2]
        parameters['scale'] = argu[-1][3]
        parameters['weight'] = argu[-1][4]
        parameters['anglestuff'] = argu[-1][5]
        parameters['tensorfunc'] = argu[-1][6]
        parameters['freq'] = argu[-1][7]
        for name in ['bgrnd']:
            if struc[name][0][0] == 1:
                parameters[name] = param[struc[name][0][1]]
            elif struc[name][0][0] == 0:
                parameters[name] = argu[struc[name][0][1]]
            else:
                altStruc = struc[name][0][1]
                if struc[altStruc[0]][altStruc[1]][0] == 1:
                    parameters[name] = altStruc[2] * allParam[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                elif struc[altStruc[0]][altStruc[1]][0] == 0:
                    parameters[name] = altStruc[2] * allArgu[altStruc[4]][struc[altStruc[0]][altStruc[1]][1]] + altStruc[3]
        for i in range(numExp):
            for name in ['pos', 'cq', 'eta', 'amp', 'lor2', 'gauss2', 'lor1', 'gauss1']:
                if struc[name][i][0] == 1:
                    parameters[name] = param[struc[name][i][1]]
                elif struc[name][i][0] == 0:
                    parameters[name] = argu[struc[name][i][1]]
                else:
                    altStruc = struc[name][i][1]
                    strucTarget = allStruc[altStruc[4]]
                    if strucTarget[altStruc[0]][altStruc[1]][0] == 1:
                        parameters[name] = altStruc[2] * allParam[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
                    elif strucTarget[altStruc[0]][altStruc[1]][0] == 0:
                        parameters[name] = altStruc[2] * allArgu[altStruc[4]][strucTarget[altStruc[0]][altStruc[1]][1]] + altStruc[3]
            testFunc += parameters['amp'] * parameters['tensorfunc'](x, parameters['I'], parameters['MQ'], parameters['shear'], parameters['scale'], parameters['pos']/axMult, parameters['cq'], parameters['eta'], [parameters['lor1'], parameters['lor2']], [parameters['gauss1'], parameters['gauss2']], parameters['anglestuff'], parameters['freq'], sw, parameters['weight'])
        testFunc = np.real(np.fft.fftshift(np.fft.fft(testFunc, axis=0), axes=0))
        testFunc += parameters['bgrnd']
        fullTestFunc.append(testFunc)
    return np.array(fullTestFunc)
