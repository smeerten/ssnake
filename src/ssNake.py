#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

EXE = False

import sys
import os
import importlib
from PyQt5 import QtGui, QtCore, QtWidgets
QT = 5

QtCore.pyqtRemoveInputHook()
import matplotlib
# First import matplotlib and Qt
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import multiprocessing

# Create splash window
if __name__ == '__main__':
    multiprocessing.freeze_support() #Fix multiprocessing for pyinstaller on windows (line does nothing otherwise)
    root = QtWidgets.QApplication(sys.argv)
    root.setWindowIcon(QtGui.QIcon(os.path.dirname(os.path.realpath(__file__)) + '/Icons/logo.gif'))
    splash_pix = QtGui.QPixmap(os.path.dirname(os.path.realpath(__file__)) + '/Icons/logo.gif')
    splash = QtWidgets.QSplashScreen(splash_pix, QtCore.Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    progressBar = QtWidgets.QProgressBar(splash)
    progressBar.setGeometry(int(2.5 * splash.width() / 10.0), int(0.89 * splash.height()), int(5 * splash.width() / 10.0), int(splash.height() / 20.0))
    splash.show()

def splashProgressStep(splashStep):  # A function to easily increase the progressbar value
    if __name__ == '__main__':
        splashStep = splashStep + 1
        progressBar.setValue(int(splashStep // splashSteps + (splashStep % splashSteps > 0)))  # Rounds up without math or numpy module
        root.processEvents()
    return splashStep

def import_lib(name, nameAs, className, splashStep):
    # Function to load a library from string names
    if className is None:
        globals()[nameAs] = importlib.import_module(name)
    else:
        mod = importlib.import_module(name)
        globals()[nameAs] = getattr(mod, className)
    return splashProgressStep(splashStep)

# List of all libs to be imported:
# [name, name to be saved as, import specific class]
importList = [['matplotlib.figure', 'Figure', 'Figure'],
              ['traceback', 'tb', None],
              ['numpy', 'np', None],
              ['copy', 'copy', None],
              ['gc', 'gc', None],
              ['multiprocessing', 'multiprocessing', None],
              ['datetime', 'datetime', None],
              ['webbrowser', 'webbrowser', None],
              ['spectrum', 'sc', None],
              ['hypercomplex', 'hc', None],
              ['fitting', 'fit', None],
              ['safeEval', 'safeEval', 'safeEval'],
              ['widgetClasses', 'wc', None],
              ['updateWindow', 'UpdateWindow', 'UpdateWindow'],
              ['saveFigure', 'SaveFigureWindow', 'SaveFigureWindow'],
              ['functions', 'func', None],
              ['specIO', 'io', None],
              ['views', 'views', None],
              ['simFunctions', 'sim', None],
              ['loadIsotopes', 'loadIsotopes', None],
              ['scipy', 'optimize', 'optimize']]

splashSteps = len(importList) / 100.0
splashStep = 0
# Import everything else
for elem in importList:
    splashStep = import_lib(elem[0], elem[1], elem[2], splashStep)
isoPath = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + "IsotopeProperties"
ISOTOPES = loadIsotopes.getIsotopeInfo(isoPath)
matplotlib.rc('font', family='DejaVu Sans')
np.set_printoptions(threshold=sys.maxsize)
QtCore.QLocale.setDefault(QtCore.QLocale('en_US'))

VERSION = 'v1.5'
# Required library version
NPVERSION = '1.11.0'
MPLVERSION = '1.5.0'
SPVERSION = '0.14.1'
PY2VERSION = '2.7'
PY3VERSION = '3.4'

def splitString(val, size):
    #Split a string at spaces, with maxsize 'size'
    total = ''
    while len(val) > size:
        tmp = val[0:size]
        space = tmp.rfind(' ')
        total = total + val[0:space] + '\n'
        val = val[space + 1::] #increment string. +1 to cut off the space
    total = total + val
    return total

#Prepare TOOLTIPS dictionary
TOOLTIPS = dict()
with open(os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Tooltips', 'r') as f:
    data = f.read().split('\n')
for line in data:
    if line:
        tmp = line.split('\t')
        if len(tmp) > 2:
            text = tmp[1] + '\n\n' + splitString(tmp[2], 80) #Header plus split main text
        else:
            text = tmp[1]
        TOOLTIPS[tmp[0]] = text


class SsnakeException(sc.SpectrumException):
    pass


class MainProgram(QtWidgets.QMainWindow):

    def __init__(self, root):
        super(MainProgram, self).__init__()
        self.root = root
        self.VERSION = VERSION
        self.errors = []
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setAcceptDrops(True)
        self.mainWindow = None
        self.workspaces = []
        self.workspaceNames = []
        self.workspaceNum = 0
        self.macros = {}
        self.macroActions = {}
        self.referenceName = []  # List with saved reference names
        self.referenceValue = []  # List with saved reference values
        self.referenceActions = {}
        self.loadDefaults()
        if self.defaultStartupBool:
            self.lastLocation = os.path.expanduser(self.defaultStartupDir)
        else:
            if EXE:
                self.lastLocation = os.path.expanduser('~')
            else:
                self.lastLocation = os.getcwd()
        if not self.defaultTooltips: #Disable tooltips by setting them to empty strings
            for elem in TOOLTIPS.keys():
                TOOLTIPS[elem] = ''
        self.initMenu()
        self.menuCheck()
        self.main_widget = QtWidgets.QSplitter(self)
        self.main_widget.setHandleWidth(10)
        self.gridWidget = QtWidgets.QWidget(self)
        self.mainFrame = QtWidgets.QGridLayout(self.gridWidget)
        self.tree = wc.SsnakeTreeWidget(self)
        self.main_widget.addWidget(self.tree)
        self.main_widget.addWidget(self.gridWidget)
        self.logo = QtWidgets.QLabel(self)
        self.logo.setPixmap(QtGui.QPixmap(os.path.dirname(os.path.realpath(__file__)) + "/Icons/logo.gif"))
        self.mainFrame.addWidget(self.logo, 0, 0, QtCore.Qt.AlignCenter)
        self.tabs = wc.SsnakeTabs(self)
        self.tabs.setMovable(True)
        self.tabs.tabBar().tabMoved.connect(self.moveWorkspace)
        self.allowChange = True
        self.tabs.setTabsClosable(True)
        self.tabs.currentChanged.connect(self.changeMainWindow)
        self.tabs.tabCloseRequested.connect(self.destroyWorkspace)
        self.mainFrame.addWidget(self.tabs, 0, 0)
        self.statusBar = QtWidgets.QStatusBar(self)
        self.setStatusBar(self.statusBar)
        self.tabs.hide()
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        self.eventFilter = wc.SsnakeEventFilter(self)
        self.root.installEventFilter(self.eventFilter)
        self.initToolbar()
        self.main_widget.setStretchFactor(1, 10)
        #Set double click filter for splitter
        self.splitterEventFilter = wc.SplitterEventFilter(self.main_widget)
        self.main_widget.handle(1).installEventFilter(self.splitterEventFilter)
        self.resize(self.defaultWidth, self.defaultHeight)
        if self.defaultMaximized:
            self.showMaximized()
        QtWidgets.QShortcut(QtGui.QKeySequence.Paste, self).activated.connect(self.handlePaste)
        QtWidgets.QShortcut(QtGui.QKeySequence.Copy, self).activated.connect(self.handleCopy)

    def dispError(self, error):
        self.dispMsg("Program error. Please report.", color="red")
        CurTime = datetime.datetime.now()
        TimeStr = '{0:02d}'.format(CurTime.hour) + ':' + '{0:02d}'.format(CurTime.minute) + ':' + '{0:02d}'.format(CurTime.second)
        self.errors.append([TimeStr, error])

    def handlePaste(self):
        self.dropEvent(QtWidgets.QApplication.instance().clipboard())

    def handleCopy(self):
        """ Makes a pixelwise copy of the currently viewed canvas """
        if self.mainWindow is None:
            return
        if issubclass(type(self.mainWindow), fit.TabFittingWindow): #If fitting, take canvas from current tab
            canvas = self.mainWindow.tabs.currentWidget().canvas
        else:
            canvas = self.mainWindow.canvas
        screen = self.root.primaryScreen()
        pixmap = screen.grabWindow(canvas.winId())
        QtWidgets.QApplication.clipboard().setPixmap(pixmap)

    def resetDefaults(self):
        self.defaultUnits = 1
        self.defaultPrecis = 4
        self.defaultShowTitle = True
        self.defaultPPM = False
        self.defaultWidth = 1
        self.defaultHeight = 1
        self.defaultMaximized = False
        self.defaultAskName = True
        self.defaultToolBar = True
        self.defaultLinewidth = 1.0
        self.defaultMinXTicks = 12
        self.defaultMinYTicks = 8
        self.defaultColor = '#1F77B4'
        self.defaultGrids = [False, False]
        self.defaultDiagonalBool = False
        self.defaultDiagonalMult = 1
        self.defaultZeroScroll = True
        self.defaultZoomStep = 1
        self.defaultColorRange = 'none'
        self.defaultColorMap = 'seismic'
        self.defaultPColorMap = 'gray'
        self.defaultWidthRatio = 3.0
        self.defaultHeightRatio = 3.0
        self.defaultContourConst = True
        self.defaultPosColor = '#1F77B4'
        self.defaultNegColor = '#FF7F0E'
        self.defaultSecondOrderPhaseDialog = False
        self.defaultStartupBool = False
        self.defaultStartupDir = '~'
        self.defaultTooltips = True
        self.defaultToolbarActionList = ['File --> Open',
                                         'File -- > Save --> Matlab',
                                         'File --> Export --> Figure',
                                         'Seperator',
                                         'Workspaces --> Duplicate',
                                         'Workspaces --> Delete',
                                         'Seperator',
                                         'Edit --> Undo',
                                         'Edit --> Redo',
                                         'Edit --> Reload',
                                         'Seperator', 'Tools --> Apodize',
                                         'Tools --> Phasing --> Phase',
                                         'Tools --> Phasing --> Autophase 0',
                                         'Seperator',
                                         'Matrix --> Sizing',
                                         'Matrix --> Shift Data',
                                         'Matrix --> Multiply',
                                         'Seperator',
                                         'Fitting --> S/N',
                                         'Fitting --> FWHM',
                                         'Fitting --> Integrals',
                                         'Fitting --> Relaxation Curve',
                                         'Fitting --> Lorentzian/Gaussian',
                                         'Seperator',
                                         'Plot --> 1D Plot',
                                         'Plot --> Stack Plot',
                                         'Plot --> Array Plot',
                                         'Plot --> Contour Plot',
                                         'Plot --> Multi Plot',
                                         'Seperator',
                                         'History --> History',
                                         'History --> Clear Undo/Redo List',
                                         'Seperator',
                                         'Utilities --> NMR Table']

    def loadDefaults(self):
        self.resetDefaults()
        QtCore.QSettings.setDefaultFormat(QtCore.QSettings.IniFormat)
        QtCore.QCoreApplication.setOrganizationName("ssNake")
        QtCore.QCoreApplication.setApplicationName("ssNake")
        settings = QtCore.QSettings()
        try:
            self.defaultUnits = settings.value("plot/units", self.defaultUnits, int)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the units")
        self.defaultPPM = settings.value("plot/ppm", self.defaultPPM, bool)
        self.defaultToolbarActionList = settings.value("toolbarList", self.defaultToolbarActionList, str)
        self.defaultShowTitle = settings.value("plot/showTitle", self.defaultShowTitle, bool)
        self.defaultColor = settings.value("plot/colour", self.defaultColor, str)
        self.defaultColorRange = settings.value("plot/colourrange", self.defaultColorRange, str)
        if not str(self.defaultColorRange) in views.COLORRANGELIST:
            self.dispMsg("Incorrect colourrange in config file")
            self.defaultColorRange = views.COLORRANGELIST[0]
        try:
            self.defaultLinewidth = settings.value("plot/linewidth", self.defaultLinewidth, float)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the plot/linewidth")
        self.defaultMinXTicks = settings.value("plot/minXTicks", self.defaultMinXTicks, int)
        self.defaultMinYTicks = settings.value("plot/minYTicks", self.defaultMinYTicks, int)
        self.defaultGrids = [settings.value("plot/xgrid", self.defaultGrids[0], bool), settings.value("plot/ygrid", self.defaultGrids[1], bool)]
        self.defaultZeroScroll = settings.value("plot/zeroscroll", self.defaultZeroScroll, bool)
        self.defaultZoomStep = settings.value("plot/zoomstep", self.defaultZoomStep, float)
        self.defaultColorMap = settings.value("contour/colourmap", self.defaultColorMap, str)
        self.defaultContourConst = settings.value("contour/constantcolours", self.defaultContourConst, bool)
        self.defaultPosColor = settings.value("contour/poscolour", self.defaultPosColor, str)
        self.defaultNegColor = settings.value("contour/negcolour", self.defaultNegColor, str)
        self.defaultPColorMap = settings.value("2Dcolor/colourmap", self.defaultPColorMap, str)
        if not str(self.defaultColorMap) in views.COLORMAPLIST:
            self.dispMsg("Incorrect colourmap in config file")
            self.defaultColorMap = views.COLORMAPLIST[0]
        if not str(self.defaultPColorMap) in views.COLORMAPLIST:
            self.dispMsg("Incorrect pcolourmap in config file")
            self.defaultPColorMap = views.COLORMAPLIST[0]
        self.defaultDiagonalBool = settings.value("contour/diagonalbool", self.defaultDiagonalBool, bool)
        try:
            self.defaultDiagonalMult = settings.value("contour/diagonalmult", self.defaultDiagonalMult, float)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the diagonal multiplier")
        self.defaultPrecis = settings.value("precision", self.defaultPrecis, int)
        self.defaultMaximized = settings.value("maximized", self.defaultMaximized, bool)
        try:
            self.defaultWidth = settings.value("width", self.defaultWidth, int)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the width")
        try:
            self.defaultHeight = settings.value("height", self.defaultHeight, int)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the height")
        self.defaultAskName = settings.value("ask_name", self.defaultAskName, bool)
        self.defaultToolBar = settings.value("toolbar", self.defaultToolBar, bool)
        self.defaultStartupBool = settings.value("startupdiron", self.defaultStartupBool, bool)
        self.defaultStartupDir = settings.value("startupdir", self.defaultStartupDir, str)
        self.defaultTooltips = settings.value("tooltips", self.defaultTooltips, bool)
        try:
            self.defaultWidthRatio = settings.value("contour/width_ratio", self.defaultWidthRatio, float)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the contour/width_ratio")
        try:
            self.defaultHeightRatio = settings.value("contour/height_ratio", self.defaultHeightRatio, float)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the contour/height_ratio")
        self.defaultSecondOrderPhaseDialog = settings.value("phasing/second_order_phase_dialog", self.defaultSecondOrderPhaseDialog, bool)

    def saveDefaults(self):
        QtCore.QSettings.setDefaultFormat(QtCore.QSettings.IniFormat)
        QtCore.QCoreApplication.setOrganizationName("ssNake")
        QtCore.QCoreApplication.setApplicationName("ssNake")
        settings = QtCore.QSettings()
        settings.setValue("plot/showTitle", self.defaultShowTitle)
        settings.setValue("plot/units", self.defaultUnits)
        settings.setValue("plot/ppm", self.defaultPPM)
        settings.setValue('toolbarList', self.defaultToolbarActionList)
        settings.setValue("plot/colour", self.defaultColor)
        settings.setValue("plot/colourrange", self.defaultColorRange)
        settings.setValue("plot/linewidth", self.defaultLinewidth)
        settings.setValue("plot/minXTicks", self.defaultMinXTicks)
        settings.setValue("plot/minYTicks", self.defaultMinYTicks)
        settings.setValue("plot/xgrid", self.defaultGrids[0])
        settings.setValue("plot/ygrid", self.defaultGrids[1])
        settings.setValue("plot/zeroscroll", self.defaultZeroScroll)
        settings.setValue("plot/zoomstep", self.defaultZoomStep)
        settings.setValue("precision", self.defaultPrecis)
        settings.setValue("maximized", self.defaultMaximized)
        settings.setValue("width", self.defaultWidth)
        settings.setValue("height", self.defaultHeight)
        settings.setValue("ask_name", self.defaultAskName)
        settings.setValue("toolbar", self.defaultToolBar)
        settings.setValue("startupdiron", self.defaultStartupBool)
        settings.setValue("startupdir", self.defaultStartupDir)
        settings.setValue("tooltips", self.defaultTooltips)
        settings.setValue("contour/colourmap", self.defaultColorMap)
        settings.setValue("contour/constantcolours", self.defaultContourConst)
        settings.setValue("contour/poscolour", self.defaultPosColor)
        settings.setValue("contour/negcolour", self.defaultNegColor)
        settings.setValue("contour/width_ratio", self.defaultWidthRatio)
        settings.setValue("contour/height_ratio", self.defaultHeightRatio)
        settings.setValue("contour/diagonalbool", self.defaultDiagonalBool)
        settings.setValue("contour/diagonalmult", self.defaultDiagonalMult)
        settings.setValue("2Dcolor/colourmap", self.defaultPColorMap)
        settings.setValue("phasing/second_order_phase_dialog", self.defaultSecondOrderPhaseDialog)

    def dispMsg(self, msg, color='black'):
        if color == 'red':
            self.statusBar.setStyleSheet("QStatusBar{padding-left:8px;color:red;}")
        else:
            self.statusBar.setStyleSheet("QStatusBar{padding-left:8px;color:blue;}")
        self.statusBar.showMessage(msg, 10000)

    def initToolbar(self):
        if self.defaultToolBar:
            self.toolbar = self.addToolBar('Toolbar')
            self.toolbar.setMovable(False)
            self.toolbar.setIconSize(QtCore.QSize(22, 22))
            self.toolbar.toggleViewAction().setEnabled(False)
            self.seperatorAction = []
            self.allActionsList = [['Seperator', None],
                                   ['File --> Open', self.openAct],
                                   ['File --> Save --> JSON', self.saveAct],
                                   ['File -- > Save --> Matlab', self.saveMatAct],
                                   ['File --> Export --> Figure', self.savefigAct],
                                   ['File --> Export --> Simpson', self.saveSimpsonAct],
                                   ['File --> Export --> ASCII (1D/2D)', self.saveASCIIAct],
                                   ['File --> Export --> CSV (1D/2D)', self.saveCSVAct],
                                   ['File --> Preferences', self.preferencesAct],
                                   ['File --> Quit', self.quitAct],
                                   ['Workspaces --> Duplicate', self.newAct],
                                   ['Workspaces --> Slice to Workspace', self.newSlice],
                                   ['Workspaces --> Delete', self.closeAct],
                                   ['Workspaces --> Rename', self.renameWorkspaceAct],
                                   ['Workspaces --> Next', self.forwardAct],
                                   ['Workspaces --> Previous', self.backAct],
                                   ['Workspaces --> Info', self.workInfoAct],
                                   ['Macro --> Start Recording', self.macrostartAct],
                                   ['Macro --> Stop Recording', self.macrostopAct],
                                   ['Macro --> Load', self.macroLoadAct],
                                   ['Edit --> Undo', self.undoAction],
                                   ['Edit --> Redo', self.redoAction],
                                   ['Edit --> Reload', self.reloadAct],
                                   ['Edit --> Monitor', self.monitorAct],
                                   ['Tools --> Real', self.realAct],
                                   ['Tools --> Imag', self.imagAct],
                                   ['Tools --> Abs', self.absAct],
                                   ['Tools --> Complex Conjugate', self.conjAct],
                                   ['Tools --> Apodize', self.apodizeAct],
                                   ['Tools --> Phasing --> Phase', self.phaseAct],
                                   ['Tools --> Phasing --> Autophase 0', self.autoPhaseAct0],
                                   ['Tools --> Phasing --> Autophase 0+1', self.autoPhaseAct1],
                                   ['Tools --> Phasing --> Autophase per trace 0', self.autoPhaseAllAct0],
                                   ['Tools --> Phasing --> Autophase per trace 0+1', self.autoPhaseAllAct1],
                                   ['Tools --> Swap Echo', self.swapEchoAct],
                                   ['Tools --> Offset Correction', self.corOffsetAct],
                                   ['Tools --> Baseline Correction', self.baselineAct],
                                   ['Tools --> Subtract Averages', self.subAvgAct],
                                   ['Tools --> Reference Deconvolution', self.refDeconvAct],
                                   ['Tools --> Correct Digital Filter', self.digitalFilterAct],
                                   ['Tools --> Scale SW', self.scaleSWAct],
                                   ['Tools --> Scale freq and ref', self.scaleFreqRefAct],
                                   ['Tools --> LPSVD', self.lpsvdAct],
                                   ['Matrix --> Sizing', self.sizingAct],
                                   ['Matrix --> Shift Data', self.shiftAct],
                                   ['Matrix --> Roll Data', self.rollAct],
                                   ['Matrix --> Align Maxima', self.alignAct],
                                   ['Matrix --> Multiply', self.multiplyAct],
                                   ['Matrix --> Normalize', self.normalizeAct],
                                   ['Matrix --> Region --> Integrate', self.intRegionAct],
                                   ['Matrix --> Region --> Sum', self.sumRegionAct],
                                   ['Matrix --> Region --> Max', self.maxRegionAct],
                                   ['Matrix --> Region --> Min', self.minRegionAct],
                                   ['Matrix --> Region --> Max Position', self.maxposRegionAct],
                                   ['Matrix --> Region --> Min Position', self.minposRegionAct],
                                   ['Matrix --> Region --> Average', self.averageRegionAct],
                                   ['Matrix --> Diff', self.diffAct],
                                   ['Matrix --> Cumsum', self.cumsumAct],
                                   ['Matrix --> Extract Part', self.extractpartAct],
                                   ['Matrix --> Flip L/R', self.fliplrAct],
                                   ['Matrix --> Delete', self.matrixdelAct],
                                   ['Matrix --> Split', self.splitAct],
                                   ['Matrix --> Reorder', self.reorderAct],
                                   ['Matrix --> Regrid', self.regridAct],
                                   ['Matrix --> Concatenate', self.concatAct],
                                   ['Matrix --> Shearing', self.shearAct],
                                   ['Transforms --> Fourier Transform', self.fourierAct],
                                   ['Transforms --> Real Fourier Transform', self.realFourierAct],
                                   ['Transforms --> Fftshift', self.fftshiftAct],
                                   ['Transforms --> Inv fftshift', self.invfftshiftAct],
                                   ['Transforms --> Hilbert Transform', self.hilbertAct],
                                   ['Transforms --> NUS --> FFM', self.ffmAct],
                                   ['Transforms --> NUS --> CLEAN', self.cleanAct],
                                   ['Transforms --> NUS --> IST', self.istAct],
                                   ['Transforms --> Hypercomplex --> States', self.statesAct],
                                   ['Transforms --> Hypercomplex --> TPPI', self.statesTPPIAct],
                                   ['Transforms --> Hypercomplex --> Echo-antiecho', self.echoantiAct],
                                   ['Fitting --> S/N', self.snrAct],
                                   ['Fitting --> FWHM', self.fwhmAct],
                                   ['Fitting --> Centre of Mass', self.massAct],
                                   ['Fitting --> Integrals', self.intfitAct],
                                   ['Fitting --> Relaxation Curve', self.relaxAct],
                                   ['Fitting --> Diffusion Curve', self.diffusionAct],
                                   ['Fitting --> Lorentzian/Gaussian', self.lorentzfitAct],
                                   ['Fitting --> CSA', self.csastaticAct],
                                   ['Fitting --> Quadrupole', self.quadAct],
                                   ['Fitting --> Quadrupole+CSA', self.quadCSAAct],
                                   ['Fitting --> Czjzek', self.czjzekAct],
                                   ['Fitting --> MQMAS', self.mqmasAct],
                                   ['Fitting --> Czjzek MQMAS', self.mqmasCzjzekAct],
                                   ['Fitting --> External', self.externalFitAct],
                                   ['Fitting --> Function', self.functionFitAct],
                                   ['Combine --> Combine Workspaces', self.combineWorkspaceAct],
                                   ['Combine --> Insert From Workspace', self.insertdatAct],
                                   ['Combine --> Add', self.adddatAct],
                                   ['Combine --> Subtract', self.subdatAct],
                                   ['Combine --> Multiply', self.multdatAct],
                                   ['Combine --> Divide', self.divdatAct],
                                   ['Plot --> 1D Plot', self.onedplotAct],
                                   ['Plot --> Scatter', self.scatterplotAct],
                                   ['Plot --> Stack Plot', self.stackplotAct],
                                   ['Plot --> Array Plot', self.arrayplotAct],
                                   ['Plot --> Contour Plot', self.contourplotAct],
                                   ['Plot --> 2D Colour Plot', self.colour2DplotAct],
                                   ['Plot --> Multi Plot', self.multiplotAct],
                                   ['Plot --> Set Reference', self.setrefAct],
                                   ['Plot --> Clear Current Reference', self.delrefAct],
                                   ['Plot --> Load Reference', self.loadrefAct],
                                   ['Plot --> User X-axis', self.userxAct],
                                   ['Plot --> Plot Settings', self.plotprefAct],
                                   ['History --> History', self.historyAct],
                                   ['History --> Clear Undo/Redo List', self.clearundoAct],
                                   ['Utilities --> Chemical Shift Conversion Tool', self.shiftconvAct],
                                   ['Utilities --> Dipolar Distance Tool', self.dipolarconvAct],
                                   ['Utilities --> Quadrupole Coupling Conversion Tool', self.quadconvAct],
                                   ['Utilities --> NMR Table', self.nmrtableAct],
                                   ['Help --> GitLab Page', self.githubAct],
                                   ['Help --> ssNake Tutorials', self.tutorialAct],
                                   ['Help --> About', self.aboutAct]]
            for element in self.defaultToolbarActionList:
                if element == 'Seperator':
                    self.seperatorAction.append(QtWidgets.QAction(self))
                    self.seperatorAction[-1].setSeparator(True)
                    self.toolbar.addAction(self.seperatorAction[-1])
                else:
                    for action in self.allActionsList:
                        if element == action[0]:
                            self.toolbar.addAction(action[1])

    def initMenu(self):
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        self.menubar = self.menuBar()
        self.filemenu = QtWidgets.QMenu('&File', self)
        self.menubar.addMenu(self.filemenu)
        self.openAct = self.filemenu.addAction(QtGui.QIcon(IconDirectory + 'open.png'), '&Open', self.loadFromMenu, QtGui.QKeySequence.Open)
        self.openAct.setToolTip('Open a File')
        self.combineLoadAct = self.filemenu.addAction(QtGui.QIcon(IconDirectory + 'combine.png'), '&Open && Combine', self.createCombineLoadWindow)
        self.combineLoadAct.setToolTip('Open and Combine Multiple Files')
        self.savemenu = QtWidgets.QMenu('&Save', self)
        self.filemenu.addMenu(self.savemenu)
        self.saveAct = self.savemenu.addAction(QtGui.QIcon(IconDirectory + 'JSON.png'), 'JSON', self.saveJSONFile, QtGui.QKeySequence.Save)
        self.saveAct.setToolTip('Save as JSON File')
        self.saveMatAct = self.savemenu.addAction(QtGui.QIcon(IconDirectory + 'Matlab.png'), 'MATLAB', self.saveMatlabFile)
        self.saveMatAct.setToolTip('Save as MATLAB File')
        self.exportmenu = QtWidgets.QMenu('&Export', self)
        self.filemenu.addMenu(self.exportmenu)
        self.savefigAct = self.exportmenu.addAction(QtGui.QIcon(IconDirectory + 'figure.png'), 'Figure', self.saveFigure, QtGui.QKeySequence.Print)
        self.savefigAct.setToolTip('Export as Figure')
        self.saveSimpsonAct = self.exportmenu.addAction(QtGui.QIcon(IconDirectory + 'simpson.png'), 'Simpson', self.saveSimpsonFile)
        self.saveSimpsonAct.setToolTip('Export as Simpson File')
        self.saveASCIIAct = self.exportmenu.addAction(QtGui.QIcon(IconDirectory + 'ASCII.png'), 'ASCII (1D/2D)', self.saveASCIIFile)
        self.saveASCIIAct.setToolTip('Save as ASCII Text File')
        self.saveCSVAct = self.exportmenu.addAction(QtGui.QIcon(IconDirectory + 'CSV.png'),'CSV (1D/2D)', self.saveCSVFile)
        self.saveCSVAct.setToolTip('Save as CSV Text File')
        self.preferencesAct = self.filemenu.addAction(QtGui.QIcon(IconDirectory + 'preferences.png'), '&Preferences', lambda: PreferenceWindow(self))
        self.preferencesAct.setToolTip('Open Preferences Window')
        self.quitAct = self.filemenu.addAction(QtGui.QIcon(IconDirectory + 'quit.png'), '&Quit', self.fileQuit, QtGui.QKeySequence.Quit)
        self.quitAct.setToolTip('Close ssNake')
        self.saveActList = [self.saveAct, self.saveMatAct]
        self.exportActList = [self.savefigAct, self.saveSimpsonAct, self.saveASCIIAct,self.saveCSVAct]
        self.fileActList = [self.openAct, self.saveAct, self.saveMatAct,
                            self.savefigAct, self.saveSimpsonAct, self.saveASCIIAct,self.saveCSVAct,
                            self.combineLoadAct, self.preferencesAct, self.quitAct]
        # Workspaces menu
        self.workspacemenu = QtWidgets.QMenu('&Workspaces', self)
        self.menubar.addMenu(self.workspacemenu)
        self.newAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'duplicate.png'), 'D&uplicate', self.duplicateWorkspace, QtGui.QKeySequence.New)
        self.newAct.setToolTip('Duplicate Workspace')
        self.newSlice = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'duplicate.png'), 'Slice to Workspace', lambda: self.duplicateWorkspace(sliceOnly=True))
        self.newSlice.setToolTip('Copy Current Slice to New Workspace')
        self.closeAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'), '&Delete', self.destroyWorkspace, QtGui.QKeySequence.Close)
        self.closeAct.setToolTip('Delete Workspace')
        self.renameWorkspaceAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'), '&Rename', self.renameWorkspace, QtCore.Qt.Key_F2)
        self.renameWorkspaceAct.setToolTip('Rename Workspace')
        self.activemenu = QtWidgets.QMenu('&Go to', self)
        self.workspacemenu.addMenu(self.activemenu)
        self.forwardAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'next.png'), '&Next', lambda: self.stepWorkspace(1), QtGui.QKeySequence.Forward)
        self.forwardAct.setToolTip('Next Workspace')
        self.backAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'previous.png'), '&Previous', lambda: self.stepWorkspace(-1), QtGui.QKeySequence.Back)
        self.backAct.setToolTip('Previous Workspace')
        self.workInfoAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'about.png'), '&Info', lambda: self.mainWindowCheck(lambda mainWindow: WorkInfoWindow(mainWindow)))
        self.workInfoAct.setToolTip('Workspace Information')
        self.workspaceActList = [self.newAct, self.newSlice, self.closeAct, self.renameWorkspaceAct,
                                 self.forwardAct, self.backAct, self.workInfoAct]
        # Macro menu
        self.macromenu = QtWidgets.QMenu('&Macros', self)
        self.menubar.addMenu(self.macromenu)
        self.macrostartAct = self.macromenu.addAction(QtGui.QIcon(IconDirectory + 'record.png'), 'St&art Recording', self.macroCreate)
        self.macrostartAct.setToolTip('Start Recording Macro')
        self.macrostopAct = self.macromenu.addAction(QtGui.QIcon(IconDirectory + 'stop.png'), 'St&op Recording', self.stopMacro)
        self.macrostopAct.setToolTip('Stop Recording Macro')
        self.macrolistmenu = QtWidgets.QMenu('&Run', self)
        self.macromenu.addMenu(self.macrolistmenu)
        self.macrorenamemenu = QtWidgets.QMenu('Re&name', self)
        self.macromenu.addMenu(self.macrorenamemenu)
        self.macrodeletemenu = QtWidgets.QMenu('&Delete', self)
        self.macromenu.addMenu(self.macrodeletemenu)
        self.macrosavemenu = QtWidgets.QMenu('&Save', self)
        self.macromenu.addMenu(self.macrosavemenu)
        self.macroLoadAct = self.macromenu.addAction(QtGui.QIcon(IconDirectory + 'open.png'), '&Load', self.loadMacro)
        self.macroLoadAct.setToolTip('Load Macro')
        self.macroActList = [self.macrostartAct, self.macrostopAct]
        self.multiDActions = []
        # the edit drop down menu
        self.editmenu = QtWidgets.QMenu("&Edit", self)
        self.menubar.addMenu(self.editmenu)
        self.undoAction = self.editmenu.addAction(QtGui.QIcon(IconDirectory + 'undo.png'), "&Undo", self.undo, QtGui.QKeySequence.Undo)
        self.undoAction.setShortcutContext(QtCore.Qt.WidgetShortcut)
        self.undoAction.setToolTip('Undo')
        self.redoAction = self.editmenu.addAction(QtGui.QIcon(IconDirectory + 'redo.png'), "&Redo", self.redo, QtGui.QKeySequence.Redo)
        self.redoAction.setShortcutContext(QtCore.Qt.WidgetShortcut)
        self.redoAction.setToolTip('Redo')
        self.noUndoAct = QtWidgets.QAction("&No Undo Mode", self.editmenu, checkable=True)
        self.noUndoAct.toggled.connect(self.noUndoMode)
        self.editmenu.addAction(self.noUndoAct)
        self.clearundoAct = self.editmenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'), "&Clear Undo/Redo List", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.clearUndo()))
        self.clearundoAct.setToolTip('Clear Undo/Redo List')
        self.reloadAct = self.editmenu.addAction(QtGui.QIcon(IconDirectory + 'reload.png'), "Re&load", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.reloadLast()), QtGui.QKeySequence.Refresh)
        self.reloadAct.setToolTip('Reload Current Data')
        self.monitorAct = self.editmenu.addAction(QtGui.QIcon(IconDirectory + 'monitor.png'), "&Monitor", lambda: self.mainWindowCheck(lambda mainWindow: MonitorWindow(mainWindow)))
        self.monitorAct.setToolTip('Monitor Current Data')
        self.editActList = [self.undoAction, self.redoAction, self.clearundoAct, self.noUndoAct, self.reloadAct, self.monitorAct]
        # the tool drop down menu
        self.toolMenu = QtWidgets.QMenu("&Tools", self)
        self.menubar.addMenu(self.toolMenu)
        self.realAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'real.png'), "&Real", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.real()))
        self.realAct.setToolTip('Take Real Part of Data')
        self.imagAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'imag.png'), "&Imag", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.imag()))
        self.imagAct.setToolTip('Take Imaginary Part of Data')
        self.absAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'abs.png'), "&Abs", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.abs()))
        self.absAct.setToolTip('Take Absolute of Data')
        self.conjAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'complexconj.png'), "&Complex Conjugate", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.conj()))
        self.conjAct.setToolTip('Take Complex Conjugate of Data')
        self.apodizeAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'apodize.png'), "Apo&dize", lambda: self.mainWindowCheck(lambda mainWindow: ApodWindow(mainWindow)))
        self.apodizeAct.setToolTip('Open Apodize Window')
        self.phasingmenu = QtWidgets.QMenu('&Phasing', self)
        self.toolMenu.addMenu(self.phasingmenu)
        self.phaseAct = self.phasingmenu.addAction(QtGui.QIcon(IconDirectory + 'phase.png'), "&Phase", lambda: self.mainWindowCheck(lambda mainWindow: PhaseWindow(mainWindow)))
        self.phaseAct.setToolTip('Open Phasing Window')
        self.autoPhaseAct0 = self.phasingmenu.addAction(QtGui.QIcon(IconDirectory + 'autophase0.png'), "Autophase 0", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.directAutoPhase(0)))
        self.autoPhaseAct0.setToolTip('Autophase 0 order')
        self.autoPhaseAct1 = self.phasingmenu.addAction(QtGui.QIcon(IconDirectory + 'autophase1.png'), "Autophase 0+1", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.directAutoPhase(1)))
        self.autoPhaseAct1.setToolTip('Autophase 0 and 1 order')
        self.autoPhaseAllAct0 = self.phasingmenu.addAction(QtGui.QIcon(IconDirectory + 'autophase0.png'), "Autophase per trace 0", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.autoPhaseAll(0)))
        self.autoPhaseAllAct0.setToolTip('Autophase per trace 0 order')
        self.autoPhaseAllAct1 = self.phasingmenu.addAction(QtGui.QIcon(IconDirectory + 'autophase1.png'), "Autophase per trace 0+1", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.autoPhaseAll(1)))
        self.autoPhaseAllAct1.setToolTip('Autophase per trace 0 and 1 order')
        self.swapEchoAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'swapecho.png'), "Swap &Echo", lambda: self.mainWindowCheck(lambda mainWindow: SwapEchoWindow(mainWindow)))
        self.swapEchoAct.setToolTip('Swap Echo')
        self.corOffsetAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'offset.png'), "&Offset Correction", lambda: self.mainWindowCheck(lambda mainWindow: DCWindow(mainWindow)))
        self.corOffsetAct.setToolTip('Offset Correction')
        self.baselineAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'baseline.png'), "&Baseline Correction", lambda: self.mainWindowCheck(lambda mainWindow: BaselineWindow(mainWindow)))
        self.baselineAct.setToolTip('Baseline Correction')
        self.subAvgAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'subaverage.png'), "S&ubtract Averages", lambda: self.mainWindowCheck(lambda mainWindow: SubtractAvgWindow(mainWindow)))
        self.subAvgAct.setToolTip('Subtract Averages')
        self.refDeconvAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'deconvolute.png'), "Re&ference Deconvolution", lambda: self.mainWindowCheck(lambda mainWindow: FiddleWindow(mainWindow)))
        self.refDeconvAct.setToolTip('Reference Deconvolution')
        self.digitalFilterAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'dFilter.png'), "&Correct Digital Filter", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.CorrectDigitalFilter()))
        self.digitalFilterAct.setToolTip("Correct Digital Filter")
        self.lpsvdAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'LPSVD.png'), "&LPSVD", lambda: self.mainWindowCheck(lambda mainWindow: LPSVDWindow(mainWindow)))
        self.lpsvdAct.setToolTip('LPSVD linear prediction')
        self.scaleSWAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'ScaleSW.png'), "Scale SW", lambda: self.mainWindowCheck(lambda mainWindow: ScaleSWWindow(mainWindow)))
        self.scaleSWAct.setToolTip('Scale the Current Spectral Width')
        self.scaleFreqRefAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'ScaleFreqRef.png'), "Scale Car/Ref freq", 
                                    lambda: self.mainWindowCheck(lambda mainWindow: ScaleFreqRefWindow(mainWindow)))
        self.scaleFreqRefAct.setToolTip('Scale the current Carrier and Reference frequencies')
        self.referencelistmenu = QtWidgets.QMenu('&Reference', self)
        self.toolMenu.addMenu(self.referencelistmenu)
        self.setrefAct = self.referencelistmenu.addAction(QtGui.QIcon(IconDirectory + 'setreference.png'), "&Set Reference", lambda: self.mainWindowCheck(lambda mainWindow: RefWindow(mainWindow)))
        self.setrefAct.setToolTip('Set Reference')
        self.delrefAct = self.referencelistmenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'), "&Clear Current Reference", self.referenceClear)
        self.delrefAct.setToolTip('Clear Current Reference')
        self.referencerunmenu = QtWidgets.QMenu('&Apply', self)
        self.referencelistmenu.addMenu(self.referencerunmenu)
        self.referencedeletemenu = QtWidgets.QMenu('&Delete', self)
        self.referencelistmenu.addMenu(self.referencedeletemenu)
        self.referencerenamemenu = QtWidgets.QMenu('Re&name', self)
        self.referencelistmenu.addMenu(self.referencerenamemenu)
        self.referencesavemenu = QtWidgets.QMenu('&Save', self)
        self.referencelistmenu.addMenu(self.referencesavemenu)
        self.loadrefAct = self.referencelistmenu.addAction(QtGui.QIcon(IconDirectory + 'open.png'), "&Load", self.referenceLoad)
        self.loadrefAct.setToolTip('Load Reference')
        self.toolsActList = [self.realAct, self.imagAct, self.absAct, self.conjAct,
                             self.apodizeAct, self.phaseAct, self.autoPhaseAct0,
                             self.autoPhaseAct1, self.autoPhaseAllAct0, self.phasingmenu,
                             self.autoPhaseAllAct1, self.swapEchoAct, self.corOffsetAct,
                             self.baselineAct, self.subAvgAct, self.refDeconvAct, self.lpsvdAct,
                             self.digitalFilterAct, self.scaleSWAct, self.scaleFreqRefAct]
        # the matrix drop down menu
        self.matrixMenu = QtWidgets.QMenu("M&atrix", self)
        self.menubar.addMenu(self.matrixMenu)
        self.sizingAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'sizing.png'), "&Sizing", lambda: self.mainWindowCheck(lambda mainWindow: SizeWindow(mainWindow)))
        self.sizingAct.setToolTip('Set Size')
        self.shiftAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'shift.png'), "S&hift Data", lambda: self.mainWindowCheck(lambda mainWindow: ShiftDataWindow(mainWindow)))
        self.shiftAct.setToolTip('Shift Data')
        self.rollAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'roll.png'), "Roll Data", lambda: self.mainWindowCheck(lambda mainWindow: RollDataWindow(mainWindow)))
        self.rollAct.setToolTip('Roll Data')
        self.alignAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'alignMax.png'), "Align Maxima", lambda: self.mainWindowCheck(lambda mainWindow: AlignDataWindow(mainWindow)))
        self.alignAct.setToolTip('Align Maxima')
        self.regionMenu = QtWidgets.QMenu("Region", self)
        self.matrixMenu.addMenu(self.regionMenu)
        self.intRegionAct = self.regionMenu.addAction(QtGui.QIcon(IconDirectory + 'int.png'), "&Integrate", lambda: self.mainWindowCheck(lambda mainWindow: integrateWindow(mainWindow)))
        self.intRegionAct.setToolTip('Integrate Region')
        self.sumRegionAct = self.regionMenu.addAction(QtGui.QIcon(IconDirectory + 'sum.png'), "S&um", lambda: self.mainWindowCheck(lambda mainWindow: sumWindow(mainWindow)))
        self.sumRegionAct.setToolTip('Sum Region')
        self.maxRegionAct = self.regionMenu.addAction(QtGui.QIcon(IconDirectory + 'max.png'), "&Max", lambda: self.mainWindowCheck(lambda mainWindow: maxWindow(mainWindow)))
        self.maxRegionAct.setToolTip('Maximum of Region')
        self.minRegionAct = self.regionMenu.addAction(QtGui.QIcon(IconDirectory + 'min.png'), "M&in", lambda: self.mainWindowCheck(lambda mainWindow: minWindow(mainWindow)))
        self.minRegionAct.setToolTip('Minimum of Region')
        self.maxposRegionAct = self.regionMenu.addAction(QtGui.QIcon(IconDirectory + 'maxpos.png'), "Ma&x position", lambda: self.mainWindowCheck(lambda mainWindow: argmaxWindow(mainWindow)))
        self.maxposRegionAct.setToolTip('Position of Maximum of Region')
        self.minposRegionAct = self.regionMenu.addAction(QtGui.QIcon(IconDirectory + 'minpos.png'), "Mi&n position", lambda: self.mainWindowCheck(lambda mainWindow: argminWindow(mainWindow)))
        self.minposRegionAct.setToolTip('Position of Minimum of Region')
        self.averageRegionAct = self.regionMenu.addAction(QtGui.QIcon(IconDirectory + 'average.png'), "&Average", lambda: self.mainWindowCheck(lambda mainWindow: avgWindow(mainWindow)))
        self.averageRegionAct.setToolTip('Average of Region')
        self.diffAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'diff.png'), "&Diff", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.diff()))
        self.diffAct.setToolTip('Difference')
        self.cumsumAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'cumsum.png'), "&Cumsum", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.cumsum()))
        self.cumsumAct.setToolTip('Cumulative sum')
        self.extractpartAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'extractpart.png'), "&Extract part", lambda: self.mainWindowCheck(lambda mainWindow: extractRegionWindow(mainWindow)))
        self.extractpartAct.setToolTip('Extract part')
        self.fliplrAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'fliplr.png'), "&Flip L/R", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.flipLR()))
        self.fliplrAct.setToolTip('Flip L/R')
        self.matrixdelAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'matrixdelete.png'), "De&lete", lambda: self.mainWindowCheck(lambda mainWindow: DeleteWindow(mainWindow)))
        self.matrixdelAct.setToolTip('Delete Points')
        self.splitAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'split.png'), "S&plit", lambda: self.mainWindowCheck(lambda mainWindow: SplitWindow(mainWindow)))
        self.splitAct.setToolTip('Split')
        self.multiplyAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'multiply.png'), "Mul&tiply", lambda: self.mainWindowCheck(lambda mainWindow: MultiplyWindow(mainWindow)))
        self.multiplyAct.setToolTip('Multiply')
        self.normalizeAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'normalize.png'), "Normalize", lambda: self.mainWindowCheck(lambda mainWindow: NormalizeWindow(mainWindow)))
        self.normalizeAct.setToolTip('Normalize')
        self.reorderAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'reorder.png'), "&Reorder", lambda: self.mainWindowCheck(lambda mainWindow: ReorderWindow(mainWindow)))
        self.reorderAct.setToolTip('Reorder')
        self.regridAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'regrid.png'), "Regrid", lambda: self.mainWindowCheck(lambda mainWindow: RegridWindow(mainWindow)))
        self.regridAct.setToolTip('Regrid')
        self.concatAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'concatenate.png'), "C&oncatenate", lambda: self.mainWindowCheck(lambda mainWindow: ConcatenateWindow(mainWindow)))
        self.concatAct.setToolTip('Concatenate')
        self.multiDActions.append(self.concatAct)
        self.shearAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'shear.png'), "Shearin&g", lambda: self.mainWindowCheck(lambda mainWindow: ShearingWindow(mainWindow)))
        self.shearAct.setToolTip('Shearing')
        self.multiDActions.append(self.shearAct)
        self.matrixActList = [self.sizingAct, self.shiftAct, self.rollAct, self.alignAct, self.intRegionAct,
                              self.sumRegionAct, self.maxRegionAct, self.minRegionAct,
                              self.maxposRegionAct, self.minposRegionAct, self.averageRegionAct,
                              self.diffAct, self.cumsumAct, self.extractpartAct,
                              self.fliplrAct, self.matrixdelAct, self.splitAct,
                              self.multiplyAct, self.normalizeAct, self.reorderAct, self.regridAct,
                              self.concatAct, self.shearAct]
        # the Transforms drop down menu
        self.transformsMenu = QtWidgets.QMenu("T&ransforms", self)
        self.menubar.addMenu(self.transformsMenu)
        self.fourierAct = self.transformsMenu.addAction(QtGui.QIcon(IconDirectory + 'fourier.png'), "&Fourier Transform", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.fourier()), QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        self.fourierAct.setToolTip('Fourier Transform')
        self.realFourierAct = self.transformsMenu.addAction(QtGui.QIcon(IconDirectory + 'realfourier.png'), "&Real Fourier Transform", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.realFourier()))
        self.realFourierAct.setToolTip('Real Fourier Transform')
        self.fftshiftAct = self.transformsMenu.addAction(QtGui.QIcon(IconDirectory + 'fftshift.png'), "Fft&shift", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.fftshift()))
        self.fftshiftAct.setToolTip('Fftshift')
        self.invfftshiftAct = self.transformsMenu.addAction(QtGui.QIcon(IconDirectory + 'ifftshift.png'), "&Inv fftshift", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.invFftshift()))
        self.invfftshiftAct.setToolTip('Inverse fftshift')
        self.hilbertAct = self.transformsMenu.addAction(QtGui.QIcon(IconDirectory + 'hilbert.png'), "&Hilbert Transform", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.hilbert()))
        self.hilbertAct.setToolTip('Hilbert Transform')
        self.nusMenu = QtWidgets.QMenu("&NUS", self)
        self.transformsMenu.addMenu(self.nusMenu)
        self.ffmAct = self.nusMenu.addAction(QtGui.QIcon(IconDirectory + 'ffm.png'), "&FFM", lambda: self.mainWindowCheck(lambda mainWindow: FFMWindow(mainWindow)))
        self.ffmAct.setToolTip('FFM')
        self.cleanAct = self.nusMenu.addAction(QtGui.QIcon(IconDirectory + 'clean.png'), "&CLEAN", lambda: self.mainWindowCheck(lambda mainWindow: CLEANWindow(mainWindow)))
        self.cleanAct.setToolTip('CLEAN')
        self.istAct = self.nusMenu.addAction(QtGui.QIcon(IconDirectory + 'ist.png'), "&IST", lambda: self.mainWindowCheck(lambda mainWindow: ISTWindow(mainWindow)))
        self.istAct.setToolTip('IST')
        self.hypercomplexMenu = QtWidgets.QMenu("Hypercomplex", self)
        self.transformsMenu.addMenu(self.hypercomplexMenu)
        self.statesAct = self.hypercomplexMenu.addAction(QtGui.QIcon(IconDirectory + 'States.png'), "&States", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.states()))
        self.statesAct.setToolTip('States Hypercomplex Data Processing')
        self.statesTPPIAct = self.hypercomplexMenu.addAction(QtGui.QIcon(IconDirectory + 'statestppi.png'), "States-&TPPI", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.statesTPPI()))
        self.statesTPPIAct.setToolTip('States-TPPI Hypercomplex Data Processing')
        self.echoantiAct = self.hypercomplexMenu.addAction(QtGui.QIcon(IconDirectory + 'echoantiecho.png'), "Ec&ho-antiecho", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.echoAntiEcho()))
        self.echoantiAct.setToolTip('Ec&ho-antiecho Hypercomplex Data Processing')
        self.transformActList = [self.fourierAct, self.realFourierAct, self.fftshiftAct,
                                 self.invfftshiftAct, self.hilbertAct, self.ffmAct,
                                 self.cleanAct, self.istAct, self.statesAct, self.statesTPPIAct, self.echoantiAct]
        # the fitting drop down menu
        self.fittingMenu = QtWidgets.QMenu("F&itting", self)
        self.menubar.addMenu(self.fittingMenu)
        self.snrAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'snr.png'), "&S/N", lambda: self.mainWindowCheck(lambda mainWindow: SNWindow(mainWindow)))
        self.snrAct.setToolTip('Signal-to-Noise Ratio')
        self.fwhmAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'fwhm.png'), "&FWHM", lambda: self.mainWindowCheck(lambda mainWindow: FWHMWindow(mainWindow)))
        self.fwhmAct.setToolTip('Full Width at Half Maximum')
        self.massAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'mass.png'), "Centre of Mass", lambda: self.mainWindowCheck(lambda mainWindow: COMWindow(mainWindow)))
        self.massAct.setToolTip('Centre of Mass')
        self.intfitAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'int.png'), "&Integrals", lambda: self.mainWindowCheck(lambda mainWindow: IntegralsWindow(mainWindow)))
        self.intfitAct.setToolTip('Get Integrals')
        self.relaxAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'relaxation.png'), "&Relaxation Curve", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createRelaxWindow()))
        self.relaxAct.setToolTip('Fit Relaxation Curve')
        self.diffusionAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'diffusion.png'), "&Diffusion Curve", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createDiffusionWindow()))
        self.diffusionAct.setToolTip('Fit Diffusion Curve')
        self.lorentzfitAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'lorentz.png'), "&Lorentzian/Gaussian", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createPeakDeconvWindow()))
        self.lorentzfitAct.setToolTip('Fit Lorentzian/Gaussian')
        self.csastaticAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'csastatic.png'), "&CSA", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createCsaDeconvWindow()))
        self.csastaticAct.setToolTip('Fit CSA')
        self.quadAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'quadconversion.png'), "&Quadrupole", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuadDeconvWindow()))
        self.quadAct.setToolTip('Fit Quadrupole')
        self.quadCSAAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'quadcsa.png'), "Q&uadrupole+CSA", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuadCSADeconvWindow()))
        self.quadCSAAct.setToolTip('Fit Quadrupole+CSA')
        self.czjzekAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'czjzekstatic.png'), "C&zjzek", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuadCzjzekWindow()))
        self.czjzekAct.setToolTip('Fit Czjzek Pattern')
        self.mqmasAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'mqmas.png'), "&MQMAS", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createMQMASWindow()))
        self.mqmasAct.setToolTip('Fit MQMAS')
        self.mqmasCzjzekAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'mqmas.png'), "Cz&jzek MQMAS", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createMQMASCzjzekWindow()))
        self.mqmasCzjzekAct.setToolTip('Fit Czjzek MQMAS')
        self.externalFitAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'simpson.png'), "&External", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createExternalFitWindow()))
        self.externalFitAct.setToolTip('Fit External')
        self.functionFitAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'function.png'), "F&unction fit", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createFunctionFitWindow()))
        self.functionFitAct.setToolTip('Fit Function')
        self.fittingActList = [self.snrAct, self.fwhmAct, self.massAct,
                               self.intfitAct, self.relaxAct, self.diffusionAct,
                               self.lorentzfitAct, self.csastaticAct, self.quadAct, self.quadCSAAct,
                               self.czjzekAct, self.externalFitAct, self.functionFitAct]
        # the combine drop down menu
        self.combineMenu = QtWidgets.QMenu("Com&bine", self)
        self.menubar.addMenu(self.combineMenu)
        self.combineWorkspaceAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'combine.png'), '&Combine Workspaces', self.createCombineWorkspaceWindow)
        self.combineWorkspaceAct.setToolTip('Combine Workspaces')
        self.insertdatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'insert.png'), "&Insert From Workspace", lambda: self.mainWindowCheck(lambda mainWindow: InsertWindow(mainWindow)))
        self.insertdatAct.setToolTip('Insert From Workspace')
        self.adddatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'add.png'), "&Add", lambda: self.mainWindowCheck(lambda mainWindow: CombineWindow(mainWindow, 0)))
        self.adddatAct.setToolTip('Add Data From Workspace')
        self.subdatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'subtract.png'), "&Subtract", lambda: self.mainWindowCheck(lambda mainWindow: CombineWindow(mainWindow, 1)))
        self.subdatAct.setToolTip('Subtract Data From Workspace')
        self.multdatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'multiplyWorkspace.png'), "&Multiply", lambda: self.mainWindowCheck(lambda mainWindow: CombineWindow(mainWindow, 2)))
        self.multdatAct.setToolTip('Multiply Data From Workspace')
        self.divdatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'divideWorkspace.png'), "&Divide", lambda: self.mainWindowCheck(lambda mainWindow: CombineWindow(mainWindow, 3)))
        self.divdatAct.setToolTip('Divide Data From Workspace')
        self.combineActList = [self.combineWorkspaceAct, self.insertdatAct, self.adddatAct,
                               self.subdatAct, self.multdatAct, self.divdatAct]
        # the plot drop down menu
        self.plotMenu = QtWidgets.QMenu("&Plot", self)
        self.menubar.addMenu(self.plotMenu)
        self.onedplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + '1dplot.png'), "&1D Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plot1D()))
        self.onedplotAct.setToolTip('1D plot')
        self.scatterplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'scatterplot.png'), "&Scatter Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotScatter()))
        self.scatterplotAct.setToolTip('Scatter Plot')
        self.stackplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'stack.png'), "S&tack Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotStack()))
        self.stackplotAct.setToolTip('Stack Plot')
        self.multiDActions.append(self.stackplotAct)
        self.arrayplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'array.png'), "&Array Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotArray()))
        self.arrayplotAct.setToolTip('Array Plot')
        self.multiDActions.append(self.arrayplotAct)
        self.contourplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'contour.png'), "&Contour Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotContour()))
        self.contourplotAct.setToolTip('Contour Plot')
        self.multiDActions.append(self.contourplotAct)
        self.multiContourplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'multicontour.png'),"Mu&lti Contour Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotMultiContour()))
        self.multiContourplotAct.setToolTip('Multi Contour Plot')
        self.multiDActions.append(self.multiContourplotAct)
        self.colour2DplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + '2DColour.png'), "2D Colour Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotColour2D()))
        self.colour2DplotAct.setToolTip('2D Colour Plot')
        self.multiDActions.append(self.colour2DplotAct)
        self.multiplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'multi.png'), "&Multi Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotMulti()))
        self.multiplotAct.setToolTip('Multi Plot')
        #==========
        self.userxAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'xaxis.png'), "&User X-axis", lambda: self.mainWindowCheck(lambda mainWindow: XaxWindow(mainWindow)))
        self.userxAct.setToolTip('User X-axis')
        self.plotprefAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'preferences.png'), "&Plot Settings", lambda: self.mainWindowCheck(lambda mainWindow: PlotSettingsWindow(mainWindow)))
        self.plotprefAct.setToolTip('Plot Settings')
        self.plotActList = [self.onedplotAct, self.scatterplotAct, self.stackplotAct,
                            self.arrayplotAct, self.contourplotAct, self.colour2DplotAct, self.multiplotAct,
                            self.multiContourplotAct, self.setrefAct, self.delrefAct, self.userxAct, self.plotprefAct]
        # the history drop down menu
        self.historyMenu = QtWidgets.QMenu("&History", self)
        self.menubar.addMenu(self.historyMenu)
        self.historyAct = self.historyMenu.addAction(QtGui.QIcon(IconDirectory + 'history.png'), "&History", lambda: self.mainWindowCheck(lambda mainWindow: HistoryWindow(mainWindow)))
        self.historyAct.setToolTip('Show Processing History')
        self.errorAct = self.historyMenu.addAction(QtGui.QIcon(IconDirectory + 'error.png'), "&Error Messages", lambda: errorWindow(self))
        self.errorAct.setToolTip('Show Error Messages')
        self.historyActList = [self.historyAct]
        # Utilities dropdown menu
        self.utilitiesMenu = QtWidgets.QMenu("&Utilities", self)
        self.menubar.addMenu(self.utilitiesMenu)
        self.shiftconvAct = self.utilitiesMenu.addAction(QtGui.QIcon(IconDirectory + 'shifttool.png'), "&Chemical Shift Conversion Tool", self.createShiftConversionWindow)
        self.shiftconvAct.setToolTip('Chemical Shift Conversion Tool')
        self.dipolarconvAct = self.utilitiesMenu.addAction(QtGui.QIcon(IconDirectory + 'dipolar.png'), "Dipolar Distance Tool", self.createDipolarDistanceWindow)
        self.dipolarconvAct.setToolTip('Dipolar Distance Tool')
        self.quadconvAct = self.utilitiesMenu.addAction(QtGui.QIcon(IconDirectory + 'quadconversion.png'), "&Quadrupole Coupling Conversion Tool", self.createQuadConversionWindow)
        self.quadconvAct.setToolTip('Quadrupole Coupling Conversion Tool')
        self.mqmasconvAct = self.utilitiesMenu.addAction(QtGui.QIcon(IconDirectory + 'mqmas.png'), "MQMAS Parameter Extraction Tool", self.createMqmasExtractWindow)
        self.mqmasconvAct.setToolTip('MQMAS Parameter Extraction Tool')
        self.tempcalAct = self.utilitiesMenu.addAction(QtGui.QIcon(IconDirectory + 'temperature.png'), "Temperature Calibration Tool", self.createTempcalWindow)
        self.tempcalAct.setToolTip('Dipolar Distance Tool')
        self.nmrtableAct = self.utilitiesMenu.addAction(QtGui.QIcon(IconDirectory + 'table.png'), "&NMR Table", self.nmrTable)
        self.nmrtableAct.setToolTip('NMR Periodic Table')
        self.utilitiesActList = [self.shiftconvAct, self.quadconvAct, self.nmrtableAct, self.dipolarconvAct, self.mqmasconvAct, self.tempcalAct]
        # the help drop down menu
        self.helpMenu = QtWidgets.QMenu("&Help", self)
        self.menubar.addMenu(self.helpMenu)
        if not EXE:
            self.updateAct = self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'update.png'), "&Update", self.updateMenu)
            self.updateAct.setToolTip('Update ssNake')
            self.helpActList = [self.updateAct]
        else:
            self.helpActList = []
        self.refmanAct = self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'manual.png'), "Reference Manual", openRefMan)
        self.refmanAct.setToolTip('Open the Reference Manual')
        self.basTutorialAct = self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'Tutorial.png'), "Basic Tutorial", openTutorial)
        self.basTutorialAct.setToolTip('Open the Tutorial Folder')
        self.tutorialAct = self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'Tutorial.png'), "Advanced Tutorials", lambda: webbrowser.open('https://gitlab.science.ru.nl/mrrc/nmrzoo/ssnake_tutorials'))
        self.tutorialAct.setToolTip('Link to ssNake Advanced Processing Tutorials')
        self.githubAct = self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'GitHub.png'), "GitLab Page", lambda: webbrowser.open('https://gitlab.science.ru.nl/mrrc/nmrzoo/ssnake'))
        self.githubAct.setToolTip('ssNake GitLab Page')
        self.aboutAct = self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'about.png'), "&About", lambda: aboutWindow(self))
        self.aboutAct.setToolTip('About Menu')
        self.helpActList = self.helpActList +  [self.shiftconvAct, self.quadconvAct, self.nmrtableAct, self.githubAct,
                                                self.tutorialAct, self.aboutAct, self.basTutorialAct]
        # Extra event lists:
        self.specOnlyList = [self.regridAct, self.csastaticAct, self.quadAct, self.quadCSAAct, self.czjzekAct]
        self.fidOnlyList = [self.relaxAct, self.diffusionAct, self.swapEchoAct]
        self.Only1DPlot = [self.snrAct, self.fwhmAct, self.massAct, self.intfitAct]
        self.notInArrayPlot = [self.userxAct, self.setrefAct, self.swapEchoAct, self.corOffsetAct, self.baselineAct, self.subAvgAct,
                               self.refDeconvAct, self.intRegionAct, self.sumRegionAct, self.maxRegionAct, self.maxRegionAct,
                               self.minRegionAct, self.maxposRegionAct, self.minposRegionAct, self.averageRegionAct,
                               self.extractpartAct, self.matrixdelAct, self.normalizeAct, self.regridAct]

    def mainWindowCheck(self, transfer):
        # checks if mainWindow exist to execute the function
        if self.mainWindow is not None:
            transfer(self.mainWindow)
        else:
            raise SsnakeException("No workspaces open")

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        fileList = [url.toLocalFile() for url in event.mimeData().urls()]
        self.loadData(fileList)

    def menuCheck(self):
        if self.mainWindow is None:
            self.savemenu.menuAction().setEnabled(False)
            self.exportmenu.menuAction().setEnabled(False)
            self.workspacemenu.menuAction().setEnabled(False)
            self.macrolistmenu.menuAction().setEnabled(False)
            self.editmenu.menuAction().setEnabled(False)
            self.toolMenu.menuAction().setEnabled(True)
            self.matrixMenu.menuAction().setEnabled(False)
            self.transformsMenu.menuAction().setEnabled(False)
            self.fittingMenu.menuAction().setEnabled(False)
            self.combineMenu.menuAction().setEnabled(False)
            self.referencerunmenu.menuAction().setEnabled(False)
            for act in self.saveActList + self.exportActList + self.workspaceActList + self.macroActList + self.editActList + self.toolsActList + self.matrixActList + self.transformActList + self.fittingActList + self.plotActList + self.combineActList + self.historyActList:
                act.setEnabled(False)
        else:
            self.editmenu.menuAction().setEnabled(True)
            self.toolMenu.menuAction().setEnabled(True)
            self.matrixMenu.menuAction().setEnabled(True)
            self.transformsMenu.menuAction().setEnabled(True)
            self.fittingMenu.menuAction().setEnabled(True)
            self.combineMenu.menuAction().setEnabled(True)
            self.referencerunmenu.menuAction().setEnabled(True)
            for act in self.editActList + self.toolsActList + self.matrixActList + self.transformActList + self.fittingActList + self.plotActList + self.historyActList + self.combineActList:
                act.setEnabled(True)
            if type(self.mainWindow) is Main1DWindow:
                self.menuEnable(True)
                for act in self.specOnlyList:
                    act.setEnabled(int(self.mainWindow.current.spec() == 1))  # Only on for spec
                for act in self.fidOnlyList:
                    act.setEnabled(int(self.mainWindow.current.spec() == 0))  # Only on for FID
                  #Limit functions based on plot type
                if type(self.mainWindow.current) == views.CurrentMulti or type(self.mainWindow.current) == views.CurrentStacked or type(self.mainWindow.current) == views.CurrentArrayed:
                    for act in self.Only1DPlot:
                        act.setEnabled(False)
                if type(self.mainWindow.current) == views.CurrentArrayed:
                    for act in self.notInArrayPlot:
                        act.setEnabled(False)
                if self.mainWindow.masterData.noUndo:  # Set menu check to the same value as in the data
                    self.noUndoAct.setChecked(True)
                else:
                    self.noUndoAct.setChecked(False)
                if len(self.mainWindow.masterData.shape()) < 2:
                    for i in self.multiDActions:
                        i.setEnabled(False)
                else:
                    for i in self.multiDActions:
                        i.setEnabled(True)
                if not self.mainWindow.masterData.undoList:
                    self.undoAction.setEnabled(False)
                else:
                    self.undoAction.setEnabled(True)
                if not self.mainWindow.masterData.redoList:
                    self.redoAction.setEnabled(False)
                else:
                    self.redoAction.setEnabled(True)
                if not self.mainWindow.masterData.undoList and not self.mainWindow.masterData.redoList:
                    self.clearundoAct.setEnabled(False)
                self.savemenu.menuAction().setEnabled(True)
                self.exportmenu.menuAction().setEnabled(True)
                self.savefigAct.setEnabled(True)
                self.macrolistmenu.menuAction().setEnabled(True)
                if self.mainWindow.currentMacro is None:
                    self.macrostopAct.setEnabled(False)
                    self.macrostartAct.setEnabled(True)
                else:
                    self.macrostopAct.setEnabled(True)
                    self.macrostartAct.setEnabled(False)
                self.savemenu.menuAction().setEnabled(True)
                self.exportmenu.menuAction().setEnabled(True)
                self.workspacemenu.menuAction().setEnabled(True)
            elif type(self.mainWindow) is SaveFigureWindow:
                self.menuEnable(False, True)
                self.savemenu.menuAction().setEnabled(False)
                self.exportmenu.menuAction().setEnabled(False)
                for act in self.workspaceActList:
                    act.setEnabled(True)
                for act in self.saveActList + self.exportActList:
                    act.setEnabled(False)
                self.workspacemenu.menuAction().setEnabled(True)
                self.macrolistmenu.menuAction().setEnabled(False)
                self.workInfoAct.setEnabled(False)
            else: #Fitting menu
                self.menuEnable(False, True)
                self.savemenu.menuAction().setEnabled(False)
                self.exportmenu.menuAction().setEnabled(True)
                self.workspacemenu.menuAction().setEnabled(True)
                self.macrolistmenu.menuAction().setEnabled(False)
                for act in self.saveActList + self.exportActList + self.workspaceActList:
                    act.setEnabled(False)
                for act in self.workspaceActList:
                    act.setEnabled(True)
                self.savefigAct.setEnabled(True)
                self.workInfoAct.setEnabled(False)

    def menuEnable(self, enable=True, internalWindow=False):
        self.menuActive = enable or internalWindow

        self.macrolistmenu.menuAction().setEnabled(enable)
        self.editmenu.menuAction().setEnabled(enable)
        self.matrixMenu.menuAction().setEnabled(enable)
        self.transformsMenu.menuAction().setEnabled(enable)
        self.fittingMenu.menuAction().setEnabled(enable)
        self.combineMenu.menuAction().setEnabled(enable)
        self.referencerunmenu.menuAction().setEnabled(enable)
        # Actions:
        for act in self.macroActList + self.editActList + self.toolsActList + self.matrixActList + self.transformActList + self.fittingActList + self.plotActList + self.combineActList + self.historyActList:
            act.setEnabled(enable)
        if not internalWindow:
            self.tree.setEnabled(enable)
            self.filemenu.menuAction().setEnabled(enable)
            self.workspacemenu.menuAction().setEnabled(enable)
            self.toolMenu.menuAction().setEnabled(enable)
            for act in self.fileActList + self.workspaceActList:
                act.setEnabled(enable)
            for i in range(self.tabs.count()):
                if i != self.workspaceNum:
                    self.tabs.setTabEnabled(i, enable)
        self.undoAction.setEnabled(enable)
        self.redoAction.setEnabled(enable)

    def askName(self, filePath=None, name=None):
        if filePath is None:
            message = 'Spectrum name'
        else:
            message = 'Spectrum name for: ' + filePath
        count = 0
        if name is None:
            name = 'spectrum' + str(count)
        while name in self.workspaceNames:
            count += 1
            name = 'spectrum' + str(count)
        givenName, ok = QtWidgets.QInputDialog.getText(self, message, 'Name:', text=name)
        if not ok:
            return None
        while (givenName in self.workspaceNames) or givenName == '':
            self.dispMsg("Workspace name '" + givenName + "' already exists")
            givenName, ok = QtWidgets.QInputDialog.getText(self, message, 'Name:', text=name)
            if not ok:
                return None
        return givenName

    def undo(self, *args):
        if isinstance(self.mainWindow, Main1DWindow):
            self.mainWindow.undo()

    def redo(self, *args):
        if self.mainWindow is not None:
            self.mainWindow.redo()

    def macroCreate(self):
        if self.mainWindow is None:
            return
        if self.mainWindow.currentMacro is not None:
            return
        count = 0
        name = 'macro' + str(count)
        while name in self.macros.keys():
            count += 1
            name = 'macro' + str(count)
        givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)
        if not ok:
            return
        while (givenName in self.macros.keys()) or (givenName == ''):
            self.dispMsg("Macro name '" + givenName + "' already exists")
            givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)
            if not ok:
                return
        self.macros[givenName] = []
        self.mainWindow.redoMacro = []
        self.mainWindow.currentMacro = givenName
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        action1 = self.macrolistmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'), givenName, lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(QtGui.QIcon(IconDirectory + 'save.png'), givenName, lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'), givenName, lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'), givenName, lambda name=givenName: self.renameMacro(name))
        self.macroActions[givenName] = [action1, action2, action3, action4]
        self.menuCheck()

    def renameMacro(self, oldName):
        if self.mainWindow is None:
            return
        count = 0
        name = 'macro' + str(count)
        while name in self.macros.keys():
            count += 1
            name = 'macro' + str(count)
        givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)
        while (givenName in self.macros.keys()) or givenName == '':
            if not ok:
                return
            self.dispMsg("Macro name '" + givenName + "' already exists")
            givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)
        self.macros[givenName] = self.macros.pop(oldName)

        if self.mainWindow.currentMacro == oldName:
            self.mainWindow.currentMacro = givenName
                    

        oldActions = self.macroActions.pop(oldName)
        self.macrolistmenu.removeAction(oldActions[0])
        self.macrosavemenu.removeAction(oldActions[1])
        self.macrodeletemenu.removeAction(oldActions[2])
        self.macrorenamemenu.removeAction(oldActions[3])
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        action1 = self.macrolistmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'), givenName, lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(QtGui.QIcon(IconDirectory + 'save.png'), givenName, lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'), givenName, lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'), givenName, lambda name=givenName: self.renameMacro(name))
        self.macroActions[givenName] = [action1, action2, action3, action4]
        self.menuCheck()

    def stopMacro(self):
        if self.mainWindow is None:
            return
        if self.mainWindow.currentMacro is None:
            return
        self.mainWindow.redoMacro = []
        self.mainWindow.currentMacro = None
        self.menuCheck()

    def macroAdd(self, name, macros):
        self.macros[name].append(macros)

    def runMacro(self, name):
        if self.mainWindow is not None:
            self.mainWindow.runMacro(self.macros[name])

    def saveMacro(self, name):
        fileName = QtWidgets.QFileDialog.getSaveFileName(self, 'Save file', self.lastLocation + os.path.sep + name + '.macro', 'MACRO (*.macro)')
        if isinstance(fileName, tuple):
            fileName = fileName[0]
        if not fileName:
            return
        self.lastLocation = os.path.dirname(fileName)
        outputMacro = self.macros[name]
        with open(fileName, 'w') as f:
            for line in outputMacro:
                f.write(line[0])
                f.write(repr(line[1]).replace('\n', '').replace(' ', ''))
                f.write('\n')

    def deleteMacro(self, name):
        self.macrolistmenu.removeAction(self.macroActions[name][0])
        self.macrosavemenu.removeAction(self.macroActions[name][1])
        self.macrodeletemenu.removeAction(self.macroActions[name][2])
        self.macrorenamemenu.removeAction(self.macroActions[name][3])
        del self.macros[name]
        del self.macroActions[name]
        for i in self.workspaces:
            if i.currentMacro == name:
                i.currentMacro = None
        self.menuCheck()

    def loadMacro(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.lastLocation)
        if isinstance(filename, tuple):
            filename = filename[0]
        if filename:  # if not cancelled
            self.lastLocation = os.path.dirname(filename)  # Save used path
        if not filename:
            return
        self.stopMacro()
        count = 0
        name = 'macro' + str(count)
        while name in self.macros.keys():
            count += 1
            name = 'macro' + str(count)
        givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)
        while (givenName in self.macros.keys()) or (givenName == ''):
            if not ok:
                return
            self.dispMsg("Macro name '" + givenName + "' already exists")
            givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)
        with open(filename, 'r') as f:
            stringList = f.readlines()
        self.macros[givenName] = []
        for line in stringList:
            splitLine = line.split("(", 1)
            splitLine[0] = splitLine[0].replace(' ', '')
            splitLine[1] = safeEval("(" + splitLine[1])
            self.macros[givenName].append(splitLine)
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        action1 = self.macrolistmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'), givenName, lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(QtGui.QIcon(IconDirectory + 'save.png'), givenName, lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'), givenName, lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'), givenName, lambda name=givenName: self.renameMacro(name))
        self.macroActions[givenName] = [action1, action2, action3, action4]
        self.menuCheck()

    def referenceAdd(self, reffreq, name):
        self.referenceName.append(name)
        self.referenceValue.append(reffreq)  # List with saved refrence values
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        action1 = self.referencerunmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'), name, lambda name=name: self.referenceRun(name))
        action2 = self.referencedeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'), name, lambda name=name: self.referenceRemove(name))
        action3 = self.referencerenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'), name, lambda name=name: self.referenceRename(name))
        action4 = self.referencesavemenu.addAction(QtGui.QIcon(IconDirectory + 'save.png'), name, lambda name=name: self.referenceSave(name))
        self.referenceActions[name] = [action1, action2, action3, action4]
        self.menuCheck()

    def referenceClear(self):
        self.mainWindow.current.setRef(None)

    def referenceRun(self, name):
        reffreq = self.referenceValue[self.referenceName.index(name)]
        self.mainWindow.current.setRef(reffreq)

    def referenceRemove(self, name):
        self.referenceValue.remove(self.referenceValue[self.referenceName.index(name)])
        self.referenceName.remove(name)
        self.referencerunmenu.removeAction(self.referenceActions[name][0])
        self.referencedeletemenu.removeAction(self.referenceActions[name][1])
        self.referencerenamemenu.removeAction(self.referenceActions[name][2])
        self.referencesavemenu.removeAction(self.referenceActions[name][3])
        del self.referenceActions[name]
        self.menuCheck()

    def referenceRename(self, oldName):
        if self.mainWindow is None:
            return
        givenName, ok = QtWidgets.QInputDialog.getText(self, 'Reference name', 'Name:', text=oldName)
        if givenName == oldName or not ok:
            return
        while (givenName in self.referenceName) or (givenName == ''):
            self.dispMsg('Name exists')
            givenName, ok = QtWidgets.QInputDialog.getText(self, 'Reference name', 'Name:', text=oldName)
            if not ok:
                return
        self.referenceName[self.referenceName.index(oldName)] = givenName
        oldActions = self.referenceActions.pop(oldName)
        self.referencerunmenu.removeAction(oldActions[0])
        self.referencedeletemenu.removeAction(oldActions[1])
        self.referencerenamemenu.removeAction(oldActions[2])
        self.referencesavemenu.removeAction(oldActions[3])
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        action1 = self.referencerunmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'), givenName, lambda name=givenName: self.referenceRun(name))
        action2 = self.referencedeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'), givenName, lambda name=givenName: self.referenceRemove(name))
        action3 = self.referencerenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'), givenName, lambda name=givenName: self.referenceRename(name))
        action4 = self.referencesavemenu.addAction(QtGui.QIcon(IconDirectory + 'save.png'), givenName, lambda name=givenName: self.referenceSave(name))
        self.referenceActions[givenName] = [action1, action2, action3, action4]
        self.menuCheck()

    def referenceSave(self, name):
        fileName = QtWidgets.QFileDialog.getSaveFileName(self, 'Save reference', self.lastLocation + os.path.sep + name + '.txt', 'txt (*.txt)')
        if isinstance(fileName, tuple):
            fileName = fileName[0]
        if not fileName:
            return
        self.lastLocation = os.path.dirname(fileName)
        reffreq = self.referenceValue[self.referenceName.index(name)]
        with open(fileName, 'w') as f:
            f.write(str(reffreq))

    def referenceLoad(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.lastLocation)
        if isinstance(filename, tuple):
            filename = filename[0]
        if filename:  # if not cancelled
            self.lastLocation = os.path.dirname(filename)  # Save used path
        if not filename:
            return
        name = os.path.basename(filename)
        if name.endswith('.txt'): #If regular extension, name becomes filename - extension
            name = name[:-4]
        count = 0
        while name in self.referenceName: #If name known, cycle trough defaults
            name = 'ref' + str(count)
            count += 1
        givenName, ok = QtWidgets.QInputDialog.getText(self, 'Reference name', 'Name:', text=name)
        while (givenName in self.macros.keys()) or (givenName == ''):
            if not ok:
                return
            self.dispMsg('Name exists')
            givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)
        with open(filename, 'r') as f:
            self.referenceName.append(givenName)
            try:
                freq = float(f.read())
            except Exception:
                raise SsnakeException("Failed loading '" + filename + "' as reference.")
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        self.referenceValue.append(freq)
        action1 = self.referencerunmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'), givenName, lambda name=givenName: self.referenceRun(name))
        action2 = self.referencedeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'), givenName, lambda name=givenName: self.referenceRemove(name))
        action3 = self.referencerenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'), givenName, lambda name=givenName: self.referenceRename(name))
        action4 = self.referencesavemenu.addAction(QtGui.QIcon(IconDirectory + 'save.png'), givenName, lambda name=givenName: self.referenceSave(name))
        self.referenceActions[givenName] = [action1, action2, action3, action4]
        self.menuCheck()

    def noUndoMode(self, val):
        self.mainWindow.current.setNoUndo(val)
        self.menuCheck()

    def changeMainWindow(self, var):
        if not self.allowChange:
            return
        self.logo.hide()
        self.tabs.show()
        if isinstance(var, int):
            num = var
        else:
            num = self.workspaceNames.index(var)
        self.workspaceNum = num
        self.mainWindow = self.workspaces[num]
        self.tabs.setCurrentIndex(num)
        self.updWorkspaceMenu()
        self.menuCheck()
        try:
            if type(self.mainWindow.current) is views.CurrentMulti:
                self.mainWindow.sideframe.checkChanged()
        except Exception:
            pass

    def moveWorkspace(self, end, start):
        self.workspaces.insert(end, self.workspaces.pop(start))
        self.workspaceNames.insert(end, self.workspaceNames.pop(start))
        if self.workspaceNum == start:
            self.workspaceNum = end
        elif self.workspaceNum > end and self.workspaceNum < start:
            self.workspaceNum += 1
        elif self.workspaceNum < end and self.workspaceNum > start:
            self.workspaceNum -= 1

    def stepWorkspace(self, step):
        if len(self.workspaces) > 1:
            self.workspaceNum += step
            self.workspaceNum = self.workspaceNum % len(self.workspaces)
            self.mainWindow = self.workspaces[self.workspaceNum]
            self.tabs.setCurrentIndex(self.workspaceNum)
            self.updWorkspaceMenu()
            self.menuCheck()
            if type(self.mainWindow) is not SaveFigureWindow:
                if type(self.mainWindow.current) is views.CurrentMulti:
                    self.mainWindow.sideframe.checkChanged()

    def duplicateWorkspace(self, sliceOnly=False, *args):
        name = self.askName()
        if sliceOnly:
            data = copy.deepcopy(self.mainWindow.get_current().data1D)
            data.setNoUndo(self.mainWindow.get_masterData().noUndo)
        else:
            data = copy.deepcopy(self.mainWindow.get_masterData())
        if name is None:
            return
        self.workspaces.append(Main1DWindow(self, data, self.mainWindow.get_current()))
        self.workspaces[-1].rename(name)
        self.tabs.addTab(self.workspaces[-1], name)
        self.workspaceNames.append(name)
        self.changeMainWindow(name)

    def renameWorkspace(self, *args):
        tmp = self.workspaceNames[self.workspaceNum]
        self.workspaceNames[self.workspaceNum] = ''
        name = self.askName(tmp, tmp)
        if name is None:
            self.workspaceNames[self.workspaceNum] = tmp
            return
        self.workspaceNames[self.workspaceNum] = name
        self.tabs.setTabText(self.workspaceNum, name)
        self.updWorkspaceMenu()
        self.workspaces[self.workspaceNum].rename(name)

    def destroyWorkspace(self, num=None):
        if self.mainWindow is None or self.menuActive is False:
            return
        if num is None:
            num = self.workspaceNum
        self.allowChange = False
        self.tabs.removeTab(num)
        self.allowChange = True
        if num == self.workspaceNum:
            self.mainWindow.kill()
            self.mainWindow = None
        else:
            self.workspaces[num].kill()
        del self.workspaceNames[num]
        del self.workspaces[num]
        if num == self.workspaceNum:
            if num == len(self.workspaces):
                self.workspaceNum = num - 1
        if num < self.workspaceNum:
            self.workspaceNum -= 1
        if self.workspaces:
            self.changeMainWindow(self.workspaceNames[self.workspaceNum])
        else:
            self.logo.show()
            self.tabs.hide()
            self.updWorkspaceMenu()

    def updWorkspaceMenu(self):
        self.activemenu.clear()
        for i in self.workspaceNames:
            self.activemenu.addAction(i, lambda i=i: self.changeMainWindow(i))
        self.menuCheck()

    def newWorkspace(self, masterData):
        name = self.askName()
        if name is None:
            raise SsnakeException("No name given")
        self.workspaces.append(Main1DWindow(self, masterData))
        self.workspaces[-1].rename(name)
        self.tabs.addTab(self.workspaces[-1], name)
        self.workspaceNames.append(name)
        self.changeMainWindow(name)

    def createCombineWorkspaceWindow(self):
        CombineWorkspaceWindow(self)

    def createCombineLoadWindow(self):
        CombineLoadWindow(self)

    def combineWorkspace(self, combineNames):
        wsname = self.askName()
        if wsname is None:
            return
        i = self.workspaceNames.index(combineNames[0])
        combineMasterData = copy.deepcopy(self.workspaces[i].get_masterData())
        shapeRequired = combineMasterData.shape()
        hyperShape = len(combineMasterData.data)
        combineMasterData.split(1, -1)
        for name in combineNames[1:]:
            i = self.workspaceNames.index(name)
            addData = self.workspaces[i].get_masterData()
            if addData.shape() != shapeRequired:
                raise SsnakeException("Not all the data has the same shape")
            if len(addData.data) != hyperShape:
                raise SsnakeException("Not all the data has the same hypercomplex shape")
            combineMasterData.insert(addData.data, combineMasterData.shape()[0], 0)
        self.workspaces.append(Main1DWindow(self, combineMasterData))
        self.workspaces[-1].rename(wsname)
        self.tabs.addTab(self.workspaces[-1], wsname)
        self.workspaceNames.append(wsname)
        self.changeMainWindow(wsname)

    def loadFromMenu(self):
        fileList = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open File', self.lastLocation)
        if isinstance(fileList, tuple):
            fileList = fileList[0]
        self.loadData(fileList)

    def loadData(self, fileList):
        for filePath in fileList:
            if filePath:  # if not cancelled
                self.lastLocation = os.path.dirname(filePath)  # Save used path
            if not filePath:
                return
            masterData = io.autoLoad(filePath)
            if masterData is None:
                raise SsnakeException("Could not load data")
            if masterData == -1:
                dialog = AsciiLoadWindow(self, filePath)
                if dialog.exec_():
                    if dialog.closed:
                        return
                asciiInfo = (dialog.dataDimension, dialog.dataOrder, dialog.dataSpec, dialog.delim, dialog.sw, dialog.axisMulti)
                masterData = io.autoLoad(filePath, [asciiInfo])
            if self.defaultAskName:
                name = self.askName(filePath, masterData.name)
                if name is None:
                    return
            else:
                name = masterData.name
                count = 0
                while name in self.workspaceNames:
                    name = 'spectrum' + str(count)
                    count += 1
            masterData.rename(name)
            if masterData is not None:
                self.workspaces.append(Main1DWindow(self, masterData))
                self.tabs.addTab(self.workspaces[-1], name)
                self.workspaceNames.append(name)
                self.changeMainWindow(name)

    def loadFitLibDir(self):
        #fileName = QtWidgets.QFileDialog.getExistingDirectory(self, 'Open Library Directory', self.lastLocation)
        fileName = QtWidgets.QFileDialog.getOpenFileName(self, 'Open Library File', self.lastLocation)
        if isinstance(fileName, tuple):
            fileName = fileName[0]
        return fileName

    def loadSIMPSONScript(self):
        fileName = QtWidgets.QFileDialog.getOpenFileName(self, 'Open SIMPSON Script', self.lastLocation)
        if isinstance(fileName, tuple):
            fileName = fileName[0]
        return fileName

    def loadAndCombine(self, filePathList):
        masterData = io.autoLoad(filePathList)
        wsname = self.askName()
        if wsname is None:
            return
        masterData.rename(wsname)
        self.workspaces.append(Main1DWindow(self, masterData))
        self.workspaces[-1].rename(wsname)
        self.tabs.addTab(self.workspaces[-1], wsname)
        self.workspaceNames.append(wsname)
        self.changeMainWindow(wsname)

    def dataFromFit(self, data, filePath, freq, sw, spec, wholeEcho, ref, xaxArray, axes):
        name = self.askName()
        if name is None:
            return
        masterData = sc.Spectrum(data,
                                 filePath,
                                 freq,
                                 sw,
                                 spec=spec,
                                 wholeEcho=wholeEcho,
                                 ref=ref,
                                 xaxArray=xaxArray,
                                 history=['Data obtained from fit'],
                                 name=name)
        masterData.resetXax(axes)
        self.workspaces.append(Main1DWindow(self, masterData))
        self.tabs.addTab(self.workspaces[-1], name)
        self.workspaceNames.append(name)
        self.changeMainWindow(name)

    def saveSimpsonFile(self):
        self.mainWindow.get_mainWindow().SaveSimpsonFile()

    def saveASCIIFile(self):
        self.mainWindow.get_mainWindow().saveASCIIFile()

    def saveCSVFile(self):
        self.mainWindow.get_mainWindow().saveCSVFile()

    def saveJSONFile(self):
        self.mainWindow.get_mainWindow().saveJSONFile()

    def saveMatlabFile(self):
        self.mainWindow.get_mainWindow().saveMatlabFile()

    def saveFigure(self):
        if self.mainWindow is None:
            return
        self.allowChange = False
        self.menuEnable(False, True)
        num = self.workspaces.index(self.mainWindow)
        self.mainWindow = SaveFigureWindow(self, self.mainWindow)
        self.tabs.removeTab(num)
        self.tabs.insertTab(num, self.mainWindow, self.workspaceNames[num])
        self.workspaces[num] = self.mainWindow
        self.tabs.setCurrentIndex(num)
        self.menuCheck()
        self.allowChange = True

    def closeSaveFigure(self, mainWindow):
        self.allowChange = False
        num = self.workspaces.index(self.mainWindow)
        self.tabs.removeTab(num)
        self.mainWindow = mainWindow
        self.workspaces[num] = self.mainWindow
        self.tabs.insertTab(num, self.mainWindow, self.workspaceNames[num])
        self.tabs.setCurrentIndex(num)
        self.menuEnable(True, True)
        self.tabs.setCurrentIndex(num)
        self.menuCheck()
        self.allowChange = True

    def createFitWindow(self, fitWindow):
        if self.mainWindow is None:
            return
        self.allowChange = False
        self.menuEnable(False, True)
        num = self.workspaces.index(self.mainWindow)
        self.tabs.removeTab(num)
        self.mainWindow = fitWindow
        self.tabs.insertTab(num, self.mainWindow, self.workspaceNames[num])
        self.workspaces[num] = self.mainWindow
        self.tabs.setCurrentIndex(num)
        self.menuCheck()
        self.allowChange = True

    def closeFitWindow(self, mainWindow):
        self.allowChange = False
        num = self.workspaces.index(self.mainWindow)
        self.tabs.removeTab(num)
        del self.mainWindow
        self.mainWindow = mainWindow
        self.workspaces[num] = self.mainWindow
        self.tabs.insertTab(num, self.mainWindow, self.workspaceNames[num])
        self.menuEnable(True, True)
        self.tabs.setCurrentIndex(num)
        self.menuCheck()
        self.allowChange = True

    def updateMenu(self):
        UpdateWindow(self)

    def createShiftConversionWindow(self):
        shiftConversionWindow(self)

    def createDipolarDistanceWindow(self):
        dipolarDistanceWindow(self)

    def createQuadConversionWindow(self):
        quadConversionWindow(self)

    def createMqmasExtractWindow(self):
        mqmasExtractWindow(self)

    def createTempcalWindow(self):
        tempCalWindow(self)

    def nmrTable(self):
        import nmrTable
        nmrTable.PeriodicTable()

    def fileQuit(self):
        self.close()

    def closeEvent(self, event):
        quit_msg = "Are you sure you want to close ssNake?"
        close = True
        if len(self.workspaces) != 0:
            close = QtWidgets.QMessageBox.Yes == QtWidgets.QMessageBox.question(self, 'Close', quit_msg, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
        if close:
            for item in fit.stopDict.keys():  # Send stop commands to all threads
                fit.stopDict[item] = True
            event.accept()
        else:
            event.ignore()

######################################################################################################


class Main1DWindow(QtWidgets.QWidget):

    def __init__(self, father, masterData, duplicateCurrent=None):
        super(Main1DWindow, self).__init__(father)
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        grid = QtWidgets.QGridLayout(self)
        grid.addWidget(self.canvas, 0, 0)
        self.currentMacro = None
        self.redoMacro = []
        self.monitor = None  # Monitor of files
        self.monitorMacros = []
        self.father = father
        self.masterData = masterData
        if duplicateCurrent is not None:
            self.current = duplicateCurrent.copyCurrent(self, self.fig, self.canvas, masterData)
        else:
            self.current = views.Current1D(self, self.fig, self.canvas, masterData)
        self.menubar = self.father.menubar
        self.sideframe = SideFrame(self)
        grid.addWidget(self.sideframe, 0, 1)
        self.bottomframe = BottomFrame(self)
        grid.addWidget(self.bottomframe, 1, 0, 1, 2)
        self.textframe = TextFrame(self)
        grid.addWidget(self.textframe, 2, 0, 1, 2)
        grid.setColumnStretch(0, 1)
        grid.setRowStretch(0, 1)
        self.grid = grid
        self.canvas.mpl_connect('button_press_event', self.buttonPress)
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.canvas.setFocus()

    def rename(self, name):
        self.current.rename(name)

    def buttonPress(self, event):
        self.current.buttonPress(event)

    def buttonRelease(self, event):
        self.current.buttonRelease(event)

    def pan(self, event):
        self.current.pan(event)

    def scroll(self, event):
        self.current.scroll(event)

    def get_mainWindow(self):
        return self

    def get_masterData(self):
        return self.masterData

    def get_current(self):
        return self.current

    def kill(self):
        self.current.kill()
        del self.current
        del self.masterData
        del self.canvas
        del self.fig  # fig is destroyed
        for i in reversed(range(self.grid.count())):
            self.grid.itemAt(i).widget().deleteLater()
        self.grid.deleteLater()
        self.bottomframe.kill()
        del self.bottomframe
        self.textframe.kill()
        del self.textframe
        self.sideframe.kill()
        del self.sideframe
        self.deleteLater()
        gc.collect()

    def rescue(self):
        self.current.kill()
        self.current = views.Current1D(self, self.current.fig, self.current.canvas, self.masterData)

    def menuEnable(self, enable=True):
        self.father.menuEnable(enable)
        self.sideframe.frameEnable(enable)
        self.bottomframe.frameEnable(enable)
        self.textframe.frameEnable(enable)
        if enable:
            self.menuCheck()

    def menuCheck(self):
        self.father.menuCheck()

    def runMacro(self, macro, display=True):
        for i, _ in enumerate(macro):
            iter1 = macro[i] # Do not loop over the macro list itself to prevent recursion if the running macro is also the one being recorded
            self.addMacro(iter1)
            try:
                getattr(self.masterData, iter1[0])(*iter1[1])
            except AttributeError:
                raise SsnakeException('unknown macro command: ' + iter1[0])
        if display:
            self.current.upd()  # get the first slice of data
            self.current.showFid()  # plot the data
            self.current.plotReset()  # reset the axes limits
            self.updAllFrames()
            self.menuCheck()

    def addMacro(self, macroStep):
        if self.currentMacro is not None:
            self.father.macroAdd(self.currentMacro, macroStep)
            self.redoMacro = []

    def saveJSONFile(self):
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum]  # Set name of file to be saved to workspace name to start
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.lastLocation + os.path.sep + WorkspaceName + '.json', 'JSON (*.json)')
        if isinstance(name, tuple):
            name = name[0]
        if not name:
            return
        self.father.lastLocation = os.path.dirname(name)  # Save used path
        io.saveJSONFile(name, self.masterData)

    def saveMatlabFile(self):
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum]  # Set name of file to be saved to workspace name to start
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.lastLocation + os.path.sep + WorkspaceName + '.mat', 'MATLAB file (*.mat)')
        if isinstance(name, tuple):
            name = name[0]
        if not name:
            return
        self.father.lastLocation = os.path.dirname(name)  # Save used path
        io.saveMatlabFile(name, self.masterData, self.father.workspaceNames[self.father.workspaceNum])

    def SaveSimpsonFile(self):
        if self.masterData.ndim() > 2:
            raise SsnakeException('Saving to Simpson format only allowed for 1D and 2D data!')
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum]  # Set name of file to be saved to workspace name to start
        if sum(self.masterData.spec) / len(self.masterData.spec) == 1:
            name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.lastLocation + os.path.sep + WorkspaceName + '.spe', 'SIMPSON file (*.spe)')
            if isinstance(name, tuple):
                name = name[0]
            if not name:
                return
        elif sum(self.masterData.spec) == 0:
            name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.lastLocation + os.path.sep + WorkspaceName + '.fid', 'SIMPSON file (*.fid)')
            if isinstance(name, tuple):
                name = name[0]
            if not name:
                return
        else:
            raise SsnakeException('Saving to Simpson format not allowed for mixed time/frequency domain data!')
        self.father.lastLocation = os.path.dirname(name)  # Save used path
        io.saveSimpsonFile(name, self.masterData)

    def saveASCIIFile(self):
        if self.masterData.ndim() > 2:
            raise SsnakeException('Saving to ASCII format only allowed for 1D and 2D data!')
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum]  # Set name of file to be saved to workspace name to start
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.lastLocation + os.path.sep + WorkspaceName + '.txt', 'ASCII file (*.txt)')
        if isinstance(name, tuple):
            name = name[0]
        if not name:
            return
        self.father.lastLocation = os.path.dirname(name)  # Save used path
        axMult = self.current.getCurrentAxMult()
        io.saveASCIIFile(name, self.masterData, axMult)

    def saveCSVFile(self):
        if self.masterData.ndim() > 2:
            raise SsnakeException('Saving to CSV format only allowed for 1D and 2D data!')
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum]  # Set name of file to be saved to workspace name to start
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.lastLocation + os.path.sep + WorkspaceName + '.txt', 'CSV file (*.csv)')
        if isinstance(name, tuple):
            name = name[0]
        if not name:
            return
        self.father.lastLocation = os.path.dirname(name)  # Save used path
        axMult = self.current.getCurrentAxMult()
        io.saveASCIIFile(name, self.masterData, axMult,delim = ',')

    def reloadLast(self):
        self.current.reload()
        self.updAllFrames()
        self.menuCheck()
        gc.collect()

    def monitorLoad(self, filePath, delay=0.5):
        self.monitor.blockSignals(True)
        if not os.path.exists(filePath):
            self.stopMonitor()
            return
        loadData = io.autoLoad(*self.masterData.filePath)
        self.masterData.restoreData(loadData, None)
        for name in self.monitorMacros:
            self.runMacro(self.father.macros[name], display=False)
        self.current.upd()
        # self.current.plotReset()
        self.current.showFid()
        self.updAllFrames()
        self.menuCheck()
        QtCore.QTimer.singleShot(delay * 1000, lambda: self.monitor.blockSignals(False))
        if filePath in self.monitor.files() or filePath in self.monitor.directories():
            return
        self.monitor.addPath(filePath)

    def startMonitor(self, macroNames, delay=0.5):
        self.monitorMacros = macroNames
        self.monitor = QtCore.QFileSystemWatcher(self.masterData.filePath[0], self)
        self.monitor.fileChanged.connect(lambda a: self.monitorLoad(a, delay))
        self.monitor.directoryChanged.connect(lambda a: self.monitorLoad(a, delay))

    def stopMonitor(self):
        self.monitorMacros = []
        if self.monitor is not None:
            for name in self.masterData.filePath[0]:
                self.monitor.removePath(name)
        self.monitor = None

    def real(self):
        self.current.real()
        self.sideframe.upd()
        self.menuCheck()

    def imag(self):
        self.current.imag()
        self.sideframe.upd()
        self.menuCheck()

    def abs(self):
        self.current.abs()
        self.sideframe.upd()
        self.menuCheck()

    def conj(self):
        self.current.conj()
        self.sideframe.upd()
        self.menuCheck()

    def fourier(self):
        self.current.complexFourier()
        self.bottomframe.upd()
        self.menuCheck()

    def realFourier(self):
        self.current.realFourier()
        self.bottomframe.upd()
        self.menuCheck()

    def fftshift(self):
        self.current.fftshift()
        self.updAllFrames()
        self.menuCheck()

    def invFftshift(self):
        self.current.fftshift(inv=True)
        self.updAllFrames()
        self.menuCheck()

    def diff(self):
        self.current.diff()
        self.updAllFrames()
        self.menuCheck()

    def cumsum(self):
        self.current.cumsum()
        self.updAllFrames()
        self.menuCheck()

    def hilbert(self):
        self.current.hilbert()
        self.menuCheck()

    def states(self):
        self.current.states()
        self.updAllFrames()
        self.menuCheck()

    def statesTPPI(self):
        self.current.statesTPPI()
        self.updAllFrames()
        self.menuCheck()

    def echoAntiEcho(self):
        self.current.echoAntiEcho()
        self.updAllFrames()
        self.menuCheck()

    def setFreq(self, freq, sw):
        self.current.setFreq(freq, sw)
        self.menuCheck()

    def flipLR(self):
        self.current.flipLR()
        self.menuCheck()

    def directAutoPhase(self, phaseNum):
        self.current.directAutoPhase(phaseNum)
        self.menuCheck()

    def autoPhaseAll(self, phaseNum):
        self.current.autoPhaseAll(phaseNum)
        self.menuCheck()

    def CorrectDigitalFilter(self):
        if self.current.data.dFilter is None:
            raise SsnakeException('Digital filter: no value defined')
        self.current.correctDFilter()
        self.menuCheck()

    def createRelaxWindow(self):
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'relax'))

    def createDiffusionWindow(self):
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'diffusion'))

    def createPeakDeconvWindow(self):
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'peakdeconv'))

    def createCsaDeconvWindow(self):
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'csadeconv'))

    def createQuadDeconvWindow(self):
        if self.current.freq() == 0.0:
            raise SsnakeException("Please set the spectrometer frequency first!")
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'quaddeconv'))

    def createQuadCSADeconvWindow(self):
        if self.current.freq() == 0.0:
            raise SsnakeException("Please set the spectrometer frequency first!")
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'quadcsadeconv'))

    def createQuadCzjzekWindow(self):
        if self.current.freq() == 0.0:
            raise SsnakeException("Please set the spectrometer frequency first!")
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'quadczjzek'))

    def createMQMASWindow(self):
        if self.masterData.ndim() < 2:
            raise SsnakeException("Data has not enough dimensions for MQMAS fitting")
        if self.current.freq() == 0.0:
            raise SsnakeException("Please set the spectrometer frequency first!")
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'mqmas'))

    def createMQMASCzjzekWindow(self):
        if self.masterData.ndim() < 2:
            raise SsnakeException("Data has not enough dimensions for MQMAS fitting")
        if self.current.freq() == 0.0:
            raise SsnakeException("Please set the spectrometer frequency first!")
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'mqmasczjzek'))

    def createExternalFitWindow(self):
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'external'))

    def createFunctionFitWindow(self):
        self.father.createFitWindow(fit.TabFittingWindow(self.father, self.father.mainWindow, 'function'))

    def plot1D(self):
        tmpcurrent = views.Current1D(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()
        self.menuCheck()

    def plotScatter(self):
        tmpcurrent = views.CurrentScatter(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()
        self.menuCheck()

    def plotStack(self):
        if len(self.masterData.shape()) < 2:
            raise SsnakeException("Data does not have enough dimensions")
        tmpcurrent = views.CurrentStacked(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()
        self.menuCheck()

    def plotArray(self):
        if len(self.masterData.shape()) < 2:
            raise SsnakeException("Data does not have enough dimensions")
        tmpcurrent = views.CurrentArrayed(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()
        self.menuCheck()

    def plotContour(self):
        if len(self.masterData.shape()) < 2:
            raise SsnakeException("Data does not have enough dimensions")
        tmpcurrent = views.CurrentContour(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()
        self.menuCheck()

    def plotMultiContour(self):
        if len(self.masterData.shape()) < 2:
            raise SsnakeException("Data does not have enough dimensions")
        tmpcurrent = views.CurrentMultiContour(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()
        self.menuCheck()

    def plotColour2D(self):
        if len(self.masterData.shape()) < 2:
            raise SsnakeException("Data does not have enough dimensions")
        tmpcurrent = views.CurrentColour2D(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()
        self.menuCheck()

    def plotMulti(self):
        tmpcurrent = views.CurrentMulti(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()
        self.menuCheck()

    def updAllFrames(self):
        self.sideframe.upd()
        self.bottomframe.upd()
        self.textframe.upd()

    def undo(self, *args):
        self.father.dispMsg(self.masterData.undo())
        self.current.upd()
        self.current.showFid()
        self.current.plotReset()
        self.updAllFrames()
        if self.currentMacro is not None:
            self.redoMacro.append(self.father.macros[self.currentMacro].pop())
        self.menuCheck()

    def redo(self, *args):
        self.masterData.redo()
        self.current.upd()
        self.current.showFid()
        self.current.plotReset()
        self.updAllFrames()
        if self.currentMacro is not None:
            self.father.macroAdd(self.currentMacro, self.redoMacro.pop())
        self.menuCheck()

    def clearUndo(self):
        self.masterData.clearUndo()
        self.menuCheck()

########################################################################################


class SideFrame(QtWidgets.QScrollArea):

    FITTING = False

    def __init__(self, parent):
        super(SideFrame, self).__init__(parent)
        self.father = parent
        self.entries = []
        self.plotIs2D = False
        content = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout(content)
        grid.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        frame1Widget = QtWidgets.QWidget()
        frame2Widget = QtWidgets.QWidget()
        grid.addWidget(frame1Widget, 0, 0)
        grid.addWidget(frame2Widget, 1, 0)
        self.frame1 = QtWidgets.QGridLayout()
        self.frame2 = QtWidgets.QGridLayout()
        frame1Widget.setLayout(self.frame1)
        frame2Widget.setLayout(self.frame2)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.grid = grid
        self.setWidget(content)
        self.upd()

    def kill(self):
        for i in reversed(range(self.grid.count())):
            self.grid.itemAt(i).widget().deleteLater()
        self.grid.deleteLater()

    def frameEnable(self, enable=True):
        self.setEnabled(enable)

    def upd(self):
        current = self.father.current
        self.shape = current.data.shape()
        self.length = len(self.shape)
        for i in reversed(range(self.frame1.count())):
            item = self.frame1.itemAt(i).widget()
            if self.FITTING:
                item.hide()
            self.frame1.removeWidget(item)
            item.deleteLater()
        for i in reversed(range(self.frame2.count())):
            item = self.frame2.itemAt(i).widget()
            if self.FITTING:
                item.hide()
            self.frame2.removeWidget(item)
            item.deleteLater()
        offset = 0
        self.plotIs2D = isinstance(current, views.CurrentStacked)
        if self.plotIs2D:
            offset = 1
        self.entries = []
        self.buttons1 = []
        self.buttons1Group = QtWidgets.QButtonGroup(self)
        self.buttons1Group.buttonClicked.connect(lambda: self.setAxes(True))
        self.buttons2 = []
        self.buttons2Group = QtWidgets.QButtonGroup(self)
        self.buttons2Group.buttonClicked.connect(lambda: self.setAxes(False))
        if self.length > 1:
            for num in range(self.length):
                if not self.FITTING:
                    self.buttons1.append(QtWidgets.QRadioButton(''))
                    self.buttons1Group.addButton(self.buttons1[num], num)
                    self.buttons1[num].setToolTip(TOOLTIPS['sideframeDimension1'])
                    self.frame1.addWidget(self.buttons1[num], num * 2 + 1, 0)
                    if self.plotIs2D:
                        self.buttons2.append(QtWidgets.QRadioButton(''))
                        self.buttons2Group.addButton(self.buttons2[num], num)
                        self.buttons2[num].setToolTip(TOOLTIPS['sideframeDimension2'])
                        self.frame1.addWidget(self.buttons2[num], num * 2 + 1, 1)
                if current.isComplex(num):
                    tmpLabel = "*D"
                else:
                    tmpLabel = "D"
                self.frame1.addWidget(wc.QLabel(tmpLabel + str(num + 1), self), num * 2, 1 + offset)
                self.entries.append(wc.SliceSpinBox(self, 0, self.shape[num] - 1))
                self.entries[-1].setToolTip(TOOLTIPS['sideFrameDimensionSlice'])
                self.frame1.addWidget(self.entries[num], num * 2 + 1, 1 + offset)
                if not self.plotIs2D:
                    self.entries[num].setValue(current.locList[num])
                else:
                    self.entries[num].setValue(current.locList[num])
                if self.FITTING and num in current.axes:
                    self.entries[num].setDisabled(True)
                self.entries[num].valueChanged.connect(lambda event, num=num: self.getSlice(num))
            if type(current) in (views.CurrentStacked, views.CurrentArrayed):
                if current.viewSettings["stackBegin"] is not None:
                    from2D = current.viewSettings["stackBegin"]
                else:
                    from2D = 0
                if current.viewSettings["stackEnd"] is not None:
                    to2D = current.viewSettings["stackEnd"]
                else:
                    to2D = self.shape[current.axes[-2]]
                if current.viewSettings["stackStep"] is not None:
                    step2D = current.viewSettings["stackStep"]
                else:
                    step2D = 1
                self.frame2.addWidget(wc.QLabel("From", self), 1, 0)
                self.fromSpin = wc.SliceSpinBox(self, 0, to2D - 1)
                self.fromSpin.setToolTip(TOOLTIPS['sideFrom'])
                self.frame2.addWidget(self.fromSpin, 2, 0)
                self.fromSpin.setValue(from2D)
                self.fromSpin.valueChanged.connect(self.setToFrom)
                self.frame2.addWidget(wc.QLabel("To", self), 3, 0)
                self.toSpin = wc.SliceSpinBox(self, from2D + 1, self.shape[current.axes[-2]])
                self.toSpin.setToolTip(TOOLTIPS['sideTo'])
                self.frame2.addWidget(self.toSpin, 4, 0)
                self.toSpin.setValue(to2D)
                self.toSpin.valueChanged.connect(self.setToFrom)
                self.frame2.addWidget(wc.QLabel("Step", self), 5, 0)
                self.stepSpin = wc.SliceSpinBox(self, 1, self.shape[current.axes[-2]])
                self.stepSpin.setToolTip(TOOLTIPS['stackStep'])
                self.frame2.addWidget(self.stepSpin, 6, 0)
                self.stepSpin.setValue(step2D)
                self.stepSpin.valueChanged.connect(self.setToFrom)
                self.frame2.addWidget(wc.QLabel("Spacing", self), 7, 0)
                self.spacingEntry = QtWidgets.QLineEdit(self)
                self.spacingEntry.setToolTip(TOOLTIPS['stackSpacing'])
                self.spacingEntry.setText('%#.3g' % current.viewSettings["spacing"])
                self.spacingEntry.returnPressed.connect(self.setSpacing)
                self.frame2.addWidget(self.spacingEntry, 8, 0)
            if isinstance(current, (views.CurrentContour)):
                if type(current) in (views.CurrentContour, views.CurrentMultiContour):
                    self.contourTypeGroup = QtWidgets.QGroupBox('Contour type:')
                    self.contourTypeFrame = QtWidgets.QGridLayout()
                    self.contourNumberLabel = wc.QLeftLabel("Number:", self)
                    self.contourTypeFrame.addWidget(self.contourNumberLabel, 0, 0)
                    self.numLEntry = wc.SsnakeSpinBox()
                    self.numLEntry.setMaximum(100000)
                    self.numLEntry.setMinimum(1)
                    self.numLEntry.setToolTip(TOOLTIPS['contourNumber'])
                    self.numLEntry.setValue(current.viewSettings["numLevels"])
                    self.numLEntry.valueChanged.connect(self.setContour)
                    self.contourTypeFrame.addWidget(self.numLEntry, 0, 1)
                    self.contourTypeFrame.addWidget(wc.QLeftLabel("Sign:", self), 1, 0)
                    self.contourSignEntry = QtWidgets.QComboBox()
                    self.contourSignEntry.setToolTip(TOOLTIPS['contourSign'])
                    self.contourSignEntry.addItems(['Both', '+ only', '- only'])
                    self.contourSignEntry.setCurrentIndex(current.viewSettings["contourSign"])
                    self.contourSignEntry.currentIndexChanged.connect(self.setContour)
                    self.contourTypeFrame.addWidget(self.contourSignEntry, 1, 1)
                    self.contourTypeLabel = wc.QLeftLabel("Type:", self)
                    self.contourTypeFrame.addWidget(self.contourTypeLabel, 2, 0)
                    self.contourTypeEntry = QtWidgets.QComboBox()
                    self.contourTypeEntry.setToolTip(TOOLTIPS['contourType'])
                    self.contourTypeEntry.addItems(['Linear', 'Multiplier'])
                    self.contourTypeEntry.setCurrentIndex(current.viewSettings["contourType"])
                    self.contourTypeEntry.currentIndexChanged.connect(self.setContour)
                    self.contourTypeFrame.addWidget(self.contourTypeEntry, 2, 1)
                    self.multiValueLabel = wc.QLeftLabel("Multiplier:", self)
                    self.contourTypeFrame.addWidget(self.multiValueLabel, 3, 0)
                    self.multiValue = wc.QLineEdit(current.viewSettings["multiValue"], self.setContour)
                    self.multiValue.setToolTip(TOOLTIPS['contourMultiplier'])
                    self.multiValue.setMaximumWidth(120)
                    self.contourTypeFrame.addWidget(self.multiValue, 3, 1)
                    if current.viewSettings["contourType"] != 1:
                        self.multiValueLabel.hide()
                        self.multiValue.hide()
                    self.contourTypeGroup.setLayout(self.contourTypeFrame)
                    self.frame2.addWidget(self.contourTypeGroup, 6, 0, 1, 3)
                    # Contour limits
                    self.contourLimitsGroup = QtWidgets.QGroupBox('Contour limits [%]:')
                    self.contourLimitsFrame = QtWidgets.QGridLayout()
                    self.maxLEntry = wc.QLineEdit(format(current.viewSettings["maxLevels"] * 100.0, '.7g'), self.setContour)
                    self.maxLEntry.setMaximumWidth(120)
                    self.maxLEntry.setToolTip(TOOLTIPS['contourMax'])
                    self.contourLimitsFrame.addWidget(self.maxLEntry, 1, 1)
                    self.minLEntry = wc.QLineEdit(format(current.viewSettings["minLevels"] * 100.0, '.7g'), self.setContour)
                    self.minLEntry.setMaximumWidth(120)
                    self.minLEntry.setToolTip(TOOLTIPS['contourMin'])
                    self.contourLimitsFrame.addWidget(self.minLEntry, 2, 1)
                    self.contourLimType = QtWidgets.QComboBox()
                    self.contourLimType.addItems(['Current 2D', 'Full data'])
                    self.contourLimType.setCurrentIndex(current.viewSettings["limitType"])
                    self.contourLimType.setToolTip(TOOLTIPS['contourLimType'])
                    self.contourLimType.currentIndexChanged.connect(self.setContour)
                    self.contourLimitsFrame.addWidget(self.contourLimType, 0, 1)
                    self.maxLabel = wc.QLeftLabel("Max:", self)
                    self.minLabel = wc.QLeftLabel("Min:", self)
                    self.relLabel = wc.QLeftLabel("Rel. to:", self)
                    self.contourLimitsFrame.addWidget(self.relLabel, 0, 0)
                    self.contourLimitsFrame.addWidget(self.maxLabel, 1, 0)
                    self.contourLimitsFrame.addWidget(self.minLabel, 2, 0)
                    self.contourLimitsGroup.setLayout(self.contourLimitsFrame)
                    self.frame2.addWidget(self.contourLimitsGroup, 7, 0, 1, 3)
                # Projections
                self.contourProjGroup = QtWidgets.QGroupBox('Projections:')
                self.contourProjFrame = QtWidgets.QGridLayout()
                self.projTopLabel = wc.QLeftLabel("Top:", self)
                self.contourProjFrame.addWidget(self.projTopLabel, 0, 0)
                self.projDropTop = QtWidgets.QComboBox()
                self.projDropTop.setToolTip(TOOLTIPS['contourTopProjection'])
                self.projDropTop.addItems(["Sum", "Max", "Min", "Off", "Slice", "Diagonal"])
                self.projDropTop.setCurrentIndex(current.viewSettings["projTop"])
                self.projDropTop.activated.connect(lambda val, self=self: self.changeProj(val, 1))
                self.contourProjFrame.addWidget(self.projDropTop, 0, 1)
                self.projTraceTop = wc.SsnakeSpinBox()
                self.projTraceTop.setMaximum(self.shape[current.axes[-2]] - 1)
                self.projTraceTop.setMinimum(0)
                self.projTraceTop.setValue(current.viewSettings["projPos"][0])
                self.projTraceTop.valueChanged.connect(lambda val, self=self: self.changeTrace(val, 0))
                self.projTraceTop.setToolTip(TOOLTIPS['contourProjTopTrac'])
                self.contourProjFrame.addWidget(self.projTraceTop, 1, 1)
                if current.viewSettings["projTop"] != 4:
                    self.projTraceTop.hide()
                self.projRightLabel = wc.QLeftLabel("Right:", self)
                self.contourProjFrame.addWidget(self.projRightLabel, 2, 0)
                self.projDropRight = QtWidgets.QComboBox()
                self.projDropRight.setToolTip(TOOLTIPS['contourRightProjection'])
                self.projDropRight.addItems(["Sum", "Max", "Min", "Off", "Slice", "Diagonal"])
                self.projDropRight.setCurrentIndex(current.viewSettings["projRight"])
                self.projDropRight.activated.connect(lambda val, self=self: self.changeProj(val, 2))
                self.contourProjFrame.addWidget(self.projDropRight, 2, 1)
                self.projTraceRight = wc.SsnakeSpinBox()
                self.projTraceRight.setMaximum(self.shape[current.axes[-1]] - 1)
                self.projTraceRight.setMinimum(0)
                self.projTraceRight.setValue(current.viewSettings["projPos"][1])
                self.projTraceRight.valueChanged.connect(lambda val, self=self: self.changeTrace(val, 1))
                self.projTraceRight.setToolTip(TOOLTIPS['contourProjRightTrac'])
                self.contourProjFrame.addWidget(self.projTraceRight, 3, 1)
                if current.viewSettings["projRight"] != 4:
                    self.projTraceRight.hide()
                self.selectTraceButton = QtWidgets.QPushButton("Select slices", self)
                self.selectTraceButton.clicked.connect(self.selectTraces)
                self.contourProjFrame.addWidget(self.selectTraceButton, 4, 1)
                if (current.viewSettings["projTop"] != 4) and (current.viewSettings["projRight"] != 4):
                    self.selectTraceButton.hide()
                # Ranges
                self.rangeCheckbox = QtWidgets.QCheckBox('Projection ranges', self)
                self.rangeCheckbox.setChecked(current.viewSettings["projLimitsBool"])
                self.rangeCheckbox.stateChanged.connect(self.activateRanges)
                self.rangeCheckbox.setToolTip(TOOLTIPS['contourProjRanges'])
                self.contourProjFrame.addWidget(self.rangeCheckbox, 5, 0, 1, 2)
                self.projTopRangeMaxLabel = wc.QLeftLabel("Top max:", self)
                self.projTopRangeMaxLabel.hide()
                self.contourProjFrame.addWidget(self.projTopRangeMaxLabel, 6, 0)
                self.projTopRangeMax = wc.SsnakeSpinBox()
                self.projTopRangeMax.setMaximum(self.shape[current.axes[-2]] - 1)
                self.projTopRangeMax.setMinimum(0)
                self.projTopRangeMax.setToolTip(TOOLTIPS['contourTopRangeMax'])
                if current.viewSettings["projLimits"][0] is None:
                    self.projTopRangeMax.setValue(self.shape[current.axes[-2]] - 1)
                else:
                    self.projTopRangeMax.setValue(current.viewSettings["projLimits"][0])
                self.projTopRangeMax.valueChanged.connect(self.changeRanges)
                self.projTopRangeMax.hide()
                self.contourProjFrame.addWidget(self.projTopRangeMax, 6, 1)
                self.projTopRangeMinLabel = wc.QLeftLabel("Top min:", self)
                self.projTopRangeMinLabel.hide()
                self.contourProjFrame.addWidget(self.projTopRangeMinLabel, 7, 0)
                self.projTopRangeMin = wc.SsnakeSpinBox()
                self.projTopRangeMin.setMaximum(self.shape[current.axes[-2]] - 1)
                self.projTopRangeMin.setMinimum(0)
                self.projTopRangeMin.setToolTip(TOOLTIPS['contourTopRangeMin'])
                if current.viewSettings["projLimits"][1] is None:
                    self.projTopRangeMin.setValue(0)
                else:
                    self.projTopRangeMin.setValue(current.viewSettings["projLimits"][1])
                self.projTopRangeMin.valueChanged.connect(self.changeRanges)
                self.projTopRangeMin.hide()
                self.contourProjFrame.addWidget(self.projTopRangeMin, 7, 1)
                self.projRightRangeMaxLabel = wc.QLeftLabel("Right max:", self)
                self.projRightRangeMaxLabel.hide()
                self.contourProjFrame.addWidget(self.projRightRangeMaxLabel, 8, 0)
                self.projRightRangeMax = wc.SsnakeSpinBox()
                self.projRightRangeMax.setMaximum(self.shape[current.axes[-1]] - 1)
                self.projRightRangeMax.setMinimum(0)
                self.projRightRangeMax.setToolTip(TOOLTIPS['contourRightRangeMax'])
                if current.viewSettings["projLimits"][2] is None:
                    self.projRightRangeMax.setValue(self.shape[current.axes[-1]] - 1)
                else:
                    self.projRightRangeMax.setValue(current.viewSettings["projLimits"][2])
                self.projRightRangeMax.valueChanged.connect(self.changeRanges)
                self.projRightRangeMax.hide()
                self.contourProjFrame.addWidget(self.projRightRangeMax, 8, 1)
                self.projRightRangeMinLabel = wc.QLeftLabel("Right min:", self)
                self.contourProjFrame.addWidget(self.projRightRangeMinLabel, 9, 0)
                self.projRightRangeMinLabel.hide()
                self.projRightRangeMin = wc.SsnakeSpinBox()
                self.projRightRangeMin.setMaximum(self.shape[current.axes[-1]] - 1)
                self.projRightRangeMin.setMinimum(0)
                self.projRightRangeMin.setToolTip(TOOLTIPS['contourRightRangeMin'])
                if current.viewSettings["projLimits"][3] is None:
                    self.projRightRangeMin.setValue(0)
                else:
                    self.projRightRangeMin.setValue(current.viewSettings["projLimits"][3])
                self.projRightRangeMin.valueChanged.connect(self.changeRanges)
                self.projRightRangeMin.hide()
                self.contourProjFrame.addWidget(self.projRightRangeMin, 9, 1)
                self.contourProjGroup.setLayout(self.contourProjFrame)
                self.frame2.addWidget(self.contourProjGroup, 9, 0, 1, 3)
                self.activateRanges(self.rangeCheckbox.checkState())
                # Diagonal group
                self.diagonalGroup = QtWidgets.QGroupBox('Diagonal:')
                self.diagonalGroup.setCheckable(True)
                self.diagonalGroup.setChecked(current.viewSettings["diagonalBool"])
                self.diagonalGroup.toggled.connect(self.switchDiagonal)
                self.diagonalGroup.setToolTip(TOOLTIPS['contourDiagonal'])
                self.diagonalFrame = QtWidgets.QGridLayout()
                self.diagMultiLabel = wc.QLeftLabel("Multiplier:", self)
                self.diagonalFrame.addWidget(self.diagMultiLabel, 0, 0)
                self.diagonalEntry = wc.QLineEdit(current.viewSettings["diagonalMult"], self.setDiagonal)
                self.diagonalEntry.setMaximumWidth(120)
                self.diagonalEntry.setToolTip(TOOLTIPS['contourDiagonalMulti'])
                self.diagonalFrame.addWidget(self.diagonalEntry, 0, 1)
                self.diagonalGroup.setLayout(self.diagonalFrame)
                self.frame2.addWidget(self.diagonalGroup, 10, 0, 1, 3)
            if not self.FITTING:
                self.buttons1Group.button(current.axes[-1]).toggle()
                if self.plotIs2D:
                    self.buttons2Group.button(current.axes[-2]).toggle()
        if isinstance(current, (views.CurrentMulti, views.CurrentMultiContour)):
            self.extraEntries = []
            self.extraButtons1 = []
            self.extraButtons1Group = []
            self.extraButtons2 = []
            self.extraButtons2Group = []
            self.nameLabels = []
            iter1 = 0
            for i in range(len(current.viewSettings["extraData"])):
                frameWidget = QtWidgets.QWidget(self)
                frame = QtWidgets.QGridLayout(frameWidget)
                self.frame2.addWidget(frameWidget, iter1, 0)
                frameWidget.setLayout(frame)
                name = current.viewSettings["extraName"][i]
                if len(name) > 20:
                    nameLabel = wc.QLabel(name[:20] + '…', self)
                    nameLabel.setToolTip(name)
                else:
                    nameLabel = wc.QLabel(name, self)
                self.nameLabels.append(nameLabel)
                frame.addWidget(self.nameLabels[i], 0, 0, 1, 3)
                self.nameLabels[i].setStyleSheet("QLabel { color: rgb" + str(current.getExtraColor(i)) + ";}")
                colorbutton = QtWidgets.QPushButton("Colour", self)
                colorbutton.clicked.connect(lambda arg, num=i: self.setExtraColor(num))
                colorbutton.setToolTip(TOOLTIPS['multiplotColour'])
                frame.addWidget(colorbutton, 1, 0)
                button = QtWidgets.QPushButton("x", self)
                button.clicked.connect(lambda arg, num=i: self.delMultiSpec(num))
                button.setToolTip(TOOLTIPS['multiplotX'])
                if isinstance(current, (views.CurrentMulti)):
                    frame.addWidget(button, 1, 1)
                    self.OOM = self.father.current.getOOM()  # Order of Magnitude
                    self.scaleLabel = wc.QLeftLabel("Scale:", self)
                    frame.addWidget(self.scaleLabel, 2, 0)
                    self.offsetLabel = wc.QLeftLabel(u"Offset (×1e" + str(self.OOM) + "):", self)
                    frame.addWidget(self.offsetLabel, 3, 0)
                    self.shiftLabel = wc.QLeftLabel("Shift:", self)
                    frame.addWidget(self.shiftLabel, 4, 0)
                    scaleEntry = wc.SsnakeDoubleSpinBox()
                    scaleEntry.setDecimals(4)
                    scaleEntry.setMaximum(1e6)
                    scaleEntry.setMinimum(-1e6)
                    scaleEntry.setSingleStep(0.1)
                    scaleEntry.setValue(self.father.current.viewSettings["extraScale"][i])
                    scaleEntry.valueChanged.connect(lambda arg, num=i: self.setScale(arg, num))
                    scaleEntry.setToolTip(TOOLTIPS['multiplotScale'])
                    frame.addWidget(scaleEntry, 2, 1)
                    offsetEntry = wc.SsnakeDoubleSpinBox()
                    offsetEntry.setDecimals(4)
                    offsetEntry.setMaximum(1e3)
                    offsetEntry.setMinimum(-1e3)
                    offsetEntry.setSingleStep(0.1)
                    offsetEntry.setValue(self.father.current.viewSettings["extraOffset"][i] / (10**self.OOM))
                    offsetEntry.valueChanged.connect(lambda arg, num=i: self.setOffset(arg, num))
                    offsetEntry.setToolTip(TOOLTIPS['multiplotOffset'])
                    frame.addWidget(offsetEntry, 3, 1)
                    shiftEntry = wc.SsnakeDoubleSpinBox()
                    shiftEntry.setDecimals(4)
                    shiftEntry.setMaximum(1e3)
                    shiftEntry.setMinimum(-1e3)
                    shiftEntry.setSingleStep(0.1)
                    shiftEntry.setValue(self.father.current.viewSettings["extraShift"][i])
                    shiftEntry.valueChanged.connect(lambda arg, num=i: self.setShift(arg, num))
                    shiftEntry.setToolTip(TOOLTIPS['multiplotShift1'])
                    frame.addWidget(shiftEntry, 4, 1)
                elif isinstance(current, (views.CurrentMultiContour)):
                    frame.addWidget(button, 1, 1, 1, 2)
                    self.OOM = self.father.current.getOOM()  # Order of Magnitude
                    self.scaleLabel = wc.QLeftLabel("Scale:", self)
                    frame.addWidget(self.scaleLabel, 2, 0)
                    self.shift1Label = wc.QLeftLabel("x Shift:", self)
                    frame.addWidget(self.shift1Label, 3, 0)
                    self.shift2Label = wc.QLeftLabel("y Shift:", self)
                    frame.addWidget(self.shift2Label, 4, 0)
                    scaleEntry = wc.SsnakeDoubleSpinBox()
                    scaleEntry.setDecimals(4)
                    scaleEntry.setMaximum(1e6)
                    scaleEntry.setMinimum(-1e6)
                    scaleEntry.setSingleStep(0.1)
                    scaleEntry.setValue(self.father.current.viewSettings["extraScale"][i])
                    scaleEntry.valueChanged.connect(lambda arg, num=i: self.setScale(arg, num))
                    scaleEntry.setToolTip(TOOLTIPS['multiplotScale'])
                    frame.addWidget(scaleEntry, 2, 1, 1, 2)
                    shiftEntry = wc.SsnakeDoubleSpinBox()
                    shiftEntry.setDecimals(4)
                    shiftEntry.setMaximum(1e3)
                    shiftEntry.setMinimum(-1e3)
                    shiftEntry.setSingleStep(0.1)
                    shiftEntry.setValue(self.father.current.viewSettings["extraShift"][i])
                    shiftEntry.valueChanged.connect(lambda arg, num=i: self.setShift(arg, num))
                    shiftEntry.setToolTip(TOOLTIPS['multiplotShift1'])
                    frame.addWidget(shiftEntry, 3, 1, 1, 2)
                    shiftEntry = wc.SsnakeDoubleSpinBox()
                    shiftEntry.setDecimals(4)
                    shiftEntry.setMaximum(1e3)
                    shiftEntry.setMinimum(-1e3)
                    shiftEntry.setSingleStep(0.1)
                    shiftEntry.setValue(self.father.current.viewSettings["extraShift2"][i])
                    shiftEntry.valueChanged.connect(lambda arg, num=i: self.setShift2(arg, num))
                    shiftEntry.setToolTip(TOOLTIPS['multiplotShift2'])
                    frame.addWidget(shiftEntry, 4, 1, 1, 2)
                entries = []
                self.extraEntries.append(entries)
                buttons1 = []
                self.extraButtons1.append(buttons1)
                self.extraButtons1Group.append(QtWidgets.QButtonGroup(self))
                self.extraButtons1Group[i].buttonClicked.connect(lambda: self.setExtraAxes(True))
                buttons2 = []
                self.extraButtons2.append(buttons1)
                self.extraButtons2Group.append(QtWidgets.QButtonGroup(self))
                self.extraButtons2Group[i].buttonClicked.connect(lambda: self.setExtraAxes(False))
                if current.viewSettings["extraData"][i].ndim() > 1:
                    for num in range(current.viewSettings["extraData"][i].ndim()):
                        offset = 0
                        buttons1.append(QtWidgets.QRadioButton(''))
                        buttons1[-1].setToolTip(TOOLTIPS['multiplotDim1'])
                        self.extraButtons1Group[i].addButton(buttons1[num], num)
                        frame.addWidget(buttons1[num], num * 3 + 6, 0)
                        if self.plotIs2D:
                            offset = 1
                            buttons2.append(QtWidgets.QRadioButton(''))
                            buttons2[-1].setToolTip(TOOLTIPS['multiplotDim2'])
                            self.extraButtons2Group[i].addButton(buttons2[num], num)
                            frame.addWidget(buttons2[num], num * 3 + 6, 1)
                        frame.addWidget(wc.QLabel("D" + str(num + 1), self), num * 3 + 5, 1 + offset)
                        entries.append(wc.SliceSpinBox(self, 0, current.viewSettings["extraData"][i].shape()[num] - 1))
                        entries[-1].setToolTip(TOOLTIPS['sideFrameDimensionSlice'])
                        frame.addWidget(entries[num], num * 3 + 6, 1 + offset)
                        entries[num].setValue(current.viewSettings["extraLoc"][i][num])
                        entries[num].valueChanged.connect(lambda event, num=num, i=i: self.getExtraSlice(num, i))
                    self.extraButtons1Group[i].button(current.viewSettings["extraAxes"][i][-1]).toggle()
                    if self.plotIs2D:
                        self.extraButtons2Group[i].button(current.viewSettings["extraAxes"][i][-2]).toggle()
                iter1 += 1
            addButton = QtWidgets.QPushButton("Add plot", self)
            addButton.setToolTip(TOOLTIPS['multiplotAddPlot'])
            addButton.clicked.connect(self.addMultiSpec)
            self.frame2.addWidget(addButton, iter1, 0, 1, 2)
        QtCore.QTimer.singleShot(100, self.resizeAll)

    def resizeAll(self):
        self.setMinimumWidth(self.grid.sizeHint().width() + self.verticalScrollBar().sizeHint().width())

    def setToFrom(self, *args):
        current = self.father.current
        if not (type(current) is views.CurrentStacked or type(current) is views.CurrentArrayed):
            return
        fromVar = self.fromSpin.value()
        toVar = self.toSpin.value()
        stepVar = self.stepSpin.value()
        current.stackSelect(fromVar, toVar, stepVar)
        self.fromSpin.setMaximum(toVar - 1)
        self.toSpin.setMinimum(fromVar + 1)

    def scrollSpacing(self, var):
        self.spacingEntry.setText('%#.3g' % var)

    def setSpacing(self, *args):
        var = safeEval(self.spacingEntry.text(), length=self.father.current.len(), Type='FI')
        self.spacingEntry.setText('%#.3g' % var)
        self.father.current.setSpacing(var)

    def setContour(self, *args):
        var1 = self.numLEntry.value()
        maxC = safeEval(self.maxLEntry.text(), length=self.father.current.len(), Type='FI')
        if maxC is None:
            maxC = self.father.current.viewSettings["maxLevels"] * 100
            self.father.father.dispMsg('Invalid value for contour maximum')
        else:
            maxC = abs(float(maxC))
        minC = safeEval(self.minLEntry.text(), length=self.father.current.len(), Type='FI')
        if minC is None:
            minC = self.father.current.viewSettings["minLevels"] * 100
            self.father.father.dispMsg('Invalid value for contour minimum')
        else:
            minC = abs(float(minC))
        if minC > maxC:  # if wrong order, interchange
            maxC, minC = (minC, maxC)
        self.maxLEntry.setText(str(maxC))
        self.minLEntry.setText(str(minC))
        cSign = self.contourSignEntry.currentIndex()
        cType = self.contourTypeEntry.currentIndex()
        if cType == 0:
            self.multiValue.hide()
            self.multiValueLabel.hide()
        else:
            self.multiValue.show()
            self.multiValueLabel.show()
        multi = safeEval(self.multiValue.text(), length=self.father.current.len(), Type='FI')
        if multi is None or multi <= 1.0:
            multi = self.father.current.viewSettings["multiValue"]
            self.father.father.dispMsg('Invalid value for contour multiplier')
        else:
            multi = abs(float(multi))
        self.multiValue.setText(str(multi))
        limitType = self.contourLimType.currentIndex()
        self.father.current.setLevels(var1, maxC / 100.0, minC / 100.0, limitType, cSign, cType, multi)

    def changeProj(self, pType, direc):
        if pType == 4:
            if direc == 1:
                self.projTraceTop.show()
            else:
                self.projTraceRight.show()
        else:
            self.selectTraceButton.hide()
            if direc == 1:
                self.projTraceTop.hide()
            else:
                self.projTraceRight.hide()
        self.father.current.setProjType(pType, direc)
        if (self.father.current.viewSettings["projTop"] == 4) or (self.father.current.viewSettings["projRight"] == 4):
            self.selectTraceButton.show()
        else:
            self.selectTraceButton.hide()
        # if not self.FITTING:
        #     self.father.current.clearProj()
        #     self.father.current.showAllProj()
        # else:
        #     self.father.current.showFid()

    def changeTrace(self, num, direc):
        self.father.current.setProjTraces(num, direc)
        # if not self.FITTING:
        #     self.father.current.clearProj()
        #     self.father.current.showAllProj()
        # else:
        #     self.father.current.showFid()

    def selectTraces(self, *args):
        self.father.current.peakPickFunc = lambda pos, self=self: self.pickedTraces(pos)
        self.father.current.peakPick = 3

    def pickedTraces(self, pos):
        if self.father.current.viewSettings["projTop"] == 4 and pos[3] != self.projTraceTop.value():
            self.projTraceTop.setValue(pos[3])
        if self.father.current.viewSettings["projRight"] == 4 and pos[3] != self.projTraceRight.value():
            self.projTraceRight.setValue(pos[0])

    def changeRanges(self):
        check = self.rangeCheckbox.isChecked()
        ranges = [self.projTopRangeMax.value(), self.projTopRangeMin.value(), self.projRightRangeMax.value(), self.projRightRangeMin.value()]
        self.father.current.setProjLimits(check, ranges)
        # if not self.FITTING:
        #     self.father.current.clearProj()
        #     self.father.current.showAllProj()
        # else:
        #     self.father.current.showFid()

    def activateRanges(self, state):
        if state:
            self.projTopRangeMaxLabel.show()
            self.projTopRangeMax.show()
            self.projTopRangeMinLabel.show()
            self.projTopRangeMin.show()
            self.projRightRangeMaxLabel.show()
            self.projRightRangeMax.show()
            self.projRightRangeMinLabel.show()
            self.projRightRangeMin.show()
        else:
            self.projTopRangeMaxLabel.hide()
            self.projTopRangeMax.hide()
            self.projTopRangeMinLabel.hide()
            self.projTopRangeMin.hide()
            self.projRightRangeMaxLabel.hide()
            self.projRightRangeMax.hide()
            self.projRightRangeMinLabel.hide()
            self.projRightRangeMin.hide()
        self.changeRanges()

    def setAxes(self, first=True):
        axes = self.buttons1Group.checkedId()
        if self.plotIs2D:
            axes2 = self.buttons2Group.checkedId()
            if axes == axes2:
                if first:
                    axes2 = self.father.current.axes[-1]
                else:
                    axes = self.father.current.axes[-2]
                if isinstance(self.father.current, (views.CurrentContour)):  # If contour
                    # Correct proj values and maxima
                    newRanges = [self.projRightRangeMax.value(), self.projRightRangeMin.value(), self.projTopRangeMax.value(), self.projTopRangeMin.value()]
                    self.father.current.setProjLimits(self.rangeCheckbox.isChecked(), newRanges)
                    self.father.current.setProjTraces(self.projTraceTop.value(), 1)
                    self.father.current.setProjTraces(self.projTraceRight.value(), 0)
                    #Flip diagonal multiplier:
                    inp = safeEval(self.diagonalEntry.text(), length=self.father.current.len(), Type='FI')
                    if inp is not None:
                        self.father.current.viewSettings["diagonalMult"] = 1.0 / inp
                    #Make sure the bottom frame nicely inverts the axis units
                    time1 = self.father.bottomframe.axisDropTime.currentIndex()
                    time2 = self.father.bottomframe.axisDropTime2.currentIndex()
                    freq1 = self.father.bottomframe.axisDropFreq.currentIndex()
                    freq2 = self.father.bottomframe.axisDropFreq2.currentIndex()
                    self.father.bottomframe.axisDropTime.setCurrentIndex(time2)
                    self.father.bottomframe.axisDropTime2.setCurrentIndex(time1)
                    self.father.bottomframe.axisDropFreq.setCurrentIndex(freq2)
                    self.father.bottomframe.axisDropFreq2.setCurrentIndex(freq1)
                    if bool(self.father.current.spec()) is True:
                        tmp1 = freq1
                    else:
                        tmp1 = time1
                    if bool(self.father.current.spec(-2)) is True:
                        tmp2 = freq2
                    else:
                        tmp2 = time2
                    self.father.bottomframe.changeAxis(tmp2, update=False)
                    self.father.bottomframe.changeAxis2(tmp1, update=False)
            self.buttons2Group.button(axes2).toggle()
        self.getSlice(axes, True)
        self.upd()

    def getSlice(self, entryNum, button=False):
        axisChange = False
        if button:
            dimNum = entryNum
            axisChange = True
        elif not self.plotIs2D:
            if entryNum == self.father.current.axes[-1]:
                if entryNum == self.length - 1:
                    dimNum = self.length - 2
                    axisChange = True
                else:
                    dimNum = self.length - 1
                    axisChange = True
            else:
                dimNum = self.father.current.axes[-1]
        else:
            dimNum = self.father.current.axes[-1]
        locList = np.array(self.father.current.locList, dtype=int)
        for num in range(self.length):
            locList[num] = self.entries[num].value()
        if not self.FITTING:
            self.buttons1Group.button(dimNum).toggle()
            axes = np.array([self.buttons2Group.checkedId(), dimNum], dtype=int)
        else:
            axes = self.father.current.axes
        if self.plotIs2D:
            self.father.current.setSlice(axes, locList)
        else:
            self.father.current.setSlice(np.array([dimNum]), locList)
        if not self.FITTING:
            self.father.bottomframe.upd()
            if axisChange:
                self.father.menuCheck()
                self.upd()

    def setExtraAxes(self, first=True):
        for i in range(len(self.extraButtons1Group)):
            axes = self.extraButtons1Group[i].checkedId()
            if self.plotIs2D:
                axes2 = self.extraButtons2Group[i].checkedId()
                if axes == axes2:
                    if first:
                        axes2 = self.father.current.viewSettings["extraAxes"][i][-1]
                    else:
                        axes = self.father.current.viewSettings["extraAxes"][i][-2]
                self.extraButtons2Group[i].button(axes2).toggle()
            self.getExtraSlice(axes, i, True)
        self.father.current.showFid()

    def getExtraSlice(self, entryNum, entryi, button=False):
        length = self.father.current.viewSettings["extraData"][entryi].ndim()
        if button:
            dimNum = entryNum
        else:
            if entryNum == self.father.current.viewSettings["extraAxes"][entryi][-1]:
                if entryNum == length - 1:
                    dimNum = length - 2
                else:
                    dimNum = length - 1
            else:
                dimNum = self.father.current.viewSettings["extraAxes"][entryi][-1]
        locList = np.array(self.father.current.viewSettings["extraLoc"][entryi])
        for num in range(length):
            locList[num] = self.extraEntries[entryi][num].value()
        self.extraButtons1Group[entryi].button(dimNum).toggle()
        axes = np.array([self.extraButtons2Group[entryi].checkedId(), dimNum], dtype=int)
        if self.plotIs2D:
            self.father.current.setExtraSlice(entryi, axes, locList)
        else:
            self.father.current.setExtraSlice(entryi, np.array([dimNum]), locList)
        if not button:
            self.father.current.showFid()
        # self.upd()

    def setScale(self, scale, num):
        self.father.current.setExtraScale(num, scale)

    def setOffset(self, offset, num):
        self.father.current.setExtraOffset(num, offset * 10**self.OOM)

    def setShift(self, shift, num):
        self.father.current.setExtraShift(num, shift)

    def setShift2(self, shift, num):
        self.father.current.setExtraShift2(num, shift)

    def switchDiagonal(self, val):
        self.father.current.setDiagonal(bool(val))

    def setDiagonal(self):
        inp = safeEval(self.diagonalEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            inp = self.father.current.viewSettings["diagonalMult"]
            self.father.father.dispMsg('Invalid value for diagonal multiplier')
        else:
            inp = float(inp)
        self.diagonalEntry.setText(str(inp))
        self.father.current.setDiagonal(None, inp)

    def checkChanged(self):
        for i in range(len(self.father.current.viewSettings["extraData"])):
            extraData = self.father.current.viewSettings["extraData"][i]
            if extraData.ndim() > 1:
                for j in range(len(self.extraEntries[i])):
                    self.extraEntries[i][j].setMaximum(extraData.data[0].shape[j] - 1)
            self.upd()
            self.father.current.showFid()

    def setExtraColor(self, num):
        color = QtWidgets.QColorDialog.getColor()
        if not color.isValid():
            return
        self.father.current.setExtraColor(num, color.getRgbF())
        self.nameLabels[num].setStyleSheet("QLabel { color: rgb" + str(self.father.current.getExtraColor(num)) + ";}")

    def addMultiSpec(self, *args):
        text = QtWidgets.QInputDialog.getItem(self, "Select spectrum to show", "Spectrum name:", self.father.father.workspaceNames, 0, False)
        if text[1]:
            self.father.current.addExtraData(self.father.father.workspaces[self.father.father.workspaceNames.index(text[0])].get_masterData(), str(text[0]))
            self.upd()

    def delMultiSpec(self, num):
        self.father.current.delExtraData(num)
        self.upd()

################################################################################


class BottomFrame(QtWidgets.QWidget):

    def __init__(self, parent):
        super(BottomFrame, self).__init__(parent)
        self.father = parent
        grid = QtWidgets.QGridLayout(self)
        self.setLayout(grid)
        fourierButton = QtWidgets.QPushButton("Fourier", parent=self)
        fourierButton.setToolTip(TOOLTIPS['Fourier'])
        fourierButton.clicked.connect(self.father.fourier)
        grid.addWidget(fourierButton, 0, 0, 2, 1)
        self.specGroup = QtWidgets.QButtonGroup(self)
        self.specGroup.buttonClicked.connect(self.changeSpec)
        timeButton = QtWidgets.QRadioButton('Time', parent=self)
        timeButton.setToolTip(TOOLTIPS['timeButton'])
        self.specGroup.addButton(timeButton, 0)
        grid.addWidget(timeButton, 0, 1)
        freqButton = QtWidgets.QRadioButton('Frequency', parent=self)
        freqButton.setToolTip(TOOLTIPS['freqButton'])
        self.specGroup.addButton(freqButton, 1)
        grid.addWidget(freqButton, 1, 1)
        self.wholeEcho = QtWidgets.QCheckBox("Whole echo", parent=self)
        self.wholeEcho.setToolTip(TOOLTIPS['wholeEcho'])
        self.wholeEcho.clicked.connect(self.setWholeEcho)
        grid.addWidget(self.wholeEcho, 0, 2, 2, 1)
        grid.addWidget(wc.QLabel("Freq [MHz]:", self), 0, 3)
        self.freqEntry = wc.QLineEdit('', self.changeFreq, parent=self)
        self.freqEntry.setToolTip(TOOLTIPS['freqEntry'])
        grid.addWidget(self.freqEntry, 1, 3)
        grid.addWidget(wc.QLabel("Sweepwidth [kHz]:", self), 0, 4)
        self.swEntry = wc.QLineEdit('', self.changeFreq, parent=self)
        self.swEntry.setToolTip(TOOLTIPS['swEntry'])
        grid.addWidget(self.swEntry, 1, 4)
        grid.addWidget(wc.QLabel("Plot:", self), 0, 5)
        self.plotDrop = QtWidgets.QComboBox(parent=self)
        self.plotDrop.addItems(["Real", "Imag", "Both", "Abs"])
        self.plotDrop.setToolTip(TOOLTIPS['plotDrop'])
        self.plotDrop.activated.connect(self.changePlot)
        grid.addWidget(self.plotDrop, 1, 5)
        grid.addWidget(wc.QLabel("Axis:", self), 0, 6)
        self.axisDropTime = QtWidgets.QComboBox(parent=self)
        self.axisDropTime.setToolTip(TOOLTIPS['axisDrop'])
        self.axisDropTime.addItems(["s", "ms", u"μs"])
        self.axisDropTime.activated.connect(self.changeAxis)
        grid.addWidget(self.axisDropTime, 1, 6)
        self.axisDropFreq = QtWidgets.QComboBox(parent=self)
        self.axisDropFreq.addItems(["Hz", "kHz", "MHz", "ppm"])
        self.axisDropFreq.setToolTip(TOOLTIPS['axisDrop'])
        self.axisDropFreq.activated.connect(self.changeAxis)
        grid.addWidget(self.axisDropFreq, 1, 6)
        self.ax2Label = wc.QLabel("Axis2:", self)
        grid.addWidget(self.ax2Label, 0, 7)
        self.axisDropTime2 = QtWidgets.QComboBox(parent=self)
        self.axisDropTime2.addItems(["s", "ms", u"μs"])
        self.axisDropTime2.setToolTip(TOOLTIPS['axis2Drop'])
        self.axisDropTime2.activated.connect(self.changeAxis2)
        grid.addWidget(self.axisDropTime2, 1, 7)
        self.axisDropFreq2 = QtWidgets.QComboBox(parent=self)
        self.axisDropFreq2.setToolTip(TOOLTIPS['axis2Drop'])
        self.axisDropFreq2.addItems(["Hz", "kHz", "MHz", "ppm"])
        self.axisDropFreq2.activated.connect(self.changeAxis2)
        grid.addWidget(self.axisDropFreq2, 1, 7)
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        self.grid = grid
        self.upd()

    def kill(self):
        for i in reversed(range(self.grid.count())):
            self.grid.itemAt(i).widget().deleteLater()
        self.grid.deleteLater()

    def frameEnable(self, enable=True):
        self.setEnabled(enable)

    def upd(self):
        self.freqEntry.setText('%.6f' % (self.father.current.freq() / 1000000.0))
        self.swEntry.setText('%.6f' % (self.father.current.sw() / 1000.0))
        self.axisDropTime2.hide()
        self.axisDropFreq2.hide()
        self.axisDropFreq.model().item(3).setEnabled(True)
        if self.father.current.spec() == 0:
            self.specGroup.button(0).toggle()
            self.axisDropFreq.hide()
            self.axisDropTime.show()
            self.ax2Label.hide()
            self.axisDropTime.setCurrentIndex(self.father.current.getAxType())
        elif self.father.current.spec() == 1:
            self.specGroup.button(1).toggle()
            self.axisDropTime.hide()
            self.axisDropFreq.show()
            if self.father.current.freq() == 0.0 or self.father.current.ref() == 0.0:
                self.axisDropFreq.model().item(3).setEnabled(False)
            self.ax2Label.hide()
            if self.father.current.getppm() and self.father.current.freq() != 0.0 and self.father.current.ref() != 0.0:
                self.axisDropFreq.setCurrentIndex(3)
            else:
                self.axisDropFreq.setCurrentIndex(self.father.current.getAxType())
        if isinstance(self.father.current, views.CurrentContour):
            self.ax2Label.show()
            self.axisDropFreq2.model().item(3).setEnabled(True)
            if self.father.current.spec(-2) == 0:
                self.axisDropTime2.show()
                self.axisDropTime2.setCurrentIndex(self.father.current.getAxType(-2))
            elif self.father.current.spec(-2) == 1:
                self.axisDropFreq2.show()
                if self.father.current.freq(-2) == 0.0 or self.father.current.ref(-2) == 0.0:
                    self.axisDropFreq2.model().item(3).setEnabled(False)
                if self.father.current.getppm(-2) and self.father.current.freq(-2) != 0.0 and self.father.current.ref(-2) != 0.0:
                    self.axisDropFreq2.setCurrentIndex(3)
                else:
                    self.axisDropFreq2.setCurrentIndex(self.father.current.getAxType(-2))
        if type(self.father.current) is views.CurrentArrayed:
            self.ax2Label.show()
            self.axisDropFreq2.model().item(3).setEnabled(True)
            if self.father.current.spec(-2) == 0:
                self.axisDropTime2.show()
                self.axisDropTime2.setCurrentIndex(self.father.current.getAxType(-2))
            elif self.father.current.spec(-2) == 1:
                self.axisDropFreq2.show()
                if self.father.current.freq(-2) == 0.0:
                    self.axisDropFreq2.model().item(3).setEnabled(False)
                if self.father.current.getppm(-2) and self.father.current.freq(-2) != 0.0 and self.father.current.ref(-2) != 0.0:
                    self.axisDropFreq2.setCurrentIndex(3)
                else:
                    self.axisDropFreq2.setCurrentIndex(self.father.current.getAxType(-2))
        if self.father.current.wholeEcho():
            self.wholeEcho.setCheckState(QtCore.Qt.Checked)
        else:
            self.wholeEcho.setCheckState(QtCore.Qt.Unchecked)

    def setWholeEcho(self, inp):
        self.father.current.setWholeEcho(inp)
        self.father.menuCheck()

    def changeSpec(self):
        self.father.current.setSpec(self.specGroup.checkedId())
        self.upd()
        self.father.menuCheck()

    def changeFreq(self):
        freq = safeEval(self.freqEntry.text(), length=self.father.current.len(), Type='FI')
        sw = safeEval(self.swEntry.text(), length=self.father.current.len(), Type='FI')
        if sw is None:
            self.father.father.dispMsg('Invalid sweepwidth')
        elif sw == 0.0:
            sw = None
            self.father.father.dispMsg('Sweepwidth cannot be 0')
        else:
            sw *= 1000
        if freq is not None:
            freq *= 1e6
        else:
            self.father.father.dispMsg('Invalid spectrum frequency')
        self.father.setFreq(freq, sw)
        self.upd()

    def changePlot(self, pType):
        self.father.current.viewSettings["plotType"] = pType
        self.father.current.showFid()

    def changeAxis(self, pType, update=True):
        self.father.current.setAxType(pType, update)

    def changeAxis2(self, pType, update=True):
        self.father.current.setAxType(pType, update, -2)

##################################################################


class TextFrame(QtWidgets.QScrollArea):

    def __init__(self, parent):
        super(TextFrame, self).__init__(parent)
        self.father = parent
        self.oldx = 0.0
        self.oldy = 0.0
        self.oldamp = 0.0
        widthScale = 0.6
        content = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout(content)
        getButton = QtWidgets.QPushButton("&Get Position")
        getButton.setToolTip(TOOLTIPS['GetPos'])
        getButton.clicked.connect(self.getPosition)
        grid.addWidget(getButton, 0, 1)
        grid.addWidget(wc.QLabel("x-Position:"), 0, 2)
        self.xpos = wc.QLineEdit("0")
        self.xpos.setReadOnly(True)
        self.xpos.setToolTip(TOOLTIPS['xPosition'])
        self.xpos.setFixedWidth(int(self.xpos.sizeHint().width() * widthScale))
        grid.addWidget(self.xpos, 0, 3)
        self.yposlabel = wc.QLabel("y-Position:")
        grid.addWidget(self.yposlabel, 0, 4)
        self.ypos = wc.QLineEdit("0")
        self.ypos.setToolTip(TOOLTIPS['yPosition'])
        self.ypos.setReadOnly(True)
        self.ypos.setFixedWidth(int(self.ypos.sizeHint().width() * widthScale))
        grid.addWidget(self.ypos, 0, 5)
        grid.addWidget(wc.QLabel("x-Value:"), 0, 6)
        self.xpoint = wc.QLineEdit("0.0")
        self.xpoint.setToolTip(TOOLTIPS['xValue'])
        self.xpoint.setReadOnly(True)
        self.xpoint.setFixedWidth(int(self.xpoint.sizeHint().width() * widthScale))
        grid.addWidget(self.xpoint, 0, 7)
        self.ylabel = wc.QLabel("y-Value:")
        grid.addWidget(self.ylabel, 0, 8)
        self.ypoint = wc.QLineEdit("0.0")
        self.ypoint.setToolTip(TOOLTIPS['yValue'])
        self.ypoint.setReadOnly(True)
        self.ypoint.setFixedWidth(int(self.ypoint.sizeHint().width() * widthScale))
        grid.addWidget(self.ypoint, 0, 9)
        grid.addWidget(wc.QLabel("Amp:"), 0, 10)
        self.amppoint = wc.QLineEdit("0.0")
        self.amppoint.setToolTip(TOOLTIPS['ampValue'])
        self.amppoint.setReadOnly(True)
        self.amppoint.setFixedWidth(int(self.amppoint.sizeHint().width() * widthScale))
        grid.addWidget(self.amppoint, 0, 11)
        grid.addWidget(wc.QLabel(u"Δx:"), 0, 12)
        self.deltaxpoint = wc.QLineEdit("0.0")
        self.deltaxpoint.setToolTip(TOOLTIPS['deltaxvalue'])
        self.deltaxpoint.setReadOnly(True)
        self.deltaxpoint.setFixedWidth(int(self.deltaxpoint.sizeHint().width() * widthScale))
        grid.addWidget(self.deltaxpoint, 0, 13)
        self.deltaylabel = wc.QLabel(u"Δy:")
        grid.addWidget(self.deltaylabel, 0, 14)
        self.deltaypoint = wc.QLineEdit("0.0")
        self.deltaypoint.setToolTip(TOOLTIPS['deltayvalue'])
        self.deltaypoint.setReadOnly(True)
        self.deltaypoint.setFixedWidth(int(self.deltaypoint.sizeHint().width() * widthScale))
        grid.addWidget(self.deltaypoint, 0, 15)
        grid.addWidget(wc.QLabel(u"Δamp:"), 0, 16)
        self.deltaamppoint = wc.QLineEdit("0.0")
        self.deltaamppoint.setToolTip(TOOLTIPS['deltaamplitude'])
        self.deltaamppoint.setReadOnly(True)
        self.deltaamppoint.setFixedWidth(int(self.deltaamppoint.sizeHint().width() * widthScale))
        grid.addWidget(self.deltaamppoint, 0, 17)
        grid.setColumnStretch(20, 1)
        self.grid = grid
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setWidget(content)
        self.setMaximumHeight(self.grid.sizeHint().height() + self.horizontalScrollBar().sizeHint().height())
        self.upd()

    def upd(self):
        if isinstance(self.father.current, views.CurrentContour):
            self.ypos.show()
            self.yposlabel.show()
            self.ypoint.show()
            self.deltaypoint.show()
            self.ylabel.show()
            self.deltaylabel.show()
        else:
            self.ypos.hide()
            self.yposlabel.hide()
            self.ypoint.hide()
            self.deltaypoint.hide()
            self.ylabel.hide()
            self.deltaylabel.hide()

    def kill(self):
        for i in reversed(range(self.grid.count())):
            self.grid.itemAt(i).widget().deleteLater()
        self.grid.deleteLater()

    def frameEnable(self, enable=True):
        for child in self.children():
            child.setEnabled(enable)

    def setLabels(self, position):
        if len(position) > 3:
            self.ypos.setText(str(position[3]))
            self.deltaypoint.setText(('%#.'+str(self.father.father.defaultPrecis)+'g') % np.abs(self.oldy - position[4]))
            self.ypoint.setText(('%#.'+str(self.father.father.defaultPrecis)+'g') % position[4])
            self.oldy = position[4]
        self.deltaxpoint.setText(('%#.'+str(self.father.father.defaultPrecis)+'g') % np.abs(self.oldx - position[1]))
        self.deltaamppoint.setText(('%#.'+str(self.father.father.defaultPrecis)+'g') % np.abs(self.oldamp - position[2]))
        self.xpos.setText(str(position[0]))
        self.xpoint.setText(('%#.'+str(self.father.father.defaultPrecis)+'g') % position[1])
        self.amppoint.setText(('%#.'+str(self.father.father.defaultPrecis)+'g') % position[2])
        self.oldx = position[1]
        self.oldamp = position[2]

    def getPosition(self, *args):
        self.father.current.peakPickFunc = lambda pos, self=self: self.setLabels(pos)
        if isinstance(self.father.current, views.CurrentContour):
            self.father.current.peakPick = 2
        else:
            self.father.current.peakPick = True

#################################################################################

class AsciiLoadWindow(QtWidgets.QDialog):

    dataOrders = ['XRI', 'XR', 'XI', 'RI', 'R']
    delimiters = ['Tab', 'Space', 'Comma']
    timeUnits = ['s','ms',u'μs']
    timeMultiVals = [1.0,1.0e-3,1.0e-6]
    freqUnits = ['Hz','kHz','MHz','ppm']
    freqMultiVals = [1.0, 1.0e3, 1.0e6, None]

    def __init__(self, parent, file):
        super(AsciiLoadWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool | QtCore.Qt.WindowContextHelpButtonHint)
        self.dataDimension = 1
        self.dataSpec = False
        self.dataOrder = 'XRI'
        self.sw = 0.0
        self.delim = 'Tab'
        self.closed = False
        self.axisMulti = 1.0
        self.setWindowTitle("Load ASCII")
        grid = QtWidgets.QGridLayout(self)
        grid.addWidget(QtWidgets.QLabel("# Dimensions:"), 1, 0)
        self.numDims = wc.SsnakeSpinBox()
        self.numDims.setMinimum(1)
        self.numDims.setValue(1)
        self.numDims.setMaximum(2)
        grid.addWidget(self.numDims, 2, 0, 1, 2)
        grid.addWidget(QtWidgets.QLabel("Data Type:"), 3, 0)
        self.specGroup = QtWidgets.QButtonGroup(self)
        self.specGroup.buttonClicked.connect(self.checkState)
        self.timeButton = QtWidgets.QRadioButton('Time', parent=self)
        self.timeButton.toggle()
        self.specGroup.addButton(self.timeButton, 0)
        grid.addWidget(self.timeButton, 4, 0)
        self.freqButton = QtWidgets.QRadioButton('Frequency', parent=self)
        self.specGroup.addButton(self.freqButton, 1)
        grid.addWidget(self.freqButton, 4, 1)
        grid.addWidget(QtWidgets.QLabel("Data Order:"), 5, 0)
        self.datOrderBox = QtWidgets.QComboBox()
        self.datOrderBox.addItems(self.dataOrders)
        grid.addWidget(self.datOrderBox, 6, 0, 1, 2)
        self.unitLabel = wc.QLeftLabel("x-axis unit:")
        grid.addWidget(self.unitLabel, 7, 0, 1, 2)
        self.timeUnitBox = QtWidgets.QComboBox()
        self.timeUnitBox.addItems(self.timeUnits)
        grid.addWidget(self.timeUnitBox, 8, 0, 1, 2)
        self.freqUnitBox = QtWidgets.QComboBox()
        self.freqUnitBox.addItems(self.freqUnits)
        self.freqUnitBox.currentIndexChanged.connect(self.checkState)
        grid.addWidget(self.freqUnitBox, 8, 0, 1, 2)
        self.freqUnitBox.hide()
        self.swLabel = wc.QLeftLabel("Spectral Width [kHz]:")
        grid.addWidget(self.swLabel, 9, 0, 1, 2)
        self.swEntry = wc.QLineEdit("0.0")
        grid.addWidget(self.swEntry, 10, 0, 1, 2)
        self.swLabel.hide()
        self.swEntry.hide()
        self.datOrderBox.currentIndexChanged.connect(self.checkState)
        self.freqLabel = wc.QLeftLabel("Spectrometer freq. [MHz]:")
        grid.addWidget(self.freqLabel, 11, 0, 1, 2)
        self.freqEntry = wc.QLineEdit("0.0")
        grid.addWidget(self.freqEntry, 12, 0, 1, 2)
        self.freqLabel.hide()
        self.freqEntry.hide()
        grid.addWidget(QtWidgets.QLabel("Data Delimiter:"), 13, 0)
        self.datDelimBox = QtWidgets.QComboBox()
        self.datDelimBox.addItems(self.delimiters)
        grid.addWidget(self.datDelimBox, 14, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        box = QtWidgets.QDialogButtonBox()
        box.addButton(cancelButton, QtWidgets.QDialogButtonBox.RejectRole)
        box.addButton(okButton, QtWidgets.QDialogButtonBox.AcceptRole)
        grid.addWidget(box, 15, 0, 1, 2)
        self.show()
        self.setFixedSize(self.size())
        self.checkType(file)

    def checkState(self):
        tmp = self.dataOrders[self.datOrderBox.currentIndex()]
        if tmp in ('RI', 'R'):
            self.swLabel.show()
            self.swEntry.show()
            self.unitLabel.hide()
            self.timeUnitBox.hide()
            self.freqUnitBox.hide()
            self.freqLabel.hide()
            self.freqEntry.hide()
        else:
            self.swLabel.hide()
            self.swEntry.hide()
            self.unitLabel.show()
            if self.timeButton.isChecked():
                self.timeUnitBox.show()
                self.freqUnitBox.hide()
                self.freqLabel.hide()
                self.freqEntry.hide()
            else:
                self.timeUnitBox.hide()
                self.freqUnitBox.show()
                self.freqLabel.hide()
                self.freqEntry.hide()
                if self.freqUnitBox.currentIndex() == 3:
                    self.freqLabel.show()
                    self.freqEntry.show()

    def checkType(self, file):
        if file.endswith('.zip'):
            return # cannot read zipped files from here
        try:
            with open(file, 'r') as f:
                line = f.readline()
                while len(line)>0 and (line[0]=='#' or line[0]=='\n'):
                    line = f.readline()
            if line.count(',') > 0:
                sep = 'Comma'
            elif line.count('\t') > 0:
                sep = 'Tab'
            else:
                sep = 'Space'
            self.datDelimBox.setCurrentIndex(self.delimiters.index(sep))
            sepList = ['\t', ' ', ',']
            data = np.fromstring(line, sep=sepList[self.delimiters.index(sep)])
            if len(data) > 3:
                self.numDims.setValue(2)
        except Exception:
            return

    def closeEvent(self, *args):
        self.closed = True
        self.accept()
        self.deleteLater()

    def applyAndClose(self):
        self.dataOrder = self.dataOrders[self.datOrderBox.currentIndex()]
        self.delim = self.delimiters[self.datDelimBox.currentIndex()]
        if self.dataOrder == 'RI' or self.dataOrder == 'R':
            self.sw = safeEval(self.swEntry.text(), Type='FI')
            if self.sw == 0 or self.sw is None:
                raise SsnakeException('Spectral Width input is not valid')
        self.dataDimension = self.numDims.value()
        if self.timeButton.isChecked():
            self.dataSpec = False
        else:
            self.dataSpec = True

        if 'X' in self.dataOrder: # If there is an x-axis
            if self.dataSpec == False:
                self.axisMulti = self.timeMultiVals[self.timeUnitBox.currentIndex()]
            else:
                if self.freqUnitBox.currentIndex() == 3: #If ppm
                    self.axisMulti = safeEval(self.freqEntry.text())
                else:
                    self.axisMulti = self.freqMultiVals[self.freqUnitBox.currentIndex()]
        self.accept()
        self.deleteLater()

#################################################################################


class WorkInfoWindow(QtWidgets.QDialog):
    #A window to view info of the workspace (size etc)
    def __init__(self, parent):
        super(WorkInfoWindow, self).__init__(parent)
        self.father = parent
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool | QtCore.Qt.WindowContextHelpButtonHint)
        self.setWindowTitle("Workspace Info")
        grid = QtWidgets.QGridLayout(self)
        workGroup = QtWidgets.QGroupBox('Data:')
        workFrame = QtWidgets.QGridLayout()
        workGroup.setLayout(workFrame)
        workFrame.addWidget(QtWidgets.QLabel("Name:"), 0, 0)
        workFrame.addWidget(wc.QSelectLabel(self.father.masterData.name), 0, 1)
        workFrame.addWidget(QtWidgets.QLabel("# Dimensions:"), 1, 0)
        workFrame.addWidget(wc.QSelectLabel(str(self.father.masterData.ndim())), 1, 1)
        sw = self.father.masterData.sw
        np = self.father.masterData.shape()
        freq = self.father.masterData.freq
        ref = [x/1e6 if x is not None else x for x in self.father.masterData.ref]
        whole = self.father.masterData.wholeEcho
        spec = ['Time' if x == 0 else 'Frequency' for x in self.father.masterData.spec]
        for x in range(self.father.masterData.ndim()):
            workFrame.addWidget(wc.QSelectLabel('D' + str(x+1)), 2, x+1)
            workFrame.addWidget(wc.QSelectLabel(str(sw[x]/1000)), 3, x+1)
            workFrame.addWidget(wc.QSelectLabel(str(freq[x]/1e6)), 4, x+1)
            workFrame.addWidget(wc.QSelectLabel(str(ref[x])), 5, x+1)
            workFrame.addWidget(wc.QSelectLabel(str(np[x])), 6, x+1)
            workFrame.addWidget(wc.QSelectLabel(spec[x]), 7, x+1)
            workFrame.addWidget(wc.QSelectLabel(str(self.father.masterData.isComplex(x))), 8, x+1)
            workFrame.addWidget(wc.QSelectLabel(str(whole[x])), 9, x+1)
        workFrame.addWidget(QtWidgets.QLabel('Spectral Width [kHz]:'), 3, 0)
        workFrame.addWidget(QtWidgets.QLabel('Frequency [MHz]:'), 4, 0)
        workFrame.addWidget(QtWidgets.QLabel('Reference [MHz]:'), 5, 0)
        workFrame.addWidget(QtWidgets.QLabel('Number of Points:'), 6, 0)
        workFrame.addWidget(QtWidgets.QLabel('Type:'), 7, 0)
        workFrame.addWidget(QtWidgets.QLabel('Complex:'), 8, 0)
        workFrame.addWidget(QtWidgets.QLabel('Whole Echo:'), 9, 0)
        grid.addWidget(workGroup, 0, 0, 1, 3)
        metaGroup = QtWidgets.QGroupBox('Metadata:')
        metaFrame = QtWidgets.QGridLayout()
        metaGroup.setLayout(metaFrame)
        for pos, key in enumerate(self.father.masterData.metaData):
            metaFrame.addWidget(QtWidgets.QLabel(key), pos, 0)
            metaFrame.addWidget(wc.QSelectLabel(self.father.masterData.metaData[key]), pos, 1)
        grid.addWidget(metaGroup, 1, 0, 1, 3)
        okButton = QtWidgets.QPushButton("&Close")
        okButton.clicked.connect(self.closeEvent)
        grid.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())

    def closeEvent(self, *args):
        self.closed = True
        self.accept()
        self.deleteLater()

#################################################################################
class PhaseWindow(wc.ToolWindow):

    NAME = "Phasing"
    SINGLESLICE = True
    RESOLUTION = 1000
    P1LIMIT = 540.0
    P2LIMIT = 1440.0
    PHASE0STEP = 1.0
    PHASE1STEP = 1.0
    PHASE2STEP = 1.0

    def __init__(self, parent):
        super(PhaseWindow, self).__init__(parent)
        self.zeroVal = 0.0
        self.firstVal = 0.0
        self.secondVal = 0.0
        self.pivotFirstVal = 0.0
        self.pivotSecondVal = 0.0
        self.available = True

        # Zero order
        self.zeroOrderGroup = QtWidgets.QGroupBox('Zero order:')
        self.zeroOrderFrame = QtWidgets.QGridLayout()
        autoZero = QtWidgets.QPushButton("Autophase 0th")
        autoZero.clicked.connect(lambda: self.autophase(0))
        self.zeroOrderFrame.addWidget(autoZero, 0, 1)
        self.zeroEntry = wc.QLineEdit("0.000", self.inputZeroOrder)
        self.zeroOrderFrame.addWidget(self.zeroEntry, 2, 1)
        self.leftZero = QtWidgets.QPushButton("<")
        self.leftZero.clicked.connect(lambda: self.stepPhase(-1, 0, 0))
        self.leftZero.setAutoRepeat(True)
        self.zeroOrderFrame.addWidget(self.leftZero, 2, 0)
        self.rightZero = QtWidgets.QPushButton(">")
        self.rightZero.clicked.connect(lambda: self.stepPhase(1, 0, 0))
        self.rightZero.setAutoRepeat(True)
        self.zeroOrderFrame.addWidget(self.rightZero, 2, 2)
        self.zeroScale = wc.SsnakeSlider(QtCore.Qt.Horizontal)
        self.zeroScale.setRange(-self.RESOLUTION, self.RESOLUTION)
        self.zeroScale.valueChanged.connect(self.setZeroOrder)
        self.zeroOrderFrame.addWidget(self.zeroScale, 3, 0, 1, 3)
        self.zeroOrderGroup.setLayout(self.zeroOrderFrame)
        self.grid.addWidget(self.zeroOrderGroup, 0, 0, 1, 3)

        # First order
        self.firstOrderGroup = QtWidgets.QGroupBox('First order:')
        self.firstOrderFrame = QtWidgets.QGridLayout()
        autoFirst = QtWidgets.QPushButton("Autophase 0th+1st")
        autoFirst.clicked.connect(lambda: self.autophase(1))
        self.firstOrderFrame.addWidget(autoFirst, 0, 1)
        self.firstEntry = wc.QLineEdit("0.000", self.inputFirstOrder)
        self.firstOrderFrame.addWidget(self.firstEntry, 1, 1)
        self.leftFirst = QtWidgets.QPushButton("<")
        self.leftFirst.clicked.connect(lambda: self.stepPhase(0, -1, 0))
        self.leftFirst.setAutoRepeat(True)
        self.firstOrderFrame.addWidget(self.leftFirst, 1, 0)
        self.rightFirst = QtWidgets.QPushButton(">")
        self.rightFirst.clicked.connect(lambda: self.stepPhase(0, 1, 0))
        self.rightFirst.setAutoRepeat(True)
        self.firstOrderFrame.addWidget(self.rightFirst, 1, 2)
        self.firstScale = wc.SsnakeSlider(QtCore.Qt.Horizontal)
        self.firstScale.setRange(-self.RESOLUTION, self.RESOLUTION)
        self.firstScale.valueChanged.connect(self.setFirstOrder)
        self.firstOrderFrame.addWidget(self.firstScale, 2, 0, 1, 3)

        if self.father.current.spec() > 0:
            self.firstOrderFrame.addWidget(wc.QLabel("Pivot point [Hz]:"), 3, 0, 1, 3)
            pickFirstRef = QtWidgets.QPushButton("Pick pivot")
            pickFirstRef.clicked.connect(lambda: self.pickRef(1))
            self.firstOrderFrame.addWidget(pickFirstRef, 4, 1)
            self.refFirstEntry = wc.QLineEdit(('%.3f' % self.pivotFirstVal), lambda: self.inputRef(1))
            self.firstOrderFrame.addWidget(self.refFirstEntry, 5, 1)
        self.firstOrderGroup.setLayout(self.firstOrderFrame)
        self.grid.addWidget(self.firstOrderGroup, 1, 0, 1, 3)

        # Second order
        self.secondOrderGroup = QtWidgets.QGroupBox('Second order:')
        self.secondOrderFrame = QtWidgets.QGridLayout()
        self.secondEntry = wc.QLineEdit("0.000", self.inputSecondOrder)
        self.secondOrderFrame.addWidget(self.secondEntry, 0, 1)
        self.leftSecond = QtWidgets.QPushButton("<")
        self.leftSecond.clicked.connect(lambda: self.stepPhase(0, 0, -1))
        self.leftSecond.setAutoRepeat(True)
        self.secondOrderFrame.addWidget(self.leftSecond, 0, 0)
        self.rightSecond = QtWidgets.QPushButton(">")
        self.rightSecond.clicked.connect(lambda: self.stepPhase(0, 0, 1))
        self.rightSecond.setAutoRepeat(True)
        self.secondOrderFrame.addWidget(self.rightSecond, 0, 2)
        self.secondScale = wc.SsnakeSlider(QtCore.Qt.Horizontal)
        self.secondScale.setRange(-self.RESOLUTION, self.RESOLUTION)
        self.secondScale.valueChanged.connect(self.setSecondOrder)
        self.secondOrderFrame.addWidget(self.secondScale, 1, 0, 1, 3)

        if self.father.current.spec() > 0:
            self.secondOrderFrame.addWidget(wc.QLabel("Pivot point [Hz]:"), 2, 0, 1, 3)
            pickSecondRef = QtWidgets.QPushButton("Pick pivot")
            pickSecondRef.clicked.connect(lambda: self.pickRef(2))
            self.secondOrderFrame.addWidget(pickSecondRef, 3, 1)
            self.refSecondEntry = wc.QLineEdit(('%.3f' % self.pivotSecondVal), lambda: self.inputRef(2))
            self.secondOrderFrame.addWidget(self.refSecondEntry, 4, 1)
        self.secondOrderGroup.setLayout(self.secondOrderFrame)
        self.grid.addWidget(self.secondOrderGroup, 2, 0, 1, 3)
        self.secondOrderGroup.setVisible(self.father.father.defaultSecondOrderPhaseDialog)

        self.secondOrderCheckBox = QtWidgets.QCheckBox("2nd order phasing")
        self.secondOrderCheckBox.setChecked(self.father.father.defaultSecondOrderPhaseDialog)
        self.layout.addWidget(self.secondOrderCheckBox, 2, 0, 1, 3)
        self.secondOrderCheckBox.stateChanged.connect(self.setSecondOrderVisible)

    def setSecondOrderVisible(self):
        self.secondOrderGroup.setVisible(self.secondOrderCheckBox.isChecked())

    def setModifierTexts(self, event):
        sign = u"\u00D7"
        if event.modifiers() & QtCore.Qt.AltModifier:
            sign = '/'
        left = [self.leftZero, self.leftFirst, self.leftSecond]
        right = [self.rightZero, self.rightFirst, self.rightSecond]
        if event.modifiers() & QtCore.Qt.ControlModifier and event.modifiers() & QtCore.Qt.ShiftModifier:
            text = ' ' + sign + '1000'
        elif event.modifiers() & QtCore.Qt.ControlModifier:
            text = ' ' + sign + '10'
        elif event.modifiers() & QtCore.Qt.ShiftModifier:
            text = ' ' + sign + '100'
        else:
            text = ''
        for widget in left:
            widget.setText('<' + text)
        for widget in right:
            widget.setText('>' + text)

    def keyPressEvent(self, event):
        self.setModifierTexts(event)

    def keyReleaseEvent(self, event):
        self.setModifierTexts(event)

    def setZeroOrder(self, value, *args):
        if self.available:
            self.zeroVal = float(value) / self.RESOLUTION * 180
            self.zeroEntry.setText('%.3f' % self.zeroVal)
            self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, np.pi * self.secondVal / 180.0)

    def inputZeroOrder(self, *args):
        inp = safeEval(self.zeroEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException('Phasing: zero order value input is not valid!')
        self.zeroVal = np.mod(inp + 180, 360) - 180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.available = False
        self.zeroScale.setValue(int(round(self.zeroVal / 180.0 * self.RESOLUTION)))
        self.available = True
        self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, np.pi * self.secondVal / 180.0)

    def setFirstOrder(self, value, *args):
        if self.available:
            value = float(value) / self.RESOLUTION * self.P1LIMIT
            newZero = (self.zeroVal - (value - self.firstVal) * self.pivotFirstVal / self.father.current.sw())
            self.zeroVal = np.mod(newZero + 180, 360) - 180
            self.zeroEntry.setText('%.3f' % self.zeroVal)
            self.firstVal = value
            self.firstEntry.setText('%.3f' % self.firstVal)
            self.available = False
            self.zeroScale.setValue(int(round(self.zeroVal / 180.0 * self.RESOLUTION)))
            self.available = True
            self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, np.pi * self.secondVal / 180.0)

    def inputFirstOrder(self, *args):
        value = safeEval(self.firstEntry.text(), length=self.father.current.len(), Type='FI')
        if value is None:
            raise SsnakeException('Phasing: first order value input is not valid!')
        newZero = (self.zeroVal - (value - self.firstVal) * self.pivotFirstVal / self.father.current.sw())
        self.zeroVal = np.mod(newZero + 180, 360) - 180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.firstVal = value
        self.firstEntry.setText('%.3f' % self.firstVal)
        self.available = False
        self.zeroScale.setValue(int(round(self.zeroVal / 180.0 * self.RESOLUTION)))
        self.firstScale.setValue(int(round(self.firstVal / self.P1LIMIT * self.RESOLUTION)))
        self.available = True
        self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, np.pi * self.secondVal / 180.0)

    def setSecondOrder(self, value, *args):
        if self.available:
            value = float(value) / self.RESOLUTION * self.P2LIMIT
            pivot1 = self.pivotFirstVal / self.father.current.sw()
            pivot2 = self.pivotSecondVal / self.father.current.sw()
            if pivot1 == pivot2:
                newFirst = (self.firstVal - 2 * (value - self.secondVal) * pivot1)
            else:
                newFirst = (self.firstVal + (value - self.secondVal) * (np.power(pivot2, 2) - np.power(pivot1, 2)) / (pivot1 -pivot2))
            newZero = (self.zeroVal - (newFirst - self.firstVal) * pivot1 - (value - self.secondVal) * np.power(pivot1, 2))
            self.zeroVal = np.mod(newZero + 180, 360) - 180
            self.zeroEntry.setText('%.3f' % self.zeroVal)
            self.firstVal = newFirst
            self.firstEntry.setText('%.3f' % newFirst)
            self.secondVal = value
            self.secondEntry.setText('%.3f' % self.secondVal)
            self.available = False
            self.zeroScale.setValue(int(round(self.zeroVal / 180.0 * self.RESOLUTION)))
            self.firstScale.setValue(int(round(self.firstVal / self.P1LIMIT * self.RESOLUTION)))
            self.available = True
            self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, np.pi * self.secondVal / 180.0)

    def inputSecondOrder(self, *args):
        value = safeEval(self.secondEntry.text(), length=self.father.current.len(), Type='FI')
        if value is None:
            raise SsnakeException('Phasing: second order value input is not valid!')
        pivot1 = self.pivotFirstVal / self.father.current.sw()
        pivot2 = self.pivotSecondVal / self.father.current.sw()
        if pivot1 == pivot2:
            newFirst = (self.firstVal - 2 * (value - self.secondVal) * pivot1)
        else:
            newFirst = (self.firstVal + (value - self.secondVal) * (np.power(pivot2, 2) - np.power(pivot1, 2)) / (pivot1 - pivot2))
        newZero = (self.zeroVal - (newFirst - self.firstVal) * pivot1 - (value - self.secondVal) * np.power(pivot1, 2))
        self.zeroVal = np.mod(newZero + 180, 360) - 180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.firstVal = newFirst
        self.firstEntry.setText('%.3f' % newFirst)
        self.secondVal = value
        self.secondEntry.setText('%.3f' % self.secondVal)
        self.available = False
        self.zeroScale.setValue(int(round(self.zeroVal / 180.0 * self.RESOLUTION)))
        self.firstScale.setValue(int(round(self.firstVal / self.P1LIMIT * self.RESOLUTION)))
        self.secondScale.setValue(int(round(self.secondVal / self.P2LIMIT * self.RESOLUTION)))
        self.available = True
        self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, np.pi * self.secondVal / 180.0)

    def autophase(self, num):
        phases = self.father.current.autoPhase(num)
        val = phases[0] / np.pi * 180.0
        self.zeroVal = (np.mod(val + 180, 360) - 180)
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.available = False
        self.zeroScale.setValue(round(self.zeroVal / 180.0 * self.RESOLUTION))
        self.available = True
        if num == 1:
            val = phases[1] / np.pi * 180.0
            self.firstVal = val
            self.firstEntry.setText('%.3f' % self.firstVal)
        self.inputFirstOrder()

    def stepPhase(self, phase0, phase1, phase2):
        step = 1
        multiplier = 1
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier and QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            multiplier *= 1000
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            multiplier *= 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            multiplier *= 100
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.AltModifier:
            step = step / multiplier
        else:
            step = step * multiplier
        phase0 = step * phase0
        phase1 = step * phase1
        phase2 = step * phase2
        inp = safeEval(self.zeroEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException('Phasing: zero order value input is not valid!')
        inp += phase0 * self.PHASE0STEP
        self.zeroVal = np.mod(inp + 180, 360) - 180
        value = safeEval(self.firstEntry.text(), length=self.father.current.len(), Type='FI')
        if value is None:
            raise SsnakeException('Phasing: first order value input is not valid!')
        value += phase1 * self.PHASE1STEP
        if self.father.current.spec() > 0:
            self.inputRef(1)
        second = safeEval(self.secondEntry.text(), length=self.father.current.len(), Type='FI')
        if second is None:
            raise SsnakeException('Phasing: second order value input is not valid!')
        second += phase2 * self.PHASE2STEP
        if self.father.current.spec() > 0:
            self.inputRef(2)
        pivot1 = self.pivotFirstVal / self.father.current.sw()
        pivot2 = self.pivotSecondVal / self.father.current.sw()
        if pivot1 == pivot2:
            newFirst = (value - 2 * (second - self.secondVal) * pivot1)
        else:
            newFirst = (value + (second - self.secondVal) * (np.power(pivot2, 2) - np.power(pivot1, 2)) / (pivot1 - pivot2))
        newZero = (self.zeroVal - (newFirst - self.firstVal) * pivot1 - (second - self.secondVal) * np.power(pivot1, 2))
        self.zeroVal = np.mod(newZero + 180, 360) - 180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.firstVal = newFirst
        self.firstEntry.setText('%.3f' % self.firstVal)
        self.secondVal = second
        self.secondEntry.setText('%.3f' % self.secondVal)
        self.available = False
        self.zeroScale.setValue(round(self.zeroVal / 180.0 * self.RESOLUTION))
        self.firstScale.setValue(round(self.firstVal / self.P1LIMIT * self.RESOLUTION))
        self.secondScale.setValue(round(self.secondVal / self.P2LIMIT * self.RESOLUTION))
        self.available = True
        self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, np.pi * self.secondVal / 180.0)

    def inputRef(self, order, *args):
        if order == 1:
            Val = safeEval(self.refFirstEntry.text(), length=self.father.current.len(), Type='FI')
            if Val is None:
                raise SsnakeException('Phasing: pivot input is not valid!')
            self.pivotFirstVal = Val
            self.refFirstEntry.setText('%.3f' % self.pivotFirstVal)
        elif order == 2:
            Val = safeEval(self.refSecondEntry.text(), length=self.father.current.len(), Type='FI')
            if Val is None:
                raise SsnakeException('Phasing: pivot input is not valid!')
            self.pivotSecondVal = Val
            self.refSecondEntry.setText('%.3f' % self.pivotSecondVal)

    def setRef(self, value, order, *args):
        if order == 1:
            self.pivotFirstVal = float(value)
            self.refFirstEntry.setText('%.3f' % self.pivotFirstVal)
        elif order == 2:
            self.pivotSecondVal = float(value)
            self.refSecondEntry.setText('%.3f' % self.pivotSecondVal)

    def pickRef(self, order, *args):
        if order == 1:
            self.father.current.peakPickFunc = lambda pos, self=self: self.setRef(self.father.current.xax()[pos[0]], 1)
        elif order == 2:
            self.father.current.peakPickFunc = lambda pos, self=self: self.setRef(self.father.current.xax()[pos[0]], 2)
        self.father.current.peakPick = True

    def applyFunc(self):
        if self.father.current.spec() > 0:
            self.inputRef(1)
            if self.secondOrderCheckBox.isChecked():
                self.inputRef(2)
        self.inputZeroOrder()
        self.inputFirstOrder()
        if self.secondOrderCheckBox.isChecked():
            self.inputSecondOrder()
        self.father.current.applyPhase(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, np.pi * self.secondVal / 180.0, (self.singleSlice.isChecked() == 1))

################################################################


class ApodWindow(wc.ToolWindow):

    RESOLUTION = 10000
    NAME = "Apodize"
    SINGLESLICE = True

    def __init__(self, parent):
        super(ApodWindow, self).__init__(parent)
        self.entries = {}
        self.ticks = {}
        boldFont = QtGui.QFont()
        boldFont.setBold(True)
        self.maximum = 100.0 * self.father.current.sw() / (self.father.current.len())
        self.lbstep = 1.0
        self.available = True
        self.lorGroup = QtWidgets.QGroupBox()
        self.lorFrame = QtWidgets.QGridLayout()
        lorTick = QtWidgets.QCheckBox("Lorentzian [Hz]:")
        lorTick.setFont(boldFont)
        lorTick.toggled.connect(lambda x: self.checkEval('lor'))
        self.lorFrame.addWidget(lorTick, 0, 0, 1, 3)
        self.ticks['lor'] = lorTick
        lorEntry = wc.QLineEdit("0.00", self.apodPreview)
        lorEntry.setMinimumWidth(150)
        lorEntry.setEnabled(False)
        self.lorFrame.addWidget(lorEntry, 1, 1)
        self.entries['lor'] = [lorEntry]
        self.leftLor = QtWidgets.QPushButton("<")
        self.leftLor.clicked.connect(lambda: self.stepLB(-0.5 * self.father.current.sw() / (self.father.current.len()), 'lor'))
        self.leftLor.setAutoRepeat(True)
        self.lorFrame.addWidget(self.leftLor, 1, 0)
        self.rightLor = QtWidgets.QPushButton(">")
        self.rightLor.clicked.connect(lambda: self.stepLB(0.5 * self.father.current.sw() / (self.father.current.len()), 'lor'))
        self.rightLor.setAutoRepeat(True)
        self.lorFrame.addWidget(self.rightLor, 1, 2)
        self.lorScale = wc.SsnakeSlider(QtCore.Qt.Horizontal)
        self.lorScale.setRange(0, self.RESOLUTION)
        self.lorScale.valueChanged.connect(lambda x: self.setLorGauss(x, 'lor'))
        self.lorFrame.addWidget(self.lorScale, 2, 0, 1, 3)
        self.lorMax = 100.0 * self.father.current.sw() / (self.father.current.len())
        self.lorGroup.setLayout(self.lorFrame)
        self.grid.addWidget(self.lorGroup, 0, 0, 1, 3)
        self.gaussGroup = QtWidgets.QGroupBox()
        self.gaussFrame = QtWidgets.QGridLayout()
        gaussTick = QtWidgets.QCheckBox("Gaussian [Hz]:")
        gaussTick.setFont(boldFont)
        gaussTick.toggled.connect(lambda: self.checkEval('gauss'))
        self.gaussFrame.addWidget(gaussTick, 3, 0, 1, 3)
        self.ticks['gauss'] = gaussTick
        gaussEntry = wc.QLineEdit("0.00", self.apodPreview)
        gaussEntry.setEnabled(False)
        gaussEntry.setMinimumWidth(150)
        self.gaussFrame.addWidget(gaussEntry, 4, 1)
        self.entries['gauss'] = [gaussEntry]
        self.leftGauss = QtWidgets.QPushButton("<")
        self.leftGauss.clicked.connect(lambda: self.stepLB(-0.5 * self.father.current.sw() / (self.father.current.len()), 'gauss'))
        self.leftGauss.setAutoRepeat(True)
        self.gaussFrame.addWidget(self.leftGauss, 4, 0)
        self.rightGauss = QtWidgets.QPushButton(">")
        self.rightGauss.clicked.connect(lambda: self.stepLB(0.5 * self.father.current.sw() / (self.father.current.len()), 'gauss'))
        self.rightGauss.setAutoRepeat(True)
        self.gaussFrame.addWidget(self.rightGauss, 4, 2)
        self.gaussScale = wc.SsnakeSlider(QtCore.Qt.Horizontal)
        self.gaussScale.setRange(0, self.RESOLUTION)
        self.gaussScale.valueChanged.connect(lambda x: self.setLorGauss(x, 'gauss'))
        self.gaussFrame.addWidget(self.gaussScale, 5, 0, 1, 3)
        self.gaussMax = 100.0 * self.father.current.sw() / (self.father.current.len())
        self.gaussGroup.setLayout(self.gaussFrame)
        self.grid.addWidget(self.gaussGroup, 1, 0, 1, 3)
        self.cos2Group = QtWidgets.QGroupBox()
        self.cos2Frame = QtWidgets.QGridLayout()
        cos2Tick = QtWidgets.QCheckBox("Cos^2:")
        cos2Tick.setFont(boldFont)
        cos2Tick.clicked.connect(lambda: self.checkEval('cos2'))
        self.cos2Frame.addWidget(cos2Tick, 0, 0, 1, 2)
        self.ticks['cos2'] = cos2Tick
        cos2Entry = wc.QLineEdit("1", self.apodPreview)
        cos2Entry.setEnabled(False)
        widthHint = cos2Entry.minimumSizeHint()
        widthHint.setWidth(widthHint.width() *4)
        cos2Entry.setMinimumSize(widthHint)
        self.cos2Frame.addWidget(cos2Entry, 1, 2)
        self.entries['cos2'] = [cos2Entry]
        cos2Label = wc.QLeftLabel("Frequency:")
        cos2Label.setEnabled(False)
        self.entries['cos2'].append(cos2Label)
        self.cos2Frame.addWidget(cos2Label, 1, 0)
        cos2DegEntry = wc.QLineEdit("0", self.apodPreview)
        cos2DegEntry.setEnabled(False)
        cos2DegEntry.setMinimumSize(widthHint)
        self.cos2Frame.addWidget(cos2DegEntry, 2, 2)
        self.entries['cos2'].append(cos2DegEntry)
        cos2PhLabel = wc.QLeftLabel("Phase [deg]:")
        cos2PhLabel.setEnabled(False)
        self.entries['cos2'].append(cos2PhLabel)
        self.cos2Frame.addWidget(QtWidgets.QWidget(), 1, 1)
        self.cos2Frame.setColumnStretch(1, 1)
        self.cos2Frame.addWidget(cos2PhLabel, 2, 0)
        self.cos2Group.setLayout(self.cos2Frame)
        self.grid.addWidget(self.cos2Group, 2, 0, 1, 3)
        self.hammingGroup = QtWidgets.QGroupBox()
        self.hammingFrame = QtWidgets.QGridLayout()
        hammingTick = QtWidgets.QCheckBox("Hamming:")
        hammingTick.setFont(boldFont)
        hammingTick.clicked.connect(lambda: self.checkEval('hamming'))
        self.hammingFrame.addWidget(hammingTick, 0, 0, 1, 2)
        self.ticks['hamming'] = hammingTick
        hammingLabel = wc.QLeftLabel("Frequency:\t")
        hammingLabel.setEnabled(False)
        self.entries['hamming'] = [hammingLabel]
        self.hammingFrame.addWidget(hammingLabel, 1, 0)
        hammingEntry = wc.QLineEdit("1", self.apodPreview)
        hammingEntry.setEnabled(False)
        hammingEntry.setMinimumSize(widthHint)
        self.hammingFrame.addWidget(hammingEntry, 1, 2)
        self.entries['hamming'].append(hammingEntry)
        self.hammingFrame.addWidget(QtWidgets.QWidget(), 1, 1)
        self.hammingFrame.setColumnStretch(1, 1)
        self.hammingGroup.setLayout(self.hammingFrame)
        self.grid.addWidget(self.hammingGroup, 3, 0, 1, 3)
        self.shiftGroup = QtWidgets.QGroupBox()
        self.shiftFrame = QtWidgets.QGridLayout()
        shiftTick = QtWidgets.QCheckBox("Shift:")
        shiftTick.setFont(boldFont)
        shiftTick.clicked.connect(lambda: self.checkEval('shift'))
        self.ticks['shift'] = shiftTick
        self.shiftFrame.addWidget(shiftTick, 0, 0)
        shiftLabel = wc.QLeftLabel("Value [s]:")
        shiftLabel.setEnabled(False)
        self.shiftFrame.addWidget(shiftLabel, 1, 0)
        shiftEntry = wc.QLineEdit("0", self.apodPreview)
        self.shiftFrame.addWidget(QtWidgets.QWidget(), 1, 1)
        shiftEntry.setMinimumSize(widthHint)
        shiftEntry.setEnabled(False)
        self.entries['shift'] = [shiftEntry, shiftLabel]
        self.shiftFrame.addWidget(shiftEntry, 1, 2)
        self.shiftGroup.setLayout(self.shiftFrame)
        self.shiftFrame.setColumnStretch(1, 1)
        self.grid.addWidget(self.shiftGroup, 4, 0, 1, 3)
        if self.father.current.data.ndim() > 1:
            self.shiftingGroup = QtWidgets.QGroupBox()
            self.shiftingFrame = QtWidgets.QGridLayout()
            shiftingTick = QtWidgets.QCheckBox("Shifting:")
            shiftingTick.clicked.connect(lambda: self.checkEval('shifting'))
            shiftingTick.setFont(boldFont)
            self.ticks['shifting'] = shiftingTick
            self.shiftingFrame.addWidget(shiftingTick, 0, 0)
            self.shiftingDropdown = QtWidgets.QComboBox()
            drop_list = ['User Defined', 
                     'Spin 3/2, 3QMAS', 'Spin 3/2, ST1MAS', 'Spin 3/2, DQ-STMAS',
                     'Spin 5/2, 3QMAS', 'Spin 5/2, 5QMAS', 'Spin 5/2, ST1MAS', 'Spin 5/2, DQ-STMAS',
                     'Spin 7/2, 3QMAS', 'Spin 5/2, 5QMAS', 'Spin 5/2, 7QMAS', 'Spin 7/2, ST1MAS', 'Spin 7/2, DQ-STMAS',
                     'Spin 9/2, 3QMAS', 'Spin 9/2, 5QMAS', 'Spin 9/2, 7QMAS', 'Spin 9/2, 9QMAS', 'Spin 9/2, ST1MAS', 'Spin 9/2, DQ-STMAS',
                    ]
            self.shiftingDropdown.addItems(drop_list)
            #['User Defined', 'Spin 3/2, -3Q (7/9)', 'Spin 5/2, 3Q (19/12)',
            #                                'Spin 5/2, -5Q (25/12)', 'Spin 7/2, 3Q (101/45)', 'Spin 7/2, 5Q (11/9)',
            #                                'Spin 7/2, -7Q (161/45)', 'Spin 9/2, 3Q (91/36)', 'Spin 9/2, 5Q (95/36)',
            #                                'Spin 9/2, 7Q (7/18)', 'Spin 9/2, -9Q (31/6)'])
            self.shiftingDropdown.activated.connect(self.dropdownChanged)
            #self.shiftingList = [0, 7.0 / 9.0, 19.0 / 12.0,
            #                     25.0 / 12.0, 101.0 / 45.0, 11.0 / 9.0,
            #                     161.0 / 45.0, 91.0 / 36.0, 95.0 / 36.0,
            #                     7.0 / 18.0, 31.0 / 6.0]
            self.shiftingList = [0,
                          (3/2, -3/2, 3/2), (3/2, -3/2, -1/2), (3/2, -3/2, 1/2),
                          (5/2, -3/2, 3/2), (5/2, -5/2, 5/2), (5/2, -3/2, -1/2), (5/2, -3/2, 1/2),
                          (7/2, -3/2, 3/2), (7/2, -5/2, 5/2), (7/2, -7/2, 7/2), (7/2, -3/2, -1/2), (7/2, -3/2, 1/2),
                          (9/2, -3/2, 3/2), (9/2, -5/2, 5/2), (9/2, -7/2, 7/2), (9/2, -9/2, 9/2), (9/2, -3/2, -1/2), (9/2, -3/2, 1/2),
                        ]
            self.shiftingDropdown.setMinimumSize(widthHint)
            self.shiftingDropdown.setEnabled(False)
            self.shiftingFrame.addWidget(self.shiftingDropdown, 1, 2)
            shiftingTypeLabel = wc.QLeftLabel("Type:")
            shiftingTypeLabel.setEnabled(False)
            self.shiftingFrame.addWidget(shiftingTypeLabel, 1, 0)
            self.shiftingEntry = wc.QLineEdit("0.00", self.apodPreview)
            self.shiftingEntry.setEnabled(False)
            self.shiftingEntry.setMinimumSize(widthHint)
            self.shiftingFrame.addWidget(self.shiftingEntry, 2, 2)
            shiftingValueLabel = wc.QLeftLabel("Value:")
            shiftingValueLabel.setEnabled(False)
            self.shiftingFrame.addWidget(shiftingValueLabel, 2, 0)
            self.shiftingAxis = QtWidgets.QComboBox()
            self.shiftingValues = list(map(str, np.delete(range(1, self.father.current.data.ndim() + 1), self.father.current.axes[-1])))
            self.shiftingAxis.addItems(self.shiftingValues)
            self.shiftingAxis.currentIndexChanged.connect(self.apodPreview)
            self.shiftingAxis.setMinimumSize(widthHint)
            self.shiftingAxis.setEnabled(False)
            self.shiftingFrame.addWidget(self.shiftingAxis, 3, 2)
            shiftingAxisLabel = wc.QLeftLabel("Axis:")
            shiftingAxisLabel.setEnabled(False)
            self.shiftingFrame.addWidget(shiftingAxisLabel, 3, 0)
            self.shiftingFrame.addWidget(QtWidgets.QWidget(), 1, 1)
            self.shiftingFrame.setColumnStretch(1, 1)
            self.entries['shifting'] = [self.shiftingDropdown, shiftingTypeLabel, self.shiftingEntry, shiftingValueLabel, self.shiftingAxis, shiftingAxisLabel]
            self.shiftingGroup.setLayout(self.shiftingFrame)
            self.grid.addWidget(self.shiftingGroup, 5, 0, 1, 3)

    def setModifierTexts(self, event):
        sign = u"\u00D7"
        if event.modifiers() & QtCore.Qt.AltModifier:
            sign = '/'
        left = [self.leftLor, self.leftGauss]
        right = [self.rightLor, self.rightGauss]
        if event.modifiers() & QtCore.Qt.ControlModifier and event.modifiers() & QtCore.Qt.ShiftModifier:
            text = ' ' + sign + '1000'
        elif event.modifiers() & QtCore.Qt.ControlModifier:
            text = ' ' + sign + '10'
        elif event.modifiers() & QtCore.Qt.ShiftModifier:
            text = ' ' + sign + '100'
        else:
            text = ''
        for widget in left:
            widget.setText('<' + text)
        for widget in right:
            widget.setText('>' + text)

    def keyPressEvent(self, event):
        self.setModifierTexts(event)

    def keyReleaseEvent(self, event):
        self.setModifierTexts(event)

    def dropdownChanged(self, update=True):
        index = self.shiftingDropdown.currentIndex()
        if index == 0:
            shifting = "0"
        else:
            shifting = f"{abs(func.R(*self.shiftingList[index]))}"
#        if index == 0:
#            self.shiftingEntry.setEnabled(True)
#        else:
#            self.shiftingEntry.setEnabled(False)
        if update:
            self.shiftingEntry.setText(shifting)
            self.apodPreview()

    def checkEval(self, key):
        if self.ticks[key].isChecked():
            for elem in self.entries[key]:
                elem.setEnabled(True)
        else:
            for elem in self.entries[key]:
                elem.setEnabled(False)
        if self.father.current.data.ndim() > 1:
            if self.ticks['shifting'].isChecked():
                self.dropdownChanged(update=False) #Check dropdown state
        if key in ('lor', 'gauss'):  # for lorentzian and gaussian
            if safeEval(self.entries[key][0].text(), length=self.father.current.len(), Type='FI') != 0.0:  # only update if value was not zero
                self.apodPreview()
        else:
            self.apodPreview()

    def setLorGauss(self,value, type, *args):
        #type: 'lor' or 'gauss'
        if self.available:
            self.entries[type][0].setText(('%.'+str(self.father.father.defaultPrecis)+'g') % (float(value) * self.maximum / self.RESOLUTION))
            if not self.ticks[type].isChecked():
                self.ticks[type].setChecked(1)
            self.apodPreview()

    def apodPreview(self, *args):
        self.available = False #turn off gauss/lorentz callbacks
        lor, gauss, cos2, cos2Ph, hamming, shift, shifting, shiftingAxis = self.checkInput()
        self.available = True
        self.father.current.apodPreview(lor, gauss, [cos2, cos2Ph], hamming, shift, shifting, shiftingAxis)

    def stepLB(self, incr, type):
        step = incr * self.lbstep
        multiplier = 1
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier and QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            multiplier = 1000
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            multiplier = 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            multiplier = 100
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.AltModifier:
            step = step / multiplier
        else:
            step = step * multiplier
        if not self.ticks[type].isChecked():
            self.ticks[type].setChecked(1)
        lor, gauss, cos2, cos2Ph, hamming, shift, shifting, shiftingAxis = self.checkInput()
        if type == 'lor':
            self.entries[type][0].setText(('%.'+str(self.father.father.defaultPrecis)+'g') % (lor + step))
        elif type == 'gauss':
            self.entries[type][0].setText(('%.'+str(self.father.father.defaultPrecis)+'g') % (gauss + step))
        self.apodPreview()

    def checkInput(self):
        lor = None
        gauss = None
        cos2 = None
        cos2Ph = None
        hamming = None
        shift = 0.0
        shifting = 0.0
        shiftingAxis = None
        if self.ticks['lor'].isChecked():
            lor = safeEval(self.entries['lor'][0].text(), length=self.father.current.len(), Type='FI')
            if lor is None:
                self.father.current.showFid()
                raise SsnakeException('Apodize: Lorentzian value is not valid!')
            self.lorScale.setValue(int(round(lor * self.RESOLUTION / self.maximum)))
        if self.ticks['gauss'].isChecked():
            gauss = safeEval(self.entries['gauss'][0].text(), length=self.father.current.len(), Type='FI')
            if gauss is None:
                self.father.current.showFid()
                raise SsnakeException('Apodize: Gaussian value is not valid!')
            self.gaussScale.setValue(int(round(gauss * self.RESOLUTION / self.maximum)))
        if self.ticks['cos2'].isChecked():
            cos2 = safeEval(self.entries['cos2'][0].text(), length=self.father.current.len(), Type='FI')
            if cos2 is None:
                self.father.current.showFid()
                raise SsnakeException('Apodize: cos^2 frequency value is not valid!')
        if self.ticks['cos2'].isChecked():
            cos2Ph = safeEval(self.entries['cos2'][2].text(), length=self.father.current.len(), Type='FI')
            if cos2Ph is None:
                self.father.current.showFid()
                raise SsnakeException('Apodize: cos^2 phase value is not valid!')
        if self.ticks['hamming'].isChecked():
            hamming = safeEval(self.entries['hamming'][1].text(), length=self.father.current.len(), Type='FI')
            if hamming is None:
                self.father.current.showFid()
                raise SsnakeException('Apodize: Hamming value is not valid!')
        if self.ticks['shift'].isChecked():
            shift = safeEval(self.entries['shift'][0].text(), length=self.father.current.len(), Type='FI')
            if shift is None:
                self.father.current.showFid()
                raise SsnakeException('Apodize: Shift value is not valid!')
        if self.father.current.data.ndim() > 1:
            if self.ticks['shifting'].isChecked():
                if self.father.current.data.ndim() > 1:
                    shifting = safeEval(self.shiftingEntry.text(), length=self.father.current.len(), Type='FI')
                    if shifting is None:
                        self.father.current.showFid()
                        raise SsnakeException('Apodize: Shifting value is not valid!')
                    shiftingAxis = int(self.shiftingValues[self.shiftingAxis.currentIndex()]) - 1
                else:
                    shiftingAxis = None
        return lor, gauss, cos2, cos2Ph, hamming, shift, shifting, shiftingAxis

    def applyFunc(self):
        lor, gauss, cos2, cos2Ph, hamming, shift, shifting, shiftingAxis = self.checkInput()
        self.father.current.applyApod(lor, gauss, [cos2, cos2Ph], hamming, shift, shifting, shiftingAxis, (self.singleSlice.isChecked()))

#######################################################################################

class SizeWindow(wc.ToolWindow):

    NAME = "Set size"

    def __init__(self, parent):
        super(SizeWindow, self).__init__(parent)
        self.sizeGroup = QtWidgets.QGroupBox('Size:')
        self.sizeFrame = QtWidgets.QGridLayout()
        self.sizeVal = parent.current.len()
        self.sizeEntry = wc.QLineEdit(self.sizeVal, self.sizePreview)
        self.sizeEntry.setMinimumWidth(100)
        self.sizeFrame.addWidget(self.sizeEntry, 0, 1)
        rightPower = QtWidgets.QPushButton("+ 2^n")
        rightPower.clicked.connect(lambda: self.stepSize(True))
        # rightZero.setAutoRepeat(True)
        self.sizeFrame.addWidget(rightPower, 0, 2)
        leftPower = QtWidgets.QPushButton("- 2^n")
        leftPower.clicked.connect(lambda: self.stepSize(False))
        self.sizeFrame.addWidget(leftPower, 0, 0)
        self.sizeGroup.setLayout(self.sizeFrame)
        self.grid.addWidget(self.sizeGroup, 0, 0, 1, 3)
        # offset
        self.offGroup = QtWidgets.QGroupBox('Offset:')
        self.offFrame = QtWidgets.QGridLayout()
        if self.father.current.wholeEcho():
            self.posVal = int(np.floor(parent.current.len() / 2.0))
        else:
            self.posVal = parent.current.len()
        self.posEntry = wc.QLineEdit(self.posVal, self.sizePreview)
        self.offFrame.addWidget(self.posEntry, 0, 1)
        self.offGroup.setLayout(self.offFrame)
        self.grid.addWidget(self.offGroup, 1, 0, 1, 3)
        if not self.father.current.spec():
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
            self.father.current.peakPick = True

    def stepSize(self, forward):
        inp = safeEval(self.sizeEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException('Sizing: \'Size\' input is not valid')
        inp = int(round(inp))
        if inp < 1:
            raise SsnakeException('Sizing: \'Size\' cannot be below 1')
        if forward:  # If + button
            new = int(np.floor(np.log2(inp)) + 1)
        else:
            new = int(np.ceil(np.log2(inp)) - 1)
        if new < 0:
            new = 0
        self.sizeEntry.setText(str(2**new))
        self.sizePreview()

    def sizePreview(self, *args):
        inp = safeEval(self.sizeEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException('Sizing: \'Size\' input is not valid')
        self.sizeVal = int(round(inp))
        if self.sizeVal < 1:
            raise SsnakeException('Sizing: \'Size\' cannot be below 1')
        self.sizeEntry.setText(str(self.sizeVal))
        inp = safeEval(self.posEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException('Sizing: \'Offset\' input is not valid')
        self.posVal = int(round(inp))
        if self.posVal < 1:
            raise SsnakeException('Sizing: \'Offset\' cannot be below 1')
        self.posEntry.setText(str(self.posVal))
        self.father.current.resizePreview(self.sizeVal, self.posVal)

    def applyFunc(self):
        inp = safeEval(self.sizeEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException('Sizing: \'Size\' input is not valid')
        self.sizeVal = int(round(inp))
        if self.sizeVal < 1:
            raise SsnakeException('Sizing: \'Size\' cannot be below 1')
        inp = safeEval(self.posEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException('Sizing: \'Offset\' input is not valid')
        self.posVal = int(round(inp))
        if self.posVal < 1:
            raise SsnakeException('Sizing: \'Offset\' cannot be below 1')
        self.father.current.resize(self.sizeVal, self.posVal)
        self.father.sideframe.upd()

    def picked(self, pos):
        self.posEntry.setText(str(pos[0]))
        self.sizePreview()
        self.father.current.peakPick = True
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)

##########################################################################################


class SwapEchoWindow(wc.ToolWindow):

    NAME = "Swap echo"

    def __init__(self, parent):
        super(SwapEchoWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Echo position:"), 0, 0)
        self.posVal = int(round(0.5 * parent.current.len()))
        self.posEntry = wc.QLineEdit(self.posVal, self.swapEchoPreview)
        self.grid.addWidget(self.posEntry, 1, 0)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def swapEchoPreview(self, *args):
        inp = safeEval(self.posEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Swap echo: not a valid index")
        self.posVal = int(round(inp))
        self.posEntry.setText(str(self.posVal))
        if self.posVal > 0 and self.posVal < self.father.current.len():
            self.father.current.swapEchoPreview(self.posVal)
            self.father.current.peakPick = False

    def applyFunc(self):
        self.father.current.peakPickReset()
        inp = safeEval(self.posEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Swap echo: not a valid index")
        self.posVal = int(round(inp))
        self.posEntry.setText(str(self.posVal))
        if self.posVal < 0 or self.posVal >= (self.father.current.len()):
            raise SsnakeException("Swap echo: not a valid index")
        self.father.current.swapEcho(self.posVal)
        self.father.bottomframe.upd()

    def picked(self, pos):
        self.father.current.swapEchoPreview(pos[0])
        self.posEntry.setText(str(pos[0]))
        self.father.current.peakPick = False

###########################################################################


class LPSVDWindow(wc.ToolWindow):

    NAME = "LPSVD"

    def __init__(self, parent):
        super(LPSVDWindow, self).__init__(parent)
        self.specGroup = QtWidgets.QButtonGroup(self)
        backwardButton = QtWidgets.QRadioButton('Backward', parent=self)
        self.specGroup.addButton(backwardButton, 1)
        forwardButton = QtWidgets.QRadioButton('Forward', parent=self)
        self.specGroup.addButton(forwardButton, 0)
        self.grid.addWidget(backwardButton, 1, 0)
        self.grid.addWidget(forwardButton, 2, 0)
        backwardButton.setChecked(True)
        self.grid.addWidget(wc.QLabel("Number of points for analysis:"), 3, 0)
        analPoints = int(np.floor(self.father.current.len() * 3 / 4.0))
        self.aPointsEntry = wc.QLineEdit(analPoints)
        self.grid.addWidget(self.aPointsEntry, 4, 0)
        self.grid.addWidget(wc.QLabel("Max number of frequencies:"), 5, 0)
        numberFreq = 20
        self.nFreqEntry = wc.QLineEdit(numberFreq)
        self.grid.addWidget(self.nFreqEntry, 6, 0)
        self.grid.addWidget(wc.QLabel("Number of points to predict:"), 7, 0)
        predictPoints = 10
        self.nPredictEntry = wc.QLineEdit(predictPoints)
        self.grid.addWidget(self.nPredictEntry, 8, 0)

    def applyFunc(self):
        analPoints = safeEval(self.aPointsEntry.text(), length=self.father.current.len(), Type='FI')
        if analPoints is None:
            raise SsnakeException('LPSVD: Number of points for analysis is not valid')
        numberFreq = safeEval(self.nFreqEntry.text(), length=self.father.current.len(), Type='FI')
        if numberFreq is None:
            raise SsnakeException('LPSVD: Number of frequencies is not valid')
        predictPoints = safeEval(self.nPredictEntry.text(), length=self.father.current.len(), Type='FI')
        if predictPoints is None:
            raise SsnakeException('LPSVD: Number of points to predict is not valid')
        if analPoints > self.father.current.len():
            raise SsnakeException('LPSVD: Number of points for analysis cannot be larger than data size')
        if analPoints < 2:
            raise SsnakeException('LPSVD: Number of points for analysis should be at least 2')
        if self.specGroup.checkedId() == 0:
            forward = True
        else:
            forward = False
        self.father.current.lpsvd(predictPoints, numberFreq, forward, analPoints)
        self.father.sideframe.upd()

###########################################################################

class ScaleSWWindow(wc.ToolWindow):

    NAME = "Scale SW"
    def __init__(self, parent):
        super(ScaleSWWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Scale Factor:"), 0, 0)
        self.scaleDropdown = QtWidgets.QComboBox()
        drop_list = ['User Defined', 
                     'Spin 3/2, 3QMAS', 'Spin 3/2, ST1MAS', 'Spin 3/2, DQ-STMAS',
                     'Spin 5/2, 3QMAS', 'Spin 5/2, 5QMAS', 'Spin 5/2, ST1MAS', 'Spin 5/2, DQ-STMAS',
                     'Spin 7/2, 3QMAS', 'Spin 5/2, 5QMAS', 'Spin 5/2, 7QMAS', 'Spin 7/2, ST1MAS', 'Spin 7/2, DQ-STMAS',
                     'Spin 9/2, 3QMAS', 'Spin 9/2, 5QMAS', 'Spin 9/2, 7QMAS', 'Spin 9/2, 9QMAS', 'Spin 9/2, ST1MAS', 'Spin 9/2, DQ-STMAS',
                    ]
        self.scaleDropdown.addItems(drop_list)
        #['User Defined','Spin 3/2, -3Q (9/34)',  'Spin 3/2, ST-1 (8/9)', 'Spin 5/2, 3Q (-12/17)', 
        #                             'Spin 5/2, -5Q (12/85)', 'Spin 5/2, ST1 (12/85)', 'Spin 7/2, 3Q (-45/34)',
        #                             'Spin 7/2, 5Q (-9/34)', 'Spin 7/2, -7Q (45/476)', 'Spin 9/2, 3Q (-36/17)', 'Spin 9/2, 5Q (-36/85)', 'Spin 9/2, 7Q (-18/117)', 'Spin 9/2, -9Q (6/85)'])
        self.scaleDropdown.activated.connect(self.dropdownChanged)
#        self.scaleList = [0, 9.0/34.0, 9/17, -12.0/17.0, 12.0/85.0, -45.0/34.0, -9/34.0, 45.0/476.0, -36.0/17.0, -36.0/85.0, -18.0/117.0, 6.0/85.0]
        self.scaleList = [1,
              (3/2, -3/2, 3/2), (3/2, -3/2, -1/2), (3/2, -3/2, 1/2),
              (5/2, -3/2, 3/2), (5/2, -5/2, 5/2), (5/2, -3/2, -1/2), (5/2, -3/2, 1/2),
              (7/2, -3/2, 3/2), (7/2, -5/2, 5/2), (7/2, -7/2, 7/2), (7/2, -3/2, -1/2), (7/2, -3/2, 1/2),
              (9/2, -3/2, 3/2), (9/2, -5/2, 5/2), (9/2, -7/2, 7/2), (9/2, -9/2, 9/2), (9/2, -3/2, -1/2), (9/2, -3/2, 1/2),
                        ]
        self.grid.addWidget(self.scaleDropdown, 1, 0)
        self.scaleEntry = wc.QLineEdit("0.0")
        self.grid.addWidget(self.scaleEntry, 3, 0)

    def dropdownChanged(self):
        index = self.scaleDropdown.currentIndex()
        if index == 0:
            scale = "1"
        else:
            scale = f"{func.scale_SW_ratio(*self.scaleList[index])}"
        self.scaleEntry.setText(scale)
        #self.scaleEntry.setText("%.9f" % self.scaleList[index])

    def applyFunc(self):
        scale = safeEval(self.scaleEntry.text(), length=self.father.current.len(), Type='FI')
        if scale is None:
            raise SsnakeException("Scale SW: Factor not a valid value")
        self.father.current.scaleSw(scale)


class ScaleFreqRefWindow(wc.ToolWindow):

    NAME = "Scale Freq and Ref for MQMAS"

    def __init__(self, parent):
        super(ScaleFreqRefWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("MQMAS Scale Factor:"), 0, 0)
        self.scaleDropdown = QtWidgets.QComboBox()
        # f_ratio = abs(ratio - mq) with mq negative if mq=S+1/2
        drop_list = ['User Defined', 
                     'Spin 3/2, 3QMAS', 'Spin 3/2, ST1MAS', 'Spin 3/2, DQ-STMAS',
                     'Spin 5/2, 3QMAS', 'Spin 5/2, 5QMAS', 'Spin 5/2, ST1MAS', 'Spin 5/2, DQ-STMAS',
                     'Spin 7/2, 3QMAS', 'Spin 5/2, 5QMAS', 'Spin 5/2, 7QMAS', 'Spin 7/2, ST1MAS', 'Spin 7/2, DQ-STMAS',
                     'Spin 9/2, 3QMAS', 'Spin 9/2, 5QMAS', 'Spin 9/2, 7QMAS', 'Spin 9/2, 9QMAS', 'Spin 9/2, ST1MAS', 'Spin 9/2, DQ-STMAS',
                    ]
        self.scaleDropdown.addItems(drop_list)
        #['User Defined', 'Spin 3/2, -3Q', 'Spin 5/2, 3Q', 'Spin 5/2, -5Q',
        #                             'Spin 7/2, 3Q', 'Spin 7/2, 5Q', 'Spin 7/2, -7Q', 
        #                             'Spin 9/2, 3Q', 'Spin 9/2, 5Q', 'Spin 9/2, 7Q', 'Spin 9/2, -9Q'])
        self.scaleDropdown.activated.connect(self.dropdownChanged)
        self.scaleList = [1,
                          (3/2, -3/2, 3/2), (3/2, -3/2, -1/2), (3/2, -3/2, 1/2),
                          (5/2, -3/2, 3/2), (5/2, -5/2, 5/2), (5/2, -3/2, -1/2), (5/2, -3/2, 1/2),
                          (7/2, -3/2, 3/2), (7/2, -5/2, 5/2), (7/2, -7/2, 7/2), (7/2, -3/2, -1/2), (7/2, -3/2, 1/2),
                          (9/2, -3/2, 3/2), (9/2, -5/2, 5/2), (9/2, -7/2, 7/2), (9/2, -9/2, 9/2), (9/2, -3/2, -1/2), (9/2, -3/2, 1/2),
                         ]
        self.grid.addWidget(self.scaleDropdown, 1, 0)
        self.scaleEntry = wc.QLineEdit("0.0")
        self.grid.addWidget(self.scaleEntry, 3, 0)

    def dropdownChanged(self):
        index = self.scaleDropdown.currentIndex()
        if index == 0:
            scale = "1"
        else:
            scale = f"{func.scale_CarRef_ratio(*self.scaleList[index])}"
        self.scaleEntry.setText(scale)

    def applyFunc(self):
        scale = safeEval(self.scaleEntry.text(), length=self.father.current.len(), Type='FI')
        if scale is None:
            raise SsnakeException("Scale: Factor not a valid value")
        freq = self.father.current.freq()
        ref = self.father.current.ref()
        sw = self.father.current.sw()
        self.father.current.setFreq(freq*scale, sw)
        self.father.current.setRef(ref*scale)


###########################################################################

class ShiftDataWindow(wc.ToolWindow):

    NAME = "Shifting data"
    SINGLESLICE = True

    def __init__(self, parent):
        super(ShiftDataWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Data points to shift:"), 0, 0, 1, 3)
        self.shiftVal = 0
        self.shiftEntry = wc.QLineEdit(self.shiftVal, self.shiftPreview)
        self.shiftEntry.setMinimumWidth(100)
        self.grid.addWidget(self.shiftEntry, 1, 1)
        leftShift = QtWidgets.QPushButton("<")
        leftShift.clicked.connect(self.stepDownShift)
        leftShift.setAutoRepeat(True)
        self.grid.addWidget(leftShift, 1, 0)
        rightShift = QtWidgets.QPushButton(">")
        rightShift.clicked.connect(self.stepUpShift)
        rightShift.setAutoRepeat(True)
        self.grid.addWidget(rightShift, 1, 2)

    def stepUpShift(self, *args):
        inp = safeEval(self.shiftEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Shift data: shift value not valid")
        self.shiftVal = int(round(inp))
        shift = 1
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier and QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            shift *= 1000
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            shift *= 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            shift *= 100
        self.shiftVal = self.shiftVal + shift
        self.shiftEntry.setText(str(self.shiftVal))
        self.shiftPreview()

    def stepDownShift(self, *args):
        inp = safeEval(self.shiftEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Shift data: shift value not valid")
        self.shiftVal = int(round(inp))
        shift = -1
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier and QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            shift *= 1000
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            shift *= 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            shift *= 100
        self.shiftVal = self.shiftVal + shift
        self.shiftEntry.setText(str(self.shiftVal))
        self.shiftPreview()

    def shiftPreview(self, *args):
        inp = safeEval(self.shiftEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Shift data: shift value not valid")
        self.shiftVal = int(round(inp))
        self.shiftEntry.setText(str(self.shiftVal))
        self.father.current.shiftPreview(self.shiftVal)

    def applyFunc(self):
        inp = safeEval(self.shiftEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Shift data: shift value not valid")
        shift = int(round(inp))
        self.father.current.shift(shift, (self.singleSlice.isChecked()))

###########################################################################

class RollDataWindow(wc.ToolWindow):

    NAME = "Roll data"
    SINGLESLICE = True

    def __init__(self, parent):
        super(RollDataWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Data points to roll:"), 0, 0, 1, 3)
        self.shiftVal = 0
        self.shiftEntry = wc.QLineEdit(self.shiftVal, self.rollPreview)
        self.shiftEntry.setMinimumWidth(100)
        self.grid.addWidget(self.shiftEntry, 1, 1)
        leftShift = QtWidgets.QPushButton("<")
        leftShift.clicked.connect(self.stepDownShift)
        leftShift.setAutoRepeat(True)
        self.grid.addWidget(leftShift, 1, 0)
        rightShift = QtWidgets.QPushButton(">")
        rightShift.clicked.connect(self.stepUpShift)
        rightShift.setAutoRepeat(True)
        self.grid.addWidget(rightShift, 1, 2)
        self.shift_axisCB = QtWidgets.QCheckBox("Shift axis")
        self.layout.addWidget(self.shift_axisCB, 2, 0)

    def stepUpShift(self, *args):
        inp = safeEval(self.shiftEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Roll data: roll value not valid")
        self.shiftVal = inp
        shift = 1
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier and QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            shift *= 1000
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            shift *= 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            shift *= 100
        self.shiftVal = self.shiftVal + shift
        self.shiftEntry.setText(str(self.shiftVal))
        self.rollPreview()

    def stepDownShift(self, *args):
        inp = safeEval(self.shiftEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Roll data: roll value not valid")
        self.shiftVal = inp
        shift = -1
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier and QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            shift *= 1000
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            shift *= 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            shift *= 100
        self.shiftVal = self.shiftVal + shift
        self.shiftEntry.setText(str(self.shiftVal))
        self.rollPreview()

    def rollPreview(self, *args):
        inp = safeEval(self.shiftEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Roll data: roll value not valid")
        self.shiftVal = inp
        self.shiftEntry.setText(str(self.shiftVal))
        self.father.current.rollPreview(self.shiftVal, shift_axis=self.shift_axisCB.isChecked())

    def applyFunc(self):
        inp = safeEval(self.shiftEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Roll data: roll value not valid")
        if self.singleSlice.isChecked() and self.shift_axisCB.isChecked():
            raise SsnakeException("Shifting axis is not allowed when single slice is checked!")
        shift = inp
        self.father.current.roll(shift, (self.singleSlice.isChecked()), shift_axis=self.shift_axisCB.isChecked())

#############################################################


class DCWindow(wc.ToolWindow):

    NAME = "Offset correction"
    SINGLESLICE = True

    def __init__(self, parent):
        super(DCWindow, self).__init__(parent)
        self.startVal = int(round(0.8 * parent.current.len()))
        self.endVal = parent.current.len()
        self.grid.addWidget(wc.QLabel("Start index:"), 0, 0)
        self.startEntry = wc.QLineEdit(self.startVal, self.offsetPreview)
        self.grid.addWidget(self.startEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End index:"), 2, 0)
        self.endEntry = wc.QLineEdit(self.endVal, self.offsetPreview)
        self.grid.addWidget(self.endEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Offset:"), 4, 0)
        val = parent.current.getdcOffset(int(round(0.8 * parent.current.len())), parent.current.len())
        self.offsetEntry = wc.QLineEdit('{:.2e}'.format(val), lambda: self.offsetPreview(True))
        self.grid.addWidget(self.offsetEntry, 5, 0)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def picked(self, pos, second=False):
        dataLength = self.father.current.len()
        if second:
            inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
            if inp is not None:
                self.startVal = int(round(inp))
            if self.startVal < 0:
                self.startVal = 0
            elif self.startVal > dataLength:
                self.startVal = dataLength
            self.endVal = pos[0]
            self.endEntry.setText(str(self.endVal))
            if inp is not None:
                self.startEntry.setText(str(self.startVal))
                val = self.father.current.getdcOffset(self.startVal, self.endVal)
                self.offsetEntry.setText('{:.2e}'.format(val))
                self.father.current.dcOffset(val)
            else:
                self.offsetEntry.setText('')
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
            self.father.current.peakPick = True
        else:
            self.startEntry.setText(str(pos[0]))
            inp = safeEval(self.endEntry.text(), length=self.father.current.len(), Type='FI')
            if inp is not None:
                self.endVal = int(round(inp))
            if self.endVal < 0:
                self.endVal = 0
            elif self.endVal > dataLength:
                self.endVal = dataLength
            self.startVal = pos[0]
            if inp is not None:
                val = self.father.current.getdcOffset(self.startVal, self.endVal)
                self.offsetEntry.setText('{:.2e}'.format(val))
            else:
                self.offsetEntry.setText('')
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, True)
            self.father.current.peakPick = True

    def offsetPreview(self, inserted=False):
        if inserted:
            dcVal = safeEval(self.offsetEntry.text(), length=self.father.current.len(), Type='C')
            if dcVal is None:
                raise SsnakeException("Offset correction: offset value not valid")
            self.father.current.dcOffset(dcVal)
        else:
            dataLength = self.father.current.len()
            inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
            if inp is None:
                raise SsnakeException("Offset correction: start value not valid")
            self.startVal = int(round(inp))
            if self.startVal < 0:
                self.startVal = 0
            elif self.startVal > dataLength:
                self.startVal = dataLength
            self.startEntry.setText(str(self.startVal))
            inp = safeEval(self.endEntry.text(), length=self.father.current.len(), Type='FI')
            if inp is None:
                raise SsnakeException("Offset correction: end value not valid")
            self.endVal = int(round(inp))
            if self.endVal < 0:
                self.endVal = 0
            elif self.endVal > dataLength:
                self.endVal = dataLength
            self.endEntry.setText(str(self.endVal))
            val = self.father.current.getdcOffset(self.startVal, self.endVal)
            self.offsetEntry.setText('{:.2e}'.format(val))
            self.father.current.dcOffset(val)

    def applyFunc(self):
        inp = safeEval(self.offsetEntry.text(), length=self.father.current.len(), Type='C')
        if inp is None:
            raise SsnakeException("Offset correction: offset value not valid")
        self.father.current.peakPickReset()
        self.father.current.subtract([inp], self.singleSlice.isChecked())

#############################################################


class BaselineWindow(wc.ToolWindow):

    NAME = "Baseline correction"
    SINGLESLICE = True
    TYPES = ['poly','sin/cos']
    TYP_NAMES = ['Polynomial','sine/cosine']

    def __init__(self, parent):
        super(BaselineWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Type:"), 0, 0, 1, 2)
        self.typeDropdown = QtWidgets.QComboBox()
        self.typeDropdown.addItems(self.TYP_NAMES)
        self.grid.addWidget(self.typeDropdown, 1, 0, 1, 2)
        self.grid.addWidget(wc.QLabel("Degree:"), 2, 0, 1, 2)
        self.removeList = []
        self.degreeEntry = wc.SsnakeSpinBox()
        self.degreeEntry.setMaximum(10000)
        self.degreeEntry.setMinimum(1)
        self.degreeEntry.setValue(3)
        self.degreeEntry.setAlignment(QtCore.Qt.AlignCenter)
        self.grid.addWidget(self.degreeEntry, 3, 0, 1, 2)
        self.invertButton = QtWidgets.QCheckBox("Invert selection")
        self.invertButton.stateChanged.connect(self.preview)
        self.grid.addWidget(self.invertButton, 4, 0, 1, 2)
        self.allFitButton = QtWidgets.QCheckBox("Fit traces separately")
        self.grid.addWidget(self.allFitButton, 5, 0, 1, 2)
        resetButton = QtWidgets.QPushButton("&Reset")
        resetButton.clicked.connect(self.reset)
        self.grid.addWidget(resetButton, 6, 0)
        fitButton = QtWidgets.QPushButton("&Fit")
        fitButton.clicked.connect(self.preview)
        self.grid.addWidget(fitButton, 6, 1)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def picked(self, pos):
        self.removeList.append(pos[0])
        self.father.current.previewRemoveList(self.removeList, invert=self.invertButton.isChecked())
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def preview(self, *args):
        type = self.TYPES[self.typeDropdown.currentIndex()]
        inp = self.degreeEntry.value()
        self.father.current.previewRemoveList(self.removeList, invert=self.invertButton.isChecked())
        self.father.current.previewBaselineCorrection(inp, self.removeList, type, invert=self.invertButton.isChecked())
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def reset(self, *args):
        self.removeList = []
        self.father.current.resetPreviewRemoveList()
        self.preview()

    def closeEvent(self, *args):
        self.father.current.removeListLines = []
        del self.father.current.removeListLines
        super(BaselineWindow, self).closeEvent(*args)

    def applyFunc(self):
        inp = self.degreeEntry.value()
        type = self.TYPES[self.typeDropdown.currentIndex()]
        if self.allFitButton.isChecked():
            self.father.current.baselineCorrectionAll(inp, self.removeList, type, invert=self.invertButton.isChecked())
        else:
            self.father.current.baselineCorrection(inp, self.removeList, type, self.singleSlice.isChecked(), invert=self.invertButton.isChecked())
        self.father.current.peakPickReset()
        self.father.current.resetPreviewRemoveList()

#############################################################


class regionWindow(wc.ToolWindow):

    def __init__(self, parent, name):
        self.NAME = name
        super(regionWindow, self).__init__(parent)
        self.startVal = [0]  # dummy variables
        self.endVal = [parent.current.len()]  # dummy variables
        self.grid.addWidget(wc.QLabel("Start index:"), 0, 0)
        self.grid.addWidget(wc.QLabel("End index:"), 0, 1)
        self.startEntry = []
        self.endEntry = []
        self.deleteButton = []
        self.partIter = 0
        self.entryCount = 1
        self.first = True
        self.startEntry.append(wc.QLineEdit(""))
        self.startEntry[0].editingFinished.connect(lambda self=self, tmp=self.startEntry[0]: self.setVal(tmp, True))
        self.grid.addWidget(self.startEntry[0], 1, 0)
        self.endEntry.append(wc.QLineEdit(""))
        self.endEntry[0].editingFinished.connect(lambda self=self, tmp=self.endEntry[0]: self.setVal(tmp, False))
        self.grid.addWidget(self.endEntry[0], 1, 1)
        self.deleteButton.append(QtWidgets.QPushButton("X"))
        self.deleteButton[0].clicked.connect(lambda extra, self=self: self.deleteEntry(self.deleteButton[0]))
        self.grid.addWidget(self.deleteButton[0], 1, 2)
        self.newSpec = QtWidgets.QCheckBox("Result in new workspace")
        self.layout.addWidget(self.newSpec, 1, 0, 1, 2)
        self.grid.setRowStretch(100, 1)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def addValue(self, pos):
        if self.first:
            self.startVal[self.partIter] = pos
            self.startEntry[self.partIter].setText(str(pos))
            self.first = False
        else:
            tmp = self.startVal[self.partIter]
            self.startVal[self.partIter] = min(pos, tmp)
            self.endVal[self.partIter] = max(pos, tmp)
            self.startVal = np.append(self.startVal, 0)
            self.endVal = np.append(self.endVal, self.father.current.len())
            self.startEntry[self.partIter].setText(str(self.startVal[self.partIter]))
            self.endEntry[self.partIter].setText(str(self.endVal[self.partIter]))
            self.partIter += 1
            self.startEntry.append(wc.QLineEdit())
            self.startEntry[self.partIter].editingFinished.connect(lambda self=self, tmp=self.startEntry[self.partIter]: self.setVal(tmp, True))
            self.grid.addWidget(self.startEntry[self.partIter], 1 + self.entryCount, 0)
            self.endEntry.append(wc.QLineEdit())
            self.endEntry[self.partIter].editingFinished.connect(lambda self=self, tmp=self.endEntry[self.partIter]: self.setVal(tmp, False))
            self.grid.addWidget(self.endEntry[self.partIter], 1 + self.entryCount, 1)
            self.deleteButton.append(QtWidgets.QPushButton("X"))
            self.deleteButton[self.partIter].clicked.connect(lambda extra, self=self, tmp=self.deleteButton[self.partIter]: self.deleteEntry(tmp))
            self.grid.addWidget(self.deleteButton[self.partIter], 1 + self.entryCount, 2)
            self.entryCount += 1
            self.first = True

    def deleteEntry(self, button=None, num=None):
        if num is None:
            num = self.deleteButton.index(button)
        if num == self.partIter:
            self.startVal[num] = 0
            self.endVal[num] = self.father.current.len()
            self.startEntry[num].setText("")
            self.endEntry[num].setText("")
            self.first = True
            return
        self.grid.removeWidget(self.endEntry[num])
        self.grid.removeWidget(self.startEntry[num])
        self.grid.removeWidget(self.deleteButton[num])
        self.endEntry[num].deleteLater()
        self.startEntry[num].deleteLater()
        self.deleteButton[num].deleteLater()
        self.endEntry.pop(num)
        self.startEntry.pop(num)
        self.deleteButton.pop(num)
        self.startVal = np.delete(self.startVal, num)
        self.endVal = np.delete(self.endVal, num)
        self.partIter -= 1

    def picked(self, pos):
        self.addValue(pos[0])
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def setVal(self, entry, isMin=False):
        inp = safeEval(entry.text(), length=self.father.current.len(), Type='FI')
        error = False
        if inp is not None:
            inp = int(inp)
            if inp < 0:
                inp = 0
            if inp > self.father.current.len():
                inp = self.father.current.len()
        if isMin:
            num = self.startEntry.index(entry)
            if inp is None:
                self.startVal[num] = -1  # If the input is wrong, use -1 as a placeholder for it in the value list
                self.father.father.dispMsg(self.NAME + ": wrong input")
                error = True
            elif self.endVal[num] == -1:
                self.startVal[num] = inp
            else:
                self.startVal[num] = min(inp, self.endVal[num])
                self.endVal[num] = max(inp, self.endVal[num])
        else:
            num = self.endEntry.index(entry)
            if inp is None:
                self.endVal[num] = -1
                self.father.father.dispMsg(self.NAME + ": wrong input")
                error = True
            elif self.startVal[num] == -1:
                self.endVal[num] = inp
            else:
                self.endVal[num] = max(inp, self.startVal[num])
                self.startVal[num] = min(inp, self.startVal[num])
        if num == self.partIter:
            self.partIter += 1
            self.startVal = np.append(self.startVal, 0)
            self.endVal = np.append(self.endVal, self.father.current.len())
            self.startEntry.append(wc.QLineEdit())
            self.startEntry[self.partIter].editingFinished.connect(lambda self=self, tmp=self.startEntry[self.partIter]: self.setVal(tmp, True))
            self.grid.addWidget(self.startEntry[self.partIter], 1 + self.entryCount, 0)
            self.endEntry.append(wc.QLineEdit())
            self.endEntry[self.partIter].editingFinished.connect(lambda self=self, tmp=self.endEntry[self.partIter]: self.setVal(tmp, False))
            self.grid.addWidget(self.endEntry[self.partIter], 1 + self.entryCount, 1)
            self.deleteButton.append(QtWidgets.QPushButton("X"))
            self.deleteButton[self.partIter].clicked.connect(lambda extra, self=self, tmp=self.deleteButton[self.partIter]: self.deleteEntry(tmp))
            self.grid.addWidget(self.deleteButton[self.partIter], 1 + self.entryCount, 2)
            self.entryCount += 1
            self.first = True
        if error:  # Return only after partIter check
            return
        if self.startVal[num] != -1:  # Only if the input is OK, reprint it
            self.startEntry[num].setText(str(self.startVal[num]))
        if self.endVal[num] != -1:
            self.endEntry[num].setText(str(self.endVal[num]))

    def apply(self, maximum, minimum, newSpec):
        pass

    def applyFunc(self):
        if self.partIter == 0:
            self.apply(np.array([0]), np.array([self.father.current.len()]), self.newSpec.isChecked())
        else:
            self.apply(self.startVal[:self.partIter], self.endVal[:self.partIter], self.newSpec.isChecked())

############################################################


class integrateWindow(regionWindow):

    def __init__(self, parent):
        super(integrateWindow, self).__init__(parent, 'Integrate')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):  # Check for errors in the inputs
            raise SsnakeException(self.NAME + ": wrong input")
        if newSpec:
            self.father.father.newWorkspace(self.father.current.integrate(minimum, maximum, newSpec))
        else:
            self.father.current.integrate(minimum, maximum, newSpec)
            self.father.updAllFrames()

############################################################


class sumWindow(regionWindow):

    def __init__(self, parent):
        super(sumWindow, self).__init__(parent, 'Sum')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            raise SsnakeException(self.NAME + ": wrong input")
        if newSpec:
            self.father.father.newWorkspace(self.father.current.sum(minimum, maximum, newSpec))
        else:
            self.father.current.sum(minimum, maximum, newSpec)
            self.father.updAllFrames()

############################################################


class maxWindow(regionWindow):

    def __init__(self, parent):
        super(maxWindow, self).__init__(parent, 'Max')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            raise SsnakeException(self.NAME + ": wrong input")
        if newSpec:
            self.father.father.newWorkspace(self.father.current.max(minimum, maximum, newSpec))
        else:
            self.father.current.max(minimum, maximum, newSpec)
            self.father.updAllFrames()

############################################################


class minWindow(regionWindow):

    def __init__(self, parent):
        super(minWindow, self).__init__(parent, 'Min')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            raise SsnakeException(self.NAME + ": wrong input")
        if newSpec:
            self.father.father.newWorkspace(self.father.current.min(minimum, maximum, newSpec))
        else:
            self.father.current.min(minimum, maximum, newSpec)
            self.father.updAllFrames()

############################################################


class argmaxWindow(regionWindow):

    def __init__(self, parent):
        super(argmaxWindow, self).__init__(parent, 'Max position')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            raise SsnakeException(self.NAME + ": wrong input")
        if newSpec:
            self.father.father.newWorkspace(self.father.current.argmax(minimum, maximum, newSpec))
        else:
            self.father.current.argmax(minimum, maximum, newSpec)
            self.father.updAllFrames()

############################################################


class argminWindow(regionWindow):

    def __init__(self, parent):
        super(argminWindow, self).__init__(parent, 'Min position')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            raise SsnakeException(self.NAME + ": wrong input")
        if newSpec:
            self.father.father.newWorkspace(self.father.current.argmin(minimum, maximum, newSpec))
        else:
            self.father.current.argmin(minimum, maximum, newSpec)
            self.father.updAllFrames()

############################################################


class avgWindow(regionWindow):

    def __init__(self, parent):
        super(avgWindow, self).__init__(parent, 'Average')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            raise SsnakeException(self.NAME + ": wrong input")
        if newSpec:
            self.father.father.newWorkspace(self.father.current.average(minimum, maximum, newSpec))
        else:
            self.father.current.average(minimum, maximum, newSpec)
            self.father.updAllFrames()

#############################################################


class regionWindow2(wc.ToolWindow):

    def __init__(self, parent, name, newSpecOption):
        self.NAME = name
        super(regionWindow2, self).__init__(parent)
        self.startVal = 0
        self.endVal = parent.current.len()
        self.grid.addWidget(wc.QLabel("Start index:"), 0, 0)
        self.startEntry = wc.QLineEdit(self.startVal, self.checkValues)
        self.grid.addWidget(self.startEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End index:"), 2, 0)
        self.endEntry = wc.QLineEdit(self.endVal, self.checkValues)
        self.grid.addWidget(self.endEntry, 3, 0)
        self.newSpec = QtWidgets.QCheckBox("Result in new workspace")
        if not newSpecOption:
            self.newSpec.hide()
        self.layout.addWidget(self.newSpec, 1, 0, 1, 2)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def preview(self, maximum, minimum):
        pass

    def picked(self, pos, second=False):
        if second:
            dataLength = self.father.current.len()
            inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
            if inp is not None:
                self.startVal = int(round(inp))
            if self.startVal < 0:
                self.startVal = 0
            elif self.startVal > dataLength:
                self.startVal = dataLength
            self.startEntry.setText(str(self.startVal))
            self.endVal = pos[0]
            self.endEntry.setText(str(self.endVal))
            self.preview(self.startVal, self.endVal)
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
            self.father.current.peakPick = True
        else:
            self.startEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, True)
            self.father.current.peakPick = True

    def checkValues(self, *args):
        dataLength = self.father.current.len()
        inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is not None:
            self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        self.startEntry.setText(str(self.startVal))
        inp = safeEval(self.endEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is not None:
            self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        self.endEntry.setText(str(self.endVal))
        self.preview(self.startVal, self.endVal)

    def applyFunc(self):
        dataLength = self.father.current.len()
        inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException(self.NAME + ": value not valid")
        inp = safeEval(self.endEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException(self.NAME + ": value not valid")
        self.checkValues(self)
        self.apply(self.startVal, self.endVal, self.newSpec.isChecked())

    def apply(self, maximum, minimum, newSpec):
        pass

#############################################################


class regionWindowStep(wc.ToolWindow):

    def __init__(self, parent, name, newSpecOption):
        self.NAME = name
        super(regionWindowStep, self).__init__(parent)
        self.startVal = 0
        self.endVal = parent.current.len()
        self.stepVal = 1;
        self.grid.addWidget(wc.QLabel("Start index:"), 0, 0)
        self.startEntry = wc.QLineEdit(self.startVal, self.checkValues)
        self.grid.addWidget(self.startEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End index:"), 2, 0)
        self.endEntry = wc.QLineEdit(self.endVal, self.checkValues)
        self.grid.addWidget(self.endEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Step size:"), 4, 0)
        self.stepEntry = wc.QLineEdit(self.stepVal, self.checkValues)
        self.grid.addWidget(self.stepEntry, 5, 0)
        self.newSpec = QtWidgets.QCheckBox("Result in new workspace")
        if not newSpecOption:
            self.newSpec.hide()
        self.layout.addWidget(self.newSpec, 1, 0, 1, 2)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def preview(self, maximum, minimum, step):
        pass

    def picked(self, pos, second=False):
        if second:
            dataLength = self.father.current.len()
            inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
            if inp is not None:
                self.startVal = int(round(inp))
            if self.startVal < 0:
                self.startVal = 0
            elif self.startVal > dataLength:
                self.startVal = dataLength
            self.startEntry.setText(str(self.startVal))
            self.endVal = pos[0]
            self.endEntry.setText(str(self.endVal))
            self.preview(self.startVal, self.endVal, self.stepVal)
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
            self.father.current.peakPick = True
        else:
            self.startEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, True)
            self.father.current.peakPick = True

    def checkValues(self, *args):
        dataLength = self.father.current.len()
        inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is not None:
            self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        self.startEntry.setText(str(self.startVal))
        inp = safeEval(self.endEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is not None:
            self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        self.endEntry.setText(str(self.endVal))
        inp = safeEval(self.stepEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is not None:
            self.stepVal = int(round(inp))
        if self.stepVal <= 0:
            self.stepVal = 1
        elif self.stepVal > dataLength:
            self.stepVal = dataLength
        self.stepEntry.setText(str(self.stepVal))
        self.preview(self.startVal, self.endVal, self.stepVal)

    def applyFunc(self):
        dataLength = self.father.current.len()
        inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException(self.NAME + ": value not valid")
        inp = safeEval(self.endEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException(self.NAME + ": value not valid")
        inp = safeEval(self.stepEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException(self.NAME + ": value not valid")
        self.checkValues(self)
        self.apply(self.startVal, self.endVal, self.newSpec.isChecked(), self.stepVal)

    def apply(self, maximum, minimum, newSpec, step):
        pass
    
############################################################


class extractRegionWindow(regionWindowStep):

    def __init__(self, parent):
        super(extractRegionWindow, self).__init__(parent, 'Extract part', True)

    def apply(self, maximum, minimum, newSpec, step):
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.extract(minimum, maximum, newSpec, step)) is None:
                return None
        else:
            self.father.current.extract(minimum, maximum, newSpec, step)
            self.father.updAllFrames()
        return 1

############################################################


class SubtractAvgWindow(regionWindow2):

    def __init__(self, parent):
        super(SubtractAvgWindow, self).__init__(parent, 'Subtract Avg', False)

    def apply(self, maximum, minimum, newSpec):
        self.father.current.subtractAvg(maximum, minimum)
        self.father.updAllFrames()
        return 1

    def preview(self, maximum, minimum):
        self.father.current.subtractAvgPreview(maximum, minimum)

############################################################


class AlignDataWindow(regionWindow2):

    def __init__(self, parent):
        super(AlignDataWindow, self).__init__(parent, 'Align Maxima', False)

    def apply(self, maximum, minimum, newSpec):
        self.father.current.align(maximum, minimum)
        self.father.updAllFrames()
        return 1

    # def preview(self, maximum, minimum):
    #     self.father.current.subtractAvgPreview(maximum, minimum)

#############################################################


class FiddleWindow(wc.ToolWindow):

    NAME = "Reference deconvolution"

    def __init__(self, parent):
        super(FiddleWindow, self).__init__(parent)
        self.startVal = 0
        self.endVal = parent.current.len()
        self.grid.addWidget(wc.QLabel("Start index:"), 0, 0)
        self.startEntry = wc.QLineEdit(self.startVal, self.checkValues)
        self.grid.addWidget(self.startEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End index:"), 2, 0)
        self.endEntry = wc.QLineEdit(self.endVal, self.checkValues)
        self.grid.addWidget(self.endEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Linebroadening [Hz]:"), 4, 0)
        self.lbEntry = wc.QLineEdit("1.0", self.checkValues)
        self.grid.addWidget(self.lbEntry, 5, 0)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def picked(self, pos, second=False):
        if second:
            dataLength = self.father.current.len()
            inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
            if inp is not None:
                self.startVal = int(round(inp))
            if self.startVal < 0:
                self.startVal = 0
            elif self.startVal > dataLength:
                self.startVal = dataLength
            self.startEntry.setText(str(self.startVal))
            self.endVal = pos[0]
            self.endEntry.setText(str(self.endVal))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
            self.father.current.peakPick = True
        else:
            self.startEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, True)
            self.father.current.peakPick = True

    def checkValues(self, *args):
        dataLength = self.father.current.len()
        inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is not None:
            self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        self.startEntry.setText(str(self.startVal))
        inp = safeEval(self.endEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is not None:
            self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        self.endEntry.setText(str(self.endVal))
        inp = safeEval(self.lbEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is not None:
            self.lbEntry.setText(str(inp))

    def applyFunc(self):
        dataLength = self.father.current.len()
        inp = safeEval(self.startEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Reference deconv: start entry not valid")
        self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        inp = safeEval(self.endEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Reference deconv: end entry not valid")
        self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        lb = safeEval(self.lbEntry.text(), length=self.father.current.len(), Type='FI')
        if lb is None:
            raise SsnakeException("Reference deconv: Linebroadening entry not valid")
        self.father.current.fiddle(self.startVal, self.endVal, lb)

##############################################################


class DeleteWindow(wc.ToolWindow):

    NAME = "Delete"

    def __init__(self, parent):
        super(DeleteWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Indexes to delete:"), 0, 0)
        self.delEntry = wc.QLineEdit('0', self.preview)
        self.grid.addWidget(self.delEntry, 1, 0)

    def preview(self, *args):
        length = int(self.father.current.len())
        pos = safeEval(self.delEntry.text(), length=self.father.current.len())
        if pos is None:
            raise SsnakeException('Delete: not all values are valid indexes to delete')
        pos = np.array(pos)
        pos[pos < 0] = pos[pos < 0] + length
        if (pos < 0).any() or (pos >= length).any():
            raise SsnakeException('Delete: not all values are valid indexes to delete')
        self.father.current.deletePreview(pos)

    def applyFunc(self):
        length = self.father.current.len()
        pos = safeEval(self.delEntry.text(), length=self.father.current.len())
        if pos is None:
            raise SsnakeException('Delete: not all values are valid indexes to delete')
        if isinstance(pos, (int, float)):
            pos = np.array([pos])
        else:
            pos = np.array(pos)
        pos[pos < 0] = pos[pos < 0] + length
        if (pos < 0).any() or (pos >= length).any():
            raise SsnakeException('Delete: not all values are valid indexes to delete')
        self.father.current.delete(pos)

##############################################################


class SplitWindow(wc.ToolWindow):

    NAME = "Split"

    def __init__(self, parent):
        super(SplitWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Sections:"), 0, 0)
        self.splitEntry = wc.QLineEdit('1', self.preview)
        self.grid.addWidget(self.splitEntry, 1, 0)

    def preview(self, *args):
        val = safeEval(self.splitEntry.text(), length=self.father.current.len(), Type='FI')
        if val is None:
            raise SsnakeException("Split: input not valid")
        self.splitEntry.setText(str(int(round(val))))

    def applyFunc(self):
        val = safeEval(self.splitEntry.text(), length=self.father.current.len(), Type='FI')
        if val is None:
            raise SsnakeException("Split: input not valid")
        val = int(val)
        if val <= 0:
            raise SsnakeException("Split: input not valid")
        self.father.current.split(int(round(val)))

##############################################################


class ConcatenateWindow(wc.ToolWindow):

    NAME = "Concatenate"

    def __init__(self, parent):
        super(ConcatenateWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Concatenation axes:"), 0, 0)
        self.axesEntry = QtWidgets.QComboBox()
        self.axesEntry.addItems(np.array(np.arange(self.father.current.data.ndim() - 1) + 1, dtype=str))
        self.grid.addWidget(self.axesEntry, 1, 0)

    def applyFunc(self):
        self.father.current.concatenate(self.axesEntry.currentIndex())

##############################################################


class InsertWindow(wc.ToolWindow):

    NAME = "Insert"

    def __init__(self, parent):
        super(InsertWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Start insert at index:"), 0, 0)
        self.posEntry = wc.QLineEdit(self.father.current.len(), self.preview)
        self.grid.addWidget(self.posEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("Workspace to insert:"), 2, 0)
        self.wsEntry = QtWidgets.QComboBox()
        self.wsEntry.addItems(self.father.father.workspaceNames)
        self.grid.addWidget(self.wsEntry, 3, 0)

    def preview(self, *args):
        pos = safeEval(self.posEntry.text(), length=self.father.current.len(), Type='FI')
        if pos is None:
            return
        pos = int(round(pos))
        if pos > self.father.current.len():
            pos = self.father.current.len()
        elif pos < 0:
            pos = 0
        self.posEntry.setText(str(pos))

    def applyFunc(self):
        pos = safeEval(self.posEntry.text(), length=self.father.current.len(), Type='FI')
        if pos is None:
            raise SsnakeException("Not a valid value")
        pos = int(round(pos))
        if pos > self.father.current.len():
            pos = self.father.current.len()
        elif pos < 0:
            pos = 0
        ws = self.wsEntry.currentIndex()
        self.father.current.insert(self.father.father.workspaces[ws].masterData.getData(), pos)

##############################################################


class CombineWindow(wc.ToolWindow):

    SINGLESLICE = True
    RESIZABLE = True

    def __init__(self, parent, combType):
        super(CombineWindow, self).__init__(parent)
        self.combType = combType  # 0 = add, 1 = subtract, 2 = multiply, 3 = divide
        if self.combType == 0:
            self.WindowTitle = "Add"
            self.grid.addWidget(wc.QLabel("Workspace to add:"), 0, 0)
        elif self.combType == 1:
            self.WindowTitle = "Subtract"
            self.grid.addWidget(wc.QLabel("Workspace to subtract:"), 0, 0)
        elif self.combType == 2:
            self.WindowTitle = "Multiply"
            self.grid.addWidget(wc.QLabel("Workspace to multiply:"), 0, 0)
        elif self.combType == 3:
            self.WindowTitle = "Divide"
            self.grid.addWidget(wc.QLabel("Workspace to divide:"), 0, 0)
        self.setWindowTitle(self.WindowTitle)
        self.wsEntry = QtWidgets.QComboBox()
        self.wsEntry.addItems(self.father.father.workspaceNames)
        self.grid.addWidget(self.wsEntry, 1, 0)

    def applyFunc(self):
        ws = self.wsEntry.currentIndex()
        if self.combType == 0:
            returnValue = self.father.current.add(self.father.father.workspaces[ws].masterData.getData(), self.singleSlice.isChecked())
        elif self.combType == 1:
            returnValue = self.father.current.subtract(self.father.father.workspaces[ws].masterData.getData(), self.singleSlice.isChecked())
        elif self.combType == 2:
            returnValue = self.father.current.multiply(self.father.father.workspaces[ws].masterData.getData(), self.singleSlice.isChecked())
        elif self.combType == 3:
            returnValue = self.father.current.divide(self.father.father.workspaces[ws].masterData.getData(), self.singleSlice.isChecked())

##############################################################


class SNWindow(wc.ToolWindow):

    NAME = "Signal to noise"
    CANCELNAME = "&Close"
    OKNAME = "C&alc"
    APPLYANDCLOSE = False

    def __init__(self, parent):
        super(SNWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Start index noise:"), 0, 0)
        self.minNoiseEntry = wc.QLineEdit('0', self.checkValues)
        self.grid.addWidget(self.minNoiseEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End index noise:"), 2, 0)
        self.maxNoiseEntry = wc.QLineEdit(parent.current.len(), self.checkValues)
        self.grid.addWidget(self.maxNoiseEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Start index signal:"), 4, 0)
        self.minEntry = wc.QLineEdit('0', self.checkValues)
        self.grid.addWidget(self.minEntry, 5, 0)
        self.grid.addWidget(wc.QLabel("End index signal:"), 6, 0)
        self.maxEntry = wc.QLineEdit(parent.current.len(), self.checkValues)
        self.grid.addWidget(self.maxEntry, 7, 0)
        self.grid.addWidget(wc.QLabel("S/N:"), 8, 0)
        self.snEntry = wc.QLineEdit("0.0")
        self.grid.addWidget(self.snEntry, 9, 0)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def picked(self, pos, num=0):
        if num == 0:
            self.minNoiseEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 1)
            self.father.current.peakPick = True
        elif num == 1:
            self.maxNoiseEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 2)
            self.father.current.peakPick = True
        elif num == 2:
            self.minEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 3)
            self.father.current.peakPick = True
        elif num == 3:
            self.maxEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 0)
            self.father.current.peakPick = True
            self.applyFunc()

    def checkValues(self, *args):
        dataLength = self.father.current.len()
        inp = safeEval(self.minNoiseEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minNoiseEntry.setText(str(minimum))
        inp = safeEval(self.maxNoiseEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxNoiseEntry.setText(str(maximum))
        inp = safeEval(self.minEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.applyFunc()

    def applyFunc(self):
        dataLength = self.father.current.len()
        inp = safeEval(self.minNoiseEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("S/N: invalid range")
        minimumNoise = int(round(inp))
        if minimumNoise < 0:
            minimumNoise = 0
        elif minimumNoise > dataLength:
            minimumNoise = dataLength
        self.minNoiseEntry.setText(str(minimumNoise))
        inp = safeEval(self.maxNoiseEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("S/N: invalid range")
        maximumNoise = int(round(inp))
        if maximumNoise < 0:
            maximumNoise = 0
        elif maximumNoise > dataLength:
            maximumNoise = dataLength
        self.maxNoiseEntry.setText(str(maximumNoise))
        inp = safeEval(self.minEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("S/N: invalid range")
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("S/N: invalid range")
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.snEntry.setText(str(self.father.current.SN(minimumNoise, maximumNoise, minimum, maximum)))

##############################################################


class FWHMWindow(wc.ToolWindow):

    NAME = "FWHM"
    CANCELNAME = "&Close"
    OKNAME = "C&alc"
    APPLYANDCLOSE = False

    def __init__(self, parent):
        super(FWHMWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Start index:"), 0, 0)
        self.minEntry = wc.QLineEdit('0', self.checkValues)
        self.grid.addWidget(self.minEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End index:"), 2, 0)
        self.maxEntry = wc.QLineEdit(parent.current.len(), self.checkValues)
        self.grid.addWidget(self.maxEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Units:"), 4, 0)
        unitSelect = self.father.current.getAxType()
        if self.father.current.spec() == 1:
            unitList = ['Hz', 'kHz', 'MHz', 'ppm']
            if self.father.current.getppm():
                unitSelect = 3
        else:
            unitList = ['s', 'ms', u'μs']
        self.unitDrop = QtWidgets.QComboBox()
        self.unitDrop.addItems(unitList)
        self.unitDrop.setCurrentIndex(unitSelect)
        self.unitDrop.currentIndexChanged.connect(self.checkValues)
        self.grid.addWidget(self.unitDrop, 5, 0)
        self.grid.addWidget(wc.QLabel(u"FWHM:"), 6, 0)
        self.fwhmEntry = wc.QLineEdit("0.0")
        self.grid.addWidget(self.fwhmEntry, 7, 0)
        self.grid.addWidget(wc.QLabel(u"0.55%:"), 8, 0)
        self.zffEntry = wc.QLineEdit("0.0")
        self.grid.addWidget(self.zffEntry, 9, 0)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def picked(self, pos, num=0):
        if num == 0:
            self.minEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 1)
            self.father.current.peakPick = True
        elif num == 1:
            self.maxEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 0)
            self.father.current.peakPick = True
            self.applyFunc()

    def checkValues(self, *args):
        dataLength = self.father.current.len()
        inp = safeEval(self.minEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.applyFunc()

    def applyFunc(self):
        dataLength = self.father.current.len()
        inp = safeEval(self.minEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("FWHM: invalid range")
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("FWHM: invalid range")
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.fwhmEntry.setText(str(self.father.current.fwhm(minimum, maximum, 0.5, self.unitDrop.currentIndex())))
        self.zffEntry.setText(str(self.father.current.fwhm(minimum, maximum, 0.0055, self.unitDrop.currentIndex())))

##############################################################


class COMWindow(wc.ToolWindow):  # Centre of Mass Window

    NAME = "Centre of Mass"
    CANCELNAME = "&Close"
    OKNAME = "C&alc"
    APPLYANDCLOSE = False

    def __init__(self, parent):
        super(COMWindow, self).__init__(parent)
        self.pickDim = 1
        if isinstance(self.father.current, views.CurrentContour):
            self.pickDim = 2
        self.grid.addWidget(wc.QLabel("X axis:"), 0, 0, 1, 2)
        self.grid.addWidget(wc.QLabel("Start:"), 1, 0)
        self.grid.addWidget(wc.QLabel("End:"), 2, 0)
        if self.pickDim == 2:
            unitSelectY = self.father.current.getAxType(-2)
            if self.father.current.spec(-2) == 1:
                unitListY = ['Hz', 'kHz', 'MHz', 'ppm']
                if self.father.current.getppm(-2):
                    unitSelectY = 3
            else:
                unitListY = ['s', 'ms', u'μs']
            self.grid.addWidget(wc.QLabel("Y axis:"), 3, 0, 1, 2)
            self.grid.addWidget(wc.QLabel("Start:"), 4, 0)
            self.grid.addWidget(wc.QLabel("End:"), 5, 0)
            self.minEntryY = wc.QLineEdit("0", lambda: self.applyFunc(False))
            self.grid.addWidget(self.minEntryY, 4, 1)
            self.maxEntryY = wc.QLineEdit(parent.current.len(-2), lambda: self.applyFunc(False))
            self.grid.addWidget(self.maxEntryY, 5, 1)
            self.grid.addWidget(wc.QLabel("Y Unit:"), 11, 0)
            self.unitDropY = QtWidgets.QComboBox()
            self.unitDropY.addItems(unitListY)
            self.unitDropY.setCurrentIndex(unitSelectY)
            self.unitDropY.currentIndexChanged.connect(lambda: self.applyFunc(True))
            self.grid.addWidget(self.unitDropY, 11, 1)
            self.comEntryY = wc.QLineEdit("0.0")
            self.grid.addWidget(wc.QLabel(u"Y COM:"), 12, 0)
            self.grid.addWidget(self.comEntryY, 12, 1)
        if self.father.current.spec() == 1:
            unitList = ['Hz', 'kHz', 'MHz', 'ppm']
            if self.father.current.getppm():
                unitSelect = 3
        else:
            unitList = ['s', 'ms', u'μs']
        self.minEntry = wc.QLineEdit("0", lambda: self.applyFunc(False))
        self.grid.addWidget(self.minEntry, 1, 1)
        self.maxEntry = wc.QLineEdit(parent.current.len(), lambda: self.applyFunc(False))
        self.grid.addWidget(self.maxEntry, 2, 1)
        unitSelect = self.father.current.getAxType()
        self.grid.addWidget(wc.QLabel(u"Centre of Mass:"), 8, 0, 1, 2)
        self.grid.addWidget(wc.QLabel("X Unit:"), 9, 0)
        self.unitDrop = QtWidgets.QComboBox()
        self.unitDrop.addItems(unitList)
        self.unitDrop.setCurrentIndex(unitSelect)
        self.unitDrop.currentIndexChanged.connect(lambda: self.applyFunc(True))
        self.grid.addWidget(self.unitDrop, 9, 1)
        self.comEntry = wc.QLineEdit("0.0")
        self.grid.addWidget(wc.QLabel(u"X COM:"), 10, 0)
        self.grid.addWidget(self.comEntry, 10, 1)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = self.pickDim

    def picked(self, pos, num=0):
        if num == 0:
            self.minEntry.setText(str(pos[0]))
            if self.pickDim == 2:
                self.minEntryY.setText(str(pos[3]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 1)
            self.father.current.peakPick = self.pickDim
        elif num == 1:
            self.maxEntry.setText(str(pos[0]))
            if self.pickDim == 2:
                self.maxEntryY.setText(str(pos[3]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 0)
            self.father.current.peakPick = self.pickDim
            self.applyFunc()

    def applyFunc(self, calc=True):
        dataLength = self.father.current.len()
        inp = safeEval(self.minEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Centre of Mass: invalid range")
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Centre of Mass: invalid range")
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        #For contour
        if self.pickDim == 2:
            dataLengthY = self.father.current.len(-2)
            inp = safeEval(self.minEntryY.text(), length=self.father.current.len(), Type='FI')
            if inp is None:
                raise SsnakeException("Centre of Mass: invalid range")
            minimumY = int(round(inp))
            if minimumY < 0:
                minimumY = 0
            elif minimumY > dataLengthY:
                minimumY = dataLengthY
            self.minEntryY.setText(str(minimumY))
            inp = safeEval(self.maxEntryY.text(), length=self.father.current.len(), Type='FI')
            if inp is None:
                raise SsnakeException("Centre of Mass: invalid range")
            maximumY = int(round(inp))
            if maximumY < 0:
                maximumY = 0
            elif maximumY > dataLengthY:
                maximumY = dataLengthY
            self.maxEntryY.setText(str(maximumY))

        if calc:
            if self.pickDim == 1:
                self.comEntry.setText(str(self.father.current.COM([minimum], [maximum], [self.unitDrop.currentIndex()])[0]))
            elif self.pickDim == 2:
                com = self.father.current.COM([minimum, minimumY], [maximum, maximumY], [self.unitDrop.currentIndex(), self.unitDropY.currentIndex()])
                self.comEntry.setText(str(com[0]))
                self.comEntryY.setText(str(com[1]))



##########################################################################################

class IntegralsWindow(wc.ToolWindow):
    NAME = "Integrals"
    CANCELNAME = "&Close"
    OKNAME = "C&alc"
    APPLYANDCLOSE = False

    def __init__(self, parent):
        super(IntegralsWindow, self).__init__(parent)
        self.pickDim = 1
        if isinstance(self.father.current, views.CurrentContour):
            self.pickDim = 2
        self.grid.addWidget(wc.QLabel("Start index X:"), 0, 0)
        self.grid.addWidget(wc.QLabel("End index X:"), 0, 1)
        if self.pickDim == 2:
            self.grid.addWidget(wc.QLabel("Start index Y:"), 0, 2)
            self.grid.addWidget(wc.QLabel("End index Y:"), 0, 3)
        self.grid.addWidget(wc.QLabel("Integral:"), 0, 4)
        self.scaling = 1
        self.num = 0
        self.pickType = 0
        self.minEntries = []
        self.maxEntries = []
        self.minEntriesY = []
        self.maxEntriesY = []
        self.intEntries = []
        self.intValues = []
        self.xValues = []
        self.yValues = []
        self.datMax = 0
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = self.pickDim

    def picked(self, pos):
        if self.pickDim == 2:
            posY = str(pos[3])
        pos = str(pos[0])
        if self.pickType == 0:
            self.minEntries.append(wc.QLineEdit(pos, self.applyFunc))
            self.maxEntries.append(wc.QLineEdit('', self.applyFunc))
            self.intEntries.append(wc.QLineEdit('', (lambda n: lambda: self.setScaling(n))(self.num)))
            self.intValues.append(None)
            self.xValues.append(None)
            self.yValues.append(None)
            self.intEntries[-1].setMinimumWidth(120)
            self.grid.addWidget(self.minEntries[-1], self.num+1, 0)
            self.grid.addWidget(self.maxEntries[-1], self.num+1, 1)
            self.grid.addWidget(self.intEntries[-1], self.num+1, 4)
            if self.pickDim == 2:
                self.minEntriesY.append(wc.QLineEdit(posY, self.applyFunc))
                self.maxEntriesY.append(wc.QLineEdit('', self.applyFunc))
                self.grid.addWidget(self.minEntriesY[-1], self.num+1, 2)
                self.grid.addWidget(self.maxEntriesY[-1], self.num+1, 3)
            self.pickType = 1
        elif self.pickType == 1:
            self.maxEntries[-1].setText(pos)
            if self.pickDim == 2:
                self.maxEntriesY[-1].setText(posY)
            self.num += 1
            self.applyFunc()
            self.pickType = 0
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = self.pickDim

    def preview(self):
        if self.pickDim == 1:
            self.father.current.integralsPreview(self.xValues, self.yValues, self.datMax)
            self.father.current.peakPick = True
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        if self.pickDim == 2:
            xMin = [int(x.text()) for x in self.minEntries]
            xMax = [int(x.text()) for x in self.maxEntries]
            yMin = [int(x.text()) for x in self.minEntriesY]
            yMax = [int(x.text()) for x in self.maxEntriesY]
            self.father.current.integralsPreview(xMin, xMax, yMin, yMax)

    def setScaling(self, num):
        inp = safeEval(self.intEntries[num].text(), length=self.father.current.len(), Type='FI')
        int = self.intValues[num]
        if inp is None:
            return
        else:
            self.scaling = int / inp
        self.applyFunc()

    def applyFunc(self):
        dataLength = [self.father.current.shape()[-1] - 1]
        Parts = [[self.minEntries], [self.maxEntries]]
        if self.pickDim == 2:
            dataLength.append(self.father.current.shape()[-2] - 1)
            Parts[0].append(self.minEntriesY)
            Parts[1].append(self.maxEntriesY)
        for num in range(len(self.minEntries)):
            results = [[], []] #The min/max results
            ok = []
            for place, _ in enumerate(Parts):
                for i, part in enumerate(Parts[place]):
                    inp = safeEval(part[num].text(), length=dataLength, Type='FI')
                    if inp is None:
                        part[num].setText('')
                        ok.append(False)
                    else:
                        ok.append(True)
                        tmp = int(round(inp))
                        tmp = min(max(tmp, 0), dataLength[i]) #makes sure that 0 < value < Length
                        results[place].append(tmp)
                        part[num].setText(str(tmp))
            if all(ok):
                self.intValues[num], self.xValues[num], self.yValues[num], self.datMax = self.father.current.Integrals(*results)
                self.intEntries[num].setText('%#.7g' % (self.intValues[num] / self.scaling))
            else:
                self.intEntries[num].setText('')
                self.intValues[num] = None
        self.preview()

##########################################################################################

class ReorderWindow(wc.ToolWindow):

    NAME = "Reorder"

    def __init__(self, parent):
        super(ReorderWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Positions of the spectra:"), 0, 0)
        self.valEntry = wc.QLineEdit('', self.preview)
        self.grid.addWidget(self.valEntry, 1, 0)
        fileButton = QtWidgets.QPushButton("&Browse")
        fileButton.clicked.connect(self.getPosFromFile)
        self.grid.addWidget(fileButton, 2, 0)
        self.grid.addWidget(wc.QLabel("Length of dimension:"), 3, 0)
        self.lengthEntry = wc.QLineEdit('', self.preview)
        self.grid.addWidget(self.lengthEntry, 4, 0)

    def preview(self, *args):
        pass

    def getPosFromFile(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.lastLocation)
        if isinstance(filename, tuple):
            filename = filename[0]
        if filename:  # if not cancelled
            self.father.father.lastLocation = os.path.dirname(filename)  # Save used path
        if not filename:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def applyFunc(self):
        newLength = self.lengthEntry.text()
        if newLength == '':
            newLength = None
        else:
            newLength = safeEval(self.lengthEntry.text(), length=self.father.current.len(), Type='FI')
            if newLength is None:
                raise SsnakeException("Reorder: `Length' input is not valid")
        val = safeEval(self.valEntry.text(), length=int(self.father.current.len()))
        if not isinstance(val, (list, np.ndarray)):
            raise SsnakeException("Reorder: `Positions' input is not a list or array")
        if len(val) != self.father.current.len():
            raise SsnakeException("Reorder: length of input does not match length of data")
        val = np.array(val, dtype=int)
        check = self.father.current.reorder(val, newLength)
        if check is False:
            raise SsnakeException("Reorder: error during applying")

##########################################################################################


class RegridWindow(wc.ToolWindow):

    NAME = "Regrid"

    def __init__(self, parent):
        super(RegridWindow, self).__init__(parent)
        self.typeDrop = QtWidgets.QComboBox(parent=self)
        self.typeDrop.addItems(["Min/max input"])
        self.grid.addWidget(self.typeDrop, 0, 0, 1, 2)
        # Get unit
        if self.father.current.spec() == 1:
            if self.father.masterData.shape()[self.father.current.axes[-1]] == 1:
                self.closeEvent()
                raise SsnakeException("Regrid: Regrid not possible with size 1")
            if self.father.current.getppm():
                self.unit = 'ppm'
            else:
                axType = self.father.current.getAxType()
                if axType == 0:
                    self.unit = 'Hz'
                elif axType == 1:
                    self.unit = 'kHz'
                elif axType == 2:
                    self.unit = 'MHz'
                elif axType == 3:
                    self.unit = 'ppm'
            maxVal = self.father.current.xax()[-1]
            minVal = self.father.current.xax()[0]
            if self.unit == 'kHz':
                maxVal /= 1e3
                minVal /= 1e3
            elif self.unit == 'MHz':
                maxVal /= 1e6
                minVal /= 1e6
            elif self.unit == 'ppm':
                maxVal /= self.father.masterData.ref[self.father.current.axes[-1]] / 1e6
                minVal /= self.father.masterData.ref[self.father.current.axes[-1]] / 1e6
            self.maxValue = wc.QLineEdit(maxVal)
            self.maxValue.setMinimumWidth(150)
            self.maxLabel = wc.QLeftLabel('Max [' + self.unit + ']:')
            self.minValue = wc.QLineEdit(minVal)
            self.minLabel = wc.QLeftLabel('Min [' + self.unit + ']:')
            self.points = wc.QLineEdit(self.father.masterData.shape()[self.father.current.axes[-1]])
            self.pointsLabel = wc.QLeftLabel('# of points:')
            self.grid.addWidget(self.minValue, 1, 1)
            self.grid.addWidget(self.minLabel, 1, 0)
            self.grid.addWidget(self.maxValue, 2, 1)
            self.grid.addWidget(self.maxLabel, 2, 0)
            self.grid.addWidget(self.pointsLabel, 3, 0)
            self.grid.addWidget(self.points, 3, 1)
        else:
            self.closeEvent()

    def applyFunc(self):
        maxVal = safeEval(self.maxValue.text(), length=self.father.current.len(), Type='FI')
        if maxVal is None:
            raise SsnakeException("Regrid: 'Max' input not valid")
        minVal = safeEval(self.minValue.text(), length=self.father.current.len(), Type='FI')
        if minVal is None:
            raise SsnakeException("Regrid: 'Min' input not valid")
        numPoints = safeEval(self.points.text(), length=self.father.current.len(), Type='FI')
        if numPoints is None or numPoints == 1:
            raise SsnakeException("Regrid: '# of points' input not valid")
        numPoints = int(numPoints)
        # Convert to Hz/s
        if self.unit == 'kHz':
            maxVal *= 1e3
            minVal *= 1e3
        elif self.unit == 'MHz':
            maxVal *= 1e6
            minVal *= 1e6
        elif self.unit == 'ppm':
            maxVal *= self.father.masterData.ref[self.father.current.axes[-1]] / 1e6
            minVal *= self.father.masterData.ref[self.father.current.axes[-1]] / 1e6
        self.father.current.regrid([minVal, maxVal], numPoints)

##########################################################################################


class FFMWindow(wc.ToolWindow):

    NAME = "FFM"

    def __init__(self, parent):
        super(FFMWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Positions of the spectra:"), 0, 0)
        self.valEntry = wc.QLineEdit('', self.preview)
        self.grid.addWidget(self.valEntry, 1, 0)
        fileButton = QtWidgets.QPushButton("&Browse")
        fileButton.clicked.connect(self.getPosFromFile)
        self.grid.addWidget(fileButton, 2, 0)
        self.grid.addWidget(wc.QLabel("Type of the position list:"), 3, 0)
        self.typeDrop = QtWidgets.QComboBox(parent=self)
        self.typeDrop.addItems(["Complex", "States/States-TPPI", "TPPI"])
        self.grid.addWidget(self.typeDrop, 4, 0)
        self.grid.addWidget(wc.QLabel("Reconstruction may take a while"), 5, 0)

    def preview(self, *args):
        pass

    def getPosFromFile(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.lastLocation)
        if isinstance(filename, tuple):
            filename = filename[0]
        if filename:  # if not cancelled
            self.father.father.lastLocation = os.path.dirname(filename)  # Save used path
        if not filename:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def applyFunc(self):
        val = safeEval(self.valEntry.text(), length=self.father.current.len())
        if not isinstance(val, (list, np.ndarray)):
            raise SsnakeException("FFM: 'Positions' is not a list or array")
        val = np.array(val, dtype=int)
        check = self.father.current.ffm(val, self.typeDrop.currentIndex())
        if check is False:
            raise SsnakeException("FFM: error")

##########################################################################################


class CLEANWindow(wc.ToolWindow):

    NAME = "CLEAN"

    def __init__(self, parent):
        super(CLEANWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Positions of the spectra:"), 0, 0)
        self.valEntry = wc.QLineEdit('', self.preview)
        self.grid.addWidget(self.valEntry, 1, 0)
        fileButton = QtWidgets.QPushButton("&Browse")
        fileButton.clicked.connect(self.getPosFromFile)
        self.grid.addWidget(fileButton, 2, 0)
        self.grid.addWidget(wc.QLabel("Type of the position list:"), 3, 0)
        self.typeDrop = QtWidgets.QComboBox(parent=self)
        self.typeDrop.addItems(["Complex", "States/States-TPPI", "TPPI"])
        self.grid.addWidget(self.typeDrop, 4, 0)
        self.grid.addWidget(wc.QLabel("Gamma:"), 5, 0)
        self.gammaEntry = wc.QLineEdit("0.2")
        self.grid.addWidget(self.gammaEntry, 6, 0)
        self.grid.addWidget(wc.QLabel("Threshold:"), 7, 0)
        self.thresholdEntry = wc.QLineEdit("2.0")
        self.grid.addWidget(self.thresholdEntry, 8, 0)
        self.grid.addWidget(wc.QLabel("Max. iterations:"), 11, 0)
        self.maxIterEntry = wc.QLineEdit("2000")
        self.grid.addWidget(self.maxIterEntry, 12, 0)

    def preview(self, *args):
        pass

    def getPosFromFile(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.lastLocation)
        if isinstance(filename, tuple):
            filename = filename[0]
        if filename:  # if not cancelled
            self.father.father.lastLocation = os.path.dirname(filename)  # Save used path
        if not filename:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def applyFunc(self):
        val = safeEval(self.valEntry.text(), length=self.father.current.len())
        if not isinstance(val, (list, np.ndarray)):
            raise SsnakeException("CLEAN: 'Positions' is not a list or array")
        val = np.array(val, dtype=int)
        gamma = safeEval(self.gammaEntry.text(), length=self.father.current.len(), Type='FI')
        if gamma is None:
            raise SsnakeException("CLEAN: 'Gamma' input is not valid")
        threshold = safeEval(self.thresholdEntry.text(), length=self.father.current.len(), Type='FI')
        if threshold is None:
            raise SsnakeException("CLEAN: 'Threshold' input is not valid")
        threshold = threshold
        maxIter = safeEval(self.maxIterEntry.text(), length=self.father.current.len(), Type='FI')
        if maxIter is None:
            raise SsnakeException("CLEAN: 'Max. iter.' is not valid")
        maxIter = int(maxIter)
        check = self.father.current.clean(val, self.typeDrop.currentIndex(), gamma, threshold, maxIter)
        if check is False:
            raise SsnakeException("CLEAN: error")

################################################################


class ISTWindow(wc.ToolWindow):

    NAME = "IST"

    def __init__(self, parent):
        super(ISTWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Positions of the spectra:"), 0, 0)
        self.valEntry = wc.QLineEdit('', self.preview)
        self.grid.addWidget(self.valEntry, 1, 0)
        fileButton = QtWidgets.QPushButton("&Browse")
        fileButton.clicked.connect(self.getPosFromFile)
        self.grid.addWidget(fileButton, 2, 0)
        self.grid.addWidget(wc.QLabel("Type of the position list:"), 3, 0)
        self.typeDrop = QtWidgets.QComboBox(parent=self)
        self.typeDrop.addItems(["Complex", "States/States-TPPI", "TPPI"])
        self.grid.addWidget(self.typeDrop, 4, 0)
        self.grid.addWidget(wc.QLabel("Threshold:"), 5, 0)
        self.thresholdEntry = wc.QLineEdit("0.9")
        self.grid.addWidget(self.thresholdEntry, 6, 0)
        self.grid.addWidget(wc.QLabel("Max. iterations:"), 7, 0)
        self.maxIterEntry = wc.QLineEdit("100")
        self.grid.addWidget(self.maxIterEntry, 8, 0)
        self.grid.addWidget(wc.QLabel("Stop when residual below (% of ND max):"), 9, 0)
        self.tracelimitEntry = wc.QLineEdit("2.0")
        self.grid.addWidget(self.tracelimitEntry, 10, 0)

    def preview(self, *args):
        pass

    def getPosFromFile(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.lastLocation)
        if isinstance(filename, tuple):
            filename = filename[0]
        if filename:  # if not cancelled
            self.father.father.lastLocation = os.path.dirname(filename)  # Save used path
        if not filename:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def applyFunc(self):
        val = safeEval(self.valEntry.text(), length=self.father.current.len())
        if not isinstance(val, (list, np.ndarray)):
            raise SsnakeException("IST: 'Positions' input is not a list or array")
        val = np.array(val, dtype=int)
        tracelimit = safeEval(self.tracelimitEntry.text(), length=self.father.current.len(), Type='FI')
        if tracelimit is None:
            raise SsnakeException("IST: 'Residual' input is not valid")
        tracelimit /= 100
        threshold = safeEval(self.thresholdEntry.text(), length=self.father.current.len(), Type='FI')
        if threshold is None:
            raise SsnakeException("IST: 'Threshold' input is not valid")
        maxIter = safeEval(self.maxIterEntry.text(), length=self.father.current.len(), Type='FI')
        if maxIter is None:
            raise SsnakeException("IST: 'Max. iter.' input is not valid")
        maxIter = int(maxIter)
        check = self.father.current.ist(val, self.typeDrop.currentIndex(), threshold, maxIter, tracelimit)
        if check is False:
            raise SsnakeException("IST: error")

################################################################


class ShearingWindow(wc.ToolWindow):

    NAME = "Shearing"

    def __init__(self, parent):
        super(ShearingWindow, self).__init__(parent)
        options = list(map(str, range(1, self.father.masterData.ndim() + 1)))
        self.grid.addWidget(wc.QLabel("Shearing constant:"), 0, 0)
        self.shearDropdown = QtWidgets.QComboBox()
        
        drop_list = ['User Defined', 
                     'Spin 3/2, 3QMAS', 'Spin 3/2, ST1MAS', 'Spin 3/2, DQ-STMAS',
                     'Spin 5/2, 3QMAS', 'Spin 5/2, 5QMAS', 'Spin 5/2, ST1MAS', 'Spin 5/2, DQ-STMAS',
                     'Spin 7/2, 3QMAS', 'Spin 5/2, 5QMAS', 'Spin 5/2, 7QMAS', 'Spin 7/2, ST1MAS', 'Spin 7/2, DQ-STMAS',
                     'Spin 9/2, 3QMAS', 'Spin 9/2, 5QMAS', 'Spin 9/2, 7QMAS', 'Spin 9/2, 9QMAS', 'Spin 9/2, ST1MAS', 'Spin 9/2, DQ-STMAS',
                    ]
        self.shearDropdown.addItems(drop_list)
        #['User Defined', 'Spin 3/2, -3Q (7/9)', 'Spin 5/2, 3Q (19/12)', 'Spin 5/2, -5Q (25/12)', 'Spin 7/2, 3Q (101/45)',
        #'Spin 7/2, 5Q (11/9)', 'Spin 7/2, -7Q (161/45)', 'Spin 9/2, 3Q (91/36)', 'Spin 9/2, 5Q (95/36)', 'Spin 9/2, 7Q (7/18)', 'Spin 9/2, -9Q (31/6)']
        self.shearDropdown.activated.connect(self.dropdownChanged)
        #self.shearList = [0, 7.0 / 9.0, 19.0 / 12.0, 25.0 / 12.0, 101.0 / 45.0, 11.0 / 9.0, 161.0 / 45.0, 91.0 / 36.0, 95.0 / 36.0, 7.0 / 18.0, 31.0 / 6.0]
        self.scaleList = [1,
              (3/2, -3/2, 3/2), (3/2, -3/2, -1/2), (3/2, -3/2, 1/2),
              (5/2, -3/2, 3/2), (5/2, -5/2, 5/2), (5/2, -3/2, -1/2), (5/2, -3/2, 1/2),
              (7/2, -3/2, 3/2), (7/2, -5/2, 5/2), (7/2, -7/2, 7/2), (7/2, -3/2, -1/2), (7/2, -3/2, 1/2),
              (9/2, -3/2, 3/2), (9/2, -5/2, 5/2), (9/2, -7/2, 7/2), (9/2, -9/2, 9/2), (9/2, -3/2, -1/2), (9/2, -3/2, 1/2),
                        ]
        self.grid.addWidget(self.shearDropdown, 1, 0)
        self.shearEntry = wc.QLineEdit("0.0", self.shearPreview)
        self.grid.addWidget(self.shearEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Shearing direction:"), 4, 0)
        self.dirEntry = QtWidgets.QComboBox()
        self.dirEntry.addItems(options)
        self.dirEntry.setCurrentIndex(self.father.masterData.ndim() - 2)
        self.grid.addWidget(self.dirEntry, 5, 0)
        self.grid.addWidget(wc.QLabel("Shearing axis:"), 6, 0)
        self.axEntry = QtWidgets.QComboBox()
        self.axEntry.addItems(options)
        self.axEntry.setCurrentIndex(self.father.masterData.ndim() - 1)
        self.grid.addWidget(self.axEntry, 7, 0)
        self.toRefCheck = QtWidgets.QCheckBox("Relative to Reference")
        self.grid.addWidget(self.toRefCheck, 8, 0)

    def dropdownChanged(self):
        index = self.shearDropdown.currentIndex()
        if index == 0:
            scale = "1"
        else:
            scale = f"{-func.R(*self.scaleList[index])}"
        self.shearEntry.setText(scale)
#        self.shearEntry.setText("%.9f" % self.shearList[index])

    def shearPreview(self, *args):
        shear = safeEval(self.shearEntry.text(), length=self.father.current.len(), Type='FI')
        if shear is not None:
            self.shearEntry.setText(str(float(shear)))

    def applyFunc(self):
        shear = safeEval(self.shearEntry.text(), length=self.father.current.len(), Type='FI')
        if shear is None:
            raise SsnakeException("Shearing: 'constant' not a valid value")
        axis = self.dirEntry.currentIndex()
        axis2 = self.axEntry.currentIndex()
        if axis == axis2:
            raise SsnakeException("Shearing: axes cannot be the same for shearing")
        self.father.current.shearing(float(shear), axis, axis2, self.toRefCheck.isChecked())

##########################################################################################


class MultiplyWindow(wc.ToolWindow):

    NAME = "Multiply"
    SINGLESLICE = True

    def __init__(self, parent):
        super(MultiplyWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Values:"), 0, 0)
        self.valEntry = wc.QLineEdit('', self.preview)
        self.grid.addWidget(self.valEntry, 1, 0)

    def preview(self, *args):
        val = safeEval(self.valEntry.text(), length=self.father.current.len())
        if val is None:
            raise SsnakeException("Multiply: input not valid")
        self.father.current.multiplyPreview(np.array(val))

    def applyFunc(self):
        val = safeEval(self.valEntry.text(), length=self.father.current.len())
        if val is None:
            raise SsnakeException("Multiply: input not valid")
        self.father.current.multiply(np.array(val), self.singleSlice.isChecked())

##########################################################################################

class NormalizeWindow(wc.ToolWindow):

    NAME = "Normalize"
    SINGLESLICE = True

    def __init__(self, parent):
        super(NormalizeWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Start index:"), 0, 0)
        self.minEntry = wc.QLineEdit("0", self.checkValues)
        self.grid.addWidget(self.minEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End index:"), 2, 0)
        self.maxEntry = wc.QLineEdit(parent.current.len(), self.checkValues)
        self.grid.addWidget(self.maxEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Type:"), 4, 0)
        self.typeDrop = QtWidgets.QComboBox()
        self.typeDrop.addItems(['Integral', 'Maximum', 'Minimum'])
        self.typeDrop.setCurrentIndex(0)
        self.typeDrop.currentIndexChanged.connect(self.checkValues)
        self.grid.addWidget(self.typeDrop, 5, 0)
        self.grid.addWidget(wc.QLabel("Multiplier:"), 6, 0)
        self.valEntry = wc.QLineEdit("1.0")
        self.grid.addWidget(self.valEntry, 7, 0)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def picked(self, pos, num=0):
        if num == 0:
            self.minEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 1)
            self.father.current.peakPick = True
        elif num == 1:
            self.maxEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 0)
            self.father.current.peakPick = True
            #self.applyFunc()

    def checkValues(self, *args):
        dataLength = self.father.current.len()
        inp = safeEval(self.minEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        #self.applyFunc()

    def applyFunc(self):
        dataLength = self.father.current.len()
        inp = safeEval(self.minEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Normalize: invalid range")
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text(), length=self.father.current.len(), Type='FI')
        if inp is None:
            raise SsnakeException("Normalize: invalid range")
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        try:
            scale = safeEval(self.valEntry.text(), length=self.father.current.len(), Type='FI')
        except Exception:
            raise SsnakeException("Normalize: invalid multiplier")
        type = self.typeDrop.currentIndex()
        if type == 0:
            val, xValues, yValues, datMax = self.father.current.Integrals([minimum], [maximum])
        elif type == 1:
            val = self.father.current.MaxMin(minimum, maximum, type='max')
        elif type == 2:
            val = self.father.current.MaxMin(minimum, maximum, type='min')
        self.father.current.normalize( 1.0 / val, scale, type, self.singleSlice.isChecked())

##########################################################################################

class XaxWindow(wc.ToolWindow):

    RESIZABLE = True
    NAME = "User defined x-axis"

    def __init__(self, parent):
        super(XaxWindow, self).__init__(parent)
        self.axisSize = self.father.current.len()
        self.grid.addWidget(wc.QLabel("Input x-axis values:"), 0, 0, 1, 2)
        self.typeDropdown = QtWidgets.QComboBox()
        self.typeDropdown.addItems(['Expression', 'Linear', 'Logarithmic'])
        self.typeDropdown.activated.connect(self.typeChanged)
        self.grid.addWidget(self.typeDropdown, 1, 0, 1, 2)
        self.exprEntry = wc.QLineEdit('', self.xaxPreview)
        self.grid.addWidget(self.exprEntry, 2, 0, 1, 2)
        # Linear
        self.linStartLabel = wc.QLeftLabel("Start [s]:")
        self.linStopLabel = wc.QLeftLabel("Stop [s]:")
        self.linStartLabel.hide()
        self.linStopLabel.hide()
        self.grid.addWidget(self.linStartLabel, 3, 0, 1, 1)
        self.grid.addWidget(self.linStopLabel, 4, 0, 1, 1)
        self.linStartEntry = wc.QLineEdit('', self.xaxPreview)
        self.linStartEntry.setMaximumWidth(120)
        self.linStartEntry.hide()
        self.grid.addWidget(self.linStartEntry, 3, 1, 1, 1)
        self.linStopEntry = wc.QLineEdit('', self.xaxPreview)
        self.linStopEntry.setMaximumWidth(120)
        self.linStopEntry.hide()
        self.grid.addWidget(self.linStopEntry, 4, 1, 1, 1)
        # Log
        self.logStartLabel = wc.QLeftLabel("Start [s]:")
        self.logStopLabel = wc.QLeftLabel("Stop [s]:")
        self.logStartLabel.hide()
        self.logStopLabel.hide()
        self.grid.addWidget(self.logStartLabel, 5, 0, 1, 1)
        self.grid.addWidget(self.logStopLabel, 6, 0, 1, 1)
        self.logStartEntry = wc.QLineEdit('', self.xaxPreview)
        self.logStartEntry.setMaximumWidth(120)
        self.logStartEntry.hide()
        self.grid.addWidget(self.logStartEntry, 5, 1, 1, 1)
        self.logStopEntry = wc.QLineEdit('', self.xaxPreview)
        self.logStopEntry.setMaximumWidth(120)
        self.logStopEntry.hide()
        self.grid.addWidget(self.logStopEntry, 6, 1, 1, 1)
        self.table = QtWidgets.QTableWidget(self.axisSize, 2)
        self.table.setHorizontalHeaderLabels(['Index', 'Value [s]'])
        self.table.verticalHeader().hide()
        for val in range(self.axisSize):
            item = QtWidgets.QTableWidgetItem(str(val))
            item.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table.setItem(int(val), 0, item)
            item2 = QtWidgets.QTableWidgetItem('')
            item2.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table.setItem(int(val), 1, item2)
#        self.table.setVerticalHeaderLabels([str(a) for a in range(self.axisSize)])
        self.grid.addWidget(self.table, 12, 0, 1, 2)
        self.resize(250, 500)

    def typeChanged(self, index):
        if index == 0:  # If expr
            self.exprEntry.show()
            self.linStartLabel.hide()
            self.linStopLabel.hide()
            self.linStartEntry.hide()
            self.linStopEntry.hide()
            self.logStartLabel.hide()
            self.logStopLabel.hide()
            self.logStartEntry.hide()
            self.logStopEntry.hide()
        elif index == 1:
            self.exprEntry.hide()
            self.linStartLabel.show()
            self.linStopLabel.show()
            self.linStartEntry.show()
            self.linStopEntry.show()
            self.logStartLabel.hide()
            self.logStopLabel.hide()
            self.logStartEntry.hide()
            self.logStopEntry.hide()
        elif index == 2:
            self.exprEntry.hide()
            self.linStartLabel.hide()
            self.linStopLabel.hide()
            self.linStartEntry.hide()
            self.linStopEntry.hide()
            self.logStartLabel.show()
            self.logStopLabel.show()
            self.logStartEntry.show()
            self.logStopEntry.show()

    def getValues(self):
        if self.typeDropdown.currentIndex() == 0:
            env = vars(np).copy()
            env['length'] = self.father.current.len()  # so length can be used to in equations
            env['euro'] = lambda fVal, num=self.axisSize: func.euro(fVal, num)
            try:
                val = np.array(eval(self.exprEntry.text(), env),dtype=float)                # find a better solution, also add catch for exceptions
            except SyntaxError:
                try:
                    val = np.fromstring(self.exprEntry.text(), sep=' ')
                    val2 = np.fromstring(self.exprEntry.text(), sep=',')
                    if len(val2) > len(val):
                        val = val2
                except Exception as e:
                    raise SsnakeException(str(e))
            except Exception as e:
                raise SsnakeException(str(e))
            if not isinstance(val, (list, np.ndarray)):
                raise SsnakeException("X-axis: Input is not a list or array")
            if len(val) != self.father.current.len():
                raise SsnakeException("X-axis: Length of input does not match length of data")
            if not all(isinstance(x, (int, float)) for x in val):
                raise SsnakeException("X-axis: Array is not all of int or float type")
        elif self.typeDropdown.currentIndex() == 1:
            start = safeEval(self.linStartEntry.text(), Type='FI')
            stop = safeEval(self.linStopEntry.text(), Type='FI')
            if start is None:
                raise SsnakeException("X-axis: linear start value is not valid")
            if stop is None:
                raise SsnakeException("X-axis: linear stop value is not valid")
            val = np.linspace(start, stop, self.axisSize)
        elif self.typeDropdown.currentIndex() == 2:
            start = safeEval(self.logStartEntry.text(), Type='FI')
            stop = safeEval(self.logStopEntry.text(), Type='FI')
            if start is None or start <= 0.0:
                raise SsnakeException("X-axis: logarithmic start value is not valid")
            if stop is None or stop <= 0.0:
                raise SsnakeException("X-axis: logarithmic stop value is not valid")
            val = np.logspace(np.log10(start), np.log10(stop), self.axisSize)
        return val

    def xaxPreview(self, *args):
        val = self.getValues()
        if val is None:  # if error return. Messages are handled by the called function
            return
        for i in range(self.axisSize):
            item = QtWidgets.QTableWidgetItem('{:.6g}'.format(val[i]))
            item.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table.setItem(i, 1, item)
        self.father.current.setXaxPreview(np.array(val))

    def applyFunc(self):
        val = self.getValues()
        if val is None:  # if error return. Messages are handled by the called function
            return
        self.father.current.setXax(np.array(val))

##########################################################################################


class RefWindow(wc.ToolWindow):

    NAME = "Reference"

    def __init__(self, parent):
        super(RefWindow, self).__init__(parent)
        # Secondary reference definitions
        file = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + "References.txt"
        with open(file) as refFile:
            refList = [line.strip().split('\t') for line in refFile]
        secRefNames = ["User Defined"]
        secRefValues = ["0.0"]
        for entry in refList:
            secRefNames.append(entry[0])
            secRefValues.append(entry[1])
        self.secRefNames = secRefNames
        self.secRefValues = secRefValues
        if parent.current.spec() == 0:
            self.closeEvent()
            raise SsnakeException('Setting ppm is only available for frequency data')
        self.grid.addWidget(wc.QLabel("Name:"), 0, 0)
        self.refName = wc.QLineEdit()
        self.grid.addWidget(self.refName, 1, 0)
        self.grid.addWidget(wc.QLabel("Frequency [MHz]:"), 2, 0)
        self.freqEntry = wc.QLineEdit(("%.7f" % (self.father.current.ref() * 1e-6)), self.preview)
        self.grid.addWidget(self.freqEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Secondary Reference:"), 4, 0)
        self.refSecond = QtWidgets.QComboBox(parent=self)
        self.refSecond.addItems(self.secRefNames)
        self.refSecond.activated.connect(self.fillSecondaryRef)
        self.grid.addWidget(self.refSecond, 5, 0)
        self.grid.addWidget(wc.QLabel("Reference [ppm]:"), 6, 0)
        self.refEntry = wc.QLineEdit("0.0", self.preview)
        self.grid.addWidget(self.refEntry, 7, 0)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def preview(self, *args):
        freq = safeEval(self.freqEntry.text(), length=self.father.current.len(), Type='FI')
        ref = safeEval(self.refEntry.text(), length=self.father.current.len(), Type='FI')
        if freq is None or ref is None:
            return
        self.freqEntry.setText("%.7f" % (freq))
        self.refEntry.setText(str(ref))

    def fillSecondaryRef(self):
        self.refEntry.setText(self.secRefValues[self.refSecond.currentIndex()])

    def applyAndClose(self):
        self.father.current.peakPickReset()
        freq = safeEval(self.freqEntry.text(), length=self.father.current.len(), Type='FI')
        ref = safeEval(self.refEntry.text(), length=self.father.current.len(), Type='FI')
        if freq is None or ref is None:
            raise SsnakeException("Not a valid value")
        freq = freq * 1e6
        reffreq = freq / (1.0 + ref * 1e-6)
        givenname = self.refName.text()
        nameOK = True
        if givenname:  # If name is filled in
            if givenname in self.father.father.referenceName:  # if exists
                self.father.father.dispMsg("Reference name '" + givenname + "' already exists")
                nameOK = False
            else:
                self.father.father.referenceAdd(reffreq, givenname)
        if nameOK:
            self.father.current.setRef(reffreq)
            self.father.bottomframe.upd()
            self.closeEvent()

    def picked(self, pos):
        self.freqEntry.setText("%.7f" % ((self.father.current.ref() + self.father.current.xax()[pos[0]]) * 1e-6))
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

##########################################################################################


class HistoryWindow(wc.ToolWindow):

    NAME = "Processing history"
    RESIZABLE = True
    MENUDISABLE = True

    def __init__(self, parent):
        super(HistoryWindow, self).__init__(parent)
        self.cancelButton.hide()
        self.valEntry = QtWidgets.QTextEdit()
        self.valEntry.setReadOnly(True)
        self.valEntry.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.valEntry.setText(self.father.masterData.getHistory())
        self.grid.addWidget(self.valEntry, 1, 0)
        self.resize(550, 700)

#########################################################################################


class OrigListWidget(QtWidgets.QListWidget):

    def __init__(self, parent=None, dest=None):
        super(OrigListWidget, self).__init__(parent)
        self.setDragDropMode(QtWidgets.QAbstractItemView.DragDrop)
        self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.setAcceptDrops(True)
        self.dest = dest

    def dropEvent(self, event):
        if event.source() == self:
            pass
        else:
            if self.dest is not None:
                for item in self.dest.selectedItems():
                    self.dest.takeItem(self.dest.row(item))

    def mouseDoubleClickEvent(self, event):
        for item in self.selectedItems():
            QtWidgets.QListWidgetItem(item.text(), self.dest)


#########################################################################################


class DestListWidget(QtWidgets.QListWidget):

    def __init__(self, parent=None):
        super(DestListWidget, self).__init__(parent)
        self.setDragDropMode(QtWidgets.QAbstractItemView.DragDrop)
        self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.setAcceptDrops(True)

    def dropEvent(self, event):
        if event.source() == self:
            event.setDropAction(QtCore.Qt.MoveAction)
            super(DestListWidget, self).dropEvent(event)
        else:
            event.setDropAction(QtCore.Qt.CopyAction)
            super(DestListWidget, self).dropEvent(event)

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Delete:
            self.deleteSelected()

    def deleteSelected(self):
        for item in self.selectedItems():
            self.takeItem(self.row(item))

    def moveSelection(self,direction = 'up'):

        #Get selected items
        index = [self.row(item) for item in self.selectedItems()]
        items = [item for item in self.selectedItems()]
        #Sort items based on index, to get move order right
        items = [x for _,x in sorted(zip(index,items))]

        if direction == 'up':
            check = 0
            step = -1
        elif direction == 'down':
            check = self.count() - 1
            step = +1
            items = items[::-1] #Invert move order

        if check in index: #If one item already at limit
            return

        #If not, move one line
        for item in items:
            row = self.row(item)
            currentItem = self.takeItem(row)
            self.insertItem(row + step, currentItem)

        #Reselect the items
        for item in items:
            item.setSelected(True)

    def mouseDoubleClickEvent(self, event):
        self.deleteSelected()

##########################################################################################


class CombineWorkspaceWindow(wc.ToolWindow):

    NAME = "Combine workspaces"
    RESIZABLE = True

    def __init__(self, parent):
        super(CombineWorkspaceWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Workspaces:"), 0, 0)
        self.grid.addWidget(wc.QLabel("Combined spectrum:"), 0, 2)
        self.listB = DestListWidget(self)
        self.listA = OrigListWidget(self, self.listB)
        for i in self.father.workspaceNames:
            QtWidgets.QListWidgetItem(i, self.listA).setToolTip(i)
        self.grid.addWidget(self.listA, 1, 0, 2, 1)
        self.grid.addWidget(self.listB, 1, 2, 2, 1)
        self.rightPush = QtWidgets.QPushButton(u"\u2192", self)
        self.leftPush = QtWidgets.QPushButton(u"\u2190", self)
        self.downPush = QtWidgets.QPushButton(u"\u2193", self)
        self.upPush = QtWidgets.QPushButton(u"\u2191", self)
        self.rightPush.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Expanding)
        self.leftPush.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Expanding)
        self.downPush.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Expanding)
        self.upPush.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Expanding)
        self.grid.addWidget(self.rightPush, 1, 1)
        self.grid.addWidget(self.leftPush, 2, 1)
        self.grid.addWidget(self.upPush, 1, 3)
        self.grid.addWidget(self.downPush, 2, 3)
        self.leftPush.clicked.connect(self.right2left)
        self.rightPush.clicked.connect(self.left2right)
        self.upPush.clicked.connect(self.moveUp)
        self.downPush.clicked.connect(self.moveDown)
        self.resize(500, 400)

    def right2left(self):
        self.listB.deleteSelected()

    def left2right(self):
        for item in self.listA.selectedItems():
            self.listB.addItem(item.text())

    def moveUp(self):
        self.listB.moveSelection('up')

    def moveDown(self):
        self.listB.moveSelection('down')

    def applyFunc(self, *args):
        items = []
        for index in range(self.listB.count()):
            items.append(self.listB.item(index).text())
        if not items:
            raise SsnakeException("Please select at least one workspace to combine")
        self.father.combineWorkspace(items)

    def closeEvent(self, *args):
        if self.MENUDISABLE:
            self.father.menuEnable(True)
        self.deleteLater()

##########################################################################################


class CombineLoadWindow(wc.ToolWindow):

    NAME = "Open & Combine"
    BROWSE = True
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent):
        super(CombineLoadWindow, self).__init__(parent)
        self.setAcceptDrops(True)
        self.grid.addWidget(wc.QLabel("Data to be Combined:"), 0, 0)
        self.specList = DestListWidget(self)
        self.grid.addWidget(self.specList, 1, 0)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            name = url.toLocalFile()
            self.specList.addItem(name)

    def browse(self):
        fileList = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open File', self.father.lastLocation)
        if isinstance(fileList, tuple):
            fileList = fileList[0]
        for filePath in fileList:
            if filePath:  # if not cancelled
                self.father.lastLocation = os.path.dirname(filePath)  # Save used path
            if not filePath:
                return
            self.specList.addItem(filePath)

    def applyFunc(self, *args):
        items = []
        for index in range(self.specList.count()):
            items.append(self.specList.item(index).text())
        if not items:
            raise SsnakeException("Please select at least one workspace to combine")
        self.father.loadAndCombine(items)

    def closeEvent(self, *args):
        self.deleteLater()

##########################################################################################


class MonitorWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        super(MonitorWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Monitor")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        fileName = self.father.masterData.filePath[0][0]
        if len(fileName) > 58:
            fileName = fileName[:55] + '...'
        fileLabel = wc.QLabel("File: " + fileName)
        fileLabel.setToolTip(self.father.masterData.filePath[0][0])
        layout.addWidget(fileLabel, 0, 0, 1, 3)
        layout.addLayout(grid, 1, 0, 1, 3)
        grid.addWidget(wc.QLabel("Macros:"), 0, 0)
        grid.addWidget(wc.QLabel("Apply after loading:"), 0, 1)
        self.listB = DestListWidget(self)
        for i in self.father.monitorMacros:
            QtWidgets.QListWidgetItem(i, self.listB).setToolTip(i)
        self.listA = OrigListWidget(self, self.listB)
        for i in self.father.father.macros.keys():
            QtWidgets.QListWidgetItem(i, self.listA).setToolTip(i)
        grid.addWidget(self.listA, 1, 0)
        grid.addWidget(self.listB, 1, 1)
        grid.addWidget(wc.QLabel("Delay [s]:"), 2, 0)
        self.delTime = wc.SsnakeDoubleSpinBox()
        self.delTime.setMaximum(10000)
        self.delTime.setMinimum(0)
        self.delTime.setSingleStep(0.1)
        self.delTime.setValue(0.5)
        grid.addWidget(self.delTime, 2, 1)
        cancelButton = QtWidgets.QPushButton("&Close")
        cancelButton.clicked.connect(self.closeEvent)
        watchButton = QtWidgets.QPushButton("&Watch")
        watchButton.clicked.connect(self.applyAndClose)
        unwatchButton = QtWidgets.QPushButton("&Unwatch")
        unwatchButton.clicked.connect(self.stopAndClose)

        box = QtWidgets.QDialogButtonBox()
        box.addButton(cancelButton, QtWidgets.QDialogButtonBox.RejectRole)
        box.addButton(watchButton, QtWidgets.QDialogButtonBox.ActionRole)
        box.addButton(unwatchButton, QtWidgets.QDialogButtonBox.ActionRole)
        layout.addWidget(box, 3, 0)
        layout.setColumnStretch(4, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuEnable(False)
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def applyAndClose(self, *args):
        self.father.stopMonitor()
        items = []
        for index in range(self.listB.count()):
            items.append(self.listB.item(index).text())
        delay = self.delTime.value()
        self.father.startMonitor(items, delay)
        self.closeEvent()

    def stopAndClose(self, *args):
        self.father.stopMonitor()
        self.closeEvent()

    def closeEvent(self, *args):
        self.father.menuEnable(True)
        self.deleteLater()

##############################################################################


class PlotSettingsWindow(wc.ToolWindow):

    NAME = "Preferences"

    def __init__(self, parent):
        super(PlotSettingsWindow, self).__init__(parent)
        tabWidget = QtWidgets.QTabWidget()
        tab1 = QtWidgets.QWidget()
        tab2 = QtWidgets.QWidget()
        tab3 = QtWidgets.QWidget()
        tabWidget.addTab(tab1, "Plot")
        tabWidget.addTab(tab2, "Contour")
        tabWidget.addTab(tab3, "2D Colour")
        grid1 = QtWidgets.QGridLayout()
        grid2 = QtWidgets.QGridLayout()
        grid3 = QtWidgets.QGridLayout()
        tab1.setLayout(grid1)
        tab2.setLayout(grid2)
        tab3.setLayout(grid3)
        grid1.setColumnStretch(10, 1)
        grid1.setRowStretch(10, 1)
        grid2.setColumnStretch(10, 1)
        grid2.setRowStretch(10, 1)
        grid3.setColumnStretch(10, 1)
        grid3.setRowStretch(10, 1)
        grid1.addWidget(QtWidgets.QLabel("Linewidth:"), 1, 0)
        self.lwSpinBox = wc.SsnakeDoubleSpinBox()
        self.lwSpinBox.setSingleStep(0.1)
        self.lwSpinBox.setValue(self.father.current.viewSettings["linewidth"])
        self.lwSpinBox.valueChanged.connect(self.preview)
        grid1.addWidget(self.lwSpinBox, 1, 1)
        self.color = self.father.current.viewSettings["color"]
        lineColorButton = QtWidgets.QPushButton("Line colour")
        lineColorButton.clicked.connect(self.setColor)
        grid1.addWidget(lineColorButton, 2, 0)
        grid1.addWidget(QtWidgets.QLabel("Colour range:"), 3, 0)
        self.crEntry = QtWidgets.QComboBox(self)
        self.crEntry.addItems(views.COLORRANGELIST)
        self.crEntry.setCurrentIndex(self.father.current.getColorRange())
        self.crEntry.currentIndexChanged.connect(self.preview)
        grid1.addWidget(self.crEntry, 3, 1)
        self.xgridCheck = QtWidgets.QCheckBox("x-grid")
        self.xgridCheck.setChecked(self.father.current.viewSettings["grids"][0])
        self.xgridCheck.stateChanged.connect(self.preview)
        grid1.addWidget(self.xgridCheck, 4, 0, 1, 2)
        self.ygridCheck = QtWidgets.QCheckBox("y-grid")
        self.ygridCheck.setChecked(self.father.current.viewSettings["grids"][1])
        grid1.addWidget(self.ygridCheck, 5, 0, 1, 2)
        self.ygridCheck.stateChanged.connect(self.preview)
        grid1.addWidget(QtWidgets.QLabel("Min X Ticks:"), 6, 0)
        self.xTicksSpinBox = wc.SsnakeSpinBox()
        self.xTicksSpinBox.setValue(self.father.current.viewSettings["minXTicks"])
        self.xTicksSpinBox.valueChanged.connect(self.preview)
        grid1.addWidget(self.xTicksSpinBox, 6, 1)
        grid1.addWidget(QtWidgets.QLabel("Min Y Ticks:"), 7, 0)
        self.yTicksSpinBox = wc.SsnakeSpinBox()
        self.yTicksSpinBox.setValue(self.father.current.viewSettings["minYTicks"])
        self.yTicksSpinBox.valueChanged.connect(self.preview)
        grid1.addWidget(self.yTicksSpinBox, 7, 1)
        grid2.addWidget(QtWidgets.QLabel("Colourmap:"), 0, 0)
        self.cmEntry = QtWidgets.QComboBox(self)
        self.cmEntry.addItems(views.COLORMAPLIST)
        self.cmEntry.setCurrentIndex(self.father.current.getColorMap())
        self.cmEntry.currentIndexChanged.connect(self.preview)
        grid2.addWidget(self.cmEntry, 0, 1)
        self.constColorCheck = QtWidgets.QCheckBox("Constant colours")
        self.constColorCheck.setChecked(self.father.current.viewSettings["contourConst"])
        grid2.addWidget(self.constColorCheck, 1, 0)
        self.constColorCheck.stateChanged.connect(self.preview)
        self.posColor = self.father.current.viewSettings["contourColors"][0]
        posColorButton = QtWidgets.QPushButton("Positive colour")
        posColorButton.clicked.connect(self.setPosColor)
        grid2.addWidget(posColorButton, 2, 0)
        self.negColor = self.father.current.viewSettings["contourColors"][1]
        negColorButton = QtWidgets.QPushButton("Negative colour")
        negColorButton.clicked.connect(self.setNegColor)
        grid2.addWidget(negColorButton, 3, 0)
        grid3.addWidget(QtWidgets.QLabel("Colourmap:"), 0, 0)
        self.cmEntry2D = QtWidgets.QComboBox(self)
        self.cmEntry2D.addItems(views.COLORMAPLIST)
        self.cmEntry2D.setCurrentIndex(self.father.current.getPColorMap())
        self.cmEntry2D.currentIndexChanged.connect(self.preview)
        grid3.addWidget(self.cmEntry2D, 0, 1)
        self.grid.addWidget(tabWidget, 0, 0)

    def preview(self, *args):
        tmpLw = self.father.current.viewSettings["linewidth"]
        self.father.current.setLw(self.lwSpinBox.value())
        tmpXTicks = self.father.current.viewSettings["minXTicks"]
        tmpYTicks = self.father.current.viewSettings["minYTicks"]
        self.father.current.setTickNum(self.xTicksSpinBox.value(), self.yTicksSpinBox.value())
        tmpColor = self.father.current.viewSettings["color"]
        self.father.current.setColor(self.color)
        tmpColorRange = self.father.current.getColorRange()
        self.father.current.setColorRange(self.crEntry.currentIndex())
        tmpColorMap = self.father.current.getColorMap()
        self.father.current.setColorMap(self.cmEntry.currentIndex())
        tmpGrids = self.father.current.viewSettings["grids"]
        self.father.current.setGrids([self.xgridCheck.isChecked(), self.ygridCheck.isChecked()])
        tmpContourConst = self.father.current.viewSettings["contourConst"]
        self.father.current.setContourConst(self.constColorCheck.isChecked())
        tmpContourColors = self.father.current.viewSettings["contourColors"]
        tmpColorMap2D = self.father.current.getPColorMap()
        self.father.current.setPColorMap(self.cmEntry2D.currentIndex())
        self.father.current.setContourColors([self.posColor, self.negColor])
        self.father.current.showFid()
        self.father.current.setLw(tmpLw)
        self.father.current.setTickNum(tmpXTicks, tmpYTicks)
        self.father.current.setColor(tmpColor)
        self.father.current.setColorRange(tmpColorRange)
        self.father.current.setColorMap(tmpColorMap)
        self.father.current.setGrids(tmpGrids)
        self.father.current.setContourConst(tmpContourConst)
        self.father.current.setContourColors(tmpContourColors)
        self.father.current.setPColorMap(tmpColorMap2D)

    def setColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtGui.QColor(self.color))
        if tmp.isValid():
            self.color = tmp.name()
        self.preview()

    def setPosColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtGui.QColor(self.posColor))
        if tmp.isValid():
            self.posColor = tmp.name()
        self.preview()

    def setNegColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtGui.QColor(self.negColor))
        if tmp.isValid():
            self.negColor = tmp.name()
        self.preview()

    def applyFunc(self, *args):
        self.father.current.setColor(self.color)
        self.father.current.setLw(self.lwSpinBox.value())
        self.father.current.setTickNum(self.xTicksSpinBox.value(), self.yTicksSpinBox.value())
        self.father.current.setGrids([self.xgridCheck.isChecked(), self.ygridCheck.isChecked()])
        self.father.current.setColorRange(self.crEntry.currentIndex())
        self.father.current.setColorMap(self.cmEntry.currentIndex())
        self.father.current.setContourConst(self.constColorCheck.isChecked())
        self.father.current.setContourColors([self.posColor, self.negColor])
        self.father.current.setPColorMap(self.cmEntry2D.currentIndex())

##############################################################################


class errorWindow(wc.ToolWindow):

    NAME = "Error Messages"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent):
        super(errorWindow, self).__init__(parent)
        self.cancelButton.hide()
        self.errorQList = QtWidgets.QListWidget(self)
        self.errorQList.currentRowChanged.connect(self.rowChange)
        for error in self.father.errors:
            if len(error[1]) == 3:
                tmp = QtWidgets.QListWidgetItem(error[0] + ': Program error. Please report.', self.errorQList)
                tmp.setForeground(QtGui.QBrush(QtGui.QColor('red')))
            elif len(error[1]) == 1:
                QtWidgets.QListWidgetItem(error[0] + ': ' + error[1][0], self.errorQList)
        self.errorEdit = QtWidgets.QTextEdit(self)
        self.errorEdit.setReadOnly(True)
        errorText = ''
        self.errorEdit.setHtml(errorText)
        self.grid.addWidget(self.errorQList, 0, 0, 1, 3)
        self.grid.addWidget(self.errorEdit, 1, 0, 1, 3)
        self.resize(550, 700)

    def rowChange(self, row):
        errorText = ''
        error = self.father.errors[row]
        if len(error[1]) == 3:
            errorText = errorText + error[0] + '<br>'
            for line in tb.format_exception(error[1][0], error[1][1], error[1][2]):
                errorText = errorText + line + '<br>'
        self.errorEdit.setHtml(errorText)

    def closeEvent(self, *args):
        self.deleteLater()

##############################################################################


class PreferenceWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        super(PreferenceWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Preferences")
        tabWidget = QtWidgets.QTabWidget()
        tab1 = QtWidgets.QWidget()
        tab2 = QtWidgets.QWidget()
        tab3 = QtWidgets.QWidget()
        tab4 = QtWidgets.QWidget()
        tab5 = QtWidgets.QWidget()
        tabWidget.addTab(tab1, "Window")
        tabWidget.addTab(tab2, "Plot")
        tabWidget.addTab(tab3, "Contour")
        tabWidget.addTab(tab4, "2D Colour")
        tabWidget.addTab(tab5, "Phasing")
        grid1 = QtWidgets.QGridLayout()
        grid2 = QtWidgets.QGridLayout()
        grid3 = QtWidgets.QGridLayout()
        grid4 = QtWidgets.QGridLayout()
        grid5 = QtWidgets.QGridLayout()
        tab1.setLayout(grid1)
        tab2.setLayout(grid2)
        tab3.setLayout(grid3)
        tab4.setLayout(grid4)
        tab5.setLayout(grid5)
        grid1.setColumnStretch(10, 1)
        grid1.setRowStretch(10, 1)
        grid2.setColumnStretch(10, 1)
        grid2.setRowStretch(10, 1)
        grid3.setColumnStretch(10, 1)
        grid3.setRowStretch(10, 1)
        grid4.setColumnStretch(10, 1)
        grid4.setRowStretch(10, 1)
        grid5.setColumnStretch(10, 1)
        grid5.setRowStretch(10, 1)
        # grid1.addWidget(wc.QLabel("Window size:"), 0, 0, 1, 2)
        grid1.addWidget(wc.QLabel("Width:"), 1, 0)
        self.widthSpinBox = wc.SsnakeSpinBox()
        self.widthSpinBox.setMaximum(100000)
        self.widthSpinBox.setMinimum(1)
        self.widthSpinBox.setValue(self.father.defaultWidth)
        grid1.addWidget(self.widthSpinBox, 1, 1)
        grid1.addWidget(wc.QLabel("Height:"), 2, 0)
        self.heightSpinBox = wc.SsnakeSpinBox()
        self.heightSpinBox.setMaximum(100000)
        self.heightSpinBox.setMinimum(1)
        self.heightSpinBox.setValue(self.father.defaultHeight)
        grid1.addWidget(self.heightSpinBox, 2, 1)
        self.maximizedCheck = QtWidgets.QCheckBox("Open maximized")
        self.maximizedCheck.setChecked(self.father.defaultMaximized)
        grid1.addWidget(self.maximizedCheck, 3, 0, 1, 2)
        self.askNameCheck = QtWidgets.QCheckBox("Ask workspace name when loading")
        self.askNameCheck.setChecked(self.father.defaultAskName)
        grid1.addWidget(self.askNameCheck, 4, 0, 1, 2)
        self.toolbarCheck = QtWidgets.QCheckBox("Show Shortcut Toolbar")
        self.toolbarCheck.setChecked(self.father.defaultToolBar)
        grid1.addWidget(self.toolbarCheck, 5, 0, 1, 2)
        self.tooltipCheck = QtWidgets.QCheckBox("Show Tooltips")
        self.tooltipCheck.setChecked(self.father.defaultTooltips)
        grid1.addWidget(self.tooltipCheck, 7, 0, 1, 2)
        editToolbarButton = QtWidgets.QPushButton("Edit Toolbar")
        editToolbarButton.clicked.connect(lambda: ToolbarWindow(self))
        grid1.addWidget(editToolbarButton, 6, 0, 1, 2)
        self.currentToolbar = self.father.defaultToolbarActionList
        self.startupgroupbox = QtWidgets.QGroupBox("Startup Directory")
        self.startupgroupbox.setCheckable(True)
        self.startupgroupbox.setChecked(self.father.defaultStartupBool)
        grid1.addWidget(self.startupgroupbox, 8, 0, 1, 2)
        startupgrid = QtWidgets.QGridLayout()
        self.startupgroupbox.setLayout(startupgrid)
        self.startupDirEntry = QtWidgets.QLineEdit(self)
        self.startupDirEntry.setText(self.father.defaultStartupDir)
        startupgrid.addWidget(self.startupDirEntry, 0, 0)
        self.startupDirButton = QtWidgets.QPushButton("Browse", self)
        self.startupDirButton.clicked.connect(self.browseStartup)
        startupgrid.addWidget(self.startupDirButton, 0, 1)
        # grid2 definitions
        grid2.addWidget(QtWidgets.QLabel("Linewidth:"), 1, 0)
        self.lwSpinBox = wc.SsnakeDoubleSpinBox()
        self.lwSpinBox.setSingleStep(0.1)
        self.lwSpinBox.setValue(self.father.defaultLinewidth)
        grid2.addWidget(self.lwSpinBox, 1, 1)
        self.color = self.father.defaultColor
        lineColorButton = QtWidgets.QPushButton("Line colour")
        lineColorButton.clicked.connect(self.setColor)
        grid2.addWidget(lineColorButton, 2, 0)
        grid2.addWidget(QtWidgets.QLabel("Colour range:"), 3, 0)
        self.crEntry = QtWidgets.QComboBox(self)
        self.crEntry.addItems(views.COLORRANGELIST)
        self.crEntry.setCurrentIndex(views.COLORRANGELIST.index(self.father.defaultColorRange))
        grid2.addWidget(self.crEntry, 3, 1)
        self.xgridCheck = QtWidgets.QCheckBox("x-grid")
        self.xgridCheck.setChecked(self.father.defaultGrids[0])
        grid2.addWidget(self.xgridCheck, 4, 0, 1, 2)
        self.ygridCheck = QtWidgets.QCheckBox("y-grid")
        self.ygridCheck.setChecked(self.father.defaultGrids[1])
        grid2.addWidget(self.ygridCheck, 5, 0, 1, 2)
        grid2.addWidget(QtWidgets.QLabel("Min X Ticks:"), 6, 0)
        self.xTicksSpinBox = wc.SsnakeSpinBox()
        self.xTicksSpinBox.setValue(self.father.defaultMinXTicks)
        grid2.addWidget(self.xTicksSpinBox, 6, 1)
        grid2.addWidget(QtWidgets.QLabel("Min Y Ticks:"), 7, 0)
        self.yTicksSpinBox = wc.SsnakeSpinBox()
        self.yTicksSpinBox.setValue(self.father.defaultMinYTicks)
        grid2.addWidget(self.yTicksSpinBox, 7, 1)
        grid2.addWidget(QtWidgets.QLabel("Units:"), 8, 0)
        self.unitGroup = QtWidgets.QButtonGroup()
        button = QtWidgets.QRadioButton("s/Hz")
        self.unitGroup.addButton(button, 0)
        grid2.addWidget(button, 9, 1)
        button = QtWidgets.QRadioButton("ms/kHz")
        self.unitGroup.addButton(button, 1)
        grid2.addWidget(button, 10, 1)
        button = QtWidgets.QRadioButton(u"μs/MHz")
        self.unitGroup.addButton(button, 2)
        grid2.addWidget(button, 11, 1)
        self.unitGroup.button(self.father.defaultUnits).setChecked(True)
        self.ppmCheck = QtWidgets.QCheckBox("ppm")
        self.ppmCheck.setChecked(self.father.defaultPPM)
        grid2.addWidget(self.ppmCheck, 12, 1)
        self.zeroScrollCheck = QtWidgets.QCheckBox("Scroll y-axis from zero")
        self.zeroScrollCheck.setChecked(self.father.defaultZeroScroll)
        grid2.addWidget(self.zeroScrollCheck, 13, 0, 1, 2)
        grid2.addWidget(QtWidgets.QLabel("Zoom step:"), 14, 0)
        self.ZoomStepSpinBox = wc.SsnakeDoubleSpinBox()
        self.ZoomStepSpinBox.setSingleStep(0.1)
        self.ZoomStepSpinBox.setValue(self.father.defaultZoomStep)
        grid2.addWidget(self.ZoomStepSpinBox, 14, 1)
        self.showTitleCheck = QtWidgets.QCheckBox("Show title in plot")
        self.showTitleCheck.setChecked(self.father.defaultShowTitle)
        grid2.addWidget(self.showTitleCheck, 15, 0, 1, 2)
        grid2.addWidget(QtWidgets.QLabel("Significant digits:"), 16, 0)
        self.precisSpinBox = wc.SsnakeSpinBox()
        self.precisSpinBox.setValue(self.father.defaultPrecis)
        grid2.addWidget(self.precisSpinBox, 16, 1)
        # grid3 definitions
        grid3.addWidget(QtWidgets.QLabel("Colourmap:"), 0, 0)
        self.cmEntry = QtWidgets.QComboBox(self)
        self.cmEntry.addItems(views.COLORMAPLIST)
        self.cmEntry.setCurrentIndex(views.COLORMAPLIST.index(self.father.defaultColorMap))
        grid3.addWidget(self.cmEntry, 0, 1)
        self.constColorCheck = QtWidgets.QCheckBox("Constant colours")
        self.constColorCheck.setChecked(self.father.defaultContourConst)
        grid3.addWidget(self.constColorCheck, 1, 0)
        self.posColor = self.father.defaultPosColor
        posColorButton = QtWidgets.QPushButton("Positive colour")
        posColorButton.clicked.connect(self.setPosColor)
        grid3.addWidget(posColorButton, 2, 0)
        self.negColor = self.father.defaultNegColor
        negColorButton = QtWidgets.QPushButton("Negative colour")
        negColorButton.clicked.connect(self.setNegColor)
        grid3.addWidget(negColorButton, 3, 0)
        grid3.addWidget(QtWidgets.QLabel("Width ratio:"), 4, 0)
        self.WRSpinBox = wc.SsnakeDoubleSpinBox()
        self.WRSpinBox.setSingleStep(0.1)
        self.WRSpinBox.setValue(self.father.defaultWidthRatio)
        grid3.addWidget(self.WRSpinBox, 4, 1)
        grid3.addWidget(QtWidgets.QLabel("Height ratio:"), 5, 0)
        self.HRSpinBox = wc.SsnakeDoubleSpinBox()
        self.HRSpinBox.setSingleStep(0.1)
        self.HRSpinBox.setValue(self.father.defaultHeightRatio)
        grid3.addWidget(self.HRSpinBox, 5, 1)
        # 2D Colour defs
        grid4.addWidget(QtWidgets.QLabel("Colourmap:"), 0, 0)
        self.cmEntry2D = QtWidgets.QComboBox(self)
        self.cmEntry2D.addItems(views.COLORMAPLIST)
        self.cmEntry2D.setCurrentIndex(views.COLORMAPLIST.index(self.father.defaultPColorMap))
        grid4.addWidget(self.cmEntry2D, 0, 1)
        # Phasing Options (if 2nd order should be available)
        self.secondOrderPhaseCheckBox = QtWidgets.QCheckBox("Always show 2nd order phase correction")
        self.secondOrderPhaseCheckBox.setChecked(self.father.defaultSecondOrderPhaseDialog)
        grid5.addWidget(self.secondOrderPhaseCheckBox, 0, 1)
        # Others
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(tabWidget, 0, 0, 1, 4)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        okButton = QtWidgets.QPushButton("&Store")
        okButton.clicked.connect(self.applyAndClose)
        resetButton = QtWidgets.QPushButton("&Reset")
        resetButton.clicked.connect(self.reset)
        box = QtWidgets.QDialogButtonBox()
        box.addButton(cancelButton, QtWidgets.QDialogButtonBox.RejectRole)
        box.addButton(okButton, QtWidgets.QDialogButtonBox.ActionRole)
        box.addButton(resetButton, QtWidgets.QDialogButtonBox.ActionRole)
        layout.addWidget(box, 1,0)
        layout.setColumnStretch(3, 1)
        self.show()

    def browseStartup(self, *args):
        newDir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select Directory', self.father.lastLocation, QtWidgets.QFileDialog.ShowDirsOnly)
        if newDir:
            self.startupDirEntry.setText(newDir)

    def setColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtGui.QColor(self.color))
        if tmp.isValid():
            self.color = tmp.name()

    def setPosColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtGui.QColor(self.posColor))
        if tmp.isValid():
            self.posColor = tmp.name()

    def setNegColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtGui.QColor(self.negColor))
        if tmp.isValid():
            self.negColor = tmp.name()

    def applyAndClose(self, *args):
        self.father.defaultUnits = self.unitGroup.checkedId()
        self.father.defaultPPM = self.ppmCheck.isChecked()
        self.father.defaultWidth = self.widthSpinBox.value()
        self.father.defaultHeight = self.heightSpinBox.value()
        self.father.defaultMaximized = self.maximizedCheck.isChecked()
        self.father.defaultAskName = self.askNameCheck.isChecked()
        self.father.defaultToolBar = self.toolbarCheck.isChecked()
        self.father.defaultTooltips = self.tooltipCheck.isChecked()
        self.father.defaultToolbarActionList = self.currentToolbar
        self.father.defaultStartupBool = self.startupgroupbox.isChecked()
        self.father.defaultStartupDir = self.startupDirEntry.text()
        self.father.defaultLinewidth = self.lwSpinBox.value()
        self.father.defaultMinXTicks = self.xTicksSpinBox.value()
        self.father.defaultMinYTicks = self.yTicksSpinBox.value()
        self.father.defaultColor = self.color
        self.father.defaultColorRange = self.crEntry.currentText()
        self.father.defaultGrids[0] = self.xgridCheck.isChecked()
        self.father.defaultGrids[1] = self.ygridCheck.isChecked()
        self.father.defaultZeroScroll = self.zeroScrollCheck.isChecked()
        self.father.defaultShowTitle = self.showTitleCheck.isChecked()
        self.father.defaultPrecis = self.precisSpinBox.value()
        self.father.defaultZoomStep = self.ZoomStepSpinBox.value()
        self.father.defaultColorMap = self.cmEntry.currentText()
        self.father.defaultContourConst = self.constColorCheck.isChecked()
        self.father.defaultPosColor = self.posColor
        self.father.defaultNegColor = self.negColor
        self.father.defaultWidthRatio = self.WRSpinBox.value()
        self.father.defaultHeightRatio = self.HRSpinBox.value()
        self.father.defaultPColorMap = self.cmEntry2D.currentText()
        self.father.defaultSecondOrderPhaseDialog = self.secondOrderPhaseCheckBox.isChecked()
        self.father.saveDefaults()
        self.closeEvent()

    def reset(self, *args):
        self.father.resetDefaults()
        self.father.saveDefaults()
        self.closeEvent()

    def closeEvent(self, *args):
        self.deleteLater()

##############################################################################


class ToolbarWindow(wc.ToolWindow):

    NAME = "Change Toolbar"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent):
        super(ToolbarWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Actions:"), 0, 0)
        self.grid.addWidget(wc.QLabel("Toolbar Actions:"), 0, 1)
        self.listB = DestListWidget(self)
        for i in self.father.father.defaultToolbarActionList:
            QtWidgets.QListWidgetItem(i, self.listB).setToolTip(i)
        self.listA = OrigListWidget(self, self.listB)
        for i in self.father.father.allActionsList:
            QtWidgets.QListWidgetItem(i[0], self.listA).setToolTip(i[0])
        self.grid.addWidget(self.listA, 1, 0, 2, 1)
        self.grid.addWidget(self.listB, 1, 2, 2, 1)
        self.rightPush = QtWidgets.QPushButton(u"\u2192", self)
        self.leftPush = QtWidgets.QPushButton(u"\u2190", self)
        self.downPush = QtWidgets.QPushButton(u"\u2193", self)
        self.upPush = QtWidgets.QPushButton(u"\u2191", self)
        self.rightPush.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Expanding)
        self.leftPush.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Expanding)
        self.downPush.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Expanding)
        self.upPush.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Expanding)
        self.grid.addWidget(self.rightPush, 1, 1)
        self.grid.addWidget(self.leftPush, 2, 1)
        self.grid.addWidget(self.upPush, 1, 3)
        self.grid.addWidget(self.downPush, 2, 3)
        self.leftPush.clicked.connect(self.right2left)
        self.rightPush.clicked.connect(self.left2right)
        self.upPush.clicked.connect(self.moveUp)
        self.downPush.clicked.connect(self.moveDown)
        self.resize(650, 500)

    def right2left(self):
        self.listB.deleteSelected()

    def left2right(self):
        for item in self.listA.selectedItems():
            self.listB.addItem(item.text())

    def moveUp(self):
        self.listB.moveSelection('up')

    def moveDown(self):
        self.listB.moveSelection('down')

    def applyAndClose(self, *args):
        items = []
        for index in range(self.listB.count()):
            items.append(self.listB.item(index).text())
        self.father.currentToolbar = items
        self.closeEvent()

    def closeEvent(self, *args):
        self.deleteLater()

##############################################################################


class aboutWindow(wc.ToolWindow):

    NAME = "About ssNake"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent):
        super(aboutWindow, self).__init__(parent)
        self.cancelButton.hide()
        self.logo = QtWidgets.QLabel(self)
        self.logo.setPixmap(QtGui.QPixmap(os.path.dirname(os.path.realpath(__file__)) + "/Icons/logo.gif"))
        self.tabs = QtWidgets.QTabWidget(self)
        self.text = QtWidgets.QTextBrowser(self)
        self.text.setOpenExternalLinks(True)
        self.license = QtWidgets.QTextBrowser(self)
        self.license.setOpenExternalLinks(True)
        licenseText = ''
        with open(os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'licenseHtml.txt') as f:
            licenseText = f.read()
        self.license.setHtml(licenseText)
        pythonVersion = sys.version
        pythonVersion = pythonVersion[:pythonVersion.index(' ')]
        from scipy import __version__ as scipyVersion
        self.text.setText('<p><b>ssNake ' + VERSION + '</b></p>' +
                          '<p>Copyright (&copy;) 2016&ndash;2024 Bas van Meerten & Wouter Franssen</p>' + '<p>Email: <a href="mailto:ssnake@science.ru.nl" >ssnake@science.ru.nl</a></p>' +
                          '<p>Publication: <a href="https://doi.org/10.1016/j.jmr.2019.02.006" >https://doi.org/10.1016/j.jmr.2019.02.006</a></p>' +
                          '<b>Library versions</b>:<br>Python ' + pythonVersion + '<br>numpy ' + np.__version__ +
                          '<br>SciPy ' + scipyVersion +
                          '<br>matplotlib ' + matplotlib.__version__ +
                          '<br>PyQt ' + QtCore.PYQT_VERSION_STR +
                          '<br>Qt ' + QtCore.QT_VERSION_STR)
        self.thanks = QtWidgets.QTextEdit(self)
        self.thanks.setReadOnly(True)
        self.thanks.setHtml('<p><b>The ssNake team wishes to thank:</b></p>Prof. Arno Kentgens<br>Koen Tijssen<br>' +
                            'Ole Brauckmann<br>Merijn Blaakmeer<br>Vincent Breukels<br>Ernst van Eck<br>Fleur van Zelst<br>' +
                            'Sander Lambregts<br>Dr. Andreas Brinkmann<br>Julien Trébosc<br>Henrik Bradtmüller')
        self.tabs.addTab(self.text, 'Version')
        self.tabs.addTab(self.thanks, 'Thanks')
        self.tabs.addTab(self.license, 'License')
        self.grid.addWidget(self.logo, 0, 0, 1, 3, QtCore.Qt.AlignHCenter)
        self.grid.addWidget(self.tabs, 1, 0, 1, 3)
        self.resize(550, 700)

    def closeEvent(self, *args):
        self.deleteLater()

##############################################################################


class shiftConversionWindow(wc.ToolWindow):

    NAME = "Chemical Shift Conversions"
    MENUDISABLE = False
    RESIZABLE = True

    def __init__(self, parent):
        super(shiftConversionWindow, self).__init__(parent)
        self.standardGroup = QtWidgets.QGroupBox('Standard Convention:')
        self.standardFrame = QtWidgets.QGridLayout()
        D11label = wc.QLabel(u'δ' + '<sub>11</sub> [ppm]')
        self.standardFrame.addWidget(D11label, 0, 1)
        D22label = wc.QLabel(u'δ' + '<sub>22</sub> [ppm]')
        self.standardFrame.addWidget(D22label, 0, 2)
        D33label = wc.QLabel(u'δ' + '<sub>33</sub> [ppm]')
        self.standardFrame.addWidget(D33label, 0, 3)
        standardGO = QtWidgets.QPushButton("Go")
        standardGO.setMinimumWidth(100)
        self.standardFrame.addWidget(standardGO, 1, 0)
        standardGO.clicked.connect(lambda: self.shiftCalc(0))
        self.D11 = wc.QLineEdit("0")
        self.D11.setMinimumWidth(100)
        self.standardFrame.addWidget(self.D11, 1, 1)
        self.D22 = wc.QLineEdit("0")
        self.D22.setMinimumWidth(100)
        self.standardFrame.addWidget(self.D22, 1, 2)
        self.D33 = wc.QLineEdit("0")
        self.D33.setMinimumWidth(100)
        self.standardFrame.addWidget(self.D33, 1, 3)
        self.standardGroup.setLayout(self.standardFrame)
        self.grid.addWidget(self.standardGroup, 0, 0, 1, 3)
        # xyz Convention
        self.xyzGroup = QtWidgets.QGroupBox('xyz Convention:')
        self.xyzFrame = QtWidgets.QGridLayout()
        dxxlabel = wc.QLabel(u'δ' + '<sub>xx</sub> [ppm]')
        self.xyzFrame.addWidget(dxxlabel, 3, 1)
        dyylabel = wc.QLabel(u'δ' + '<sub>yy</sub> [ppm]')
        self.xyzFrame.addWidget(dyylabel, 3, 2)
        dzzlabel = wc.QLabel(u'δ' + '<sub>zz</sub> [ppm]')
        self.xyzFrame.addWidget(dzzlabel, 3, 3)
        xyzGO = QtWidgets.QPushButton("Go")
        xyzGO.setMinimumWidth(100)
        self.xyzFrame.addWidget(xyzGO, 4, 0)
        xyzGO.clicked.connect(lambda: self.shiftCalc(1))
        self.dxx = wc.QLineEdit("0")
        self.dxx.setMinimumWidth(100)
        self.xyzFrame.addWidget(self.dxx, 4, 1)
        self.dyy = wc.QLineEdit("0")
        self.dyy.setMinimumWidth(100)
        self.xyzFrame.addWidget(self.dyy, 4, 2)
        self.dzz = wc.QLineEdit("0")
        self.dzz.setMinimumWidth(100)
        self.xyzFrame.addWidget(self.dzz, 4, 3)
        self.xyzGroup.setLayout(self.xyzFrame)
        self.grid.addWidget(self.xyzGroup, 1, 0, 1, 3)
        # Haeberlen Convention
        self.haebGroup = QtWidgets.QGroupBox('Haeberlen Convention')
        self.haebFrame = QtWidgets.QGridLayout()
        disolabel = wc.QLabel(u'δ' + '<sub>iso</sub> [ppm]')
        self.haebFrame.addWidget(disolabel, 6, 1)
        danisolabel = wc.QLabel(u'δ' + '<sub>aniso</sub> [ppm]')
        self.haebFrame.addWidget(danisolabel, 6, 2)
        etalabel = wc.QLabel(u'η')
        self.haebFrame.addWidget(etalabel, 6, 3)
        haeberGO = QtWidgets.QPushButton("Go")
        haeberGO.setMinimumWidth(100)
        self.haebFrame.addWidget(haeberGO, 7, 0)
        haeberGO.clicked.connect(lambda: self.shiftCalc(2))
        self.diso = wc.QLineEdit("0")
        self.diso.setMinimumWidth(100)
        self.haebFrame.addWidget(self.diso, 7, 1)
        self.daniso = wc.QLineEdit("0")
        self.daniso.setMinimumWidth(100)
        self.haebFrame.addWidget(self.daniso, 7, 2)
        self.eta = wc.QLineEdit("0")
        self.eta.setMinimumWidth(100)
        self.haebFrame.addWidget(self.eta, 7, 3)
        self.haebGroup.setLayout(self.haebFrame)
        self.grid.addWidget(self.haebGroup, 2, 0, 1, 3)
        # Hertzfeld berger
        self.hbGroup = QtWidgets.QGroupBox('Hertzfeld-Berger Convention')
        self.hbFrame = QtWidgets.QGridLayout()
        hbdisolabel = wc.QLabel(u'δ' + '<sub>iso</sub> [ppm]')
        self.hbFrame.addWidget(hbdisolabel, 9, 1)
        omegalabel = wc.QLabel(u'Ω [ppm]')
        self.hbFrame.addWidget(omegalabel, 9, 2)
        skewlabel = wc.QLabel(u'κ')
        self.hbFrame.addWidget(skewlabel, 9, 3)
        hbGO = QtWidgets.QPushButton("Go")
        hbGO.setMinimumWidth(100)
        self.hbFrame.addWidget(hbGO, 10, 0)
        hbGO.clicked.connect(lambda: self.shiftCalc(3))
        self.hbdiso = wc.QLineEdit("0")
        self.hbdiso.setMinimumWidth(100)
        self.hbFrame.addWidget(self.hbdiso, 10, 1)
        self.hbdaniso = wc.QLineEdit("0")
        self.hbdaniso.setMinimumWidth(100)
        self.hbFrame.addWidget(self.hbdaniso, 10, 2)
        self.hbskew = wc.QLineEdit("0")
        self.hbskew.setMinimumWidth(100)
        self.hbFrame.addWidget(self.hbskew, 10, 3)
        self.hbGroup.setLayout(self.hbFrame)
        self.grid.addWidget(self.hbGroup, 3, 0, 1, 3)
        # Reset
        self.cancelButton.setText("Close")
        self.cancelButton.clicked.disconnect()
        self.cancelButton.clicked.connect(self.closeEvent)
        self.okButton.setText("Reset")
        self.okButton.clicked.disconnect()
        self.okButton.clicked.connect(self.valueReset)

    def shiftCalc(self, Type):
        if Type == 0:  # If from standard
            try:
                delta11 = float(safeEval(self.D11.text(), Type='FI'))
                delta22 = float(safeEval(self.D22.text(), Type='FI'))
                delta33 = float(safeEval(self.D33.text(), Type='FI'))
                Values = [delta11, delta22, delta33]
            except Exception:
                raise SsnakeException("Shift Conversion: Invalid input in Standard Convention")
        if Type == 1:  # If from xyz
            try:
                delta11 = float(safeEval(self.dxx.text(), Type='FI'))  # Treat xyz as 123, as it reorders them anyway
                delta22 = float(safeEval(self.dyy.text(), Type='FI'))
                delta33 = float(safeEval(self.dzz.text(), Type='FI'))
                Values = [delta11, delta22, delta33]
            except Exception:
                raise SsnakeException("Shift Conversion: Invalid input in xyz Convention")
        if Type == 2:  # From haeberlen
            try:
                eta = float(safeEval(self.eta.text(), Type='FI'))
                delta = float(safeEval(self.daniso.text(), Type='FI'))
                iso = float(safeEval(self.diso.text(), Type='FI'))
                Values = [iso, delta, eta]
            except Exception:
                raise SsnakeException("Shift Conversion: Invalid input in Haeberlen Convention")
        if Type == 3:  # From Hertzfeld-Berger
            try:
                iso = float(safeEval(self.hbdiso.text(), Type='FI'))
                span = float(safeEval(self.hbdaniso.text(), Type='FI'))
                skew = float(safeEval(self.hbskew.text(), Type='FI'))
                Values = [iso, span, skew]
            except Exception:
                raise SsnakeException("Shift Conversion: Invalid input in Hertzfeld-Berger Convention")
        Results = func.shiftConversion(Values, Type)  # Do the actual conversion
        # Standard convention
        self.D11.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[0][0])
        self.D22.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[0][1])
        self.D33.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[0][2])
        # Convert to haeberlen convention and xxyyzz
        self.dxx.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[1][0])
        self.dyy.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[1][1])
        self.dzz.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[1][2])
        # Haeberlen def
        self.diso.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[2][0])
        self.daniso.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[2][1])
        try:  # If a number
            self.eta.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[2][2])
        except Exception:
            self.eta.setText('ND')
        # Convert to Herzfeld-Berger Convention
        self.hbdiso.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[3][0])
        self.hbdaniso.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[3][1])
        try:
            self.hbskew.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Results[3][2])
        except Exception:
            self.hbskew.setText('ND')

    def valueReset(self):  # Resets all the boxes to 0
        self.D11.setText('0')
        self.D22.setText('0')
        self.D33.setText('0')
        self.dxx.setText('0')
        self.dyy.setText('0')
        self.dzz.setText('0')
        self.eta.setText('0')
        self.diso.setText('0')
        self.daniso.setText('0')
        self.hbskew.setText('0')
        self.hbdiso.setText('0')
        self.hbdaniso.setText('0')

    def closeEvent(self, *args):
        self.deleteLater()

##############################################################################

class quadMqStToolWindow(wc.ToolWindow):
    Ioptions = ['1', '3/2', '2', '5/2', '3', '7/2', '4', '9/2', '5', '6', '7']
    Ivalues = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0]
    NAME = "Factor calculation for MQMAS/STMAS experiments"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent):
        super(quadMqStToolWindow, self).__init__(parent)

class quadConversionWindow(wc.ToolWindow):

    Ioptions = ['1', '3/2', '2', '5/2', '3', '7/2', '4', '9/2', '5', '6', '7']
    Ivalues = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0]
    NAME = "Quadrupolar Coupling Conversions"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent):
        super(quadConversionWindow, self).__init__(parent)
        self.comGroup = QtWidgets.QGroupBox("Common Parameters:")
        self.comFrame = QtWidgets.QGridLayout()
        nucLabel = wc.QLabel("Nucleus:")
        self.comFrame.addWidget(nucLabel, 0, 0)
        qindex = [x for (x, val) in enumerate(ISOTOPES['q']) if val is not None and ISOTOPES['spin'][x] != 0.5]
        self.names = ['User'] +  [val for (x, val) in enumerate(ISOTOPES['formatName']) if x in qindex]
        self.I = [0.0] + [val for (x, val) in enumerate(ISOTOPES['spin']) if x in qindex]
        self.Qvalues = [0.0] + [val for (x, val) in enumerate(ISOTOPES['q']) if x in qindex]
        self.nucDrop = QtWidgets.QComboBox()
        self.nucDrop.addItems(self.names)
        self.nucDrop.currentIndexChanged.connect(self.setNuc)
        self.comFrame.addWidget(self.nucDrop, 1, 0)
        Itext = wc.QLabel("I:")
        self.comFrame.addWidget(Itext, 0, 1)
        self.IEntry = QtWidgets.QComboBox()
        self.IEntry.addItems(self.Ioptions)
        self.IEntry.setCurrentIndex(0)
        self.IEntry.activated.connect(self.userChange)
        self.comFrame.addWidget(self.IEntry, 1, 1)
        etalabel = wc.QLabel(u'η:')
        self.comFrame.addWidget(etalabel, 0, 3)
        self.Eta = wc.QLineEdit("0")
        self.Eta.setMinimumWidth(100)
        self.comFrame.addWidget(self.Eta, 1, 3)
        momentlabel = wc.QLabel('Q [fm<sup>2</sup>]:')
        self.comFrame.addWidget(momentlabel, 0, 2)
        self.Moment = wc.QLineEdit("ND")
        self.Moment.setMinimumWidth(100)
        self.Moment.textEdited.connect(self.userChange)
        self.comFrame.addWidget(self.Moment, 1, 2)
        self.comGroup.setLayout(self.comFrame)
        self.grid.addWidget(self.comGroup, 1, 0, 1, 3)
        self.CqGroup = QtWidgets.QGroupBox("C_Q Convention:")
        self.CqFrame = QtWidgets.QGridLayout()
        Cqlabel = wc.QLabel(u'C' + u'<sub>Q</sub>/2π [MHz:]')
        self.CqFrame.addWidget(Cqlabel, 3, 1)
        CqGO = QtWidgets.QPushButton("Go")
        self.CqFrame.addWidget(CqGO, 4, 0)
        CqGO.clicked.connect(lambda: self.quadCalc(0))
        self.Cq = wc.QLineEdit("0")
        self.Cq.setMinimumWidth(100)
        self.CqFrame.addWidget(self.Cq, 4, 1)
        self.CqGroup.setLayout(self.CqFrame)
        self.grid.addWidget(self.CqGroup, 3, 0, 1, 2)
        self.WqGroup = QtWidgets.QGroupBox(u"ω_Q Convention:")
        self.WqFrame = QtWidgets.QGridLayout()
        Wqlabel = wc.QLabel(u'ω' + u'<sub>Q</sub>/2π [MHz]:')
        self.WqFrame.addWidget(Wqlabel, 6, 1)
        WqGO = QtWidgets.QPushButton("Go")
        self.WqFrame.addWidget(WqGO, 7, 0)
        WqGO.clicked.connect(lambda: self.quadCalc(1))
        self.Wq = wc.QLineEdit("0")
        self.Wq.setMinimumWidth(100)
        self.WqFrame.addWidget(self.Wq, 7, 1)
        self.WqGroup.setLayout(self.WqFrame)
        self.grid.addWidget(self.WqGroup, 7, 0, 1, 2)
        self.fieldGroup = QtWidgets.QGroupBox('Field Gradients:')
        self.fieldFrame = QtWidgets.QGridLayout()
        # Vxx and Vyy labels interchanged to follow the Quad+CSa fitting deifnitions
        # WF: 2021-01
        Vxxlabel = wc.QLabel('V<sub>yy</sub> [V/m<sup>2</sup>]:')
        self.fieldFrame.addWidget(Vxxlabel, 9, 1)
        VGO = QtWidgets.QPushButton("Go")
        self.fieldFrame.addWidget(VGO, 10, 0)
        VGO.clicked.connect(lambda: self.quadCalc(2))
        self.Vxx = wc.QLineEdit("ND")
        self.Vxx.setMinimumWidth(100)
        self.fieldFrame.addWidget(self.Vxx, 10, 1)
        Vyylabel = wc.QLabel('V<sub>xx</sub> [V/m<sup>2</sup>]:')
        self.fieldFrame.addWidget(Vyylabel, 9, 2)
        self.Vyy = wc.QLineEdit("ND")
        self.Vyy.setMinimumWidth(100)
        self.fieldFrame.addWidget(self.Vyy, 10, 2)
        Vzzlabel = wc.QLabel('V<sub>zz</sub> [V/m<sup>2</sup>]:')
        self.fieldFrame.addWidget(Vzzlabel, 9, 3)
        self.Vzz = wc.QLineEdit("ND")
        self.Vzz.setMinimumWidth(100)
        self.fieldFrame.addWidget(self.Vzz, 10, 3)
        self.fieldGroup.setLayout(self.fieldFrame)
        self.grid.addWidget(self.fieldGroup, 8, 0, 1, 4)
        # Reset
        self.cancelButton.setText("Close")
        self.cancelButton.clicked.disconnect()
        self.cancelButton.clicked.connect(self.closeEvent)
        self.okButton.setText("Reset")
        self.okButton.clicked.disconnect()
        self.okButton.clicked.connect(self.valueReset)

    def setNuc(self, index):
        I = self.I[index]
        if int(I * 2) > 0:
            self.IEntry.setCurrentIndex(int(I * 2 - 2))
            self.Moment.setText(str(self.Qvalues[index]))

    def userChange(self, index=None):
        self.nucDrop.setCurrentIndex(0)

    def quadCalc(self, Type):
        I = self.Ivalues[self.IEntry.currentIndex()]
        if Type == 0:  # Cq as input
            # Czz is equal to Cq, via same definition (scale) Cxx and Cyy can be found
            try:
                Cq = float(safeEval(self.Cq.text(), Type='FI'))
                Eta = float(safeEval(self.Eta.text(), Type='FI'))
                Values = [Cq, Eta]
            except Exception:
                raise SsnakeException("Quad Conversion: Invalid input in Cq definition")
        if Type == 1:
            try:
                Wq = float(safeEval(self.Wq.text(), Type='FI'))
                Eta = float(safeEval(self.Eta.text(), Type='FI'))
                Values = [Wq, Eta]
            except Exception:
                raise SsnakeException("Quad Conversion: Invalid input in Wq definition")
        if Type == 2:
            try:
                Vxx = float(safeEval(self.Vxx.text(), Type='FI'))
                Vyy = float(safeEval(self.Vyy.text(), Type='FI'))
                Vzz = float(safeEval(self.Vzz.text(), Type='FI'))
                Values = [Vxx, Vyy, Vzz]
            except Exception:
                raise SsnakeException("Quad Conversion: Invalid input in field gradients")
        try:
            Q = float(safeEval(self.Moment.text(), Type='FI')) * 1e-30  # get moment and convert from fm^2
        except Exception:
            if Type in (0, 1):
                Q = None
            else:
                raise SsnakeException("Quad Conversion: Invalid input in quadrupole moment Q")
        #Do conversion
        Result = func.quadConversion(Values, I, Type, Q)
        if Result[0][1] is None:
            self.Eta.setText('ND')
        else:
            self.Eta.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Result[0][1])
        self.Cq.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Result[0][0])
        self.Wq.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Result[1][0])
        if Result[2][0] is None:
            self.Moment.setText('ND')
            self.Vxx.setText('ND')
            self.Vyy.setText('ND')
            self.Vzz.setText('ND')
        else:
            self.Vxx.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Result[2][0])
            self.Vyy.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Result[2][1])
            self.Vzz.setText(('%#.'+str(self.father.defaultPrecis)+'g') % Result[2][2])

    def valueReset(self):  # Resets all the boxes to 0
        self.Cq.setText('0')
        self.Eta.setText('0')
        self.Wq.setText('0')
        self.Moment.setText('ND')
        self.Vxx.setText('ND')
        self.Vyy.setText('ND')
        self.Vzz.setText('ND')

    def closeEvent(self, *args):
        self.deleteLater()

##############################################################################

class dipolarDistanceWindow(wc.ToolWindow):

    NAME = "Dipolar Distance Calculation"
    RESIZABLE = True
    MENUDISABLE = False
    Ioptions = ['1/2','1', '3/2', '2', '5/2', '3', '7/2', '4', '9/2', '5', '6', '7']
    Ivalues = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0]

    def __init__(self, parent):
        super(dipolarDistanceWindow, self).__init__(parent)
        # Nuclei group
        self.comGroup = QtWidgets.QGroupBox("Gyromagnetic Ratios and Spin Quantum Numbers:")
        self.comFrame = QtWidgets.QGridLayout()
        gamma1label = wc.QLabel(u'γ<sub>1</sub> [10<sup>7</sup> rad/s/T]:')
        self.comFrame.addWidget(gamma1label, 0, 0)

        gammaindex = [x for (x, val) in enumerate(ISOTOPES['gamma']) if val is not None]
        self.gammaValues = [0.0] + [val for (x, val) in enumerate(ISOTOPES['gamma']) if x in gammaindex]
        self.spinValues = [0.0] + [val for (x, val) in enumerate(ISOTOPES['spin']) if x in gammaindex]
        self.names = ['User'] + [val for (x, val) in enumerate(ISOTOPES['formatName']) if x in gammaindex]
        self.gamma1Drop = QtWidgets.QComboBox()
        self.gamma1Drop.addItems(self.names)
        self.gamma1Drop.currentIndexChanged.connect(self.setGamma1)
        self.comFrame.addWidget(self.gamma1Drop, 1, 0)
        self.gamma2Drop = QtWidgets.QComboBox()
        self.gamma2Drop.addItems(self.names)
        self.gamma2Drop.currentIndexChanged.connect(self.setGamma2)
        self.comFrame.addWidget(self.gamma2Drop, 1, 1)
        self.gamma1 = wc.QLineEdit("0.0")
        self.gamma1.setMinimumWidth(100)
        self.gamma1.textEdited.connect(self.gamma1Changed)
        self.comFrame.addWidget(self.gamma1, 2, 0)
        spin1label = wc.QLabel(u'<i>I</i><sub>1</sub>:')
        self.comFrame.addWidget(spin1label, 3, 0)

        self.spin1 = QtWidgets.QComboBox()
        self.spin1.addItems(self.Ioptions)
        self.spin1.currentIndexChanged.connect(self.spin1Changed)
        self.comFrame.addWidget(self.spin1, 4, 0)
        gamma2label = wc.QLabel(u'γ<sub>2</sub> [10<sup>7</sup> rad/s/T]:')
        self.comFrame.addWidget(gamma2label, 0, 1)
        self.gamma2 = wc.QLineEdit("0.0")
        self.gamma2.setMinimumWidth(100)
        self.gamma2.textEdited.connect(self.gamma2Changed)
        self.comFrame.addWidget(self.gamma2, 2, 1)
        spin2label = wc.QLabel(u'<i>I</i><sub>2</sub>:')
        self.comFrame.addWidget(spin2label, 3, 1)
        self.spin2 = QtWidgets.QComboBox()
        self.spin2.addItems(self.Ioptions)
        self.spin2.currentIndexChanged.connect(self.spin2Changed)
        self.comFrame.addWidget(self.spin2, 4, 1)
        self.comGroup.setLayout(self.comFrame)
        #addWidget - fromRow, fromColumn, rowSpan, columnSpan
        self.grid.addWidget(self.comGroup, 0, 0, 5, 2)
        # Distance group
        self.distanceGroup = QtWidgets.QGroupBox("Distance:")
        self.distanceFrame = QtWidgets.QGridLayout()
        distancelabel = wc.QLabel(u'r [Å]')
        self.distanceFrame.addWidget(distancelabel, 0, 1)
        distanceGO = QtWidgets.QPushButton("Go")
        self.distanceFrame.addWidget(distanceGO, 1, 0)
        distanceGO.clicked.connect(lambda: self.Calc(0))
        self.distance = wc.QLineEdit("0")
        self.distance.setMinimumWidth(100)
        self.distanceFrame.addWidget(self.distance, 1, 1)
        self.distanceGroup.setLayout(self.distanceFrame)
        self.grid.addWidget(self.distanceGroup, 5, 0, 2, 2)
        # Dipolar coupling group
        self.dipolarGroup = QtWidgets.QGroupBox("Dipolar Coupling:")
        self.dipolarFrame = QtWidgets.QGridLayout()
        dipolarlabel = wc.QLabel(u'D [kHz]')
        self.dipolarFrame.addWidget(dipolarlabel, 0, 1)
        dipolarGO = QtWidgets.QPushButton("Go")
        self.dipolarFrame.addWidget(dipolarGO, 1, 0)
        dipolarGO.clicked.connect(lambda: self.Calc(1))
        self.dipolar = wc.QLineEdit("0")
        self.dipolar.setMinimumWidth(100)
        self.dipolarFrame.addWidget(self.dipolar, 1, 1)
        self.grid.addWidget(self.dipolarGroup, 7, 0, 2, 2)
        self.dipolarGroup.setLayout(self.dipolarFrame)
        # Second-moment group
        self.M2Group = QtWidgets.QGroupBox(u"Second Moments [10⁶ rad²/s²]:")
        self.M2Frame = QtWidgets.QGridLayout()
        
        self.M2hetGO = QtWidgets.QPushButton("Go")
        self.M2Frame.addWidget(self.M2hetGO, 1, 0)
        self.M2hetGO.clicked.connect(lambda: self.Calc(2))
        self.M2label = wc.QLabel(u'Heteronuclear')
        self.M2Frame.addWidget(self.M2label, 0, 0)
        self.M2 = wc.QLineEdit("0")
        self.M2.setMinimumWidth(100)
        self.M2Frame.addWidget(self.M2, 1, 1)
        
        
        self.M2SEDGO = QtWidgets.QPushButton("Go")
        self.M2Frame.addWidget(self.M2SEDGO, 6, 0)
        self.M2SEDGO.clicked.connect(lambda: self.Calc(3))
        self.M2SEDlabel = wc.QLeftLabel(u'Spin Echo Decay Envelope (C<sub>Q</sub> > 0)')
        self.M2Frame.addWidget(self.M2SEDlabel, 5, 0,1,2)
        self.M2SED = wc.QLineEdit("0")
        self.M2SED.setMinimumWidth(100)
        self.M2Frame.addWidget(self.M2SED, 6, 1)
        
        self.grid.addWidget(self.M2Group, 9, 0, 2, 2)
        self.M2Group.setLayout(self.M2Frame)

        # Reset
        self.cancelButton.setText("Close")
        self.cancelButton.clicked.disconnect()
        self.cancelButton.clicked.connect(self.closeEvent)
        self.okButton.setText("Reset")
        self.okButton.clicked.disconnect()
        self.okButton.clicked.connect(self.valueReset)

        self.checkHomoHetero()


    def checkHomoHetero(self):
        gamma1 = float(safeEval(self.gamma1.text(), Type='FI')) * 1e7
        gamma2 = float(safeEval(self.gamma2.text(), Type='FI')) * 1e7
        I1 = self.Ivalues[self.spin1.currentIndex()]
        I2 = self.Ivalues[self.spin2.currentIndex()]
        if gamma1 == gamma2 and I1 == I2:
            self.M2label.setText(u'M₂ Homonuclear')
            if I1 > 0.5:
                self.M2SEDGO.show()
                self.M2SEDlabel.show()
                self.M2SED.show()
            else:
                self.M2SEDGO.hide()
                self.M2SEDlabel.hide()
                self.M2SED.hide()
        else:
            self.M2label.setText(u'M₂ Heteronuclear')
            self.M2SEDGO.hide()
            self.M2SEDlabel.hide()
            self.M2SED.hide()
    
    def spin1Changed(self):
        self.gamma1Drop.setCurrentIndex(0)
        self.checkHomoHetero()
        
    def spin2Changed(self):
        self.gamma2Drop.setCurrentIndex(0)
        self.checkHomoHetero()
            
    def gamma1Changed(self):
        self.gamma1Drop.setCurrentIndex(0)
        self.checkHomoHetero()

    def gamma2Changed(self):
        self.gamma2Drop.setCurrentIndex(0)
        self.checkHomoHetero()

    def setGamma1(self, index):
        if index != 0:
            self.spin1.setCurrentIndex(self.Ivalues.index(self.spinValues[index]))
            self.gamma1.setText(str(self.gammaValues[index]))
            self.gamma1Drop.setCurrentIndex(index) #Needed to avoid resets by other boxes
            self.checkHomoHetero()

    def setGamma2(self, index):
        if index != 0:
            self.spin2.setCurrentIndex(self.Ivalues.index(self.spinValues[index]))
            self.gamma2.setText(str(self.gammaValues[index]))
            self.gamma2Drop.setCurrentIndex(index) #Needed to avoid resets by other boxes
            self.checkHomoHetero()

    def Calc(self, Type):
        try:
            gamma1 = float(safeEval(self.gamma1.text(), Type='FI')) * 1e7
            gamma2 = float(safeEval(self.gamma2.text(), Type='FI')) * 1e7
            I1 = self.Ivalues[self.spin1.currentIndex()]
            I2 = self.Ivalues[self.spin2.currentIndex()]
            #spin1 = float(safeEval(self.spin1.text(), Type='FI'))
            #spin2 = float(safeEval(self.spin2.text(), Type='FI'))
        except Exception:
            raise SsnakeException("Dipolar Distance: Invalid input in gamma values")
        if Type == 0:  # Distance as input
            try:
                r = abs(float(safeEval(self.distance.text(), Type='FI')))
            except Exception:
                raise SsnakeException("Dipolar Distance: Invalid input in r")
        if Type == 1:  # Dipolar coupling as input
            try:
                D = abs(float(safeEval(self.dipolar.text(), Type='FI')))
            except Exception:
                raise SsnakeException("Dipolar Distance: Invalid input in D")
        if Type == 2:  # Heteronuclear second-moment as input
            try:
                M2 = abs(float(safeEval(self.M2.text(), Type='FI'))) 
            except Exception:
                raise SsnakeException("Dipolar Distance: Invalid input in M2")

        if Type == 3: # Homonuclear second-moment from spin echo decay
            try:
                M2SED = abs(float(safeEval(self.M2SED.text(), Type='FI'))) 
            except Exception:
                raise SsnakeException("Dipolar Distance: Invalid input in M2")
                
        hbar = 1.054573e-34
        def wm(m):
            return (1/2) * np.sqrt(I2 * (I2 + 1) - m * (m + 1))
        
        M2_Factor = (2/(9*(2*I2+1))) * (1 + 4*wm(-1/2)**2 + 4*wm(-1/2)**4 + wm(1/2)**4)

        if Type == 0:
            if r == 0.0:
                D = np.inf
                M2 = np.inf
                M2SED = np.inf
            else:
                D = abs(- 1e-7 * gamma1 * gamma2 * hbar / (r * 10**-10) **3 / (2 * np.pi))
                D /= 1000
                if gamma1 == gamma2:
                    M2SED = abs(M2_Factor * (4/5) * (9/4) * 1e-14 * gamma1**4 * hbar**2 / ((r * 10**-10) **6))
                    # factor 4/5 results after powder averaging
                    M2SED /= 1e6
                    M2 = abs((3/5) * 1e-14 * gamma1**4 * I2 * (I2 + 1)* hbar**2 / ((r * 10**-10) **6))
                    M2 /= 1e6
                    
                else:
                    M2 = abs((4/15) * 1e-14 * gamma1**2 * gamma2**2 * I2 * (I2 + 1) * hbar**2 / ((r * 10**-10) **6))
                    M2 /= 1e6
                    M2SED = 0
                    
        if Type == 1:
            if D == 0.0:
                r = np.inf
                M2 = 0.0
                M2SED = 0.0
            else:
                r = 1 / abs(D * 1000 / gamma1 / gamma2 / hbar / 1e-7 * (2 * np.pi))**(1.0/3)
                r *= 1e10
                if gamma1 == gamma2:
                    M2SED = abs(M2_Factor * (4/5) * (9/4) * 1e-14 * gamma1**4 * hbar**2 / ((r * 10**-10) **6))
                    M2SED /= 1e6
                    M2 = abs((3/5) * 1e-14 * gamma1**4 * I2 * (I2 + 1)* hbar**2 / ((r * 10**-10) **6))
                    M2 /= 1e6
                    #factor 4/5 results after powder averaging
                    
                else:
                    M2 = abs((4/15) * 1e-14 * gamma1**2 * gamma2**2 * I2 * (I2 + 1) * hbar**2 / ((r * 10**-10) **6))
                    M2 /= 1e6
                    M2SED = 0.0
                    
        if Type == 2:
            if M2 == 0:
                r = np.inf
                D = 0.0
                M2 = 0.0
            else:
                if gamma1 == gamma2:
                    r = abs((3/5) * 1e-14 * gamma1**4 * I2 * (I2 + 1) * hbar**2 / (M2 * 10**6))
                    r = r**(1/6)
                    r *= 10**10
                    D = abs(- 1e-7 * gamma1 * gamma2 * hbar / (r * 10**-10) **3 / (2 * np.pi))
                    D /= 1000
                    M2SED = abs(M2_Factor * (4/5) * (9/4) * 1e-14 * gamma1**4 * hbar**2 / ((r * 10**-10) **6))
                    M2SED /= 1e6
                else:
                    r = abs((4/15) * 1e-14 * gamma1**2 * gamma2**2 * I2 * (I2 + 1) * hbar**2 / (M2 * 10**6))
                    r = r**(1/6)
                    r *= 10**10
                    D = abs(- 1e-7 * gamma1 * gamma2 * hbar / (r * 10**-10) **3 / (2 * np.pi))
                    D /= 1000
                    M2SED = 0.0
                
        if Type == 3:
            if M2SED == 0:
                r = np.inf
                D = 0.0
                M2 = 0.0
            else:
                if gamma1 == gamma2:
                    r = abs(M2_Factor * (4/5) * (9/4) * 1e-14 * gamma1**4 * hbar**2 / (M2SED * 10**6))
                    # factors 4/5 and 9/4 result from powder averaging
                    r = r**(1/6)
                    r *= 10**10
                    D = abs(- 1e-7 * gamma1 * gamma2 * hbar / (r * 10**-10) **3 / (2 * np.pi))
                    D /= 1000
                    M2 = abs((3/5) * 1e-14 * gamma1**4 * I2 * (I2 + 1)* hbar**2 / ((r * 10**-10) **6))
                    M2 /= 1e6
                else:
                    raise SsnakeException("Choose two identical nuclei.")
            # else:
                
        self.dipolar.setText('%#.5g' % D)
        self.distance.setText('%#.5g' % r)
        self.M2.setText('%#.5g' % M2)
        self.M2SED.setText('%#.5g' % M2SED)
        # Implement FWHM of dipolar line as input

    def valueReset(self):  # Resets all the boxes to 0
        self.dipolar.setText('0.0')
        self.distance.setText('0.0')
        self.M2het.setText('0.0')
        self.M2hom.setText('0.0')
        self.M2SED.setText('0.0')
        self.gamma1.setText('0.0')
        self.gamma2.setText('0.0')
        self.spin1.setCurrentIndex(0)
        self.spin1.setCurrentIndex(0)
        self.gamma1Drop.setCurrentIndex(0)
        self.gamma2Drop.setCurrentIndex(0)
        self.checkHomoHetero()

    def closeEvent(self, *args):
        self.deleteLater()

##############################################################################



class tempCalWindow(QtWidgets.QWidget):
    #[minTemp, maxTemp, shiftToTemp, tempToShift, Delta0, T0, RefName]
    METHANOL = [178,330,
                lambda Delta: 409.0 - 36.54 * Delta - 21.85 * Delta**2,
                lambda Temp: (36.54 - np.sqrt(36.54**2 - 4 * -21.85 * (409.0 - Temp))) / (2 * -21.85),
                'absShift',None,None,'Ammann et al., JMR, 46, 319 (1982)']
    ETH_GLYCOL = [273,416,
                  lambda Delta: 466.5 - 102.00 * Delta,
                  lambda Temp: (Temp - 466.5) / -102.0,
                  'absShift',None,None,'Ammann et al., JMR, 46, 319 (1982)']
    PBNO3 = [143, 423,
             lambda Delta, Delta0, T0: (Delta0 - Delta) / 0.753 + T0 , 
             lambda T, Delta0, T0: (T0 - T) * 0.753 + Delta0 ,
             'relShift','-3473','293','Bielecki et al., JMR, 116, 215 (1995)']
    KBR = [170, 320,
             lambda Delta, Delta0, T0: (Delta0 - Delta) / 0.0250 + T0 , 
             lambda T, Delta0, T0: (T0 - T) * 0.0250 + Delta0 ,
             'relShift','0','293','Thurber et al., JMR, 196, 84 (2009)']
    DEFINITIONS = [METHANOL, ETH_GLYCOL, PBNO3 ,KBR]
    TEXTLIST = ['1H: Methanol (178 K < T < 330 K)', '1H: Ethylene Glycol (273 K < T < 416 K)',
            '207Pb: Lead Nitrate (143 K < T < 423 K)','79Br: KBr (170 K < T < 320 K)']
    T1_KBr = [20,296,
             lambda Relax: optimize.brentq(lambda T,T1: 0.0145 + 5330/T**2 + 1.42e7/T**4 + 2.48e9/T**6 - T1, 20, 296 ,args=(Relax,)),
             lambda T: 0.0145 + 5330/T**2 + 1.42e7/T**4 + 2.48e9/T**6,
             'Thurber et al., JMR, 196, 84 (2009)']
    T1_CsI = [8,104,
              lambda Relax: optimize.brentq(lambda T,T1: -1.6e-3 + 1.52e3/T**2 + 0.387e6/T**4 + 0.121e9/T**6 - T1, 8, 104 ,args=(Relax,)),
              lambda T: -1.6e-3 + 1.52e3/T**2 + 0.387e6/T**4 + 0.121e9/T**6,
              'Sarkar et al., JMR, 212, 460 (2011)']
    T1_DEFINITIONS = [T1_KBr, T1_CsI]
    T1_TEXTLIST = ['79Br: KBr (20 K < T < 296 K, 9.4 T)', '127I: CsI (8 K < T < 104 K, 9.4 T)']

    def __init__(self, parent):
        super(tempCalWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Temperature Calibration")
        tabWidget = QtWidgets.QTabWidget()
        tab1 = QtWidgets.QWidget()
        tabWidget.addTab(tab1, "Chemical Shift Based")
        grid1 = QtWidgets.QGridLayout()
        tab1.setLayout(grid1)
        # Shift based
        self.typeDrop = QtWidgets.QComboBox()
        self.typeDrop.addItems(self.TEXTLIST)
        grid1.addWidget(self.typeDrop, 0, 0)
        self.typeDrop.currentIndexChanged.connect(self.changeType)
        self.RefGroup = QtWidgets.QGroupBox("Relative to:")
        self.RefFrame = QtWidgets.QGridLayout()
        self.Delta0Label = wc.QLabel(u'δ [ppm]')
        self.RefFrame.addWidget(self.Delta0Label, 0, 0)
        self.Delta0 = wc.QLineEdit("")
        self.Delta0.setMinimumWidth(100)
        self.RefFrame.addWidget(self.Delta0, 1, 0)
        self.T0Label = wc.QLabel(u'T [K]')
        self.RefFrame.addWidget(self.T0Label, 0, 1)
        self.T0 = wc.QLineEdit("")
        self.T0.setMinimumWidth(100)
        self.RefFrame.addWidget(self.T0, 1, 1)
        self.RefGroup.setLayout(self.RefFrame)
        grid1.addWidget(self.RefGroup, 1, 0)
        self.RefGroup.hide()
        self.DeltaGroup = QtWidgets.QGroupBox("Shift to Temperature:")
        self.DeltaFrame = QtWidgets.QGridLayout()
        self.DeltaLabel = wc.QLabel(u'Δδ [ppm]')
        self.DeltaFrame.addWidget(self.DeltaLabel, 0, 1)
        DeltaGO = QtWidgets.QPushButton("Go")
        self.DeltaFrame.addWidget(DeltaGO, 1, 0)
        DeltaGO.clicked.connect(self.shiftToTemp)
        self.Delta = wc.QLineEdit("")
        self.Delta.setMinimumWidth(100)
        self.DeltaFrame.addWidget(self.Delta, 1, 1)
        self.DeltaGroup.setLayout(self.DeltaFrame)
        grid1.addWidget(self.DeltaGroup, 2, 0)
        self.TempGroup = QtWidgets.QGroupBox("Temperature to Shift:")
        self.TempFrame = QtWidgets.QGridLayout()
        TempLabel = wc.QLabel(u'Temperature [K]')
        self.TempFrame.addWidget(TempLabel, 0, 1)
        TempGO = QtWidgets.QPushButton("Go")
        TempGO.clicked.connect(self.tempToShift)
        self.TempFrame.addWidget(TempGO, 1, 0)
        self.Temp = wc.QLineEdit("")
        self.Temp.setMinimumWidth(100)
        self.TempFrame.addWidget(self.Temp, 1, 1)
        self.TempGroup.setLayout(self.TempFrame)
        grid1.addWidget(self.TempGroup, 3, 0)
        self.refname = wc.QLabel(self.DEFINITIONS[0][7])
        grid1.addWidget(self.refname, 4, 0)
        # T1 based
        tab2 = QtWidgets.QWidget()
        tabWidget.addTab(tab2, "T1 Based")
        grid2 = QtWidgets.QGridLayout()
        tab2.setLayout(grid2)
        self.T1typeDrop = QtWidgets.QComboBox()
        self.T1typeDrop.addItems(self.T1_TEXTLIST)
        grid2.addWidget(self.T1typeDrop, 0, 0)
        self.T1typeDrop.currentIndexChanged.connect(self.changeTypeT1)
        self.T1Group = QtWidgets.QGroupBox("T1 to Temperature:")
        self.T1Frame = QtWidgets.QGridLayout()
        self.T1Label = wc.QLabel(u'T1 [s]')
        self.T1Frame.addWidget(self.T1Label, 0, 1)
        T1GO = QtWidgets.QPushButton("Go")
        self.T1Frame.addWidget(T1GO, 1, 0)
        T1GO.clicked.connect(self.t1ToTemp)
        self.T1 = wc.QLineEdit("")
        self.T1.setMinimumWidth(100)
        self.T1Frame.addWidget(self.T1, 1, 1)
        self.T1Group.setLayout(self.T1Frame)
        grid2.addWidget(self.T1Group, 2, 0)
        self.TempT1Group = QtWidgets.QGroupBox("Temperature to T1:")
        self.TempT1Frame = QtWidgets.QGridLayout()
        TempT1Label = wc.QLabel(u'Temperature [K]')
        self.TempT1Frame.addWidget(TempT1Label, 0, 1)
        TempT1GO = QtWidgets.QPushButton("Go")
        TempT1GO.clicked.connect(self.tempToT1)
        self.TempT1Frame.addWidget(TempT1GO, 1, 0)
        self.TempT1 = wc.QLineEdit("")
        self.TempT1.setMinimumWidth(100)
        self.TempT1Frame.addWidget(self.TempT1, 1, 1)
        self.TempT1Group.setLayout(self.TempT1Frame)
        grid2.addWidget(self.TempT1Group, 3, 0)
        self.refnameT1 = wc.QLabel(self.T1_DEFINITIONS[0][4])
        grid2.addWidget(self.refnameT1, 4, 0)
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(tabWidget, 0, 0, 1, 4)
        cancelButton = QtWidgets.QPushButton("&Close")
        cancelButton.clicked.connect(self.closeEvent)
        box = QtWidgets.QDialogButtonBox()
        box.addButton(cancelButton, QtWidgets.QDialogButtonBox.RejectRole)
        layout.addWidget(box, 1, 0, 1, 4)
        layout.setColumnStretch(3, 1)
#        layout.setRowStretch(3, 1)
        self.show()
        self.setFixedSize(self.size())

    def changeType(self, index):
        self.Temp.setText('')
        self.Delta.setText('')
        self.refname.setText(self.DEFINITIONS[index][7])
        if self.DEFINITIONS[self.typeDrop.currentIndex()][4] == 'absShift':
            self.DeltaLabel.setText(u'Δδ [ppm]')
            self.RefGroup.hide()
            self.Delta0.setText('')
            self.T0.setText('')
        elif self.DEFINITIONS[self.typeDrop.currentIndex()][4] == 'relShift':
            self.DeltaLabel.setText(u'δ [ppm]')
            self.RefGroup.show()
            self.Delta0.setText(self.DEFINITIONS[self.typeDrop.currentIndex()][5])
            self.T0.setText(self.DEFINITIONS[self.typeDrop.currentIndex()][6])
        self.father.root.processEvents()
        self.setFixedSize(self.sizeHint())

    def changeTypeT1(self, index):
        self.refnameT1.setText(self.T1_DEFINITIONS[index][4])
        self.father.root.processEvents()
        self.setFixedSize(self.sizeHint())

    def tempToT1(self):
        Data = self.T1_DEFINITIONS[self.T1typeDrop.currentIndex()]
        try:
            Temp = float(safeEval(self.TempT1.text(), Type='FI'))
        except Exception:
            self.T1.setText('?')
            raise SsnakeException("Temperature Calibration: Invalid input in Temp value")
        T1 = Data[3](Temp)
        if Temp < Data[0] or Temp > Data[1]:
            self.T1.setText('?')
            raise SsnakeException("Temperature Calibration: Temperature outside calibration range")
        self.T1.setText('%#.6g' % T1)

    def t1ToTemp(self):
        Data = self.T1_DEFINITIONS[self.T1typeDrop.currentIndex()]
        try:
            T1 = float(safeEval(self.T1.text(), Type='FI'))
        except Exception:
            self.TempT1.setText('?')
            raise SsnakeException("Temperature Calibration: Invalid input in Temp value")
        try:
            Temp = Data[2](T1)
        except Exception:
            self.TempT1.setText('?')
            raise SsnakeException("Temperature Calibration: Temperature outside calibration range")
        self.TempT1.setText('%#.6g' % Temp)

    def shiftToTemp(self):
        Data = self.DEFINITIONS[self.typeDrop.currentIndex()]
        try:
            Delta = float(safeEval(self.Delta.text(), Type='FI'))
        except Exception:
            self.Temp.setText('?')
            raise SsnakeException("Temperature Calibration: Invalid input in Delta value")
        if Data[4] == 'relShift':
            try:
                Delta0 = float(safeEval(self.Delta0.text(), Type='FI'))
                T0 = float(safeEval(self.T0.text(), Type='FI'))
            except Exception:
                self.Temp.setText('?')
                raise SsnakeException("Temperature Calibration: Invalid input in References values")
            Temp = Data[2](Delta, Delta0, T0)
        else:
            Temp = Data[2](Delta)
        if Temp < Data[0] or Temp > Data[1]:
            self.Temp.setText('?')
            raise SsnakeException("Temperature Calibration: Temperature outside calibration range")
        self.Temp.setText('%#.6g' % Temp)

    def tempToShift(self):
        Data = self.DEFINITIONS[self.typeDrop.currentIndex()]
        try:
            Temp = float(safeEval(self.Temp.text(), Type='FI'))
        except Exception:
            self.Delta.setText('?')
            raise SsnakeException("Temperature Calibration: Invalid input in Temp value")
        if Data[4] == 'relShift':
            try:
                Delta0 = float(safeEval(self.Delta0.text(), Type='FI'))
                T0 = float(safeEval(self.T0.text(), Type='FI'))
            except Exception:
                self.Delta.setText('?')
                raise SsnakeException("Temperature Calibration: Invalid input in References values")
            Delta = Data[3](Temp, Delta0, T0)
        else:
            Delta = Data[3](Temp)
        if Temp < Data[0] or Temp > Data[1]:
            self.Delta.setText('?')
            raise SsnakeException("Temperature Calibration: Temperature outside calibration range")
        self.Delta.setText('%#.6g' % Delta)

    def closeEvent(self, *args):
        self.deleteLater()

##############################################################################
class mqmasExtractWindow(wc.ToolWindow):

    Ioptions = ['3/2','5/2', '7/2', '9/2']
    Ivalues = [1.5, 2.5, 3.5, 4.5]
    z = [680.0/27.0, 8500.0/81.0, 6664.0/27.0, 1360.0/3.0]
    BdevA = [1.0/68.0, 3.0/850.0, 5.0/3332.0, 1.0/1224.0]
    NAME = "MQMAS"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent):
        super(mqmasExtractWindow, self).__init__(parent)
        self.comGroup = QtWidgets.QGroupBox()
        self.comFrame = QtWidgets.QGridLayout()
        self.comFrame.addWidget(wc.QLabel("I:"), 0, 0)
        self.IEntry = QtWidgets.QComboBox()
        self.IEntry.addItems(self.Ioptions)
        self.IEntry.setCurrentIndex(0)
        self.comFrame.addWidget(self.IEntry, 0, 1)
        self.comFrame.addWidget(wc.QLabel(u'ν' + '<sub>0</sub> [MHz]'), 1, 0)
        self.nu0 = wc.QLineEdit("0.0")
        self.comFrame.addWidget(self.nu0, 1, 1)
        self.comGroup.setLayout(self.comFrame)
        self.grid.addWidget(self.comGroup, 0, 0, 2, 2)
        self.onetwoGroup = QtWidgets.QGroupBox("δ1/δ2:")
        self.onetwoFrame = QtWidgets.QGridLayout()
        self.onetwoFrame.addWidget(wc.QLabel(u'δ' + '<sub>1</sub> [ppm]'), 2, 0)
        self.onetwoFrame.addWidget(wc.QLabel(u'δ' + '<sub>2</sub> [ppm]'), 3, 0)
        self.delta1 = wc.QLineEdit("0.0")
        self.onetwoFrame.addWidget(self.delta1, 2, 1)
        self.delta2 = wc.QLineEdit("0.0")
        self.onetwoFrame.addWidget(self.delta2, 3, 1)
        self.delta1.setMinimumWidth(200)
        self.calcIsoPqButton = QtWidgets.QPushButton("Calc δiso/PQ", self)
        self.calcIsoPqButton.clicked.connect(self.calcIsoPq)
        self.onetwoFrame.addWidget(self.calcIsoPqButton, 4, 0, 1, 2)
        self.onetwoGroup.setLayout(self.onetwoFrame)
        self.grid.addWidget(self.onetwoGroup, 2, 0, 4, 2)
        self.isopqGroup = QtWidgets.QGroupBox("δiso/PQ:")
        self.isopqFrame = QtWidgets.QGridLayout()
        self.isopqFrame.addWidget(wc.QLabel(u'δ' + '<sub>iso</sub> [ppm]'), 6, 0)
        self.deltaIso = wc.QLineEdit("0.0")
        self.isopqFrame.addWidget(self.deltaIso, 6, 1)
        self.isopqFrame.addWidget(wc.QLabel('P<sub>Q</sub> [MHz]'), 7, 0)
        self.pq = wc.QLineEdit("0.0")
        self.isopqFrame.addWidget(self.pq, 7, 1)
        self.calc12Button = QtWidgets.QPushButton("Calc δ1/δ2", self)
        self.calc12Button.clicked.connect(self.calc12)
        self.isopqFrame.addWidget(self.calc12Button, 8, 0, 1, 2)
        self.isopqGroup.setLayout(self.isopqFrame)
        self.grid.addWidget(self.isopqGroup, 6, 0, 4, 2)
        self.cancelButton.setText("Close")
        self.cancelButton.clicked.disconnect()
        self.cancelButton.clicked.connect(self.closeEvent)
        self.okButton.setText("Reset")
        self.okButton.clicked.disconnect()
        self.okButton.clicked.connect(self.valueReset)

    def calcIsoPq(self):
        nu0 = safeEval(self.nu0.text(), Type='FI')
        wrong = False
        if nu0 is None:
            self.father.dispMsg("MQMAS Extract: Invalid input in V0")
            wrong = True
        delta1 = safeEval(self.delta1.text(), Type='FI')
        if delta1 is None and wrong is False:
            self.father.dispMsg("MQMAS Extract: Invalid input in Delta1")
            wrong = True
        delta2 = safeEval(self.delta2.text(), Type='FI')
        if delta2 is None and wrong is False:
            self.father.dispMsg("MQMAS Extract: Invalid input in Delta2")
            wrong = True
        zval = self.z[self.IEntry.currentIndex()]
        if wrong is False and delta1 is not None and delta2 is not None:
            if delta1 - delta2 < 0.0:
                self.father.dispMsg("MQMAS Extract: Delta1 should be larger than Delta2!")
                wrong = True
        if wrong:
            self.deltaIso.setText('-')
            self.pq.setText('-')
            return
        iso = (17.0 * delta1 + 10.0 * delta2) / 27
        self.deltaIso.setText(str(iso))
        pq = np.sqrt(zval * 1e-6 * (nu0 * 1e6)**2 * (delta1 - delta2)) / 1e6
        self.pq.setText(str(pq))

    def calc12(self):
        nu0 = safeEval(self.nu0.text(), Type='FI')
        wrong = False
        if nu0 is None or nu0 == 0.0:
            self.father.dispMsg("MQMAS Extract: Invalid input in V0")
            wrong = True
        iso = safeEval(self.deltaIso.text(), Type='FI')
        if iso is None and wrong is False:
            self.father.dispMsg("MQMAS Extract: Invalid input in DeltaIso")
            wrong = True
        pq = safeEval(self.pq.text(), Type='FI')
        if pq is None and wrong is False:
            self.father.dispMsg("MQMAS Extract: Invalid input in PQ")
            wrong = True
        if wrong:
            self.delta1.setText('-')
            self.delta2.setText('-')
            return
        BdevA = self.BdevA[self.IEntry.currentIndex()]
        delta1 = iso + BdevA * pq**2/nu0**2 * 1e6
        self.delta1.setText(str(delta1))
        delta2 = (27 * iso -  17 * delta1) / 10
        self.delta2.setText(str(delta2))

    def valueReset(self):  # Resets all the boxes to 0
        self.nu0.setText('0')
        self.delta1.setText('0')
        self.delta2.setText('0')
        self.deltaIso.setText('0')
        self.pq.setText('0')
        self.IEntry.setCurrentIndex(0)

    def closeEvent(self, *args):
        self.deleteLater()


def libVersionChecker(version,needed):
    """Compares a two library version strings ('1.2.3' format)

        First compares major, then minor, etc.
        If version is lower than needed, False is returned
    """
    current = [int(x) for x in version.split('.')]
    required = [int(x) for x in needed.split('.')]
    check = True
    if current[0] < required[0]:
        check = False
    elif  current[0] == required[0]:
        if current[1] < required[1]:
            check = False
        elif  current[1] == required[1]:
            if len(current) > 2 and len(required) > 2:
                if current[2] < required[2]:
                    check = False
    return check


def checkVersions():
    """Checks versions of relevant python libraries

        Compares specified versions of libraries against the loaded version.
        If the values are to low, an error message is returned.
    """
    from scipy import __version__ as scipyVersion # Scipy is not fully imported, so only load version
    libs = [['numpy', np.__version__, NPVERSION],
            ['matplotlib', matplotlib.__version__, MPLVERSION],
            ['scipy', scipyVersion, SPVERSION]]
    libs.append(['python', str(sys.version_info.major) + '.' + str(sys.version_info.minor), PY3VERSION])
    messages = []
    error = False
    for elem in libs:
        check = libVersionChecker(elem[1], elem[2])
        if not check:
            error = True
            messages.append('"' + elem[0] + '" version is too low (need "' + elem[2] + '" have "' + elem[1] +'")')
    return error, messages

def popupVersionError(messages):
    """Gives a message window displaying version issues

        Input is a list of strings
    """
    msg = ""
    for elem in messages:
        msg = msg + elem + '\n'
    reply = QtWidgets.QMessageBox.warning(QtWidgets.QWidget(), 'Invalid software version', msg, QtWidgets.QMessageBox.Ignore, QtWidgets.QMessageBox.Abort)
    quit = False
    if reply == QtWidgets.QMessageBox.Abort:
        quit = True
    return quit

def openRefMan():
    file = os.path.dirname(os.path.realpath(__file__))  + os.path.sep + '..' + os.path.sep + 'ReferenceManual.pdf'
    if sys.platform.startswith('linux'):
        os.system("xdg-open " + '"' + file + '"')
    elif sys.platform.startswith('darwin'):
        os.system("open " + '"' + file + '"')
    elif sys.platform.startswith('win'):
        os.startfile(file)

def openTutorial():
    path = os.path.dirname(os.path.realpath(__file__))  + os.path.sep + '..' + os.path.sep + '/Tutorial'
    if sys.platform.startswith('linux'):
        os.system("xdg-open " + '"' + path + '"')
    elif sys.platform.startswith('darwin'):
        os.system("open " + '"' + path + '"')
    elif sys.platform.startswith('win'):
        os.startfile(path)

if __name__ == '__main__':
    error, messages = checkVersions()
    quit = False
    if error:
        splash.close()
        quit = popupVersionError(messages)
    if not quit:
        mainProgram = MainProgram(root)
        mainProgram.setWindowTitle("ssNake - " + VERSION)
        mainProgram.show()
        if not error:
            splash.finish(mainProgram)
        sys._excepthook = sys.excepthook
        if len(sys.argv) > 1:
            mainProgram.loadData(sys.argv[1:])

        def exception_hook(exctype, value, traceback):
            if not isinstance(value, Exception): # Do not catch keyboard interrupts
                sys._excepthook(exctype, value, traceback)
            elif isinstance(value, (sc.SpectrumException, hc.HComplexException, sim.SimException)):
                mainProgram.dispMsg(str(value))
            else:
                mainProgram.dispError([exctype, value, traceback])
        sys.excepthook = exception_hook
        sys.exit(root.exec_())

