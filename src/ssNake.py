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


import sip
import sys
import os
sip.setapi('QString', 2)
try:
    from PyQt4 import QtGui, QtCore
    from PyQt4 import QtGui as QtWidgets
    QT = 4
except ImportError:
    from PyQt5 import QtGui, QtCore, QtWidgets
    QT = 5


#Create splash window
if __name__ == '__main__':
    root = QtWidgets.QApplication(sys.argv)
    root.setWindowIcon(QtGui.QIcon(os.path.dirname(os.path.realpath(__file__)) + '/logo.gif'))
    splash_pix = QtGui.QPixmap(os.path.dirname(os.path.realpath(__file__)) + '/logo.gif')
    splash = QtWidgets.QSplashScreen(splash_pix, QtCore.Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    progressBar = QtWidgets.QProgressBar(splash)
    progressBar.setGeometry(2.5*splash.width()/10, 0.89*splash.height(),5*splash.width()/10, splash.height()/20)
    splash.show()

    
splashSteps=14.0/100
splashStep = 0.0
def splashProgressStep(splashStep): #A function to easily increase the progressbar value
    if __name__ == '__main__':
        splashStep=splashStep+1
        progressBar.setValue(splashStep // splashSteps + (splashStep % splashSteps > 0)) #Rounds up without math or numpy module
        root.processEvents()   
    return splashStep

  
import matplotlib  
splashStep = splashProgressStep(splashStep)
if QT ==4:
    matplotlib.use('Qt4Agg')
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
else:
    matplotlib.use('Qt5Agg')
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
splashStep = splashProgressStep(splashStep)
from matplotlib.figure import Figure
splashStep = splashProgressStep(splashStep)
import numpy as np
splashStep = splashProgressStep(splashStep)
import re
splashStep = splashProgressStep(splashStep)
import copy
splashStep = splashProgressStep(splashStep)
import spectrum_classes as sc
splashStep = splashProgressStep(splashStep)
import fitting as fit
splashStep = splashProgressStep(splashStep)
from safeEval import safeEval
splashStep = splashProgressStep(splashStep)
import widgetClasses as wc
splashStep = splashProgressStep(splashStep)
from updateWindow import UpdateWindow
splashStep = splashProgressStep(splashStep)
from plotWindow import MainPlotWindow
splashStep = splashProgressStep(splashStep)
from euro import euro
splashStep = splashProgressStep(splashStep)
import scipy.constants as SC
splashStep = splashProgressStep(splashStep)


matplotlib.rc('font', family='DejaVu Sans')
np.set_printoptions(threshold=np.nan)
QtCore.QLocale.setDefault(QtCore.QLocale('en_US'))

pi = np.pi

VERSION = 'v0.6b'


class MainProgram(QtWidgets.QMainWindow):

    def __init__(self, root):
        QtWidgets.QMainWindow.__init__(self)
        self.root = root
        self.VERSION = VERSION
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
        self.LastLocation = ''
        self.initMenu()
        self.menuCheck()
#        self.initToolbar()
        self.main_widget = QtWidgets.QWidget(self)
        self.mainFrame = QtWidgets.QGridLayout(self.main_widget)
        self.logo = QtWidgets.QLabel(self)
        self.logo.setPixmap(QtGui.QPixmap(os.path.dirname(os.path.realpath(__file__)) + "/logo.gif"))
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
        self.eventFilter = wc.MyEventFilter(self)
        self.root.installEventFilter(self.eventFilter)
        self.loadDefaults()
        self.resize(self.defaultWidth, self.defaultHeight)
        if self.defaultMaximized:
            self.showMaximized()
        QtWidgets.QShortcut(QtGui.QKeySequence.Paste, self).activated.connect(self.handlePaste)
        QtWidgets.QShortcut(QtGui.QKeySequence.Copy, self).activated.connect(self.handleCopy)

    def handlePaste(self):
        self.dropEvent(QtWidgets.QApplication.instance().clipboard())

    def handleCopy(self):
        if self.mainWindow is None:
            return
        pixmap = QtGui.QPixmap.grabWidget(self.mainWindow.canvas)
        QtWidgets.QApplication.clipboard().setPixmap(pixmap)
        
    def resetDefaults(self):
        self.defaultWidth = 1
        self.defaultHeight = 1
        self.defaultMaximized = False
        self.defaultAskName = True
        self.defaultLinewidth = 1.0
        self.defaultColor = '#0000FF'
        self.defaultGrids = [False, False]
        self.defaultColorMap = 'seismic'
        self.defaultWidthRatio = 3.0
        self.defaultHeightRatio = 3.0
        self.defaultContourConst = True
        self.defaultPosColor = '#FF0000'
        self.defaultNegColor = '#0000FF'

    def loadDefaults(self):
        self.resetDefaults()
        QtCore.QSettings.setDefaultFormat(QtCore.QSettings.IniFormat)
        QtCore.QCoreApplication.setOrganizationName("ssNake")
        QtCore.QCoreApplication.setApplicationName("ssNake")
        settings = QtCore.QSettings()
        self.defaultColor = settings.value("plot/color", self.defaultColor, str)
        try:
            self.defaultLinewidth = settings.value("plot/linewidth", self.defaultLinewidth, float)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the plot/linewidth")
        self.defaultGrids = [settings.value("plot/xgrid", self.defaultGrids[0], bool), settings.value("plot/ygrid", self.defaultGrids[1], bool)]
        self.defaultColorMap = settings.value("contour/colormap", self.defaultColorMap, str)
        self.defaultContourConst = settings.value("contour/constantcolors", self.defaultContourConst, bool)
        self.defaultPosColor = settings.value("contour/poscolor", self.defaultPosColor, str)
        self.defaultNegColor = settings.value("contour/negcolor", self.defaultNegColor, str)
        if not str(self.defaultColorMap) in sc.COLORMAPLIST:
            self.dispMsg("Incorrect colormap in config file")
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
        try:
            self.defaultWidthRatio = settings.value("contour/width_ratio", self.defaultWidthRatio, float)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the contour/width_ratio")
        try:
            self.defaultHeightRatio = settings.value("contour/height_ratio", self.defaultHeightRatio, float)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the contour/height_ratio")

    def saveDefaults(self):
        QtCore.QSettings.setDefaultFormat(QtCore.QSettings.IniFormat)
        QtCore.QCoreApplication.setOrganizationName("ssNake")
        QtCore.QCoreApplication.setApplicationName("ssNake")
        settings = QtCore.QSettings()
        settings.setValue("plot/color", self.defaultColor)
        settings.setValue("plot/linewidth", self.defaultLinewidth)
        settings.setValue("plot/xgrid", self.defaultGrids[0])
        settings.setValue("plot/ygrid", self.defaultGrids[1])
        settings.setValue("maximized", self.defaultMaximized)
        settings.setValue("width", self.defaultWidth)
        settings.setValue("height", self.defaultHeight)
        settings.setValue("ask_name", self.defaultAskName)
        settings.setValue("contour/colormap", self.defaultColorMap)
        settings.setValue("contour/constantcolors", self.defaultContourConst)
        settings.setValue("contour/poscolor", self.defaultPosColor)
        settings.setValue("contour/negcolor", self.defaultNegColor)
        settings.setValue("contour/width_ratio", self.defaultWidthRatio)
        settings.setValue("contour/height_ratio", self.defaultHeightRatio)

    def dispMsg(self, msg):
        self.statusBar.showMessage(msg, 10000)
    
    
    def initToolbar(self):
        self.toolbar = self.addToolBar('Toolbar')
        self.toolbar.setMovable(False)
        self.toolbar.setIconSize(QtCore.QSize(22,22))
        

        
        self.toobarActionList = [self.openAct,self.saveAct,self.saveMatAct,self.savefigAct,
                                 self.saveSimpsonAct,self.saveASCIIAct,self.saveASCIIAct,self.preferencesAct,
                                 self.quitAct,self.newAct,self.closeAct,self.renameWorkspaceAct,self.forwardAct,
                                 self.backAct,self.combineWorkspaceAct,self.macrostartAct,self.macrostopAct,
                                 self.macroLoadAct,self.undoAction,self.redoAction,self.reloadAct,self.monitorAct,
                                 self.realAct,self.imagAct,self.absAct,self.apodizeAct,self.phaseAct,
                                 self.swapEchoAct,self.corOffsetAct,self.baselineAct,self.subAvgAct,self.refDeconvAct,
                                 self.statesAct,self.statesTPPIAct,self.echoantiAct,self.brukDigitalAct,
                                 self.lpsvdAct,self.sizingAct,self.shiftAct,self.intRegionAct,
                                 self.sumRegionAct,self.maxRegionAct,self.minRegionAct,self.maxposRegionAct,
                                 self.minposRegionAct,self.averageRegionAct,self.diffAct,self.cumsumAct,
                                 self.extractpartAct,self.fliplrAct, self.matrixdelAct,self.splitAct,
                                 self.multiplyAct,self.reorderAct,self.concatAct,self.shearAct,
                                 self.fourierAct,self.realFourierAct,self.fftshiftAct,self.invfftshiftAct,
                                 self.hilbertAct,self.ffmAct,self.cleanAct,self.snrAct,self.fwhmAct,
                                 self.massAct,self.intfitAct,self.relaxAct,self.diffusionAct,self.lorentzfitAct,
                                 self.csastaticAct,self.csamasAct,self.firstquadstatAct,self.firstquadmasAct,
                                 self.secondquadstatAct,self.secondquadmasAct,self.czjzekstatAct,self.czjzekmasAct,
                                 self.insertdatAct,self.adddatAct,self.subdatAct,self.onedplotAct,self.scatterplotAct,
                                 self.stackplotAct,self.arrayplotAct,self.contourplotAct,self.skewplotAct,self.multiplotAct,
                                 self.setrefAct, self.delrefAct,self.loadrefAct,self.userxAct,self.plotprefAct]
       
      
#        self.emptyAction = QtGui.QAction('', self)
#        self.emptyAction.setEnabled(False)
#        self.toolbar.addAction(self.emptyAction)
#        
#        self.seperatorAction = QtGui.QAction(self)
#        self.seperatorAction.setSeparator(True)
#        self.toolbar.addAction(self.seperatorAction)
        
       
        for entry in self.toobarActionList:
                    self.toolbar.addAction(entry)
    
        
        
        
    def initMenu(self):
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        self.menubar = self.menuBar()
        self.filemenu = QtWidgets.QMenu('&File', self)
        self.menubar.addMenu(self.filemenu)
        self.openAct = self.filemenu.addAction(QtGui.QIcon(IconDirectory + 'open.png'), '&Open', self.loadFromMenu, QtGui.QKeySequence.Open)
        self.openAct.setToolTip('Open a File')
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
        self.saveASCIIAct = self.exportmenu.addAction(QtGui.QIcon(IconDirectory + 'ssnake.png'), 'ASCII (1D/2D)', self.saveASCIIFile)
        self.saveASCIIAct.setToolTip('Save as ASCII Text File')        
        self.preferencesAct = self.filemenu.addAction(QtGui.QIcon(IconDirectory + 'preferences.png'), '&Preferences', lambda: PreferenceWindow(self))
        self.preferencesAct.setToolTip('Open Preferences Window')
        self.quitAct = self.filemenu.addAction(QtGui.QIcon(IconDirectory + 'quit.png'), '&Quit', self.fileQuit, QtGui.QKeySequence.Quit)
        self.quitAct.setToolTip('Close ssNake')

        # Workspaces menu
        self.workspacemenu = QtWidgets.QMenu('&Workspaces', self)
        self.menubar.addMenu(self.workspacemenu)
        self.newAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'duplicate.png'), 'D&uplicate', self.duplicateWorkspace, QtGui.QKeySequence.New)
        self.newAct.setToolTip('Duplicate Workspace')
        self.closeAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'), '&Delete', self.destroyWorkspace, QtGui.QKeySequence.Close)
        self.closeAct.setToolTip('Delete Workspace')
        self.renameWorkspaceAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'), '&Rename', self.renameWorkspace, QtCore.Qt.Key_F2)
        self.renameWorkspaceAct.setToolTip('Rename Workspace')
        self.activemenu = QtWidgets.QMenu('&Active', self)
        self.workspacemenu.addMenu(self.activemenu)
        self.forwardAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'next.png'), '&Next', lambda: self.stepWorkspace(1), QtGui.QKeySequence.Forward)
        self.forwardAct.setToolTip('Next Workspace')
        self.backAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'previous.png'), '&Previous', lambda: self.stepWorkspace(-1), QtGui.QKeySequence.Back)
        self.backAct.setToolTip('Previous Workspace')
        self.combineWorkspaceAct = self.workspacemenu.addAction(QtGui.QIcon(IconDirectory + 'combine.png'), '&Combine', self.createCombineWorkspaceWindow)
        self.combineWorkspaceAct.setToolTip('Combine Workspaces')

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
        self.reloadAct = self.editmenu.addAction(QtGui.QIcon(IconDirectory + 'reload.png'), "Re&load", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.reloadLast()), QtGui.QKeySequence.Refresh)
        self.reloadAct.setToolTip('Reload Current Data')
        self.monitorAct = self.editmenu.addAction(QtGui.QIcon(IconDirectory + 'monitor.png'),"&Monitor", lambda: self.mainWindowCheck(lambda mainWindow: MonitorWindow(mainWindow)))
        self.monitorAct.setToolTip('Monitor Current Data')
        
        # the tool drop down menu
        self.toolMenu = QtWidgets.QMenu("&Tools", self)
        self.menubar.addMenu(self.toolMenu)
        self.realAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'real.png'), "&Real", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.real()))
        self.realAct.setToolTip('Take Real Part of Data')
        self.imagAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'imag.png'), "&Imag", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.imag()))
        self.imagAct.setToolTip('Take Imaginary Part of Data')
        self.absAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'abs.png'), "&Abs", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.abs()))
        self.absAct.setToolTip('Take Absolute of Data')
        self.apodizeAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'apodize.png'),"Apo&dize", lambda: self.mainWindowCheck(lambda mainWindow: ApodWindow(mainWindow)))
        self.apodizeAct.setToolTip('Open Apodize Window')
        self.phaseAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'phase.png'), "&Phasing", lambda: self.mainWindowCheck(lambda mainWindow: PhaseWindow(mainWindow)))
        self.phaseAct.setToolTip('Open Phasing Window')
        self.swapEchoAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'fliplr.png'), "Swap &Echo", lambda: self.mainWindowCheck(lambda mainWindow: SwapEchoWindow(mainWindow)))
        self.swapEchoAct.setToolTip('Swap Echo')
        self.corOffsetAct = self.toolMenu.addAction("&Offset Correction", lambda: self.mainWindowCheck(lambda mainWindow: DCWindow(mainWindow)))
        self.corOffsetAct.setToolTip('Offset Correction')
        self.baselineAct = self.toolMenu.addAction("&Baseline Correction", lambda: self.mainWindowCheck(lambda mainWindow: BaselineWindow(mainWindow)))
        self.baselineAct.setToolTip('Baseline Correction')
        self.subAvgAct = self.toolMenu.addAction("S&ubtract Averages", lambda: self.mainWindowCheck(lambda mainWindow: SubtractAvgWindow(mainWindow)))
        self.subAvgAct.setToolTip('Subtract Averages')
        self.refDeconvAct = self.toolMenu.addAction("Re&ference Deconvolution", lambda: self.mainWindowCheck(lambda mainWindow: FiddleWindow(mainWindow)))
        self.refDeconvAct.setToolTip('Reference Deconvolution')
        self.statesAct = self.toolMenu.addAction("&States", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.states()))
        self.statesAct.setToolTip('States Hypercomplex Data Processing')
        self.statesTPPIAct = self.toolMenu.addAction("States-&TPPI", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.statesTPPI()))
        self.statesTPPIAct.setToolTip('States-TPPI Hypercomplex Data Processing')
        self.echoantiAct = self.toolMenu.addAction("Ec&ho-antiecho", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.echoAntiEcho()))
        self.echoantiAct.setToolTip('Ec&ho-antiecho Hypercomplex Data Processing')
        self.brukDigitalAct = self.toolMenu.addAction("&Correct Bruker Digital Filter", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.BrukerDigital()))
        self.brukDigitalAct.setToolTip("Correct Bruker Digital Filter")
        self.lpsvdAct = self.toolMenu.addAction("&LPSVD", lambda: self.mainWindowCheck(lambda mainWindow: LPSVDWindow(mainWindow)))
        self.lpsvdAct.setToolTip('LPSVD linear prediction')

        # the matrix drop down menu
        self.matrixMenu = QtWidgets.QMenu("M&atrix", self)
        self.menubar.addMenu(self.matrixMenu)
        self.sizingAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'sizing.png'),"&Sizing", lambda: self.mainWindowCheck(lambda mainWindow: SizeWindow(mainWindow)))
        self.sizingAct.setToolTip('Set Size')
        self.shiftAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'shift.png'), "S&hift Data", lambda: self.mainWindowCheck(lambda mainWindow: ShiftDataWindow(mainWindow)))
        self.shiftAct.setToolTip('Shift Data')
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
        self.maxposRegionAct = self.regionMenu.addAction("Ma&x position", lambda: self.mainWindowCheck(lambda mainWindow: argmaxWindow(mainWindow)))
        self.maxposRegionAct.setToolTip('Position of Maximum of Region')
        self.minposRegionAct = self.regionMenu.addAction("Mi&n position", lambda: self.mainWindowCheck(lambda mainWindow: argminWindow(mainWindow)))
        self.minposRegionAct.setToolTip('Position of Minimum of Region')
        self.averageRegionAct = self.regionMenu.addAction(QtGui.QIcon(IconDirectory + 'average.png'), "&Average", lambda: self.mainWindowCheck(lambda mainWindow: avgWindow(mainWindow)))
        self.averageRegionAct.setToolTip('Average of Region')
        self.diffAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'diff.png'), "&Diff", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.diff()))
        self.diffAct.setToolTip('Difference')
        self.cumsumAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'cumsum.png'), "&Cumsum", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.cumsum()))
        self.cumsumAct.setToolTip('Cumulative sum')
        self.extractpartAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'extractpart.png'),"&Extract part", lambda: self.mainWindowCheck(lambda mainWindow: extractRegionWindow(mainWindow)))
        self.extractpartAct.setToolTip('Extract part')
        self.fliplrAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'fliplr.png'), "&Flip L/R", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.flipLR()))
        self.fliplrAct.setToolTip('Flip L/R')
        self.matrixdelAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'matrixdelete.png'), "De&lete", lambda: self.mainWindowCheck(lambda mainWindow: DeleteWindow(mainWindow)))
        self.matrixdelAct.setToolTip('Delete Points')
        self.splitAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'split.png'),"S&plit", lambda: self.mainWindowCheck(lambda mainWindow: SplitWindow(mainWindow)))
        self.splitAct.setToolTip('Split')
        self.multiplyAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'multiply.png'), "Mul&tiply", lambda: self.mainWindowCheck(lambda mainWindow: MultiplyWindow(mainWindow)))
        self.multiplyAct.setToolTip('Multiply')
        self.reorderAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'reorder.png'), "&Reorder", lambda: self.mainWindowCheck(lambda mainWindow: ReorderWindow(mainWindow)))
        self.reorderAct.setToolTip('Reorder')
        self.concatAct = self.matrixMenu.addAction("C&oncatenate", lambda: self.mainWindowCheck(lambda mainWindow: ConcatenateWindow(mainWindow)))
        self.concatAct.setToolTip('Concatenate')
        self.multiDActions.append(self.concatAct)
        self.shearAct = self.matrixMenu.addAction("Shearin&g", lambda: self.mainWindowCheck(lambda mainWindow: ShearingWindow(mainWindow)))
        self.shearAct.setToolTip('Shearing')
        self.multiDActions.append(self.shearAct)


        # the fft drop down menu
        self.fftMenu = QtWidgets.QMenu("T&ransforms", self)
        self.menubar.addMenu(self.fftMenu)
        self.fourierAct = self.fftMenu.addAction("&Fourier Transform", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.fourier()), QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        self.fourierAct.setToolTip('Fourier Transform')
        self.realFourierAct = self.fftMenu.addAction("&Real Fourier Transform", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.realFourier()))
        self.realFourierAct.setToolTip('Real Fourier Transform')
        self.fftshiftAct = self.fftMenu.addAction("Fft&shift", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.fftshift()))
        self.fftshiftAct.setToolTip('Fftshift')
        self.invfftshiftAct = self.fftMenu.addAction("&Inv fftshift", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.invFftshift()))
        self.invfftshiftAct.setToolTip('Inverse fftshift')
        self.hilbertAct = self.fftMenu.addAction("&Hilbert Transform", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.hilbert()))
        self.hilbertAct.setToolTip('Hilbert Transform') 
        self.nusMenu = QtWidgets.QMenu("&NUS", self)
        self.fftMenu.addMenu(self.nusMenu)
        self.ffmAct = self.nusMenu.addAction("&FFM", lambda: self.mainWindowCheck(lambda mainWindow: FFMWindow(mainWindow)))
        self.ffmAct.setToolTip('FFM') 
        self.cleanAct = self.nusMenu.addAction("&CLEAN", lambda: self.mainWindowCheck(lambda mainWindow: CLEANWindow(mainWindow)))
        self.cleanAct.setToolTip('CLEAN') 

        # the fitting drop down menu
        self.fittingMenu = QtWidgets.QMenu("F&itting", self)
        self.menubar.addMenu(self.fittingMenu)
        self.snrAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'snr.png'),"&S/N", lambda: self.mainWindowCheck(lambda mainWindow: SNWindow(mainWindow)))
        self.snrAct.setToolTip('Signal-to-Noise Ratio')
        self.fwhmAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'fwhm.png'),"&FWHM", lambda: self.mainWindowCheck(lambda mainWindow: FWHMWindow(mainWindow)))
        self.fwhmAct.setToolTip('Full Width at Half Maximum')
        self.massAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'mass.png'),"Centre of Mass", lambda: self.mainWindowCheck(lambda mainWindow: COMWindow(mainWindow)))
        self.massAct.setToolTip('Centre of Mass')
        self.intfitAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'int.png'),"&Integrals", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createIntegralsWindow()))
        self.intfitAct.setToolTip('Get Integrals')
        self.relaxAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'relaxation.png'),"&Relaxation Curve", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createRelaxWindow()))
        self.relaxAct.setToolTip('Fit Relaxation Curve')
        self.diffusionAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'diffusion.png'),"&Diffusion Curve", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createDiffusionWindow()))
        self.diffusionAct.setToolTip('Fit Diffusion Curve')
        self.lorentzfitAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'lorentz.png'),"&Lorentzian/Gaussian", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createPeakDeconvWindow()))
        self.lorentzfitAct.setToolTip('Fit Lorentzian/Gaussian')
        self.csastaticAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'csastatic.png'),"&CSA Static", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createTensorDeconvWindow()))
        self.csastaticAct.setToolTip('Fit CSA Static')
        self.csamasAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'csaMAS.png'),"CSA MAS", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createHerzfeldBergerWindow()))
        self.csamasAct.setToolTip('Fit CSA MAS')
        self.firstquadstatAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'firstquadstatic.png'),"First Order &Quadrupole Static", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuad1DeconvWindow()))
        self.firstquadstatAct.setToolTip('Fit First Order Quadrupole Static')
        self.firstquadmasAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'firstquadMAS.png'),"First Order &Quadrupole MAS", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuad1MASDeconvWindow()))
        self.firstquadmasAct.setToolTip('Fit First Order Quadrupole MAS')
        self.secondquadstatAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'secondquadstatic.png'),"S&econd Order Quadrupole Static", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuad2StaticDeconvWindow()))
        self.secondquadstatAct.setToolTip('Fit Second Order Quadrupole Static')
        self.secondquadmasAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'secondquadMAS.png'),"Se&cond Order Quadrupole MAS", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuad2MASDeconvWindow()))
        self.secondquadmasAct.setToolTip('Fit Second Order Quadrupole MAS')
        self.czjzekstatAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'czjzekstatic.png'),"Czjzek S&tatic", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuad2StaticCzjzekWindow()))
        self.czjzekstatAct.setToolTip('Fit Static Czjzek Pattern')
        self.czjzekmasAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'czjzekMAS.png'),"Czjzek &MAS", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuad2MASCzjzekWindow()))
        self.czjzekmasAct.setToolTip('Fit MAS Czjzek Pattern')

        # the combine drop down menu
        self.combineMenu = QtWidgets.QMenu("Com&bine", self)
        self.menubar.addMenu(self.combineMenu)
        self.insertdatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'insert.png'),"&Insert From Workspace", lambda: self.mainWindowCheck(lambda mainWindow: InsertWindow(mainWindow)))
        self.insertdatAct.setToolTip('Insert From Workspace')
        self.adddatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'add.png'), "&Add", lambda: self.mainWindowCheck(lambda mainWindow: AddWindow(mainWindow)))
        self.adddatAct.setToolTip('Add Data From Workspace')
        self.subdatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'subtract.png'), "&Subtract", lambda: self.mainWindowCheck(lambda mainWindow: SubtractWindow(mainWindow)))
        self.subdatAct.setToolTip('Subtract Data From Workspace')

        # the plot drop down menu
        self.plotMenu = QtWidgets.QMenu("&Plot", self)
        self.menubar.addMenu(self.plotMenu)
        self.onedplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + '1dplot.png'), "&1D Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plot1D()))
        self.onedplotAct.setToolTip('1D plot')
        self.scatterplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'scatterplot.png'), "&Scatter Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotScatter()))
        self.scatterplotAct.setToolTip('Scatter Plot')
        self.stackplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'stack.png'),"S&tack Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotStack()))
        self.stackplotAct.setToolTip('Stack Plot')        
        self.multiDActions.append(self.stackplotAct)
        self.arrayplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'array.png'),"&Array Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotArray()))
        self.arrayplotAct.setToolTip('Array Plot')        
        self.multiDActions.append(self.arrayplotAct)
        self.contourplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'contour.png'), "&Contour Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotContour()))      
        self.contourplotAct.setToolTip('Contour Plot')
        self.multiDActions.append(self.contourplotAct)
        self.skewplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'skewed.png'),"S&kewed Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotSkewed()))
        self.skewplotAct.setToolTip('Skew Plot')
        self.multiDActions.append(self.skewplotAct)
        self.multiplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'multi.png'),"&Multi Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotMulti()))
        self.multiplotAct.setToolTip('Multi Plot')

        self.referencelistmenu = QtWidgets.QMenu('&Reference', self)
        self.plotMenu.addMenu(self.referencelistmenu)
        self.setrefAct = self.referencelistmenu.addAction(QtGui.QIcon(IconDirectory + 'setreference.png'),"&Set Reference", lambda: self.mainWindowCheck(lambda mainWindow: RefWindow(mainWindow)))
        self.setrefAct.setToolTip('Set Reference')
        self.delrefAct = self.referencelistmenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'),"&Clear Current Reference", self.referenceClear)
        self.delrefAct.setToolTip('Clear Current Reference')
        self.referencerunmenu = QtWidgets.QMenu('&Run', self)
        self.referencelistmenu.addMenu(self.referencerunmenu)
        self.referencedeletemenu = QtWidgets.QMenu('&Delete', self)
        self.referencelistmenu.addMenu(self.referencedeletemenu)
        self.referencerenamemenu = QtWidgets.QMenu('Re&name', self)
        self.referencelistmenu.addMenu(self.referencerenamemenu)
        self.referencesavemenu = QtWidgets.QMenu('&Save', self)
        self.referencelistmenu.addMenu(self.referencesavemenu)
        self.loadrefAct = self.referencelistmenu.addAction(QtGui.QIcon(IconDirectory + 'open.png'), "&Load", self.referenceLoad)
        self.loadrefAct.setToolTip('Load Reference')

        self.userxAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'xaxis.png'),"&User X-axis", lambda: self.mainWindowCheck(lambda mainWindow: XaxWindow(mainWindow)))
        self.userxAct.setToolTip('User X-axis')
        self.plotprefAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'preferences.png'),"&Plot Settings", lambda: self.mainWindowCheck(lambda mainWindow: PlotSettingsWindow(mainWindow)))
        self.plotprefAct.setToolTip('Plot Settings')

        # the history drop down menu
        self.historyMenu = QtWidgets.QMenu("&History", self)
        self.menubar.addMenu(self.historyMenu)
        self.historyMenu.addAction(QtGui.QIcon(IconDirectory + 'history.png'), "&History", lambda: self.mainWindowCheck(lambda mainWindow: HistoryWindow(mainWindow)))
        self.historyMenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'),"&Clear Undo/Redo List", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.clearUndo()))
        
        
        # the help drop down menu
        self.helpMenu = QtWidgets.QMenu("&Help", self)
        self.menubar.addMenu(self.helpMenu)
        self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'update.png'),"&Update", self.updateMenu)
        self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'shifttool.png'),"&Chemical Shift Conversion Tool", self.createShiftConversionWindow)
        self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'quadconversion.png'),"&Quadrupole Coupling Conversion Tool", self.createQuadConversionWindow)
        self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'table.png'),"&NMR Table", self.nmrTable)
        
        self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'about.png'),"&About", self.about)

    def mainWindowCheck(self, transfer):
        # checks if mainWindow exist to execute the function
        if self.mainWindow is not None:
            transfer(self.mainWindow)
        else:
            self.dispMsg("No workspaces open")

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if path.endswith('.zip'):
                import tempfile
                import shutil
                import zipfile
                try:
                    temp_dir = tempfile.mkdtemp()
                    zipfile.ZipFile(path).extractall(temp_dir)
                    for i in os.listdir(temp_dir): #Send the original path too,  for the workspace name
                        if self.autoLoad(os.path.join(temp_dir, i),realpath=path):
                            break
                finally:
                    shutil.rmtree(temp_dir)
            else:
                self.autoLoad(path)
            if path != '':  # if not cancelled
                self.LastLocation = os.path.dirname(path)  # Save used path

    def menuCheck(self):
        if self.mainWindow is None:
            self.savemenu.menuAction().setEnabled(False)
            self.exportmenu.menuAction().setEnabled(False)
            self.workspacemenu.menuAction().setEnabled(False)
            self.macromenu.menuAction().setEnabled(False)
            self.editmenu.menuAction().setVisible(False)
            self.toolMenu.menuAction().setVisible(False)
            self.matrixMenu.menuAction().setVisible(False)
            self.fftMenu.menuAction().setVisible(False)
            self.fittingMenu.menuAction().setVisible(False)
            self.combineMenu.menuAction().setVisible(False)
            self.plotMenu.menuAction().setVisible(False)
            self.historyMenu.menuAction().setVisible(False)
        else:
            self.editmenu.menuAction().setVisible(True)
            self.toolMenu.menuAction().setVisible(True)
            self.matrixMenu.menuAction().setVisible(True)
            self.fftMenu.menuAction().setVisible(True)
            self.fittingMenu.menuAction().setVisible(True)
            self.combineMenu.menuAction().setVisible(True)
            self.plotMenu.menuAction().setVisible(True)
            self.historyMenu.menuAction().setVisible(True)
            if isinstance(self.mainWindow, Main1DWindow):
                self.menuEnable()
                if (len(self.mainWindow.masterData.data.shape) < 2):
                    for i in self.multiDActions:
                        i.setEnabled(False)
                else:
                    for i in self.multiDActions:
                        i.setEnabled(True)
                if not self.mainWindow.undoList:
                    self.undoAction.setEnabled(False)
                else:
                    self.undoAction.setEnabled(True)
                if not self.mainWindow.redoList:
                    self.redoAction.setEnabled(False)
                else:
                    self.redoAction.setEnabled(True)
                self.savemenu.menuAction().setEnabled(True)
                self.exportmenu.menuAction().setEnabled(True)
                self.savefigAct.setEnabled(True)
                self.macromenu.menuAction().setEnabled(True)
                if self.mainWindow.currentMacro is None:
                    self.macrostopAct.setEnabled(False)
                    self.macrostartAct.setEnabled(True)
                else:
                    self.macrostopAct.setEnabled(True)
                    self.macrostartAct.setEnabled(False)
                self.savemenu.menuAction().setEnabled(True)
                self.exportmenu.menuAction().setEnabled(True)
                self.workspacemenu.menuAction().setEnabled(True)
            elif isinstance(self.mainWindow, MainPlotWindow):
                self.menuDisable(True)
                self.savemenu.menuAction().setEnabled(True)
                self.exportmenu.menuAction().setEnabled(True)
                self.savefigAct.setEnabled(False)
                self.workspacemenu.menuAction().setEnabled(True)
                self.macromenu.menuAction().setEnabled(False)
            else:
                self.menuDisable(True)
                self.savemenu.menuAction().setEnabled(True)
                self.exportmenu.menuAction().setEnabled(True)
                self.savefigAct.setEnabled(True)
                self.workspacemenu.menuAction().setEnabled(True)
                self.macromenu.menuAction().setEnabled(False)

    def menuEnable(self, internalWindow=False):
        self.macromenu.menuAction().setEnabled(True)
        self.editmenu.menuAction().setEnabled(True)
        self.toolMenu.menuAction().setEnabled(True)
        self.matrixMenu.menuAction().setEnabled(True)
        self.fftMenu.menuAction().setEnabled(True)
        self.fittingMenu.menuAction().setEnabled(True)
        self.combineMenu.menuAction().setEnabled(True)
        self.plotMenu.menuAction().setEnabled(True)
        self.historyMenu.menuAction().setEnabled(True)
        if not internalWindow:
            self.filemenu.menuAction().setEnabled(True)
            self.workspacemenu.menuAction().setEnabled(True)
            self.openAct.setEnabled(True)
            self.savefigAct.setEnabled(True)
            self.saveAct.setEnabled(True)
            self.newAct.setEnabled(True)
            self.closeAct.setEnabled(True)
            self.forwardAct.setEnabled(True)
            self.backAct.setEnabled(True)
            for i in range(self.tabs.count()):
                self.tabs.setTabEnabled(i, True)
        self.undoAction.setEnabled(True)
        self.redoAction.setEnabled(True)

    def menuDisable(self, internalWindow=False):
        self.macromenu.menuAction().setEnabled(False)
        self.editmenu.menuAction().setEnabled(False)
        self.toolMenu.menuAction().setEnabled(False)
        self.matrixMenu.menuAction().setEnabled(False)
        self.fftMenu.menuAction().setEnabled(False)
        self.fittingMenu.menuAction().setEnabled(False)
        self.combineMenu.menuAction().setEnabled(False)
        self.plotMenu.menuAction().setEnabled(False)
        self.historyMenu.menuAction().setEnabled(False)
        if not internalWindow:
            self.filemenu.menuAction().setEnabled(False)
            self.workspacemenu.menuAction().setEnabled(False)
            self.openAct.setEnabled(False)
            self.savefigAct.setEnabled(False)
            self.saveAct.setEnabled(False)
            self.newAct.setEnabled(False)
            self.closeAct.setEnabled(False)
            self.forwardAct.setEnabled(False)
            self.backAct.setEnabled(False)
            for i in range(self.tabs.count()):
                if i != self.workspaceNum:
                    self.tabs.setTabEnabled(i, False)
        self.undoAction.setEnabled(False)
        self.redoAction.setEnabled(False)

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
            return
        while (givenName in self.workspaceNames) or givenName == '':
            self.dispMsg("Workspace name '" + givenName + "' already exists")
            givenName, ok = QtWidgets.QInputDialog.getText(self, message, 'Name:', text=name)
            if not ok:
                return
        return givenName

    def undo(self, *args):
        if self.mainWindow is not None:
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
        while (givenName in self.macros.keys()) or givenName is '':
            self.dispMsg("Macro name '" + givenName + "' already exists")
            givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)
            if not ok:
                return
        self.macros[givenName] = []
        self.mainWindow.redoMacro = []
        self.mainWindow.currentMacro = givenName
        action1 = self.macrolistmenu.addAction(givenName, lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(givenName, lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(givenName, lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(givenName, lambda name=givenName: self.renameMacro(name))
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
        while (givenName in self.macros.keys()) or givenName is '':
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
        action1 = self.macrolistmenu.addAction(givenName, lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(givenName, lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(givenName, lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(givenName, lambda name=givenName: self.renameMacro(name))
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
        import json
        fileName = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.LastLocation + os.path.sep + name + '.json', 'JSON (*.json)')
        if type(fileName) is tuple:
            fileName = fileName[0]
        if fileName:  # if not cancelled
            self.LastLocation = os.path.dirname(fileName)
        if not fileName:
            return
        with open(fileName, 'w') as f:
            json.dump(self.macros[name], f, indent=4)

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
        import json
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.LastLocation)
        if type(filename) is tuple:
            filename = filename[0]
        if filename:  # if not cancelled
            self.LastLocation = os.path.dirname(filename)  # Save used path
        if len(filename) == 0:
            return
        self.stopMacro()
        count = 0
        name = 'macro' + str(count)
        while name in self.macros.keys():
            count += 1
            name = 'macro' + str(count)
        givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)
        while (givenName in self.macros.keys()) or givenName is '':
            if not ok:
                return
            self.dispMsg("Macro name '" + givenName + "' already exists")
            givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)
        with open(filename, 'r') as f:
            self.macros[givenName] = json.load(f)
        action1 = self.macrolistmenu.addAction(givenName, lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(givenName, lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(givenName, lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(givenName, lambda name=givenName: self.renameMacro(name))
        self.macroActions[givenName] = [action1, action2, action3, action4]
        self.menuCheck()

    def referenceAdd(self, reffreq, name):
        self.referenceName.append(name)
        self.referenceValue.append(reffreq)  # List with saved refrence values
        action1 = self.referencerunmenu.addAction(name, lambda name=name: self.referenceRun(name))
        action2 = self.referencedeletemenu.addAction(name, lambda name=name: self.referenceRemove(name))
        action3 = self.referencerenamemenu.addAction(name, lambda name=name: self.referenceRename(name))
        action4 = self.referencesavemenu.addAction(name, lambda name=name: self.referenceSave(name))
        self.referenceActions[name] = [action1, action2, action3, action4]
        self.menuCheck()

    def referenceClear(self):
        self.mainWindow.undoList.append(self.mainWindow.current.setRef(None))

    def referenceRun(self, name):
        reffreq = self.referenceValue[self.referenceName.index(name)]
        self.mainWindow.undoList.append(self.mainWindow.current.setRef(reffreq))

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
        while (givenName in self.referenceName) or givenName is '':
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
#        self.macrorenamemenu.removeAction(oldActions[3])
        action1 = self.referencerunmenu.addAction(givenName, lambda name=givenName: self.referenceRun(name))
        action2 = self.referencedeletemenu.addAction(givenName, lambda name=givenName: self.referenceRemove(name))
        action3 = self.referencerenamemenu.addAction(givenName, lambda name=givenName: self.referenceRename(name))
        action4 = self.referencesavemenu.addAction(givenName, lambda name=givenName: self.referenceSave(name))

#        action4 = self.macrorenamemenu.addAction(givenName, lambda name=givenName: self.renameMacro(name))
        self.referenceActions[givenName] = [action1, action2, action3, action4]
        self.menuCheck()

    def referenceSave(self, name):
        fileName = QtWidgets.QFileDialog.getSaveFileName(self, 'Save reference', self.LastLocation + os.path.sep + name + '.txt', 'txt (*.json)')
        if type(fileName) is tuple:
            fileName = fileName[0]
        if not fileName:
            return
        else:
            self.LastLocation = os.path.dirname(fileName)
        reffreq = self.referenceValue[self.referenceName.index(name)]
        with open(fileName, 'w') as f:
            f.write(str(reffreq))

    def referenceLoad(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.LastLocation)
        if type(filename) is tuple:
            filename = filename[0]
        if filename:  # if not cancelled
            self.LastLocation = os.path.dirname(filename)  # Save used path
        if len(filename) == 0:
            return
        count = 0
        name = 'ref' + str(count)
        while name in self.referenceName:
            count += 1
            name = 'ref' + str(count)

        givenName, ok = QtWidgets.QInputDialog.getText(self, 'Reference name', 'Name:', text=name)

        while (givenName in self.macros.keys()) or givenName is '':
            if not ok:
                return
            self.dispMsg('Name exists')
            givenName, ok = QtWidgets.QInputDialog.getText(self, 'Macro name', 'Name:', text=name)

        with open(filename, 'r') as f:
            self.referenceName.append(givenName)
            try:
                freq = float(f.read())

            except:
                self.dispMsg("Failed loading '" + filename + "' as reference.")
                return

        self.referenceValue.append(freq)
        action1 = self.referencerunmenu.addAction(givenName, lambda name=givenName: self.referenceRun(name))
        action2 = self.referencedeletemenu.addAction(givenName, lambda name=givenName: self.referenceRemove(name))
        action3 = self.referencerenamemenu.addAction(givenName, lambda name=givenName: self.referenceRename(name))
        action4 = self.referencesavemenu.addAction(givenName, lambda name=givenName: self.referenceSave(name))
        self.referenceActions[givenName] = [action1, action2, action3, action4]
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
        self.updWorkspaceMenu(var)
        self.menuCheck()
        if isinstance(self.mainWindow.current, (sc.CurrentMulti)):
            self.mainWindow.sideframe.checkChanged()

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
            self.updWorkspaceMenu(self.workspaceNames[self.workspaceNum])
            self.menuCheck()
            if isinstance(self.mainWindow.current, (sc.CurrentMulti)):
                self.mainWindow.sideframe.checkChanged()

    def duplicateWorkspace(self, *args):
        name = self.askName()
        if name is None:
            return
        self.workspaces.append(Main1DWindow(self, copy.deepcopy(self.mainWindow.get_masterData()), self.mainWindow.get_current()))
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
        self.updWorkspaceMenu(name)
        self.workspaces[self.workspaceNum].rename(name)

    def destroyWorkspace(self, num=None):
        if self.mainWindow is None:
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
            self.workspaceNum = 0
        if num < self.workspaceNum:
            self.workspaceNum -= 1
        if len(self.workspaces) > 0:
            self.changeMainWindow(self.workspaceNames[self.workspaceNum])
        else:
            self.logo.show()
            self.tabs.hide()
            self.updWorkspaceMenu(None)

    def updWorkspaceMenu(self, var):
        self.activemenu.clear()
        for i in self.workspaceNames:
            self.activemenu.addAction(i, lambda i=i: self.changeMainWindow(i))
        self.menuCheck()

    def newWorkspace(self, masterData):
        name = self.askName()
        if name is None:
            return
        self.workspaces.append(Main1DWindow(self, masterData))
        self.workspaces[-1].rename(name)
        self.tabs.addTab(self.workspaces[-1], name)
        self.workspaceNames.append(name)
        self.changeMainWindow(name)
        return 1

    def createCombineWorkspaceWindow(self):
        CombineWorkspaceWindow(self)

    def combineWorkspace(self, combineNames):
        wsname = self.askName()
        if wsname is None:
            return
        i = self.workspaceNames.index(combineNames[0])
        combineMasterData = copy.deepcopy(self.workspaces[i].get_masterData())
        shapeRequired = combineMasterData.data.shape
        combineMasterData.split(1, -1)
        for name in combineNames[1:]:
            i = self.workspaceNames.index(name)
            addData = self.workspaces[i].get_masterData()
            if addData.data.shape != shapeRequired:
                self.dispMsg("Not all the data has the required shape")
                return False
            combineMasterData.insert(addData.data, combineMasterData.data.shape[0], 0)
        self.workspaces.append(Main1DWindow(self, combineMasterData))
        self.workspaces[-1].rename(wsname)
        self.tabs.addTab(self.workspaces[-1], wsname)
        self.workspaceNames.append(wsname)
        self.changeMainWindow(wsname)
        return True

    def loadFromMenu(self):
        fileList = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open File', self.LastLocation)
        if type(fileList) is tuple:
            fileList = fileList[0]
        for filePath in fileList:
            if filePath:  # if not cancelled
                self.LastLocation = os.path.dirname(filePath)  # Save used path
            if len(filePath) == 0:
                return
            if filePath.endswith('.zip'):
                import tempfile
                import shutil
                import zipfile
                try:
                    temp_dir = tempfile.mkdtemp()
                    zipfile.ZipFile(filePath).extractall(temp_dir)
                    for i in os.listdir(temp_dir):
                        if self.autoLoad(os.path.join(temp_dir, i),realpath=filePath):
                            break
                finally:
                    shutil.rmtree(temp_dir)
            else:
                self.autoLoad(filePath)

    def autoLoad(self, filePath,realpath=False):
        returnVal = 0
        if os.path.isfile(filePath):
            filename = os.path.basename(filePath)
            if filename.endswith('.fid') or filename.endswith('.spe'):
                with open(filePath, 'r') as f:
                    check = int(np.fromfile(f, np.float32, 1))
                if check == 0:
                    self.loading(8, filePath)  # Suspected NMRpipe format
                else:  # SIMPSON
                    self.loading(4, filePath)
                return returnVal
            elif filename.endswith('.json') or filename.endswith('.JSON'):
                self.loading(5, filePath)
                return returnVal
            elif filename.endswith('.mat') or filename.endswith('.MAT'):
                self.loading(6, filePath)
                return returnVal
            elif filename.endswith('.jdf'):#JEOL delta format
                self.loading(9, filePath)
                return returnVal
            filePath = os.path.dirname(filePath)
            returnVal = 1
        direc = filePath
        if os.path.exists(direc + os.path.sep + 'procpar') and os.path.exists(direc + os.path.sep + 'fid'):
            self.loading(0, filePath)
            return returnVal
            # And for varian processed data
        if (os.path.exists(direc + os.path.sep + '..' + os.path.sep + 'procpar') or os.path.exists(direc + os.path.sep + 'procpar')) and os.path.exists(direc + os.path.sep + 'data'):
            self.loading(0, filePath)
            return returnVal
        elif os.path.exists(direc + os.path.sep + 'acqus') and (os.path.exists(direc + os.path.sep + 'fid') or os.path.exists(direc + os.path.sep + 'ser')):
            self.loading(1, filePath)
            return returnVal
        elif os.path.exists(direc + os.path.sep + 'procs') and (os.path.exists(direc + os.path.sep + '1r') or os.path.exists(direc + os.path.sep + '2rr')):
            self.loading(7, filePath)
            return returnVal
        elif os.path.exists(direc + os.path.sep + 'acq') and os.path.exists(direc + os.path.sep + 'data'):
            self.loading(2, filePath)
            return returnVal
        elif os.path.exists(direc + os.path.sep + 'acqu.par'):
            dirFiles = os.listdir(direc)
            files2D = [x for x in dirFiles if '.2d' in x]
            files1D = [x for x in dirFiles if '.1d' in x]
            if len(files2D) != 0 or len(files1D) != 0:
                self.loading(3, filePath,realpath=realpath)
                return returnVal

    def dataFromFit(self, data, filePath, freq, sw, spec, wholeEcho, ref, xaxArray, axes):
        name = self.askName()
        if name is None:
            return
        masterData = sc.Spectrum(name,
                                 data,
                                 filePath,
                                 freq,
                                 sw,
                                 spec,
                                 wholeEcho,
                                 ref,
                                 xaxArray,
                                 msgHandler=lambda msg: self.dispMsg(msg),
                                 history=['Data obtained from fit'])
        masterData.resetXax(axes)
        self.workspaces.append(Main1DWindow(self, masterData))
        self.tabs.addTab(self.workspaces[-1], name)
        self.workspaceNames.append(name)
        self.changeMainWindow(name)

    def loading(self, num, filePath, returnBool=False,realpath=False):
        if returnBool:
            name = None
        else:
            if realpath: #If there is a temp file, use the real path for name
                name = os.path.splitext(os.path.basename(realpath))[0]
            else:
                name = os.path.splitext(os.path.basename(filePath))[0]
            if self.defaultAskName:
                if realpath: #If there is a temperary directory
                    name = self.askName(realpath, name)
                else:
                    name = self.askName(filePath, name)
                if name is None:
                    return
            else:
                count = 0
                while name in self.workspaceNames:
                    name = 'spectrum' + str(count)
                    count += 1
        if num == 0:
            masterData = self.LoadVarianFile(filePath, name)
        elif num == 1:
            masterData = self.LoadBrukerTopspin(filePath, name)
        elif num == 2:
            masterData = self.LoadChemFile(filePath, name)
        elif num == 3:
            masterData = self.LoadMagritek(filePath, name,realpath)
        elif num == 4:
            masterData = self.LoadSimpsonFile(filePath, name)
        elif num == 5:
            masterData = self.loadJSONFile(filePath, name)
        elif num == 6:
            masterData = self.loadMatlabFile(filePath, name)
        elif num == 7:
            masterData = self.LoadBrukerSpectrum(filePath, name)
        elif num == 8:
            masterData = self.LoadPipe(filePath, name)
        elif num == 9:
            masterData = self.LoadJEOLDelta(filePath, name)
        if returnBool:
            return masterData
        else:
            if masterData is not None:
                self.workspaces.append(Main1DWindow(self, masterData))
                self.tabs.addTab(self.workspaces[-1], name)
                self.workspaceNames.append(name)
                self.changeMainWindow(name)

    def LoadVarianFile(self, filePath, name=''):
        from struct import unpack
        if os.path.isfile(filePath):
            Dir = os.path.dirname(filePath)
        else:
            Dir = filePath
        freq = 300e6
        sw = 50e3
        sw1 = 50e3
        freq1 = 0  # intialize second dimension freqency as 0

        if os.path.exists(Dir + os.path.sep + 'procpar'):
            file = Dir + os.path.sep + 'procpar'
        elif os.path.exists(Dir + os.path.sep + '..' + os.path.sep + 'procpar'):
            file = Dir + os.path.sep + '..' + os.path.sep + 'procpar'
        else:
            self.dispMsg(Dir + os.path.sep + 'procpar does not exits, used standard sw and freq')
            file = None
        if file is not None:
            with open(file, 'r') as f:
                data = f.read().split('\n')
            for s in range(0, len(data)):
                if data[s].startswith('sfrq '):
                    freq = float(data[s + 1].split()[1]) * 1e6
                elif data[s].startswith('sw '):
                    sw = float(data[s + 1].split()[1])
                elif data[s].startswith('sw1 '):
                    sw1 = float(data[s + 1].split()[1])
                elif data[s].startswith('dfrq '):
                    freq1 = float(data[s + 1].split()[1]) * 1e6
                elif data[s].startswith('reffrq '):
                    reffreq = float(data[s + 1].split()[1]) * 1e6
                elif data[s].startswith('reffrq1 '):
                    reffreq1 = float(data[s + 1].split()[1]) * 1e6
                elif data[s].startswith('phfid '):
                    if int(data[s].split()[-2]): #if on
                        phfid = float(data[s + 1].split()[1]) 
                    else:
                        phfid = 0
                elif data[s].startswith('rp '):
                    rp = float(data[s + 1].split()[1])  
                    
                

        if os.path.exists(Dir + os.path.sep + 'fid'):
            datafile = Dir + os.path.sep + 'fid'
            filePath = datafile
        elif os.path.exists(Dir + os.path.sep + 'data'):
            datafile = Dir + os.path.sep + 'data'
            filePath = datafile
        else:
            self.dispMsg('No valid data file found')
            return None
        with open(datafile, "rb") as f:
            raw = np.fromfile(f, np.int32, 6)
            nblocks = unpack('>l', raw[0])[0]
            ntraces = unpack('>l', raw[1])[0]
            npoints = unpack('>l', raw[2])[0]
            ebytes = unpack('>l', raw[3])[0]
            tbytes = unpack('>l', raw[4])[0]
            bbytes = unpack('>l', raw[5])[0]
            raw = np.fromfile(f, np.int16, 2)
            vers_id = unpack('>h', raw[0])[0]
            status = unpack('>h', raw[1])[0]
            spec = bool(int(bin(status)[-2]))
            raw = np.fromfile(f, np.int32, 1)
            nbheaders = unpack('>l', raw[0])[0]
            SizeTD2 = npoints
            SizeTD1 = nblocks * ntraces
            a = []
            fid32 = int(bin(status)[-3])
            fidfloat = int(bin(status)[-4])
            hypercomplex = bool(bin(status)[-5])
            flipped = bool(bin(status)[-10])

            if not fid32 and fidfloat:  # only for `newest' format, use fast routine
            
                totalpoints = (ntraces * npoints + nbheaders**2 * 7)*nblocks
                raw = np.fromfile(f, np.float32, totalpoints)
                a = raw.newbyteorder('>f')
#                print(bin(status)[-10])
                if not spec or (spec and not hypercomplex):
                    a = a.reshape(nblocks, int(totalpoints / nblocks))
                    a = a[:, 7::]
                    fid = a[:, ::2] - 1j * a[:, 1::2]
                else:
                    fid = a[nbheaders*7::4] - 1j * a[nbheaders*7+1::4]
                    fid = fid.reshape(int(SizeTD1/4),SizeTD2)
                    if flipped:
                        fid = np.fliplr(fid)
 
            elif fid32 and not fidfloat:  # for VNMRJ 2 data
                totalpoints = (ntraces * npoints + nbheaders**2 * 7)*nblocks
                raw = np.fromfile(f, np.int32, totalpoints)
                a = raw.newbyteorder('>l')
                a = a.reshape(nblocks, int(totalpoints / nblocks))
                a = a[:, 7::]
                fid = a[:, ::2] - 1j * a[:, 1::2]
            else:  # use slow, but robust routine
                for iter1 in range(0, nblocks):
                    b = []
                    for iter2 in range(0, nbheaders):
                        raw = np.fromfile(f, np.int16, nbheaders * 14)
                    if not fid32 and not fidfloat:
                        raw = np.fromfile(f, np.int16, ntraces * npoints)
                        for iter3 in raw:
                            b.append(unpack('>h', iter3)[0])
                    elif fid32 and not fidfloat:
                        raw = np.fromfile(f, np.int32, ntraces * npoints)
                        for iter3 in raw:
                            b.append(unpack('>l', iter3)[0])
                    else:
                        raw = np.fromfile(f, np.float32, ntraces * npoints)
                        for iter3 in raw:
                            b.append(unpack('>f', iter3)[0])
                    b = np.array(b)
                    if(len(b) != ntraces * npoints):
                        b.append(np.zeros(ntraces * npoints - len(b)))
                    a.append(b)
                a = np.complex128(a)
                fid = a[:, ::2] - 1j * a[:, 1::2]
        
        
        fid = fid * np.exp((rp + phfid) / 180 * np.pi * 1j) #apply zero order phase
        if SizeTD1 is 1:
            fid = fid[0][:]
            if spec:  # flip if spectrum
                fid = np.flipud(fid)
            masterData = sc.Spectrum(name, fid, (0, filePath), [freq], [sw], [bool(int(spec))],ref = [reffreq], msgHandler=lambda msg: self.dispMsg(msg))
        else:
            masterData = sc.Spectrum(name, fid, (0, filePath), [freq1, freq], [sw1, sw], [bool(int(spec))] * 2,ref = [reffreq1,reffreq], msgHandler=lambda msg: self.dispMsg(msg))
        masterData.addHistory("Varian data loaded from " + filePath)
        return masterData

    def LoadPipe(self, filePath, name=''):
        with open(filePath, 'r') as f:
            header = np.fromfile(f, np.float32, 512)

            NumberofPoints = int(header[99])
            data = np.fromfile(f, np.float32, NumberofPoints)
            if int(header[106]) == 0:  # if complex
                data = data + 1j * np.fromfile(f, np.float32, NumberofPoints)

            spec = int(header[220])  # 1 if ft, 0 if time
            freq = header[119] * 1e6
            sw = header[100]
            reference = header[101]  # frequency of last point in Hz

        sidefreq = -np.floor(NumberofPoints / 2) / NumberofPoints * sw  # freqeuency of last point on axis
        ref = sidefreq + freq - reference
        if spec == 1:
            data = np.flipud(data)

        masterData = sc.Spectrum(name, data, (8, filePath), [freq], [sw], [spec], ref=[ref], msgHandler=lambda msg: self.dispMsg(msg))
        masterData.addHistory("NMR pipe data loaded from " + filePath)
        return masterData

    def LoadJEOLDelta(self, filePath, name=''):
        from struct import unpack
        with open(filePath, "rb") as f:
            file_identifier = f.read(8)
            endian = unpack('>B',f.read(1))[0]
            f.read(3) #placeholder to get rid of unused data
            data_dimension_number = unpack('>B',f.read(1))[0]
            data_dimension_exist = unpack('>B',f.read(1))[0]
            data_type = unpack('>B',f.read(1))[0]
            f.read(1) #placeholder to get rid of unused data
            translate = np.fromfile(f, '>B', 8)
            data_axis_type = np.fromfile(f, '>B', 8)
            data_units = np.fromfile(f, '>B', 16).reshape(8,2)#Reshape for later unit extraction
            title = f.read(76)
            f.read(52)
            #176
            data_points = np.fromfile(f, '>I', 8)
            data_offset_start = np.fromfile(f, '>I', 8)
            data_offset_stop = np.fromfile(f, '>I', 8)
            data_axis_start = np.fromfile(f, '>d', 8)
            data_axis_stop = np.fromfile(f, '>d', 8)
            #400
            f.read(664)
            base_freq = np.fromfile(f, '>d', 8)
            #1128
            zero_point = np.fromfile(f, '>d', 8)
            reverse = np.fromfile(f, '>B', 8)
            #1200
            f.read(12)
            param_start = np.fromfile(f, '>I', 1)[0]
            param_length = np.fromfile(f, '>I', 1)[0]
            #1220
            f.read(64)
            data_start = np.fromfile(f, '>I', 1)[0]
            data_length = np.fromfile(f, '>Q', 1)[0]
            #1296
            f.read(data_start-1296) #skip to data_start
            #start reading the data
            if endian:
                dataendian = '<d'
            else:
                dataendian = '>d'
            if data_dimension_number == 1:
                if data_axis_type[0] == 1: #if real
                    data = np.fromfile(f, dataendian, data_points[0])
                elif data_axis_type[0] == 3: #if complex
                    data = np.fromfile(f, dataendian, data_points[0]) - 1j*np.fromfile(f, dataendian, data_points[0])
                    data = data[0:data_offset_stop[0]+1]
            elif data_dimension_number == 2:
                if data_axis_type[0] == 4: #if real-complex (no hypercomplex)
                    Step = 4 #Step size of block
                    pointsD2 = data_points[0]
                    pointsD1 = data_points[1]
                    datalength = pointsD2 * pointsD1
                    datareal = np.fromfile(f, dataendian, datalength)
                    dataimag = np.fromfile(f, dataendian, datalength)
                    data = datareal - 1j*dataimag
                    data = np.reshape(data,[pointsD1/Step,datalength/pointsD1/Step,Step,Step])
                    data = np.concatenate(np.concatenate(data,1),1)
                    data = data[0:data_offset_stop[1]+1,0:data_offset_stop[0]+1] #cut back to real size
                if data_axis_type[0] == 3: #if complex (i.e. hypercomplex)
                    Step = 32 #Step size of block
                    pointsD2 = data_points[0]
                    pointsD1 = data_points[1]
                    datalength = pointsD2 * pointsD1
                    datareal1 = np.fromfile(f, dataendian, datalength)
                    dataimag1 = np.fromfile(f, dataendian, datalength)
                    datareal2 = np.fromfile(f, dataendian, datalength)
                    dataimag2 = np.fromfile(f, dataendian, datalength)
                    data1 = datareal1 - 1j*dataimag1
                    data1 = np.reshape(data1,[pointsD1/Step,datalength/pointsD1/Step,Step,Step])
                    data1 = np.concatenate(np.concatenate(data1,1),1)
                    data2 = datareal2 - 1j*dataimag2
                    data2 = np.reshape(data2,[pointsD1/Step,datalength/pointsD1/Step,Step,Step])
                    data2 = np.concatenate(np.concatenate(data2,1),1)
                    data = np.zeros([data2.shape[0]*2,data2.shape[1]],dtype=complex)
                    if reverse[0] == 0:
                        data[::2,:] = data1 #Interleave both types in D1
                        data[1::2,:] = data2
                    else:
                        data[::2,:] = data2 #Interleave both types in D1
                        data[1::2,:] = data1           
                    data = data[0:(data_offset_stop[1] + 1)*2,0:data_offset_stop[0]+1] #Cut back to real size
                    if reverse[1] == 1:
                        data = np.conjugate(data)
            else:
                    self.dispMsg('No valid data file found')
                    return None
	#unitExp = data_units[0][0] &15#unit_exp: if (unitExp > 7) >unitExp -= 16;
	#scaleType = (data_units[0][1]>>4) &15 #scaleType: if (scaleType > 7) scaleType -= 16;
        sw = []
        spec = []
        ref = []
        freq = []
        for axisNum in reversed(range(data_dimension_number)):
            freq.append(base_freq[0]*1e6)
            axisType = data_units[axisNum][1] #Sec = 28, Hz = 13, PPM = 26
            axisScale = data_units[axisNum][0]
            if axisType == 28: #Sec
                scale = (axisScale >> 4) & 15
                if scale > 7:
                    scale = scale -16
                dw = (data_axis_stop[axisNum]-data_axis_start[axisNum])/(data_offset_stop[axisNum]+1 -1)*10**(-scale*3)#minus one to give same axis as spectrum???
                #scale for SI prefix
                sw.append(1.0 / dw)
                spec.append(False)
                sidefreq = -np.floor((data_offset_stop[axisNum]+1) / 2) / data_offset_stop[axisNum]+1 * sw[-1]  # frequency of last point on axis
                ref.append(base_freq[axisNum]*1e6)
            if axisType == 13: #Hz
                sw.append(np.abs(data_axis_start[axisNum]-data_axis_stop[0]))
                spec.append(True)
                if data_dimension_number == 1:
                    data = np.flipud(data)
                sidefreq = -np.floor(data_points[axisNum] / 2) / (data_offset_stop[axisNum]+1) * sw[-1]  # frequency of last point on axis
                ref.append(sidefreq + base_freq[axisNum]*1e6 - data_axis_stop[axisNum])
            if axisType == 26: #ppm
                sw.append(np.abs(data_axis_start[axisNum]-data_axis_stop[axisNum])*base_freq[axisNum])
                spec.append(True)
                if data_dimension_number == 1:
                    data = np.flipud(data)
                sidefreq = -np.floor((data_offset_stop[axisNum]+1) / 2) / (data_offset_stop[axisNum]+1) * sw[-1]  # frequency of last point on axis
                ref.append(sidefreq + base_freq[axisNum]*1e6 - data_axis_stop[axisNum]*base_freq[axisNum])
        masterData = sc.Spectrum(name, data, (9, filePath), freq, sw, spec,ref=ref, msgHandler=lambda msg: self.dispMsg(msg))
        masterData.addHistory("JEOL delta data loaded from " + filePath)
        return masterData


    def loadJSONFile(self, filePath, name=''):
        import json
        with open(filePath, 'r') as inputfile:
            struct = json.load(inputfile)
        data = np.array(struct['dataReal']) + 1j * np.array(struct['dataImag'])
        ref = np.where(np.isnan(struct['ref']), None, struct['ref'])
        if 'history' in struct.keys():
            history = struct['history']
        else:
            history = None
        xaxA = []
        for i in struct['xaxArray']:
            xaxA.append(np.array(i))
        masterData = sc.Spectrum(name,
                                 data,
                                 (5, filePath),
                                 list(struct['freq']),
                                 list(struct['sw']),
                                 list(struct['spec']),
                                 list(np.array(struct['wholeEcho'], dtype=bool)),
                                 list(ref),
                                 xaxA,
                                 msgHandler=lambda msg: self.dispMsg(msg),
                                 history=history)
        masterData.addHistory("JSON data loaded from " + filePath)
        return masterData

    def loadMatlabFile(self, filePath, name=''):
        import scipy.io
        with open(filePath, 'rb') as inputfile:  # read first several bytes the check .mat version
            teststring = inputfile.read(13)
        version = float(teststring.decode("utf-8")[7:10])  # extract version from the binary array
        if version < 7.3:  # all versions below 7.3 are supported
            matlabStruct = scipy.io.loadmat(filePath)
            var = [k for k in matlabStruct.keys() if not k.startswith('__')][0]
            mat = matlabStruct[var]
            if mat['dim'] == 1:
                data = mat['data'][0, 0][0]
                xaxA = [k[0] for k in (mat['xaxArray'][0])]
            else:
                data = mat['data'][0, 0]
                if all(x == data.shape[0] for x in data.shape):
                    xaxA = [k for k in (mat['xaxArray'][0, 0])]
                else:
                    xaxA = [k[0] for k in (mat['xaxArray'][0, 0][0])]
            # insert some checks for data type
            ref = mat['ref'][0, 0][0]
            ref = np.where(np.isnan(ref), None, ref)
            if 'history' in mat.dtype.names:
                history = list(np.array(mat['history'][0, 0], dtype=str))
            else:
                history = None
            masterData = sc.Spectrum(name,
                                     data,
                                     (6, filePath),
                                     list(mat['freq'][0, 0][0]),
                                     list(mat['sw'][0, 0][0]),
                                     list(mat['spec'][0, 0][0]),
                                     list(np.array(mat['wholeEcho'][0, 0][0]) > 0),
                                     list(ref),
                                     xaxA,
                                     msgHandler=lambda msg: self.dispMsg(msg),
                                     history=history)
            masterData.addHistory("Matlab data loaded from " + filePath)
            return masterData
        else:  # If the version is 7.3, use HDF5 type loading
            import h5py  # For .mat v 7.3 support
            f = h5py.File(filePath, 'r')
            Groups = []
            for name in f:
                if name != '#refs#':
                    Groups.append(name)
            DataGroup = Groups[0]  # get the groupo name
            mat = f[DataGroup]
            if np.array(mat['dim'])[0][0] == 1:
                xaxA = list([np.array(mat['xaxArray'])[:, 0]])
                data = np.array(mat['data'])
                data = (data['real'] + data['imag'] * 1j)[:, 0]  # split and use real and imag part
            else:
                data = np.transpose(np.array(mat['data']))
                data = data['real'] + data['imag'] * 1j
                if all(x == data.shape[0] for x in data.shape):
                    xaxA = [np.array(mat[k]) for k in (mat['xaxArray'])]
                else:
                    xaxA = [np.array(mat[k[0]]) for k in (mat['xaxArray'])]
            ref = np.array(mat['ref'])[:, 0]
            ref = np.where(np.isnan(ref), None, ref)
            if 'history' in mat.keys():
                history = list()
                history.append([item.astype(np.int8).tostring().decode("ascii") for item in np.array(mat['history']).transpose()])
                history = history[0]
            else:
                history = None
            masterData = sc.Spectrum(name,
                                     data,
                                     (6, filePath),
                                     list(np.array(mat['freq'])[:, 0]),
                                     list(np.array(mat['sw'])[:, 0]),
                                     list(np.array(mat['spec'])[:, 0]),
                                     list(np.array(mat['wholeEcho'])[:, 0] > 0),
                                     list(ref),
                                     xaxA,
                                     msgHandler=lambda msg: self.dispMsg(msg),
                                     history=history)
            masterData.addHistory("Matlab data loaded from " + filePath)
            return masterData

    def LoadBrukerTopspin(self, filePath, name=''):
        if os.path.isfile(filePath):
            Dir = os.path.dirname(filePath)
        else:
            Dir = filePath
        if os.path.exists(Dir + os.path.sep + 'acqus'):
            with open(Dir + os.path.sep + 'acqus', 'r') as f:
                data = f.read().split('\n')
            for s in range(0, len(data)):
                if data[s].startswith('##$TD='):
                    sizeTD2 = int(data[s][6:])
                if data[s].startswith('##$SFO1='):
                    freq2 = float(data[s][8:]) * 1e6
                if data[s].startswith('##$SW_h='):
                    SW2 = float(data[s][8:])
                if data[s].startswith('##$BYTORDA='):
                    ByteOrder = int(data[s][11:])
        sizeTD1 = 1
        if os.path.exists(Dir + os.path.sep + 'acqu2s'):
            with open(Dir + os.path.sep + 'acqu2s', 'r') as f:
                data2 = f.read().split('\n')
            for s in range(0, len(data2)):
                if data2[s].startswith('##$TD='):
                    sizeTD1 = int(data2[s][6:])
                if data2[s].startswith('##$SFO1='):
                    freq1 = float(data2[s][8:]) * 1e6
                if data2[s].startswith('##$SW_h='):
                    SW1 = float(data2[s][8:])
        if os.path.exists(Dir + os.path.sep + 'fid'):
            filePath = Dir + os.path.sep + 'fid'
            with open(Dir + os.path.sep + 'fid', "rb") as f:
                raw = np.fromfile(f, np.int32, sizeTD1* sizeTD2)
        elif os.path.exists(Dir + os.path.sep + 'ser'):
            filePath = Dir + os.path.sep + 'ser'
            with open(Dir + os.path.sep + 'ser', "rb") as f:
                raw = np.fromfile(f, np.int32, sizeTD1 * int(np.ceil(sizeTD2 / 256))*256) #Always load full 1024 byte blocks (256 data points)
        if ByteOrder:
            RawInt = raw.newbyteorder('b')
        else:
            RawInt = raw.newbyteorder('l')
        ComplexData = np.array(RawInt[0:len(RawInt):2]) + 1j * np.array(RawInt[1:len(RawInt):2])
        spec = [False]
        if sizeTD1 is 1:
            masterData = sc.Spectrum(name, ComplexData, (1, filePath), [freq2], [SW2], spec, msgHandler=lambda msg: self.dispMsg(msg))
        else:
            ComplexData = ComplexData.reshape(sizeTD1, int(np.ceil(sizeTD2 / 256) * 256 / 2))
            ComplexData = ComplexData[:,0:int(sizeTD2/2)] #Cut off placeholder data
            masterData = sc.Spectrum(name, ComplexData, (1, filePath), [freq1, freq2], [SW1, SW2], spec * 2, msgHandler=lambda msg: self.dispMsg(msg))
        masterData.addHistory("Bruker data loaded from " + filePath)
        return masterData

    def LoadBrukerSpectrum(self, filePath, name=''):
        if os.path.isfile(filePath):
            Dir = os.path.dirname(filePath)
        else:
            Dir = filePath
        if os.path.exists(Dir + os.path.sep + 'procs'):  # Get D2 parameters
            with open(Dir + os.path.sep + 'procs', 'r') as f:
                data = f.read().split('\n')
            for s in range(0, len(data)):
                if data[s].startswith('##$SI='):
                    sizeTD2 = int(re.findall("\#\#\$SI= (.*.)", data[s])[0])
#                if data[s].startswith('##$XDIM='):
#                    blockingD2 = int(data[s][8:])
                if data[s].startswith('##$BYTORDP='):
                    ByteOrder = int(data[s][11:])
                if data[s].startswith('##$SW_p='):
                    SW2 = float(data[s][8:])
                if data[s].startswith('##$SF='):
                    Ref2 = float(data[s][6:])*1e6 
                    
        freq2 = 0
        if os.path.exists(Dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'acqus'):  # Get D2 parameters from fid directory, if available
            with open(Dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'acqus', 'r') as f:
                data = f.read().split('\n')
            for s in range(0, len(data)):
                if data[s].startswith('##$SFO1='):
                    freq2 = float(data[s][8:]) * 1e6
        sizeTD1 = 1
        if os.path.exists(Dir + os.path.sep + 'proc2s'):  # Get D1 parameters
            with open(Dir + os.path.sep + 'proc2s', 'r') as f:
                data2 = f.read().split('\n')
            for s in range(0, len(data2)):
                if data2[s].startswith('##$SI='):
                    sizeTD1 = int(data2[s][6:])
#                if data2[s].startswith('##$XDIM='):
#                    blockingD1 = int(data[s][8:])
                if data2[s].startswith('##$SW_p='):
                    SW1 = float(data2[s][8:])
                if data2[s].startswith('##$SF='):
                    Ref1 = float(data2[s][6:])*1e6
        freq1 = 0
        if os.path.exists(Dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'acqu2s'):  # Get D1 parameters from fid directory, if available
            with open(Dir + os.path.sep + '..' + os.path.sep + '..' + os.path.sep + 'acqu2s', 'r') as f:
                data = f.read().split('\n')
            for s in range(0, len(data)):
                if data[s].startswith('##$SFO1='):
                    freq1 = float(data[s][8:]) * 1e6
        if os.path.exists(Dir + os.path.sep + '1r'):  # Get D2 data
            filePath = Dir + os.path.sep + '1r'
            with open(Dir + os.path.sep + '1r', "rb") as f:
                RawReal = np.fromfile(f, np.int32, sizeTD1 * sizeTD2)
            RawImag = np.zeros([sizeTD1 * sizeTD2])
            if os.path.exists(Dir + os.path.sep + '1i'):
                with open(Dir + os.path.sep + '1i', "rb") as f:
                    RawImag = np.fromfile(f, np.int32, sizeTD1 * sizeTD2)
        elif os.path.exists(Dir + os.path.sep + '2rr'):  # Get D1 data
            filePath = Dir + os.path.sep + '2rr'
            with open(Dir + os.path.sep + '2rr', "rb") as f:
                RawReal = np.fromfile(f, np.int32, sizeTD1 * sizeTD2)
            RawImag = np.zeros([sizeTD1 * sizeTD2])
            if os.path.exists(Dir + os.path.sep + '2ir'):  # If hypercomplex
                with open(Dir + os.path.sep + '2ir', "rb") as f:
                    RawImag = np.fromfile(f, np.int32, sizeTD1 * sizeTD2)
            elif os.path.exists(Dir + os.path.sep + '2ii'):
                with open(Dir + os.path.sep + '2ii', "rb") as f:
                    RawImag = np.fromfile(f, np.int32, sizeTD1 * sizeTD2)
        if ByteOrder:
            RawReal = RawReal.newbyteorder('b')
            RawImag = RawImag.newbyteorder('b')
        else:
            RawReal = RawReal.newbyteorder('l')
            RawImag = RawImag.newbyteorder('l')
        Data = np.flipud(RawReal) - 1j * np.flipud(RawImag)
        spec = [True]
        if sizeTD1 is 1:
            masterData = sc.Spectrum(name, Data, (7, filePath), [freq2], [SW2], spec,ref=[Ref2], msgHandler=lambda msg: self.dispMsg(msg))
        else:
            Data = Data.reshape(sizeTD1, sizeTD2)
            masterData = sc.Spectrum(name, Data, (7, filePath), [freq1, freq2], [SW1, SW2], spec * 2,ref=[Ref1,Ref2], msgHandler=lambda msg: self.dispMsg(msg))
        masterData.addHistory("Bruker spectrum data loaded from " + filePath)
        return masterData

    def LoadChemFile(self, filePath, name=''):
        if os.path.isfile(filePath):
            Dir = os.path.dirname(filePath)
        else:
            Dir = filePath
        sizeTD1 = 1
        sw1 = 50e3
        H = dict(line.strip().split('=') for line in open(Dir + os.path.sep + 'acq', 'r'))
        sizeTD2 = int(H['al'])
        freq = float(H['sf' + H['ch1']])
        sw = 1 / float(H['dw'][:-1])
        if any('array_num_values_' in s for s in H.keys()):
            if 'use_array=1' in open(Dir + '/acq_2').read():
                for s in H.keys():
                    if ('array_num_values_' in s):
                        sizeTD1 = sizeTD1 * int(H[s])
            else:
                if 'al2' in H:
                    sizeTD1 = int(float(H['al2']))
                    if 'dw2' in H:
                        sw1 = 1 / float(H['dw2'][:-1])
        else:
            if 'al2' in H:
                sizeTD1 = int(float(H['al2']))
                if 'dw2' in H:
                    sw1 = 1 / float(H['dw2'][:-1])
        with open(Dir + os.path.sep + 'data', 'rb') as f:
            raw = np.fromfile(f, np.int32)
            b = np.complex128(raw.byteswap())
        filePath = Dir + os.path.sep + 'data'
        fid = b[:len(b) / 2] + 1j * b[len(b) / 2:]
        fid = np.reshape(fid, (sizeTD1, sizeTD2))
        data = np.array(fid)
        spec = [False]
        if sizeTD1 is 1:
            data = data[0][:]
            masterData = sc.Spectrum(name, data, (2, filePath), [freq * 1e6], [sw], spec, msgHandler=lambda msg: self.dispMsg(msg))
        else:
            data = data.reshape((sizeTD1, sizeTD2))
            masterData = sc.Spectrum(name, data, (2, filePath), [freq * 1e6] * 2, [sw1, sw], spec * 2, msgHandler=lambda msg: self.dispMsg(msg))
        masterData.addHistory("Chemagnetics data loaded from " + filePath)
        return masterData

    def LoadMagritek(self, filePath, name='',realPath=''):
        # Magritek load script based on some Matlab files by Ole Brauckman
        if os.path.isfile(filePath):
            Dir = os.path.dirname(filePath)
        else:
            Dir = filePath
            
        if realPath:
            rememberPath = realPath
        else:
            rememberPath = filePath
        DirFiles = os.listdir(Dir)
        Files2D = [x for x in DirFiles if '.2d' in x]
        Files1D = [x for x in DirFiles if '.1d' in x]
        H = dict(line.strip().split('=') for line in open(Dir + os.path.sep + 'acqu.par', 'r'))
        sw = float(H['bandwidth                 ']) * 1000
        sizeTD2 = int(H['nrPnts                    '])
        freq = float(H['b1Freq                    '])*1e6
        lastfreq = float(H['lowestFrequency           '])
        
        sidefreq = -np.floor(sizeTD2 / 2) / sizeTD2 * sw  # freqeuency of last point on axis
        ref = sidefreq + freq - lastfreq
        if len(Files2D) == 1:
            File = Files2D[0]
            sizeTD1 = int(H['nrSteps                   '])
            if 'bandwidth2                ' in H:
                sw1 = float(H['bandwidth2                ']) *1000
            else:
                sw1 = 50e3
            if  'lowestFrequency2          ' in H:   
                lastfreq1 = float(H['lowestFrequency2          '])
            
                sidefreq1 = -np.floor(sizeTD1 / 2) / sizeTD1 * sw1  # freqeuency of last point on axis
                ref1 = sidefreq1 + freq - lastfreq1
            else:
                ref1 = None
            
            with open(Dir + os.path.sep + File, 'rb') as f:
                raw = np.fromfile(f, np.float32)
            Data = raw[-2 * sizeTD2 * sizeTD1::]
            ComplexData = Data[0:Data.shape[0]:2] - 1j * Data[1:Data.shape[0]:2]
            ComplexData = ComplexData.reshape((sizeTD1, sizeTD2))
            masterData = sc.Spectrum(name, ComplexData, (3, rememberPath), [freq] * 2, [sw1,sw], [False] * 2,ref=[ref1,ref], msgHandler=lambda msg: self.dispMsg(msg))
        elif len(Files1D) != 0:
            File = 'data.1d'
            with open(Dir + os.path.sep + File, 'rb') as f:
                raw = np.fromfile(f, np.float32)
            Data = raw[-2 * sizeTD2::]
            ComplexData = Data[0:Data.shape[0]:2] - 1j * Data[1:Data.shape[0]:2]
            masterData = sc.Spectrum(name, ComplexData, (3, rememberPath), [freq], [sw], [False],ref=[ref], msgHandler=lambda msg: self.dispMsg(msg))
        masterData.addHistory("Magritek data loaded from " + rememberPath)
        return masterData

    def LoadSimpsonFile(self, filePath, name=''):
        with open(filePath, 'r') as f:
            Lines = f.read().split('\n')
        NP, NI, SW, SW1, TYPE, FORMAT = 0, 1, 0, 0, '', 'Normal'
        DataStart = Lines.index('DATA')
        DataEnd = Lines.index('END')
        for s in range(0, DataStart):
            if Lines[s].startswith('NP='):
                NP = int(re.sub('NP=', '', Lines[s]))
            elif Lines[s].startswith('NI='):
                NI = int(re.sub('NI=', '', Lines[s]))
            elif Lines[s].startswith('SW='):
                SW = float(re.sub('SW=', '', Lines[s]))
            elif Lines[s].startswith('SW1='):
                SW1 = float(re.sub('SW1=', '', Lines[s]))
            elif Lines[s].startswith('TYPE='):
                TYPE = re.sub('TYPE=', '', Lines[s])
            elif Lines[s].startswith('FORMAT='):
                FORMAT = re.sub('FORMAT=', '', Lines[s])
        if 'Normal' in FORMAT:
            length = DataEnd - DataStart - 1
            data = np.zeros(length, dtype=complex)
            for i in range(length):
                temp = Lines[DataStart + 1 + i].split()
                data[i] = float(temp[0]) + 1j * float(temp[1])
        elif 'BINARY' in FORMAT:
            # Binary code based on:
            # pysimpson: Python module for reading SIMPSON files
            # By: Jonathan J. Helmus (jjhelmus@gmail.com)
            # Version: 0.1 (2012-04-13)
            # License: GPL
            chardata = ''.join(Lines[DataStart + 1:DataEnd])
            nquads, mod = divmod(len(chardata), 4)
            assert mod == 0     # character should be in blocks of 4
            BASE = 33
            charst =  np.fromstring(chardata, dtype=np.uint8)
            charst = charst.reshape(nquads,4) - BASE
            FIRST = lambda f, x: ((x) & ~(~0 << f))
            LAST = lambda f, x: ((x) & (~0 << (8 - f)))
            
            first = FIRST(6, charst[:,0]) | LAST(2, charst[:,1] << 2)
            second  = FIRST(4, charst[:,1]) | LAST(4, charst[:,2] << 2)
            third = FIRST(2, charst[:,2]) | LAST(6, charst[:,3] << 2)
            
            Bytes = np.ravel(np.transpose(np.array([first,second,third]))).astype('int64')
            
            # convert every 4 'bytes' to a float
            num_points, num_pad = divmod(len(Bytes), 4)
            Bytes = np.array(Bytes)
            Bytes=Bytes[:-num_pad]
            Bytes=Bytes.reshape(num_points,4)
            mantissa = ((Bytes[:,2] % 128) << 16) + (Bytes[:,1] << 8) + Bytes[:,0]
            exponent = (Bytes[:,3] % 128) * 2 + (Bytes[:,2] >= 128) * 1
            negative = Bytes[:,3] >= 128
            e = exponent - 127
            m = np.abs(mantissa) / np.float64(1 << 23)
            data = np.float32((-1)**negative*np.ldexp(m,e))
            data = data.view('complex64')
        if NI != 1:  # 2D data, reshape to NI, NP
            data = data.reshape(int(NI), -1)
        if 'FID' in TYPE:
            spec = [False]
        elif 'SPE' in TYPE:
            spec = [True]
        if NI is 1:
            masterData = sc.Spectrum(name, data, (4, filePath), [0], [SW], spec, msgHandler=lambda msg: self.dispMsg(msg))
        else:
            masterData = sc.Spectrum(name, data, (4, filePath), [0, 0], [SW1, SW], spec * 2, msgHandler=lambda msg: self.dispMsg(msg))
        masterData.addHistory("SIMPSON data loaded from " + filePath)
        return masterData

    def saveSimpsonFile(self):
        self.mainWindow.get_mainWindow().SaveSimpsonFile()

    def saveASCIIFile(self):
        self.mainWindow.get_mainWindow().saveASCIIFile()

    def saveJSONFile(self):
        self.mainWindow.get_mainWindow().saveJSONFile()

    def saveMatlabFile(self):
        self.mainWindow.get_mainWindow().saveMatlabFile()

    def saveFigure(self):
        if self.mainWindow is None:
            return
        self.allowChange = False
        self.menuDisable(True)
        num = self.workspaces.index(self.mainWindow)
        self.mainWindow = MainPlotWindow(self, self.mainWindow)
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
        self.menuEnable(True)
        self.tabs.setCurrentIndex(num)
        self.menuCheck()
        self.allowChange = True

    def createFitWindow(self, fitWindow):
        if self.mainWindow is None:
            return
        self.allowChange = False
        self.menuDisable(True)
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
        self.menuEnable(True)
        self.tabs.setCurrentIndex(num)
        self.menuCheck()
        self.allowChange = True

    def updateMenu(self):
        UpdateWindow(self)
        
    def createShiftConversionWindow(self):
        shiftConversionWindow(self)
        
    def createQuadConversionWindow(self):
        quadConversionWindow(self)

    def nmrTable(self):
        import subprocess
        subprocess.Popen([sys.executable, os.path.dirname(os.path.realpath(__file__)) + '/nmrTable.py'])
        
    def about(self):
        message = "ssNake " + VERSION
        QtWidgets.QMessageBox.about(self, 'About', message)

    def fileQuit(self):
        self.close()

    def closeEvent(self, event):
        quit_msg = "Are you sure you want to close ssNake?"
        reply = QtWidgets.QMessageBox.question(self, 'Close', quit_msg, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
        if reply == QtWidgets.QMessageBox.Yes:
            for item in fit.stopDict.keys(): #Send stop commands to all threads
                fit.stopDict[item] = True
            
            
            event.accept()
        else:
            event.ignore()

######################################################################################################


class Main1DWindow(QtWidgets.QWidget):

    def __init__(self, father, masterData, duplicateCurrent=None):
        QtWidgets.QWidget.__init__(self, father)
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        grid = QtWidgets.QGridLayout(self)
        grid.addWidget(self.canvas, 0, 0)
        self.undoList = []
        self.redoList = []
        self.currentMacro = None
        self.redoMacro = []
        self.monitor = None # Monitor of files
        self.monitorMacros = []
        self.father = father
        self.mainProgram = self.father  # remove all references to mainprogram to father
        self.masterData = masterData
        if duplicateCurrent is not None:
            self.current = duplicateCurrent.copyCurrent(self, self.fig, self.canvas, masterData)
        else:
            self.current = sc.Current1D(self, self.fig, self.canvas, masterData)
        self.menubar = self.father.menubar
        self.sideframe = SideFrame(self)
        grid.addWidget(self.sideframe, 0, 1)
        self.bottomframe = BottomFrame(self)
        grid.addWidget(self.bottomframe, 1, 0)
        self.textframe = TextFrame(self)
        grid.addWidget(self.textframe, 2, 0)
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

    def rescue(self):
        self.current.kill()
        self.current = sc.Current1D(self, self.current.fig, self.current.canvas, self.masterData)

    def menuEnable(self):
        self.father.menuEnable()
        self.sideframe.frameEnable()
        self.bottomframe.frameEnable()
        self.textframe.frameEnable()
        self.menuCheck()

    def menuDisable(self):
        self.father.menuDisable()
        self.sideframe.frameDisable()
        self.bottomframe.frameDisable()
        self.textframe.frameDisable()

    def menuCheck(self):
        self.father.menuCheck()

    def runMacro(self, macro, display=True):
        self.redoList = []
        for iter1 in macro:
            if iter1[0] == 'reload':
                loadData = self.father.loading(self.masterData.filePath[0], self.masterData.filePath[1], True) 
                self.undoList.append(self.masterData.restoreData(loadData, None))
            elif iter1[0] == 'real':
                self.undoList.append(self.masterData.real())
            elif iter1[0] == 'imag':
                self.undoList.append(self.masterData.imag())
            elif iter1[0] == 'abs':
                self.undoList.append(self.masterData.abs())
            elif iter1[0] == 'phase':
                self.undoList.append(self.masterData.setPhase(*iter1[1]))
            elif iter1[0] == 'fourier':
                self.undoList.append(self.masterData.fourier(*iter1[1]))
            elif iter1[0] == 'realFourier':
                self.undoList.append(self.masterData.realFourier(*iter1[1]))
            elif iter1[0] == 'fftshift':
                self.undoList.append(self.masterData.fftshift(*iter1[1]))
            elif iter1[0] == 'diff':
                self.undoList.append(self.masterData.diff(*iter1[1]))
            elif iter1[0] == 'cumsum':
                self.undoList.append(self.masterData.cumsum(*iter1[1]))
            elif iter1[0] == 'apodize':
                self.undoList.append(self.masterData.apodize(*iter1[1]))
            elif iter1[0] == 'freq':
                self.undoList.append(self.masterData.setFreq(*iter1[1]))
            elif iter1[0] == 'ref':
                self.undoList.append(self.masterData.setRef(*iter1[1]))
            elif iter1[0] == 'size':
                self.undoList.append(self.masterData.setSize(*iter1[1]))
            elif iter1[0] == 'spec':
                self.undoList.append(self.masterData.changeSpec(*iter1[1]))
            elif iter1[0] == 'swapecho':
                self.undoList.append(self.masterData.swapEcho(*iter1[1]))
            elif iter1[0] == 'wholeEcho':
                self.undoList.append(self.masterData.wholeEcho(*iter1[1]))
            elif iter1[0] == 'shift':
                self.undoList.append(self.masterData.shiftData(*iter1[1]))
            elif iter1[0] == 'states':
                self.undoList.append(self.masterData.states(*iter1[1]))
            elif iter1[0] == 'statesTPPI':
                self.undoList.append(self.masterData.statesTPPI(*iter1[1]))
            elif iter1[0] == 'echoAntiEcho':
                self.undoList.append(self.masterData.echoAntiEcho(*iter1[1]))
            elif iter1[0] == 'baselineCorrection':
                self.undoList.append(self.masterData.baselineCorrection(*iter1[1]))
            elif iter1[0] == 'integrate':
                self.undoList.append(self.masterData.matrixManip(*iter1[1], which=0))
            elif iter1[0] == 'sum':
                self.undoList.append(self.masterData.matrixManip(*iter1[1], which=5))
            elif iter1[0] == 'max':
                self.undoList.append(self.masterData.matrixManip(*iter1[1], which=1))
            elif iter1[0] == 'min':
                self.undoList.append(self.masterData.matrixManip(*iter1[1], which=2))
            elif iter1[0] == 'argmax':
                self.undoList.append(self.masterData.matrixManip(*iter1[1], which=3))
            elif iter1[0] == 'argmin':
                self.undoList.append(self.masterData.matrixManip(*iter1[1], which=4))
            elif iter1[0] == 'average':
                self.undoList.append(self.masterData.matrixManip(*iter1[1], which=6))
            elif iter1[0] == 'fliplr':
                self.undoList.append(self.masterData.flipLR(*iter1[1]))
            elif iter1[0] == 'concatenate':
                self.undoList.append(self.masterData.concatenate(*iter1[1]))
            elif iter1[0] == 'split':
                self.undoList.append(self.masterData.split(*iter1[1]))
            elif iter1[0] == 'insert':
                self.undoList.append(self.masterData.insert(*iter1[1]))
            elif iter1[0] == 'delete':
                self.undoList.append(self.masterData.remove(*iter1[1]))
            elif iter1[0] == 'add':
                self.undoList.append(self.masterData.add(*iter1[1]))
            elif iter1[0] == 'subtract':
                self.undoList.append(self.masterData.subtract(*iter1[1]))
            elif iter1[0] == 'multiply':
                self.undoList.append(self.masterData.multiply(*iter1[1]))
            elif iter1[0] == 'subtractAvg':
                self.undoList.append(self.masterData.subtractAvg(*iter1[1]))
            elif iter1[0] == 'FIDDLE':
                self.undoList.append(self.masterData.fiddle(*iter1[1]))
            elif iter1[0] == 'reorder':
                self.undoList.append(self.masterData.reorder(*iter1[1]))
            elif iter1[0] == 'ffm':
                self.undoList.append(self.masterData.ffm_1d(*iter1[1]))
            elif iter1[0] == 'clean':
                self.undoList.append(self.masterData.clean(*iter1[1]))
            elif iter1[0] == 'shear':
                self.undoList.append(self.masterData.shear(*iter1[1]))
            elif iter1[0] == 'extract':
                self.undoList.append(self.masterData.getRegion(*iter1[1]))
            elif iter1[0] == 'setxax':
                self.undoList.append(self.masterData.setXax(*iter1[1]))
            elif iter1[0] == 'hilbert':
                self.undoList.append(self.masterData.hilbert(*iter1[1]))
            else:
                self.father.dispMsg('unknown macro command: ' + iter1[0])
        if display:
            self.current.upd()  # get the first slice of data
            self.current.plotReset()  # reset the axes limits
            self.current.showFid()  # plot the data
            self.updAllFrames()
            self.menuCheck()

    def addMacro(self, macroStep):
        if self.currentMacro is not None:
            self.father.macroAdd(self.currentMacro, macroStep)
            self.redoMacro = []

    def saveJSONFile(self):
        import json
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum]  # Set name of file to be saved to workspace name to start
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.LastLocation + os.path.sep + WorkspaceName + '.json', 'JSON (*.json)')
        if type(name) is tuple:
            name = name[0]        
        if not name:
            return
        self.father.LastLocation = os.path.dirname(name)  # Save used path
        struct = {}
        struct['dataReal'] = np.real(self.masterData.data).tolist()
        struct['dataImag'] = np.imag(self.masterData.data).tolist()
        struct['freq'] = self.masterData.freq.tolist()
        struct['sw'] = list(self.masterData.sw)
        struct['spec'] = list(1.0 * np.array(self.masterData.spec))
        struct['wholeEcho'] = list(1.0 * np.array(self.masterData.wholeEcho))
        struct['ref'] = np.array(self.masterData.ref, dtype=np.float).tolist()
        struct['history'] = self.masterData.history
        tmpXax = []
        for i in self.masterData.xaxArray:
            tmpXax.append(i.tolist())
        struct['xaxArray'] = tmpXax
        with open(name, 'w') as outfile:
            json.dump(struct, outfile)

    def saveMatlabFile(self):
        import scipy.io
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum]  # Set name of file to be saved to workspace name to start
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.LastLocation + os.path.sep + WorkspaceName + '.mat', 'MATLAB file (*.mat)')
        if type(name) is tuple:
            name = name[0]        
        if not name:
            return
        self.father.LastLocation = os.path.dirname(name)  # Save used path
        struct = {}
        struct['dim'] = self.masterData.data.ndim
        struct['data'] = self.masterData.data
        struct['freq'] = self.masterData.freq
        struct['sw'] = self.masterData.sw
        struct['spec'] = self.masterData.spec
        struct['wholeEcho'] = self.masterData.wholeEcho
        struct['ref'] = np.array(self.masterData.ref, dtype=np.float)
        struct['history'] = self.masterData.history
        struct['xaxArray'] = self.masterData.xaxArray
        matlabStruct = {self.father.workspaceNames[self.father.workspaceNum]: struct}
        scipy.io.savemat(name, matlabStruct)

    def SaveSimpsonFile(self):
        if self.masterData.data.ndim > 2:
            self.father.dispMsg('Saving to Simpson format only allowed for 1D and 2D data!')
            return
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum]  # Set name of file to be saved to workspace name to start
        if sum(self.masterData.spec) / len(self.masterData.spec) == 1:
            name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.LastLocation + os.path.sep + WorkspaceName + '.spe', 'SIMPSON file (*.spe)')
            if type(name) is tuple:
                name = name[0]       
            if not name:
                return
        elif sum(self.masterData.spec) == 0:
            name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.LastLocation + os.path.sep + WorkspaceName + '.fid', 'SIMPSON file (*.fid)')
            if type(name) is tuple:
                name = name[0]        
            if not name:
                return
        else:
            self.father.dispMsg('Saving to Simpson format not allowed for mixed time/frequency domain data!')
            return
        self.father.LastLocation = os.path.dirname(name)  # Save used path
        with open(name, 'w') as f:
            f.write('SIMP\n')
            if self.masterData.data.ndim is 2:
                f.write('NP=' + str(self.masterData.data.shape[1]) + '\n')
                f.write('NI=' + str(self.masterData.data.shape[0]) + '\n')
                f.write('SW=' + str(self.masterData.sw[1]) + '\n')
                f.write('SW1=' + str(self.masterData.sw[0]) + '\n')
            else:
                f.write('NP=' + str(self.masterData.data.shape[0]) + '\n')
                f.write('SW=' + str(self.masterData.sw[0]) + '\n')
            if self.masterData.spec[0]:
                f.write('TYPE=SPE' + '\n')
            else:
                f.write('TYPE=FID' + '\n')
            f.write('DATA' + '\n')
            if self.masterData.data.ndim is 1:
                for Line in self.masterData.data:
                    f.write(str(Line.real) + ' ' + str(Line.imag) + '\n')
            if self.masterData.data.ndim is 2:
                Points = self.masterData.data.shape
                for iii in range(0, Points[0]):
                    for jjj in range(0, Points[1]):
                        f.write(str(self.masterData.data[iii][jjj].real) + ' ' + str(self.masterData.data[iii][jjj].imag) + '\n')
            f.write('END')

    def saveASCIIFile(self):
        if self.masterData.data.ndim > 2:
            self.father.dispMsg('Saving to ASCII format only allowed for 1D and 2D data!')
            return
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum]  # Set name of file to be saved to workspace name to start
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.LastLocation + os.path.sep + WorkspaceName + '.txt', 'ASCII file (*.txt)')
        if type(name) is tuple:
            name = name[0]        
        if not name:
            return
        self.father.LastLocation = os.path.dirname(name)  # Save used path
        axis = np.array([self.masterData.xaxArray[-1]]).transpose()
        if self.masterData.data.ndim == 1:  # create nx1 matrix if it is a 1d data set
            data = np.array([self.masterData.data]).transpose()
        else:
            data = self.masterData.data.transpose()
        splitdata = np.zeros([data.shape[0], data.shape[1] * 2])
        for line in np.arange(data.shape[1]):
            splitdata[:, line * 2] = np.real(data[:, line])
            splitdata[:, line * 2 + 1] = np.imag(data[:, line])
        data = np.concatenate((axis, splitdata), axis=1)
        np.savetxt(name, data, delimiter='\t')

    def reloadLast(self):
        self.redoList = []
        path = self.masterData.filePath[1]
        if path.endswith('.zip'):
            import tempfile
            import shutil
            import zipfile
            temp_dir = tempfile.mkdtemp()
            zipfile.ZipFile(path).extractall(temp_dir)
            loadData = self.father.loading(self.masterData.filePath[0], temp_dir, True,realpath=path) 
        else:
            loadData = self.father.loading(self.masterData.filePath[0], self.masterData.filePath[1], True) 
        self.undoList.append(self.masterData.restoreData(loadData, None))
        self.current.upd()
        self.current.plotReset()
        self.current.showFid()
        self.updAllFrames()
        self.addMacro(['reload'])
        self.menuCheck()

    def monitorLoad(self, filePath):
        if not os.path.exists(filePath):
            self.stopMonitor()
            return
        self.redoList = []
        self.undoList = []
        loadData = self.father.loading(self.masterData.filePath[0], self.masterData.filePath[1], True)
        self.masterData.restoreData(loadData, None)
        for name in self.monitorMacros:
            self.runMacro(self.father.macros[name], display=False)
        self.current.upd() 
        #self.current.plotReset()  
        self.current.showFid()  
        self.updAllFrames()
        self.menuCheck()
        if filePath in self.monitor.files() or filePath in self.monitor.directories():
            return
        self.monitor.addPath(filePath)

    def startMonitor(self, macroNames):
        self.monitorMacros = macroNames
        self.monitor = QtCore.QFileSystemWatcher([self.masterData.filePath[1]], self)
        self.monitor.fileChanged.connect(self.monitorLoad)
        self.monitor.directoryChanged.connect(self.monitorLoad)

    def stopMonitor(self):
        self.monitorMacros = []
        if self.monitor is not None:
            self.monitor.removePath(self.masterData.filePath[1])
        self.monitor = None
        
    def real(self):
        self.redoList = []
        self.undoList.append(self.masterData.real())
        self.current.upd()
        self.current.showFid()
        self.addMacro(['real'])
        self.menuCheck()

    def imag(self):
        self.redoList = []
        self.undoList.append(self.masterData.imag())
        self.current.upd()
        self.current.showFid()
        self.addMacro(['imag'])
        self.menuCheck()

    def abs(self):
        self.redoList = []
        self.undoList.append(self.masterData.abs())
        self.current.upd()
        self.current.showFid()
        self.addMacro(['abs'])
        self.menuCheck()

    def fourier(self):
        self.redoList = []
        self.undoList.append(self.current.fourier())
        self.bottomframe.upd()
        self.menuCheck()

    def realFourier(self):
        self.redoList = []
        self.undoList.append(self.current.realFourier())
        self.bottomframe.upd()
        self.menuCheck()

    def fftshift(self):
        self.redoList = []
        self.undoList.append(self.current.fftshift())
        self.updAllFrames()
        self.menuCheck()

    def invFftshift(self):
        self.redoList = []
        self.undoList.append(self.current.fftshift(inv=True))
        self.updAllFrames()
        self.menuCheck()

    def diff(self):
        self.redoList = []
        self.undoList.append(self.current.diff())
        self.updAllFrames()
        self.menuCheck()

    def cumsum(self):
        self.redoList = []
        self.undoList.append(self.current.cumsum())
        self.updAllFrames()
        self.menuCheck()

    def hilbert(self):
        self.redoList = []
        self.undoList.append(self.current.hilbert())
        self.menuCheck()

    def states(self):
        self.redoList = []
        self.undoList.append(self.current.states())
        self.updAllFrames()
        self.menuCheck()

    def statesTPPI(self):
        self.redoList = []
        self.undoList.append(self.current.statesTPPI())
        self.updAllFrames()
        self.menuCheck()

    def echoAntiEcho(self):
        self.redoList = []
        self.undoList.append(self.current.echoAntiEcho())
        self.updAllFrames()
        self.menuCheck()

    def setFreq(self, freq, sw):
        self.redoList = []
        self.undoList.append(self.current.setFreq(freq, sw))
        self.menuCheck()

    def flipLR(self):
        self.redoList = []
        self.undoList.append(self.current.flipLR())
        self.menuCheck()

    def BrukerDigital(self):
        FilePath = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.LastLocation)
        if type(FilePath) is tuple:
            FilePath = FilePath[0]
        self.father.LastLocation = os.path.dirname(FilePath)  # Save used path
        if FilePath is '':
            return
        Dir = os.path.dirname(FilePath)
        if not os.path.exists(Dir + os.path.sep + 'acqus'):
            self.father.dispMsg("acqus file does not exist")
            return
        with open(Dir + os.path.sep + 'acqus', 'r') as f:
            data = f.read().split('\n')
        FilterCorrection = -1.0
        for s in range(0, len(data)):
            if data[s].startswith('##$GRPDLY='):
                FilterCorrection = float(data[s][10:])
            if data[s].startswith('##$DECIM='):
                DECIM = int(float(data[s][9:]))
            if data[s].startswith('##$DSPFVS='):
                DSPFVS = int(float(data[s][10:]))
        if DSPFVS == 10 or DSPFVS == 11 or DSPFVS == 12:  # get from table
            CorrectionList = [{'2': 44.7500, '3': 33.5000, '4': 66.6250, '6': 59.0833, '8': 68.5625, '12': 60.3750,
                               '16': 69.5313, '24': 61.0208, '32': 70.0156, '48': 61.3438, '64': 70.2578, '96': 61.5052,
                               '128': 70.3789, '192': 61.5859, '256': 70.4395, '384': 61.6263, '512': 70.4697, '768': 61.6465,
                               '1024': 70.4849, '1536': 61.6566, '2048': 70.4924},
                              {'2': 46.0000, '3': 36.5000, '4': 48.0000, '6': 50.1667, '8': 53.2500, '12': 69.5000,
                               '16': 72.2500, '24': 70.1667, '32': 72.7500, '48': 70.5000, '64': 73.0000, '96': 70.6667,
                               '128': 72.5000, '192': 71.3333, '256': 72.2500, '384': 71.6667, '512': 72.1250, '768': 71.8333,
                               '1024': 72.0625, '1536': 71.9167, '2048': 72.0313},
                              {'2': 46.311, '3': 36.530, '4': 47.870, '6': 50.229, '8': 53.289, '12': 69.551, '16': 71.600,
                               '24': 70.184, '32': 72.138, '48': 70.528, '64': 72.348, '96': 70.700, '128': 72.524}]
            # Take correction from database. Based on matNMR routine (Jacco van Beek), which is itself based
            # on a text by W. M. Westler and F. Abildgaard.
            FilterCorrection = CorrectionList[10 - DSPFVS][str(DECIM)]
        if FilterCorrection == -1.0:
            self.father.dispMsg('DSPFVS value not recognized (Bruker hardware version not known)')
            return
        if FilterCorrection != -1.0:  # If changed
            self.redoList = []
            self.undoList.append(self.current.applyPhase(0, FilterCorrection * 2 * np.pi))
            self.menuCheck()

    def createIntegralsWindow(self):
        self.father.createFitWindow(fit.IntegralsWindow(self.father, self.father.mainWindow))

    def createRelaxWindow(self):
        self.father.createFitWindow(fit.RelaxWindow(self.father, self.father.mainWindow))

    def createDiffusionWindow(self):
        self.father.createFitWindow(fit.DiffusionWindow(self.father, self.father.mainWindow))

    def createPeakDeconvWindow(self):
        self.father.createFitWindow(fit.PeakDeconvWindow(self.father, self.father.mainWindow))

    def createTensorDeconvWindow(self):
        self.father.createFitWindow(fit.TensorDeconvWindow(self.father, self.father.mainWindow))

    def createHerzfeldBergerWindow(self):
        self.father.createFitWindow(fit.HerzfeldBergerWindow(self.father, self.father.mainWindow))

    def createQuad1MASDeconvWindow(self):
        self.father.createFitWindow(fit.Quad1MASDeconvWindow(self.father, self.father.mainWindow))

    def createQuad1DeconvWindow(self):
        self.father.createFitWindow(fit.Quad1DeconvWindow(self.father, self.father.mainWindow))

    def createQuad2StaticDeconvWindow(self):
        if self.current.freq == 0.0:
            self.father.dispMsg("Please set the spectrometer frequency first!")
            return
        self.father.createFitWindow(fit.Quad2DeconvWindow(self.father, self.father.mainWindow))

    def createQuad2MASDeconvWindow(self):
        if self.current.freq == 0.0:
            self.father.dispMsg("Please set the spectrometer frequency first!")
            return
        self.father.createFitWindow(fit.Quad2DeconvWindow(self.father, self.father.mainWindow, True))

    def createQuad2StaticCzjzekWindow(self):
        if self.current.freq == 0.0:
            self.father.dispMsg("Please set the spectrometer frequency first!")
            return
        self.father.createFitWindow(fit.Quad2CzjzekWindow(self.father, self.father.mainWindow))

    def createQuad2MASCzjzekWindow(self):
        if self.current.freq == 0.0:
            self.father.dispMsg("Please set the spectrometer frequency first!")
            return
        self.father.createFitWindow(fit.Quad2CzjzekWindow(self.father, self.father.mainWindow, True))

    def plot1D(self):
        tmpcurrent = sc.Current1D(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()

    def plotScatter(self):
        tmpcurrent = sc.CurrentScatter(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()

    def plotStack(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentStacked(self, self.fig, self.canvas, self.masterData, self.current)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            self.father.dispMsg("Data does not have enough dimensions")

    def plotArray(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentArrayed(self, self.fig, self.canvas, self.masterData, self.current)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            self.father.dispMsg("Data does not have enough dimensions")

    def plotContour(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentContour(self, self.fig, self.canvas, self.masterData, self.current)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            self.father.dispMsg("Data does not have enough dimensions")

    def plotSkewed(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentSkewed(self, self.fig, self.canvas, self.masterData, self.current)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            self.father.dispMsg("Data does not have enough dimensions")

    def plotMulti(self):
        tmpcurrent = sc.CurrentMulti(self, self.fig, self.canvas, self.masterData, self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()

    def updAllFrames(self):
        self.sideframe.upd()
        self.bottomframe.upd()
        self.textframe.upd()

    def undo(self, *args):
        undoFunc = None
        while undoFunc is None and self.undoList:
            undoFunc = self.undoList.pop()
        if undoFunc is None:
            self.father.dispMsg("no undo information")
            return
        self.redoList.append(undoFunc(self.masterData))
        self.masterData.removeFromHistory(2)
        self.current.upd()
        self.current.plotReset()
        self.current.showFid()
        self.updAllFrames()
        if self.currentMacro is not None:
            self.redoMacro.append(self.father.macros[self.currentMacro].pop())
        self.menuCheck()

    def redo(self, *args):
        if self.redoList:
            self.undoList.append(self.redoList.pop()(self.masterData))
            self.current.upd()
            self.current.plotReset()
            self.current.showFid()
            self.updAllFrames()
            if self.currentMacro is not None:
                self.father.macroAdd(self.currentMacro, self.redoMacro.pop())
            self.menuCheck()
        else:
            self.father.dispMsg("no redo information")
            
    def clearUndo(self):
        self.undoList = []
        self.redoList = []
        self.menuCheck()
        

########################################################################################


class SideFrame(QtWidgets.QScrollArea):

    def __init__(self, parent):
        QtWidgets.QScrollArea.__init__(self, parent)
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

    def frameEnable(self):
        self.setEnabled(True)

    def frameDisable(self):
        self.setEnabled(False)

    def upd(self):
        current = self.father.current
        self.shape = current.data.data.shape
        self.length = len(self.shape)
        for i in reversed(range(self.frame1.count())):
            item = self.frame1.itemAt(i).widget()
            self.frame1.removeWidget(item)
            item.deleteLater()
        for i in reversed(range(self.frame2.count())):
            item = self.frame2.itemAt(i).widget()
            self.frame2.removeWidget(item)
            item.deleteLater()
        offset = 0
        self.plotIs2D = isinstance(current, (sc.CurrentStacked, sc.CurrentArrayed, sc.CurrentContour, sc.CurrentSkewed))
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
                self.buttons1.append(QtWidgets.QRadioButton(''))
                self.buttons1Group.addButton(self.buttons1[num], num)
                self.frame1.addWidget(self.buttons1[num], num * 2 + 1, 0)
                if self.plotIs2D:
                    self.buttons2.append(QtWidgets.QRadioButton(''))
                    self.buttons2Group.addButton(self.buttons2[num], num)
                    self.frame1.addWidget(self.buttons2[num], num * 2 + 1, 1)
                self.frame1.addWidget(wc.QLabel("D" + str(num + 1), self), num * 2, 1 + offset)
                self.entries.append(wc.SliceSpinBox(self, 0, self.shape[num] - 1))
                self.frame1.addWidget(self.entries[num], num * 2 + 1, 1 + offset)
                if not self.plotIs2D:
                    if num < current.axes:
                        self.entries[num].setValue(current.locList[num])
                    elif num == current.axes:
                        self.entries[num].setValue(0)
                    else:
                        self.entries[num].setValue(current.locList[num - 1])
                else:
                    if (num < current.axes) and (num < current.axes2):
                        self.entries[num].setValue(current.locList[num])
                    elif (num == current.axes) or (num == current.axes2):
                        self.entries[num].setValue(0)
                    elif (num > current.axes) and (num > current.axes2):
                        self.entries[num].setValue(current.locList[num - 2])
                    else:
                        self.entries[num].setValue(current.locList[num - 1])
                self.entries[num].valueChanged.connect(lambda event=None, num=num: self.getSlice(event, num))
            if isinstance(current, (sc.CurrentStacked, sc.CurrentArrayed, sc.CurrentSkewed)):
                if current.stackBegin is not None:
                    from2D = current.stackBegin
                else:
                    from2D = 0
                if current.stackEnd is not None:
                    to2D = current.stackEnd
                else:
                    to2D = self.shape[current.axes2]
                if current.stackStep is not None:
                    step2D = current.stackStep
                else:
                    step2D = 1
                self.frame2.addWidget(wc.QLabel("From", self), 1, 0)
                self.fromSpin = wc.SliceSpinBox(self, 0, to2D - 1)
                self.frame2.addWidget(self.fromSpin, 2, 0)
                self.fromSpin.setValue(from2D)
                self.fromSpin.valueChanged.connect(self.setToFrom)
                self.frame2.addWidget(wc.QLabel("To", self), 3, 0)
                self.toSpin = wc.SliceSpinBox(self, from2D + 1, self.shape[current.axes2])
                self.frame2.addWidget(self.toSpin, 4, 0)
                self.toSpin.setValue(to2D)
                self.toSpin.valueChanged.connect(self.setToFrom)
                self.frame2.addWidget(wc.QLabel("Step", self), 5, 0)
                self.stepSpin = wc.SliceSpinBox(self, 1, self.shape[current.axes2])
                self.frame2.addWidget(self.stepSpin, 6, 0)
                self.stepSpin.setValue(step2D)
                self.stepSpin.valueChanged.connect(self.setToFrom)
                if isinstance(current, (sc.CurrentStacked, sc.CurrentArrayed)):
                    self.frame2.addWidget(wc.QLabel("Spacing", self), 7, 0)
                    self.spacingEntry = QtWidgets.QLineEdit(self)
                    self.spacingEntry.setText('%#.3g' % current.spacing)
                    self.spacingEntry.returnPressed.connect(self.setSpacing)
                    self.frame2.addWidget(self.spacingEntry, 8, 0)
                elif isinstance(current, (sc.CurrentSkewed)):
                    self.frame2.addWidget(wc.QLabel("Skew", self), 7, 0)
                    self.skewEntry = QtWidgets.QLineEdit(self)
                    self.skewEntry.setText('%.2f' % current.skewed)
                    self.skewEntry.returnPressed.connect(self.setSkew)
                    self.frame2.addWidget(self.skewEntry, 8, 0)
                    self.frame2.addWidget(wc.QLabel("Elevation", self), 9, 0)
                    self.elevEntry = QtWidgets.QLineEdit(self)
                    self.elevEntry.setText('%.1f' % current.elevation)
                    self.elevEntry.returnPressed.connect(self.setSkew)
                    self.frame2.addWidget(self.elevEntry, 10, 0)
            if isinstance(current, (sc.CurrentContour)):
                self.frame2.addWidget(wc.QLabel("Number of contours:", self), 1, 0)
                self.numLEntry = QtWidgets.QLineEdit(self)
                self.numLEntry.setText(str(current.numLevels))
                self.numLEntry.returnPressed.connect(self.setContour)
                self.frame2.addWidget(self.numLEntry, 2, 0)
                self.contourTypeEntry = QtWidgets.QComboBox()
                self.contourTypeEntry.addItems(['Linear','Multiplier'])
                self.contourTypeEntry.currentIndexChanged.connect(self.setContour)
                self.frame2.addWidget(wc.QLabel("Contour scale:", self), 3, 0)
                self.frame2.addWidget(self.contourTypeEntry, 4, 0)
                self.multiValueLabel = wc.QLabel("Multiplier value:", self)
                self.frame2.addWidget(self.multiValueLabel, 5, 0)
                self.multiValueLabel.hide()
                self.multiValue = QtWidgets.QLineEdit(self)
                self.multiValue.setText(str(2.0))
                self.multiValue.returnPressed.connect(self.setContour)
                self.frame2.addWidget(self.multiValue, 6, 0)
                self.multiValue.hide()
                self.frame2.addWidget(wc.QLabel("Highest contour [%]:", self), 7, 0)
                self.maxLEntry = QtWidgets.QLineEdit(self)
                self.maxLEntry.setText(str(current.maxLevels * 100.0))
                self.maxLEntry.returnPressed.connect(self.setContour)
                self.frame2.addWidget(self.maxLEntry, 8, 0)
                self.frame2.addWidget(wc.QLabel("Lowest contour [%]:", self), 9, 0)
                self.minLEntry = QtWidgets.QLineEdit(self)
                self.minLEntry.setText(str(current.minLevels * 100.0))
                self.minLEntry.returnPressed.connect(self.setContour)
                self.frame2.addWidget(self.minLEntry, 10, 0)
            self.buttons1Group.button(current.axes).toggle()
            if self.plotIs2D:
                self.buttons2Group.button(current.axes2).toggle()
        if isinstance(current, (sc.CurrentMulti)):
            self.extraEntries = []
            self.extraButtons1 = []
            self.extraButtons1Group = []
            self.nameLabels = []
            iter1 = 0
            for i in range(len(current.extraData)):
                frameWidget = QtWidgets.QWidget(self)
                frame = QtWidgets.QGridLayout(frameWidget)
                self.frame2.addWidget(frameWidget, iter1, 0)
                frameWidget.setLayout(frame)
                name = current.extraName[i]
                if len(name) > 20:
                    name = name[:20]
                self.nameLabels.append(wc.QLabel(name, self))
                frame.addWidget(self.nameLabels[i], 0, 0, 1, 2)
                self.nameLabels[i].setStyleSheet("QLabel { color: rgb" + str(current.getExtraColor(i)) + ";}")
                colorbutton = QtWidgets.QPushButton("Color", self)
                colorbutton.clicked.connect(lambda arg, num=i: self.setExtraColor(num))
                frame.addWidget(colorbutton, 1, 0)
                button = QtWidgets.QPushButton("x", self)
                button.clicked.connect(lambda arg, num=i: self.delMultiSpec(num))
                frame.addWidget(button, 1, 1)
                self.OOM = self.father.current.getOOM()  # Order of Magnitude
                frame.addWidget(wc.QLabel("Scale", self), 2, 0)
                frame.addWidget(wc.QLabel("Offset (X1e" + str(self.OOM) + ")", self), 2, 1)
                scaleEntry = QtWidgets.QDoubleSpinBox()
                scaleEntry.setMaximum(1e3)
                scaleEntry.setMinimum(-1e3)
                scaleEntry.setSingleStep(0.1)
                scaleEntry.setValue(self.father.current.extraScale[i])
                scaleEntry.valueChanged.connect(lambda arg, num=i: self.setScale(arg, num))
                frame.addWidget(scaleEntry, 3, 0)
                offsetEntry = QtWidgets.QDoubleSpinBox()
                offsetEntry.setMaximum(1e3)
                offsetEntry.setMinimum(-1e3)
                offsetEntry.setSingleStep(0.1)
                offsetEntry.setValue(self.father.current.extraOffset[i]/(10**self.OOM))
                offsetEntry.valueChanged.connect(lambda arg, num=i: self.setOffset(arg, num))
                frame.addWidget(offsetEntry, 3, 1)
                entries = []
                self.extraEntries.append(entries)
                buttons1 = []
                self.extraButtons1.append(buttons1)
                self.extraButtons1Group.append(QtWidgets.QButtonGroup(self))
                self.extraButtons1Group[i].buttonClicked.connect(lambda: self.setExtraAxes(True))
                if current.extraData[i].data.ndim > 1:
                    for num in range(current.extraData[i].data.ndim):
                        buttons1.append(QtWidgets.QRadioButton(''))
                        self.extraButtons1Group[i].addButton(buttons1[num], num)
                        frame.addWidget(buttons1[num], num * 2 + 5, 0)
                        frame.addWidget(wc.QLabel("D" + str(num + 1), self), num * 2 + 4, 1)
                        entries.append(wc.SliceSpinBox(self, 0, current.extraData[i].data.shape[num] - 1))
                        frame.addWidget(entries[num], num * 2 + 5, 1)
                        if num < current.extraAxes[i]:
                            entries[num].setValue(current.extraLoc[i][num])
                        elif num == current.extraAxes[i]:
                            entries[num].setValue(0)
                        else:
                            entries[num].setValue(current.extraLoc[i][num - 1])
                        entries[num].valueChanged.connect(lambda event=None, num=num, i=i: self.getExtraSlice(event, num, i))
                    self.extraButtons1Group[i].button(current.extraAxes[i]).toggle()
                iter1 += 1
            addButton = QtWidgets.QPushButton("Add spectrum", self)
            addButton.clicked.connect(self.addMultiSpec)
            self.frame2.addWidget(addButton, iter1, 0, 1, 2)
        QtCore.QTimer.singleShot(100, self.resizeAll)

    def resizeAll(self):
        self.setMinimumWidth(self.grid.sizeHint().width() + self.verticalScrollBar().sizeHint().width())

    def setToFrom(self, *args):
        current = self.father.current
        if not isinstance(current, (sc.CurrentStacked, sc.CurrentArrayed, sc.CurrentSkewed)):
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
        var = float(safeEval(self.spacingEntry.text()))
        self.spacingEntry.setText('%#.3g' % var)
        self.father.current.setSpacing(var)

    def setSkew(self, *args):
        var = float(safeEval(self.skewEntry.text()))
        self.skewEntry.setText('%.2f' % var)
        var2 = float(safeEval(self.elevEntry.text()))
        self.elevEntry.setText('%.1f' % var2)
        self.father.current.setSkewed(var, var2)

    def setContour(self, *args):
        var1 = int(round(safeEval(self.numLEntry.text())))
        self.numLEntry.setText(str(var1))
        var2 = abs(float(safeEval(self.maxLEntry.text())))
        var3 = abs(float(safeEval(self.minLEntry.text())))
        if var3 > var2: #if wrong order, interchange
            var2, var3 = (var3 , var2)
        self.maxLEntry.setText('%.1f' % var2)
        self.minLEntry.setText('%.1f' % var3)
        var4 = self.contourTypeEntry.currentIndex()
        if var4 == 0:
            self.multiValue.hide()
            self.multiValueLabel.hide()
        else:
            self.multiValue.show()
            self.multiValueLabel.show()
        var5 = abs(float(safeEval(self.multiValue.text())))
        self.father.current.setLevels(var1, var2 / 100.0, var3 / 100.0, var4, var5)

    def setAxes(self, first=True):
        axes = self.buttons1Group.checkedId()
        if self.plotIs2D:
            axes2 = self.buttons2Group.checkedId()
            if axes == axes2:
                if first:
                    axes2 = self.father.current.axes
                else:
                    axes = self.father.current.axes2
            self.buttons2Group.button(axes2).toggle()
        self.getSlice(None, axes, True)
        self.upd()

    def getSlice(self, event, entryNum, button=False):
        if button:
            dimNum = entryNum
        elif not self.plotIs2D:
            if entryNum == self.father.current.axes:
                if entryNum == self.length - 1:
                    dimNum = self.length - 2
                else:
                    dimNum = self.length - 1
            else:
                dimNum = self.father.current.axes
        else:
            dimNum = self.father.current.axes
        locList = []
        for num in range(self.length):
            appendLoc = True
            if self.plotIs2D and (num == self.buttons2Group.checkedId()):
                appendLoc = False
            inp = self.entries[num].value()
            if num == dimNum:
                pass
            else:
                if appendLoc:
                    locList.append(inp)
        self.buttons1Group.button(dimNum).toggle()
        if self.plotIs2D:
            self.father.current.setBlock(dimNum, self.buttons2Group.checkedId(), locList)
        else:
            self.father.current.setSlice(dimNum, locList)
        self.father.bottomframe.upd()
        # self.upd()

    def setExtraAxes(self, first=True):
        for i in range(len(self.extraButtons1Group)):
            axes = self.extraButtons1Group[i].checkedId()
            self.getExtraSlice(None, axes, i, True)
        self.father.current.showFid()

    def getExtraSlice(self, event, entryNum, entryi, button=False):
        length = self.father.current.extraData[entryi].data.ndim
        if button:
            dimNum = entryNum
        else:
            if entryNum == self.father.current.extraAxes[entryi]:
                if entryNum == length - 1:
                    dimNum = length - 2
                else:
                    dimNum = length - 1
            else:
                dimNum = self.father.current.extraAxes[entryi]
        locList = []
        for num in range(length):
            inp = self.extraEntries[entryi][num].value()
            if num == dimNum:
                pass
            else:
                locList.append(inp)
        self.extraButtons1Group[entryi].button(dimNum).toggle()
        self.father.current.setExtraSlice(entryi, dimNum, locList)
        if not button:
            self.father.current.showFid()
        # self.upd()

    def setScale(self, scale, num):
        self.father.current.setExtraScale(num, scale)

    def setOffset(self, offset, num):
        self.father.current.setExtraOffset(num, offset*10**self.OOM)

    def checkChanged(self):
        for i in range(len(self.father.current.extraData)):
            extraData = self.father.current.extraData[i]
            if extraData.data.ndim > 1:
                for j in range(len(self.extraEntries[i])):
                    self.extraEntries[i][j].setMaximum(extraData.data.shape[j] - 1)
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
        QtWidgets.QWidget.__init__(self, parent)
        self.father = parent
        grid = QtWidgets.QGridLayout(self)
        self.setLayout(grid)
        fourierButton = QtWidgets.QPushButton("Fourier", parent=self)
        fourierButton.clicked.connect(self.father.fourier)
        grid.addWidget(fourierButton, 0, 0, 2, 1)
        self.specGroup = QtWidgets.QButtonGroup(self)
        self.specGroup.buttonClicked.connect(self.changeSpec)
        timeButton = QtWidgets.QRadioButton('Time', parent=self)
        self.specGroup.addButton(timeButton, 0)
        grid.addWidget(timeButton, 0, 1)
        freqButton = QtWidgets.QRadioButton('Frequency', parent=self)
        self.specGroup.addButton(freqButton, 1)
        grid.addWidget(freqButton, 1, 1)
        self.wholeEcho = QtWidgets.QCheckBox("Whole echo", parent=self)
        self.wholeEcho.clicked.connect(self.setWholeEcho)
        grid.addWidget(self.wholeEcho, 0, 2, 2, 1)
        grid.addWidget(wc.QLabel("Freq [MHz]:", self), 0, 3)
        self.freqEntry = QtWidgets.QLineEdit(parent=self)
        self.freqEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.freqEntry.returnPressed.connect(self.changeFreq)
        grid.addWidget(self.freqEntry, 1, 3)
        grid.addWidget(wc.QLabel("Sweepwidth [kHz]:", self), 0, 4)
        self.swEntry = QtWidgets.QLineEdit(parent=self)
        self.swEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.swEntry.returnPressed.connect(self.changeFreq)
        grid.addWidget(self.swEntry, 1, 4)
        grid.addWidget(wc.QLabel("Plot:", self), 0, 5)
        self.plotDrop = QtWidgets.QComboBox(parent=self)
        self.plotDrop.addItems(["Real", "Imag", "Both", "Abs"])
        self.plotDrop.activated.connect(self.changePlot)
        grid.addWidget(self.plotDrop, 1, 5)
        grid.addWidget(wc.QLabel("Axis:", self), 0, 6)
        self.axisDropTime = QtWidgets.QComboBox(parent=self)
        self.axisDropTime.addItems(["s", "ms", u"\u03bcs"])
        self.axisDropTime.activated.connect(self.changeAxis)
        grid.addWidget(self.axisDropTime, 1, 6)
        self.axisDropFreq = QtWidgets.QComboBox(parent=self)
        self.axisDropFreq.addItems(["Hz", "kHz", "MHz", "ppm"])
        self.axisDropFreq.activated.connect(self.changeAxis)
        grid.addWidget(self.axisDropFreq, 1, 6)
        self.ax2Label = wc.QLabel("Axis2:", self)
        grid.addWidget(self.ax2Label, 0, 7)
        self.axisDropTime2 = QtWidgets.QComboBox(parent=self)
        self.axisDropTime2.addItems(["s", "ms", u"\u03bcs"])
        self.axisDropTime2.activated.connect(self.changeAxis2)
        grid.addWidget(self.axisDropTime2, 1, 7)
        self.axisDropFreq2 = QtWidgets.QComboBox(parent=self)
        self.axisDropFreq2.addItems(["Hz", "kHz", "MHz", "ppm"])
        self.axisDropFreq2.activated.connect(self.changeAxis2)
        grid.addWidget(self.axisDropFreq2, 1, 7)
        self.proj1Label = wc.QLabel("Proj top:", self)
        grid.addWidget(self.proj1Label, 0, 8)
        self.projDrop1 = QtWidgets.QComboBox(parent=self)
        self.projDrop1.addItems(["sum", "max", "min"])
        self.projDrop1.activated.connect(lambda val, self=self: self.changeProj(val, 1))
        grid.addWidget(self.projDrop1, 1, 8)
        self.proj2Label = wc.QLabel("Proj right:", self)
        grid.addWidget(self.proj2Label, 0, 9)
        self.projDrop2 = QtWidgets.QComboBox(parent=self)
        self.projDrop2.addItems(["sum", "max", "min"])
        self.projDrop2.activated.connect(lambda val, self=self: self.changeProj(val, 2))
        grid.addWidget(self.projDrop2, 1, 9)
        grid.setColumnStretch(10, 1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        self.grid = grid
        self.upd()

    def kill(self):
        for i in reversed(range(self.grid.count())):
            self.grid.itemAt(i).widget().deleteLater()
        self.grid.deleteLater()

    def frameEnable(self):
        self.setEnabled(True)

    def frameDisable(self):
        self.setEnabled(False)

    def upd(self):
        self.freqEntry.setText('%.6f' % (self.father.current.freq / 1000000.0))
        self.swEntry.setText('%.6f' % (self.father.current.sw / 1000.0))
        self.axisDropTime2.hide()
        self.axisDropFreq2.hide()
        self.proj1Label.hide()
        self.proj2Label.hide()
        self.projDrop1.hide()
        self.projDrop2.hide()
        self.axisDropFreq.model().item(3).setEnabled(True)
        if self.father.current.spec == 0:
            self.specGroup.button(0).toggle()
            self.axisDropFreq.hide()
            self.axisDropTime.show()
            self.ax2Label.hide()
            self.axisDropTime.setCurrentIndex(self.father.current.axType)
        elif self.father.current.spec == 1:
            self.specGroup.button(1).toggle()
            self.axisDropTime.hide()
            self.axisDropFreq.show()
            if self.father.current.freq == 0.0:
                self.axisDropFreq.model().item(3).setEnabled(False)
            self.ax2Label.hide()
            if self.father.current.ppm:
                self.axisDropFreq.setCurrentIndex(3)
            else:
                self.axisDropFreq.setCurrentIndex(self.father.current.axType)
        if isinstance(self.father.current, sc.CurrentContour):
            self.ax2Label.show()
            self.proj1Label.show()
            self.proj2Label.show()
            self.projDrop1.show()
            self.projDrop2.show()
            self.axisDropFreq2.model().item(3).setEnabled(True)
            if self.father.current.spec2 == 0:
                self.axisDropTime2.show()
                self.axisDropTime2.setCurrentIndex(self.father.current.axType2)
            elif self.father.current.spec2 == 1:
                self.axisDropFreq2.show()
                if self.father.current.freq2 == 0.0:
                    self.axisDropFreq2.model().item(3).setEnabled(False)
                self.axisDropFreq2.setCurrentIndex(self.father.current.axType2)
        if self.father.current.wholeEcho:
            self.wholeEcho.setCheckState(QtCore.Qt.Checked)
        else:
            self.wholeEcho.setCheckState(QtCore.Qt.Unchecked)

    def setWholeEcho(self, inp):
        self.father.undoList.append(self.father.current.setWholeEcho(inp))
        self.father.menuCheck()

    def changeSpec(self):
        self.father.redoList = []
        self.father.undoList.append(self.father.current.changeSpec(self.specGroup.checkedId()))
        self.upd()
        self.father.menuCheck()

    def changeFreq(self):
        freq = safeEval(self.freqEntry.text()) * 1e6
        sw = safeEval(self.swEntry.text()) * 1e3
        if freq != 0 and sw != 0:
            self.father.setFreq(freq, sw)
        self.upd()

    def changePlot(self, pType):
        self.father.current.plotType = pType
        self.father.current.showFid()

    def changeAxis(self, pType):
        self.father.current.setAxType(pType)
        self.father.current.showFid()

    def changeAxis2(self, pType):
        self.father.current.setAxType2(pType)
        self.father.current.showFid()

    def changeProj(self, pType, direc):
        self.father.current.setProjType(pType, direc)
        self.father.current.showProj()

##################################################################


class TextFrame(QtWidgets.QScrollArea):

    def __init__(self, parent):
        QtWidgets.QScrollArea.__init__(self, parent)
        self.father = parent
        self.oldx = 0.0
        self.oldy = 0.0
        self.oldamp = 0.0
        widthScale = 0.6
        content = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout(content)
        getButton = QtWidgets.QPushButton("&Get Position")
        getButton.clicked.connect(self.getPosition)
        grid.addWidget(getButton, 0, 1)
        grid.addWidget(wc.QLabel("x-Position:"), 0, 2)
        self.xpos = QtWidgets.QLineEdit()
        self.xpos.setAlignment(QtCore.Qt.AlignHCenter)
        self.xpos.setFixedWidth(self.xpos.sizeHint().width() * widthScale)
        self.xpos.setText("0")
        grid.addWidget(self.xpos, 0, 3)
        self.yposlabel = wc.QLabel("y-Position:")
        grid.addWidget(self.yposlabel, 0, 4)
        self.ypos = QtWidgets.QLineEdit()
        self.ypos.setAlignment(QtCore.Qt.AlignHCenter)
        self.ypos.setFixedWidth(self.ypos.sizeHint().width() * widthScale)
        self.ypos.setText("0")
        grid.addWidget(self.ypos, 0, 5)
        grid.addWidget(wc.QLabel("x-Value:"), 0, 6)
        self.xpoint = QtWidgets.QLineEdit()
        self.xpoint.setAlignment(QtCore.Qt.AlignHCenter)
        self.xpoint.setFixedWidth(self.xpoint.sizeHint().width() * widthScale)
        self.xpoint.setText("0.0")
        grid.addWidget(self.xpoint, 0, 7)
        self.ylabel = wc.QLabel("y-Value:")
        grid.addWidget(self.ylabel, 0, 8)
        self.ypoint = QtWidgets.QLineEdit()
        self.ypoint.setAlignment(QtCore.Qt.AlignHCenter)
        self.ypoint.setFixedWidth(self.ypoint.sizeHint().width() * widthScale)
        self.ypoint.setText("0.0")
        grid.addWidget(self.ypoint, 0, 9)
        grid.addWidget(wc.QLabel("Amp:"), 0, 10)
        self.amppoint = QtWidgets.QLineEdit()
        self.amppoint.setAlignment(QtCore.Qt.AlignHCenter)
        self.amppoint.setFixedWidth(self.amppoint.sizeHint().width() * widthScale)
        self.amppoint.setText("0.0")
        grid.addWidget(self.amppoint, 0, 11)
        grid.addWidget(wc.QLabel(u"\u0394x:"), 0, 12)
        self.deltaxpoint = QtWidgets.QLineEdit()
        self.deltaxpoint.setAlignment(QtCore.Qt.AlignHCenter)
        self.deltaxpoint.setFixedWidth(self.deltaxpoint.sizeHint().width() * widthScale)
        self.deltaxpoint.setText("0.0")
        grid.addWidget(self.deltaxpoint, 0, 13)
        self.deltaylabel = wc.QLabel(u"\u0394y:")
        grid.addWidget(self.deltaylabel, 0, 14)
        self.deltaypoint = QtWidgets.QLineEdit()
        self.deltaypoint.setAlignment(QtCore.Qt.AlignHCenter)
        self.deltaypoint.setFixedWidth(self.deltaypoint.sizeHint().width() * widthScale)
        self.deltaypoint.setText("0.0")
        grid.addWidget(self.deltaypoint, 0, 15)
        grid.addWidget(wc.QLabel(u"\u0394amp:"), 0, 16)
        self.deltaamppoint = QtWidgets.QLineEdit()
        self.deltaamppoint.setAlignment(QtCore.Qt.AlignHCenter)
        self.deltaamppoint.setFixedWidth(self.deltaamppoint.sizeHint().width() * widthScale)
        self.deltaamppoint.setText("0.0")
        grid.addWidget(self.deltaamppoint, 0, 17)
        grid.setColumnStretch(20, 1)
        self.grid = grid
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setWidget(content)
        self.setMaximumHeight(self.grid.sizeHint().height() + self.horizontalScrollBar().sizeHint().height())
        self.upd()

    def upd(self):
        if isinstance(self.father.current, sc.CurrentContour):
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

    def frameEnable(self):
        for child in self.children():
            child.setEnabled(True)

    def frameDisable(self):
        for child in self.children():
            child.setEnabled(False)

    def setLabels(self, position):
        if len(position) > 3:
            self.ypos.setText(str(position[3]))
            self.deltaypoint.setText('%#.3g' % np.abs(self.oldy - position[4]))
            self.ypoint.setText('%#.3g' % position[4])
            self.oldy = position[3]
        self.deltaxpoint.setText('%#.3g' % np.abs(self.oldx - position[1]))
        self.deltaamppoint.setText('%#.3g' % np.abs(self.oldamp - position[2]))
        self.xpos.setText(str(position[0]))
        self.xpoint.setText('%#.3g' % position[1])
        self.amppoint.setText('%#.3g' % position[2])
        self.oldx = position[1]
        self.oldamp = position[2]

    def getPosition(self, *args):
        self.father.current.peakPickFunc = lambda pos, self=self: self.setLabels(pos)
        if isinstance(self.father.current, sc.CurrentContour):
            self.father.current.peakPick = 2
        else:
            self.father.current.peakPick = True

#################################################################################


class PhaseWindow(QtWidgets.QWidget):

    RESOLUTION = 10000.0
    P1LIMIT = 540.0
    PHASE0STEP = 1.0
    PHASE1STEP = 1.0

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent.father)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.zeroVal = 0.0
        self.firstVal = 0.0
        self.refVal = 0.0
        self.available = True
        self.setWindowTitle("Phasing")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Zero order phasing:"), 0, 0, 1, 3)
        autoZero = QtWidgets.QPushButton("Autophase 0th")
        autoZero.clicked.connect(lambda: self.autophase(0))
        grid.addWidget(autoZero, 1, 1)
        self.zeroEntry = QtWidgets.QLineEdit()
        self.zeroEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.zeroEntry.returnPressed.connect(self.inputZeroOrder)
        self.zeroEntry.setText("0.000")
        grid.addWidget(self.zeroEntry, 2, 1)
        leftZero = QtWidgets.QPushButton("<")
        leftZero.clicked.connect(lambda: self.stepPhase(-1, 0))
        leftZero.setAutoRepeat(True)
        grid.addWidget(leftZero, 2, 0)
        rightZero = QtWidgets.QPushButton(">")
        rightZero.clicked.connect(lambda: self.stepPhase(1, 0))
        rightZero.setAutoRepeat(True)
        grid.addWidget(rightZero, 2, 2)
        self.zeroScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.zeroScale.setRange(-self.RESOLUTION, self.RESOLUTION)
        self.zeroScale.valueChanged.connect(self.setZeroOrder)
        grid.addWidget(self.zeroScale, 3, 0, 1, 3)
        grid.addWidget(wc.QLabel("First order phasing:"), 4, 0, 1, 3)
        autoFirst = QtWidgets.QPushButton("Autophase 0th+1st")
        autoFirst.clicked.connect(lambda: self.autophase(1))
        grid.addWidget(autoFirst, 5, 1)
        self.firstEntry = QtWidgets.QLineEdit()
        self.firstEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.firstEntry.returnPressed.connect(self.inputFirstOrder)
        self.firstEntry.setText("0.000")
        grid.addWidget(self.firstEntry, 6, 1)
        leftFirst = QtWidgets.QPushButton("<")
        leftFirst.clicked.connect(lambda: self.stepPhase(0, -1))
        leftFirst.setAutoRepeat(True)
        grid.addWidget(leftFirst, 6, 0)
        rightFirst = QtWidgets.QPushButton(">")
        rightFirst.clicked.connect(lambda: self.stepPhase(0, 1))
        rightFirst.setAutoRepeat(True)
        grid.addWidget(rightFirst, 6, 2)
        self.firstScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.firstScale.setRange(-self.RESOLUTION, self.RESOLUTION)
        self.firstScale.valueChanged.connect(self.setFirstOrder)
        grid.addWidget(self.firstScale, 7, 0, 1, 3)
        if self.father.current.spec > 0:
            grid.addWidget(wc.QLabel("Reference:"), 8, 0, 1, 3)
            pickRef = QtWidgets.QPushButton("Pick reference")
            pickRef.clicked.connect(self.pickRef)
            grid.addWidget(pickRef, 9, 1)
            self.refEntry = QtWidgets.QLineEdit()
            self.refEntry.setAlignment(QtCore.Qt.AlignHCenter)
            self.refEntry.setText('%.3f' % self.refVal)
            self.refEntry.returnPressed.connect(self.inputRef)
            grid.addWidget(self.refEntry, 10, 1)

        self.singleSlice = QtWidgets.QCheckBox("Single slice")
        grid.addWidget(self.singleSlice, 11, 0, 1, 3)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 1, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def setZeroOrder(self, value, *args):
        if self.available:
            self.zeroVal = float(value) / self.RESOLUTION * 180
            self.zeroEntry.setText('%.3f' % self.zeroVal)
            self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0)

    def inputZeroOrder(self, *args):
        inp = safeEval(self.zeroEntry.text())
        self.zeroVal = np.mod(inp + 180, 360) - 180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.available = False
        self.zeroScale.setValue(round(self.zeroVal / 180.0 * self.RESOLUTION))
        self.available = True
        self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0)

    def setFirstOrder(self, value, *args):
        if self.available:
            value = float(value) / self.RESOLUTION * self.P1LIMIT
            newZero = (self.zeroVal - (value - self.firstVal) * self.refVal / self.father.current.sw)
            self.zeroVal = np.mod(newZero + 180, 360) - 180
            self.zeroEntry.setText('%.3f' % self.zeroVal)
            self.firstVal = value
            self.firstEntry.setText('%.3f' % self.firstVal)
            self.available = False
            self.zeroScale.setValue(round(self.zeroVal / 180.0 * self.RESOLUTION))
            self.available = True
            self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0)

    def inputFirstOrder(self, *args):
        value = float(safeEval(self.firstEntry.text()))
        newZero = (self.zeroVal - (value - self.firstVal) * self.refVal / self.father.current.sw)
        self.zeroVal = np.mod(newZero + 180, 360) - 180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.firstVal = value
        self.firstEntry.setText('%.3f' % self.firstVal)
        self.available = False
        self.zeroScale.setValue(round(self.zeroVal / 180.0 * self.RESOLUTION))
        self.firstScale.setValue(round(self.firstVal / self.P1LIMIT * self.RESOLUTION))
        self.available = True
        self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0)

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

    def stepPhase(self, phase0, phase1):
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            multiplier = 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            multiplier = 100
        else:
            multiplier = 1
        phase0 = multiplier * phase0
        phase1 = multiplier * phase1
        inp = safeEval(self.zeroEntry.text()) + phase0 * self.PHASE0STEP
        self.zeroVal = np.mod(inp + 180, 360) - 180
        value = safeEval(self.firstEntry.text()) + phase1 * self.PHASE1STEP
        newZero = (self.zeroVal - (value - self.firstVal) * self.refVal / self.father.current.sw)
        self.zeroVal = np.mod(newZero + 180, 360) - 180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.firstVal = value
        self.firstEntry.setText('%.3f' % self.firstVal)
        self.available = False
        self.zeroScale.setValue(round(self.zeroVal / 180.0 * self.RESOLUTION))
        self.firstScale.setValue(round(self.firstVal / self.P1LIMIT * self.RESOLUTION))
        self.available = True
        self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0)

    def inputRef(self, *args):
        self.refVal = safeEval(self.refEntry.text())
        self.refEntry.setText('%.3f' % self.refVal)

    def setRef(self, value, *args):
        self.refVal = float(value)
        self.refEntry.setText('%.3f' % self.refVal)

    def pickRef(self, *args):
        self.father.current.peakPickFunc = lambda pos, self=self: self.setRef(self.father.current.xax[pos[0]])
        self.father.current.peakPick = True

    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        self.inputZeroOrder()
        self.inputFirstOrder()
        self.father.redoList = []
        self.father.undoList.append(self.father.current.applyPhase(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, (self.singleSlice.isChecked() == 1)))
        self.father.menuEnable()
        self.deleteLater()

################################################################


class ApodWindow(QtWidgets.QWidget):

    RESOLUTION = 10000

    def __init__(self, parent):
        parent.menuDisable()
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.entries = []
        self.ticks = []
        self.maximum = 100.0 * self.father.current.sw / (self.father.current.data1D.shape[-1])
        self.lorstep = 1.0
        self.gaussstep = 1.0
        self.available = True
        self.setWindowTitle("Apodize")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        lorTick = QtWidgets.QCheckBox("Lorentzian:")
        lorTick.toggled.connect(lambda: self.checkEval(0))
        grid.addWidget(lorTick, 0, 0, 1, 3)
        self.ticks.append(lorTick)
        lorEntry = QtWidgets.QLineEdit()
        lorEntry.setEnabled(False)
        lorEntry.setAlignment(QtCore.Qt.AlignHCenter)
        lorEntry.setText("0.00")
        lorEntry.returnPressed.connect(self.apodPreview)
        grid.addWidget(lorEntry, 1, 1)
        self.entries.append(lorEntry)
        leftLor = QtWidgets.QPushButton("<")
        leftLor.clicked.connect(lambda: self.stepLB(-0.5 * self.father.current.sw / (self.father.current.data1D.shape[-1]), 0))
        leftLor.setAutoRepeat(True)
        grid.addWidget(leftLor, 1, 0)
        rightLor = QtWidgets.QPushButton(">")
        rightLor.clicked.connect(lambda: self.stepLB(0.5 * self.father.current.sw / (self.father.current.data1D.shape[-1]), 0))
        rightLor.setAutoRepeat(True)
        grid.addWidget(rightLor, 1, 2)
        self.lorScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.lorScale.setRange(0, self.RESOLUTION)
        self.lorScale.valueChanged.connect(self.setLor)
        grid.addWidget(self.lorScale, 2, 0, 1, 3)
        self.lorMax = 100.0 * self.father.current.sw / (self.father.current.data1D.shape[-1])

        gaussTick = QtWidgets.QCheckBox("Gaussian:")
        gaussTick.toggled.connect(lambda: self.checkEval(1))
        grid.addWidget(gaussTick, 3, 0, 1, 3)
        self.ticks.append(gaussTick)
        gaussEntry = QtWidgets.QLineEdit()
        gaussEntry.setEnabled(False)
        gaussEntry.setAlignment(QtCore.Qt.AlignHCenter)
        gaussEntry.setText("0.00")
        gaussEntry.returnPressed.connect(self.apodPreview)
        grid.addWidget(gaussEntry, 4, 1)
        self.entries.append(gaussEntry)
        leftGauss = QtWidgets.QPushButton("<")
        leftGauss.clicked.connect(lambda: self.stepLB(0, -0.5 * self.father.current.sw / (self.father.current.data1D.shape[-1])))
        leftGauss.setAutoRepeat(True)
        grid.addWidget(leftGauss, 4, 0)
        rightGauss = QtWidgets.QPushButton(">")
        rightGauss.clicked.connect(lambda: self.stepLB(0, 0.5 * self.father.current.sw / (self.father.current.data1D.shape[-1])))
        rightGauss.setAutoRepeat(True)
        grid.addWidget(rightGauss, 4, 2)
        self.gaussScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.gaussScale.setRange(0, self.RESOLUTION)
        self.gaussScale.valueChanged.connect(self.setGauss)
        grid.addWidget(self.gaussScale, 5, 0, 1, 3)
        self.gaussMax = 100.0 * self.father.current.sw / (self.father.current.data1D.shape[-1])

        cos2Tick = QtWidgets.QCheckBox("Cos^2:")
        cos2Tick.clicked.connect(lambda: self.checkEval(2))
        grid.addWidget(cos2Tick, 6, 0, 1, 3)
        self.ticks.append(cos2Tick)
        cos2Entry = QtWidgets.QLineEdit()
        cos2Entry.setEnabled(False)
        cos2Entry.setAlignment(QtCore.Qt.AlignHCenter)
        cos2Entry.setText("1.00")
        cos2Entry.returnPressed.connect(self.apodPreview)
        grid.addWidget(cos2Entry, 7, 1)
        self.entries.append(cos2Entry)

        hammingTick = QtWidgets.QCheckBox("Hamming:")
        hammingTick.clicked.connect(lambda: self.checkEval(3))
        grid.addWidget(hammingTick, 8, 0, 1, 3)
        self.ticks.append(hammingTick)
        hammingEntry = QtWidgets.QLineEdit()
        hammingEntry.setEnabled(False)
        hammingEntry.setAlignment(QtCore.Qt.AlignHCenter)
        hammingEntry.setText("1.00")
        hammingEntry.returnPressed.connect(self.apodPreview)
        grid.addWidget(hammingEntry, 9, 1)
        self.entries.append(hammingEntry)

        grid.addWidget(wc.QLabel("Shift:"), 10, 0, 1, 3)
        self.shiftEntry = QtWidgets.QLineEdit()
        self.shiftEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.shiftEntry.setText("0.00")
        self.shiftEntry.returnPressed.connect(self.apodPreview)
        grid.addWidget(self.shiftEntry, 11, 1)

        if self.father.current.data.data.ndim > 1:
            grid.addWidget(wc.QLabel("Shifting:"), 12, 0, 1, 3)
            self.shiftingEntry = QtWidgets.QLineEdit()
            self.shiftingEntry.setAlignment(QtCore.Qt.AlignHCenter)
            self.shiftingEntry.setText("0.00")
            self.shiftingEntry.returnPressed.connect(self.apodPreview)
            grid.addWidget(self.shiftingEntry, 13, 1)
            self.shiftingAxes = QtWidgets.QComboBox()
            self.shiftingValues = list(map(str, np.delete(range(1, self.father.current.data.data.ndim + 1), self.father.current.axes)))
            self.shiftingAxes.addItems(self.shiftingValues)
            self.shiftingAxes.currentIndexChanged.connect(self.apodPreview)
            grid.addWidget(self.shiftingAxes, 14, 1)

        self.singleSlice = QtWidgets.QCheckBox("Single slice")
        grid.addWidget(self.singleSlice, 15, 0, 1, 3)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 1, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def checkEval(self, num):
        if self.ticks[num].isChecked():
            self.entries[num].setEnabled(True)
        else:
            self.entries[num].setEnabled(False)
        if num == 0 or num == 1: #for lorentzian and gaussian
            if safeEval(self.entries[num].text()) != 0.0: #only update if value was not zero
               self.apodPreview()
        else:
            self.apodPreview()

    def setLor(self, value, *args):
        if self.available:
            self.entries[0].setText('%.4g' % (float(value) * self.maximum / self.RESOLUTION))
            if not self.ticks[0].isChecked():
                self.ticks[0].setChecked(1)
            self.apodPreview()

    def setGauss(self, value, *args):
        if self.available:
            self.entries[1].setText('%.4g' % (float(value) * self.maximum / self.RESOLUTION))
            if not self.ticks[1].isChecked():
                self.ticks[1].setChecked(1)
            self.apodPreview()

    def apodPreview(self, *args):
        self.available = False
        lor = None
        gauss = None
        cos2 = None
        hamming = None
        shifting = None
        shiftingAxes = 0
        if self.ticks[0].isChecked():
            lor = safeEval(self.entries[0].text())
            if lor is None:
                self.father.father.dispMsg('Lorentzian value is not valid!')
            self.entries[0].setText('%.4g' % lor)
            self.lorScale.setValue(round(lor * self.RESOLUTION / self.maximum))
        if self.ticks[1].isChecked():
            gauss = safeEval(self.entries[1].text())
            if gauss is None:
                self.father.father.dispMsg('Gaussian value is not valid!')
            self.entries[1].setText('%.4g' % gauss)
            self.gaussScale.setValue(round(gauss * self.RESOLUTION / self.maximum))
        if self.ticks[2].isChecked():
            cos2 = safeEval(self.entries[2].text())
            if cos2 is None:
                self.father.father.dispMsg('cos^2 value is not valid!')
            self.entries[2].setText('%.4g' % cos2)
        if self.ticks[3].isChecked():
            hamming = safeEval(self.entries[3].text())
            if hamming is None:
                self.father.father.dispMsg('Hamming value is not valid!')
            self.entries[3].setText('%.4g' % hamming)
        shift = safeEval(self.shiftEntry.text())
        if shift is None:
            self.father.father.dispMsg('Shift value is not valid!')
        self.shiftEntry.setText('%.4g' % shift)
        if self.father.current.data.data.ndim > 1:
            shifting = safeEval(self.shiftingEntry.text())
            if shifting is None:
                self.father.father.dispMsg('Shifting value is not valid!')
            self.shiftingEntry.setText('%.4g' % shifting)
            shiftingAxes = int(self.shiftingValues[self.shiftingAxes.currentIndex()]) - 1
        else:
            shiftingAxes = None
        self.available = True
        self.father.current.apodPreview(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes)

    def stepLB(self, lorincr, gaussincr):
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            multiplier = 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            multiplier = 100
        else:
            multiplier = 1
        self.entries[0].setText('%.4g' % (float(self.entries[0].text()) + multiplier*lorincr * self.lorstep))
        self.entries[1].setText('%.4g' % (float(self.entries[1].text()) + multiplier*gaussincr * self.gaussstep))
        if (lorincr != 0) and (not self.ticks[0].isChecked()):
            self.ticks[0].setChecked(1)
        if (gaussincr != 0) and (not self.ticks[1].isChecked()):
            self.ticks[1].setChecked(1)
        self.apodPreview()

    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        lor = None
        gauss = None
        cos2 = None
        hamming = None
        shifting = None
        shiftingAxes = 0
        if self.ticks[0].isChecked():
            lor = safeEval(self.entries[0].text())
            if lor is None:
                self.father.father.dispMsg('Lorentzian value is not valid!')
                return
        if self.ticks[1].isChecked():
            gauss = safeEval(self.entries[1].text())
            if gauss is None:
                self.father.father.dispMsg('Gaussian value is not valid!')
                return
        if self.ticks[2].isChecked():
            cos2 = safeEval(self.entries[2].text())
            if cos2 is None:
                self.father.father.dispMsg('cos^2 value is not valid!')
                return
        if self.ticks[3].isChecked():
            hamming = safeEval(self.entries[3].text())
            if hamming is None:
                self.father.father.dispMsg('Hamming value is not valid!')
                return
        shift = safeEval(self.shiftEntry.text())
        if shift is None:
            self.father.father.dispMsg('Shift value is not valid!')
            return
        if self.father.current.data.data.ndim > 1:
            shifting = safeEval(self.shiftingEntry.text())
            if shifting is None:
                self.father.father.dispMsg('Shifting value is not valid!')
                return
            shiftingAxes = int(self.shiftingValues[self.shiftingAxes.currentIndex()]) - 1
        else:
            shiftingAxes = None
        self.father.redoList = []
        self.father.undoList.append(self.father.current.applyApod(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, (self.singleSlice.isChecked())))
        self.father.menuEnable()
        self.deleteLater()

#######################################################################################


class SizeWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        parent.menuDisable()
        self.father = parent
        self.setWindowTitle("Set size")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Size:"), 0, 0)
        self.sizeVal = parent.current.data1D.shape[-1]
        self.sizeEntry = QtWidgets.QLineEdit()
        self.sizeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.sizeEntry.setText(str(self.sizeVal))
        self.sizeEntry.returnPressed.connect(self.sizePreview)

        grid.addWidget(self.sizeEntry, 1, 0)
        grid.addWidget(wc.QLabel("Offset:"), 2, 0)
        if self.father.current.wholeEcho:
            self.posVal = int(np.floor(parent.current.data1D.shape[-1] / 2.0))
        else:
            self.posVal = parent.current.data1D.shape[-1]
        self.posEntry = QtWidgets.QLineEdit()
        self.posEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.posEntry.setText(str(self.posVal))
        self.posEntry.returnPressed.connect(self.sizePreview)
        grid.addWidget(self.posEntry, 3, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 1, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        if not self.father.current.spec:
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
            self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def sizePreview(self, *args):
        inp = safeEval(self.sizeEntry.text())
        if inp is not None:
            self.sizeVal = int(round(inp))
        if self.sizeVal < 1:
            self.father.father.dispMsg('Value is not valid')
            return
        self.sizeEntry.setText(str(self.sizeVal))
        inp = safeEval(self.posEntry.text())
        if inp is not None:
            self.posVal = int(round(inp))
        self.posEntry.setText(str(self.posVal))
        self.father.current.setSizePreview(self.sizeVal, self.posVal)

    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        inp = safeEval(self.sizeEntry.text())
        if inp is not None:
            self.sizeVal = int(round(inp))
        if self.sizeVal < 1:
            self.father.father.dispMsg('Value is not valid')
            return
        inp = safeEval(self.posEntry.text())
        if inp is not None:
            self.posVal = int(round(inp))
        self.father.redoList = []
        self.father.undoList.append(self.father.current.applySize(self.sizeVal, self.posVal))
        self.father.sideframe.upd()
        self.father.menuEnable()
        self.deleteLater()

    def picked(self, pos):
        self.posEntry.setText(str(pos[0]))
        self.sizePreview()
        self.father.current.peakPick = True
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)

##########################################################################################


class SwapEchoWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Swap echo")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Echo position:"), 0, 0)
        self.posVal = int(round(0.5 * len(parent.current.data1D)))
        self.posEntry = QtWidgets.QLineEdit()
        self.posEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.posEntry.setText(str(self.posVal))
        self.posEntry.returnPressed.connect(self.swapEchoPreview)
        grid.addWidget(self.posEntry, 1, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 1, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def swapEchoPreview(self, *args):
        inp = safeEval(self.posEntry.text())
        if inp is not None:
            self.posVal = int(round(inp))
        if self.posVal > 0 and self.posVal < (self.father.current.data1D.shape[-1]):
            self.father.current.setSwapEchoPreview(self.posVal)
            self.father.current.peakPick = False

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        self.father.current.peakPickReset()
        inp = safeEval(self.posEntry.text())
        if inp is not None:
            self.posVal = int(round(inp))
        if self.posVal > 0 and self.posVal < (self.father.current.data1D.shape[-1]):
            self.father.redoList = []
            self.father.undoList.append(self.father.current.applySwapEcho(self.posVal))
            self.father.bottomframe.upd()
            self.father.menuEnable()
            self.deleteLater()
        else:
            self.father.father.dispMsg("not a valid index for swap echo")

    def picked(self, pos):
        self.father.current.setSwapEchoPreview(pos[0])
        self.posEntry.setText(str(pos[0]))
        self.father.current.peakPick = False

###########################################################################


class LPSVDWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        parent.menuDisable()
        self.father = parent
        self.setWindowTitle("LPSVD")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        # grid.addWidget(wc.QLabel("# points for analysis:"), 2, 0)
        self.specGroup = QtWidgets.QButtonGroup(self)
        # self.specGroup.buttonClicked.connect(self.changeSpec)
        backwardButton = QtWidgets.QRadioButton('Backward', parent=self)
        self.specGroup.addButton(backwardButton, 1)
        forwardButton = QtWidgets.QRadioButton('Forward', parent=self)
        self.specGroup.addButton(forwardButton, 0)
        grid.addWidget(backwardButton, 1, 0)
        grid.addWidget(forwardButton, 2, 0)
        backwardButton.setChecked(True)

        grid.addWidget(wc.QLabel("# points for analysis:"), 3, 0)
        self.analPoints = 200
        self.aPointsEntry = QtWidgets.QLineEdit()
        self.aPointsEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.aPointsEntry.setText(str(self.analPoints))
        grid.addWidget(self.aPointsEntry, 4, 0)
        grid.addWidget(wc.QLabel("Number of frequencies:"), 5, 0)
        self.numberFreq = 1
        self.nFreqEntry = QtWidgets.QLineEdit()
        self.nFreqEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.nFreqEntry.setText(str(self.numberFreq))
        grid.addWidget(self.nFreqEntry, 6, 0)

        grid.addWidget(wc.QLabel("Number prediction points:"), 7, 0)
        self.predictPoints = 10
        self.nPredictEntry = QtWidgets.QLineEdit()
        self.nPredictEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.nPredictEntry.setText(str(self.predictPoints))
        grid.addWidget(self.nPredictEntry, 8, 0)

        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 1, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        self.analPoints = safeEval(self.aPointsEntry.text())
        self.numberFreq = safeEval(self.nFreqEntry.text())
        self.predictPoints = safeEval(self.nPredictEntry.text())

        if self.analPoints > len(self.father.current.data1D):
            self.father.father.dispMsg('Number of points for analysis cannot be more than data size!')
            return
        if self.analPoints <= self.numberFreq * 4:
            self.father.father.dispMsg('Number of points for analysis must be more than 4 times the number of frequencies!')
            return

        self.father.redoList = []
        self.father.undoList.append(self.father.current.applyLPSVD(self.analPoints, self.numberFreq, self.predictPoints, self.specGroup.checkedId()))
        self.father.sideframe.upd()
        self.father.menuEnable()
        self.deleteLater()

###########################################################################


class ShiftDataWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Shifting data")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Data points to shift:"), 0, 0, 1, 3)
        self.shiftVal = 0
        self.shiftEntry = QtWidgets.QLineEdit()
        self.shiftEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.shiftEntry.setText(str(self.shiftVal))
        self.shiftEntry.returnPressed.connect(self.shiftPreview)
        grid.addWidget(self.shiftEntry, 1, 1)
        leftShift = QtWidgets.QPushButton("<")
        leftShift.clicked.connect(self.stepDownShift)
        leftShift.setAutoRepeat(True)
        grid.addWidget(leftShift, 1, 0)
        rightShift = QtWidgets.QPushButton(">")
        rightShift.clicked.connect(self.stepUpShift)
        rightShift.setAutoRepeat(True)
        grid.addWidget(rightShift, 1, 2)
        self.singleSlice = QtWidgets.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice, 1, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def stepUpShift(self, *args):
        inp = safeEval(self.shiftEntry.text())
        if inp is not None:
            self.shiftVal = int(round(inp))
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            shift = +10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            shift = +100
        else:
            shift = +1
        self.shiftVal = self.shiftVal + shift
        self.shiftEntry.setText(str(self.shiftVal))
        self.shiftPreview()

    def stepDownShift(self, *args):
        inp = safeEval(self.shiftEntry.text())
        if inp is not None:
            self.shiftVal = int(round(inp))
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            shift = -10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            shift = -100
        else:
            shift = -1
        self.shiftVal = self.shiftVal + shift
        self.shiftEntry.setText(str(self.shiftVal))
        self.shiftPreview()

    def shiftPreview(self, *args):
        inp = safeEval(self.shiftEntry.text())
        if inp is not None:
            self.shiftVal = int(round(inp))
        self.shiftEntry.setText(str(self.shiftVal))
        self.father.current.setShiftPreview(self.shiftVal)

    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        inp = safeEval(self.shiftEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid input")
            return
        shift = int(round(inp))
        self.father.redoList = []
        self.father.undoList.append(self.father.current.applyShift(shift, (self.singleSlice.isChecked())))
        self.father.menuEnable()
        self.deleteLater()

#############################################################


class DCWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Offset correction")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        self.startVal = int(round(0.8 * parent.current.data1D.shape[-1]))
        self.endVal = parent.current.data1D.shape[-1]
        grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.startEntry = QtWidgets.QLineEdit()
        self.startEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.startEntry.setText(str(self.startVal))
        self.startEntry.returnPressed.connect(self.offsetPreview)
        grid.addWidget(self.startEntry, 1, 0)
        grid.addWidget(wc.QLabel("End point:"), 2, 0)
        self.endEntry = QtWidgets.QLineEdit()
        self.endEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.endEntry.setText(str(self.endVal))
        self.endEntry.returnPressed.connect(self.offsetPreview)
        grid.addWidget(self.endEntry, 3, 0)
        grid.addWidget(wc.QLabel("Offset:"), 4, 0)
        self.offsetEntry = QtWidgets.QLineEdit()
        self.offsetEntry.setAlignment(QtCore.Qt.AlignHCenter)
        val = parent.current.getdcOffset(int(round(0.8 * parent.current.data1D.shape[-1])), parent.current.data1D.shape[-1])
        self.offsetEntry.setText('{:.2e}'.format(val))
        self.offsetEntry.returnPressed.connect(lambda: self.offsetPreview(True))
        grid.addWidget(self.offsetEntry, 5, 0)
        self.singleSlice = QtWidgets.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice, 1, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def picked(self, pos, second=False):
        dataLength = self.father.current.data1D.shape[-1]
        if second:
            inp = safeEval(self.startEntry.text())
            if inp is not None:
                self.startVal = int(round(inp))
            if self.startVal < 0:
                self.startVal = 0
            elif self.startVal > dataLength:
                self.startVal = dataLength
            self.startEntry.setText(str(self.startVal))
            self.endVal = pos[0]
            self.endEntry.setText(str(self.endVal))
            val = self.father.current.getdcOffset(self.startVal, self.endVal)
            self.offsetEntry.setText('{:.2e}'.format(val))
            self.father.current.dcOffset(val)
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
            self.father.current.peakPick = True
        else:
            self.startEntry.setText(str(pos[0]))
            inp = safeEval(self.endEntry.text())
            if inp is not None:
                self.endVal = int(round(inp))
            if self.endVal < 0:
                self.endVal = 0
            elif self.endVal > dataLength:
                self.endVal = dataLength
            self.startVal = pos[0]
            val = self.father.current.getdcOffset(self.startVal, self.endVal)
            self.offsetEntry.setText('{:.2e}'.format(val))
            self.father.current.peakPickFunc = lambda pos, self= self: self.picked(pos, True)
            self.father.current.peakPick = True

    def offsetPreview(self, inserted=False):
        if inserted:
            self.father.current.dcOffset(safeEval(self.offsetEntry.text()))
        else:
            dataLength = self.father.current.data1D.shape[-1]
            inp = safeEval(self.startEntry.text())
            if inp is not None:
                self.startVal = int(round(inp))
            if self.startVal < 0:
                self.startVal = 0
            elif self.startVal > dataLength:
                self.startVal = dataLength
            self.startEntry.setText(str(self.startVal))
            inp = safeEval(self.endEntry.text())
            if inp is not None:
                self.endVal = int(round(inp))
            if self.endVal < 0:
                self.endVal = 0
            elif self.endVal > dataLength:
                self.endVal = dataLength
            self.endEntry.setText(str(self.endVal))
            val = self.father.current.getdcOffset(self.startVal, self.endVal)
            self.offsetEntry.setText('{:.2e}'.format(val))
            self.father.current.dcOffset(val)

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        inp = safeEval(self.offsetEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        self.father.current.peakPickReset()
        self.father.redoList = []
        self.father.undoList.append(self.father.current.subtract(inp, self.singleSlice.isChecked()))
        self.father.menuEnable()
        self.deleteLater()

#############################################################


class BaselineWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Baseline correction")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Polynomial Degree:"), 0, 0, 1, 2)
        self.removeList = []
        self.degree = 3
        self.degreeEntry = QtWidgets.QLineEdit()
        self.degreeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.degreeEntry.setText(str(self.degree))
        self.degreeEntry.returnPressed.connect(self.setDegree)
        grid.addWidget(self.degreeEntry, 1, 0, 1, 2)
        resetButton = QtWidgets.QPushButton("&Reset")
        resetButton.clicked.connect(self.reset)
        grid.addWidget(resetButton, 2, 0)
        fitButton = QtWidgets.QPushButton("&Fit")
        fitButton.clicked.connect(self.preview)
        grid.addWidget(fitButton, 2, 1)
        self.singleSlice = QtWidgets.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice, 1, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def picked(self, pos):
        self.removeList.append(pos[0])
        self.father.current.previewRemoveList(self.removeList)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def setDegree(self):
        inp = safeEval(self.degreeEntry.text())
        if inp is not None:
            self.degree = inp
        self.degreeEntry.setText(str(self.degree))

    def preview(self, *args):
        inp = safeEval(self.degreeEntry.text())
        if inp is not None:
            self.degree = inp
        self.degreeEntry.setText(str(self.degree))
        self.father.current.previewBaseline(self.degree, self.removeList)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def reset(self, *args):
        self.removeList = []
        self.father.current.resetPreviewRemoveList()

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.current.resetPreviewRemoveList()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        inp = safeEval(self.degreeEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        returnValue = self.father.current.applyBaseline(self.degree, self.removeList, self.singleSlice.isChecked())
        if returnValue is None:
            return
        self.father.current.peakPickReset()
        self.father.current.resetPreviewRemoveList()
        self.father.redoList = []
        self.father.undoList.append(returnValue)
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

#############################################################


class regionWindow(QtWidgets.QWidget):

    def __init__(self, parent, name):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle(name)
        layout = QtWidgets.QGridLayout(self)
        self.grid = QtWidgets.QGridLayout()
        layout.addLayout(self.grid, 0, 0, 1, 2)
        self.startVal = [0]  # dummy variables
        self.endVal = [parent.current.data1D.shape[-1]]  # dummy variables
        self.grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.grid.addWidget(wc.QLabel("End point:"), 0, 1)
        self.startEntry = []
        self.endEntry = []
        self.deleteButton = []
        self.partIter = 0
        self.entryCount = 1
        self.first = True
        self.startEntry.append(QtWidgets.QLineEdit())
        self.startEntry[0].setAlignment(QtCore.Qt.AlignHCenter)
        self.startEntry[0].setText("")
        self.startEntry[0].editingFinished.connect(lambda self=self, tmp=self.startEntry[0]: self.setVal(tmp, True))
        self.grid.addWidget(self.startEntry[0], 1, 0)
        self.endEntry.append(QtWidgets.QLineEdit())
        self.endEntry[0].setAlignment(QtCore.Qt.AlignHCenter)
        self.endEntry[0].setText("")
        self.endEntry[0].editingFinished.connect(lambda self=self, tmp=self.endEntry[0]: self.setVal(tmp, False))
        self.grid.addWidget(self.endEntry[0], 1, 1)
        self.deleteButton.append(QtWidgets.QPushButton("X"))
        self.deleteButton[0].clicked.connect(lambda extra, self=self: self.deleteEntry(self.deleteButton[0]))
        self.grid.addWidget(self.deleteButton[0], 1, 2)
        self.newSpec = QtWidgets.QCheckBox("Result in new workspace")
        layout.addWidget(self.newSpec, 1, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok", self)
        okButton.clicked.connect(self.applyAndClose)
        okButton.setFocus()
        layout.addWidget(okButton, 2, 1)
        self.grid.setRowStretch(100, 1)
        self.show()
        # self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

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
            self.endVal = np.append(self.endVal, self.father.current.data1D.shape[-1])
            self.startEntry[self.partIter].setText(str(self.startVal[self.partIter]))
            self.endEntry[self.partIter].setText(str(self.endVal[self.partIter]))
            self.partIter += 1
            self.startEntry.append(QtWidgets.QLineEdit())
            self.startEntry[self.partIter].setAlignment(QtCore.Qt.AlignHCenter)
            self.startEntry[self.partIter].editingFinished.connect(lambda self=self, tmp=self.startEntry[self.partIter]: self.setVal(tmp, True))
            self.grid.addWidget(self.startEntry[self.partIter], 1 + self.entryCount, 0)
            self.endEntry.append(QtWidgets.QLineEdit())
            self.endEntry[self.partIter].setAlignment(QtCore.Qt.AlignHCenter)
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
            self.endVal[num] = self.father.current.data1D.shape[-1]
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
        inp = safeEval(entry.text())
        if inp is None:
            return
        inp = int(inp)
        if inp < 0:
            inp = 0
        if inp > self.father.current.data1D.shape[-1]:
            inp = self.father.current.data1D.shape[-1]
        if isMin:
            num = self.startEntry.index(entry)
            self.startVal[num] = min(inp, self.endVal[num])
            self.endVal[num] = max(inp, self.endVal[num])
        else:
            num = self.endEntry.index(entry)
            self.endVal[num] = max(inp, self.startVal[num])
            self.startVal[num] = min(inp, self.startVal[num])
        if num == self.partIter:
            self.partIter += 1
            self.startVal = np.append(self.startVal, 0)
            self.endVal = np.append(self.endVal, self.father.current.data1D.shape[-1])
            self.startEntry.append(QtWidgets.QLineEdit())
            self.startEntry[self.partIter].setAlignment(QtCore.Qt.AlignHCenter)
            self.startEntry[self.partIter].editingFinished.connect(lambda self=self, tmp=self.startEntry[self.partIter]: self.setVal(tmp, True))
            self.grid.addWidget(self.startEntry[self.partIter], 1 + self.entryCount, 0)
            self.endEntry.append(QtWidgets.QLineEdit())
            self.endEntry[self.partIter].setAlignment(QtCore.Qt.AlignHCenter)
            self.endEntry[self.partIter].editingFinished.connect(lambda self=self, tmp=self.endEntry[self.partIter]: self.setVal(tmp, False))
            self.grid.addWidget(self.endEntry[self.partIter], 1 + self.entryCount, 1)
            self.deleteButton.append(QtWidgets.QPushButton("X"))
            self.deleteButton[self.partIter].clicked.connect(lambda extra, self=self, tmp=self.deleteButton[self.partIter]: self.deleteEntry(tmp))
            self.grid.addWidget(self.deleteButton[self.partIter], 1 + self.entryCount, 2)
            self.entryCount += 1
            self.first = True
        self.startEntry[num].setText(str(self.startVal[num]))
        self.endEntry[num].setText(str(self.endVal[num]))

    def apply(self, maximum, minimum, newSpec):
        pass

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.updAllFrames()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        if self.partIter == 0:
            if self.apply(np.array([0]), np.array([self.father.current.data1D.shape[-1]]), self.newSpec.isChecked()) is None:
                return
        else:
            if self.apply(self.startVal[:self.partIter], self.endVal[:self.partIter], self.newSpec.isChecked()) is None:
                return
        self.father.current.peakPickReset()
        self.father.menuEnable()
        self.deleteLater()

############################################################


class integrateWindow(regionWindow):

    def __init__(self, parent):
        regionWindow.__init__(self, parent, 'Integrate')

    def apply(self, maximum, minimum, newSpec):
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.integrate(minimum, maximum, newSpec)) is None:
                return None
        else:
            returnValue = self.father.current.integrate(minimum, maximum, newSpec)
            if returnValue is None:
                return None
            self.father.redoList = []
            self.father.undoList.append(returnValue)
            self.father.updAllFrames()
        return 1

############################################################


class sumWindow(regionWindow):

    def __init__(self, parent):
        regionWindow.__init__(self, parent, 'Sum')

    def apply(self, maximum, minimum, newSpec):
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.sum(minimum, maximum, newSpec)) is None:
                return None
        else:
            returnValue = self.father.current.sum(minimum, maximum, newSpec)
            if returnValue is None:
                return None
            self.father.redoList = []
            self.father.undoList.append(returnValue)
            self.father.updAllFrames()
        return 1

############################################################


class maxWindow(regionWindow):

    def __init__(self, parent):
        regionWindow.__init__(self, parent, 'Max')

    def apply(self, maximum, minimum, newSpec):
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.maxMatrix(minimum, maximum, newSpec)) is None:
                return None
        else:
            returnValue = self.father.current.maxMatrix(minimum, maximum, newSpec)
            if returnValue is None:
                return None
            self.father.redoList = []
            self.father.undoList.append(returnValue)
            self.father.updAllFrames()
        return 1

############################################################


class minWindow(regionWindow):

    def __init__(self, parent):
        regionWindow.__init__(self, parent, 'Min')

    def apply(self, maximum, minimum, newSpec):
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.minMatrix(minimum, maximum, newSpec)) is None:
                return None
        else:
            returnValue = self.father.current.minMatrix(minimum, maximum, newSpec)
            if returnValue is None:
                return None
            self.father.redoList = []
            self.father.undoList.append(returnValue)
            self.father.updAllFrames()
        return 1

############################################################


class argmaxWindow(regionWindow):

    def __init__(self, parent):
        regionWindow.__init__(self, parent, 'Max position')

    def apply(self, maximum, minimum, newSpec):
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.argmaxMatrix(minimum, maximum, newSpec)) is None:
                return None
        else:
            returnValue = self.father.current.argmaxMatrix(minimum, maximum, newSpec)
            if returnValue is None:
                return None
            self.father.redoList = []
            self.father.undoList.append(returnValue)
            self.father.updAllFrames()
        return 1

############################################################


class argminWindow(regionWindow):

    def __init__(self, parent):
        regionWindow.__init__(self, parent, 'Min position')

    def apply(self, maximum, minimum, newSpec):
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.argminMatrix(minimum, maximum, newSpec)) is None:
                return None
        else:
            returnValue = self.father.current.argminMatrix(minimum, maximum, newSpec)
            if returnValue is None:
                return None
            self.father.redoList = []
            self.father.undoList.append(returnValue)
            self.father.updAllFrames()
        return 1

############################################################


class avgWindow(regionWindow):

    def __init__(self, parent):
        regionWindow.__init__(self, parent, 'Average')

    def apply(self, maximum, minimum, newSpec):
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.average(minimum, maximum, newSpec)) is None:
                return None
        else:
            returnValue = self.father.current.average(minimum, maximum, newSpec)
            if returnValue is None:
                return None
            self.father.redoList = []
            self.father.undoList.append(returnValue)
            self.father.updAllFrames()
        return 1

#############################################################


class regionWindow2(QtWidgets.QWidget):

    def __init__(self, parent, name, newSpecOption):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle(name)
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        self.startVal = 0
        self.endVal = parent.current.data1D.shape[-1]
        grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.startEntry = QtWidgets.QLineEdit()
        self.startEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.startEntry.setText(str(self.startVal))
        self.startEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.startEntry, 1, 0)
        grid.addWidget(wc.QLabel("End point:"), 2, 0)
        self.endEntry = QtWidgets.QLineEdit()
        self.endEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.endEntry.setText(str(self.endVal))
        self.endEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.endEntry, 3, 0)
        self.newSpec = QtWidgets.QCheckBox("Result in new workspace")
        if not newSpecOption:
            self.newSpec.hide()
        layout.addWidget(self.newSpec, 1, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def preview(self, maximum, minimum):
        pass

    def picked(self, pos, second=False):
        if second:
            dataLength = self.father.current.data1D.shape[-1]
            inp = safeEval(self.startEntry.text())
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
            self.preview(self.startVal, self.endVal)
        else:
            self.startEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, True)
            self.father.current.peakPick = True

    def checkValues(self, *args):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.startEntry.text())
        if inp is not None:
            self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        self.startEntry.setText(str(self.startVal))
        inp = safeEval(self.endEntry.text())
        if inp is not None:
            self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        self.endEntry.setText(str(self.endVal))
        self.preview(self.startVal, self.endVal)

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.current.upd()
        self.father.current.showFid()
        self.father.updAllFrames()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.startEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        inp = safeEval(self.endEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        if self.apply(self.startVal, self.endVal, self.newSpec.isChecked()) is None:
            return
        self.father.current.peakPickReset()
        self.father.menuEnable()
        self.deleteLater()

    def apply(self, maximum, minimum, newSpec):
        pass

############################################################


class extractRegionWindow(regionWindow2):

    def __init__(self, parent):
        regionWindow2.__init__(self, parent, 'Extract part', True)

    def apply(self, maximum, minimum, newSpec):
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.getRegion(minimum, maximum, newSpec)) is None:
                return None
        else:
            returnValue = self.father.current.getRegion(minimum, maximum, newSpec)
            if returnValue is None:
                return None
            self.father.redoList = []
            self.father.undoList.append(returnValue)
            self.father.updAllFrames()
        return 1

############################################################


class SubtractAvgWindow(regionWindow2):

    def __init__(self, parent):
        regionWindow2.__init__(self, parent, 'Subtract Avg', False)

    def apply(self, maximum, minimum, newSpec):
        returnValue = self.father.current.subtractAvg(maximum, minimum)
        if returnValue is None:
            return None
        self.father.redoList = []
        self.father.undoList.append(returnValue)
        self.father.updAllFrames()
        return 1

    def preview(self, maximum, minimum):
        self.father.current.subtractAvgPreview(maximum, minimum)

#############################################################


class FiddleWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Reference deconvolution")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        self.startVal = 0
        self.endVal = parent.current.data1D.shape[-1]
        grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.startEntry = QtWidgets.QLineEdit()
        self.startEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.startEntry.setText(str(self.startVal))
        self.startEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.startEntry, 1, 0)
        grid.addWidget(wc.QLabel("End point:"), 2, 0)
        self.endEntry = QtWidgets.QLineEdit()
        self.endEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.endEntry.setText(str(self.endVal))
        self.endEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.endEntry, 3, 0)
        grid.addWidget(wc.QLabel("Linebroadening [Hz]:"), 4, 0)
        self.lbEntry = QtWidgets.QLineEdit()
        self.lbEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.lbEntry.setText("1.0")
        self.lbEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.lbEntry, 5, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def picked(self, pos, second=False):
        if second:
            dataLength = self.father.current.data1D.shape[-1]
            inp = safeEval(self.startEntry.text())
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
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.startEntry.text())
        if inp is not None:
            self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        self.startEntry.setText(str(self.startVal))
        inp = safeEval(self.endEntry.text())
        if inp is not None:
            self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        self.endEntry.setText(str(self.endVal))
        inp = safeEval(self.lbEntry.text())
        if inp is not None:
            self.lbEntry.setText(str(inp))

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.current.upd()
        self.father.current.showFid()
        self.father.updAllFrames()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.startEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        inp = safeEval(self.endEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        lb = safeEval(self.lbEntry.text())
        if lb is None:
            self.father.father.dispMsg("Not a valid value")
            return
        returnValue = self.father.current.fiddle(self.startVal, self.endVal, lb)
        if returnValue is None:
            return None
        self.father.redoList = []
        self.father.undoList.append(returnValue)
        self.father.updAllFrames()
        self.father.current.peakPickReset()
        self.father.menuEnable()
        self.deleteLater()

##############################################################


class DeleteWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Delete")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Indexes to delete:"), 0, 0)
        self.delEntry = QtWidgets.QLineEdit()
        self.delEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.delEntry.setText(str('0'))
        self.delEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.delEntry, 1, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def preview(self, *args):
        env = vars(np).copy()
        length = int(self.father.current.data1D.shape[-1])
        env['length'] = length  # so length can be used to in equations
        pos = np.array(eval(self.delEntry.text(), env)).flatten()                # find a better solution, also add catch for exceptions
        pos[pos < 0] = pos[pos < 0] + length
        if (pos > -1).all() and (pos < length).all():
            self.father.current.deletePreview(pos)
        else:
            self.father.father.dispMsg('Not all values are valid indexes to delete')

    def applyAndClose(self):
        env = vars(np).copy()
        length = int(self.father.current.data1D.shape[-1])
        env['length'] = length  # so length can be used to in equations
        pos = np.array(eval(self.delEntry.text(), env)).flatten()                # find a better solution, also add catch for exceptions
        pos[pos < 0] = pos[pos < 0] + length
        if (pos > -1).all() and (pos < length).all():
            self.father.redoList = []
            self.father.undoList.append(self.father.current.delete(pos))
            self.father.menuEnable()
            self.father.sideframe.upd()
            self.deleteLater()
        else:
            self.father.father.dispMsg('Not all values are valid indexes to delete')

    def closeEvent(self, *args):
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

##############################################################


class SplitWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Split")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Sections:"), 0, 0)
        self.splitEntry = QtWidgets.QLineEdit()
        self.splitEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.splitEntry.setText('1')
        self.splitEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.splitEntry, 1, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def preview(self, *args):
        val = safeEval(self.splitEntry.text(), self.father.current.data1D.shape[-1])
        if val is not None:
            self.splitEntry.setText(str(int(round(val))))

    def applyAndClose(self):
        val = safeEval(self.splitEntry.text(), self.father.current.data1D.shape[-1])
        if val is None:
            self.father.father.dispMsg("Not a valid value")
            return
        returnValue = self.father.current.split(int(round(val)))
        if returnValue is None:
            return
        self.father.redoList = []
        self.father.undoList.append(returnValue)
        self.father.menuEnable()
        self.father.updAllFrames()
        self.deleteLater()

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

##############################################################


class ConcatenateWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Concatenate")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Concatenation axes:"), 0, 0)
        self.axesEntry = QtWidgets.QComboBox()
        self.axesEntry.addItems(np.array(np.arange(self.father.current.data.data.ndim - 1) + 1, dtype=str))
        grid.addWidget(self.axesEntry, 1, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def applyAndClose(self):
        returnValue = self.father.current.concatenate(self.axesEntry.currentIndex())
        if returnValue is None:
            return
        self.father.redoList = []
        self.father.undoList.append(returnValue)
        self.father.menuEnable()
        self.father.updAllFrames()
        self.deleteLater()

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

##############################################################


class InsertWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Insert")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Start insert at index:"), 0, 0)
        self.posEntry = QtWidgets.QLineEdit()
        self.posEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.posEntry.setText(str(self.father.current.data1D.shape[-1]))
        self.posEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.posEntry, 1, 0)
        grid.addWidget(wc.QLabel("Workspace to insert:"), 2, 0)
        self.wsEntry = QtWidgets.QComboBox()
        self.wsEntry.addItems(self.father.father.workspaceNames)
        grid.addWidget(self.wsEntry, 3, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def preview(self, *args):
        pos = safeEval(self.posEntry.text())
        if pos is None:
            return
        pos = int(round(pos))
        if pos > self.father.current.data1D.shape[-1]:
            pos = self.father.current.data1D.shape[-1]
        elif pos < 0:
            pos = 0
        self.posEntry.setText(str(pos))

    def applyAndClose(self):
        pos = safeEval(self.posEntry.text())
        if pos is None:
            self.father.father.dispMsg("Not a valid value")
            return
        pos = int(round(pos))
        if pos > self.father.current.data1D.shape[-1]:
            pos = self.father.current.data1D.shape[-1]
        elif pos < 0:
            pos = 0
        ws = self.wsEntry.currentIndex()
        self.father.redoList = []
        self.father.undoList.append(self.father.current.insert(self.father.father.workspaces[ws].masterData.data, pos))
        self.father.menuEnable()
        self.father.sideframe.upd()
        self.deleteLater()

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

##############################################################


class AddWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Add")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Workspace to add:"), 0, 0)
        self.wsEntry = QtWidgets.QComboBox()
        self.wsEntry.addItems(self.father.father.workspaceNames)
        grid.addWidget(self.wsEntry, 1, 0)
        self.singleSlice = QtWidgets.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice, 1, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def applyAndClose(self):
        ws = self.wsEntry.currentIndex()
        returnValue = self.father.current.add(self.father.father.workspaces[ws].masterData.data, self.singleSlice.isChecked())
        if returnValue is None:
            return
        self.father.redoList = []
        self.father.undoList.append(returnValue)
        self.father.menuEnable()
        self.father.sideframe.upd()
        self.deleteLater()

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

##############################################################


class SubtractWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Subtract")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Workspace to subtract:"), 0, 0)
        self.wsEntry = QtWidgets.QComboBox()
        self.wsEntry.addItems(self.father.father.workspaceNames)
        grid.addWidget(self.wsEntry, 1, 0)
        self.singleSlice = QtWidgets.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice, 1, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def applyAndClose(self):
        ws = self.wsEntry.currentIndex()
        returnValue = self.father.current.subtract(self.father.father.workspaces[ws].masterData.data, self.singleSlice.isChecked())
        if returnValue is None:
            return
        self.father.redoList = []
        self.father.undoList.append(returnValue)
        self.father.menuEnable()
        self.father.sideframe.upd()
        self.deleteLater()

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

##############################################################


class SNWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Signal to noise")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Start point noise:"), 0, 0)
        self.minNoiseEntry = QtWidgets.QLineEdit()
        self.minNoiseEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.minNoiseEntry.setText("0")
        self.minNoiseEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.minNoiseEntry, 1, 0)
        grid.addWidget(wc.QLabel("End point noise:"), 2, 0)
        self.maxNoiseEntry = QtWidgets.QLineEdit()
        self.maxNoiseEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.maxNoiseEntry.setText(str(parent.current.data1D.shape[-1]))
        self.maxNoiseEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.maxNoiseEntry, 3, 0)
        grid.addWidget(wc.QLabel("Start point signal:"), 4, 0)
        self.minEntry = QtWidgets.QLineEdit()
        self.minEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.minEntry.setText("0")
        self.minEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.minEntry, 5, 0)
        grid.addWidget(wc.QLabel("End point signal:"), 6, 0)
        self.maxEntry = QtWidgets.QLineEdit()
        self.maxEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.maxEntry.setText(str(parent.current.data1D.shape[-1]))
        self.maxEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.maxEntry, 7, 0)
        grid.addWidget(wc.QLabel("S/N:"), 8, 0)
        self.snEntry = QtWidgets.QLineEdit()
        self.snEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.snEntry.setText('0.0')
        grid.addWidget(self.snEntry, 9, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Fit")
        okButton.clicked.connect(self.apply)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

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
            self.apply()

    def checkValues(self, *args):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.minNoiseEntry.text())
        if inp is None:
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minNoiseEntry.setText(str(minimum))
        inp = safeEval(self.maxNoiseEntry.text())
        if inp is None:
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxNoiseEntry.setText(str(maximum))
        inp = safeEval(self.minEntry.text())
        if inp is None:
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.apply()

    def apply(self):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.minNoiseEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        minimumNoise = int(round(inp))
        if minimumNoise < 0:
            minimumNoise = 0
        elif minimumNoise > dataLength:
            minimumNoise = dataLength
        self.minNoiseEntry.setText(str(minimumNoise))
        inp = safeEval(self.maxNoiseEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        maximumNoise = int(round(inp))
        if maximumNoise < 0:
            maximumNoise = 0
        elif maximumNoise > dataLength:
            maximumNoise = dataLength
        self.maxNoiseEntry.setText(str(maximumNoise))
        inp = safeEval(self.minEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.snEntry.setText(str(self.father.current.SN(minimumNoise, maximumNoise, minimum, maximum)))

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.menuEnable()
        self.deleteLater()

##############################################################


class FWHMWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("FWHM")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.minEntry = QtWidgets.QLineEdit()
        self.minEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.minEntry.setText("0")
        self.minEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.minEntry, 1, 0)
        grid.addWidget(wc.QLabel("End point:"), 2, 0)
        self.maxEntry = QtWidgets.QLineEdit()
        self.maxEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.maxEntry.setText(str(parent.current.data1D.shape[-1]))
        self.maxEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.maxEntry, 3, 0)
        if self.father.current.spec == 1:
            if self.father.current.ppm:
                grid.addWidget(wc.QLabel("FWHM [ppm]:"), 4, 0)
            else:
                if self.father.current.axType == 0:
                    grid.addWidget(wc.QLabel("FWHM [Hz]:"), 4, 0)
                elif self.father.current.axType == 1:
                    grid.addWidget(wc.QLabel("FWHM [kHz]:"), 4, 0)
                elif self.father.current.axType == 2:
                    grid.addWidget(wc.QLabel("FWHM [MHz]:"), 4, 0)
                elif self.father.current.axType == 3:
                    grid.addWidget(wc.QLabel("FWHM [ppm]:"), 4, 0)
        else:
            if self.father.current.axType == 0:
                grid.addWidget(wc.QLabel("FWHM [s]:"), 4, 0)
            elif self.father.current.axType == 1:
                grid.addWidget(wc.QLabel("FWHM [ms]:"), 4, 0)
            elif self.father.current.axType == 2:
                grid.addWidget(wc.QLabel(u"FWHM [\u03bcs]:"), 4, 0)
        self.fwhmEntry = QtWidgets.QLineEdit()
        self.fwhmEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.fwhmEntry.setText('0.0')
        grid.addWidget(self.fwhmEntry, 5, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Fit")
        okButton.clicked.connect(self.apply)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def picked(self, pos, num=0):
        if num == 0:
            self.minEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 1)
            self.father.current.peakPick = True
        elif num == 1:
            self.maxEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 0)
            self.father.current.peakPick = True
            self.apply()

    def checkValues(self, *args):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.minEntry.text())
        if inp is None:
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.apply()

    def apply(self):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.minEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.fwhmEntry.setText(str(self.father.current.fwhm(minimum, maximum)))

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.menuEnable()
        self.deleteLater()
##############################################################


class COMWindow(QtWidgets.QWidget):  # Centre of Mass Window

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Centre of Mass")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.minEntry = QtWidgets.QLineEdit()
        self.minEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.minEntry.setText("0")
        self.minEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.minEntry, 1, 0)
        grid.addWidget(wc.QLabel("End point:"), 2, 0)
        self.maxEntry = QtWidgets.QLineEdit()
        self.maxEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.maxEntry.setText(str(parent.current.data1D.shape[-1]))
        self.maxEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.maxEntry, 3, 0)
        if self.father.current.spec == 1:
            if self.father.current.ppm:
                grid.addWidget(wc.QLabel("Centre of Mass [ppm]:"), 4, 0)
            else:
                if self.father.current.axType == 0:
                    grid.addWidget(wc.QLabel("Centre of Mass [Hz]:"), 4, 0)
                elif self.father.current.axType == 1:
                    grid.addWidget(wc.QLabel("Centre of Mass [kHz]:"), 4, 0)
                elif self.father.current.axType == 2:
                    grid.addWidget(wc.QLabel("Centre of Mass [MHz]:"), 4, 0)
                elif self.father.current.axType == 3:
                    grid.addWidget(wc.QLabel("Centre of Mass [ppm]:"), 4, 0)
        else:
            if self.father.current.axType == 0:
                grid.addWidget(wc.QLabel("Centre of Mass [s]:"), 4, 0)
            elif self.father.current.axType == 1:
                grid.addWidget(wc.QLabel("Centre of Mass [ms]:"), 4, 0)
            elif self.father.current.axType == 2:
                grid.addWidget(wc.QLabel(u"Centre of Mass [\u03bcs]:"), 4, 0)
        self.comEntry = QtWidgets.QLineEdit()
        self.comEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.comEntry.setText('0.0')
        grid.addWidget(self.comEntry, 5, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Fit")
        okButton.clicked.connect(self.apply)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def picked(self, pos, num=0):
        if num == 0:
            self.minEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 1)
            self.father.current.peakPick = True
        elif num == 1:
            self.maxEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos, 0)
            self.father.current.peakPick = True
            self.apply()

    def checkValues(self, *args):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.minEntry.text())
        if inp is None:
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.apply()

    def apply(self):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.minEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            self.father.father.dispMsg("Not a valid value")
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.comEntry.setText(str(self.father.current.COM(minimum, maximum)))

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.menuEnable()
        self.deleteLater()

##########################################################################################


class ReorderWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Reorder")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Positions of the spectra:"), 0, 0)
        self.valEntry = QtWidgets.QLineEdit()
        self.valEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.valEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.valEntry, 1, 0)
        fileButton = QtWidgets.QPushButton("&Browse")
        fileButton.clicked.connect(self.getPosFromFile)
        grid.addWidget(fileButton, 2, 0)
        grid.addWidget(wc.QLabel("Length of dimension:"), 3, 0)
        self.lengthEntry = QtWidgets.QLineEdit()
        self.lengthEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.lengthEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.lengthEntry, 4, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def preview(self, *args):
        pass

    def getPosFromFile(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.LastLocation)
        if type(filename) is tuple:
            filename = filename[0]        
        if filename:  # if not cancelled
            self.father.father.LastLocation = os.path.dirname(filename)  # Save used path
        if len(filename) == 0:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        newLength = safeEval(self.lengthEntry.text())
        env = vars(np).copy()
        env['length'] = int(self.father.current.data1D.shape[-1])  # so length can be used to in equations
        env['euro'] = lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal, num)
        val = eval(self.valEntry.text(), env)                # find a better solution, also add catch for exceptions
        if not isinstance(val, (list, np.ndarray)):
            self.father.father.dispMsg("Input is not a list or array")
            return
        if len(val) != self.father.current.data1D.shape[-1]:
            self.father.father.dispMsg("Length of input does not match length of data")
            return
        val = np.array(val, dtype=int)
        self.father.redoList = []
        self.father.undoList.append(self.father.current.reorder(val, newLength))
        self.father.menuEnable()
        self.father.updAllFrames()
        self.deleteLater()

##########################################################################################


class FFMWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("FFM")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Positions of the spectra:"), 0, 0)
        self.valEntry = QtWidgets.QLineEdit()
        self.valEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.valEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.valEntry, 1, 0)
        fileButton = QtWidgets.QPushButton("&Browse")
        fileButton.clicked.connect(self.getPosFromFile)
        grid.addWidget(fileButton, 2, 0)
        grid.addWidget(wc.QLabel("Type of the position list:"), 3, 0)
        self.typeDrop = QtWidgets.QComboBox(parent=self)
        self.typeDrop.addItems(["Complex", "States/States-TPPI", "TPPI"])
        grid.addWidget(self.typeDrop, 4, 0)
        grid.addWidget(wc.QLabel("Reconstruction may take a while"), 5, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def preview(self, *args):
        pass

    def getPosFromFile(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.LastLocation)
        if type(filename) is tuple:
            filename = filename[0]        
        if filename:  # if not cancelled
            self.father.father.LastLocation = os.path.dirname(filename)  # Save used path
        if len(filename) == 0:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        env = vars(np).copy()
        env['length'] = int(self.father.current.data1D.shape[-1])  # so length can be used to in equations
        env['euro'] = lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal, num)
        val = eval(self.valEntry.text(), env)                # find a better solution, also add catch for exceptions
        if not isinstance(val, (list, np.ndarray)):
            self.father.father.dispMsg("Input is not a list or array")
            return
        val = np.array(val, dtype=int)
        self.father.redoList = []
        self.father.undoList.append(self.father.current.ffm(val, self.typeDrop.currentIndex()))
        self.father.updAllFrames()
        self.father.menuEnable()
        self.deleteLater()

##########################################################################################


class CLEANWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("CLEAN")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Positions of the spectra:"), 0, 0)
        self.valEntry = QtWidgets.QLineEdit()
        self.valEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.valEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.valEntry, 1, 0)
        fileButton = QtWidgets.QPushButton("&Browse")
        fileButton.clicked.connect(self.getPosFromFile)
        grid.addWidget(fileButton, 2, 0)
        grid.addWidget(wc.QLabel("Type of the position list:"), 3, 0)
        self.typeDrop = QtWidgets.QComboBox(parent=self)
        self.typeDrop.addItems(["Complex", "States/States-TPPI", "TPPI"])
        grid.addWidget(self.typeDrop, 4, 0)
        grid.addWidget(wc.QLabel("Gamma:"), 5, 0)
        self.gammaEntry = QtWidgets.QLineEdit()
        self.gammaEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.gammaEntry.setText("0.2")
        grid.addWidget(self.gammaEntry, 6, 0)
        grid.addWidget(wc.QLabel("Threshold:"), 7, 0)
        self.thresholdEntry = QtWidgets.QLineEdit()
        self.thresholdEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.thresholdEntry.setText("2.0")
        grid.addWidget(self.thresholdEntry, 8, 0)
        # grid.addWidget(wc.QLabel("Linewidth [Hz]:"), 9, 0)
        # self.lbEntry = QtWidgets.QLineEdit()
        # self.lbEntry.setAlignment(QtCore.Qt.AlignHCenter)
        # self.lbEntry.setText("1.0")
        # grid.addWidget(self.lbEntry, 10, 0)
        grid.addWidget(wc.QLabel("Max. iterations:"), 11, 0)
        self.maxIterEntry = QtWidgets.QLineEdit()
        self.maxIterEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.maxIterEntry.setText("2000")
        grid.addWidget(self.maxIterEntry, 12, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def preview(self, *args):
        pass

    def getPosFromFile(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.LastLocation)
        if type(filename) is tuple:
            filename = filename[0]        
        if filename:  # if not cancelled
            self.father.father.LastLocation = os.path.dirname(filename)  # Save used path
        if len(filename) == 0:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        env = vars(np).copy()
        env['length'] = int(self.father.current.data1D.shape[-1])  # so length can be used to in equations
        env['euro'] = lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal, num)
        val = eval(self.valEntry.text(), env)                # find a better solution, also add catch for exceptions
        if not isinstance(val, (list, np.ndarray)):
            self.father.father.dispMsg("Input is not a list or array")
            return
        val = np.array(val, dtype=int)
        gamma = safeEval(self.gammaEntry.text())
        if gamma is None:
            self.father.dispMsg("One of the inputs is not valid")
            return
        threshold = safeEval(self.thresholdEntry.text())
        if threshold is None:
            self.father.dispMsg("One of the inputs is not valid")
            return
        threshold = threshold
        # lb = safeEval(self.lbEntry.text())
        # if lb is None:
        #    self.father.dispMsg("One of the inputs is not valid")
        #    return
        maxIter = safeEval(self.maxIterEntry.text())
        if maxIter is None:
            self.father.dispMsg("One of the inputs is not valid")
            return
        maxIter = int(maxIter)
        self.father.redoList = []
        self.father.undoList.append(self.father.current.clean(val, self.typeDrop.currentIndex(), gamma, threshold, maxIter))
        self.father.updAllFrames()
        self.father.menuEnable()
        self.deleteLater()

################################################################


class ShearingWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Shearing")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        options = list(map(str, range(1, self.father.masterData.data.ndim + 1)))
        grid.addWidget(wc.QLabel("Shearing constant:"), 0, 0)
        self.shearEntry = QtWidgets.QLineEdit()
        self.shearEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.shearEntry.setText("0.0")
        self.shearEntry.returnPressed.connect(self.shearPreview)
        grid.addWidget(self.shearEntry, 1, 0)
        grid.addWidget(wc.QLabel("Shearing direction:"), 2, 0)
        self.dirEntry = QtWidgets.QComboBox()
        self.dirEntry.addItems(options)
        self.dirEntry.setCurrentIndex(self.father.masterData.data.ndim - 2)
        grid.addWidget(self.dirEntry, 3, 0)
        grid.addWidget(wc.QLabel("Shearing axis:"), 4, 0)
        self.axEntry = QtWidgets.QComboBox()
        self.axEntry.addItems(options)
        self.axEntry.setCurrentIndex(self.father.masterData.data.ndim - 1)
        grid.addWidget(self.axEntry, 5, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def shearPreview(self, *args):
        shear = safeEval(self.shearEntry.text())
        if shear is not None:
            self.shearEntry.setText(str(float(shear)))

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        shear = safeEval(self.shearEntry.text())
        if shear is None:
            self.father.father.dispMsg("Not a valid value")
            return
        axes = self.dirEntry.currentIndex()
        axes2 = self.axEntry.currentIndex()
        if axes == axes2:
            self.father.father.dispMsg("Axes can't be the same for shearing")
            return
        else:
            self.father.redoList = []
            self.father.undoList.append(self.father.current.shearing(float(shear), axes, axes2))
            self.father.menuEnable()
            self.deleteLater()

##########################################################################################


class MultiplyWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Multiply")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Values:"), 0, 0)
        self.valEntry = QtWidgets.QLineEdit()
        self.valEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.valEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.valEntry, 1, 0)
        self.singleSlice = QtWidgets.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice, 1, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def preview(self, *args):
        env = vars(np).copy()
        env['length'] = int(self.father.current.data1D.shape[-1])  # so length can be used to in equations
        env['euro'] = lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal, num)
        val = eval(self.valEntry.text(), env)                # find a better solution, also add catch for exceptions
        self.father.current.multiplyPreview(np.array(val))

    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        env = vars(np).copy()
        env['length'] = int(self.father.current.data1D.shape[-1])  # so length can be used to in equations
        env['euro'] = lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal, num)
        val = eval(self.valEntry.text(), env)                # find a better solution, also add catch for exceptions
        returnValue = self.father.current.multiply(np.array(val), self.singleSlice.isChecked())
        if returnValue is None:
            return
        self.father.redoList = []
        self.father.undoList.append(returnValue)
        self.father.menuEnable()
        self.deleteLater()

##########################################################################################


class XaxWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("User defined x-axis")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Expression for x-axis values:"), 0, 0)
        self.valEntry = QtWidgets.QLineEdit()
        self.valEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.valEntry.returnPressed.connect(self.xaxPreview)
        grid.addWidget(self.valEntry, 1, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def xaxPreview(self, *args):
        env = vars(np).copy()
        env['length'] = int(self.father.current.data1D.shape[-1])  # so length can be used to in equations
        env['euro'] = lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal, num)
        val = eval(self.valEntry.text(), env)                # find a better solution, also add catch for exceptions
        if not isinstance(val, (list, np.ndarray)):
            self.father.father.dispMsg("Input is not a list or array")
            return
        if len(val) != self.father.current.data1D.shape[-1]:
            self.father.father.dispMsg("Length of input does not match length of data")
            return
        if not all(isinstance(x, (int, float)) for x in val):
            self.father.father.dispMsg("Array is not all of int or float type")
            return
        self.father.current.setXaxPreview(np.array(val))

    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        env = vars(np).copy()
        env['length'] = int(self.father.current.data1D.shape[-1])  # so length can be used to in equations
        env['euro'] = lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal, num)
        val = eval(self.valEntry.text(), env)                # find a better solution, also add catch for exceptions
        if not isinstance(val, (list, np.ndarray)):
            self.father.father.dispMsg("Input is not a list or array")
            return
        if len(val) != self.father.current.data1D.shape[-1]:
            self.father.father.dispMsg("Length of input does not match length of data")
            return
        if not all(isinstance(x, (int, float)) for x in val):
            self.father.father.dispMsg("Array is not all of int or float type")
            return
        self.father.redoList = []
        self.father.undoList.append(self.father.current.setXax(np.array(val)))
        self.father.menuEnable()
        self.deleteLater()

##########################################################################################


class RefWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Reference")

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
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        if parent.current.spec == 0:
            self.father.father.dispMsg('Setting ppm is only available for frequency data')
            self.deleteLater()
            return
        grid.addWidget(wc.QLabel("Name:"), 0, 0)
        self.refName = QtWidgets.QLineEdit()
        self.refName.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.refName, 1, 0)
        grid.addWidget(wc.QLabel("Frequency [MHz]:"), 2, 0)
        self.freqEntry = QtWidgets.QLineEdit()
        self.freqEntry.setText("%.7f" % (self.father.current.ref * 1e-6))
        self.freqEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.freqEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.freqEntry, 3, 0)
        grid.addWidget(wc.QLabel("Secondary Reference:"), 4, 0)
        self.refSecond = QtWidgets.QComboBox(parent=self)
        self.refSecond.addItems(self.secRefNames)
        self.refSecond.activated.connect(self.fillSecondaryRef)
        grid.addWidget(self.refSecond, 5, 0)

        grid.addWidget(wc.QLabel("Reference [ppm]:"), 6, 0)
        self.refEntry = QtWidgets.QLineEdit()
        self.refEntry.setText("0.0")
        self.refEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.refEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.refEntry, 7, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def preview(self, *args):
        freq = safeEval(self.freqEntry.text())
        ref = safeEval(self.refEntry.text())
        if freq is None or ref is None:
            return
        self.freqEntry.setText("%.7f" % (freq))
        self.refEntry.setText(str(ref))

    def fillSecondaryRef(self):
        self.refEntry.setText(self.secRefValues[self.refSecond.currentIndex()])

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):

        self.father.current.peakPickReset()
        freq = safeEval(self.freqEntry.text())
        ref = safeEval(self.refEntry.text())
        if freq is None or ref is None:
            self.father.father.dispMsg("Not a valid value")
            return
        freq = freq * 1e6
        reffreq = freq / (1.0 + ref * 1e-6)
        givenname = self.refName.text()
        nameOK = True
        if givenname:  # If name is filled in
            if self.father.father.referenceName.__contains__(givenname):  # if exists
                self.father.father.dispMsg("Reference name '" + givenname + "' already exists")
                nameOK = False
            else:
                self.father.mainProgram.referenceAdd(reffreq, givenname)

        if nameOK:
            self.father.redoList = []
            self.father.undoList.append(self.father.current.setRef(reffreq))
            self.father.menuEnable()
            self.deleteLater()

    def picked(self, pos):
        self.freqEntry.setText("%.7f" % ((self.father.current.ref + self.father.current.xax[pos[0]]) * 1e-6))
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

##########################################################################################


class HistoryWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Processing history")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        # grid.addWidget(wc.QLabel("History:"), 0, 0)
        self.valEntry = QtWidgets.QTextEdit()
        self.valEntry.setReadOnly(True)
        self.valEntry.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)
        self.valEntry.setText(self.father.masterData.getHistory())
        grid.addWidget(self.valEntry, 1, 0)
        cancelButton = QtWidgets.QPushButton("&Close")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        layout.setColumnStretch(1, 1)
        self.show()
        self.father.menuDisable()

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

#########################################################################################


class OrigListWidget(QtWidgets.QListWidget):

    def __init__(self, type, parent=None):
        super(OrigListWidget, self).__init__(parent)
        self.setDragDropMode(QtWidgets.QAbstractItemView.DragDrop)
        self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.setAcceptDrops(True)

    def dropEvent(self, event):
        pass

#########################################################################################


class DestListWidget(QtWidgets.QListWidget):

    def __init__(self, type, parent=None):
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
            for item in self.selectedItems():
                self.takeItem(self.row(item))

##########################################################################################


class CombineWorkspaceWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Combine workspaces")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 3)
        grid.addWidget(wc.QLabel("Workspaces:"), 0, 0)
        grid.addWidget(wc.QLabel("Combined spectrum:"), 0, 1)
        self.listA = OrigListWidget(self)
        for i in self.father.workspaceNames:
            QtWidgets.QListWidgetItem(i, self.listA).setToolTip(i)
        self.listB = DestListWidget(self)
        grid.addWidget(self.listA, 1, 0)
        grid.addWidget(self.listB, 1, 1)
        cancelButton = QtWidgets.QPushButton("&Close")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 2, 1)
        layout.setColumnStretch(2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def applyAndClose(self, *args):
        items = []
        for index in range(self.listB.count()):
            items.append(self.listB.item(index).text())
        if self.father.combineWorkspace(items):
            self.father.menuEnable()
            self.deleteLater()

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

##########################################################################################


class MonitorWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Monitor")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        fileName = self.father.masterData.filePath[1]
        if len(fileName) > 58:
            fileName = fileName[:55] + '...'
        fileLabel = wc.QLabel("File: " + fileName)
        fileLabel.setToolTip(self.father.masterData.filePath[1])
        layout.addWidget(fileLabel, 0, 0, 1, 3)
        layout.addLayout(grid, 1, 0, 1, 3)
        grid.addWidget(wc.QLabel("Macros:"), 0, 0)
        grid.addWidget(wc.QLabel("Apply after loading:"), 0, 1)
        self.listA = OrigListWidget(self)
        for i in self.father.father.macros.keys():
            QtWidgets.QListWidgetItem(i, self.listA).setToolTip(i)
        self.listB = DestListWidget(self)
        for i in self.father.monitorMacros:
            QtWidgets.QListWidgetItem(i, self.listB).setToolTip(i)
        grid.addWidget(self.listA, 1, 0)
        grid.addWidget(self.listB, 1, 1)
        cancelButton = QtWidgets.QPushButton("&Close")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 2, 0)
        watchButton = QtWidgets.QPushButton("&Watch")
        watchButton.clicked.connect(self.applyAndClose)
        layout.addWidget(watchButton, 2, 1)
        unwatchButton = QtWidgets.QPushButton("&Unwatch")
        unwatchButton.clicked.connect(self.stopAndClose)
        layout.addWidget(unwatchButton, 2, 2)
        layout.setColumnStretch(3, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def applyAndClose(self, *args):
        self.father.stopMonitor()
        items = []
        for index in range(self.listB.count()):
            items.append(self.listB.item(index).text())
        self.father.startMonitor(items)
        self.closeEvent()
        
    def stopAndClose(self, *args):
        self.father.stopMonitor()
        self.closeEvent()
    
    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()
        
##############################################################################


class PlotSettingsWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Preferences")
        tabWidget = QtWidgets.QTabWidget()
        tab1 = QtWidgets.QWidget()
        tab2 = QtWidgets.QWidget()
        # tab3 = QtWidgets.QWidget()
        tabWidget.addTab(tab1, "Plot")
        tabWidget.addTab(tab2, "Contour")
        grid1 = QtWidgets.QGridLayout()
        grid2 = QtWidgets.QGridLayout()
        # grid3 = QtWidgets.QGridLayout()
        tab1.setLayout(grid1)
        tab2.setLayout(grid2)
        # tab3.setLayout(grid3)
        grid1.setColumnStretch(10, 1)
        grid1.setRowStretch(10, 1)
        grid2.setColumnStretch(10, 1)
        grid2.setRowStretch(10, 1)
        # grid3.setColumnStretch(10, 1)
        # grid3.setRowStretch(10, 1)

        grid1.addWidget(QtWidgets.QLabel("Linewidth:"), 1, 0)
        self.lwSpinBox = QtWidgets.QDoubleSpinBox()
        self.lwSpinBox.setSingleStep(0.1)
        self.lwSpinBox.setValue(self.father.current.linewidth)
        self.lwSpinBox.valueChanged.connect(self.preview)
        grid1.addWidget(self.lwSpinBox, 1, 1)
        self.color = self.father.current.color
        lineColorButton = QtWidgets.QPushButton("Line color")
        lineColorButton.clicked.connect(self.setColor)
        grid1.addWidget(lineColorButton, 2, 0)
        self.xgridCheck = QtWidgets.QCheckBox("x-grid")
        self.xgridCheck.setChecked(self.father.current.grids[0])
        self.xgridCheck.stateChanged.connect(self.preview)
        grid1.addWidget(self.xgridCheck, 3, 0, 1, 2)
        self.ygridCheck = QtWidgets.QCheckBox("y-grid")
        self.ygridCheck.setChecked(self.father.current.grids[1])
        grid1.addWidget(self.ygridCheck, 4, 0, 1, 2)
        self.ygridCheck.stateChanged.connect(self.preview)

        grid2.addWidget(QtWidgets.QLabel("Colormap:"), 0, 0)
        self.cmEntry = QtWidgets.QComboBox(self)
        self.cmEntry.addItems(sc.COLORMAPLIST)
        self.cmEntry.setCurrentIndex(sc.COLORMAPLIST.index(self.father.current.colorMap))
        self.cmEntry.currentIndexChanged.connect(self.preview)
        grid2.addWidget(self.cmEntry, 0, 1)
        self.constColorCheck = QtWidgets.QCheckBox("Constant colors")
        self.constColorCheck.setChecked(self.father.current.contourConst)
        grid2.addWidget(self.constColorCheck, 1, 0)
        self.constColorCheck.stateChanged.connect(self.preview)
        self.posColor = self.father.current.contourColors[0]
        posColorButton = QtWidgets.QPushButton("Positive color")
        posColorButton.clicked.connect(self.setPosColor)
        grid2.addWidget(posColorButton, 2, 0)
        self.negColor = self.father.current.contourColors[1]
        negColorButton = QtWidgets.QPushButton("Negative color")
        negColorButton.clicked.connect(self.setNegColor)
        grid2.addWidget(negColorButton, 3, 0)

        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(tabWidget, 0, 0, 1, 4)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 1, 1)
        self.show()
        self.father.menuDisable()

    def preview(self, *args):
        tmpLw = self.father.current.linewidth
        self.father.current.setLw(self.lwSpinBox.value())
        tmpColor = self.father.current.color
        self.father.current.setColor(self.color)
        tmpColorMap = self.father.current.getColorMap()
        self.father.current.setColorMap(self.cmEntry.currentIndex())
        tmpGrids = self.father.current.grids
        self.father.current.setGrids([self.xgridCheck.isChecked(), self.ygridCheck.isChecked()])
        tmpContourConst = self.father.current.contourConst
        self.father.current.setContourConst(self.constColorCheck.isChecked())
        tmpContourColors = self.father.current.contourColors
        self.father.current.setContourColors([self.posColor, self.negColor])
        self.father.current.showFid()
        self.father.current.setLw(tmpLw)
        self.father.current.setColor(tmpColor)
        self.father.current.setColorMap(tmpColorMap)
        self.father.current.setGrids(tmpGrids)
        self.father.current.setContourConst(tmpContourConst)
        self.father.current.setContourColors(tmpContourColors)

    def setColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtWidgets.QColor(self.color))
        if tmp.isValid():
            self.color = tmp.name()
        self.preview()

    def setPosColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtWidgets.QColor(self.posColor))
        if tmp.isValid():
            self.posColor = tmp.name()
        self.preview()

    def setNegColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtWidgets.QColor(self.negColor))
        if tmp.isValid():
            self.negColor = tmp.name()
        self.preview()

    def applyAndClose(self, *args):
        self.father.current.setColor(self.color)
        self.father.current.setLw(self.lwSpinBox.value())
        self.father.current.setGrids([self.xgridCheck.isChecked(), self.ygridCheck.isChecked()])
        self.father.current.setColorMap(self.cmEntry.currentIndex())
        self.father.current.setContourConst(self.constColorCheck.isChecked())
        self.father.current.setContourColors([self.posColor, self.negColor])
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def closeEvent(self, *args):
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

##############################################################################


class PreferenceWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Preferences")
        tabWidget = QtWidgets.QTabWidget()
        tab1 = QtWidgets.QWidget()
        tab2 = QtWidgets.QWidget()
        tab3 = QtWidgets.QWidget()
        tabWidget.addTab(tab1, "Window")
        tabWidget.addTab(tab2, "Plot")
        tabWidget.addTab(tab3, "Contour")
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
        # grid1.addWidget(wc.QLabel("Window size:"), 0, 0, 1, 2)
        grid1.addWidget(wc.QLabel("Width:"), 1, 0)
        self.widthSpinBox = QtWidgets.QSpinBox()
        self.widthSpinBox.setMaximum(100000)
        self.widthSpinBox.setMinimum(1)
        self.widthSpinBox.setValue(self.father.defaultWidth)
        grid1.addWidget(self.widthSpinBox, 1, 1)
        grid1.addWidget(wc.QLabel("Height:"), 2, 0)
        self.heightSpinBox = QtWidgets.QSpinBox()
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

        grid2.addWidget(QtWidgets.QLabel("Linewidth:"), 1, 0)
        self.lwSpinBox = QtWidgets.QDoubleSpinBox()
        self.lwSpinBox.setSingleStep(0.1)
        self.lwSpinBox.setValue(self.father.defaultLinewidth)
        grid2.addWidget(self.lwSpinBox, 1, 1)
        self.color = self.father.defaultColor
        lineColorButton = QtWidgets.QPushButton("Line color")
        lineColorButton.clicked.connect(self.setColor)
        grid2.addWidget(lineColorButton, 2, 0)
        self.xgridCheck = QtWidgets.QCheckBox("x-grid")
        self.xgridCheck.setChecked(self.father.defaultGrids[0])
        grid2.addWidget(self.xgridCheck, 3, 0, 1, 2)
        self.ygridCheck = QtWidgets.QCheckBox("y-grid")
        self.ygridCheck.setChecked(self.father.defaultGrids[1])
        grid2.addWidget(self.ygridCheck, 4, 0, 1, 2)

        grid3.addWidget(QtWidgets.QLabel("Colormap:"), 0, 0)
        self.cmEntry = QtWidgets.QComboBox(self)
        self.cmEntry.addItems(sc.COLORMAPLIST)
        self.cmEntry.setCurrentIndex(sc.COLORMAPLIST.index(self.father.defaultColorMap))
        grid3.addWidget(self.cmEntry, 0, 1)
        self.constColorCheck = QtWidgets.QCheckBox("Constant colors")
        self.constColorCheck.setChecked(self.father.defaultContourConst)
        grid3.addWidget(self.constColorCheck, 1, 0)
        self.posColor = self.father.defaultPosColor
        posColorButton = QtWidgets.QPushButton("Positive color")
        posColorButton.clicked.connect(self.setPosColor)
        grid3.addWidget(posColorButton, 2, 0)
        self.negColor = self.father.defaultNegColor
        negColorButton = QtWidgets.QPushButton("Negative color")
        negColorButton.clicked.connect(self.setNegColor)
        grid3.addWidget(negColorButton, 3, 0)
        grid3.addWidget(QtWidgets.QLabel("Width ratio:"), 4, 0)
        self.WRSpinBox = QtWidgets.QDoubleSpinBox()
        self.WRSpinBox.setSingleStep(0.1)
        self.WRSpinBox.setValue(self.father.defaultWidthRatio)
        grid3.addWidget(self.WRSpinBox, 4, 1)
        grid3.addWidget(QtWidgets.QLabel("Height ratio:"), 5, 0)
        self.HRSpinBox = QtWidgets.QDoubleSpinBox()
        self.HRSpinBox.setSingleStep(0.1)
        self.HRSpinBox.setValue(self.father.defaultHeightRatio)
        grid3.addWidget(self.HRSpinBox, 5, 1)

        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(tabWidget, 0, 0, 1, 4)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtWidgets.QPushButton("&Store")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 1, 1)
        resetButton = QtWidgets.QPushButton("&Reset")
        resetButton.clicked.connect(self.reset)
        layout.addWidget(resetButton, 1, 2)
        layout.setColumnStretch(3, 1)
        self.show()

    def setColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtWidgets.QColor(self.color))
        if tmp.isValid():
            self.color = tmp.name()

    def setPosColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtWidgets.QColor(self.posColor))
        if tmp.isValid():
            self.posColor = tmp.name()

    def setNegColor(self, *args):
        tmp = QtWidgets.QColorDialog.getColor(QtWidgets.QColor(self.negColor))
        if tmp.isValid():
            self.negColor = tmp.name()

    def applyAndClose(self, *args):
        self.father.defaultWidth = self.widthSpinBox.value()
        self.father.defaultHeight = self.heightSpinBox.value()
        self.father.defaultMaximized = self.maximizedCheck.isChecked()
        self.father.defaultAskName = self.askNameCheck.isChecked()
        self.father.defaultLinewidth = self.lwSpinBox.value()
        self.father.defaultColor = self.color
        self.father.defaultGrids[0] = self.xgridCheck.isChecked()
        self.father.defaultGrids[1] = self.ygridCheck.isChecked()
        self.father.defaultColorMap = self.cmEntry.currentText()
        self.father.defaultContourConst = self.constColorCheck.isChecked()
        self.father.defaultPosColor = self.posColor
        self.father.defaultNegColor = self.negColor
        self.father.defaultWidthRatio = self.WRSpinBox.value()
        self.father.defaultHeightRatio = self.HRSpinBox.value()
        self.father.saveDefaults()
        self.closeEvent()

    def reset(self, *args):
        self.father.resetDefaults()
        self.father.saveDefaults()
        self.closeEvent()

    def closeEvent(self, *args):
        self.deleteLater()

##############################################################################


class shiftConversionWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent

        self.setWindowTitle("Chemical Shift Conversions")

        grid = QtWidgets.QGridLayout()
        grid.setColumnStretch(10, 1)
        grid.setRowStretch(14, 1)

        StConv = QtWidgets.QLabel("Standard Convention:")
        StConv.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(StConv, 0, 0)
        D11label = QtWidgets.QLabel(u'\u03b4' + '<sub>11</sub> (ppm)')
        D11label.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(D11label, 0, 1)
        D22label = QtWidgets.QLabel(u'\u03b4' + '<sub>22</sub> (ppm)')
        D22label.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(D22label, 0, 2)
        D33label = QtWidgets.QLabel(u'\u03b4' + '<sub>33</sub> (ppm)')
        D33label.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(D33label, 0, 3)
        standardGO = QtWidgets.QPushButton("Go")
        grid.addWidget(standardGO, 1, 0)
        standardGO.clicked.connect(lambda: self.shiftCalc(0))
        self.D11 = QtWidgets.QLineEdit()
        self.D11.setText("0")
        self.D11.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.D11, 1, 1)
        self.D22 = QtWidgets.QLineEdit()
        self.D22.setText("0")
        self.D22.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.D22, 1, 2)
        self.D33 = QtWidgets.QLineEdit()
        self.D33.setText("0")
        self.D33.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.D33, 1, 3)

        # xyz Convention
        grid.addWidget(QtWidgets.QLabel(""), 2, 0)
        StConv = QtWidgets.QLabel("xyz Convention:")
        StConv.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(StConv, 3, 0)
        dxxlabel = QtWidgets.QLabel(u'\u03b4' + '<sub>xx</sub> (ppm)')
        dxxlabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(dxxlabel, 3, 1)
        dyylabel = QtWidgets.QLabel(u'\u03b4' + '<sub>yy</sub> (ppm)')
        dyylabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(dyylabel, 3, 2)
        dzzlabel = QtWidgets.QLabel(u'\u03b4' + '<sub>zz</sub> (ppm)')
        dzzlabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(dzzlabel, 3, 3)

        xyzGO = QtWidgets.QPushButton("Go")
        grid.addWidget(xyzGO, 4, 0)
        xyzGO.clicked.connect(lambda: self.shiftCalc(1))
        self.dxx = QtWidgets.QLineEdit()
        self.dxx.setText("0")
        self.dxx.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.dxx, 4, 1)
        self.dyy = QtWidgets.QLineEdit()
        self.dyy.setText("0")
        self.dyy.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.dyy, 4, 2)
        self.dzz = QtWidgets.QLineEdit()
        self.dzz.setText("0")
        self.dzz.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.dzz, 4, 3)

        # Haeberlen Convention
        grid.addWidget(QtWidgets.QLabel(""), 5, 0)
        StConv = QtWidgets.QLabel("Haeberlen Convention:")
        StConv.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(StConv, 6, 0)
        disolabel = QtWidgets.QLabel(u'\u03b4' + '<sub>iso</sub> (ppm)')
        disolabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(disolabel, 6, 1)
        danisolabel = QtWidgets.QLabel(u'\u03b4' + '<sub>aniso</sub> (ppm)')
        danisolabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(danisolabel, 6, 2)
        etalabel = QtWidgets.QLabel(u'\u03b7')
        etalabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(etalabel, 6, 3)

        haeberGO = QtWidgets.QPushButton("Go")
        grid.addWidget(haeberGO, 7, 0)
        haeberGO.clicked.connect(lambda: self.shiftCalc(2))
        self.diso = QtWidgets.QLineEdit()
        self.diso.setText("0")
        self.diso.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.diso, 7, 1)
        self.daniso = QtWidgets.QLineEdit()
        self.daniso.setText("0")
        self.daniso.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.daniso, 7, 2)
        self.eta = QtWidgets.QLineEdit()
        self.eta.setText("0")
        self.eta.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.eta, 7, 3)

        # Hertzfeld berger
        grid.addWidget(QtWidgets.QLabel(""), 8, 0)
        HbConv = QtWidgets.QLabel("Hertzfeld-Berger Convention:")
        HbConv.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(HbConv, 9, 0)
        hbdisolabel = QtWidgets.QLabel(u'\u03b4' + '<sub>iso</sub> (ppm)')
        hbdisolabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(hbdisolabel, 9, 1)
        omegalabel = QtWidgets.QLabel(u'\u03a9 (ppm)')
        omegalabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(omegalabel, 9, 2)
        skewlabel = QtWidgets.QLabel(u'\u03ba')
        skewlabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(skewlabel, 9, 3)

        hbGO = QtWidgets.QPushButton("Go")
        grid.addWidget(hbGO, 10, 0)
        hbGO.clicked.connect(lambda: self.shiftCalc(3))
        self.hbdiso = QtWidgets.QLineEdit()
        self.hbdiso.setText("0")
        self.hbdiso.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.hbdiso, 10, 1)
        self.hbdaniso = QtWidgets.QLineEdit()
        self.hbdaniso.setText("0")
        self.hbdaniso.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.hbdaniso, 10, 2)
        self.hbskew = QtWidgets.QLineEdit()
        self.hbskew.setText("0")
        self.hbskew.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.hbskew, 10, 3)

        # Reset
        grid.addWidget(QtWidgets.QLabel(""), 11, 0)
        resetbutton = QtWidgets.QPushButton("Reset")
        grid.addWidget(resetbutton, 12, 0)
        resetbutton.clicked.connect(self.valueReset)

        closebutton = QtWidgets.QPushButton("Close")
        grid.addWidget(closebutton, 12, 3)
        closebutton.clicked.connect(self.closeEvent)

        self.setLayout(grid)
        self.show()

    def shiftCalc(self, Type):
        if Type == 0:  # If from standard
            try:
                delta11 = float(safeEval(self.D11.text()))
                delta22 = float(safeEval(self.D22.text()))
                delta33 = float(safeEval(self.D33.text()))
                Values = [delta11,delta22,delta33]
            except:
                self.father.dispMsg("Invalid input in Standard Convention")
                return
        if Type == 1:  # If from xyz
            try:
                delta11 = float(safeEval(self.dxx.text()))  # Treat xyz as 123, as it reorders them anyway
                delta22 = float(safeEval(self.dyy.text()))
                delta33 = float(safeEval(self.dzz.text()))
                Values = [delta11,delta22,delta33]
            except:
                self.father.dispMsg("Invalid input in xyz Convention")
                return
        if Type == 2:  # From haeberlen
            try:
                eta = float(safeEval(self.eta.text()))
                delta = float(safeEval(self.daniso.text()))
                iso = float(safeEval(self.diso.text()))
                Values = [iso,delta,eta]
            except:
                self.father.dispMsg("Invalid input in Haeberlen Convention")
                return                
        if Type == 3:  # From Hertzfeld-Berger
            try:
                iso = float(safeEval(self.hbdiso.text()))
                span = float(safeEval(self.hbdaniso.text()))
                skew = float(safeEval(self.hbskew.text()))
                Values = [iso,span,skew]
            except:
                self.father.dispMsg("Invalid input in Hertzfeld-Berger Convention")
                return    

        Results = fit.shiftConversion(Values,Type) #Do the actual conversion

        #Standard convention
        self.D11.setText('%#.4g' % Results[0][0])
        self.D22.setText('%#.4g' % Results[0][1])
        self.D33.setText('%#.4g' % Results[0][2])

        # Convert to haeberlen convention and xxyyzz
        self.dxx.setText('%#.4g' % Results[1][0])
        self.dyy.setText('%#.4g' % Results[1][1])
        self.dzz.setText('%#.4g' % Results[1][2])


        #Haeberlen def
        self.diso.setText('%#.4g' % Results[2][0])
        self.daniso.setText('%#.4g' % Results[2][1])
        try: #If a number
            self.eta.setText('%#.4g' % Results[2][2])
        except:
            self.eta.setText('ND')

        # Convert to Herzfeld-Berger Convention
        self.hbdiso.setText('%#.4g' % Results[3][0])
        self.hbdaniso.setText('%#.4g' % Results[3][1])
        try:
            self.hbskew.setText('%#.4g' % Results[3][2])
        except:
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

    def closeEvent(self):
        self.deleteLater()

class quadConversionWindow(QtWidgets.QWidget):
    
    Ioptions = ['1', '3/2', '2', '5/2', '3', '7/2', '4', '9/2','5','6','7']
    Ivalues = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,5.0,6.0,7.0]
   
    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Quadrupolar Coupling Conversions")
        grid = QtWidgets.QGridLayout()
        grid.setColumnStretch(10, 1)
        grid.setRowStretch(14, 1)
        Itext = QtWidgets.QLabel("I:")
        Itext.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(Itext, 0, 0)
        self.IEntry = QtWidgets.QComboBox()
        self.IEntry.addItems(self.Ioptions)
        self.IEntry.setCurrentIndex(0)
        grid.addWidget(self.IEntry, 1, 0)
        etalabel = QtWidgets.QLabel(u'\u03b7')
        etalabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(etalabel, 0, 1)
        self.Eta = QtWidgets.QLineEdit()
        self.Eta.setText("0")
        self.Eta.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.Eta, 1, 1)
        momentlabel = QtWidgets.QLabel('Q (fm<sup>2</sup>)')
        momentlabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(momentlabel, 0, 2)
        self.Moment = QtWidgets.QLineEdit()
        self.Moment.setText("ND")
        self.Moment.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.Moment, 1, 2)
        grid.addWidget(QtWidgets.QLabel(""),2, 0)
        CqConv = QtWidgets.QLabel("C<sub>Q</sub> Convention:")
        CqConv.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(CqConv, 3, 0)
        Cqlabel = QtWidgets.QLabel(u'C' + u'<sub>Q</sub>/2\u03c0 (MHz)')
        Cqlabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(Cqlabel, 3, 1)
        CqGO = QtWidgets.QPushButton("Go")
        grid.addWidget(CqGO, 4, 0)
        CqGO.clicked.connect(lambda: self.quadCalc(0))
        self.Cq = QtWidgets.QLineEdit()
        self.Cq.setText("0")
        self.Cq.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.Cq, 4, 1)
        grid.addWidget(QtWidgets.QLabel(""),5, 0)
        WqConv = QtWidgets.QLabel(u"\u03c9<sub>Q</sub> Convention:")
        WqConv.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(WqConv, 6, 0)
        Wqlabel = QtWidgets.QLabel(u'\u03c9' + u'<sub>Q</sub>/2\u03c0 (MHz)')
        Wqlabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(Wqlabel, 6, 1)
        WqGO = QtWidgets.QPushButton("Go")
        grid.addWidget(WqGO, 7, 0)
        WqGO.clicked.connect(lambda: self.quadCalc(1))
        self.Wq = QtWidgets.QLineEdit()
        self.Wq.setText("0")
        self.Wq.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.Wq, 7, 1)        
        grid.addWidget(QtWidgets.QLabel(""),8, 0)        
        VConv = QtWidgets.QLabel("Field gradients:")
        VConv.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(VConv, 9, 0)
        Vxxlabel = QtWidgets.QLabel('V<sub>xx</sub> (V/m<sup>2</sup>)')
        Vxxlabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(Vxxlabel, 9, 1)
        VGO = QtWidgets.QPushButton("Go")
        grid.addWidget(VGO, 10, 0)
        VGO.clicked.connect(lambda: self.quadCalc(2))
        self.Vxx = QtWidgets.QLineEdit()
        self.Vxx.setText("ND")
        self.Vxx.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.Vxx, 10, 1)        
        Vyylabel = QtWidgets.QLabel('V<sub>yy</sub> (V/m<sup>2</sup>)')
        Vyylabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(Vyylabel, 9, 2)
        self.Vyy = QtWidgets.QLineEdit()
        self.Vyy.setText("ND")
        self.Vyy.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.Vyy, 10, 2)        
        Vzzlabel = QtWidgets.QLabel('V<sub>zz</sub> (V/m<sup>2</sup>)')
        Vzzlabel.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(Vzzlabel, 9, 3)
        self.Vzz = QtWidgets.QLineEdit()
        self.Vzz.setText("ND")
        self.Vzz.setAlignment(QtCore.Qt.AlignHCenter)
        grid.addWidget(self.Vzz, 10, 3)

        # Reset
        grid.addWidget(QtWidgets.QLabel(""), 11, 0)
        resetbutton = QtWidgets.QPushButton("Reset")
        grid.addWidget(resetbutton, 12, 0)
        resetbutton.clicked.connect(self.valueReset)
        closebutton = QtWidgets.QPushButton("Close")
        grid.addWidget(closebutton, 12, 3)
        closebutton.clicked.connect(self.closeEvent)        
        self.setLayout(grid)
        self.show()

    def quadCalc(self, Type):
        I = self.Ivalues[self.IEntry.currentIndex()]
        if Type == 0: #Cq as input
            #Czz is equal to Cq, via same definition (scale) Cxx and Cyy can be found
            try:
                Czz = float(safeEval(self.Cq.text())) 
                
                Eta = float(safeEval(self.Eta.text())) 
                Cxx = Czz*(Eta-1)/2
                Cyy = -Cxx-Czz
            except:
                self.father.dispMsg("Invalid input in Cq definition")
                return 
        if Type == 1:
            try:
                Vmax = float(safeEval(self.Wq.text()))                 
                Eta = float(safeEval(self.Eta.text())) 
                Czz = Vmax*(2.0*I*(2*I-1))/3.0
                Cxx = Czz*(Eta-1)/2
                Cyy = -Cxx-Czz
            except:
                self.father.dispMsg("Invalid input in Wq definition")
                return                 
        if Type ==2:
             try:
                Vxx = float(safeEval(self.Vxx.text())) 
                Vyy = float(safeEval(self.Vyy.text())) 
                Vzz = float(safeEval(self.Vzz.text())) 
                Q = float(safeEval(self.Moment.text()))*1e-30 #get moment and convert from fm^2                
                #Force traceless
                if not np.isclose(Vxx+Vyy+Vzz,0.0):
                    Diff = (Vxx+Vyy+Vzz)/3.0
                    Vxx = Vxx - Diff
                    Vyy = Vyy - Diff
                    Vzz = Vzz - Diff                    
                Scaling = SC.elementary_charge*Q/SC.Planck 
                Czz = Vzz * Scaling/1e6 #scale for Cq definition in MHz
                Cxx = Vxx * Scaling/1e6
                Cyy = Vyy * Scaling/1e6
             except:
                self.father.dispMsg("Invalid input in field gradients")
                return             
        #sort    
        CArray = np.array([Cxx, Cyy, Czz])
        Cindex = np.argsort(np.abs(CArray))
        Csort = CArray[Cindex]
        if Csort[2]<0: #If Czz negative due to weird input, make it positive
            Csort=-Csort
        CqNew = Csort[2]
        if CqNew == 0.0:
            self.Eta.setText('ND')
        else:
            EtaNew = np.abs((Csort[0]-Csort[1])/Csort[2]) #Abs to avoid -0.0 rounding error
            self.Eta.setText('%#.4g' % EtaNew)
        WqNew = CqNew*3.0/(2.0*I*(2*I-1))
        self.Cq.setText('%#.4g' % CqNew)
        self.Wq.setText('%#.4g' % WqNew)
        try:
            Q = float(safeEval(self.Moment.text()))*1e-30 #get moment and convert from fm^2
            Scaling = SC.elementary_charge*Q/SC.Planck 
            Vxx = Csort[0]/Scaling * 1e6
            Vyy = Csort[1]/Scaling * 1e6
            Vzz = Csort[2]/Scaling * 1e6
            self.Vxx.setText('%#.4g' % Vxx)
            self.Vyy.setText('%#.4g' % Vyy)
            self.Vzz.setText('%#.4g' % Vzz)
        except:
            self.Moment.setText('ND')
            self.Vxx.setText('ND')
            self.Vyy.setText('ND')
            self.Vzz.setText('ND')

    def valueReset(self):  # Resets all the boxes to 0
        self.Cq.setText('0')
        self.Eta.setText('0')
        self.Wq.setText('0')
        self.Moment.setText('ND')
        self.Vxx.setText('ND')
        self.Vyy.setText('ND')
        self.Vzz.setText('ND')

    def closeEvent(self):
        self.deleteLater()


if __name__ == '__main__':
    mainProgram = MainProgram(root)
    mainProgram.setWindowTitle("ssNake - " + VERSION)
    mainProgram.show()
    splash.finish(mainProgram)
    sys.exit(root.exec_())
    
