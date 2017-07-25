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
import traceback as tb
import datetime
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

    
splashSteps=15.0/100
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
import functions as func
splashStep = splashProgressStep(splashStep)
import scipy.constants as SC
splashStep = splashProgressStep(splashStep)
import loadFiles as LF
splashStep = splashProgressStep(splashStep)

matplotlib.rc('font', family='DejaVu Sans')
np.set_printoptions(threshold=np.nan)
QtCore.QLocale.setDefault(QtCore.QLocale('en_US'))

pi = np.pi

VERSION = 'v0.7b'


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
        self.LastLocation = ''
        self.initMenu()
        self.menuCheck()
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
        self.initToolbar()
        self.resize(self.defaultWidth, self.defaultHeight)
        if self.defaultMaximized:
            self.showMaximized()
        QtWidgets.QShortcut(QtGui.QKeySequence.Paste, self).activated.connect(self.handlePaste)
        QtWidgets.QShortcut(QtGui.QKeySequence.Copy, self).activated.connect(self.handleCopy)

    def dispError(self,error):
        CurTime = datetime.datetime.now()
        TimeStr = '{0:02d}'.format(CurTime.hour) + ':' + '{0:02d}'.format(CurTime.minute) + ':' + '{0:02d}'.format(CurTime.second)
        self.errors.append([TimeStr,error])
        
    def handlePaste(self):
        self.dropEvent(QtWidgets.QApplication.instance().clipboard())

    def handleCopy(self):
        if self.mainWindow is None:
            return
        pixmap = QtGui.QPixmap.grabWidget(self.mainWindow.canvas)
        QtWidgets.QApplication.clipboard().setPixmap(pixmap)
        
    def resetDefaults(self):
        self.defaultUnits = 1
        self.defaultPPM = False
        self.defaultWidth = 1
        self.defaultHeight = 1
        self.defaultMaximized = False
        self.defaultAskName = True
        self.defaultToolBar = True
        self.defaultLinewidth = 1.0
        self.defaultColor = '#0000FF'
        self.defaultGrids = [False, False]
        self.defaultDiagonalBool = False
        self.defaultDiagonalMult = 1
        self.defaultZeroScroll = True
        self.defaultColorMap = 'seismic'
        self.defaultWidthRatio = 3.0
        self.defaultHeightRatio = 3.0
        self.defaultContourConst = True
        self.defaultPosColor = '#FF0000'
        self.defaultNegColor = '#0000FF'
        self.defaultToolbarActionList = ['File --> Open','File -- > Save --> Matlab','File --> Export --> Figure','Seperator',
                                     'Workspaces --> Duplicate','Workspaces --> Delete','Seperator','Edit --> Undo','Edit --> Redo',
                                     'Edit --> Reload','Seperator','Tools --> Apodize','Tools --> Phase','Tools --> Autophase 0','Seperator',
                                     'Matrix --> Sizing','Matrix --> Shift Data','Matrix --> Multiply','Seperator','Fitting --> S/N','Fitting --> FWHM',
                                     'Fitting --> Integrals','Fitting --> Relaxation Curve','Fitting --> Lorentzian/Gaussian','Seperator',
                                     'Plot --> 1D Plot','Plot --> Stack Plot','Plot --> Array Plot','Plot --> Contour Plot',
                                     'Plot --> Multi Plot','Seperator','History --> History','History --> Clear Undo/Redo List',
                                     'Seperator','Utilities --> NMR Table'] 

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

        self.defaultColor = settings.value("plot/colour", self.defaultColor, str)
        try:
            self.defaultLinewidth = settings.value("plot/linewidth", self.defaultLinewidth, float)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the plot/linewidth")
        self.defaultGrids = [settings.value("plot/xgrid", self.defaultGrids[0], bool), settings.value("plot/ygrid", self.defaultGrids[1], bool)]
        self.defaultZeroScroll = settings.value("plot/zeroscroll", self.defaultZeroScroll, bool)
        self.defaultColorMap = settings.value("contour/colourmap", self.defaultColorMap, str)
        self.defaultContourConst = settings.value("contour/constantcolours", self.defaultContourConst, bool)
        self.defaultPosColor = settings.value("contour/poscolour", self.defaultPosColor, str)
        self.defaultNegColor = settings.value("contour/negcolour", self.defaultNegColor, str)
        if not str(self.defaultColorMap) in sc.COLORMAPLIST:
            self.dispMsg("Incorrect colourmap in config file")
        self.defaultDiagonalBool = settings.value("contour/diagonalbool", self.defaultDiagonalBool, bool)
        try:
            self.defaultDiagonalMult = settings.value("contour/diagonalmult", self.defaultDiagonalMult, float)
        except TypeError:
            self.dispMsg("Incorrect value in the config file for the diagonal multiplier")
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
        settings.setValue("plot/units", self.defaultUnits)
        settings.setValue("plot/ppm", self.defaultPPM)
        settings.setValue('toolbarList', self.defaultToolbarActionList)
        settings.setValue("plot/colour", self.defaultColor)
        settings.setValue("plot/linewidth", self.defaultLinewidth)
        settings.setValue("plot/xgrid", self.defaultGrids[0])
        settings.setValue("plot/ygrid", self.defaultGrids[1])
        settings.setValue("plot/zeroscroll", self.defaultZeroScroll)
        settings.setValue("maximized", self.defaultMaximized)
        settings.setValue("width", self.defaultWidth)
        settings.setValue("height", self.defaultHeight)
        settings.setValue("ask_name", self.defaultAskName)
        settings.setValue("toolbar", self.defaultToolBar)
        settings.setValue("contour/colourmap", self.defaultColorMap)
        settings.setValue("contour/constantcolours", self.defaultContourConst)
        settings.setValue("contour/poscolour", self.defaultPosColor)
        settings.setValue("contour/negcolour", self.defaultNegColor)
        settings.setValue("contour/width_ratio", self.defaultWidthRatio)
        settings.setValue("contour/height_ratio", self.defaultHeightRatio)
        settings.setValue("contour/diagonalbool", self.defaultDiagonalBool)
        settings.setValue("contour/diagonalmult", self.defaultDiagonalMult)

    def dispMsg(self, msg, color = 'black', error = True):
        if color == 'red':
            self.statusBar.setStyleSheet("QStatusBar{padding-left:8px;color:red;}")
        else:
            self.statusBar.setStyleSheet("QStatusBar{padding-left:8px;color:black;}")
        if error:
            self.dispError([msg])
        self.statusBar.showMessage(msg, 10000)
    
    def initToolbar(self):
        if self.defaultToolBar:
            self.toolbar = self.addToolBar('Toolbar')
            self.toolbar.setMovable(False)
            self.toolbar.setIconSize(QtCore.QSize(22,22))
            
            self.seperatorAction = []
            
            self.allActionsList = [['Seperator',None],['File --> Open',self.openAct],['File --> Save --> JSON',self.saveAct],['File -- > Save --> Matlab',self.saveMatAct],
                                   ['File --> Export --> Figure',self.savefigAct],['File --> Export --> Simpson',self.saveSimpsonAct],['File --> Export --> ASCII (1D/2D)',self.saveASCIIAct],
                                    ['File --> Preferences',self.preferencesAct],['File --> Quit',self.quitAct],
                                    ['Workspaces --> Duplicate',self.newAct],['Workspaces --> Delete',self.closeAct],['Workspaces --> Rename',self.renameWorkspaceAct],
                                    ['Workspaces --> Next',self.forwardAct],['Workspaces --> Previous',self.backAct],
                                    ['Macro --> Start Recording',self.macrostartAct],['Macro --> Stop Recording',self.macrostopAct],['Macro --> Load',self.macroLoadAct],
                                    ['Edit --> Undo',self.undoAction],['Edit --> Redo',self.redoAction],['Edit --> Reload',self.reloadAct],['Edit --> Monitor',self.monitorAct],
                                    ['Tools --> Real',self.realAct],['Tools --> Imag',self.imagAct],['Tools --> Abs',self.absAct],['Tools --> Complex Conjugate',self.conjAct],['Tools --> Apodize',self.apodizeAct],
                                    ['Tools --> Phase',self.phaseAct],['Tools --> Autophase 0',self.autoPhaseAct0],['Tools --> Autophase 0+1',self.autoPhaseAct1],['Tools --> Swap Echo',self.swapEchoAct],['Tools --> Offset Correction',self.corOffsetAct],
                                    ['Tools --> Baseline Correction',self.baselineAct],['Tools --> Subtract Averages',self.subAvgAct],['Tools --> Reference Deconvolution',self.refDeconvAct],
                                    ['Tools --> Correct Bruker Digital Filter',self.brukDigitalAct],['Tools --> Hypercomplex --> States',self.statesAct],['Tools --> Hypercomplex --> TPPI',self.statesTPPIAct],['Tools --> Hypercomplex --> Echo-antiecho',self.echoantiAct],
                                    ['Tools --> LPSVD',self.lpsvdAct],
                                    ['Matrix --> Sizing',self.sizingAct],['Matrix --> Shift Data',self.shiftAct],['Matrix --> Multiply',self.multiplyAct],['Matrix --> Region --> Integrate',self.intRegionAct],
                                    ['Matrix --> Region --> Sum',self.sumRegionAct],['Matrix --> Region --> Max',self.maxRegionAct],['Matrix --> Region --> Min',self.minRegionAct],
                                    ['Matrix --> Region --> Max Position',self.maxposRegionAct],['Matrix --> Region --> Min Position',self.minposRegionAct],['Matrix --> Region --> Average',self.averageRegionAct],
                                    ['Matrix --> Diff',self.diffAct],['Matrix --> Cumsum',self.cumsumAct],['Matrix --> Extract Part',self.extractpartAct],['Matrix --> Flip L/R',self.fliplrAct],
                                    ['Matrix --> Delete',self.matrixdelAct],['Matrix --> Split',self.splitAct],['Matrix --> Multiply',self.multiplyAct],['Matrix --> Reorder',self.reorderAct], ['Matrix --> Regrid',self.regridAct],
                                    ['Matrix --> Concatenate',self.concatAct],['Matrix --> Shearing',self.shearAct],
                                    ['Transforms --> Fourier Transform',self.fourierAct],['Transforms --> Real Fourier Transform',self.realFourierAct],['Transforms --> Fftshift',self.fftshiftAct],
                                    ['Transforms --> Inv fftshift',self.invfftshiftAct],['Transforms --> Hilbert Transform',self.hilbertAct],['Transforms --> NUS --> FFM',self.ffmAct],
                                    ['Transforms --> NUS --> CLEAN',self.cleanAct],['Transforms --> NUS --> IST',self.istAct],
                                    ['Fitting --> S/N',self.snrAct],['Fitting --> FWHM',self.fwhmAct],['Fitting --> Centre of Mass',self.massAct],
                                    ['Fitting --> Integrals',self.intfitAct],['Fitting --> Relaxation Curve',self.relaxAct],['Fitting --> Diffusion Curve',self.diffusionAct],
                                    ['Fitting --> Lorentzian/Gaussian',self.lorentzfitAct],['Fitting --> CSA',self.csastaticAct],
                                    ['Fitting --> First Order Quadrupole',self.firstquadstatAct],['Fitting --> Second Order Quadrupole',self.secondquadstatAct],
                                    ['Fitting --> Czjzek',self.czjzekstatAct],
                                    ['Combine --> Combine Workspaces',self.combineWorkspaceAct],['Combine --> Insert From Workspace',self.insertdatAct],['Combine --> Add',self.adddatAct],['Combine --> Subtract',self.subdatAct],['Combine --> Multiply',self.multdatAct],
                                    ['Combine --> Divide',self.divdatAct],
                                    ['Plot --> 1D Plot',self.onedplotAct],['Plot --> Scatter',self.scatterplotAct],['Plot --> Stack Plot',self.stackplotAct],
                                    ['Plot --> Array Plot',self.arrayplotAct],['Plot --> Contour Plot',self.contourplotAct],['Plot --> Multi Plot',self.multiplotAct],
                                    ['Plot --> Set Reference',self.setrefAct],['Plot --> Clear Current Reference',self.delrefAct],['Plot --> Load Reference',self.loadrefAct],['Plot --> User X-axis',self.userxAct],
                                    ['Plot --> Plot Settings',self.plotprefAct],
                                    ['History --> History',self.historyAct],['History --> Clear Undo/Redo List',self.clearundoAct],
                                    ['Utilities --> Chemical Shift Conversion Tool',self.shiftconvAct],['Utilities --> Quadrupole Coupling Conversion Tool',self.quadconvAct],['Utilities --> NMR Table',self.nmrtableAct],
                                    ['Help --> Update',self.updateAct],['Help --> About',self.aboutAct]]
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
        self.combineLoadAct = self.filemenu.addAction(QtGui.QIcon(IconDirectory + 'combine.png'),'&Open && Combine', self.createCombineLoadWindow)
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
        self.saveASCIIAct = self.exportmenu.addAction(QtGui.QIcon(IconDirectory + 'ssnake.png'), 'ASCII (1D/2D)', self.saveASCIIFile)
        self.saveASCIIAct.setToolTip('Save as ASCII Text File')
        self.preferencesAct = self.filemenu.addAction(QtGui.QIcon(IconDirectory + 'preferences.png'), '&Preferences', lambda: PreferenceWindow(self))
        self.preferencesAct.setToolTip('Open Preferences Window')
        self.quitAct = self.filemenu.addAction(QtGui.QIcon(IconDirectory + 'quit.png'), '&Quit', self.fileQuit, QtGui.QKeySequence.Quit)
        self.quitAct.setToolTip('Close ssNake')
        
        self.saveActList = [self.saveAct,self.saveMatAct]
        self.exportActList = [self.savefigAct,self.saveSimpsonAct,self.saveASCIIAct]        
        self.fileActList = [self.openAct,self.saveAct,self.saveMatAct,self.savefigAct,
                            self.saveSimpsonAct,self.saveASCIIAct,self.combineLoadAct,self.preferencesAct,self.quitAct]
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

        self.workspaceActList = [self.newAct,self.closeAct,self.renameWorkspaceAct,self.forwardAct,
                                 self.backAct]
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
        
        self.macroActList = [self.macrostartAct,self.macrostopAct]

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
        self.noUndoAct = QtWidgets.QAction("&No Undo Mode", self.editmenu,checkable = True)
        self.noUndoAct.toggled.connect(self.noUndoMode)
        self.editmenu.addAction(self.noUndoAct)
        self.clearundoAct = self.editmenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'),"&Clear Undo/Redo List", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.clearUndo()))
        self.clearundoAct.setToolTip('Clear Undo/Redo List')
        self.reloadAct = self.editmenu.addAction(QtGui.QIcon(IconDirectory + 'reload.png'), "Re&load", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.reloadLast()), QtGui.QKeySequence.Refresh)
        self.reloadAct.setToolTip('Reload Current Data')
        self.monitorAct = self.editmenu.addAction(QtGui.QIcon(IconDirectory + 'monitor.png'),"&Monitor", lambda: self.mainWindowCheck(lambda mainWindow: MonitorWindow(mainWindow)))
        self.monitorAct.setToolTip('Monitor Current Data')
        self.editActList = [self.undoAction,self.redoAction,self.clearundoAct, self.noUndoAct ,self.reloadAct,self.monitorAct]        
        
        # the tool drop down menu
        self.toolMenu = QtWidgets.QMenu("&Tools", self)
        self.menubar.addMenu(self.toolMenu)
        self.realAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'real.png'), "&Real", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.real()))
        self.realAct.setToolTip('Take Real Part of Data')
        self.imagAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'imag.png'), "&Imag", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.imag()))
        self.imagAct.setToolTip('Take Imaginary Part of Data')
        self.absAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'abs.png'), "&Abs", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.abs()))
        self.absAct.setToolTip('Take Absolute of Data')
        self.conjAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'complexconj.png'),"&Complex Conjugate", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.conj()))
        self.conjAct.setToolTip('Take Complex Conjugate of Data')
        self.apodizeAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'apodize.png'),"Apo&dize", lambda: self.mainWindowCheck(lambda mainWindow: ApodWindow(mainWindow)))
        self.apodizeAct.setToolTip('Open Apodize Window')
        self.phaseAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'phase.png'), "&Phasing", lambda: self.mainWindowCheck(lambda mainWindow: PhaseWindow(mainWindow)))
        self.phaseAct.setToolTip('Open Phasing Window')
        self.autoPhaseAct0 = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'autophase0.png'),"Autophase 0", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.directAutoPhase(0)))
        self.autoPhaseAct0.setToolTip('Autophase 0 order')
        self.autoPhaseAct1 = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'autophase1.png'),"Autophase 0+1", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.directAutoPhase(1)))
        self.autoPhaseAct1.setToolTip('Autophase 0 and 1 order')
        self.swapEchoAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'swapecho.png'), "Swap &Echo", lambda: self.mainWindowCheck(lambda mainWindow: SwapEchoWindow(mainWindow)))
        self.swapEchoAct.setToolTip('Swap Echo')
        self.corOffsetAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'offset.png'),"&Offset Correction", lambda: self.mainWindowCheck(lambda mainWindow: DCWindow(mainWindow)))
        self.corOffsetAct.setToolTip('Offset Correction')
        self.baselineAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'baseline.png'),"&Baseline Correction", lambda: self.mainWindowCheck(lambda mainWindow: BaselineWindow(mainWindow)))
        self.baselineAct.setToolTip('Baseline Correction')
        self.subAvgAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'subaverage.png'),"S&ubtract Averages", lambda: self.mainWindowCheck(lambda mainWindow: SubtractAvgWindow(mainWindow)))
        self.subAvgAct.setToolTip('Subtract Averages')
        self.refDeconvAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'deconvolute.png'),"Re&ference Deconvolution", lambda: self.mainWindowCheck(lambda mainWindow: FiddleWindow(mainWindow)))
        self.refDeconvAct.setToolTip('Reference Deconvolution')
        self.brukDigitalAct = self.toolMenu.addAction(QtGui.QIcon(IconDirectory + 'bruker.png'),"&Correct Bruker Digital Filter", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.BrukerDigital()))
        self.brukDigitalAct.setToolTip("Correct Bruker Digital Filter")
        self.lpsvdAct = self.toolMenu.addAction("&LPSVD", lambda: self.mainWindowCheck(lambda mainWindow: LPSVDWindow(mainWindow)))
        self.lpsvdAct.setToolTip('LPSVD linear prediction')
        
        self.hypercomplexMenu = QtWidgets.QMenu("Hypercomplex", self)
        self.toolMenu.addMenu(self.hypercomplexMenu)
        self.statesAct = self.hypercomplexMenu.addAction(QtGui.QIcon(IconDirectory + 'States.png'),"&States", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.states()))
        self.statesAct.setToolTip('States Hypercomplex Data Processing')
        self.statesTPPIAct = self.hypercomplexMenu.addAction(QtGui.QIcon(IconDirectory + 'statestppi.png'),"States-&TPPI", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.statesTPPI()))
        self.statesTPPIAct.setToolTip('States-TPPI Hypercomplex Data Processing')
        self.echoantiAct = self.hypercomplexMenu.addAction(QtGui.QIcon(IconDirectory + 'echoantiecho.png'),"Ec&ho-antiecho", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.echoAntiEcho()))
        self.echoantiAct.setToolTip('Ec&ho-antiecho Hypercomplex Data Processing')
        

        self.toolsActList = [self.realAct,self.imagAct,self.absAct,self.apodizeAct,self.phaseAct,self.autoPhaseAct0,self.autoPhaseAct1,
                             self.swapEchoAct,self.corOffsetAct,self.baselineAct,self.subAvgAct,self.refDeconvAct,self.statesAct,
                             self.statesTPPIAct,self.echoantiAct,self.brukDigitalAct,self.lpsvdAct]
        
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
        self.regridAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'regrid.png'),"Regrid", lambda: self.mainWindowCheck(lambda mainWindow: RegridWindow(mainWindow)))
        self.regridAct.setToolTip('Regrid')
        self.concatAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'concatenate.png'),"C&oncatenate", lambda: self.mainWindowCheck(lambda mainWindow: ConcatenateWindow(mainWindow)))
        self.concatAct.setToolTip('Concatenate')
        self.multiDActions.append(self.concatAct)
        self.shearAct = self.matrixMenu.addAction(QtGui.QIcon(IconDirectory + 'shear.png'),"Shearin&g", lambda: self.mainWindowCheck(lambda mainWindow: ShearingWindow(mainWindow)))
        self.shearAct.setToolTip('Shearing')
        self.multiDActions.append(self.shearAct)
        
        self.matrixActList = [self.sizingAct,self.shiftAct,self.intRegionAct,self.sumRegionAct,self.maxRegionAct,
                              self.minRegionAct,self.maxposRegionAct,self.minposRegionAct,self.averageRegionAct,
                              self.diffAct,self.cumsumAct,self.extractpartAct,self.fliplrAct,self.matrixdelAct,
                              self.splitAct,self.multiplyAct,self.reorderAct,self.regridAct,self.concatAct,self.shearAct]

        # the fft drop down menu
        self.fftMenu = QtWidgets.QMenu("T&ransforms", self)
        self.menubar.addMenu(self.fftMenu)
        self.fourierAct = self.fftMenu.addAction(QtGui.QIcon(IconDirectory + 'fourier.png'),"&Fourier Transform", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.fourier()), QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        self.fourierAct.setToolTip('Fourier Transform')
        self.realFourierAct = self.fftMenu.addAction(QtGui.QIcon(IconDirectory + 'realfourier.png'),"&Real Fourier Transform", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.realFourier()))
        self.realFourierAct.setToolTip('Real Fourier Transform')
        self.fftshiftAct = self.fftMenu.addAction(QtGui.QIcon(IconDirectory + 'fftshift.png'),"Fft&shift", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.fftshift()))
        self.fftshiftAct.setToolTip('Fftshift')
        self.invfftshiftAct = self.fftMenu.addAction(QtGui.QIcon(IconDirectory + 'ifftshift.png'),"&Inv fftshift", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.invFftshift()))
        self.invfftshiftAct.setToolTip('Inverse fftshift')
        self.hilbertAct = self.fftMenu.addAction(QtGui.QIcon(IconDirectory + 'hilbert.png'),"&Hilbert Transform", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.hilbert()))
        self.hilbertAct.setToolTip('Hilbert Transform') 
        self.nusMenu = QtWidgets.QMenu("&NUS", self)
        self.fftMenu.addMenu(self.nusMenu)
        self.ffmAct = self.nusMenu.addAction(QtGui.QIcon(IconDirectory + 'ffm.png'),"&FFM", lambda: self.mainWindowCheck(lambda mainWindow: FFMWindow(mainWindow)))
        self.ffmAct.setToolTip('FFM') 
        self.cleanAct = self.nusMenu.addAction(QtGui.QIcon(IconDirectory + 'clean.png'),"&CLEAN", lambda: self.mainWindowCheck(lambda mainWindow: CLEANWindow(mainWindow)))
        self.cleanAct.setToolTip('CLEAN') 
        self.istAct = self.nusMenu.addAction(QtGui.QIcon(IconDirectory + 'ist.png'),"&IST", lambda: self.mainWindowCheck(lambda mainWindow: ISTWindow(mainWindow)))
        self.istAct.setToolTip('IST') 
        
        self.fftActList = [self.fourierAct,self.realFourierAct,self.fftshiftAct,self.invfftshiftAct,
                           self.hilbertAct,self.ffmAct,self.cleanAct,self.istAct]

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
        self.csastaticAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'csastatic.png'),"&CSA", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createTensorDeconvWindow()))
        self.csastaticAct.setToolTip('Fit CSA')
        self.firstquadstatAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'firstquadstatic.png'),"First Order &Quadrupole", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuad1DeconvWindow()))
        self.firstquadstatAct.setToolTip('Fit First Order Quadrupole')
        self.secondquadstatAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'secondquadstatic.png'),"S&econd Order Quadrupole", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuad2DeconvWindow()))
        self.secondquadstatAct.setToolTip('Fit Second Order Quadrupole')
        self.czjzekstatAct = self.fittingMenu.addAction(QtGui.QIcon(IconDirectory + 'czjzekstatic.png'),"C&zjzek", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.createQuad2CzjzekWindow()))
        self.czjzekstatAct.setToolTip('Fit Czjzek Pattern')
        
        self.fittingActList = [self.snrAct,self.fwhmAct,self.massAct,self.intfitAct,self.relaxAct,
                               self.diffusionAct,self.lorentzfitAct,self.csastaticAct,
                               self.firstquadstatAct,self.secondquadstatAct,
                               self.czjzekstatAct]
        # the combine drop down menu
        self.combineMenu = QtWidgets.QMenu("Com&bine", self)
        self.menubar.addMenu(self.combineMenu)
        self.combineWorkspaceAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'combine.png'), '&Combine Workspaces', self.createCombineWorkspaceWindow)
        self.combineWorkspaceAct.setToolTip('Combine Workspaces')
        self.insertdatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'insert.png'),"&Insert From Workspace", lambda: self.mainWindowCheck(lambda mainWindow: InsertWindow(mainWindow)))
        self.insertdatAct.setToolTip('Insert From Workspace')
        self.adddatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'add.png'), "&Add", lambda: self.mainWindowCheck(lambda mainWindow: CombineWindow(mainWindow, 0)))
        self.adddatAct.setToolTip('Add Data From Workspace')
        self.subdatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'subtract.png'), "&Subtract", lambda: self.mainWindowCheck(lambda mainWindow: CombineWindow(mainWindow, 1)))
        self.subdatAct.setToolTip('Subtract Data From Workspace')
        self.multdatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'multiplyWorkspace.png'),"&Multiply", lambda: self.mainWindowCheck(lambda mainWindow: CombineWindow(mainWindow, 2)))
        self.multdatAct.setToolTip('Multiply Data From Workspace')
        self.divdatAct = self.combineMenu.addAction(QtGui.QIcon(IconDirectory + 'divideWorkspace.png'),"&Divide", lambda: self.mainWindowCheck(lambda mainWindow: CombineWindow(mainWindow, 3)))
        self.divdatAct.setToolTip('Divide Data From Workspace')
        
        self.combineActList = [self.combineWorkspaceAct,self.insertdatAct,self.adddatAct,self.subdatAct,self.multdatAct,self.divdatAct]
        
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
#        self.skewplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'skewed.png'),"S&kewed Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotSkewed()))
#        self.skewplotAct.setToolTip('Skew Plot')
#        self.multiDActions.append(self.skewplotAct)
        self.multiplotAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'multi.png'),"&Multi Plot", lambda: self.mainWindowCheck(lambda mainWindow: mainWindow.plotMulti()))
        self.multiplotAct.setToolTip('Multi Plot')

        self.referencelistmenu = QtWidgets.QMenu('&Reference', self)
        self.plotMenu.addMenu(self.referencelistmenu)
        self.setrefAct = self.referencelistmenu.addAction(QtGui.QIcon(IconDirectory + 'setreference.png'),"&Set Reference", lambda: self.mainWindowCheck(lambda mainWindow: RefWindow(mainWindow)))
        self.setrefAct.setToolTip('Set Reference')
        self.delrefAct = self.referencelistmenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'),"&Clear Current Reference", self.referenceClear)
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

        self.userxAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'xaxis.png'),"&User X-axis", lambda: self.mainWindowCheck(lambda mainWindow: XaxWindow(mainWindow)))
        self.userxAct.setToolTip('User X-axis')
        self.plotprefAct = self.plotMenu.addAction(QtGui.QIcon(IconDirectory + 'preferences.png'),"&Plot Settings", lambda: self.mainWindowCheck(lambda mainWindow: PlotSettingsWindow(mainWindow)))
        self.plotprefAct.setToolTip('Plot Settings')
        
        self.plotActList = [self.onedplotAct,self.scatterplotAct,self.stackplotAct,self.arrayplotAct,
                            self.contourplotAct,self.multiplotAct,self.setrefAct,
                            self.delrefAct,self.userxAct,self.plotprefAct]
        
        # the history drop down menu
        self.historyMenu = QtWidgets.QMenu("&History", self)
        self.menubar.addMenu(self.historyMenu)
        self.historyAct = self.historyMenu.addAction(QtGui.QIcon(IconDirectory + 'history.png'), "&History", lambda: self.mainWindowCheck(lambda mainWindow: HistoryWindow(mainWindow)))
        self.historyAct.setToolTip('Show Processing History')
        self.errorAct = self.historyMenu.addAction(QtGui.QIcon(IconDirectory + 'error.png'),"&Error Messages", lambda: errorWindow(self))
        self.errorAct.setToolTip('Show Error Messages')
        
        self.historyActList = [self.historyAct]
        
        #Utilities dropdown menu
        self.utilitiesMenu = QtWidgets.QMenu("&Utilities", self)
        self.menubar.addMenu(self.utilitiesMenu)
        self.shiftconvAct = self.utilitiesMenu.addAction(QtGui.QIcon(IconDirectory + 'shifttool.png'),"&Chemical Shift Conversion Tool", self.createShiftConversionWindow)
        self.shiftconvAct.setToolTip('Chemical Shift Conversion Tool')
        self.quadconvAct = self.utilitiesMenu.addAction(QtGui.QIcon(IconDirectory + 'quadconversion.png'),"&Quadrupole Coupling Conversion Tool", self.createQuadConversionWindow)
        self.quadconvAct.setToolTip('Quadrupole Coupling Conversion Tool')
        self.nmrtableAct = self.utilitiesMenu.addAction(QtGui.QIcon(IconDirectory + 'table.png'),"&NMR Table", self.nmrTable)
        self.nmrtableAct.setToolTip('NMR Periodic Table')
        self.utilitiesActList = [self.shiftconvAct,self.quadconvAct,self.nmrtableAct]
        
        # the help drop down menu
        self.helpMenu = QtWidgets.QMenu("&Help", self)
        self.menubar.addMenu(self.helpMenu)
        self.updateAct = self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'update.png'),"&Update", self.updateMenu)
        self.updateAct.setToolTip('Update ssNake')
        self.aboutAct = self.helpMenu.addAction(QtGui.QIcon(IconDirectory + 'about.png'),"&About", lambda: aboutWindow(self))
        self.aboutAct.setToolTip('About Menu') 

        self.helpActList = [self.updateAct,self.shiftconvAct,self.quadconvAct,self.nmrtableAct,self.aboutAct]
    

        #Extra event lists:
        self.specOnlyList = [self.regridAct,self.csastaticAct,self.firstquadstatAct,self.secondquadstatAct,self.czjzekstatAct]
        self.fidOnlyList = [self.relaxAct,self.diffusionAct]


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

            self.macrolistmenu.menuAction().setEnabled(False)
            #self.macromenu.menuAction().setEnabled(False)
            self.editmenu.menuAction().setEnabled(False)
            self.toolMenu.menuAction().setEnabled(False)
            self.matrixMenu.menuAction().setEnabled(False)
            self.fftMenu.menuAction().setEnabled(False)
            self.fittingMenu.menuAction().setEnabled(False)
            self.combineMenu.menuAction().setEnabled(False)
            #self.plotMenu.menuAction().setEnabled(False)
            self.referencerunmenu.menuAction().setEnabled(False)
            #self.historyMenu.menuAction().setEnabled(False)
            for act in self.saveActList + self.exportActList + self.workspaceActList + self.macroActList + self.editActList + self.toolsActList + self.matrixActList + self.fftActList + self.fittingActList + self.plotActList + self.combineActList + self.historyActList:
                act.setEnabled(False)
 
        else:
            self.editmenu.menuAction().setEnabled(True)
            self.toolMenu.menuAction().setEnabled(True)
            self.matrixMenu.menuAction().setEnabled(True)
            self.fftMenu.menuAction().setEnabled(True)
            self.fittingMenu.menuAction().setEnabled(True)
            self.combineMenu.menuAction().setEnabled(True)
            self.referencerunmenu.menuAction().setEnabled(True)

            for act in self.editActList + self.toolsActList + self.matrixActList + self.fftActList + self.fittingActList + self.plotActList + self.historyActList + self.combineActList:
                act.setEnabled(True)
            if isinstance(self.mainWindow, Main1DWindow):
                self.menuEnable()

                for act in self.specOnlyList:
                    act.setEnabled(self.mainWindow.current.spec == 1) #Only on for spec
                for act in self.fidOnlyList:
                    act.setEnabled(self.mainWindow.current.spec == 0) #Only on for FID
                
                if self.mainWindow.masterData.noUndo: #Set menu check to the same value as in the data
                    self.noUndoAct.setChecked(True)
                else:
                    self.noUndoAct.setChecked(False)
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
                if not self.mainWindow.undoList and not self.mainWindow.redoList:
                    self.clearundoAct.setEnabled(False)
                self.savemenu.menuAction().setEnabled(True)
                self.exportmenu.menuAction().setEnabled(True)
                self.savefigAct.setEnabled(True)
                #self.macromenu.menuAction().setEnabled(True)
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

            elif isinstance(self.mainWindow, MainPlotWindow):
                self.menuDisable(True)
                self.savemenu.menuAction().setEnabled(True)
                self.exportmenu.menuAction().setEnabled(True)
                for act in self.saveActList + self.exportActList + self.workspaceActList:
                    act.setEnabled(True)
                self.savefigAct.setEnabled(False)
                self.workspacemenu.menuAction().setEnabled(True)
                #self.macromenu.menuAction().setEnabled(False)
                self.macrolistmenu.menuAction().setEnabled(False)
            else:
                self.menuDisable(True)
                self.savemenu.menuAction().setEnabled(True)
                self.exportmenu.menuAction().setEnabled(True)
                self.savefigAct.setEnabled(True)
                self.workspacemenu.menuAction().setEnabled(True)
                #self.macromenu.menuAction().setEnabled(False)
                self.macrolistmenu.menuAction().setEnabled(False)
                for act in self.saveActList + self.exportActList + self.workspaceActList:
                    act.setEnabled(True)

    def menuEnable(self, internalWindow=False):
        #self.macromenu.menuAction().setEnabled(True)
        self.macrolistmenu.menuAction().setEnabled(True)
        self.editmenu.menuAction().setEnabled(True)
        self.toolMenu.menuAction().setEnabled(True)
        self.matrixMenu.menuAction().setEnabled(True)
        self.fftMenu.menuAction().setEnabled(True)
        self.fittingMenu.menuAction().setEnabled(True)
        self.combineMenu.menuAction().setEnabled(True)
        #self.plotMenu.menuAction().setEnabled(True)
        self.referencerunmenu.menuAction().setEnabled(True)
        #self.historyMenu.menuAction().setEnabled(True)
        
        #Actions:
        for act in self.macroActList + self.editActList + self.toolsActList + self.matrixActList + self.fftActList + self.fittingActList + self.plotActList + self.combineActList + self.historyActList:
            act.setEnabled(True)
        
        
        if not internalWindow:
            self.filemenu.menuAction().setEnabled(True)
            self.workspacemenu.menuAction().setEnabled(True)
            for act in self.fileActList + self.workspaceActList:
                act.setEnabled(True)
            for i in range(self.tabs.count()):
                self.tabs.setTabEnabled(i, True)
        self.undoAction.setEnabled(True)
        self.redoAction.setEnabled(True)

    def menuDisable(self, internalWindow=False):
        #self.macromenu.menuAction().setEnabled(False)
        self.macrolistmenu.menuAction().setEnabled(True)
        self.editmenu.menuAction().setEnabled(False)
        self.toolMenu.menuAction().setEnabled(False)
        self.matrixMenu.menuAction().setEnabled(False)
        self.fftMenu.menuAction().setEnabled(False)
        self.fittingMenu.menuAction().setEnabled(False)
        self.combineMenu.menuAction().setEnabled(False)
        #self.plotMenu.menuAction().setEnabled(False)
        self.referencerunmenu.menuAction().setEnabled(False)
        #self.historyMenu.menuAction().setEnabled(False)
        #Actions:
        for act in self.macroActList + self.editActList + self.toolsActList + self.matrixActList + self.fftActList + self.fittingActList + self.plotActList + self.combineActList + self.historyActList:
            act.setEnabled(False)
            
            
        if not internalWindow:
            self.filemenu.menuAction().setEnabled(False)
            self.workspacemenu.menuAction().setEnabled(False)
            for act in self.fileActList + self.workspaceActList:
                act.setEnabled(False)
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
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        action1 = self.macrolistmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'),givenName, lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(givenName, lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'),givenName, lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'),givenName, lambda name=givenName: self.renameMacro(name))
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
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        action1 = self.macrolistmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'),givenName, lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(givenName, lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'),givenName, lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'),givenName, lambda name=givenName: self.renameMacro(name))
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
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        action1 = self.macrolistmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'),givenName, lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(givenName, lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'),givenName, lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'),givenName, lambda name=givenName: self.renameMacro(name))
        self.macroActions[givenName] = [action1, action2, action3, action4]
        self.menuCheck()

    def referenceAdd(self, reffreq, name):
        self.referenceName.append(name)
        self.referenceValue.append(reffreq)  # List with saved refrence values
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        action1 = self.referencerunmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'),name, lambda name=name: self.referenceRun(name))
        action2 = self.referencedeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'),name, lambda name=name: self.referenceRemove(name))
        action3 = self.referencerenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'),name, lambda name=name: self.referenceRename(name))
        action4 = self.referencesavemenu.addAction(name, lambda name=name: self.referenceSave(name))
        self.referenceActions[name] = [action1, action2, action3, action4]
        self.menuCheck()

    def referenceClear(self):
        if self.mainWindow.masterData.noUndo:
            self.mainWindow.current.setRef(None)
        else:
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

        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        action1 = self.referencerunmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'),givenName, lambda name=givenName: self.referenceRun(name))
        action2 = self.referencedeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'),givenName, lambda name=givenName: self.referenceRemove(name))
        action3 = self.referencerenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'),givenName, lambda name=givenName: self.referenceRename(name))
        action4 = self.referencesavemenu.addAction(givenName, lambda name=givenName: self.referenceSave(name))


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
        IconDirectory = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'Icons' + os.path.sep
        self.referenceValue.append(freq)
        action1 = self.referencerunmenu.addAction(QtGui.QIcon(IconDirectory + 'run.png'),givenName, lambda name=givenName: self.referenceRun(name))
        action2 = self.referencedeletemenu.addAction(QtGui.QIcon(IconDirectory + 'delete.png'),givenName, lambda name=givenName: self.referenceRemove(name))
        action3 = self.referencerenamemenu.addAction(QtGui.QIcon(IconDirectory + 'rename.png'),givenName, lambda name=givenName: self.referenceRename(name))
        action4 = self.referencesavemenu.addAction(givenName, lambda name=givenName: self.referenceSave(name))
        self.referenceActions[givenName] = [action1, action2, action3, action4]
        self.menuCheck()

    def noUndoMode(self,val):
        if val:
            self.mainWindow.undoList = []
            self.mainWindow.redoList = []
            self.mainWindow.masterData.noUndo = True
        else:
            self.mainWindow.masterData.noUndo = False
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
        try:
            if isinstance(self.mainWindow.current, (sc.CurrentMulti)):
                self.mainWindow.sideframe.checkChanged()
        except:
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
            if num == len(self.workspaces):
                self.workspaceNum = num - 1
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

    def createCombineLoadWindow(self):
        CombineLoadWindow(self)

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
                self.dispMsg("Not all the data has the same shape")
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
                        if self.autoLoad(os.path.join(temp_dir, i), realpath=filePath):
                            break
                finally:
                    shutil.rmtree(temp_dir)
            else:
                self.autoLoad(filePath)

    def fileTypeCheck(self, filePath):
        returnVal = 0
        
        
        fileBase = '' 
        direc = filePath 
        if os.path.isfile(filePath):
            filename = os.path.basename(filePath)
            fileBase = os.path.splitext(filename)[0]
            direc = os.path.dirname(filePath)
            if filename.endswith('.fid') or filename.endswith('.spe'):
                with open(filePath, 'r') as f:
                    check = int(np.fromfile(f, np.float32, 1))
                if check == 0:
                    return (8, filePath, returnVal)  # Suspected NMRpipe format
                else:  # SIMPSON
                    return (4, filePath, returnVal)
            elif filename.endswith('.json') or filename.endswith('.JSON'):
                return (5, filePath, returnVal)
            elif filename.endswith('.mat') or filename.endswith('.MAT'):
                return (6, filePath, returnVal)
            elif filename.endswith('.jdf'):#JEOL delta format
                return (9, filePath, returnVal)
            elif filename.endswith('.dx') or filename.endswith('.jdx') or filename.endswith('.jcamp'):#JCAMP format
                return (10, filePath, returnVal)
            elif filename.endswith('.sig'): #Bruker minispec    
                return (12, filePath, returnVal)
            returnVal = 1
            direc = os.path.dirname(filePath)

        if os.path.exists(direc + os.path.sep + 'procpar') and os.path.exists(direc + os.path.sep + 'fid'):
            return (0, direc, returnVal)
            # And for varian processed data
        if (os.path.exists(direc + os.path.sep + '..' + os.path.sep + 'procpar') or os.path.exists(direc + os.path.sep + 'procpar')) and os.path.exists(direc + os.path.sep + 'data'):
            return (0, direc, returnVal)
        elif os.path.exists(direc + os.path.sep + 'acqus') and (os.path.exists(direc + os.path.sep + 'fid') or os.path.exists(direc + os.path.sep + 'ser')):
            return (1, direc, returnVal)
        elif os.path.exists(direc + os.path.sep + 'procs') and (os.path.exists(direc + os.path.sep + '1r') or os.path.exists(direc + os.path.sep + '2rr')):
            return (7, direc, returnVal)
        elif os.path.exists(direc + os.path.sep + 'acq') and os.path.exists(direc + os.path.sep + 'data'):
            return (2, direc, returnVal)
        elif os.path.exists(direc + os.path.sep + 'acqu.par'):
            dirFiles = os.listdir(direc)
            files2D = [x for x in dirFiles if '.2d' in x]
            files1D = [x for x in dirFiles if '.1d' in x]
            if len(files2D) != 0 or len(files1D) != 0:
                return (3, direc, returnVal)
        elif os.path.exists(direc + os.path.sep + fileBase + '.spc') and os.path.exists(direc + os.path.sep + fileBase + '.par'):
            return (13, direc + os.path.sep + fileBase, returnVal)
        elif os.path.isfile(filePath): #If not recognised, load as ascii
            return (11, filePath, returnVal)
        return (None,filePath, 2)
                
    def autoLoad(self, filePath, realpath=False):
        val = self.fileTypeCheck(filePath)
        if val[0] is not None:
            self.loading(val[0], val[1], realpath=realpath)
        return val[2]

    def loadAndCombine(self, filePathList):
        filePath = filePathList.pop(0)
        combineMasterData = None
        if filePath.endswith('.zip'):
            import tempfile
            import shutil
            import zipfile
            try:
                temp_dir = tempfile.mkdtemp()
                zipfile.ZipFile(filePath).extractall(temp_dir)
                val = self.fileTypeCheck(os.path.join(temp_dir, os.listdir(temp_dir)[0]))
                combineMasterData = self.loading(val[0], val[1], returnBool=True, realpath=filePath)
            finally:
                shutil.rmtree(temp_dir)
        else:
            val = self.fileTypeCheck(filePath)
            combineMasterData = self.loading(val[0], val[1], returnBool=True)
        if combineMasterData is None:
            self.dispMsg("Data could not be loaded")
            return False
        shapeRequired = combineMasterData.data.shape
        combineMasterData.split(1, -1)
        for filePath in filePathList:
            if filePath.endswith('.zip'):
                import tempfile
                import shutil
                import zipfile
                try:
                    temp_dir = tempfile.mkdtemp()
                    zipfile.ZipFile(filePath).extractall(temp_dir)
                    val = self.fileTypeCheck(os.path.join(temp_dir, os.listdir(temp_dir)[0]))
                    addData = self.loading(val[0], val[1], returnBool=True, realpath=filePath)
                finally:
                    shutil.rmtree(temp_dir)
            else:
                val = self.fileTypeCheck(filePath)
                addData = self.loading(val[0], val[1], returnBool=True)
            if addData.data.shape != shapeRequired:
                self.dispMsg("Not all the data has the required shape")
                return False
            combineMasterData.insert(addData.data, combineMasterData.data.shape[0], 0)
        wsname = self.askName()
        self.workspaces.append(Main1DWindow(self, combineMasterData))
        self.workspaces[-1].rename(wsname)
        self.tabs.addTab(self.workspaces[-1], wsname)
        self.workspaceNames.append(wsname)
        self.changeMainWindow(wsname)
    
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

    def loading(self, num, filePath, returnBool=False, realpath=False):
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
        elif num == 10:
            masterData = self.LoadJCAMP(filePath, name)   
        elif num == 11:
            masterData = self.LoadAscii(filePath, name) 
        elif num == 12:
            masterData = self.LoadMinispec(filePath, name)
        elif num == 13:
            masterData = self.LoadBrukerEPR(filePath, name)
        if returnBool:
            return masterData
        else:
            if masterData is not None:
                self.workspaces.append(Main1DWindow(self, masterData))
                self.tabs.addTab(self.workspaces[-1], name)
                self.workspaceNames.append(name)
                self.changeMainWindow(name)

    def LoadVarianFile(self, filePath, name=''):
        try:
            masterData = LF.LoadVarianFile(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            return masterData
        except:
            self.dispMsg("Error on loading Varian data",'red')
            return None 

    def LoadPipe(self, filePath, name=''):
        try:
            masterData = LF.LoadPipe(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            return masterData
        except:
            self.dispMsg("Error on loading NMRpipe data",'red')
            return None 

    def LoadJEOLDelta(self, filePath, name=''):
        try:
            masterData = LF.LoadJEOLDelta(filePath,name)
        except:
            self.dispMsg("Error on loading JEOL Delta data",'red')
            return None
        if masterData ==  'ND error':
            self.dispMsg("Error: JEOL Delta data of this type is not supported",'red')
        masterData.msgHandler = lambda msg: self.dispMsg(msg)
        return masterData


    def loadJSONFile(self, filePath, name=''):
        try:
            masterData = LF.loadJSONFile(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            return masterData
        except:
            self.dispMsg("Error on loading JSON data",'red')
            return None 

    def loadMatlabFile(self, filePath, name=''):
        try:
            masterData = LF.loadMatlabFile(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            return masterData
        except:
            self.dispMsg("Error on loading MATLAB data",'red')
            return None 

    def LoadBrukerTopspin(self, filePath, name=''):
        try:
            masterData = LF.LoadBrukerTopspin(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            return masterData
        except:
            self.dispMsg("Error on loading Bruker data",'red')
            return None 
        
        
    def LoadBrukerSpectrum(self, filePath, name=''):
        try:
            masterData = LF.LoadBrukerSpectrum(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            return masterData     
        except:
            self.dispMsg("Error on loading Bruker Spectrum data",'red')
            return None 

    def LoadChemFile(self, filePath, name=''):
        try:
            masterData = LF.LoadChemFile(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            return masterData   
        except:
            self.dispMsg("Error on loading Chemagnetic data",'red')
            return None

    def LoadMagritek(self, filePath, name='',realPath=''):
        try:
            masterData = LF.LoadMagritek(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            return masterData  
        except:
            self.dispMsg("Error on loading Magritek data",'red')
            return None 

    def LoadSimpsonFile(self, filePath, name=''):
        try:
            masterData = LF.LoadSimpsonFile(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            return masterData  
        except:
            self.dispMsg("Error on loading SIMPSON data",'red')
            return None 
    
    def LoadJCAMP(self, filePath, name=''):
        try:
            masterData = LF.loadJCAMP(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            masterData.addHistory("JCAMP data loaded from " + filePath)
            return masterData
        except:
            self.dispMsg("Error on loading JCAMP data",'red')
            return None 
        
    def LoadAscii(self, filePath, name=''):
        dialog = AsciiLoadWindow(self, filePath)
        if dialog.exec_():
            if dialog.closed:
                return
            else:
                try:
                    masterData = LF.LoadAscii(filePath, name, dialog.dataDimension, dialog.dataSpec, dialog.dataOrder, dialog.delim, dialog.sw)
                    masterData.msgHandler = lambda msg: self.dispMsg(msg)
                    return masterData
                except:
                    self.dispMsg("Error on loading ASCII data",'red')
                    return None
                    
    def LoadMinispec(self, filePath, name=''):
        try:
            masterData = LF.LoadMinispec(filePath,name)

            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            masterData.addHistory("Minispec data loaded from " + filePath)
            return masterData 
        except:
            self.dispMsg("Error on loading Minispec data",'red')
            return None 
          
    def LoadBrukerEPR(self, filePath, name=''):
        try:
            masterData = LF.LoadBrukerEPR(filePath,name)
            masterData.msgHandler = lambda msg: self.dispMsg(msg)
            return masterData     
        except:
            self.dispMsg("Error on loading Bruker EPR data",'red')
            return None 

   
        
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
#        self.mainWindow.current.showFid()

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
        super(Main1DWindow, self).__init__(father)
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
        grid.addWidget(self.bottomframe, 1, 0 , 1, 2)
        self.textframe = TextFrame(self)
        grid.addWidget(self.textframe, 2, 0 , 1 , 2)
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
                returnValue = self.masterData.restoreData(loadData, None)
            elif iter1[0] == 'real':
                returnValue = self.masterData.real()
            elif iter1[0] == 'imag':
                returnValue = self.masterData.imag()
            elif iter1[0] == 'abs':
                returnValue = self.masterData.abs()
            elif iter1[0] == 'conj':
                returnValue = self.masterData.conj()
            elif iter1[0] == 'phase':
                returnValue = self.masterData.setPhase(*iter1[1])
            elif iter1[0] == 'autoPhase':
                returnValue = self.masterData.autoPhase(*iter1[1])
            elif iter1[0] == 'fourier':
                returnValue = self.masterData.fourier(*iter1[1])
            elif iter1[0] == 'realFourier':
                returnValue = self.masterData.realFourier(*iter1[1])
            elif iter1[0] == 'fftshift':
                returnValue = self.masterData.fftshift(*iter1[1])
            elif iter1[0] == 'diff':
                returnValue = self.masterData.diff(*iter1[1])
            elif iter1[0] == 'cumsum':
                returnValue = self.masterData.cumsum(*iter1[1])
            elif iter1[0] == 'apodize':
                returnValue = self.masterData.apodize(*iter1[1])
            elif iter1[0] == 'freq':
                returnValue = self.masterData.setFreq(*iter1[1])
            elif iter1[0] == 'ref':
                returnValue = self.masterData.setRef(*iter1[1])
            elif iter1[0] == 'size':
                returnValue = self.masterData.setSize(*iter1[1])
            elif iter1[0] == 'spec':
                returnValue = self.masterData.changeSpec(*iter1[1])
            elif iter1[0] == 'swapecho':
                returnValue = self.masterData.swapEcho(*iter1[1])
            elif iter1[0] == 'wholeEcho':
                returnValue = self.masterData.wholeEcho(*iter1[1])
            elif iter1[0] == 'shift':
                returnValue = self.masterData.shiftData(*iter1[1])
            elif iter1[0] == 'states':
                returnValue = self.masterData.states(*iter1[1])
            elif iter1[0] == 'statesTPPI':
                returnValue = self.masterData.statesTPPI(*iter1[1])
            elif iter1[0] == 'echoAntiEcho':
                returnValue = self.masterData.echoAntiEcho(*iter1[1])
            elif iter1[0] == 'baselineCorrection':
                returnValue = self.masterData.baselineCorrection(*iter1[1])
            elif iter1[0] == 'integrate':
                returnValue = self.masterData.matrixManip(*iter1[1], which=0)
            elif iter1[0] == 'sum':
                returnValue = self.masterData.matrixManip(*iter1[1], which=5)
            elif iter1[0] == 'max':
                returnValue = self.masterData.matrixManip(*iter1[1], which=1)
            elif iter1[0] == 'min':
                returnValue = self.masterData.matrixManip(*iter1[1], which=2)
            elif iter1[0] == 'argmax':
                returnValue = self.masterData.matrixManip(*iter1[1], which=3)
            elif iter1[0] == 'argmin':
                returnValue = self.masterData.matrixManip(*iter1[1], which=4)
            elif iter1[0] == 'average':
                returnValue = self.masterData.matrixManip(*iter1[1], which=6)
            elif iter1[0] == 'regrid':
                returnValue = self.masterData.regrid(*iter1[1])
            elif iter1[0] == 'fliplr':
                returnValue = self.masterData.flipLR(*iter1[1])
            elif iter1[0] == 'concatenate':
                returnValue = self.masterData.concatenate(*iter1[1])
            elif iter1[0] == 'split':
                returnValue = self.masterData.split(*iter1[1])
            elif iter1[0] == 'insert':
                returnValue = self.masterData.insert(*iter1[1])
            elif iter1[0] == 'delete':
                returnValue = self.masterData.remove(*iter1[1])
            elif iter1[0] == 'add':
                returnValue = self.masterData.add(*iter1[1])
            elif iter1[0] == 'subtract':
                returnValue = self.masterData.subtract(*iter1[1])
            elif iter1[0] == 'multiplySpec':
                returnValue = self.masterData.multiplySpec(*iter1[1])
            elif iter1[0] == 'divideSpec':
                returnValue = self.masterData.divideSpec(*iter1[1])
            elif iter1[0] == 'multiply':
                returnValue = self.masterData.multiply(*iter1[1])
            elif iter1[0] == 'subtractAvg':
                returnValue = self.masterData.subtractAvg(*iter1[1])
            elif iter1[0] == 'FIDDLE':
                returnValue = self.masterData.fiddle(*iter1[1])
            elif iter1[0] == 'reorder':
                returnValue = self.masterData.reorder(*iter1[1])
            elif iter1[0] == 'ffm':
                returnValue = self.masterData.ffm_1d(*iter1[1])
            elif iter1[0] == 'clean':
                returnValue = self.masterData.clean(*iter1[1])
            elif iter1[0] == 'shear':
                returnValue = self.masterData.shear(*iter1[1])
            elif iter1[0] == 'extract':
                returnValue = self.masterData.getRegion(*iter1[1])
            elif iter1[0] == 'setxax':
                returnValue = self.masterData.setXax(*iter1[1])
            elif iter1[0] == 'hilbert':
                returnValue = self.masterData.hilbert(*iter1[1])
            else:
                self.father.dispMsg('unknown macro command: ' + iter1[0])
                returnValue = None

            if not self.masterData.noUndo or returnValue is not None:
                self.undoList.append(returnValue)
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
        #axis = np.array([self.masterData.xaxArray[-1]]).transpose()
        axType = self.current.axType
        if self.masterData.spec[-1] == 1:
            if self.current.ppm:
                if self.current.ref is not None:
                    axMult = 1e6 / self.masterData.ref[-1]
                else:                    
                    axMult = 1e6 / self.masterData.freq[-1]
            else:
                axMult = 1.0 / (1000.0**axType)
        elif self.masterData.spec[-1] == 0:
            axMult = 1000.0**axType

        axis = np.array([self.masterData.xaxArray[-1] * axMult]).transpose()


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
        if self.masterData.noUndo:
            self.masterData.restoreData(loadData, None)
        else:
            self.undoList.append(self.masterData.restoreData(loadData, None))
        self.current.upd()
        self.current.plotReset()
        self.current.showFid()
        self.updAllFrames()
        self.addMacro(['reload'])
        self.menuCheck()

    def monitorLoad(self, filePath, delay = 0.5):
        self.monitor.blockSignals(True)
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
        QtCore.QTimer.singleShot(delay * 1000, lambda: self.monitor.blockSignals(False))
        if filePath in self.monitor.files() or filePath in self.monitor.directories():
            return
        self.monitor.addPath(filePath)

    def startMonitor(self, macroNames, delay = 0.5):
        self.monitorMacros = macroNames
        self.monitor = QtCore.QFileSystemWatcher([self.masterData.filePath[1]], self)
        self.monitor.fileChanged.connect(lambda a: self.monitorLoad(a,delay))
        self.monitor.directoryChanged.connect(lambda a: self.monitorLoad(a,delay))

    def stopMonitor(self):
        self.monitorMacros = []
        if self.monitor is not None:
            self.monitor.removePath(self.masterData.filePath[1])
        self.monitor = None
        
    def real(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.masterData.real()
        else:
            self.undoList.append(self.masterData.real())
        self.current.upd()
        self.current.showFid()
        self.addMacro(['real'])
        self.menuCheck()

    def imag(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.masterData.imag()
        else:
            self.undoList.append(self.masterData.imag())
        self.current.upd()
        self.current.showFid()
        self.addMacro(['imag'])
        self.menuCheck()

    def abs(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.masterData.abs()
        else:
            self.undoList.append(self.masterData.abs())
        self.current.upd()
        self.current.showFid()
        self.addMacro(['abs'])
        self.menuCheck()
        
    def conj(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.masterData.conj()
        else:
            self.undoList.append(self.masterData.conj())
        self.current.upd()
        self.current.showFid()
        self.addMacro(['conj'])
        self.menuCheck()

    def fourier(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.fourier()
        else:
            self.undoList.append(self.current.fourier())
        self.bottomframe.upd()
        self.menuCheck()

    def realFourier(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.realFourier()
        else:
            self.undoList.append(self.current.realFourier())
        self.bottomframe.upd()
        self.menuCheck()

    def fftshift(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.fftshift()
        else:
            self.undoList.append(self.current.fftshift())
        self.updAllFrames()
        self.menuCheck()

    def invFftshift(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.fftshift(inv=True)
        else:
            self.undoList.append(self.current.fftshift(inv=True))
        self.updAllFrames()
        self.menuCheck()

    def diff(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.diff()
        else:
            self.undoList.append(self.current.diff())
        self.updAllFrames()
        self.menuCheck()

    def cumsum(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.cumsum()
        else:
            self.undoList.append(self.current.cumsum())
        self.updAllFrames()
        self.menuCheck()

    def hilbert(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.hilbert()
        else:
            self.undoList.append(self.current.hilbert())
        self.menuCheck()

    def states(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.states()
        else:
            self.undoList.append(self.current.states())
        self.updAllFrames()
        self.menuCheck()

    def statesTPPI(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.statesTPPI()
        else:
            self.undoList.append(self.current.statesTPPI())
        self.updAllFrames()
        self.menuCheck()

    def echoAntiEcho(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.echoAntiEcho()
        else:
            self.undoList.append(self.current.echoAntiEcho())
        self.updAllFrames()
        self.menuCheck()

    def setFreq(self, freq, sw):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.setFreq(freq, sw)
        else:
            self.undoList.append(self.current.setFreq(freq, sw))
        self.menuCheck()

    def flipLR(self):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.flipLR()
        else:
            self.undoList.append(self.current.flipLR())
        self.menuCheck()

    def directAutoPhase(self, phaseNum):
        self.redoList = []
        if self.masterData.noUndo:
            self.current.directAutoPhase(phaseNum)
        else:
            self.undoList.append(self.current.directAutoPhase(phaseNum))
        self.menuCheck()

    def BrukerDigital(self):
#        FilePath = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.LastLocation)
#        if type(FilePath) is tuple:
#            FilePath = FilePath[0]
#        self.father.LastLocation = os.path.dirname(FilePath)  # Save used path
        FilePath = self.masterData.filePath[1]
        if FilePath is '':
            return
        Dir = os.path.dirname(FilePath)
        if not os.path.exists(Dir + os.path.sep + 'acqus'):
            self.father.dispMsg("acqus file does not exist, specify load path")
            FilePath = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.LastLocation)[0]
            if FilePath == '':
                return
            self.father.LastLocation = os.path.dirname(FilePath)  # Save used path
            Dir = os.path.dirname(FilePath)
            if type(FilePath) is tuple:
                FilePath = FilePath[0]
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
            if self.masterData.noUndo:
                self.current.applyPhase(0, FilterCorrection * 2 * np.pi)
            else:   
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

    def createQuad1DeconvWindow(self):
        self.father.createFitWindow(fit.Quad1DeconvWindow(self.father, self.father.mainWindow))

    def createQuad2DeconvWindow(self):
        if self.current.freq == 0.0:
            self.father.dispMsg("Please set the spectrometer frequency first!")
            return
        self.father.createFitWindow(fit.Quad2DeconvWindow(self.father, self.father.mainWindow))

    def createQuad2CzjzekWindow(self):
        if self.current.freq == 0.0:
            self.father.dispMsg("Please set the spectrometer frequency first!")
            return
        self.father.createFitWindow(fit.Quad2CzjzekWindow(self.father, self.father.mainWindow))

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
        message = self.masterData.removeFromHistory(2)
        self.father.dispMsg("Undo: " + message,error = False)
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
                
                self.contourTypeGroup = QtWidgets.QGroupBox('Contour type:')
                self.contourTypeFrame = QtWidgets.QGridLayout()
                self.contourNumberLabel = wc.QLeftLabel("Number:", self)
                
                self.contourTypeFrame.addWidget(self.contourNumberLabel, 0, 0)
                
                self.numLEntry = QtWidgets.QSpinBox()
                self.numLEntry.setMaximum(100000)
                self.numLEntry.setMinimum(1)
                self.numLEntry.setValue(current.numLevels)
                self.numLEntry.valueChanged.connect(self.setContour)
                self.contourTypeFrame.addWidget(self.numLEntry, 0, 1)
                
                self.contourTypeFrame.addWidget(wc.QLeftLabel("Sign:", self), 1, 0)
                self.contourSignEntry = QtWidgets.QComboBox()
                self.contourSignEntry.addItems(['Both','+ only','- only'])
                self.contourSignEntry.setCurrentIndex(current.contourSign)
                self.contourSignEntry.currentIndexChanged.connect(self.setContour)
                self.contourTypeFrame.addWidget(self.contourSignEntry, 1, 1)
                
                self.contourTypeLabel = wc.QLeftLabel("Type:", self)
                self.contourTypeFrame.addWidget(self.contourTypeLabel, 2, 0)
                
                self.contourTypeEntry = QtWidgets.QComboBox()
                self.contourTypeEntry.addItems(['Linear','Multiplier'])
                self.contourTypeEntry.setCurrentIndex(current.contourType)
                self.contourTypeEntry.currentIndexChanged.connect(self.setContour)
                self.contourTypeFrame.addWidget(self.contourTypeEntry, 2, 1)
                
                self.multiValueLabel = wc.QLeftLabel("Multiplier:", self)
                self.contourTypeFrame.addWidget(self.multiValueLabel, 3, 0)
                
                self.multiValue = wc.QLineEdit(current.multiValue, self.setContour)
                self.multiValue.setMaximumWidth(120)
                self.contourTypeFrame.addWidget(self.multiValue, 3, 1)
                
                if current.contourType != 1:
                    self.multiValueLabel.hide()
                    self.multiValue.hide()
                self.contourTypeGroup.setLayout(self.contourTypeFrame)
                self.frame2.addWidget(self.contourTypeGroup, 6, 0, 1, 3)
                    
                #Contour limits    
                self.contourLimitsGroup = QtWidgets.QGroupBox('Contour limits [%]:')
                self.contourLimitsFrame = QtWidgets.QGridLayout()
                self.maxLEntry = wc.QLineEdit(format(current.maxLevels * 100.0, '.7g'), self.setContour)
                self.maxLEntry.setMaximumWidth(120)
                self.contourLimitsFrame.addWidget(self.maxLEntry, 0, 1)
                self.minLEntry = wc.QLineEdit(format(current.minLevels * 100.0, '.7g'), self.setContour)
                self.minLEntry.setMaximumWidth(120)
                self.contourLimitsFrame.addWidget(self.minLEntry, 1, 1)
                self.maxLabel = wc.QLeftLabel("Max:", self)
                self.minLabel = wc.QLeftLabel("Min:", self)
                self.contourLimitsFrame.addWidget(self.maxLabel, 0, 0)
                self.contourLimitsFrame.addWidget(self.minLabel, 1, 0)
                self.contourLimitsGroup.setLayout(self.contourLimitsFrame)
                self.frame2.addWidget(self.contourLimitsGroup, 7, 0, 1, 3)
                
                #Projections
                self.contourProjGroup = QtWidgets.QGroupBox('Projections:')
                self.contourProjFrame = QtWidgets.QGridLayout()
                self.projTopLabel = wc.QLeftLabel("Top:", self)
                self.contourProjFrame.addWidget(self.projTopLabel, 0, 0)
                self.projDropTop = QtWidgets.QComboBox()
                self.projDropTop.addItems(["sum", "max", "min" , "off"])
                self.projDropTop.setCurrentIndex(current.projTop)
                self.projDropTop.activated.connect(lambda val, self=self: self.changeProj(val, 1))
                self.contourProjFrame.addWidget(self.projDropTop, 0, 1,)
           
                self.projRightLabel = wc.QLeftLabel("Right:", self)
                self.contourProjFrame.addWidget(self.projRightLabel, 1, 0)
                self.projDropRight = QtWidgets.QComboBox()
                self.projDropRight.addItems(["sum", "max", "min" ,"off"])
                self.projDropRight.setCurrentIndex(current.projRight)
                self.projDropRight.activated.connect(lambda val, self=self: self.changeProj(val, 2))
                self.contourProjFrame.addWidget(self.projDropRight, 1, 1)
                
                #Ranges
                self.rangeCheckbox = QtWidgets.QCheckBox('Projection ranges',self)
                self.rangeCheckbox.setChecked(current.projLimitsBool)
                self.rangeCheckbox.stateChanged.connect(self.activateRanges)
                self.contourProjFrame.addWidget(self.rangeCheckbox,2,0,1,2)
                
                self.projTopRangeMaxLabel = wc.QLeftLabel("Top max:", self)
                self.projTopRangeMaxLabel.hide()
                self.contourProjFrame.addWidget(self.projTopRangeMaxLabel, 3, 0)
                self.projTopRangeMax = QtWidgets.QSpinBox()
                self.projTopRangeMax.setMaximum(self.shape[current.axes2] - 1)
                self.projTopRangeMax.setMinimum(0)
                if current.projLimits[0] is None:
                    self.projTopRangeMax.setValue(self.shape[current.axes2] - 1)
                else:
                    self.projTopRangeMax.setValue(current.projLimits[0])
                self.projTopRangeMax.valueChanged.connect(self.changeRanges)
                self.projTopRangeMax.hide()
                self.contourProjFrame.addWidget(self.projTopRangeMax, 3, 1)
                
                self.projTopRangeMinLabel = wc.QLeftLabel("Top min:", self)
                self.projTopRangeMinLabel.hide()
                self.contourProjFrame.addWidget(self.projTopRangeMinLabel, 4, 0)
                self.projTopRangeMin = QtWidgets.QSpinBox()
                self.projTopRangeMin.setMaximum(self.shape[current.axes2] - 1)
                self.projTopRangeMin.setMinimum(0)
                if current.projLimits[1] is None:
                    self.projTopRangeMin.setValue(0)
                else:
                    self.projTopRangeMin.setValue(current.projLimits[1])
                self.projTopRangeMin.valueChanged.connect(self.changeRanges)
                self.projTopRangeMin.hide()
                self.contourProjFrame.addWidget(self.projTopRangeMin, 4, 1)
                
                self.projRightRangeMaxLabel = wc.QLeftLabel("Right max:", self)
                self.projRightRangeMaxLabel.hide()
                self.contourProjFrame.addWidget(self.projRightRangeMaxLabel, 5, 0)
                self.projRightRangeMax = QtWidgets.QSpinBox()
                self.projRightRangeMax.setMaximum(self.shape[current.axes] - 1)
                self.projRightRangeMax.setMinimum(0)
                if current.projLimits[2] is None:
                    self.projRightRangeMax.setValue(self.shape[current.axes] - 1)
                else:
                    self.projRightRangeMax.setValue(current.projLimits[2])
                self.projRightRangeMax.valueChanged.connect(self.changeRanges)
                self.projRightRangeMax.hide()
                self.contourProjFrame.addWidget(self.projRightRangeMax, 5, 1)
                
                self.projRightRangeMinLabel = wc.QLeftLabel("Right min:", self)
                self.contourProjFrame.addWidget(self.projRightRangeMinLabel, 6, 0)
                self.projRightRangeMinLabel.hide()
                self.projRightRangeMin = QtWidgets.QSpinBox()
                self.projRightRangeMin.setMaximum(self.shape[current.axes] - 1)
                self.projRightRangeMin.setMinimum(0)
                if current.projLimits[3] is None:
                    self.projRightRangeMin.setValue(0)
                else:
                    self.projRightRangeMin.setValue(current.projLimits[3])
                self.projRightRangeMin.valueChanged.connect(self.changeRanges)
                self.projRightRangeMin.hide()
                self.contourProjFrame.addWidget(self.projRightRangeMin, 6, 1)
                self.contourProjGroup.setLayout(self.contourProjFrame)
                self.frame2.addWidget(self.contourProjGroup, 8, 0, 1, 3)
                self.activateRanges(self.rangeCheckbox.checkState())
                
                #Diagonal group
                self.diagonalGroup = QtWidgets.QGroupBox('Diagonal:')
                self.diagonalGroup.setCheckable(True)
                self.diagonalGroup.setChecked(current.diagonalBool)
                self.diagonalGroup.toggled.connect(self.switchDiagonal)
                self.diagonalFrame = QtWidgets.QGridLayout()
                self.diagMultiLabel = wc.QLeftLabel("Multiplier:", self)
                self.diagonalFrame.addWidget(self.diagMultiLabel, 0, 0)
                self.diagonalEntry = wc.QLineEdit(current.diagonalMult, self.setDiagonal)
                self.diagonalEntry.setMaximumWidth(120)
                self.diagonalFrame.addWidget(self.diagonalEntry, 0, 1)
                self.diagonalGroup.setLayout(self.diagonalFrame)
                self.frame2.addWidget(self.diagonalGroup, 9, 0, 1, 3)
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
                colorbutton = QtWidgets.QPushButton("Colour", self)
                colorbutton.clicked.connect(lambda arg, num=i: self.setExtraColor(num))
                frame.addWidget(colorbutton, 1, 0)
                button = QtWidgets.QPushButton("x", self)
                button.clicked.connect(lambda arg, num=i: self.delMultiSpec(num))
                frame.addWidget(button, 1, 1)
                self.OOM = self.father.current.getOOM()  # Order of Magnitude                
                self.scaleLabel = wc.QLeftLabel("Scale:", self)
                frame.addWidget(self.scaleLabel, 2, 0)
                self.offsetLabel = wc.QLeftLabel(u"Offset (\u00D71e" + str(self.OOM) + "):", self)
                frame.addWidget(self.offsetLabel, 3, 0)
                self.shiftLabel = wc.QLeftLabel("Shift:", self)
                frame.addWidget(self.shiftLabel, 4, 0)
                
                scaleEntry = QtWidgets.QDoubleSpinBox()
                scaleEntry.setMaximum(1e3)
                scaleEntry.setMinimum(-1e3)
                scaleEntry.setSingleStep(0.1)
                scaleEntry.setValue(self.father.current.extraScale[i])
                scaleEntry.valueChanged.connect(lambda arg, num=i: self.setScale(arg, num))
                frame.addWidget(scaleEntry, 2, 1)
                offsetEntry = QtWidgets.QDoubleSpinBox()
                offsetEntry.setMaximum(1e3)
                offsetEntry.setMinimum(-1e3)
                offsetEntry.setSingleStep(0.1)
                offsetEntry.setValue(self.father.current.extraOffset[i]/(10**self.OOM))
                offsetEntry.valueChanged.connect(lambda arg, num=i: self.setOffset(arg, num))
                frame.addWidget(offsetEntry, 3, 1)
                shiftEntry = QtWidgets.QDoubleSpinBox()
                shiftEntry.setMaximum(1e3)
                shiftEntry.setMinimum(-1e3)
                shiftEntry.setSingleStep(0.1)
                shiftEntry.setValue(self.father.current.extraShift[i])
                shiftEntry.valueChanged.connect(lambda arg, num=i: self.setShift(arg, num))
                frame.addWidget(shiftEntry, 4, 1)
                
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
                        frame.addWidget(buttons1[num], num * 3 + 6, 0)
                        frame.addWidget(wc.QLabel("D" + str(num + 1), self), num * 3 + 5, 1)
                        entries.append(wc.SliceSpinBox(self, 0, current.extraData[i].data.shape[num] - 1))
                        frame.addWidget(entries[num], num * 3 + 6, 1)
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
        var1 = self.numLEntry.value()
        maxC =safeEval(self.maxLEntry.text())
        if maxC is None:
            maxC = self.father.current.maxLevels * 100
            self.father.father.dispMsg('Invalid value for contour maximum')
        else:
            maxC = abs(float(maxC))
        minC =safeEval(self.minLEntry.text())
        if minC is None:
            minC = self.father.current.minLevels * 100
            self.father.father.dispMsg('Invalid value for contour minimum')
        else:
            minC = abs(float(minC))
        if minC > maxC: #if wrong order, interchange
            maxC, minC = (minC , maxC)
        self.maxLEntry.setText(str(maxC))
        self.minLEntry.setText(str(minC))
        cSign = self.contourSignEntry.currentIndex()
        cType =self.contourTypeEntry.currentIndex() 
        if cType == 0:
            self.multiValue.hide()
            self.multiValueLabel.hide()
        else:
            self.multiValue.show()
            self.multiValueLabel.show()
        multi = safeEval(self.multiValue.text())
        if multi is None:
            multi = self.father.current.multiValue
            self.father.father.dispMsg('Invalid value for contour multiplier')
        else:
            multi = abs(float(multi))
        self.multiValue.setText(str(multi))
        self.father.current.setLevels(var1, maxC / 100.0, minC / 100.0, cSign, cType, multi)

    def changeProj(self, pType, direc):
        self.father.current.setProjType(pType, direc)
        self.father.current.showProj()
    
    def changeRanges(self):
        Check = self.rangeCheckbox.isChecked()
        Ranges = [self.projTopRangeMax.value(),self.projTopRangeMin.value(),self.projRightRangeMax.value(),self.projRightRangeMin.value()]
        self.father.current.setProjLimits(Check, Ranges)
        self.father.current.showProj()
    
    def activateRanges(self,state):
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
                    axes2 = self.father.current.axes
                else:
                    axes = self.father.current.axes2
            self.buttons2Group.button(axes2).toggle()
        self.getSlice(None, axes, True)
        self.upd()
        self.father.menuCheck()

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
       
    def setShift(self, shift, num):
        self.father.current.setExtraShift(num, shift)    

    def switchDiagonal(self, val):
        self.father.current.setDiagonal(bool(val))

    def setDiagonal(self):
        inp = safeEval(self.diagonalEntry.text())
        if inp is None:
            inp = self.father.current.diagonalMult
            self.father.father.dispMsg('Invalid value for diagonal multiplier')
        else:
            inp = float(inp)
        self.diagonalEntry.setText(str(inp))
        self.father.current.setDiagonal(None, inp)

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
        super(BottomFrame, self).__init__(parent)
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
        self.freqEntry = wc.QLineEdit('', self.changeFreq, parent=self)
        grid.addWidget(self.freqEntry, 1, 3)
        grid.addWidget(wc.QLabel("Sweepwidth [kHz]:", self), 0, 4)
        self.swEntry = wc.QLineEdit('', self.changeFreq, parent=self)
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
        if self.father.masterData.noUndo:
            self.father.current.setWholeEcho(inp)
        else:
            self.father.undoList.append(self.father.current.setWholeEcho(inp))
        self.father.menuCheck()

    def changeSpec(self):
        self.father.redoList = []
        if self.father.masterData.noUndo:
            self.father.current.changeSpec(self.specGroup.checkedId())
        else:
            self.father.undoList.append(self.father.current.changeSpec(self.specGroup.checkedId()))
        self.upd()
        self.father.menuCheck()

    def changeFreq(self):
        freq = safeEval(self.freqEntry.text())
        sw = safeEval(self.swEntry.text())
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
        self.father.current.plotType = pType
        self.father.current.showFid()

    def changeAxis(self, pType):
        self.father.current.setAxType(pType)
        self.father.current.showFid()

    def changeAxis2(self, pType):
        self.father.current.setAxType2(pType)
        self.father.current.showFid()

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
        getButton.clicked.connect(self.getPosition)
        grid.addWidget(getButton, 0, 1)
        grid.addWidget(wc.QLabel("x-Position:"), 0, 2)
        self.xpos = wc.QLineEdit("0")
        self.xpos.setFixedWidth(self.xpos.sizeHint().width() * widthScale)
        grid.addWidget(self.xpos, 0, 3)
        self.yposlabel = wc.QLabel("y-Position:")
        grid.addWidget(self.yposlabel, 0, 4)
        self.ypos = wc.QLineEdit("0")
        self.ypos.setFixedWidth(self.ypos.sizeHint().width() * widthScale)
        grid.addWidget(self.ypos, 0, 5)
        grid.addWidget(wc.QLabel("x-Value:"), 0, 6)
        self.xpoint = wc.QLineEdit("0.0")
        self.xpoint.setFixedWidth(self.xpoint.sizeHint().width() * widthScale)
        grid.addWidget(self.xpoint, 0, 7)
        self.ylabel = wc.QLabel("y-Value:")
        grid.addWidget(self.ylabel, 0, 8)
        self.ypoint = wc.QLineEdit("0.0")
        self.ypoint.setFixedWidth(self.ypoint.sizeHint().width() * widthScale)
        grid.addWidget(self.ypoint, 0, 9)
        grid.addWidget(wc.QLabel("Amp:"), 0, 10)
        self.amppoint = wc.QLineEdit("0.0")
        self.amppoint.setFixedWidth(self.amppoint.sizeHint().width() * widthScale)
        grid.addWidget(self.amppoint, 0, 11)
        grid.addWidget(wc.QLabel(u"\u0394x:"), 0, 12)
        self.deltaxpoint = wc.QLineEdit("0.0")
        self.deltaxpoint.setFixedWidth(self.deltaxpoint.sizeHint().width() * widthScale)
        grid.addWidget(self.deltaxpoint, 0, 13)
        self.deltaylabel = wc.QLabel(u"\u0394y:")
        grid.addWidget(self.deltaylabel, 0, 14)
        self.deltaypoint = wc.QLineEdit("0.0")
        self.deltaypoint.setFixedWidth(self.deltaypoint.sizeHint().width() * widthScale)
        grid.addWidget(self.deltaypoint, 0, 15)
        grid.addWidget(wc.QLabel(u"\u0394amp:"), 0, 16)
        self.deltaamppoint = wc.QLineEdit("0.0")
        self.deltaamppoint.setFixedWidth(self.deltaamppoint.sizeHint().width() * widthScale)
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
            self.oldy = position[4]
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


class AsciiLoadWindow(QtWidgets.QDialog):
    
    dataOrders = ['XRI','XR','XI','RI','R']
    delimiters = ['Tab','Space','Comma']
    
    def __init__(self, parent, file):
        super(AsciiLoadWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool | QtCore.Qt.WindowContextHelpButtonHint  )
        self.father = parent
        self.dataDimension = 1
        self.dataSpec = False
        self.dataOrder = 'XRI'
        self.sw = 0.0
        self.delim = 'Tab'
        self.closed = False
        self.setWindowTitle("Load ASCII")
        grid = QtWidgets.QGridLayout(self)
        
        grid.addWidget(QtWidgets.QLabel("# Dimensions:"), 1, 0)
        self.numDims = QtWidgets.QSpinBox()
        self.numDims.setMinimum(1)
        self.numDims.setValue(1)
        self.numDims.setMaximum(2)
        grid.addWidget(self.numDims, 2, 0, 1, 2)
        
        grid.addWidget(QtWidgets.QLabel("Data Type:"), 3, 0)
        
        self.specGroup = QtWidgets.QButtonGroup(self)
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
        
        self.swLabel = wc.QLeftLabel("Spectral Width [kHz]:")
        grid.addWidget(self.swLabel, 7, 0, 1, 2)
        self.swEntry = wc.QLineEdit("0.0")
        grid.addWidget(self.swEntry, 8, 0, 1, 2)
        self.swLabel.hide()
        self.swEntry.hide()
        self.datOrderBox.currentIndexChanged.connect(self.checkDatOrder)
        
        grid.addWidget(QtWidgets.QLabel("Data Delimiter:"), 9, 0)
        self.datDelimBox = QtWidgets.QComboBox()
        self.datDelimBox.addItems(self.delimiters)
        grid.addWidget(self.datDelimBox, 10, 0, 1, 2)
        
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        grid.addWidget(cancelButton, 13, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        grid.addWidget(okButton, 13, 1)
        
        self.show()
        self.setFixedSize(self.size())
        self.checkType(file)
        
    def checkDatOrder(self):
        tmp = self.dataOrders[ self.datOrderBox.currentIndex() ]
        if tmp == 'RI' or tmp == 'R':
            self.swLabel.show()
            self.swEntry.show()
        else:
            self.swLabel.hide()
            self.swEntry.hide()
        
    def checkType(self, file):
        try:
            with open(file, 'r') as f:
                line = f.readline()
            if line.count(',') > 0:
                sep = 'Comma'
            elif line.count('\t') > 0:
                sep = 'Tab'
            else:
                sep = 'Space'
            self.datDelimBox.setCurrentIndex(self.delimiters.index(sep)) 
            sepList = ['\t',' ',',']           
            data = np.fromstring(line, sep = sepList[self.delimiters.index(sep)])
            if len(data) > 3:
               self.numDims.setValue(2) 
#            if len(data) % 2 == 0: #if odd size
#                self.datOrderBox.setCurrentIndex(1)
        except:
            return
    
    def closeEvent(self, *args):
        self.closed = True
        self.accept()
        self.deleteLater()

    def applyAndClose(self):
        self.dataOrder = self.dataOrders[ self.datOrderBox.currentIndex() ]
        self.delim = self.delimiters[ self.datDelimBox.currentIndex() ]
        if self.dataOrder == 'RI' or self.dataOrder == 'R':
            self.sw = safeEval(self.swEntry.text())
            if self.sw == 0 or self.sw == None:
                self.father.dispMsg('Spectral Width input is not valid')
                return
            
        self.dataDimension = self.numDims.value()
        if self.timeButton.isChecked():
           self.dataSpec = False
        else:
           self.dataSpec = True

        self.accept()
        self.deleteLater()

#################################################################################


class PhaseWindow(wc.ToolWindows):

    NAME = "Phasing"
    SINGLESLICE = True
    RESOLUTION = 1000.0
    P1LIMIT = 540.0
    PHASE0STEP = 1.0
    PHASE1STEP = 1.0

    def __init__(self, parent):
        super(PhaseWindow, self).__init__(parent)
        self.zeroVal = 0.0
        self.firstVal = 0.0
        self.refVal = 0.0
        self.available = True
        #Zero order
        self.zeroOrderGroup = QtWidgets.QGroupBox('Zero order:')
        self.zeroOrderFrame = QtWidgets.QGridLayout()
        autoZero = QtWidgets.QPushButton("Autophase 0th")
        autoZero.clicked.connect(lambda: self.autophase(0))
        self.zeroOrderFrame.addWidget(autoZero, 0, 1)
        self.zeroEntry = wc.QLineEdit("0.000", self.inputZeroOrder)
        self.zeroOrderFrame.addWidget(self.zeroEntry, 2, 1)
        leftZero = QtWidgets.QPushButton("<")
        leftZero.clicked.connect(lambda: self.stepPhase(-1, 0))
        leftZero.setAutoRepeat(True)
        self.zeroOrderFrame.addWidget(leftZero, 2, 0)
        rightZero = QtWidgets.QPushButton(">")
        rightZero.clicked.connect(lambda: self.stepPhase(1, 0))
        rightZero.setAutoRepeat(True)
        self.zeroOrderFrame.addWidget(rightZero, 2, 2)
        self.zeroScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.zeroScale.setRange(-self.RESOLUTION, self.RESOLUTION)
        self.zeroScale.valueChanged.connect(self.setZeroOrder)
        self.zeroOrderFrame.addWidget(self.zeroScale, 3, 0, 1, 3)
        self.zeroOrderGroup.setLayout(self.zeroOrderFrame)
        self.grid.addWidget(self.zeroOrderGroup,0,0,1,3)
        #First order
        self.firstOrderGroup = QtWidgets.QGroupBox('First order:')
        self.firstOrderFrame = QtWidgets.QGridLayout()
        autoFirst = QtWidgets.QPushButton("Autophase 0th+1st")
        autoFirst.clicked.connect(lambda: self.autophase(1))
        self.firstOrderFrame.addWidget(autoFirst, 5, 1)
        self.firstEntry = wc.QLineEdit("0.000", self.inputFirstOrder)
        self.firstOrderFrame.addWidget(self.firstEntry, 6, 1)
        leftFirst = QtWidgets.QPushButton("<")
        leftFirst.clicked.connect(lambda: self.stepPhase(0, -1))
        leftFirst.setAutoRepeat(True)
        self.firstOrderFrame.addWidget(leftFirst, 6, 0)
        rightFirst = QtWidgets.QPushButton(">")
        rightFirst.clicked.connect(lambda: self.stepPhase(0, 1))
        rightFirst.setAutoRepeat(True)
        self.firstOrderFrame.addWidget(rightFirst, 6, 2)
        self.firstScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.firstScale.setRange(-self.RESOLUTION, self.RESOLUTION)
        self.firstScale.valueChanged.connect(self.setFirstOrder)
        self.firstOrderFrame.addWidget(self.firstScale, 7, 0, 1, 3)
        if self.father.current.spec > 0:
            self.firstOrderFrame.addWidget(wc.QLabel("Reference:"), 8, 0, 1, 3)
            pickRef = QtWidgets.QPushButton("Pick reference")
            pickRef.clicked.connect(self.pickRef)
            self.firstOrderFrame.addWidget(pickRef, 9, 1)
            self.refEntry = wc.QLineEdit(('%.3f' % self.refVal), self.inputRef)
            self.firstOrderFrame.addWidget(self.refEntry, 10, 1)
        self.firstOrderGroup.setLayout(self.firstOrderFrame)
        self.grid.addWidget(self.firstOrderGroup,1,0,1,3)

    def setZeroOrder(self, value, *args):
        if self.available:
            self.zeroVal = float(value) / self.RESOLUTION * 180
            self.zeroEntry.setText('%.3f' % self.zeroVal)
            self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0)

    def inputZeroOrder(self, *args):
        inp = safeEval(self.zeroEntry.text())
        if inp is None:
            self.father.father.dispMsg('Phasing: zero order value input is not valid!')
            return None
        self.zeroVal = np.mod(inp + 180, 360) - 180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.available = False
        self.zeroScale.setValue(round(self.zeroVal / 180.0 * self.RESOLUTION))
        self.available = True
        self.father.current.setPhaseInter(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0)
        return 1

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
        value = safeEval(self.firstEntry.text())
        if value == None:
            self.father.father.dispMsg('Phasing: first order value input is not valid!')
            return None
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
        return 1

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
        inp = safeEval(self.zeroEntry.text())
        if inp == None:
            self.father.father.dispMsg('Phasing: zero order value input is not valid!')
            return None
        inp +=  phase0 * self.PHASE0STEP
        self.zeroVal = np.mod(inp + 180, 360) - 180
        value = safeEval(self.firstEntry.text())
        if value == None:
            self.father.father.dispMsg('Phasing: first order value input is not valid!')
            return None
        value += phase1 * self.PHASE1STEP

        if self.father.current.spec > 0:
            refCheck = self.inputRef()
            if refCheck == None:
                return

        value += phase1 * self.PHASE1STEP
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
        Val = safeEval(self.refEntry.text())
        if Val == None:
            self.father.father.dispMsg('Phasing: reference input is not valid!')
            return None
        self.refVal = Val
        self.refEntry.setText('%.3f' % self.refVal)
        return 1

    def setRef(self, value, *args):
        self.refVal = float(value)
        self.refEntry.setText('%.3f' % self.refVal)

    def pickRef(self, *args):
        self.father.current.peakPickFunc = lambda pos, self=self: self.setRef(self.father.current.xax[pos[0]])
        self.father.current.peakPick = True

    def applyFunc(self):
        refCheck = 1
        if self.father.current.spec > 0:
            refCheck = self.inputRef()
        zeroCheck = self.inputZeroOrder()
        if zeroCheck is not None:
            firstCheck = self.inputFirstOrder()
        if refCheck == None or zeroCheck == None or firstCheck == None: #If error. Messages are handled by functions
            return False
        self.father.redoList = []
        if self.father.current.data.noUndo:
            self.father.current.applyPhase(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, (self.singleSlice.isChecked() == 1))
        else:
            self.father.undoList.append(self.father.current.applyPhase(np.pi * self.zeroVal / 180.0, np.pi * self.firstVal / 180.0, (self.singleSlice.isChecked() == 1)))

################################################################


class ApodWindow(wc.ToolWindows):

    RESOLUTION = 10000
    NAME = "Apodize"
    SINGLESLICE = True

    def __init__(self, parent):
        super(ApodWindow, self).__init__(parent)
        self.entries = []
        self.ticks = []
        self.maximum = 100.0 * self.father.current.sw / (self.father.current.data1D.shape[-1])
        self.lorstep = 1.0
        self.gaussstep = 1.0
        self.available = True
        lorTick = QtWidgets.QCheckBox("Lorentzian:")
        lorTick.toggled.connect(lambda: self.checkEval(0))
        self.grid.addWidget(lorTick, 0, 0, 1, 3)
        self.ticks.append(lorTick)
        lorEntry = wc.QLineEdit("0.00", self.apodPreview)
        lorEntry.setMinimumWidth(150)
        lorEntry.setEnabled(False)
        self.grid.addWidget(lorEntry, 1, 1)
        self.entries.append(lorEntry)
        leftLor = QtWidgets.QPushButton("<")
        leftLor.clicked.connect(lambda: self.stepLB(-0.5 * self.father.current.sw / (self.father.current.data1D.shape[-1]), 0))
        leftLor.setAutoRepeat(True)
        self.grid.addWidget(leftLor, 1, 0)
        rightLor = QtWidgets.QPushButton(">")
        rightLor.clicked.connect(lambda: self.stepLB(0.5 * self.father.current.sw / (self.father.current.data1D.shape[-1]), 0))
        rightLor.setAutoRepeat(True)
        self.grid.addWidget(rightLor, 1, 2)
        self.lorScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.lorScale.setRange(0, self.RESOLUTION)
        self.lorScale.valueChanged.connect(self.setLor)
        self.grid.addWidget(self.lorScale, 2, 0, 1, 3)
        self.lorMax = 100.0 * self.father.current.sw / (self.father.current.data1D.shape[-1])

        gaussTick = QtWidgets.QCheckBox("Gaussian:")
        gaussTick.toggled.connect(lambda: self.checkEval(1))
        self.grid.addWidget(gaussTick, 3, 0, 1, 3)
        self.ticks.append(gaussTick)
        gaussEntry = wc.QLineEdit("0.00", self.apodPreview)
        gaussEntry.setEnabled(False)
        gaussEntry.setMinimumWidth(150)
        self.grid.addWidget(gaussEntry, 4, 1)
        self.entries.append(gaussEntry)
        leftGauss = QtWidgets.QPushButton("<")
        leftGauss.clicked.connect(lambda: self.stepLB(0, -0.5 * self.father.current.sw / (self.father.current.data1D.shape[-1])))
        leftGauss.setAutoRepeat(True)
        self.grid.addWidget(leftGauss, 4, 0)
        rightGauss = QtWidgets.QPushButton(">")
        rightGauss.clicked.connect(lambda: self.stepLB(0, 0.5 * self.father.current.sw / (self.father.current.data1D.shape[-1])))
        rightGauss.setAutoRepeat(True)
        self.grid.addWidget(rightGauss, 4, 2)
        self.gaussScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.gaussScale.setRange(0, self.RESOLUTION)
        self.gaussScale.valueChanged.connect(self.setGauss)
        self.grid.addWidget(self.gaussScale, 5, 0, 1, 3)
        self.gaussMax = 100.0 * self.father.current.sw / (self.father.current.data1D.shape[-1])

        cos2Tick = QtWidgets.QCheckBox("Cos^2:")
        cos2Tick.clicked.connect(lambda: self.checkEval(2))
        self.grid.addWidget(cos2Tick, 6, 0, 1, 3)
        self.ticks.append(cos2Tick)
        cos2Entry = wc.QLineEdit("1.00", self.apodPreview)
        cos2Entry.setEnabled(False)
        self.grid.addWidget(cos2Entry, 7, 1)
        self.entries.append(cos2Entry)

        hammingTick = QtWidgets.QCheckBox("Hamming:")
        hammingTick.clicked.connect(lambda: self.checkEval(3))
        self.grid.addWidget(hammingTick, 8, 0, 1, 3)
        self.ticks.append(hammingTick)
        hammingEntry = wc.QLineEdit("1.00", self.apodPreview)
        hammingEntry.setEnabled(False)
        self.grid.addWidget(hammingEntry, 9, 1)
        self.entries.append(hammingEntry)

        self.grid.addWidget(wc.QLabel("Shift:"), 10, 0, 1, 3)
        self.shiftEntry = wc.QLineEdit("0.00", self.apodPreview)
        self.grid.addWidget(self.shiftEntry, 11, 1)

        if self.father.current.data.data.ndim > 1:
            self.grid.addWidget(wc.QLabel("Shifting:"), 12, 0, 1, 3)
            
            self.shiftingDropdown = QtWidgets.QComboBox()
            self.shiftingDropdown.addItems(['User Defined','Spin 3/2, -3Q (7/9)','Spin 5/2, 3Q (19/12)','Spin 5/2, -5Q (25/12)','Spin 7/2, 3Q (101/45)',
                                         'Spin 7/2, 5Q (11/9)','Spin 7/2, -7Q (161/45)','Spin 9/2, 3Q (91/36)','Spin 9/2, 5Q (95/36)','Spin 9/2, 7Q (7/18)','Spin 9/2, -9Q (31/6)'])
            self.shiftingDropdown.activated.connect(self.dropdownChanged)
            self.shiftingList = [0,7.0/9.0,19.0/12.0,25.0/12.0,101.0/45.0,11.0/9.0,161.0/45.0,91.0/36.0,95.0/36.0,7.0/18.0,31.0/6.0]
        
            self.grid.addWidget(self.shiftingDropdown, 13, 1)
            
            self.shiftingEntry = wc.QLineEdit("0.00", self.apodPreview)
            self.grid.addWidget(self.shiftingEntry, 14, 1)
            self.shiftingAxes = QtWidgets.QComboBox()
            self.shiftingValues = list(map(str, np.delete(range(1, self.father.current.data.data.ndim + 1), self.father.current.axes)))
            self.shiftingAxes.addItems(self.shiftingValues)
            self.shiftingAxes.currentIndexChanged.connect(self.apodPreview)
            self.grid.addWidget(self.shiftingAxes, 15, 1)

    def dropdownChanged(self):
        index =  self.shiftingDropdown.currentIndex()
        self.shiftingEntry.setText("%.9f" % self.shiftingList[index])
        self.apodPreview()
        
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
                self.father.father.dispMsg('Apodize: Lorentzian value is not valid!')
                self.father.current.showFid()
                return False
            self.entries[0].setText('%.4g' % lor)
            self.lorScale.setValue(round(lor * self.RESOLUTION / self.maximum))
        if self.ticks[1].isChecked():
            gauss = safeEval(self.entries[1].text())
            if gauss is None:
                self.father.father.dispMsg('Apodize: Gaussian value is not valid!')
                self.father.current.showFid()
                return False
            self.entries[1].setText('%.4g' % gauss)
            self.gaussScale.setValue(round(gauss * self.RESOLUTION / self.maximum))
        if self.ticks[2].isChecked():
            cos2 = safeEval(self.entries[2].text())
            if cos2 is None:
                self.father.father.dispMsg('Apodize: cos^2 value is not valid!')
                self.father.current.showFid()
                return False
            self.entries[2].setText('%.4g' % cos2)
        if self.ticks[3].isChecked():
            hamming = safeEval(self.entries[3].text())
            if hamming is None:
                self.father.father.dispMsg('Apodize: Hamming value is not valid!')
                self.father.current.showFid()
                return False
            self.entries[3].setText('%.4g' % hamming)
        shift = safeEval(self.shiftEntry.text())
        if shift is None:
            self.father.father.dispMsg('Apodize: Shift value is not valid!')
            self.father.current.showFid()
            return False
        self.shiftEntry.setText('%.4g' % shift)
        if self.father.current.data.data.ndim > 1:
            shifting = safeEval(self.shiftingEntry.text())
            if shifting is None:
                self.father.father.dispMsg('Apodize: Shifting value is not valid!')
                self.father.current.showFid()
                return False
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
        lor = safeEval(self.entries[0].text())
        if lor is None:
            self.father.father.dispMsg('Apodize: Lorentzian value is not valid!')
            self.father.current.showFid()
            return False
        self.entries[0].setText('%.4g' % (lor + multiplier*lorincr * self.lorstep))
        gauss = safeEval(self.entries[1].text())
        if gauss is None:
            self.father.father.dispMsg('Apodize: Gaussian value is not valid!')
            self.father.current.showFid()
            return False
        self.entries[1].setText('%.4g' % (gauss + multiplier*gaussincr * self.gaussstep))
        if (lorincr != 0) and (not self.ticks[0].isChecked()):
            self.ticks[0].setChecked(1)
        if (gaussincr != 0) and (not self.ticks[1].isChecked()):
            self.ticks[1].setChecked(1)
        self.apodPreview()

    def applyFunc(self):
        lor = None
        gauss = None
        cos2 = None
        hamming = None
        shifting = None
        shiftingAxes = 0
        if self.ticks[0].isChecked():
            lor = safeEval(self.entries[0].text())
            if lor is None:
                self.father.father.dispMsg('Apodize: Lorentzian value is not valid!')
                self.father.current.showFid()
                return False
        if self.ticks[1].isChecked():
            gauss = safeEval(self.entries[1].text())
            if gauss is None:
                self.father.father.dispMsg('Apodize: Gaussian value is not valid!')
                self.father.current.showFid()
                return False
        if self.ticks[2].isChecked():
            cos2 = safeEval(self.entries[2].text())
            if cos2 is None:
                self.father.father.dispMsg('Apodize: cos^2 value is not valid!')
                self.father.current.showFid()
                return False
        if self.ticks[3].isChecked():
            hamming = safeEval(self.entries[3].text())
            if hamming is None:
                self.father.father.dispMsg('Apodize: Hamming value is not valid!')
                self.father.current.showFid()
                return False
        shift = safeEval(self.shiftEntry.text())
        if shift is None:
            self.father.father.dispMsg('Apodize: Shift value is not valid!')
            self.father.current.showFid()
            return False
        if self.father.current.data.data.ndim > 1:
            shifting = safeEval(self.shiftingEntry.text())
            if shifting is None:
                self.father.father.dispMsg('Apodize: Shifting value is not valid!')
                self.father.current.showFid()
                return False
            shiftingAxes = int(self.shiftingValues[self.shiftingAxes.currentIndex()]) - 1
        else:
            shiftingAxes = None
        self.father.redoList = []
        if self.father.current.data.noUndo:
            self.father.current.applyApod(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, (self.singleSlice.isChecked()))
        else:
            self.father.undoList.append(self.father.current.applyApod(lor, gauss, cos2, hamming, shift, shifting, shiftingAxes, (self.singleSlice.isChecked())))

#######################################################################################


class SizeWindow(wc.ToolWindows):

    NAME = "Set size"

    def __init__(self, parent):
        super(SizeWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Size:"), 0, 0)
        self.sizeVal = parent.current.data1D.shape[-1]
        self.sizeEntry = wc.QLineEdit(self.sizeVal, self.sizePreview)
        self.grid.addWidget(self.sizeEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("Offset:"), 2, 0)
        if self.father.current.wholeEcho:
            self.posVal = int(np.floor(parent.current.data1D.shape[-1] / 2.0))
        else:
            self.posVal = parent.current.data1D.shape[-1]
        self.posEntry = wc.QLineEdit(self.posVal, self.sizePreview)
        self.grid.addWidget(self.posEntry, 3, 0)
        if not self.father.current.spec:
            self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
            self.father.current.peakPick = True

    def sizePreview(self, *args):
        inp = safeEval(self.sizeEntry.text())
        if inp is not None:
            self.sizeVal = int(round(inp))
        if self.sizeVal < 1 or inp is None:
            self.father.father.dispMsg('Sizing: \'Size\' input is not valid')
            return False
        self.sizeEntry.setText(str(self.sizeVal))
        inp = safeEval(self.posEntry.text())
        if inp is not None:
            self.posVal = int(round(inp))
        if self.posVal < 1 or inp is None:
            self.father.father.dispMsg('Sizing: \'Offset\' input is not valid')
            return False
        self.posEntry.setText(str(self.posVal))
        self.father.current.setSizePreview(self.sizeVal, self.posVal)

    def applyFunc(self):
        inp = safeEval(self.sizeEntry.text())
        if inp is not None:
            self.sizeVal = int(round(inp))
        if self.sizeVal < 1 or inp is None:
            self.father.father.dispMsg('Sizing: \'Size\' input is not valid')
            return False
        inp = safeEval(self.posEntry.text())
        if inp is not None:
            self.posVal = int(round(inp))
        if self.posVal < 1 or inp is None:
            self.father.father.dispMsg('Sizing: \'Offset\' input is not valid')
            return False
        self.father.redoList = []
        if self.father.current.data.noUndo:
            self.father.current.applySize(self.sizeVal, self.posVal)
        else:
            self.father.undoList.append(self.father.current.applySize(self.sizeVal, self.posVal))
        self.father.sideframe.upd()

    def picked(self, pos):
        self.posEntry.setText(str(pos[0]))
        self.sizePreview()
        self.father.current.peakPick = True
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)

##########################################################################################


class SwapEchoWindow(wc.ToolWindows):

    NAME = "Swap echo"
    
    def __init__(self, parent):
        super(SwapEchoWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Echo position:"), 0, 0)
        self.posVal = int(round(0.5 * len(parent.current.data1D)))
        self.posEntry = wc.QLineEdit(self.posVal, self.swapEchoPreview)
        self.grid.addWidget(self.posEntry, 1, 0)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def swapEchoPreview(self, *args):
        inp = safeEval(self.posEntry.text())
        if inp is not None:
            self.posVal = int(round(inp))
        else:
            self.father.father.dispMsg("Swap echo: not a valid index")
            return False
        self.posEntry.setText(str(self.posVal))
        if self.posVal > 0 and self.posVal < (self.father.current.data1D.shape[-1]):
            self.father.current.setSwapEchoPreview(self.posVal)
            self.father.current.peakPick = False

    def applyFunc(self):
        self.father.current.peakPickReset()
        inp = safeEval(self.posEntry.text())
        if inp is not None:
            self.posVal = int(round(inp))
        else:
            self.father.father.dispMsg("Swap echo: not a valid index")
            return False
        self.posEntry.setText(str(self.posVal))
        if self.posVal > 0 and self.posVal < (self.father.current.data1D.shape[-1]):
            self.father.redoList = []
            if self.father.current.data.noUndo:
                self.father.current.applySwapEcho(self.posVal)
            else:
                self.father.undoList.append(self.father.current.applySwapEcho(self.posVal))
            self.father.bottomframe.upd()
        else:
            self.father.father.dispMsg("Swap echo: not a valid index")

    def picked(self, pos):
        self.father.current.setSwapEchoPreview(pos[0])
        self.posEntry.setText(str(pos[0]))
        self.father.current.peakPick = False

###########################################################################


class LPSVDWindow(wc.ToolWindows):

    NAME = "LPSVD"

    def __init__(self, parent):
        super(LPSVDWindow, self).__init__(parent)
        # self.grid.addWidget(wc.QLabel("# points for analysis:"), 2, 0)
        self.specGroup = QtWidgets.QButtonGroup(self)
        # self.specGroup.buttonClicked.connect(self.changeSpec)
        backwardButton = QtWidgets.QRadioButton('Backward', parent=self)
        self.specGroup.addButton(backwardButton, 1)
        forwardButton = QtWidgets.QRadioButton('Forward', parent=self)
        self.specGroup.addButton(forwardButton, 0)
        self.grid.addWidget(backwardButton, 1, 0)
        self.grid.addWidget(forwardButton, 2, 0)
        backwardButton.setChecked(True)

        self.grid.addWidget(wc.QLabel("# points for analysis:"), 3, 0)
        self.analPoints = 200
        self.aPointsEntry = wc.QLineEdit(self.analPoints)
        self.grid.addWidget(self.aPointsEntry, 4, 0)
        self.grid.addWidget(wc.QLabel("Number of frequencies:"), 5, 0)
        self.numberFreq = 1
        self.nFreqEntry = wc.QLineEdit(self.numberFreq)
        self.grid.addWidget(self.nFreqEntry, 6, 0)

        self.grid.addWidget(wc.QLabel("Number prediction points:"), 7, 0)
        self.predictPoints = 10
        self.nPredictEntry = wc.QLineEdit(self.predictPoints)
        self.grid.addWidget(self.nPredictEntry, 8, 0)

    def applyFunc(self):
        analPoints = safeEval(self.aPointsEntry.text())
        if analPoints is None:
            self.father.father.dispMsg('LPSVD: Number of points for analysis is not valid')
            return False
        numberFreq = safeEval(self.nFreqEntry.text())
        if numberFreq is None:
            self.father.father.dispMsg('LPSVD: Number of frequencies is not valid')
            return False
        predictPoints = safeEval(self.nPredictEntry.text())
        if predictPoints is None:
            self.father.father.dispMsg('LPSVD: Number of predication points is not valid')
            return False
        if self.analPoints > len(self.father.current.data1D):
            self.father.father.dispMsg('LPSVD: number of points for analysis cannot be more than data size')
            return False
        if self.analPoints <= self.numberFreq * 4:
            self.father.father.dispMsg('LPSVD: number of points for analysis must be more than 4 times the number of frequencies')
            return False
        self.father.redoList = []
        if self.father.current.data.noUndo:
            self.father.current.applyLPSVD(self.analPoints, self.numberFreq, self.predictPoints, self.specGroup.checkedId())
        else:
            self.father.undoList.append(self.father.current.applyLPSVD(self.analPoints, self.numberFreq, self.predictPoints, self.specGroup.checkedId()))
        self.father.sideframe.upd()

###########################################################################


class ShiftDataWindow(wc.ToolWindows):

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
        inp = safeEval(self.shiftEntry.text())
        if inp is not None:
            self.shiftVal = int(round(inp))
        else:
            self.father.father.dispMsg("Shift data: shift value not valid")
            return 
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
        else:
            self.father.father.dispMsg("Shift data: shift value not valid")
            return 
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
        else:
            self.father.father.dispMsg("Shift data: shift value not valid")
            return 
        self.shiftEntry.setText(str(self.shiftVal))
        self.father.current.setShiftPreview(self.shiftVal)

    def applyFunc(self):
        inp = safeEval(self.shiftEntry.text())
        if inp is None:
            self.father.father.dispMsg("Shift data: shift value not valid")
            return False
        shift = int(round(inp))
        self.father.redoList = []
        if self.father.current.data.noUndo:
            self.father.current.applyShift(shift, (self.singleSlice.isChecked()))
        else:
            self.father.undoList.append(self.father.current.applyShift(shift, (self.singleSlice.isChecked())))

#############################################################


class DCWindow(wc.ToolWindows):

    NAME = "Offset correction"
    SINGLESLICE = True

    def __init__(self, parent):
        super(DCWindow, self).__init__(parent)
        self.startVal = int(round(0.8 * parent.current.data1D.shape[-1]))
        self.endVal = parent.current.data1D.shape[-1]
        self.grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.startEntry = wc.QLineEdit(self.startVal, self.offsetPreview)
        self.grid.addWidget(self.startEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End point:"), 2, 0)
        self.endEntry = wc.QLineEdit(self.endVal, self.offsetPreview)
        self.grid.addWidget(self.endEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Offset:"), 4, 0)
        val = parent.current.getdcOffset(int(round(0.8 * parent.current.data1D.shape[-1])), parent.current.data1D.shape[-1])
        self.offsetEntry = wc.QLineEdit('{:.2e}'.format(val), lambda: self.offsetPreview(True))
        self.grid.addWidget(self.offsetEntry, 5, 0)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

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
            inp = safeEval(self.endEntry.text())
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
            self.father.current.peakPickFunc = lambda pos, self= self: self.picked(pos, True)
            self.father.current.peakPick = True

    def offsetPreview(self, inserted=False):
        if inserted:
            dcVal = safeEval(self.offsetEntry.text())
            if dcVal is None:
                self.father.father.dispMsg("Offset correction: offset value not valid")
                return
            self.father.current.dcOffset(dcVal)
        else:
            dataLength = self.father.current.data1D.shape[-1]
            inp = safeEval(self.startEntry.text())
            if inp is not None:
                self.startVal = int(round(inp))
            else:
                self.father.father.dispMsg("Offset correction: start value not valid")
                return
            if self.startVal < 0:
                self.startVal = 0
            elif self.startVal > dataLength:
                self.startVal = dataLength
            self.startEntry.setText(str(self.startVal))
            inp = safeEval(self.endEntry.text())
            if inp is not None:
                self.endVal = int(round(inp))
            else:
                self.father.father.dispMsg("Offset correction: end value not valid")
                return
            if self.endVal < 0:
                self.endVal = 0
            elif self.endVal > dataLength:
                self.endVal = dataLength
            self.endEntry.setText(str(self.endVal))
            val = self.father.current.getdcOffset(self.startVal, self.endVal)
            self.offsetEntry.setText('{:.2e}'.format(val))
            self.father.current.dcOffset(val)

    def applyFunc(self):
        inp = safeEval(self.offsetEntry.text())
        if inp is None:
            self.father.father.dispMsg("Offset correction: offset value not valid")
            return False
        self.father.current.peakPickReset()
        self.father.redoList = []
        if self.father.current.data.noUndo:
            self.father.current.subtract(inp, self.singleSlice.isChecked())
        else:
            self.father.undoList.append(self.father.current.subtract(inp, self.singleSlice.isChecked()))

#############################################################


class BaselineWindow(wc.ToolWindows):

    NAME = "Baseline correction"
    SINGLESLICE = True

    def __init__(self, parent):
        super(BaselineWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Polynomial Degree:"), 0, 0, 1, 2)
        self.removeList = []
        self.degreeEntry = QtWidgets.QSpinBox()
        self.degreeEntry.setMaximum(100)
        self.degreeEntry.setMinimum(1)
        self.degreeEntry.setValue(3)
        self.degreeEntry.setAlignment(QtCore.Qt.AlignCenter)
        self.grid.addWidget(self.degreeEntry, 1, 0, 1, 2)
        resetButton = QtWidgets.QPushButton("&Reset")
        resetButton.clicked.connect(self.reset)
        self.grid.addWidget(resetButton, 2, 0)
        fitButton = QtWidgets.QPushButton("&Fit")
        fitButton.clicked.connect(self.preview)
        self.grid.addWidget(fitButton, 2, 1)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def picked(self, pos):
        self.removeList.append(pos[0])
        self.father.current.previewRemoveList(self.removeList)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def preview(self, *args):
        inp = self.degreeEntry.value()
        check = self.father.current.previewBaseline(inp, self.removeList)
        if check == False:
            self.father.father.dispMsg("Baseline correct: error in polynomial fit",'red')
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

    def reset(self, *args):
        self.removeList = []
        self.father.current.resetPreviewRemoveList()
        self.preview()

    def closeEvent(self, *args):
        self.father.current.removeListLines = []
        del self.father.current.removeListLines
        super(BaselineWindow,self).closeEvent(*args)

    def applyFunc(self):
        inp = self.degreeEntry.value()
        returnValue = self.father.current.applyBaseline(inp, self.removeList, self.singleSlice.isChecked())
        if returnValue is None:
            self.father.father.dispMsg("Baseline correct: error in polynomial fit",'red')
            return False
        if not self.father.current.data.noUndo:
            self.father.undoList.append(returnValue)
        self.father.current.peakPickReset()
        self.father.current.resetPreviewRemoveList()
        self.father.redoList = []

#############################################################


class regionWindow(wc.ToolWindows):

    def __init__(self, parent, name):
        self.NAME = name
        super(regionWindow, self).__init__(parent)
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
            self.endVal = np.append(self.endVal, self.father.current.data1D.shape[-1])
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
        error = False
        if inp is not None:
            inp = int(inp)
            if inp < 0:
                inp = 0
            if inp > self.father.current.data1D.shape[-1]:
                inp = self.father.current.data1D.shape[-1]
        if isMin:
            num = self.startEntry.index(entry)
            if inp is None:
                self.startVal[num] = -1 #If the input is wrong, use -1 as a placeholder for it in the value list
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
            self.endVal = np.append(self.endVal, self.father.current.data1D.shape[-1])
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

        if error: #Return only after partIter check
            return

        if self.startVal[num] != -1: #Only if the input is OK, reprint it
            self.startEntry[num].setText(str(self.startVal[num]))
        if self.endVal[num] != -1:
            self.endEntry[num].setText(str(self.endVal[num]))

    def apply(self, maximum, minimum, newSpec):
        pass

    def applyFunc(self):
        if self.partIter == 0:
            if self.apply(np.array([0]), np.array([self.father.current.data1D.shape[-1]]), self.newSpec.isChecked()) is None:
                return False
        else:
            if self.apply(self.startVal[:self.partIter], self.endVal[:self.partIter], self.newSpec.isChecked()) is None:
                return False

############################################################


class integrateWindow(regionWindow):

    def __init__(self, parent):
        super(integrateWindow, self).__init__(parent, 'Integrate')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0): #Check for errors in the inputs
            self.father.father.dispMsg(self.NAME + ": wrong input")
            return None
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.integrate(minimum, maximum, newSpec)) is None:
                return None
        else:
            if self.father.current.data.noUndo:
                self.father.current.integrate(minimum, maximum, newSpec)
            else:
                returnValue = self.father.current.integrate(minimum, maximum, newSpec)
                if returnValue is None:
                    return None
                self.father.undoList.append(returnValue)
            self.father.redoList = []
            self.father.updAllFrames()
        return 1

############################################################


class sumWindow(regionWindow):

    def __init__(self, parent):
        super(sumWindow, self).__init__(parent, 'Sum')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            self.father.father.dispMsg(self.NAME + ": wrong input")
            return None
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.sum(minimum, maximum, newSpec)) is None:
                return None
        else:
            if self.father.current.data.noUndo:
                self.father.current.sum(minimum, maximum, newSpec)
            else:
                returnValue = self.father.current.sum(minimum, maximum, newSpec)
                if returnValue is None:
                    return None
                self.father.undoList.append(returnValue)
            self.father.redoList = []
            self.father.updAllFrames()
        return 1

############################################################


class maxWindow(regionWindow):

    def __init__(self, parent):
        super(maxWindow, self).__init__(parent, 'Max')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            self.father.father.dispMsg(self.NAME + ": wrong input")
            return None
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.maxMatrix(minimum, maximum, newSpec)) is None:
                return None
        else:
            if self.father.current.data.noUndo:
                self.father.current.maxMatrix(minimum, maximum, newSpec)
            else:
                returnValue = self.father.current.maxMatrix(minimum, maximum, newSpec)
                if returnValue is None:
                    return None
                self.father.undoList.append(returnValue)
            self.father.redoList = []
            self.father.updAllFrames()
        return 1

############################################################


class minWindow(regionWindow):

    def __init__(self, parent):
        super(minWindow, self).__init__(parent, 'Min')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            self.father.father.dispMsg(self.NAME + ": wrong input")
            return None
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.minMatrix(minimum, maximum, newSpec)) is None:
                return None
        else:
            if self.father.current.data.noUndo:
                self.father.current.minMatrix(minimum, maximum, newSpec)
            else:
                returnValue = self.father.current.minMatrix(minimum, maximum, newSpec)
                if returnValue is None:
                    return None
                self.father.undoList.append(returnValue)
            self.father.redoList = []
            self.father.updAllFrames()
        return 1

############################################################


class argmaxWindow(regionWindow):

    def __init__(self, parent):
        super(argmaxWindow, self).__init__(parent, 'Max position')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            self.father.father.dispMsg(self.NAME + ": wrong input")
            return None
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.argmaxMatrix(minimum, maximum, newSpec)) is None:
                return None
        else:
            if self.father.current.data.noUndo:
                returnValue = self.father.current.argmaxMatrix(minimum, maximum, newSpec)
            else:
                returnValue = self.father.current.argmaxMatrix(minimum, maximum, newSpec)
                if returnValue is None:
                    return None
                self.father.undoList.append(returnValue)
            self.father.redoList = []
            self.father.updAllFrames()
        return 1

############################################################


class argminWindow(regionWindow):

    def __init__(self, parent):
        super(argminWindow, self).__init__(parent, 'Min position')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            self.father.father.dispMsg(self.NAME + ": wrong input")
            return None
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.argminMatrix(minimum, maximum, newSpec)) is None:
                return None
        else:
            if self.father.current.data.noUndo:
                returnValue = self.father.current.argminMatrix(minimum, maximum, newSpec)
            else:
                returnValue = self.father.current.argminMatrix(minimum, maximum, newSpec)
                if returnValue is None:
                    return None
                self.father.undoList.append(returnValue)
            self.father.redoList = []
            self.father.updAllFrames()
        return 1

############################################################


class avgWindow(regionWindow):

    def __init__(self, parent):
        super(avgWindow, self).__init__(parent, 'Average')

    def apply(self, maximum, minimum, newSpec):
        if np.any(maximum < 0) or np.any(minimum < 0):
            self.father.father.dispMsg(self.NAME + ": wrong input")
            return None
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.average(minimum, maximum, newSpec)) is None:
                return None
        else:
            if self.father.current.data.noUndo:
                returnValue = self.father.current.average(minimum, maximum, newSpec)
            else:
                returnValue = self.father.current.average(minimum, maximum, newSpec)
                if returnValue is None:
                    return None
                self.father.undoList.append(returnValue)
            self.father.redoList = []
            self.father.updAllFrames()
        return 1

#############################################################


class regionWindow2(wc.ToolWindows):

    def __init__(self, parent, name, newSpecOption):
        self.NAME = name
        super(regionWindow2, self).__init__(parent)
        self.startVal = 0
        self.endVal = parent.current.data1D.shape[-1]
        self.grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.startEntry = wc.QLineEdit(self.startVal, self.checkValues)
        self.grid.addWidget(self.startEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End point:"), 2, 0)
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
            self.preview(self.startVal, self.endVal)
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
        self.preview(self.startVal, self.endVal)

    def applyFunc(self):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.startEntry.text())
        if inp is None:
            self.father.father.dispMsg(self.NAME + ": value not valid")
            return False
        self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        inp = safeEval(self.endEntry.text())
        if inp is None:
            self.father.father.dispMsg(self.NAME + ": value not valid")
            return False
        self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        if self.apply(self.startVal, self.endVal, self.newSpec.isChecked()) is None:
            return False

    def apply(self, maximum, minimum, newSpec):
        pass

############################################################


class extractRegionWindow(regionWindow2):

    def __init__(self, parent):
        super(extractRegionWindow, self).__init__(parent, 'Extract part', True)

    def apply(self, maximum, minimum, newSpec):
        if newSpec:
            if self.father.father.newWorkspace(self.father.current.getRegion(minimum, maximum, newSpec)) is None:
                return None
        else:
            if self.father.current.data.noUndo:
                self.father.current.getRegion(minimum, maximum, newSpec)
            else:
                returnValue = self.father.current.getRegion(minimum, maximum, newSpec)
                if returnValue is None:
                    return None
                    self.father.undoList.append(returnValue)
            self.father.redoList = []
            self.father.updAllFrames()
        return 1

############################################################


class SubtractAvgWindow(regionWindow2):

    def __init__(self, parent):
        super(SubtractAvgWindow, self).__init__(parent, 'Subtract Avg', False)

    def apply(self, maximum, minimum, newSpec):
        if self.father.current.data.noUndo:
            self.father.current.subtractAvg(maximum, minimum)
        else:
            returnValue = self.father.current.subtractAvg(maximum, minimum)
            if returnValue is None:
                return None
            if not self.father.current.data.noUndo:
                self.father.undoList.append(returnValue)
        self.father.redoList = []
        self.father.updAllFrames()
        return 1

    def preview(self, maximum, minimum):
        self.father.current.subtractAvgPreview(maximum, minimum)

#############################################################


class FiddleWindow(wc.ToolWindows):

    NAME = "Reference deconvolution"

    def __init__(self, parent):
        super(FiddleWindow, self).__init__(parent)
        self.startVal = 0
        self.endVal = parent.current.data1D.shape[-1]
        self.grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.startEntry = wc.QLineEdit(self.startVal, self.checkValues)
        self.grid.addWidget(self.startEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End point:"), 2, 0)
        self.endEntry = wc.QLineEdit(self.endVal, self.checkValues)
        self.grid.addWidget(self.endEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Linebroadening [Hz]:"), 4, 0)
        self.lbEntry = wc.QLineEdit("1.0", self.checkValues)
        self.grid.addWidget(self.lbEntry, 5, 0)
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

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

    def applyFunc(self):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.startEntry.text())
        if inp is None:
            self.father.father.dispMsg("Reference deconv: start entry not valid")
            return False
        self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        inp = safeEval(self.endEntry.text())
        if inp is None:
            self.father.father.dispMsg("Reference deconv: end entry not valid")
            return False
        self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        lb = safeEval(self.lbEntry.text())
        if lb is None:
            self.father.father.dispMsg("Reference deconv: Linebroadening entry not valid")
            return False
        if self.father.current.data.noUndo:
            self.father.current.fiddle(self.startVal, self.endVal, lb)
        else:
            returnValue = self.father.current.fiddle(self.startVal, self.endVal, lb)
            if returnValue is None:
                return None
            self.father.undoList.append(returnValue)
        self.father.redoList = []

##############################################################


class DeleteWindow(wc.ToolWindows):

    NAME = "Delete"

    def __init__(self, parent):
        super(DeleteWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Indexes to delete:"), 0, 0)
        self.delEntry = wc.QLineEdit('0', self.preview)
        self.grid.addWidget(self.delEntry, 1, 0)

    def preview(self, *args):
        length = int(self.father.current.data1D.shape[-1])
        pos = safeEval(self.delEntry.text())
        if pos == None:
            self.father.father.dispMsg('Delete: not all values are valid indexes to delete')
            return False
        pos = np.array(pos)
        pos[pos < 0] = pos[pos < 0] + length
        if (pos > -1).all() and (pos < length).all():
            self.father.current.deletePreview(pos)
        else:
            self.father.father.dispMsg('Delete: not all values are valid indexes to delete')

    def applyFunc(self):
        length = int(self.father.current.data1D.shape[-1])
        pos = safeEval(self.delEntry.text())
        if pos == None:
            self.father.father.dispMsg('Delete: not all values are valid indexes to delete')
            return False
        pos = np.array(pos)
        pos[pos < 0] = pos[pos < 0] + length
        if (pos > -1).all() and (pos < length).all():
            self.father.redoList = []
            if self.father.current.data.noUndo:
                self.father.current.delete(pos)
            else:
                self.father.undoList.append(self.father.current.delete(pos))
        else:
            self.father.father.dispMsg('Delete: not all values are valid indexes to delete')
            return False

##############################################################


class SplitWindow(wc.ToolWindows):

    NAME = "Split"

    def __init__(self, parent):
        super(SplitWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Sections:"), 0, 0)
        self.splitEntry = wc.QLineEdit('1', self.preview)
        self.grid.addWidget(self.splitEntry, 1, 0)

    def preview(self, *args):
        val = safeEval(self.splitEntry.text(), self.father.current.data1D.shape[-1])
        if val is None:
            self.father.father.dispMsg("Split: input not valid")
            return False
        else:
            self.splitEntry.setText(str(int(round(val))))

    def applyFunc(self):
        val = safeEval(self.splitEntry.text(), self.father.current.data1D.shape[-1])
        if val is None:
            self.father.father.dispMsg("Split: input not valid")
            return False
        val = int(val)
        if val <= 0:
            self.father.father.dispMsg("Split: input not valid")
            return False
        if self.father.current.data.noUndo:
            self.father.current.split(int(round(val)))
        else:
            returnValue = self.father.current.split(int(round(val)))
            if returnValue is None:
                return False
            self.father.undoList.append(returnValue)
        self.father.redoList = []

##############################################################


class ConcatenateWindow(wc.ToolWindows):

    NAME = "Concatenate"

    def __init__(self, parent):
        super(ConcatenateWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Concatenation axes:"), 0, 0)
        self.axesEntry = QtWidgets.QComboBox()
        self.axesEntry.addItems(np.array(np.arange(self.father.current.data.data.ndim - 1) + 1, dtype=str))
        self.grid.addWidget(self.axesEntry, 1, 0)

    def applyFunc(self):
        if self.father.current.data.noUndo:
            self.father.current.concatenate(self.axesEntry.currentIndex())
        else:
            returnValue = self.father.current.concatenate(self.axesEntry.currentIndex())
            if returnValue is None:
                return
            self.father.undoList.append(returnValue)
        self.father.redoList = []

##############################################################


class InsertWindow(wc.ToolWindows):

    NAME = "Insert"

    def __init__(self, parent):
        super(InsertWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Start insert at index:"), 0, 0)
        self.posEntry = wc.QLineEdit(self.father.current.data1D.shape[-1], self.preview)
        self.grid.addWidget(self.posEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("Workspace to insert:"), 2, 0)
        self.wsEntry = QtWidgets.QComboBox()
        self.wsEntry.addItems(self.father.father.workspaceNames)
        self.grid.addWidget(self.wsEntry, 3, 0)

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

    def applyFunc(self):
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
        if self.father.current.data.noUndo:
            self.father.current.insert(self.father.father.workspaces[ws].masterData.data, pos)
        else:
            self.father.undoList.append(self.father.current.insert(self.father.father.workspaces[ws].masterData.data, pos))

##############################################################


class CombineWindow(wc.ToolWindows):

    SINGLESLICE = True
    RESIZABLE = True
    
    def __init__(self, parent, combType):
        super(CombineWindow, self).__init__(parent)
        self.combType = combType # 0 = add, 1 = subtract, 2 = multiply, 3 = divide
        if self.combType is 0:
            self.setWindowTitle("Add")        
            self.grid.addWidget(wc.QLabel("Workspace to add:"), 0, 0)
        elif self.combType is 1:
            self.setWindowTitle("Subtract")        
            self.grid.addWidget(wc.QLabel("Workspace to subtract:"), 0, 0)
        elif self.combType is 2:
            self.setWindowTitle("Multiply")        
            self.grid.addWidget(wc.QLabel("Workspace to multiply:"), 0, 0)
        elif self.combType is 3:
            self.setWindowTitle("Divide")        
            self.grid.addWidget(wc.QLabel("Workspace to divide:"), 0, 0)
        self.wsEntry = QtWidgets.QComboBox()
        self.wsEntry.addItems(self.father.father.workspaceNames)
        self.grid.addWidget(self.wsEntry, 1, 0)

    def applyFunc(self):
        ws = self.wsEntry.currentIndex()
        if self.combType is 0:
            returnValue = self.father.current.add(self.father.father.workspaces[ws].masterData.data, self.singleSlice.isChecked())
        elif self.combType is 1:
            returnValue = self.father.current.subtract(self.father.father.workspaces[ws].masterData.data, self.singleSlice.isChecked())
        elif self.combType is 2:
            returnValue = self.father.current.multiplySpec(self.father.father.workspaces[ws].masterData.data, self.singleSlice.isChecked())
        elif self.combType is 3:
            returnValue = self.father.current.divideSpec(self.father.father.workspaces[ws].masterData.data, self.singleSlice.isChecked())
        if returnValue is None and not self.father.current.data.noUndo:
            return
        self.father.redoList = []
        if not self.father.current.data.noUndo:
            self.father.undoList.append(returnValue)
            
##############################################################


class SNWindow(wc.ToolWindows):

    NAME = "Signal to noise"
    CANCELNAME = "&Close"
    OKNAME = "C&alc"

    def __init__(self, parent):
        super(SNWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Start point noise:"), 0, 0)
        self.minNoiseEntry = wc.QLineEdit('0', self.checkValues)
        self.grid.addWidget(self.minNoiseEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End point noise:"), 2, 0)
        self.maxNoiseEntry = wc.QLineEdit(parent.current.data1D.shape[-1], self.checkValues)
        self.grid.addWidget(self.maxNoiseEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Start point signal:"), 4, 0)
        self.minEntry = wc.QLineEdit('0', self.checkValues)
        self.grid.addWidget(self.minEntry, 5, 0)
        self.grid.addWidget(wc.QLabel("End point signal:"), 6, 0)
        self.maxEntry = wc.QLineEdit(parent.current.data1D.shape[-1], self.checkValues)
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
        self.applyFunc()

    def applyFunc(self):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.minNoiseEntry.text())
        if inp is None:
            self.father.father.dispMsg("S/N: invalid range")
            return False
        minimumNoise = int(round(inp))
        if minimumNoise < 0:
            minimumNoise = 0
        elif minimumNoise > dataLength:
            minimumNoise = dataLength
        self.minNoiseEntry.setText(str(minimumNoise))
        inp = safeEval(self.maxNoiseEntry.text())
        if inp is None:
            self.father.father.dispMsg("S/N: invalid range")
            return False
        maximumNoise = int(round(inp))
        if maximumNoise < 0:
            maximumNoise = 0
        elif maximumNoise > dataLength:
            maximumNoise = dataLength
        self.maxNoiseEntry.setText(str(maximumNoise))
        inp = safeEval(self.minEntry.text())
        if inp is None:
            self.father.father.dispMsg("S/N: invalid range")
            return False
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            self.father.father.dispMsg("S/N: invalid range")
            return False
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.snEntry.setText(str(self.father.current.SN(minimumNoise, maximumNoise, minimum, maximum)))
        return False #Return to keep window

##############################################################


class FWHMWindow(wc.ToolWindows):

    NAME = "FWHM"
    CANCELNAME = "&Close"
    OKNAME = "C&alc"

    def __init__(self, parent):
        super(FWHMWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.minEntry = wc.QLineEdit('0', self.checkValues)
        self.grid.addWidget(self.minEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End point:"), 2, 0)
        self.maxEntry = wc.QLineEdit(parent.current.data1D.shape[-1], self.checkValues)
        self.grid.addWidget(self.maxEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Units:"), 4, 0)
        unitSelect = self.father.current.axType
        if self.father.current.spec == 1:
            unitList = ['Hz', 'kHz', 'MHz', 'ppm']
            if self.father.current.ppm:
                unitSelect = 3
        else:
            unitList = ['s', 'ms', u'\u03BCs']
        self.unitDrop = QtWidgets.QComboBox()
        self.unitDrop.addItems(unitList)
        self.unitDrop.setCurrentIndex(unitSelect)
        self.unitDrop.currentIndexChanged.connect(self.checkValues)
        self.grid.addWidget(self.unitDrop, 5, 0)
        self.grid.addWidget(wc.QLabel(u"FWHM:"), 6, 0)
        self.fwhmEntry = wc.QLineEdit("0.0")
        self.grid.addWidget(self.fwhmEntry, 7, 0)
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
        self.applyFunc()

    def applyFunc(self):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.minEntry.text())
        if inp is None:
            self.father.father.dispMsg("FWHM: invalid range")
            return False
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            self.father.father.dispMsg("FWHM: invalid range")
            return False
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.fwhmEntry.setText(str(self.father.current.fwhm(minimum, maximum, self.unitDrop.currentIndex())))
        return False #Return to keep window

##############################################################


class COMWindow(wc.ToolWindows):  # Centre of Mass Window

    NAME = "Centre of Mass"
    CANCELNAME = "&Close"
    OKNAME = "C&alc"

    def __init__(self, parent):
        super(COMWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Start point:"), 0, 0)
        self.minEntry = wc.QLineEdit("0", self.checkValues)
        self.grid.addWidget(self.minEntry, 1, 0)
        self.grid.addWidget(wc.QLabel("End point:"), 2, 0)
        self.maxEntry = wc.QLineEdit(parent.current.data1D.shape[-1], self.checkValues)
        self.grid.addWidget(self.maxEntry, 3, 0)
        if self.father.current.spec == 1:
            if self.father.current.ppm:
                self.grid.addWidget(wc.QLabel("Centre of Mass [ppm]:"), 4, 0)
            else:
                if self.father.current.axType == 0:
                    self.grid.addWidget(wc.QLabel("Centre of Mass [Hz]:"), 4, 0)
                elif self.father.current.axType == 1:
                    self.grid.addWidget(wc.QLabel("Centre of Mass [kHz]:"), 4, 0)
                elif self.father.current.axType == 2:
                    self.grid.addWidget(wc.QLabel("Centre of Mass [MHz]:"), 4, 0)
                elif self.father.current.axType == 3:
                    self.grid.addWidget(wc.QLabel("Centre of Mass [ppm]:"), 4, 0)
        else:
            if self.father.current.axType == 0:
                self.grid.addWidget(wc.QLabel("Centre of Mass [s]:"), 4, 0)
            elif self.father.current.axType == 1:
                self.grid.addWidget(wc.QLabel("Centre of Mass [ms]:"), 4, 0)
            elif self.father.current.axType == 2:
                self.grid.addWidget(wc.QLabel(u"Centre of Mass [\u03bcs]:"), 4, 0)
        self.comEntry = wc.QLineEdit("0.0")
        self.grid.addWidget(self.comEntry, 5, 0)
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
        self.applyFunc()

    def applyFunc(self):
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.minEntry.text())
        if inp is None:
            self.father.father.dispMsg("Centre of Mass: invalid range")
            return False
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            self.father.father.dispMsg("Centre of Mass: invalid range")
            return False
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.comEntry.setText(str(self.father.current.COM(minimum, maximum)))
        return False #Return to keep window

##########################################################################################


class ReorderWindow(wc.ToolWindows):

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
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.LastLocation)
        if type(filename) is tuple:
            filename = filename[0]        
        if filename:  # if not cancelled
            self.father.father.LastLocation = os.path.dirname(filename)  # Save used path
        if len(filename) == 0:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def applyFunc(self):
        newLength = self.lengthEntry.text()
        if newLength == '':
            newLength = None
        else:
            newLength = safeEval(self.lengthEntry.text())
            if newLength is None:
                self.father.father.dispMsg("Reorder: `Length' input is not valid")
                return False
        val = safeEval(self.valEntry.text(), int(self.father.current.data1D.shape[-1]))
        if not isinstance(val, (list, np.ndarray)):
            self.father.father.dispMsg("Reorder: `Positions' input is not a list or array")
            return False
        if len(val) != self.father.current.data1D.shape[-1]:
            self.father.father.dispMsg("Reorder: length of input does not match length of data")
            return False
        val = np.array(val, dtype=int)
        self.father.redoList = []
        check = self.father.current.reorder(val, newLength)
        if check is None:
            return False

        if not self.father.masterData.noUndo:
            self.father.undoList.append(check)
        return

##########################################################################################


class RegridWindow(wc.ToolWindows):

    NAME = "Regrid"

    def __init__(self, parent ):
        super(RegridWindow, self).__init__(parent)
        self.typeDrop = QtWidgets.QComboBox(parent=self)
        self.typeDrop.addItems(["Min/max input"])
        self.grid.addWidget(self.typeDrop, 0, 0, 1, 2)
        self.maxValue = wc.QLineEdit(10)
        #Get unit
        if self.father.current.spec == 1:
            if self.father.current.ppm:
                self.unit = 'ppm'
            else:
                if self.father.current.axType == 0:
                    self.unit = 'Hz'
                elif self.father.current.axType == 1:
                    self.unit = 'kHz'
                elif self.father.current.axType == 2:
                    self.unit = 'MHz'
                elif self.father.current.axType == 3:
                    self.unit = 'ppm'
            self.comEntry = wc.QLineEdit("0.0")
            self.maxLabel = wc.QLeftLabel('Max [' + self.unit + ']:')
            self.minValue = wc.QLineEdit(0)
            self.minLabel = wc.QLeftLabel('Min [' + self.unit + ']:')
            self.points = wc.QLineEdit(1000)
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
        maxVal = safeEval(self.maxValue.text(), type = 'FI')
        if maxVal is None:
            self.father.father.dispMsg("Regrid: 'Max' input not valid")
            return False
        minVal = safeEval(self.minValue.text(), type = 'FI')
        if minVal is None:
            self.father.father.dispMsg("Regrid: 'Min' input not valid")
            return False
        numPoints = safeEval(self.points.text(), type = 'FI')
        if numPoints is None:
            self.father.father.dispMsg("Regrid: '# of points' input not valid")
            return False
        numPoints = int(numPoints)

        #Convert to Hz/s

        if self.unit == 'kHz':
            maxVal *= 1e3
            minVal *= 1e3
        elif self.unit == 'MHz':
            maxVal *= 1e6
            minVal *= 1e6
        elif self.unit == 'ppm':
            maxVal *= self.father.masterData.ref[self.father.current.axes] / 1e6
            minVal *= self.father.masterData.ref[self.father.current.axes] / 1e6
   
        if self.father.current.data.noUndo:
            self.father.current.regrid([minVal,maxVal],numPoints)
        else:
            self.father.undoList.append(self.father.current.regrid([minVal,maxVal],numPoints))
        return

##########################################################################################


class FFMWindow(wc.ToolWindows):

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
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.LastLocation)
        if type(filename) is tuple:
            filename = filename[0]        
        if filename:  # if not cancelled
            self.father.father.LastLocation = os.path.dirname(filename)  # Save used path
        if len(filename) == 0:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def applyFunc(self):
        val = safeEval(self.valEntry.text())
        if not isinstance(val, (list, np.ndarray)):
            self.father.father.dispMsg("FFM: 'Positions' is not a list or array")
            return False
        val = np.array(val, dtype=int)
        self.father.redoList = []
        check = self.father.current.ffm(val, self.typeDrop.currentIndex())
        if check is None:
            self.father.father.dispMsg("FFM: error",color = 'red')
            return False
        if not self.father.masterData.noUndo:
            self.father.undoList.append(check)

##########################################################################################


class CLEANWindow(wc.ToolWindows):

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
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.LastLocation)
        if type(filename) is tuple:
            filename = filename[0]        
        if filename:  # if not cancelled
            self.father.father.LastLocation = os.path.dirname(filename)  # Save used path
        if len(filename) == 0:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def applyFunc(self):
        val = safeEval(self.valEntry.text())
        if not isinstance(val, (list, np.ndarray)):
            self.father.father.dispMsg("CLEAN: 'Positions' is not a list or array")
            return False
        val = np.array(val, dtype=int)
        gamma = safeEval(self.gammaEntry.text())
        if gamma is None:
            self.father.father.dispMsg("CLEAN: 'Gamma' input is not valid")
            return False
        threshold = safeEval(self.thresholdEntry.text())
        if threshold is None:
            self.father.father.dispMsg("CLEAN: 'Threshold' input is not valid")
            return False
        threshold = threshold
        maxIter = safeEval(self.maxIterEntry.text())
        if maxIter is None:
            self.father.father.dispMsg("CLEAN: 'Max. iter.' is not valid")
            return False
        maxIter = int(maxIter)
        self.father.redoList = []
        check = self.father.current.clean(val, self.typeDrop.currentIndex(), gamma, threshold, maxIter)
        if check is None:
            self.father.father.father.dispMsg("CLEAN: error",color = 'red')
            return False
        if not self.father.masterData.noUndo:
            self.father.undoList.append(check)

################################################################


class ISTWindow(wc.ToolWindows):

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
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File', self.father.father.LastLocation)
        if type(filename) is tuple:
            filename = filename[0]        
        if filename:  # if not cancelled
            self.father.father.LastLocation = os.path.dirname(filename)  # Save used path
        if len(filename) == 0:
            return
        self.valEntry.setText(repr(np.loadtxt(filename, dtype=int)))

    def applyFunc(self):
        val = safeEval(self.valEntry.text())
        if not isinstance(val, (list, np.ndarray)):
            self.father.father.dispMsg("IST: 'Positions' input is not a list or array")
            return False
        val = np.array(val, dtype=int)
        tracelimit = safeEval(self.tracelimitEntry.text()) 
        if tracelimit is None:
            self.father.father.dispMsg("IST: 'Residual' input is not valid")
            return False
        tracelimit /= 100
        threshold = safeEval(self.thresholdEntry.text())
        if threshold is None:
            self.father.father.dispMsg("IST: 'Threshold' input is not valid")
            return False
        maxIter = safeEval(self.maxIterEntry.text())
        if maxIter is None:
            self.father.father.dispMsg("IST: 'Max. iter.' input is not valid")
            return False
        maxIter = int(maxIter)
        self.father.redoList = []
        check = self.father.current.ist(val, self.typeDrop.currentIndex(), threshold, maxIter,tracelimit)
        if check is None:
            self.father.father.father.dispMsg("IST: error",color = 'red')
            return False
        if not self.father.masterData.noUndo:
            self.father.undoList.append(check)
        
################################################################


class ShearingWindow(wc.ToolWindows):

    NAME = "Shearing"

    def __init__(self, parent):
        super(ShearingWindow, self).__init__(parent)
        options = list(map(str, range(1, self.father.masterData.data.ndim + 1)))
        self.grid.addWidget(wc.QLabel("Shearing constant:"), 0, 0)
        self.shearDropdown = QtWidgets.QComboBox()
        self.shearDropdown.addItems(['User Defined','Spin 3/2, -3Q (7/9)','Spin 5/2, 3Q (19/12)','Spin 5/2, -5Q (25/12)','Spin 7/2, 3Q (101/45)',
                                     'Spin 7/2, 5Q (11/9)','Spin 7/2, -7Q (161/45)','Spin 9/2, 3Q (91/36)','Spin 9/2, 5Q (95/36)','Spin 9/2, 7Q (7/18)','Spin 9/2, -9Q (31/6)'])
        self.shearDropdown.activated.connect(self.dropdownChanged)
        self.shearList = [0,7.0/9.0,19.0/12.0,25.0/12.0,101.0/45.0,11.0/9.0,161.0/45.0,91.0/36.0,95.0/36.0,7.0/18.0,31.0/6.0]
        self.grid.addWidget(self.shearDropdown, 1, 0)
        self.shearEntry = wc.QLineEdit("0.0", self.shearPreview)
        self.grid.addWidget(self.shearEntry, 3, 0)
        self.grid.addWidget(wc.QLabel("Shearing direction:"), 4, 0)
        self.dirEntry = QtWidgets.QComboBox()
        self.dirEntry.addItems(options)
        self.dirEntry.setCurrentIndex(self.father.masterData.data.ndim - 2)
        self.grid.addWidget(self.dirEntry, 5, 0)
        self.grid.addWidget(wc.QLabel("Shearing axis:"), 6, 0)
        self.axEntry = QtWidgets.QComboBox()
        self.axEntry.addItems(options)
        self.axEntry.setCurrentIndex(self.father.masterData.data.ndim - 1)
        self.grid.addWidget(self.axEntry,7, 0)

    def dropdownChanged(self):
        index =  self.shearDropdown.currentIndex()
        self.shearEntry.setText("%.9f" % self.shearList[index])
        
    def shearPreview(self, *args):
        shear = safeEval(self.shearEntry.text())
        if shear is not None:
            self.shearEntry.setText(str(float(shear)))

    def applyFunc(self):
        shear = safeEval(self.shearEntry.text())
        if shear is None:
            self.father.father.dispMsg("Shearing: 'constant' not a valid value")
            return False
        axes = self.dirEntry.currentIndex()
        axes2 = self.axEntry.currentIndex()
        if axes == axes2:
            self.father.father.dispMsg("Shearing: axes cannot be the same for shearing")
            return False
        else:
            self.father.redoList = []
            if self.father.masterData.noUndo:
               self.father.current.shearing(float(shear), axes, axes2)
            else:
                self.father.undoList.append(self.father.current.shearing(float(shear), axes, axes2))

##########################################################################################


class MultiplyWindow(wc.ToolWindows):

    NAME = "Multiply"
    SINGLESLICE = True

    def __init__(self, parent):
        super(MultiplyWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Values:"), 0, 0)
        self.valEntry = wc.QLineEdit('', self.preview)
        self.grid.addWidget(self.valEntry, 1, 0)

    def preview(self, *args):
        val = safeEval(self.valEntry.text())
        if val is None:
            self.father.father.dispMsg("Multiply: input not valid")
            return False
        check = self.father.current.multiplyPreview(np.array(val))
        if check is not True:
            self.father.father.dispMsg("Multiply: " + check)

    def applyFunc(self):
        val = safeEval(self.valEntry.text())
        if val is None:
            self.father.father.dispMsg("Multiply: input not valid")
            return False
        if self.father.current.data.noUndo:
            self.father.current.multiply(np.array(val), self.singleSlice.isChecked())
        else:
            returnValue = self.father.current.multiply(np.array(val), self.singleSlice.isChecked())
            if returnValue is None:
                return False
            self.father.undoList.append(returnValue)
        self.father.redoList = []

##########################################################################################


class XaxWindow(wc.ToolWindows):

    RESIZABLE = True
    NAME = "User defined x-axis"

    def __init__(self, parent):
        super(XaxWindow, self).__init__(parent)
        self.axisSize = int(self.father.current.data1D.shape[-1])        
        self.grid.addWidget(wc.QLabel("Input x-axis values:"), 0, 0, 1, 2) 
        self.typeDropdown = QtWidgets.QComboBox()
        self.typeDropdown.addItems(['Expression','Linear','Logarithmic'])
        self.typeDropdown.activated.connect(self.typeChanged)
        self.grid.addWidget(self.typeDropdown, 1, 0,1 ,2)    
        self.exprEntry = wc.QLineEdit('', self.xaxPreview)
        self.grid.addWidget(self.exprEntry, 2, 0, 1, 2)
        
        #Linear 
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
        
        self.linStopEntry = wc.QLineEdit('',self.xaxPreview)
        self.linStopEntry.setMaximumWidth(120)
        self.linStopEntry.hide()
        self.grid.addWidget(self.linStopEntry, 4, 1, 1, 1)
        
        #Log
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
        
        self.table = QtWidgets.QTableWidget(self.axisSize,2)
        self.table.setHorizontalHeaderLabels(['Index','Value [s]'])
        self.table.verticalHeader().hide()
        for val in range(self.axisSize):
            item = QtWidgets.QTableWidgetItem(str(val))
            item.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table.setItem(int(val),0,item)
            item2 = QtWidgets.QTableWidgetItem('')
            item2.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table.setItem(int(val),1,item2)
            
#        self.table.setVerticalHeaderLabels([str(a) for a in range(self.axisSize)])
        self.grid.addWidget(self.table, 12, 0, 1, 2)
        self.resize(250, 500)

    def typeChanged(self,index):
        if index == 0: #If expr
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
            env['length'] = int(self.father.current.data1D.shape[-1])  # so length can be used to in equations
            env['euro'] = lambda fVal, num=self.axisSize: func.euro(fVal, num)
            try:
                val = np.array(eval(self.exprEntry.text(), env))                # find a better solution, also add catch for exceptions
            except:
                try:
                    val = np.fromstring(self.exprEntry.text(),sep=' ')
                    val2 = np.fromstring(self.exprEntry.text(),sep=',')
                    if len(val2) > len(val):
                        val = val2
                except:
                    val = None
            if not isinstance(val, (list, np.ndarray)):
                self.father.father.dispMsg("X-axis: Input is not a list or array")
                return
            if len(val) != self.father.current.data1D.shape[-1]:
                self.father.father.dispMsg("X-axis: Length of input does not match length of data")
                return
            if not all(isinstance(x, (int, float)) for x in val):
                self.father.father.dispMsg("X-axis: Array is not all of int or float type")
                return
        elif self.typeDropdown.currentIndex() == 1:
            start = safeEval(self.linStartEntry.text())
            stop = safeEval(self.linStopEntry.text())
            if start is None:
                self.father.father.dispMsg("X-axis: linear start value is not valid")
                return
            if stop is None:
                self.father.father.dispMsg("X-axis: linear stop value is not valid")
                return
            val = np.linspace(start,stop,self.axisSize)
        elif self.typeDropdown.currentIndex() == 2:
            start = safeEval(self.logStartEntry.text())
            stop = safeEval(self.logStopEntry.text())
            if start is None or start <= 0.0:
                self.father.father.dispMsg("X-axis: logarithmic start value is not valid")
                return
            if stop is None or stop <= 0.0:
                self.father.father.dispMsg("X-axis: logarithmic stop value is not valid")
                return
            val = np.logspace(np.log10(start),np.log10(stop),self.axisSize)
        return val
    
    def xaxPreview(self, *args):
        val = self.getValues()
        if val is None: #if error return. Messages are handled by the called function
            return
        for i in range(self.axisSize):
            item = QtWidgets.QTableWidgetItem('{:.6g}'.format(val[i]))
            item.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table.setItem(i,1,item)
        self.father.current.setXaxPreview(np.array(val))

    def applyFunc(self):
        val = self.getValues()
        if val is None: #if error return. Messages are handled by the called function
            return  
        self.father.redoList = []
        if self.father.current.data.noUndo:
            self.father.current.setXax(np.array(val))
        else:
            self.father.undoList.append(self.father.current.setXax(np.array(val)))

##########################################################################################


class RefWindow(wc.ToolWindows):

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
        if parent.current.spec == 0:
            self.father.father.dispMsg('Setting ppm is only available for frequency data')
            self.deleteLater()
            return
        self.grid.addWidget(wc.QLabel("Name:"), 0, 0)
        self.refName = wc.QLineEdit()
        self.grid.addWidget(self.refName, 1, 0)
        self.grid.addWidget(wc.QLabel("Frequency [MHz]:"), 2, 0)
        self.freqEntry = wc.QLineEdit(("%.7f" % (self.father.current.ref * 1e-6)), self.preview)
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
        freq = safeEval(self.freqEntry.text())
        ref = safeEval(self.refEntry.text())
        if freq is None or ref is None:
            return
        self.freqEntry.setText("%.7f" % (freq))
        self.refEntry.setText(str(ref))

    def fillSecondaryRef(self):
        self.refEntry.setText(self.secRefValues[self.refSecond.currentIndex()])

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
            if givenname in self.father.father.referenceName:  # if exists
                self.father.father.dispMsg("Reference name '" + givenname + "' already exists")
                nameOK = False
            else:
                self.father.mainProgram.referenceAdd(reffreq, givenname)
        if nameOK:
            self.father.redoList = []
            if self.father.current.data.noUndo:
                self.father.current.setRef(reffreq)
            else:
                self.father.undoList.append(self.father.current.setRef(reffreq))
            self.closeEvent()

    def picked(self, pos):
        self.freqEntry.setText("%.7f" % ((self.father.current.ref + self.father.current.xax[pos[0]]) * 1e-6))
        self.father.current.peakPickFunc = lambda pos, self=self: self.picked(pos)
        self.father.current.peakPick = True

##########################################################################################


class HistoryWindow(wc.ToolWindows):

    NAME = "Processing history"
    RESIZABLE = True
    MENUDISABLE = False

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


class CombineWorkspaceWindow(wc.ToolWindows):

    NAME = "Combine workspaces"
    RESIZABLE = True

    def __init__(self, parent):
        super(CombineWorkspaceWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Workspaces:"), 0, 0)
        self.grid.addWidget(wc.QLabel("Combined spectrum:"), 0, 1)
        self.listA = OrigListWidget(self)
        for i in self.father.workspaceNames:
            QtWidgets.QListWidgetItem(i, self.listA).setToolTip(i)
        self.listB = DestListWidget(self)
        self.grid.addWidget(self.listA, 1, 0)
        self.grid.addWidget(self.listB, 1, 1)
        self.layout.setColumnStretch(2, 1)
        self.resize(500, 400)

    def applyFunc(self, *args):
        items = []
        for index in range(self.listB.count()):
            items.append(self.listB.item(index).text())
        if len(items) == 0:
            self.father.dispMsg("Please select at least one workspace to combine")
        else:
            self.father.combineWorkspace(items)

    def closeEvent(self, *args):
        if self.MENUDISABLE:
            self.father.menuEnable()
        self.deleteLater()

##########################################################################################


class CombineLoadWindow(wc.ToolWindows):

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
        fileList = QtWidgets.QFileDialog.getOpenFileNames(self, 'Open File', self.father.LastLocation)
        if type(fileList) is tuple:
            fileList = fileList[0]
        for filePath in fileList:
            if filePath:  # if not cancelled
                self.father.LastLocation = os.path.dirname(filePath)  # Save used path
            if len(filePath) == 0:
                return
            self.specList.addItem(filePath)

    def applyFunc(self, *args):
        items = []
        for index in range(self.specList.count()):
            items.append(self.specList.item(index).text())
        if len(items) == 0:
            self.father.dispMsg("Please select at least one workspace to combine")
        else:
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
        
        grid.addWidget(wc.QLabel("Delay [s]:"), 2, 0)
        self.delTime = QtWidgets.QDoubleSpinBox()
        self.delTime.setMaximum(10000)
        self.delTime.setMinimum(0)
        self.delTime.setSingleStep(0.1)
        self.delTime.setValue(0.5)
        grid.addWidget(self.delTime, 2, 1)
        cancelButton = QtWidgets.QPushButton("&Close")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 3, 0)
        watchButton = QtWidgets.QPushButton("&Watch")
        watchButton.clicked.connect(self.applyAndClose)
        layout.addWidget(watchButton, 3, 1)
        unwatchButton = QtWidgets.QPushButton("&Unwatch")
        unwatchButton.clicked.connect(self.stopAndClose)
        layout.addWidget(unwatchButton, 3, 2)
        layout.setColumnStretch(4, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def applyAndClose(self, *args):
        self.father.stopMonitor()
        items = []
        for index in range(self.listB.count()):
            items.append(self.listB.item(index).text())
        delay = self.delTime.value()
        self.father.startMonitor(items,delay)
        self.closeEvent()
        
    def stopAndClose(self, *args):
        self.father.stopMonitor()
        self.closeEvent()
    
    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()
        
##############################################################################


class PlotSettingsWindow(wc.ToolWindows):

    NAME = "Preferences"

    def __init__(self, parent):
        super(PlotSettingsWindow, self).__init__(parent)
        tabWidget = QtWidgets.QTabWidget()
        tab1 = QtWidgets.QWidget()
        tab2 = QtWidgets.QWidget()
        tabWidget.addTab(tab1, "Plot")
        tabWidget.addTab(tab2, "Contour")
        grid1 = QtWidgets.QGridLayout()
        grid2 = QtWidgets.QGridLayout()
        tab1.setLayout(grid1)
        tab2.setLayout(grid2)
        grid1.setColumnStretch(10, 1)
        grid1.setRowStretch(10, 1)
        grid2.setColumnStretch(10, 1)
        grid2.setRowStretch(10, 1)

        grid1.addWidget(QtWidgets.QLabel("Linewidth:"), 1, 0)
        self.lwSpinBox = QtWidgets.QDoubleSpinBox()
        self.lwSpinBox.setSingleStep(0.1)
        self.lwSpinBox.setValue(self.father.current.linewidth)
        self.lwSpinBox.valueChanged.connect(self.preview)
        grid1.addWidget(self.lwSpinBox, 1, 1)
        self.color = self.father.current.color
        lineColorButton = QtWidgets.QPushButton("Line colour")
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

        grid2.addWidget(QtWidgets.QLabel("Colourmap:"), 0, 0)
        self.cmEntry = QtWidgets.QComboBox(self)
        self.cmEntry.addItems(sc.COLORMAPLIST)
        self.cmEntry.setCurrentIndex(sc.COLORMAPLIST.index(self.father.current.colorMap))
        self.cmEntry.currentIndexChanged.connect(self.preview)
        grid2.addWidget(self.cmEntry, 0, 1)
        self.constColorCheck = QtWidgets.QCheckBox("Constant colours")
        self.constColorCheck.setChecked(self.father.current.contourConst)
        grid2.addWidget(self.constColorCheck, 1, 0)
        self.constColorCheck.stateChanged.connect(self.preview)
        self.posColor = self.father.current.contourColors[0]
        posColorButton = QtWidgets.QPushButton("Positive colour")
        posColorButton.clicked.connect(self.setPosColor)
        grid2.addWidget(posColorButton, 2, 0)
        self.negColor = self.father.current.contourColors[1]
        negColorButton = QtWidgets.QPushButton("Negative colour")
        negColorButton.clicked.connect(self.setNegColor)
        grid2.addWidget(negColorButton, 3, 0)
        self.grid.addWidget(tabWidget, 0, 0)

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
        self.father.current.setGrids([self.xgridCheck.isChecked(), self.ygridCheck.isChecked()])
        self.father.current.setColorMap(self.cmEntry.currentIndex())
        self.father.current.setContourConst(self.constColorCheck.isChecked())
        self.father.current.setContourColors([self.posColor, self.negColor])

##############################################################################


class errorWindow(wc.ToolWindows):

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
                tmp = QtWidgets.QListWidgetItem(error[0] + ': Python error', self.errorQList)
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
        
    def rowChange(self,row):
        errorText = ''
        error = self.father.errors[row]
        if len(error[1]) == 3 :
            errorText = errorText + error[0] + '<br>'
            for line in tb.format_exception(error[1][0],error[1][1],error[1][2]):
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
        self.toolbarCheck = QtWidgets.QCheckBox("Show Shortcut Toolbar")
        self.toolbarCheck.setChecked(self.father.defaultToolBar)
        grid1.addWidget(self.toolbarCheck, 5, 0, 1, 2)
        editToolbarButton = QtWidgets.QPushButton("Edit Toolbar")
        editToolbarButton.clicked.connect(lambda: ToolbarWindow(self))
        grid1.addWidget(editToolbarButton, 6,0,1,2)
        self.currentToolbar = self.father.defaultToolbarActionList

        grid2.addWidget(QtWidgets.QLabel("Linewidth:"), 1, 0)
        self.lwSpinBox = QtWidgets.QDoubleSpinBox()
        self.lwSpinBox.setSingleStep(0.1)
        self.lwSpinBox.setValue(self.father.defaultLinewidth)
        grid2.addWidget(self.lwSpinBox, 1, 1)
        self.color = self.father.defaultColor
        lineColorButton = QtWidgets.QPushButton("Line colour")
        lineColorButton.clicked.connect(self.setColor)
        grid2.addWidget(lineColorButton, 2, 0)
        self.xgridCheck = QtWidgets.QCheckBox("x-grid")
        self.xgridCheck.setChecked(self.father.defaultGrids[0])
        grid2.addWidget(self.xgridCheck, 3, 0, 1, 2)
        self.ygridCheck = QtWidgets.QCheckBox("y-grid")
        self.ygridCheck.setChecked(self.father.defaultGrids[1])
        grid2.addWidget(self.ygridCheck, 4, 0, 1, 2)
        grid2.addWidget(QtWidgets.QLabel("Units:"), 5, 0)
        self.unitGroup=QtWidgets.QButtonGroup()
        button=QtWidgets.QRadioButton("s/Hz")
        self.unitGroup.addButton(button, 0)
        grid2.addWidget(button, 5, 1)
        button=QtWidgets.QRadioButton("ms/kHz")
        self.unitGroup.addButton(button, 1)
        grid2.addWidget(button, 6, 1)
        button=QtWidgets.QRadioButton(u"\u03bcs/MHz")
        self.unitGroup.addButton(button, 2)
        grid2.addWidget(button, 7, 1)
        self.unitGroup.button(self.father.defaultUnits).setChecked(True)
        self.ppmCheck = QtWidgets.QCheckBox("ppm")
        self.ppmCheck.setChecked(self.father.defaultPPM)
        grid2.addWidget(self.ppmCheck, 8, 1)
        self.zeroScrollCheck = QtWidgets.QCheckBox("Scroll y-axis from zero")
        self.zeroScrollCheck.setChecked(self.father.defaultZeroScroll)
        grid2.addWidget(self.zeroScrollCheck, 9, 0, 1, 2)

        grid3.addWidget(QtWidgets.QLabel("Colourmap:"), 0, 0)
        self.cmEntry = QtWidgets.QComboBox(self)
        self.cmEntry.addItems(sc.COLORMAPLIST)
        self.cmEntry.setCurrentIndex(sc.COLORMAPLIST.index(self.father.defaultColorMap))
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
        self.father.defaultToolbarActionList  = self.currentToolbar
        self.father.defaultLinewidth = self.lwSpinBox.value()
        self.father.defaultColor = self.color
        self.father.defaultGrids[0] = self.xgridCheck.isChecked()
        self.father.defaultGrids[1] = self.ygridCheck.isChecked()
        self.father.defaultZeroScroll = self.zeroScrollCheck.isChecked()
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


class ToolbarWindow(wc.ToolWindows):

    NAME = "Change Toolbar"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self,parent):
        super(ToolbarWindow, self).__init__(parent)
        self.grid.addWidget(wc.QLabel("Actions:"), 0, 0)
        self.grid.addWidget(wc.QLabel("Toolbar Actions:"), 0, 1)
        self.listA = OrigListWidget(self)
        for i in self.father.father.allActionsList:
            QtWidgets.QListWidgetItem(i[0], self.listA).setToolTip(i[0])
        self.listB = DestListWidget(self)
        for i in self.father.father.defaultToolbarActionList:
            QtWidgets.QListWidgetItem(i, self.listB).setToolTip(i)
        self.grid.addWidget(self.listA, 1, 0)
        self.grid.addWidget(self.listB, 1, 1)
        self.resize(650, 500)

    def applyAndClose(self, *args):
        items = []
        for index in range(self.listB.count()):
            items.append(self.listB.item(index).text())
        self.father.currentToolbar = items
        self.closeEvent()
        
    def closeEvent(self, *args):
        self.deleteLater()
        
##############################################################################


class aboutWindow(wc.ToolWindows):

    NAME = "About ssNake"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent):
        super(aboutWindow, self).__init__(parent)
        self.cancelButton.hide()
        self.logo = QtWidgets.QLabel(self)
        self.logo.setPixmap(QtGui.QPixmap(os.path.dirname(os.path.realpath(__file__)) + "/logo.gif"))
        self.tabs = QtWidgets.QTabWidget(self)
        self.text = QtWidgets.QTextEdit(self)
        self.text.setReadOnly(True)
        self.license = QtWidgets.QTextEdit(self)
        self.license.setReadOnly(True)
        licenseText = ''
        with open(os.path.dirname(os.path.realpath(__file__)) + os.path.sep + 'licenseHtml.txt') as f:
            licenseTextTemp = f.read().split('\n')
        licenseText = ' '.join(licenseTextTemp)
        self.license.setHtml(licenseText)
        pythonVersion = sys.version
        pythonVersion = pythonVersion[:pythonVersion.index(' ')]
        try:
            from PyQt4.Qt import PYQT_VERSION_STR
            from PyQt4.QtCore import QT_VERSION_STR
        except:
            from PyQt5.Qt import PYQT_VERSION_STR
            from PyQt5.QtCore import QT_VERSION_STR
        from scipy import __version__ as scipyVersion    
        self.text.setText('<p><b>ssNake ' + VERSION + '</b></p>' + 
                '<p>Copyright (&copy;) 2016&ndash;2017 Bas van Meerten & Wouter Franssen<\p>' + '<p>Email: <a href="mailto:ssnake@science.ru.nl" >ssnake@science.ru.nl</a></p>' +
        '<b>Library versions</b>:<br>Python ' + pythonVersion + '<br>numpy ' + np.__version__ +
        '<br>SciPy ' + scipyVersion +  
        '<br>matplotlib ' + matplotlib.__version__ + 
        '<br>PyQt ' + PYQT_VERSION_STR + 
        '<br>Qt ' + QT_VERSION_STR )
        self.thanks = QtWidgets.QTextEdit(self)
        self.thanks.setReadOnly(True)
        self.thanks.setHtml('<p><b>The ssNake team wishes to thank:</b></p>Koen Tijssen<br>Ole Brauckmann')
        self.tabs.addTab(self.text, 'Version') 
        self.tabs.addTab(self.thanks, 'Thanks') 
        self.tabs.addTab(self.license, 'License') 
        self.grid.addWidget(self.logo, 0, 0, 1, 3, QtCore.Qt.AlignHCenter)
        self.grid.addWidget(self.tabs, 1, 0, 1, 3)
        self.resize(550, 700)

    def closeEvent(self, *args):
        self.deleteLater()

##############################################################################


class shiftConversionWindow(wc.ToolWindows):

    NAME = "Chemical Shift Conversions"
    MENUDISABLE = False
    RESIZABLE = True

    def __init__(self, parent):
        super(shiftConversionWindow, self).__init__(parent)

        self.standardGroup = QtWidgets.QGroupBox('Standard Convention:')
        self.standardFrame = QtWidgets.QGridLayout()

        D11label = wc.QLabel(u'\u03b4' + '<sub>11</sub> [ppm]')
        self.standardFrame.addWidget(D11label, 0, 1)
        D22label = wc.QLabel(u'\u03b4' + '<sub>22</sub> [ppm]')
        self.standardFrame.addWidget(D22label, 0, 2)
        D33label = wc.QLabel(u'\u03b4' + '<sub>33</sub> [ppm]')
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
        self.grid.addWidget(self.standardGroup,0,0,1,3)



        # xyz Convention

        self.xyzGroup = QtWidgets.QGroupBox('xyz Convention:')
        self.xyzFrame = QtWidgets.QGridLayout()
        dxxlabel = wc.QLabel(u'\u03b4' + '<sub>xx</sub> [ppm]')
        self.xyzFrame.addWidget(dxxlabel, 3, 1)
        dyylabel = wc.QLabel(u'\u03b4' + '<sub>yy</sub> [ppm]')
        self.xyzFrame.addWidget(dyylabel, 3, 2)
        dzzlabel = wc.QLabel(u'\u03b4' + '<sub>zz</sub> [ppm]')
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
        self.grid.addWidget(self.xyzGroup,1,0,1,3)


        # Haeberlen Convention
        self.haebGroup = QtWidgets.QGroupBox('Haeberlen Convention')
        self.haebFrame = QtWidgets.QGridLayout()
        disolabel = wc.QLabel(u'\u03b4' + '<sub>iso</sub> [ppm]')
        self.haebFrame.addWidget(disolabel, 6, 1)
        danisolabel = wc.QLabel(u'\u03b4' + '<sub>aniso</sub> [ppm]')
        self.haebFrame.addWidget(danisolabel, 6, 2)
        etalabel = wc.QLabel(u'\u03b7')
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
        self.grid.addWidget(self.haebGroup,2,0,1,3)

        # Hertzfeld berger
        self.hbGroup = QtWidgets.QGroupBox('Hertzfeld-Berger Convention')
        self.hbFrame = QtWidgets.QGridLayout()
        hbdisolabel = wc.QLabel(u'\u03b4' + '<sub>iso</sub> [ppm]')
        self.hbFrame.addWidget(hbdisolabel, 9, 1)
        omegalabel = wc.QLabel(u'\u03a9 [ppm]')
        self.hbFrame.addWidget(omegalabel, 9, 2)
        skewlabel = wc.QLabel(u'\u03ba')
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
        self.grid.addWidget(self.hbGroup,3,0,1,3)

        # Reset
        self.cancelButton.setText("Reset")
        self.cancelButton.clicked.disconnect()
        self.cancelButton.clicked.connect(self.valueReset)
        self.okButton.setText("Close")
        self.okButton.clicked.disconnect()
        self.okButton.clicked.connect(self.closeEvent)

    def shiftCalc(self, Type):
        if Type == 0:  # If from standard
            try:
                delta11 = float(safeEval(self.D11.text()))
                delta22 = float(safeEval(self.D22.text()))
                delta33 = float(safeEval(self.D33.text()))
                Values = [delta11,delta22,delta33]
            except:
                self.father.dispMsg("Shift Conversion: Invalid input in Standard Convention")
                return
        if Type == 1:  # If from xyz
            try:
                delta11 = float(safeEval(self.dxx.text()))  # Treat xyz as 123, as it reorders them anyway
                delta22 = float(safeEval(self.dyy.text()))
                delta33 = float(safeEval(self.dzz.text()))
                Values = [delta11,delta22,delta33]
            except:
                self.father.dispMsg("Shift Conversion: Invalid input in xyz Convention")
                return
        if Type == 2:  # From haeberlen
            try:
                eta = float(safeEval(self.eta.text()))
                delta = float(safeEval(self.daniso.text()))
                iso = float(safeEval(self.diso.text()))
                Values = [iso,delta,eta]
            except:
                self.father.dispMsg("Shift Conversion: Invalid input in Haeberlen Convention")
                return                
        if Type == 3:  # From Hertzfeld-Berger
            try:
                iso = float(safeEval(self.hbdiso.text()))
                span = float(safeEval(self.hbdaniso.text()))
                skew = float(safeEval(self.hbskew.text()))
                Values = [iso,span,skew]
            except:
                self.father.dispMsg("Shift Conversion: Invalid input in Hertzfeld-Berger Convention")
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

##############################################################################


class quadConversionWindow(wc.ToolWindows):
    
    Ioptions = ['1', '3/2', '2', '5/2', '3', '7/2', '4', '9/2','5','6','7']
    Ivalues = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,5.0,6.0,7.0]
    
    NAME = "Quadrupolar Coupling Conversions"
    RESIZABLE = True
    MENUDISABLE = False

    def __init__(self, parent):
        super(quadConversionWindow, self).__init__(parent)

        self.comGroup = QtWidgets.QGroupBox("Common Parameters:")
        self.comFrame = QtWidgets.QGridLayout()
        Itext = wc.QLabel("I:")
        self.comFrame.addWidget(Itext, 0, 0)
        self.IEntry = QtWidgets.QComboBox()
        self.IEntry.addItems(self.Ioptions)
        self.IEntry.setCurrentIndex(0)
        self.comFrame.addWidget(self.IEntry, 1, 0)
        etalabel = wc.QLabel(u'\u03b7:')
        self.comFrame.addWidget(etalabel, 0, 1)
        self.Eta = wc.QLineEdit("0")
        self.Eta.setMinimumWidth(100)
        self.comFrame.addWidget(self.Eta, 1, 1)
        momentlabel = wc.QLabel('Q [fm<sup>2</sup>]:')
        self.comFrame.addWidget(momentlabel, 0, 2)
        self.Moment = wc.QLineEdit("ND")
        self.Moment.setMinimumWidth(100)
        self.comFrame.addWidget(self.Moment, 1, 2)
        self.comGroup.setLayout(self.comFrame)
        self.grid.addWidget(self.comGroup,1,0,1,3)

        self.CqGroup = QtWidgets.QGroupBox("C_Q Convention:")
        self.CqFrame = QtWidgets.QGridLayout()
        Cqlabel = wc.QLabel(u'C' + u'<sub>Q</sub>/2\u03c0 [MHz:]')
        self.CqFrame.addWidget(Cqlabel, 3, 1)
        CqGO = QtWidgets.QPushButton("Go")
        self.CqFrame.addWidget(CqGO, 4, 0)
        CqGO.clicked.connect(lambda: self.quadCalc(0))
        self.Cq = wc.QLineEdit("0")
        self.Cq.setMinimumWidth(100)
        self.CqFrame.addWidget(self.Cq, 4, 1)
        self.CqGroup.setLayout(self.CqFrame)
        self.grid.addWidget(self.CqGroup,3,0,1,2)

        self.WqGroup = QtWidgets.QGroupBox(u"\u03c9_Q Convention:")
        self.WqFrame = QtWidgets.QGridLayout()
        Wqlabel = wc.QLabel(u'\u03c9' + u'<sub>Q</sub>/2\u03c0 [MHz]:')
        self.WqFrame.addWidget(Wqlabel, 6, 1)
        WqGO = QtWidgets.QPushButton("Go")
        self.WqFrame.addWidget(WqGO, 7, 0)
        WqGO.clicked.connect(lambda: self.quadCalc(1))
        self.Wq = wc.QLineEdit("0")
        self.Wq.setMinimumWidth(100)
        self.WqFrame.addWidget(self.Wq, 7, 1)        
        self.WqGroup.setLayout(self.WqFrame)
        self.grid.addWidget(self.WqGroup,7,0,1,2)

        self.fieldGroup = QtWidgets.QGroupBox('Field Gradients:')
        self.fieldFrame = QtWidgets.QGridLayout()
        Vxxlabel = wc.QLabel('V<sub>xx</sub> [V/m<sup>2</sup>]:')
        self.fieldFrame.addWidget(Vxxlabel, 9, 1)
        VGO = QtWidgets.QPushButton("Go")
        self.fieldFrame.addWidget(VGO, 10, 0)
        VGO.clicked.connect(lambda: self.quadCalc(2))
        self.Vxx = wc.QLineEdit("ND")
        self.Vxx.setMinimumWidth(100)
        self.fieldFrame.addWidget(self.Vxx, 10, 1)        
        Vyylabel = wc.QLabel('V<sub>yy</sub> [V/m<sup>2</sup>]:')
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
        self.grid.addWidget(self.fieldGroup,8,0,1,4)

        # Reset
        self.cancelButton.setText("Reset")
        self.cancelButton.clicked.disconnect()
        self.cancelButton.clicked.connect(self.valueReset)
        self.okButton.setText("Close")
        self.okButton.clicked.disconnect()
        self.okButton.clicked.connect(self.closeEvent)

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
                self.father.dispMsg("Quad Conversion: Invalid input in Cq definition")
                return 
        if Type == 1:
            try:
                Vmax = float(safeEval(self.Wq.text()))                 
                Eta = float(safeEval(self.Eta.text())) 
                Czz = Vmax*(2.0*I*(2*I-1))/3.0
                Cxx = Czz*(Eta-1)/2
                Cyy = -Cxx-Czz
            except:
                self.father.dispMsg("Quad Conversion: Invalid input in Wq definition")
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
                self.father.dispMsg("Quad Conversion: Invalid input in field gradients")
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
    
    sys._excepthook = sys.excepthook
    def exception_hook(exctype, value, traceback):
        sys._excepthook(exctype, value, traceback)
        mainProgram.dispError([exctype, value, traceback])
    sys.excepthook = exception_hook
    sys.exit(root.exec_())
    
