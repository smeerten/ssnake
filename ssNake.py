#!/usr/bin/env python

# Copyright 2015 Bas van Meerten and Wouter Franssen

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

import sip
sip.setapi('QString', 2)

from PyQt4 import QtGui, QtCore
from numpy import arange, sin, pi
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import numpy as np
import sys
import os
import scipy.io
import json
import copy
import math
from struct import unpack

import spectrum_classes as sc
import fitting as fit
from safeEval import *

pi=math.pi

class MainProgram(QtGui.QMainWindow):
    def __init__(self,root):
        QtGui.QMainWindow.__init__(self)
        self.root = root
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setAcceptDrops(True)
        self.mainWindow = None
        self.workspaces = []
        self.workspaceNames = []
        self.workspaceNum = 0
        self.macros = {}
        self.macroActions = {} 
        
        self.menubar = self.menuBar()
        self.filemenu = QtGui.QMenu('File', self)
        self.menubar.addMenu(self.filemenu)
        self.filemenu.addAction('Open', self.loadFromMenu,QtGui.QKeySequence.Open)
        self.savemenu = QtGui.QMenu('Save',self)
        self.filemenu.addMenu(self.savemenu)
        self.savefigAct = self.savemenu.addAction('Save figure', self.saveFigure,QtGui.QKeySequence.Print)
        self.savemenu.addAction('Save JSON', self.saveJSONFile,QtGui.QKeySequence.Save)
        self.savemenu.addAction('Save MATLAB', self.saveMatlabFile)
        self.savemenu.addAction('Save as Simpson data', self.saveSimpsonFile)
        self.filemenu.addAction('Quit', self.fileQuit, QtGui.QKeySequence.Quit)
        
        self.workspacemenu = QtGui.QMenu('Workspaces',self)
        self.menubar.addMenu(self.workspacemenu)
        self.workspacemenu.addAction('Duplicate', self.duplicateWorkspace,QtGui.QKeySequence.New)
        self.workspacemenu.addAction('Delete', self.destroyWorkspace,QtGui.QKeySequence.Close)
        self.workspacemenu.addAction('Rename', self.renameWorkspace)
        self.activemenu = QtGui.QMenu('Active',self)
        self.workspacemenu.addMenu(self.activemenu)
        self.workspacemenu.addAction('Next', lambda: self.stepWorkspace(1),QtGui.QKeySequence.Forward)
        self.workspacemenu.addAction('Previous', lambda: self.stepWorkspace(-1),QtGui.QKeySequence.Back)
        self.macromenu = QtGui.QMenu('Macros',self)
        self.menubar.addMenu(self.macromenu)
        self.macrostartAct = self.macromenu.addAction('Start recording', self.macroCreate)
        self.macrostopAct = self.macromenu.addAction('Stop recording', self.stopMacro)
        self.macrolistmenu = QtGui.QMenu('Run',self)
        self.macromenu.addMenu(self.macrolistmenu)
        self.macrorenamemenu = QtGui.QMenu('Rename',self)
        self.macromenu.addMenu(self.macrorenamemenu)
        self.macrodeletemenu = QtGui.QMenu('Delete',self)
        self.macromenu.addMenu(self.macrodeletemenu)
        self.macrosavemenu = QtGui.QMenu('Save',self)
        self.macromenu.addMenu(self.macrosavemenu)
        self.macromenu.addAction('Load', self.loadMacro)

        self.menuCheck()
        
        self.main_widget = QtGui.QWidget(self)
        self.mainFrame = QtGui.QGridLayout(self.main_widget)

        self.logo = QtGui.QLabel(self)
        self.logo.setPixmap(QtGui.QPixmap(os.path.dirname(os.path.realpath(__file__)) + "/logo.gif"))
        self.mainFrame.addWidget(self.logo,0,0,QtCore.Qt.AlignCenter)
        
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        self.eventFilter = MyEventFilter(self)
        self.root.installEventFilter(self.eventFilter)
        
    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            self.autoLoad(path)
   
    def menuCheck(self):
        if self.mainWindow is None:
            self.savemenu.menuAction().setEnabled(False)
            self.workspacemenu.menuAction().setEnabled(False)
            self.macromenu.menuAction().setEnabled(False)
        elif isinstance(self.mainWindow, Main1DWindow):
            self.savemenu.menuAction().setEnabled(True)
            self.savefigAct.setEnabled(True)
            self.macromenu.menuAction().setEnabled(True)
            if self.mainWindow.currentMacro is None:
                self.macrostopAct.setEnabled(False)
                self.macrostartAct.setEnabled(True)
            else:
                self.macrostopAct.setEnabled(True)
                self.macrostartAct.setEnabled(False)
            self.savemenu.menuAction().setEnabled(True)
            self.workspacemenu.menuAction().setEnabled(True)
        elif isinstance(self.mainWindow, fit.MainPlotWindow):
            self.savemenu.menuAction().setEnabled(True)
            self.savefigAct.setEnabled(False)
            self.workspacemenu.menuAction().setEnabled(True)
            self.macromenu.menuAction().setEnabled(False)
        else:
            self.savemenu.menuAction().setEnabled(True)
            self.savefigAct.setEnabled(True)
            self.workspacemenu.menuAction().setEnabled(True)
            self.macromenu.menuAction().setEnabled(False)

    def askName(self,filePath=None):
        if filePath is None:
            message = 'Spectrum name'
        else:
            message = 'Spectrum name for: ' + filePath
        count = 0
        name = 'spectrum'+str(count)
        while name in self.workspaceNames:
            count += 1
            name = 'spectrum'+str(count)
        givenName, ok = QtGui.QInputDialog.getText(self, message, 'Name:',text=name)
        if not ok:
            return
        while (givenName in self.workspaceNames) or givenName=='':
            print('Name exists')
            givenName, ok = QtGui.QInputDialog.getText(self, message, 'Name:',text=name)
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
        name = 'macro'+str(count)
        while name in self.macros.keys():
            count += 1
            name = 'macro'+str(count)
        givenName, ok = QtGui.QInputDialog.getText(self, 'Macro name', 'Name:',text=name)
        while (givenName in self.macros.keys()) or givenName is '':
            print('Name exists')
            givenName, ok = QtGui.QInputDialog.getText(self, 'Macro name', 'Name:',text=name)
        self.macros[givenName] = []
        self.mainWindow.redoMacro = []
        self.mainWindow.currentMacro = givenName
        action1 = self.macrolistmenu.addAction(givenName,lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(givenName,lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(givenName,lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(givenName,lambda name=givenName: self.renameMacro(name))
        self.macroActions[givenName] = [action1,action2,action3,action4]
        self.menuCheck()

    def renameMacro(self,oldName):
        if self.mainWindow is None:
            return
        count = 0
        name = 'macro'+str(count)
        while name in self.macros.keys():
            count += 1
            name = 'macro'+str(count)
        givenName, ok = QtGui.QInputDialog.getText(self, 'Macro name', 'Name:',text=name)
        while (givenName in self.macros.keys()) or givenName is '':
            print('Name exists')
            givenName, ok = QtGui.QInputDialog.getText(self, 'Macro name', 'Name:',text=name)
        self.macros[givenName] = self.macros.pop(oldName)
        if self.mainWindow.currentMacro == oldName:
            self.mainWindow.currentMacro = givenName
        oldActions = self.macroActions.pop(oldName)
        self.macrolistmenu.removeAction(oldActions[0])
        self.macrosavemenu.removeAction(oldActions[1])
        self.macrodeletemenu.removeAction(oldActions[2])
        self.macrorenamemenu.removeAction(oldActions[3])
        action1 = self.macrolistmenu.addAction(givenName,lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.addAction(givenName,lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.addAction(givenName,lambda name=givenName: self.deleteMacro(name))
        action4 = self.macrorenamemenu.addAction(givenName,lambda name=givenName: self.renameMacro(name))
        self.macroActions[givenName] = [action1,action2,action3,action4]
        self.menuCheck()
        
    def stopMacro(self):
        if self.mainWindow is None:
            return
        if self.mainWindow.currentMacro is None:
            return
        self.mainWindow.redoMacro = []
        self.mainWindow.currentMacro = None
        self.menuCheck()

    def macroAdd(self,name,macros):
        self.macros[name].append(macros)

    def runMacro(self,name):
        if self.mainWindow is not None:
            self.mainWindow.runMacro(self.macros[name])

    def saveMacro(self,name):
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save File')
        if not fileName:
            return
        with open(fileName,'w') as f:
            json.dump(self.macros[name],f)

    def deleteMacro(self,name):
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
        filename = QtGui.QFileDialog.getOpenFileName(self, 'Open File')
        if len(filePath)==0:
            return
        count = 0
        name = 'macro'+str(count)
        while name in self.macros.keys():
            count += 1
            name = 'macro'+str(count)
        givenName, ok = QtGui.QInputDialog.getText(self, 'Macro name', 'Name:',text=name)
        while (givenName in self.macros.keys()) or givenName is '':
            print('Name exists')
            givenName, ok = QtGui.QInputDialog.getText(self, 'Macro name', 'Name:',text=name)
        with open(filePath,'r') as f:
            self.macros[givenName] = json.load(f)
        action1 = self.macrolistmenu.add_command(label=givenName,command=lambda name=givenName: self.runMacro(name))
        action2 = self.macrosavemenu.add_command(label=givenName,command=lambda name=givenName: self.saveMacro(name))
        action3 = self.macrodeletemenu.add_command(label=givenName,command=lambda name=givenName: self.deleteMacro(name))
        self.macroActions[givenName] = [action1,action2,action3]
        self.menuCheck()

    def changeMainWindow(self, var):
        self.logo.hide()
        if self.mainWindow is not None:
            self.mainWindow.removeFromView()
        num = self.workspaceNames.index(var)
        self.workspaceNum = num
        self.mainWindow = self.workspaces[num]
        self.mainWindow.addToView()
        self.updWorkspaceMenu(var)
        self.menuCheck()
        
    def stepWorkspace(self, step):
        if len(self.workspaces) > 1:
            if self.mainWindow is not None:
                self.mainWindow.removeFromView()
            self.workspaceNum += step
            self.workspaceNum = self.workspaceNum % len(self.workspaces)
            self.mainWindow = self.workspaces[self.workspaceNum]
            self.mainWindow.addToView()
            self.updWorkspaceMenu(self.workspaceNames[self.workspaceNum])
            self.menuCheck()

    def duplicateWorkspace(self, *args):
        name = self.askName()
        if name is None:
            return
        self.workspaces.append(Main1DWindow(self,copy.deepcopy(self.mainWindow.get_masterData()),self.mainWindow.get_current(),name=name))
        self.mainFrame.addWidget(self.workspaces[-1])
        self.mainWindow.removeFromView()
        self.workspaces[-1].initUI()
        self.workspaceNames.append(name)
        self.changeMainWindow(name)

    def renameWorkspace(self, *args):
        name = self.askName()
        if name is None:
            return
        self.workspaceNames[self.workspaceNum] = name
        self.updWorkspaceMenu(name)
        self.workspaces[self.workspaceNum].rename(name)

    def destroyWorkspace(self, *args):
        self.mainWindow.removeFromView()
        self.mainWindow.kill()
        self.mainWindow = None
        del self.workspaceNames[self.workspaceNum]
        del self.workspaces[self.workspaceNum]
        if self.workspaceNum == len(self.workspaces):
            self.workspaceNum = 0
        if len(self.workspaces) > 0:
            self.changeMainWindow(self.workspaceNames[self.workspaceNum])
        else:
            self.logo.show()
            self.updWorkspaceMenu(None)

    def updWorkspaceMenu(self,var):
        self.activemenu.clear()
        for i in self.workspaceNames:
            self.activemenu.addAction(i, lambda i=i: self.changeMainWindow(i))
        self.menuCheck()

    def loadFromMenu(self):
        filePath = str(QtGui.QFileDialog.getOpenFileName(self,'Open File'))
        if len(filePath)==0:
            return
        self.autoLoad(filePath)
        
    def autoLoad(self,filePath):
        if os.path.isfile(filePath):
            filename = os.path.basename(filePath)
            if filename.endswith('.fid') or filename.endswith('.spe'): 
                self.loading(4,filePath)
                return
            elif filename.endswith('.json') or filename.endswith('.JSON'):
                self.loading(5,filePath)
                return
            elif filename.endswith('.mat') or filename.endswith('.MAT'):
                self.loading(6,filePath)
                return
            filePath = os.path.dirname(filePath)
        direc = filePath
        if os.path.exists(direc+os.path.sep+'procpar') and os.path.exists(direc+os.path.sep+'fid'):
            self.loading(0,filePath)
            return
        elif os.path.exists(direc+os.path.sep+'acqus') and (os.path.exists(direc+os.path.sep+'fid') or os.path.exists(direc+os.path.sep+'ser')):
            self.loading(1,filePath)
            return
        elif os.path.exists(direc+os.path.sep+'acq') and os.path.exists(direc+os.path.sep+'data'):
            self.loading(2,filePath)
            return
        elif os.path.exists(direc+os.path.sep+'acqu.par'):
            dirFiles = os.listdir(direc)
            files2D = [x for x in dirFiles if '.2d' in x]
            files1D = [x for x in dirFiles if '.1d' in x]
            if len(files2D)!=0 or len(files1D)!=0:
                self.loading(3,filePath)
                return

    def dataFromFit(self, data, freq , sw , spec, wholeEcho, ref, xaxArray,axes):
        name = self.askName()
        if name is None:
            return
        masterData=sc.Spectrum(data,lambda self :self.dataFromFit(data, freq , sw , spec, wholeEcho, ref, xaxArray), freq , sw , spec, wholeEcho, ref, xaxArray)
        masterData.resetXax(axes)
        self.workspaces.append(Main1DWindow(self,masterData,name=name))
        self.mainFrame.addWidget(self.workspaces[-1])
        self.workspaces[-1].initUI()
        self.workspaceNames.append(name)
        self.changeMainWindow(name)

    def loading(self,num,filePath):
        name = self.askName(filePath)
        if name is None:
            return
        if num==0:
            masterData = self.LoadVarianFile(filePath)
        elif num==1:
            masterData = self.LoadBrukerTopspin(filePath)
        elif num==2:
            masterData = self.LoadChemFile(filePath)
        elif num==3:
            masterData = self.LoadMagritek(filePath)
        elif num==4:
            masterData = self.LoadSimpsonFile(filePath)
        elif num==5:
            masterData = self.loadJSONFile(filePath)
        elif num==6:
            masterData = self.loadMatlabFile(filePath)
        if masterData is not None:
            self.workspaces.append(Main1DWindow(self,masterData,name=name))
            self.mainFrame.addWidget(self.workspaces[-1])
            self.workspaces[-1].initUI()
            self.workspaceNames.append(name)
            self.changeMainWindow(name)

    def LoadVarianFile(self,filePath):
        Dir = filePath 
        freq = 300e6
        sw   = 50e3
        sw1  = 50e3
        if os.path.exists(Dir+os.path.sep+'procpar'):
            with open(Dir+os.path.sep+'procpar', 'r') as f: 
                data = f.read().split('\n')
            for s in range(0,len(data)): 
                if data[s].startswith('sfrq '):
                    freq=float(data[s+1].split()[1])*1e6 
                elif data[s].startswith('sw '):
                    sw=float(data[s+1].split()[1])
                elif data[s].startswith('sw1 '):
                    sw1=float(data[s+1].split()[1])
        else:
            print(Dir+os.path.sep+'procpar does not exits, used standard sw and freq')
        if os.path.exists(Dir+os.path.sep+'fid'):    
            with open(Dir+os.path.sep+'fid', "rb") as f:
                raw = np.fromfile(f, np.int32,6) 
                nblocks = unpack('>l', raw[0])[0]
                ntraces = unpack('>l', raw[1])[0]
                npoints = unpack('>l', raw[2])[0]
                ebytes = unpack('>l', raw[3])[0]
                tbytes = unpack('>l', raw[4])[0]
                bbytes = unpack('>l', raw[5])[0]
                raw = np.fromfile(f, np.int16,2)
                vers_id = unpack('>h', raw[0])[0] 
                status = unpack('>h', raw[1])[0]
                raw = np.fromfile(f, np.int32,1) 
                nbheaders = unpack('>l', raw[0])[0]
                SizeTD2 = npoints
                SizeTD1 = nblocks*ntraces
                a = []
                fid32 = int(bin(status)[-3]) 
                fidfloat = int(bin(status)[-4])
                for iter1 in range(0,nblocks): 
                    b = []
                    for iter2 in range(0,nbheaders):
                        raw = np.fromfile(f, np.int16,nbheaders*14)
                    if not fid32 and not fidfloat:
                        raw = np.fromfile(f, np.int16,ntraces*npoints)
                        for iter3 in raw:
                            b.append(unpack('>h', iter3)[0])
                    elif fid32 and not fidfloat:
                        raw = np.fromfile(f, np.int32,ntraces*npoints)
                        for iter3 in raw:
                            b.append(unpack('>l', iter3)[0])
                    else:
                        raw = np.fromfile(f, np.float32,ntraces*npoints)
                        for iter3 in raw:
                            b.append(unpack('>f', iter3)[0])
                    b=np.array(b)
                    if(len(b) != ntraces*npoints):
                        b.append(np.zeros(ntraces*npoints-len(b)))
                    a.append(b)
        a=np.complex128(a)
        fid = a[:,::2]-1j*a[:,1::2]
        if SizeTD1 is 1: 
            fid = fid[0][:]
            masterData=sc.Spectrum(fid,lambda self :self.LoadVarianFile(filePath),[freq],[sw])
        else: 
            masterData=sc.Spectrum(fid,lambda self :self.LoadVarianFile(filePath),[freq]*2,[sw]*2)
        return masterData

    def loadJSONFile(self,filePath):
        with open(filePath, 'r') as inputfile:
            struct = json.load(inputfile)
        data = np.array(struct['dataReal']) + 1j * np.array(struct['dataImag'])
        ref = np.where(np.isnan(struct['ref']), None, struct['ref'])
        xaxA = []
        for i in struct['xaxArray']:
            xaxA.append(np.array(i))
        masterData=sc.Spectrum(data,lambda self :self.loadJSONFile(filePath),list(struct['freq']),list(struct['sw']),list(struct['spec']),list(np.array(struct['wholeEcho'],dtype=bool)),list(ref),xaxA)
        return masterData
        
    def loadMatlabFile(self,filePath):
        matlabStruct = scipy.io.loadmat(filePath)
        var = [k for k in matlabStruct.keys() if not k.startswith('__')][0]
        mat = matlabStruct[var]
        if mat['dim']==1:
            xaxA = [k[0] for k in (mat['xaxArray'][0])]
            data = mat['data'][0,0][0]
        else:
            xaxA = [k[0] for k in (mat['xaxArray'][0,0][0])]
            data = mat['data'][0,0]
        #insert some checks for data type
        ref = mat['ref'][0,0][0]
        ref = np.where(np.isnan(ref), None, ref)
        masterData=sc.Spectrum(data,lambda self :self.loadMatlabFile(filePath),list(mat['freq'][0,0][0]),list(mat['sw'][0,0][0]),list(mat['spec'][0,0][0]),list(np.array(mat['wholeEcho'][0,0][0])>0),list(ref),xaxA)
        return masterData
        
    def LoadBrukerTopspin(self,filePath):
        Dir = filePath 
        if os.path.exists(Dir+os.path.sep+'acqus'):
            with open(Dir+os.path.sep+'acqus', 'r') as f: 
                data = f.read().split('\n')
            for s in range(0,len(data)):
                if data[s].startswith('##$TD='):
                    sizeTD2 = int(data[s][6:])
                if data[s].startswith('##$SFO1='):
                    freq2 = float(data[s][8:])*1e6
                if data[s].startswith('##$SW_h='):
                    SW2 = float(data[s][8:])
                if data[s].startswith('##$BYTORDA='):
                    ByteOrder = int(data[s][11:]) 
        sizeTD1=1 
        if os.path.exists(Dir+os.path.sep+'acqu2s'): 
            with open(Dir+os.path.sep+'acqu2s', 'r') as f: 
                data2 = f.read().split('\n')
            for s in range(0,len(data2)):
                if data2[s].startswith('##$TD='):
                    sizeTD1 = int(data2[s][6:])
                if data2[s].startswith('##$SFO1='):
                    freq1 = float(data2[s][8:])*1e6
                if data2[s].startswith('##$SW_h='):
                    SW1 = float(data2[s][8:])
        if os.path.exists(Dir+os.path.sep+'fid'):
            with open(Dir+os.path.sep+'fid', "rb") as f:            
                raw = np.fromfile(f, np.int32,sizeTD1*sizeTD2)
        elif os.path.exists(Dir+os.path.sep+'ser'):
            with open(Dir+os.path.sep+'ser', "rb") as f:            
                raw = np.fromfile(f, np.int32,sizeTD1*sizeTD2)
        if ByteOrder: 
            RawInt=raw.newbyteorder('b')
        else:
            RawInt=raw.newbyteorder('l')
        ComplexData = np.array(RawInt[0:len(RawInt):2])+1j*np.array(RawInt[1:len(RawInt):2])
        spec = [False]
        if sizeTD1 is 1:
            masterData=sc.Spectrum(ComplexData,lambda self :self.LoadBrukerTopspin(filePath),[freq2],[SW2],spec)
        else:
            data = ComplexData.reshape(sizeTD1,sizeTD2/2)
            masterData=sc.Spectrum(data,lambda self :self.LoadBrukerTopspin(filePath),[freq1,freq2],[SW1,SW2],spec*2)
        return masterData
                
    def LoadChemFile(self,filePath):
        Dir = filePath
        sizeTD1=1
        sw1=50e3
        H = dict(line.strip().split('=') for line in open(Dir+os.path.sep+'acq','r'))
        sizeTD2 = int(H['al'])
        freq = float(H['sf'+H['ch1']])
        sw=1/float(H['dw'][:-1])
        if any('array_num_values_' in s for s in H.keys()):
            if 'use_array=1' in open(Dir+'/acq_2').read():
                for s in H.keys():
                    if ('array_num_values_' in s):
                        sizeTD1 = sizeTD1*int(H[s])
            else:
                if 'al2' in H:
                    sizeTD1 = int(H['al2'])
                    if 'dw2' in H:
                        sw1 = 1/float(H['dw2'][:-1])
        else:
            if 'al2' in H:
                sizeTD1 = int(H['al2'])
                if 'dw2' in H:
                    sw1 = 1/float(H['dw2'][:-1])        
        with open(Dir+os.path.sep+'data','rb') as f:
            raw = np.fromfile(f, np.int32)
            b=np.complex128(raw.byteswap())
        fid = b[:len(b)/2]+1j*b[len(b)/2:]
        fid = np.reshape(fid,(sizeTD1,sizeTD2))
        data = np.array(fid) 
        spec = [False]                    
        if sizeTD1 is 1:
            data = data[0][:]
            masterData=sc.Spectrum(data,lambda self :self.LoadChemFile(filePath),[freq*1e6],[sw],spec)
        else:
            data = data.reshape((sizeTD1,sizeTD2))
            masterData=sc.Spectrum(data,lambda self :self.LoadChemFile(filePath),[freq*1e6]*2,[sw1,sw],spec*2)
        return masterData

    def LoadMagritek(self,filePath):
        #Magritek load script based on some Matlab files by Ole Brauckman
        Dir = filePath
        DirFiles = os.listdir(Dir)
        Files2D = [x for x in DirFiles if '.2d' in x]
        Files1D = [x for x in DirFiles if '.1d' in x]
        H = dict(line.strip().split('=') for line in open(Dir+os.path.sep+'acqu.par','r'))
        sw = float(H['bandwidth                 '])*1000 
        sizeTD2 = int(H['nrPnts                    '])
        freq = float(H['b1Freq                    '])
        if len(Files2D)==1:
            File=Files2D[0]
            sizeTD1 = int(H['nrSteps                   '])
            if 'bandwidth2                ' in H:
                sw1 = float(H['bandwidth2                ']) 
            else:
                sw1 = 50e3 
            with open(Dir+os.path.sep+File,'rb') as f:
                raw = np.fromfile(f, np.float32)
            Data = raw[-2*sizeTD2*sizeTD1::]
            ComplexData = Data[0:Data.shape[0]:2]-1j*Data[1:Data.shape[0]:2]
            ComplexData = ComplexData.reshape((sizeTD1,sizeTD2))
            masterData=sc.Spectrum(ComplexData,lambda self :self.LoadMagritek(filePath),[freq*1e6]*2,[sw,sw1],[False]*2)
        elif len(Files1D)!=0:
            File = 'data.1d'
            with open(Dir+os.path.sep+File,'rb') as f:
                raw = np.fromfile(f, np.float32)
            Data = raw[-2*sizeTD2::]
            ComplexData = Data[0:Data.shape[0]:2]-1j*Data[1:Data.shape[0]:2]
            masterData=sc.Spectrum(ComplexData,lambda self :self.LoadMagritek(filePath),[freq*1e6],[sw],[False])
        return masterData
            
    def LoadSimpsonFile(self,filePath):
        with open(filePath, 'r') as f: 
            Lines = f.read().split('\n')
        NP, NI, SW, SW1, TYPE, FORMAT = 0,1,0,0,'','Normal'
        DataStart = Lines.index('DATA')
        DataEnd = Lines.index('END')
        for s in range(0,DataStart):
            if Lines[s].startswith('NP='):
                NP = int(re.sub('NP=','',Lines[s]))
            elif Lines[s].startswith('NI='):
                NI = int(re.sub('NI=','',Lines[s]))
            elif Lines[s].startswith('SW='):
                SW = float(re.sub('SW=','',Lines[s]))
            elif Lines[s].startswith('SW1='):
                SW1 = float(re.sub('SW1=','',Lines[s]))
            elif Lines[s].startswith('TYPE='):
                TYPE = re.sub('TYPE=','',Lines[s])
            elif Lines[s].startswith('FORMAT='):
                FORMAT = re.sub('FORMAT=','',Lines[s])
        if 'Normal' in FORMAT: 
            data = []
            for iii in range(DataStart+1,DataEnd):
                temp = Lines[iii].split()
                data.append(float(temp[0])+1j*float(temp[1]))
        elif 'BINARY' in FORMAT: 
            RawData = np.array(Lines[DataStart+1:DataEnd])
            Ascii=[]
            for i in range(0,len(RawData)):
                for j in range(0,len(RawData[i])):
                    Ascii.append(ord(RawData[i][j]))
            Values = np.array(Ascii)
            i=0
            idx=0
            TempData=[]
            while i < 2*NP*NI:
                pts=[]
                for j in range(0,4):
                    C=[]
                    for h in range(4):
                        try: 
                            C.append(Values[idx]-33)
                        except: 
                            C.append(0)
                        idx+=1
                    pts.append(C[0]%64 + C[1]*4- C[1]*4 % 64)
                    pts.append(C[1]%16 + C[2]*4- C[2]*4%16)
                    pts.append(C[2]%4 + C[3]*4- C[3]*4%4)
                for k in range(0,3):
                    p=0
                    if i < 2*NP*NI:
                        for j in range(0,4):
                            p = np.int32(p*256);
                            p = np.int32(p | pts[4*k+j])
                        a1 = np.int32(math.floor(p)%256 * 16777216)
                        a2 = np.int32(math.floor(p/256)%256 * 65536)
                        a3 = np.int32(math.floor(p/65536)%256 * 256)
                        a4 = np.int32(math.floor(p/16777216)%256)
                        rdl = a1 | a2 | a3 | a4
                        sign = math.floor(rdl/2**31)
                        e = math.floor(rdl/8388608)%256
                        m  = rdl% 8388608
                        Value = (2.0*sign+1)*m*2.0**(e-150);   
                        #----------------
                        TempData.append(Value)
                        i+=1
            real = TempData[0:len(TempData):2]
            imag = TempData[1:len(TempData):2]
            data=[]
            for number in range(0,int(len(TempData)/2)):
                data.append(real[number]+1j*imag[number])
        data = np.array(data) 
        if 'FID' in TYPE:
            axis=0
            spec = [False]
        elif 'SPE' in TYPE:
            axis=1
            spec = [True]                    
        if NI is 1:
            masterData=sc.Spectrum(data,lambda self :self.LoadSimpsonFile(filePath),[0],[SW],spec)
        else:
            data = data.reshape((NI,NP))
            masterData=sc.Spectrum(data,lambda self :self.LoadSimpsonFile(filePath),[0,0],[SW1,SW],spec*2)
        return masterData
    
    def saveSimpsonFile(self):
        self.mainWindow.get_mainWindow().SaveSimpsonFile()
        
    def saveJSONFile(self):
        self.mainWindow.get_mainWindow().saveJSONFile()
        
    def saveMatlabFile(self):
        self.mainWindow.get_mainWindow().saveMatlabFile()
        
    def saveFigure(self):
        if self.mainWindow is not None:
            self.mainWindow.removeFromView()
        num = self.workspaces.index(self.mainWindow)
        self.mainWindow = fit.MainPlotWindow(self,self.mainWindow)
        self.mainFrame.addWidget(self.mainWindow)
        self.workspaces[num] = self.mainWindow
        self.mainWindow.addToView()
        self.menuCheck()

    def closeSaveFigure(self, mainWindow):
        self.mainWindow.removeFromView()
        num = self.workspaces.index(self.mainWindow)
        self.mainWindow = mainWindow
        self.workspaces[num] = self.mainWindow
        self.mainWindow.addToView()
        self.menuCheck()

    def createFitWindow(self,fitWindow):
        if self.mainWindow is not None:
            self.mainWindow.removeFromView()
        num = self.workspaces.index(self.mainWindow)
        self.mainWindow = fitWindow
        self.mainFrame.addWidget(self.mainWindow)
        self.workspaces[num] = self.mainWindow
        self.mainWindow.addToView()
        self.menuCheck()
        
    def closeFitWindow(self, mainWindow):
        self.mainWindow.removeFromView()
        num = self.workspaces.index(self.mainWindow)
        self.mainWindow = mainWindow
        self.workspaces[num] = self.mainWindow
        self.mainWindow.addToView()
        self.menuCheck()
        
    def fileQuit(self):
        self.close()

    def closeEvent(self, event):
        quit_msg = "Are you sure you want to exit the program?"
        reply = QtGui.QMessageBox.question(self, 'Close', 
                                           quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()
        
class Main1DWindow(QtGui.QWidget):
    def __init__(self,father,masterData,duplicateCurrent=None,name=''):
        QtGui.QWidget.__init__(self,father)
        self.name = name
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        grid = QtGui.QGridLayout(self)
        grid.addWidget(self.canvas,0,0)
        self.undoList = [] 
        self.redoList = []
        self.currentMacro = None
        self.redoMacro = []
        self.father = father
        self.mainProgram = self.father #remove all references to mainprogram to father                  
        self.masterData = masterData
        if duplicateCurrent is not None:
            self.current = duplicateCurrent.copyCurrent(self,self.fig,self.canvas,masterData)
        else:
            self.current = sc.Current1D(self,self.fig,self.canvas,masterData)
        self.menubar = self.mainProgram.menubar
        self.sideframe=SideFrame(self)
        grid.addWidget(self.sideframe,0,1)
        self.bottomframe=BottomFrame(self)
        grid.addWidget(self.bottomframe,1,0)
        self.textframe=TextFrame(self)
        grid.addWidget(self.textframe,2,0)
        grid.setColumnStretch(0,1)
        grid.setRowStretch(0,1)
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)

    def rename(self,name):
        self.name = name
        self.fig.suptitle(name)
        self.canvas.draw()
        
    def buttonPress(self,event):
        self.current.buttonPress(event)

    def buttonRelease(self,event):
        self.current.buttonRelease(event)

    def pan(self,event):
        self.current.pan(event)

    def scroll(self,event):
        self.current.scroll(event)
        
    def get_mainWindow(self):
        return self
        
    def get_masterData(self):
        return self.masterData

    def get_current(self):
        return self.current

    def kill(self):
        self.destroy()
        self.current.kill()
        del self.current
        del self.masterData
        del self.canvas
        del self.fig    #fig is destroyed

    def rescue(self):
        self.current.kill()
        self.current = sc.Current1D(self,self.current.fig,self.current.canvas,self.masterData)
        
    def initUI(self):
        self.multiDActions = []
        #the edit drop down menu
        self.editmenu = QtGui.QMenu("Edit",self)
        self.menubar.addMenu(self.editmenu)
        self.undoAction = self.editmenu.addAction("Undo",self.undo,QtGui.QKeySequence.Undo)
        self.undoAction.setShortcutContext(QtCore.Qt.WidgetShortcut)
        self.redoAction = self.editmenu.addAction("Redo",self.redo,QtGui.QKeySequence.Redo)
        self.redoAction.setShortcutContext(QtCore.Qt.WidgetShortcut)
        self.editmenu.addAction("Reload", self.reloadLast,QtGui.QKeySequence.Refresh)

	#the tool drop down menu
        self.toolMenu = QtGui.QMenu("Tools",self)
        self.menubar.addMenu(self.toolMenu)
        self.toolMenu.addAction("Real", self.real)
        self.toolMenu.addAction("Imag", self.imag)
        self.toolMenu.addAction("Abs", self.abs) 
        self.toolMenu.addAction("Apodize", self.createApodWindow)
        self.toolMenu.addAction("Phasing", self.createPhaseWindow)
        self.toolMenu.addAction("Swap Echo", self.createSwapEchoWindow)
        self.toolMenu.addAction("Offset correction", self.createDCWindow)
        self.toolMenu.addAction("Baseline correction", self.createBaselineWindow)
        self.toolMenu.addAction("States", self.states)
        self.toolMenu.addAction("States-TPPI", self.statesTPPI)
        self.toolMenu.addAction("Correct Bruker digital filter", self.BrukerDigital)

        #the matrix drop down menu
        self.matrixMenu = QtGui.QMenu("Matrix",self)
        self.menubar.addMenu(self.matrixMenu)
        self.matrixMenu.addAction("Sizing", self.createSizeWindow)
        self.matrixMenu.addAction("Shift Data", self.createShiftDataWindow)
        self.matrixMenu.addAction("Integrate", self.createIntegrateWindow)
        self.matrixMenu.addAction("Sum", self.createSumWindow)
        self.matrixMenu.addAction("Max", self.createMaxWindow)
        self.matrixMenu.addAction("Min", self.createMinWindow)
        self.matrixMenu.addAction("Max position", self.createArgMaxWindow)
        self.matrixMenu.addAction("Min position", self.createArgMinWindow)
        self.matrixMenu.addAction("Diff", self.diff)
        self.matrixMenu.addAction("Cumsum", self.cumsum)
        self.matrixMenu.addAction("Extract part", self.createRegionWindow)
        self.matrixMenu.addAction("Flip L/R", self.flipLR)
        self.matrixMenu.addAction("Delete", self.createDeleteWindow)
        self.matrixMenu.addAction("Split", self.createSplitWindow)
        self.matrixMenu.addAction("Multiply", self.createMultiplyWindow)
        self.multiDActions.append(self.matrixMenu.addAction("Concatenate", self.createConcatenateWindow))
        self.multiDActions.append(self.matrixMenu.addAction("Shearing", self.createShearingWindow))
        
        #the fft drop down menu
        self.fftMenu = QtGui.QMenu("Fourier",self)
        self.menubar.addMenu(self.fftMenu)
        self.fftMenu.addAction("Fourier transform", self.fourier, QtCore.Qt.CTRL + QtCore.Qt.Key_F)
        self.fftMenu.addAction("Real Fourier transform", self.realFourier)
        self.fftMenu.addAction("Fftshift", self.fftshift)
        self.fftMenu.addAction("Inv fftshift", self.invFftshift)
        self.fftMenu.addAction("Hilbert transform", self.hilbert)

	#the fitting drop down menu
        self.fittingMenu = QtGui.QMenu("Fitting",self)
        self.menubar.addMenu(self.fittingMenu)
        self.fittingMenu.addAction("S/N", self.createSNWindow)
        self.fittingMenu.addAction("FWHM", self.createFWHMWindow)
        self.fittingMenu.addAction("Relaxation Curve", self.createRelaxWindow)
        self.fittingMenu.addAction("Diffusion Curve", self.createDiffusionWindow)
        self.fittingMenu.addAction("Peak Deconvolution", self.createPeakDeconvWindow)
        self.fittingMenu.addAction("CSA tensor", self.createTensorDeconvWindow)
        self.fittingMenu.addAction("First order quadrupole", self.createQuad1DeconvWindow)
        self.fittingMenu.addAction("Second order quadrupole static", self.createQuad2StaticDeconvWindow)
        self.fittingMenu.addAction("Second order quadrupole MAS", self.createQuad2MASDeconvWindow)
        self.fittingMenu.addAction("Czjzek static", self.createQuad2StaticCzjzekWindow)
        self.fittingMenu.addAction("Czjzek MAS", self.createQuad2MASCzjzekWindow)
        
        #the combine drop down menu
        self.combineMenu = QtGui.QMenu("Combine",self)
        self.menubar.addMenu(self.combineMenu)
        self.combineMenu.addAction("Insert from workspace", self.createInsertWindow)
        self.combineMenu.addAction("Add", self.createAddWindow)
        self.combineMenu.addAction("Subtract", self.createSubtractWindow)

	#the plot drop down menu
        self.plotMenu = QtGui.QMenu("Plot",self)
        self.menubar.addMenu(self.plotMenu)
        self.plotMenu.addAction("1D plot", self.plot1D)
        self.plotMenu.addAction("Scatter plot", self.plotScatter)
        self.multiDActions.append(self.plotMenu.addAction("Stack plot", self.plotStack))
        self.multiDActions.append(self.plotMenu.addAction("Array plot", self.plotArray))
        self.multiDActions.append(self.plotMenu.addAction("Contour plot", self.plotContour))
        self.multiDActions.append(self.plotMenu.addAction("Skewed plot", self.plotSkewed))
        self.plotMenu.addAction("Set reference", self.createRefWindow)
        self.plotMenu.addAction("User x-axis", self.createXaxWindow)
        self.show()
        self.menuCheck()

    def addToView(self):
        self.editmenu.menuAction().setVisible(True)
        self.matrixMenu.menuAction().setVisible(True)
        self.plotMenu.menuAction().setVisible(True)
        self.toolMenu.menuAction().setVisible(True)
        self.fftMenu.menuAction().setVisible(True)
        self.fittingMenu.menuAction().setVisible(True)
        self.combineMenu.menuAction().setVisible(True)
        self.show()
        self.menuCheck()
        
    def removeFromView(self):
        self.editmenu.menuAction().setVisible(False)
        self.matrixMenu.menuAction().setVisible(False)
        self.plotMenu.menuAction().setVisible(False)
        self.toolMenu.menuAction().setVisible(False)
        self.fftMenu.menuAction().setVisible(False)
        self.fittingMenu.menuAction().setVisible(False)
        self.combineMenu.menuAction().setVisible(False)
        self.hide()

    def menuEnable(self):
        self.mainProgram.filemenu.menuAction().setEnabled(True)
        self.mainProgram.workspacemenu.menuAction().setEnabled(True)
        self.mainProgram.macromenu.menuAction().setEnabled(True)
        self.editmenu.menuAction().setEnabled(True)
        self.matrixMenu.menuAction().setEnabled(True)
        self.plotMenu.menuAction().setEnabled(True)
        self.toolMenu.menuAction().setEnabled(True)
        self.fftMenu.menuAction().setEnabled(True)
        self.fittingMenu.menuAction().setEnabled(True)
        self.combineMenu.menuAction().setEnabled(True)
        self.sideframe.frameEnable()
        self.bottomframe.frameEnable()
        self.textframe.frameEnable()
        self.menuCheck()

    def menuDisable(self):
        self.mainProgram.filemenu.menuAction().setEnabled(False)
        self.mainProgram.workspacemenu.menuAction().setEnabled(False)
        self.mainProgram.macromenu.menuAction().setEnabled(False)
        self.editmenu.menuAction().setEnabled(False)
        self.matrixMenu.menuAction().setEnabled(False)
        self.plotMenu.menuAction().setEnabled(False)
        self.toolMenu.menuAction().setEnabled(False)
        self.fftMenu.menuAction().setEnabled(False)
        self.fittingMenu.menuAction().setEnabled(False)
        self.combineMenu.menuAction().setEnabled(False)
        self.sideframe.frameDisable()
        self.bottomframe.frameDisable()
        self.textframe.frameDisable()

    def menuCheck(self):
        if (len(self.masterData.data.shape) < 2):
            for i in self.multiDActions:
                i.setEnabled(False)
        else:
            for i in self.multiDActions:
                i.setEnabled(True)
        if not self.undoList:
            self.undoAction.setEnabled(False)
        else:
            self.undoAction.setEnabled(True)
        if not self.redoList:
            self.redoAction.setEnabled(False)
        else:
            self.redoAction.setEnabled(True)

    def runMacro(self,macro):
        self.redoList = []
        for iter1 in macro:
            if iter1[0] == 'reload':
                self.undoList.append(self.masterData.reload(self))
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
            elif iter1[0] == 'baselineCorrection':
                self.undoList.append(self.masterData.baselineCorrection(*iter1[1]))
            elif iter1[0] == 'integrate':
                self.undoList.append(self.masterData.matrixManip(*iter1[1],which=0))
            elif iter1[0] == 'sum':
                self.undoList.append(self.masterData.matrixManip(*iter1[1],which=5))
            elif iter1[0] == 'max':
                self.undoList.append(self.masterData.matrixManip(*iter1[1],which=1))
            elif iter1[0] == 'min':
                self.undoList.append(self.masterData.matrixManip(*iter1[1],which=2))
            elif iter1[0] == 'argmax':
                self.undoList.append(self.masterData.matrixManip(*iter1[1],which=3))
            elif iter1[0] == 'argmin':
                self.undoList.append(self.masterData.matrixManip(*iter1[1],which=4))
            elif iter1[0] == 'fliplr':
                self.undoList.append(self.masterData.flipLR(*iter1[1]))
            elif iter1[0] == 'concatenate':
                self.undoList.append(self.masterData.concatenate(*iter1[1]))
            elif iter1[0] == 'split':
                self.undoList.append(self.masterData.split(*iter1[1]))
            elif iter1[0] == 'insert':
                self.undoList.append(self.masterData.insert(*iter1[1]))
            elif iter1[0] == 'delete':
                self.undoList.append(self.masterData.delete(*iter1[1]))
            elif iter1[0] == 'add':
                self.undoList.append(self.masterData.add(*iter1[1]))
            elif iter1[0] == 'subtract':
                self.undoList.append(self.masterData.subtract(*iter1[1]))
            elif iter1[0] == 'multiply':
                self.undoList.append(self.masterData.multiply(*iter1[1]))
            elif iter1[0] == 'shear':
                self.undoList.append(self.masterData.shear(*iter1[1]))
            elif iter1[0] == 'setxax':
                self.undoList.append(self.masterData.setXax(*iter1[1]))
            elif iter1[0] == 'hilbert':
                self.undoList.append(self.masterData.hilbert(*iter1[1]))
            else:
                print('unknown macro command: '+iter1[0])
        self.current.upd()   #get the first slice of data
        self.current.plotReset() #reset the axes limits
        self.current.showFid() #plot the data
        self.updAllFrames()
        self.menuCheck()

    def addMacro(self,macroStep):
        if self.currentMacro is not None:
            self.mainProgram.macroAdd(self.currentMacro,macroStep)
            self.redoMacro = []

    def saveJSONFile(self):
        name = QtGui.QFileDialog.getSaveFileName(self, 'Save File','','JSON (*.json)')
        if not name:
            return
        struct = {}
        struct['dataReal'] = np.real(self.masterData.data).tolist()
        struct['dataImag'] = np.imag(self.masterData.data).tolist()
        struct['freq'] = self.masterData.freq.tolist()
        struct['sw'] = list(self.masterData.sw)
        struct['spec'] = list(self.masterData.spec)
        struct['wholeEcho'] = list(np.array(self.masterData.wholeEcho,dtype=int))
        struct['ref'] = np.array(self.masterData.ref,dtype=np.float).tolist()
        tmpXax = []
        for i in self.masterData.xaxArray:
            tmpXax.append(i.tolist())
        struct['xaxArray'] = tmpXax
        with open(name, 'w') as outfile:
            json.dump(struct, outfile)

    def saveMatlabFile(self):
        name = QtGui.QFileDialog.getSaveFileName(self, 'Save File','','MATLAB file (*.mat)')
        if not name:
            return
        struct = {}
        struct['dim'] = self.masterData.dim
        struct['data'] = self.masterData.data
        struct['freq'] = self.masterData.freq
        struct['sw'] = self.masterData.sw
        struct['spec'] = self.masterData.spec
        struct['wholeEcho'] = self.masterData.wholeEcho
        struct['ref'] = np.array(self.masterData.ref,dtype=np.float)
        struct['xaxArray'] = self.masterData.xaxArray
        matlabStruct = {self.mainProgram.workspaceNames[self.mainProgram.workspaceNum]:struct}
        scipy.io.savemat(name,matlabStruct)          

    def SaveSimpsonFile(self):
        if self.masterData.dim   > 2:
            print('Saving to Simpson format only allowed for 1D and 2D data!')
            return
        if sum(self.masterData.spec)/len(self.masterData.spec)==1:
            name = QtGui.QFileDialog.getSaveFileName(self, 'Save File','','SIMPSON file (*.spe)')
        elif sum(self.masterData.spec) == 0: 
            name = QtGui.QFileDialog.getSaveFileName(self, 'Save File','','SIMPSON file (*.fid)')
        with open(name) as f: 
            f.write('SIMP\n')
            if self.masterData.dim  is 2:
                f.write('NP='+str(self.masterData.data.shape[1])+'\n')
                f.write('NI='+str(self.masterData.data.shape[0])+'\n')
                f.write('SW='+str(self.masterData.sw[1])+'\n')
                f.write('SW1='+str(self.masterData.sw[0])+'\n')
            else:
                f.write('NP='+str(self.masterData.data.shape[0])+'\n')
                f.write('SW='+str(self.masterData.sw[0])+'\n')
            if self.masterData.spec[0]:
                f.write('TYPE=SPE'+'\n') 
            else:
                f.write('TYPE=FID'+'\n') 
            f.write('DATA'+'\n')
            if self.masterData.dim  is 1:
                for Line in self.masterData.data:
                    f.write(str(Line.real)+' '+ str(Line.imag)+'\n')
            if self.masterData.dim  is 2:
                Points= self.masterData.data.shape
                for iii in range(0,Points[0]):
                    for jjj in range(0,Points[1]):
                        f.write(str(self.masterData.data[iii][jjj].real)+' '+ str(self.masterData.data[iii][jjj].imag)+'\n')
            f.write('END')

    def reloadLast(self):
        self.redoList = []
        self.undoList.append(self.masterData.reload(self.mainProgram))
        self.current.upd()
        self.current.plotReset()
        self.current.showFid()
        self.updAllFrames()
        self.addMacro(['reload'])
        self.menuCheck()

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
        
    def setFreq(self,freq,sw):
        self.redoList = []
        self.undoList.append(self.current.setFreq(freq,sw))
        self.menuCheck()

    def flipLR(self):
        self.redoList = []
        self.undoList.append(self.current.flipLR())
        self.menuCheck()
        
    def createPhaseWindow(self):
        self.extraWindow = PhaseWindow(self)
        
    def createApodWindow(self):
        self.extraWindow = ApodWindow(self)

    def createSizeWindow(self):
        self.extraWindow = SizeWindow(self)

    def createSwapEchoWindow(self):
        self.extraWindow = SwapEchoWindow(self)

    def createShiftDataWindow(self):
        self.extraWindow = ShiftDataWindow(self)

    def createDCWindow(self):
        self.extraWindow = DCWindow(self)
        
    def createBaselineWindow(self):
        self.extraWindow = BaselineWindow(self)
        
    def createRefWindow(self):
        self.extraWindow = RefWindow(self)

    def createIntegrateWindow(self):
        self.extraWindow = integrateWindow(self)
        
    def createSumWindow(self):
        self.extraWindow = sumWindow(self)
        
    def createMultiplyWindow(self):
        self.extraWindow = MultiplyWindow(self)
        
    def createMaxWindow(self):
        self.extraWindow = maxWindow(self)
        
    def createMinWindow(self):
        self.extraWindow = minWindow(self)

    def createArgMaxWindow(self):
        self.extraWindow = argmaxWindow(self)
        
    def createArgMinWindow(self):
        self.extraWindow = argminWindow(self)
        
    def createRegionWindow(self):
        self.extraWindow = extractRegionWindow(self)

    def createDeleteWindow(self):
        self.extraWindow = DeleteWindow(self)
        
    def createSplitWindow(self):
        self.extraWindow = SplitWindow(self)
        
    def createConcatenateWindow(self):
        self.extraWindow = ConcatenateWindow(self)
        
    def createInsertWindow(self):
        self.extraWindow = InsertWindow(self)
        
    def createAddWindow(self):
        self.extraWindow = AddWindow(self)
        
    def createSubtractWindow(self):
        self.extraWindow = SubtractWindow(self)
        
    def createShearingWindow(self):
        if self.masterData.dim > 1:
            self.extraWindow = ShearingWindow(self)
        else:
            print('Data has too little dimensions for shearing transform')

    def BrukerDigital(self):
        FilePath = QtGui.QFileDialog.getOpenFileName(self, 'Open File')
        if FilePath is '':
            return
        Dir = os.path.dirname(FilePath)
        if not os.path.exists(Dir+os.path.sep+'acqus'):
            print("acqus file does not exist")
            return
        with open(Dir+os.path.sep+'acqus', 'r') as f: 
            data = f.read().split('\n')
        FilterCorrection = -1.0
        for s in range(0,len(data)):
            if data[s].startswith('##$GRPDLY='):
                FilterCorrection = float(data[s][10:])
            if data[s].startswith('##$DECIM='):
                DECIM = int(float(data[s][9:]))
            if data[s].startswith('##$DSPFVS='):
                DSPFVS = int(float(data[s][10:]))
        if DSPFVS == 10 or DSPFVS == 11 or DSPFVS == 12:#get from table
            CorrectionList = [{'2':44.7500,'3':33.5000,'4':66.6250,'6':59.0833
                               ,'8':68.5625,'12':60.3750,'16':69.5313,'24':61.0208,'32':70.0156
                               ,'48':61.3438,'64':70.2578,'96':61.5052,'128':70.3789,'192':61.5859
                               ,'256':70.4395,'384':61.6263,'512':70.4697,'768':61.6465,'1024':70.4849,'1536':61.6566,'2048':70.4924},
                              {'2':46.0000,'3':36.5000,'4':48.0000,'6':50.1667,'8':53.2500,'12':69.5000,
                               '16':72.2500,'24':70.1667,'32':72.7500,'48':70.5000,'64':73.0000,'96':70.6667,
                               '128':72.5000,'192':71.3333,'256':72.2500,'384':71.6667,'512':72.1250,'768':71.8333,
                               '1024':72.0625,'1536':71.9167,'2048':72.0313},
                              {'2':46.311,'3':36.530,'4':47.870,'6':50.229,'8':53.289,'12':69.551,'16':71.600,
                               '24':70.184,'32':72.138,'48':70.528,'64':72.348,'96':70.700,'128':72.524}]
            #Take correction from database. Based on matNMR routine (Jacco van Beek), which is itself based 
            #on a text by W. M. Westler and F. Abildgaard.
            FilterCorrection = CorrectionList[10-DSPFVS][str(DECIM)]
        if FilterCorrection == -1.0:
            print('DSPFVS value not recognized (Bruker hardware version not known)')
            return
        if FilterCorrection != -1.0: #If changed
            self.redoList = []
            self.undoList.append(self.current.applyPhase(0, FilterCorrection*2*np.pi))
            self.menuCheck()
                    
    def createSNWindow(self):
        self.extraWindow = SNWindow(self)
        
    def createFWHMWindow(self):
        self.extraWindow = FWHMWindow(self)

    def createXaxWindow(self):
        self.extraWindow = XaxWindow(self)
        
    def createRelaxWindow(self):
        self.mainProgram.createFitWindow(fit.RelaxWindow(self.father,self.mainProgram,self.mainProgram.mainWindow))
        
    def createDiffusionWindow(self):
        self.mainProgram.createFitWindow(fit.DiffusionWindow(self.father,self.mainProgram,self.mainProgram.mainWindow))
        
    def createPeakDeconvWindow(self):
        self.mainProgram.createFitWindow(fit.PeakDeconvWindow(self.father,self.mainProgram,self.mainProgram.mainWindow))
        
    def createTensorDeconvWindow(self):
        self.mainProgram.createFitWindow(fit.TensorDeconvWindow(self.father,self.mainProgram,self.mainProgram.mainWindow))

    def createQuad1DeconvWindow(self):
        self.mainProgram.createFitWindow(fit.Quad1DeconvWindow(self.father,self.mainProgram,self.mainProgram.mainWindow))
        
    def createQuad2StaticDeconvWindow(self):
        self.mainProgram.createFitWindow(fit.Quad2DeconvWindow(self.father,self.mainProgram,self.mainProgram.mainWindow))
        
    def createQuad2MASDeconvWindow(self):
        self.mainProgram.createFitWindow(fit.Quad2DeconvWindow(self.father,self.mainProgram,self.mainProgram.mainWindow,True))

    def createQuad2StaticCzjzekWindow(self):
        self.mainProgram.createFitWindow(fit.Quad2CzjzekWindow(self.father,self.mainProgram,self.mainProgram.mainWindow))
        
    def createQuad2MASCzjzekWindow(self):
        self.mainProgram.createFitWindow(fit.Quad2CzjzekWindow(self.father,self.mainProgram,self.mainProgram.mainWindow,True))
        
    def plot1D(self):
        tmpcurrent = sc.Current1D(self,self.fig,self.canvas,self.masterData,self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()

    def plotScatter(self):
        tmpcurrent = sc.CurrentScatter(self,self.fig,self.canvas,self.masterData,self.current)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()

    def plotStack(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentStacked(self,self.fig,self.canvas,self.masterData,self.current)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")

    def plotArray(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentArrayed(self,self.fig,self.canvas,self.masterData,self.current)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")
            
    def plotContour(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentContour(self,self.fig,self.canvas,self.masterData,self.current)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")
            
    def plotSkewed(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentSkewed(self,self.fig,self.canvas,self.masterData,self.current)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")

    def updAllFrames(self):
        self.sideframe.upd()
        self.bottomframe.upd()

    def undo(self, *args):
        undoFunc = None
        while undoFunc is None and self.undoList:
            undoFunc = self.undoList.pop()
        if undoFunc is None:
            print("no undo information")
            return
        self.redoList.append(undoFunc(self.masterData))
        self.current.upd()
        self.current.plotReset()
        self.current.showFid()
        self.updAllFrames()
        if self.currentMacro is not None:
            self.redoMacro.append(self.mainProgram.macros[self.currentMacro].pop())
        self.menuCheck()

    def redo(self, *args):
        if self.redoList:
            self.undoList.append(self.redoList.pop()(self.masterData))
            self.current.upd()
            self.current.plotReset()
            self.current.showFid()
            self.updAllFrames()
            if self.currentMacro is not None:
                self.mainProgram.macroAdd(self.currentMacro,self.redoMacro.pop())
            self.menuCheck()
        else:
            print("no redo information")

########################################################################################
class SideFrame(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.father = parent
        self.entries=[]
        self.plotIs2D = False
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        splitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        frame1Widget = QtGui.QWidget()
        frame2Widget = QtGui.QWidget()
        splitter.addWidget(frame1Widget)
        splitter.addWidget(frame2Widget)
        splitter.setStretchFactor(1,1)
        grid.addWidget(splitter)
        self.frame1 = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        frame1Widget.setLayout(self.frame1)
        frame2Widget.setLayout(self.frame2)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.upd()

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
        self.plotIs2D = isinstance(current, (sc.CurrentStacked,sc.CurrentArrayed,sc.CurrentContour,sc.CurrentSkewed))
        if self.plotIs2D:
            offset = 1
        self.entries=[]
        self.buttons1=[]
        self.buttons1Group = QtGui.QButtonGroup(self)
        self.buttons1Group.buttonClicked.connect(lambda: self.setAxes(True))
        self.buttons2=[]
        self.buttons2Group = QtGui.QButtonGroup(self)
        self.buttons2Group.buttonClicked.connect(lambda: self.setAxes(False))
        if self.length > 1:
            for num in range(self.length):
                self.buttons1.append(QtGui.QRadioButton(''))
                self.buttons1Group.addButton(self.buttons1[num],num)
                self.frame1.addWidget(self.buttons1[num],num*2+1,0)
                if self.plotIs2D:
                    self.buttons2.append(QtGui.QRadioButton(''))
                    self.buttons2Group.addButton(self.buttons2[num],num)
                    self.frame1.addWidget(self.buttons2[num],num*2+1,1)
                self.frame1.addWidget(QLabel("D"+str(num+1),self),num*2,1+offset)
                self.entries.append(SliceSpinBox(self,0,self.shape[num]-1))
                self.frame1.addWidget(self.entries[num],num*2+1,1+offset)
                if not self.plotIs2D:
                    if num < current.axes:
                        self.entries[num].setValue(current.locList[num])
                    elif num == current.axes:
                        self.entries[num].setValue(0)
                    else:
                        self.entries[num].setValue(current.locList[num-1])
                else:
                    if (num < current.axes) and (num < current.axes2):
                        self.entries[num].setValue(current.locList[num])
                    elif (num == current.axes) or (num == current.axes2):
                        self.entries[num].setValue(0)
                    elif (num > current.axes) and (num > current.axes2):
                        self.entries[num].setValue(current.locList[num-2])
                    else:
                        self.entries[num].setValue(current.locList[num-1])
                self.entries[num].valueChanged.connect(lambda event=None,num=num: self.getSlice(event,num))
            if isinstance(current, (sc.CurrentStacked,sc.CurrentArrayed,sc.CurrentSkewed)):
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
                self.frame2.addWidget(QLabel("From",self),1,0)
                self.fromSpin = SliceSpinBox(self,0,to2D-1)
                self.frame2.addWidget(self.fromSpin,2,0)
                self.fromSpin.setValue(from2D)
                self.fromSpin.valueChanged.connect(self.setToFrom)
                self.frame2.addWidget(QLabel("To",self),3,0)
                self.toSpin = SliceSpinBox(self,from2D+1,self.shape[current.axes2])
                self.frame2.addWidget(self.toSpin,4,0)
                self.toSpin.setValue(to2D)
                self.toSpin.valueChanged.connect(self.setToFrom)
                self.frame2.addWidget(QLabel("Step",self),5,0)
                self.stepSpin = SliceSpinBox(self,1,self.shape[current.axes2])
                self.frame2.addWidget(self.stepSpin,6,0)
                self.stepSpin.setValue(step2D)
                self.stepSpin.valueChanged.connect(self.setToFrom)
                if isinstance(current, (sc.CurrentStacked,sc.CurrentArrayed)):
                    self.frame2.addWidget(QLabel("Spacing",self),7,0)
                    self.spacingEntry = QtGui.QLineEdit(self)
                    self.spacingEntry.setText('%#.3g' % current.spacing)
                    self.spacingEntry.editingFinished.connect(self.setSpacing)
                    self.frame2.addWidget(self.spacingEntry,8,0)
                elif isinstance(current, (sc.CurrentSkewed)):
                    self.frame2.addWidget(QLabel("Skew",self),7,0)
                    self.skewEntry = QtGui.QLineEdit(self)
                    self.skewEntry.setText('%.2f' % current.skewed)
                    self.skewEntry.editingFinished.connect(self.setSkew)
                    self.frame2.addWidget(self.skewEntry,8,0)
                    self.frame2.addWidget(QLabel("Elevation",self),9,0)
                    self.elevEntry = QtGui.QLineEdit(self)
                    self.elevEntry.setText('%.1f' % current.elevation)
                    self.elevEntry.editingFinished.connect(self.setSkew)
                    self.frame2.addWidget(self.elevEntry,10,0)
            if isinstance(current, (sc.CurrentContour)):
                self.frame2.addWidget(QLabel("Number of contours",self),1,0)
                self.numLEntry = QtGui.QLineEdit(self)
                self.numLEntry.setText(str(current.numLevels))
                self.numLEntry.editingFinished.connect(self.setContour)
                self.frame2.addWidget(self.numLEntry,2,0)
                self.frame2.addWidget(QLabel("Highest contour [%]",self),3,0)
                self.maxLEntry = QtGui.QLineEdit(self)
                self.maxLEntry.setText(str(current.maxLevels*100.0))
                self.maxLEntry.editingFinished.connect(self.setContour)
                self.frame2.addWidget(self.maxLEntry,4,0)
                self.frame2.addWidget(QLabel("Lowest contour [%]",self),5,0)
                self.minLEntry = QtGui.QLineEdit(self)
                self.minLEntry.setText(str(current.minLevels*100.0))
                self.minLEntry.editingFinished.connect(self.setContour)
                self.frame2.addWidget(self.minLEntry,6,0)
            self.buttons1Group.button(current.axes).toggle()
            if self.plotIs2D:
                self.buttons2Group.button(current.axes2).toggle()
                    
    def setToFrom(self, *args):
        current = self.father.current
        if not isinstance(current, (sc.CurrentStacked,sc.CurrentArrayed,sc.CurrentSkewed)):
            return
        fromVar = self.fromSpin.value()
        toVar = self.toSpin.value()
        stepVar = self.stepSpin.value()
        current.stackSelect(fromVar,toVar,stepVar)
        self.fromSpin.setMaximum(toVar-1)
        self.toSpin.setMinimum(fromVar+1)

    def scrollSpacing(self, var):
        self.spacingEntry.setText('%#.3g' % var)
            
    def setSpacing(self, *args):
        var = float(safeEval(self.spacingEntry.text()))
        self.spacingEntry.setText('%#.3g' % var)
        self.father.current.setSpacing(var)

    def setSkew(self, *args):
        var =float(safeEval(self.skewEntry.text()))
        self.skewEntry.setText('%.2f' % var)
        var2 =float(safeEval(self.elevEntry.text()))
        self.elevEntry.setText('%.1f' % var2)
        self.father.current.setSkewed(var,var2)
        
    def setContour(self, *args):
        var1 = int(round(safeEval(self.numLEntry.text())))
        self.numLEntry.setText(str(var1))
        var2 =float(safeEval(self.maxLEntry.text()))
        self.maxLEntry.setText('%.1f' % var2)
        var3 =float(safeEval(self.minLEntry.text()))
        self.minLEntry.setText('%.1f' % var3)
        self.father.current.setLevels(var1,var2/100.0,var3/100.0)
        
    def setAxes(self,first=True):
        axes = self.buttons1Group.checkedId()
        if self.plotIs2D:
            axes2 = self.buttons2Group.checkedId()
            if axes==axes2:
                if first:
                    axes2 = self.father.current.axes
                else:
                    axes = self.father.current.axes2
            self.buttons2Group.button(axes2).toggle()
        self.getSlice(None, axes,True)

    def getSlice(self, event, entryNum, button=False):
        if button:
            dimNum = entryNum
        elif not self.plotIs2D:
            if entryNum == self.father.current.axes:
                if entryNum == self.length-1:
                    dimNum = self.length-2
                else:
                    dimNum = self.length-1
            else:
                dimNum = self.father.current.axes
        else:
            dimNum = self.father.current.axes

        locList=[]
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
            self.father.current.setBlock(dimNum,self.buttons2Group.checkedId(), locList)
        else:
            self.father.current.setSlice(dimNum,locList)
        self.father.bottomframe.upd()
            
################################################################################  
class BottomFrame(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.father = parent
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        fourierButton = QtGui.QPushButton("Fourier")
        fourierButton.clicked.connect(self.father.fourier)
        grid.addWidget(fourierButton,0,0,2,1)
        self.specGroup = QtGui.QButtonGroup(self)
        self.specGroup.buttonClicked.connect(self.changeSpec)
        timeButton = QtGui.QRadioButton('Time')
        self.specGroup.addButton(timeButton,0)
        grid.addWidget(timeButton,0,1)
        freqButton = QtGui.QRadioButton('Frequency')
        self.specGroup.addButton(freqButton,1)
        grid.addWidget(freqButton,1,1)
        self.wholeEcho = QtGui.QCheckBox("Whole echo")
        self.wholeEcho.clicked.connect(self.setWholeEcho)
        grid.addWidget(self.wholeEcho,0,2,2,1)
        grid.addWidget(QLabel("Freq [MHz]"),0,3)
        self.freqEntry = QtGui.QLineEdit(self)
        self.freqEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.freqEntry.editingFinished.connect(self.changeFreq)
        grid.addWidget(self.freqEntry,1,3)
        grid.addWidget(QLabel("Sweepwidth [kHz]"),0,4)
        self.swEntry = QtGui.QLineEdit()
        self.swEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.swEntry.editingFinished.connect(self.changeFreq)
        grid.addWidget(self.swEntry,1,4)
        grid.addWidget(QLabel("Plot"),0,5)
        self.plotDrop = QtGui.QComboBox()
        self.plotDrop.addItems(["Real", "Imag", "Both","Abs"])
        self.plotDrop.activated.connect(self.changePlot)
        grid.addWidget(self.plotDrop,1,5)
        grid.addWidget(QLabel("Axis"),0,6)
        self.axisDropTime = QtGui.QComboBox()
        self.axisDropTime.addItems(["s", "ms", u"\u03bcs"])
        self.axisDropTime.activated.connect(self.changeAxis)
        grid.addWidget(self.axisDropTime,1,6)
        self.axisDropFreq = QtGui.QComboBox()
        self.axisDropFreq.addItems(["Hz", "kHz", "MHz","ppm"])
        self.axisDropFreq.activated.connect(self.changeAxis)
        grid.addWidget(self.axisDropFreq,1,6)
        self.ax2Label = QLabel("Axis2")
        grid.addWidget(self.ax2Label,0,7)
        self.axisDropTime2 = QtGui.QComboBox()
        self.axisDropTime2.addItems(["s", "ms", u"\u03bcs"])
        self.axisDropTime2.activated.connect(self.changeAxis2)
        grid.addWidget(self.axisDropTime2,1,7)
        self.axisDropFreq2 = QtGui.QComboBox()
        self.axisDropFreq2.addItems(["Hz", "kHz", "MHz","ppm"])
        self.axisDropFreq2.activated.connect(self.changeAxis2)
        grid.addWidget(self.axisDropFreq2,1,7)
        grid.setColumnStretch(10,1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        self.upd()
        
    def frameEnable(self):
        self.setEnabled(True)
            
    def frameDisable(self):
        self.setEnabled(False)
        
    def upd(self): 
        self.freqEntry.setText('%.6f' %(self.father.current.freq/1000000.0))
        self.swEntry.setText('%.6f' %(self.father.current.sw/1000.0)) 
        self.axisDropTime2.hide()
        self.axisDropFreq2.hide()
        if self.father.current.spec==0:
            self.specGroup.button(0).toggle()
            self.axisDropFreq.hide()
            self.axisDropTime.show()
            self.ax2Label.hide()
            self.axisDropTime.setCurrentIndex(self.father.current.axType)
        elif self.father.current.spec==1:
            self.specGroup.button(1).toggle()
            self.axisDropTime.hide()
            self.axisDropFreq.show()
            self.ax2Label.hide()
            self.axisDropFreq.setCurrentIndex(self.father.current.axType)
        if isinstance(self.father.current,sc.CurrentContour):
            self.ax2Label.show()
            if self.father.current.spec2==0:
                self.axisDropTime2.show()
                self.axisDropTime2.setCurrentIndex(self.father.current.axType2)
            elif self.father.current.spec2==1:
                self.axisDropFreq2.show()
                self.axisDropFreq2.setCurrentIndex(self.father.current.axType2)
        if self.father.current.wholeEcho:
            self.wholeEcho.setCheckState(QtCore.Qt.Checked)
        else:
            self.wholeEcho.setCheckState(QtCore.Qt.Unchecked)

    def setWholeEcho(self, inp):
        self.father.undoList.append(self.father.current.setWholeEcho(inp/2))
        self.father.menuCheck()

    def changeSpec(self):
        self.father.redoList = []
        self.father.undoList.append(self.father.current.changeSpec(self.specGroup.checkedId()))
        self.upd()
        self.father.menuCheck()

    def changeFreq(self):
        freq = safeEval(self.freqEntry.text())*1e6
        sw = safeEval(self.swEntry.text())*1e3
        if freq != 0 and sw != 0:
            self.father.setFreq(freq,sw)
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
class TextFrame(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.father = parent
        self.oldx = 0.0
        self.oldy = 0.0
        grid = QtGui.QGridLayout(self)
        getButton = QtGui.QPushButton("Get Position")
        getButton.clicked.connect(self.getPosition)
        grid.addWidget(getButton,0,1)
        grid.addWidget(QLabel("Position:"),0,2)
        self.pos = QtGui.QLineEdit()
        self.pos.setAlignment(QtCore.Qt.AlignHCenter)
        self.pos.setText("0")
        grid.addWidget(self.pos,0,3)
        grid.addWidget(QLabel("x-value:"),0,4)
        self.xpoint = QtGui.QLineEdit()
        self.xpoint.setAlignment(QtCore.Qt.AlignHCenter)
        self.xpoint.setText("0.0")
        grid.addWidget(self.xpoint,0,5)
        grid.addWidget(QLabel("y-value:"),0,6)
        self.ypoint = QtGui.QLineEdit()
        self.ypoint.setAlignment(QtCore.Qt.AlignHCenter)
        self.ypoint.setText("0.0")
        grid.addWidget(self.ypoint,0,7)
        grid.addWidget(QLabel(u"\u0394x:"),0,8)
        self.deltaxpoint = QtGui.QLineEdit()
        self.deltaxpoint.setAlignment(QtCore.Qt.AlignHCenter)
        self.deltaxpoint.setText("0.0")
        grid.addWidget(self.deltaxpoint,0,9)
        grid.addWidget(QLabel(u"\u0394y:"),0,10)
        self.deltaypoint = QtGui.QLineEdit()
        self.deltaypoint.setAlignment(QtCore.Qt.AlignHCenter)
        self.deltaypoint.setText("0.0")
        grid.addWidget(self.deltaypoint,0,11)
        grid.setColumnStretch(20,1)
        
    def frameEnable(self):
        for child in self.children():
            child.setEnabled(True)
            
    def frameDisable(self):
        for child in self.children():
            child.setEnabled(False)
        
    def setLabels(self,position):
        self.deltaxpoint.setText('%#.3g' % np.abs(self.oldx-position[1]))
        self.deltaypoint.setText('%#.3g' % np.abs(self.oldy-position[2]))
        self.pos.setText(str(position[0]))
        self.xpoint.setText('%#.3g' % position[1])
        self.ypoint.setText('%#.3g' % position[2])
        self.oldx = position[1]
        self.oldy = position[2]

    def getPosition(self, *args):
        self.father.current.peakPickFunc = lambda pos,self=self: self.setLabels(pos) 
        self.father.current.peakPick = True

#################################################################################   
class PhaseWindow(QtGui.QWidget): 

    RESOLUTION = 10000.0
    P1LIMIT = 540.0
    PHASE0STEP = 1.0
    PHASE1STEP = 1.0
    
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent.father)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.zeroVal = 0.0
        self.firstVal = 0.0
        self.refVal = 0.0
        self.available = True
        self.setWindowTitle("Phasing")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Zero order phasing:"),0,0,1,3)
        autoZero = QtGui.QPushButton("Autophase 0th")
        autoZero.clicked.connect(lambda: self.autophase(0))
        grid.addWidget(autoZero,1,1)
        self.zeroEntry = QtGui.QLineEdit()
        self.zeroEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.zeroEntry.returnPressed.connect(self.inputZeroOrder)
        self.zeroEntry.setText("0.000")
        grid.addWidget(self.zeroEntry,2,1)
        leftZero = QtGui.QPushButton("<")
        leftZero.clicked.connect(lambda:self.stepPhase(-1,0))
        leftZero.setAutoRepeat(True)
        grid.addWidget(leftZero,2,0)
        rightZero = QtGui.QPushButton(">")
        rightZero.clicked.connect(lambda:self.stepPhase(1,0))
        rightZero.setAutoRepeat(True)
        grid.addWidget(rightZero,2,2)
        self.zeroScale = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.zeroScale.setRange(-self.RESOLUTION,self.RESOLUTION)
        self.zeroScale.valueChanged.connect(self.setZeroOrder)
        grid.addWidget(self.zeroScale,3,0,1,3)
        grid.addWidget(QLabel("First order phasing:"),4,0,1,3)
        autoFirst = QtGui.QPushButton("Autophase 0th+1st")
        autoFirst.clicked.connect(lambda: self.autophase(1))
        grid.addWidget(autoFirst,5,1)
        self.firstEntry = QtGui.QLineEdit()
        self.firstEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.firstEntry.returnPressed.connect(self.inputFirstOrder)
        self.firstEntry.setText("0.000")
        grid.addWidget(self.firstEntry,6,1)
        leftFirst = QtGui.QPushButton("<")
        leftFirst.clicked.connect(lambda:self.stepPhase(0,-1))
        leftFirst.setAutoRepeat(True)
        grid.addWidget(leftFirst,6,0)
        rightFirst = QtGui.QPushButton(">")
        rightFirst.clicked.connect(lambda:self.stepPhase(0,1))
        rightFirst.setAutoRepeat(True)
        grid.addWidget(rightFirst,6,2)
        self.firstScale = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.firstScale.setRange(-self.RESOLUTION,self.RESOLUTION)
        self.firstScale.valueChanged.connect(self.setFirstOrder)
        grid.addWidget(self.firstScale,7,0,1,3)
        if self.father.current.spec > 0:
            grid.addWidget(QLabel("Reference:"),8,0,1,3)
            pickRef = QtGui.QPushButton("Pick reference")
            pickRef.clicked.connect(self.pickRef)
            grid.addWidget(pickRef,9,1)
            self.refEntry = QtGui.QLineEdit()
            self.refEntry.setAlignment(QtCore.Qt.AlignHCenter)
            self.refEntry.setText('%.3f' % self.refVal)
            self.refEntry.editingFinished.connect(self.inputRef)
            grid.addWidget(self.refEntry,10,1)

        self.singleSlice = QtGui.QCheckBox("Single slice")
        grid.addWidget(self.singleSlice,11,0,1,3)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,1,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,1,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)
        
    def setZeroOrder(self,value, *args):
        if self.available:
            self.zeroVal = float(value)/self.RESOLUTION*180
            self.zeroEntry.setText('%.3f' % self.zeroVal)
            self.father.current.setPhaseInter(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0)
        
    def inputZeroOrder(self, *args):
        inp = safeEval(self.zeroEntry.text())
        self.zeroVal = np.mod(inp+180,360)-180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.available = False
        self.zeroScale.setValue(round(self.zeroVal/180.0*self.RESOLUTION))
        self.available = True
        self.father.current.setPhaseInter(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0)

    def setFirstOrder(self,value, *args):
        if self.available:
            value = float(value)/self.RESOLUTION*self.P1LIMIT
            newZero = (self.zeroVal-(value-self.firstVal)*self.refVal/self.father.current.sw) 
            self.zeroVal = np.mod(newZero+180,360)-180
            self.zeroEntry.setText('%.3f' % self.zeroVal)
            self.firstVal = value
            self.firstEntry.setText('%.3f' % self.firstVal)
            self.available = False
            self.zeroScale.setValue(round(self.zeroVal/180.0*self.RESOLUTION))
            self.available = True
            self.father.current.setPhaseInter(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0)

    def inputFirstOrder(self, *args): 
        value = float(safeEval(self.firstEntry.text()))
        newZero = (self.zeroVal-(value-self.firstVal)*self.refVal/self.father.current.sw)
        self.zeroVal = np.mod(newZero+180,360)-180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.firstVal = value
        self.firstEntry.setText('%.3f' % self.firstVal)
        self.available = False
        self.zeroScale.setValue(round(self.zeroVal/180.0*self.RESOLUTION))
        self.firstScale.setValue(round(self.firstVal/self.P1LIMIT*self.RESOLUTION))
        self.available = True
        self.father.current.setPhaseInter(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0)

    def autophase(self, num):
        phases = self.father.current.autoPhase(num)
        val = phases[0]/np.pi*180.0
        self.zeroVal=(np.mod(val+180,360)-180)
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.available = False
        self.zeroScale.setValue(round(self.zeroVal/180.0*self.RESOLUTION))
        self.available = True
        if num == 1:
            val = phases[1]/np.pi*180.0
            self.firstVal = val
            self.firstEntry.setText('%.3f' % self.firstVal)
        self.inputFirstOrder()

    def stepPhase(self,phase0,phase1): 
        inp = safeEval(self.zeroEntry.text())+phase0*self.PHASE0STEP
        self.zeroVal = np.mod(inp+180,360)-180
        value = safeEval(self.firstEntry.text())+phase1*self.PHASE1STEP
        newZero = (self.zeroVal-(value-self.firstVal)*self.refVal/self.father.current.sw)
        self.zeroVal = np.mod(newZero+180,360)-180
        self.zeroEntry.setText('%.3f' % self.zeroVal)
        self.firstVal = value
        self.firstEntry.setText('%.3f' % self.firstVal)
        self.available = False
        self.zeroScale.setValue(round(self.zeroVal/180.0*self.RESOLUTION))
        self.firstScale.setValue(round(self.firstVal/self.P1LIMIT*self.RESOLUTION))
        self.available = True
        self.father.current.setPhaseInter(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0)

    def inputRef(self, *args): 
        self.refVal = safeEval(self.refEntry.text())
        self.refEntry.setText('%.3f' % self.refVal)

    def setRef(self,value,*args):
        self.refVal = float(value)
        self.refEntry.setText('%.3f' % self.refVal)

    def pickRef(self, *args): 
        self.father.current.peakPickFunc = lambda pos,self=self: self.setRef(self.father.current.xax[pos[0]])
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
        self.father.undoList.append(self.father.current.applyPhase(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0,(self.singleSlice.isChecked()==1)))
        self.father.menuEnable()
        self.deleteLater()

################################################################
class ApodWindow(QtGui.QWidget):
    
    RESOLUTION = 10000
    
    def __init__(self, parent):
        parent.menuDisable()
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.entries = []
        self.ticks = []
        self.maximum = 100.0*self.father.current.sw/(self.father.current.data1D.shape[-1])
        self.lorstep = 1.0
        self.gaussstep = 1.0
        self.available = True
        self.setWindowTitle("Apodize")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        lorTick = QtGui.QCheckBox("Lorentzian:")
        lorTick.toggled.connect(lambda: self.checkEval(0))
        grid.addWidget(lorTick,0,0,1,3)
        self.ticks.append(lorTick)
        lorEntry = QtGui.QLineEdit()
        lorEntry.setEnabled(False)
        lorEntry.setAlignment(QtCore.Qt.AlignHCenter)
        lorEntry.setText("0.00")
        lorEntry.returnPressed.connect(self.apodPreview)
        grid.addWidget(lorEntry,1,1)
        self.entries.append(lorEntry)
        leftLor = QtGui.QPushButton("<")
        leftLor.clicked.connect(lambda:self.stepLB(-0.5*self.father.current.sw/(self.father.current.data1D.shape[-1]),0))
        leftLor.setAutoRepeat(True)
        grid.addWidget(leftLor,1,0)
        rightLor = QtGui.QPushButton(">")
        rightLor.clicked.connect(lambda:self.stepLB(0.5*self.father.current.sw/(self.father.current.data1D.shape[-1]),0))
        rightLor.setAutoRepeat(True)
        grid.addWidget(rightLor,1,2)
        self.lorScale = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.lorScale.setRange(0,self.RESOLUTION)
        self.lorScale.valueChanged.connect(self.setLor)
        grid.addWidget(self.lorScale,2,0,1,3)
        self.lorMax = 100.0*self.father.current.sw/(self.father.current.data1D.shape[-1])

        gaussTick = QtGui.QCheckBox("Gaussian:")
        gaussTick.toggled.connect(lambda: self.checkEval(1))
        grid.addWidget(gaussTick,3,0,1,3)
        self.ticks.append(gaussTick)
        gaussEntry = QtGui.QLineEdit()
        gaussEntry.setEnabled(False)
        gaussEntry.setAlignment(QtCore.Qt.AlignHCenter)
        gaussEntry.setText("0.00")
        gaussEntry.returnPressed.connect(self.apodPreview)
        grid.addWidget(gaussEntry,4,1)
        self.entries.append(gaussEntry)
        leftGauss = QtGui.QPushButton("<")
        leftGauss.clicked.connect(lambda:self.stepLB(0,-0.5*self.father.current.sw/(self.father.current.data1D.shape[-1])))
        leftGauss.setAutoRepeat(True)
        grid.addWidget(leftGauss,4,0)
        rightGauss = QtGui.QPushButton(">")
        rightGauss.clicked.connect(lambda:self.stepLB(0,0.5*self.father.current.sw/(self.father.current.data1D.shape[-1])))
        rightGauss.setAutoRepeat(True)
        grid.addWidget(rightGauss,4,2)
        self.gaussScale = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.gaussScale.setRange(0,self.RESOLUTION)
        self.gaussScale.valueChanged.connect(self.setGauss)
        grid.addWidget(self.gaussScale,5,0,1,3)
        self.gaussMax = 100.0*self.father.current.sw/(self.father.current.data1D.shape[-1])

        cos2Tick = QtGui.QCheckBox("Cos^2:")
        cos2Tick.clicked.connect(lambda: self.checkEval(2))
        grid.addWidget(cos2Tick,6,0,1,3)
        self.ticks.append(cos2Tick)
        cos2Entry = QtGui.QLineEdit()
        cos2Entry.setEnabled(False)
        cos2Entry.setAlignment(QtCore.Qt.AlignHCenter)
        cos2Entry.setText("1.00")
        cos2Entry.returnPressed.connect(self.apodPreview)
        grid.addWidget(cos2Entry,7,1)
        self.entries.append(cos2Entry)

        hammingTick = QtGui.QCheckBox("Hamming:")
        hammingTick.clicked.connect(lambda: self.checkEval(3))
        grid.addWidget(hammingTick,8,0,1,3)
        self.ticks.append(hammingTick)
        hammingEntry = QtGui.QLineEdit()
        hammingEntry.setEnabled(False)
        hammingEntry.setAlignment(QtCore.Qt.AlignHCenter)
        hammingEntry.setText("1.00")
        hammingEntry.returnPressed.connect(self.apodPreview)
        grid.addWidget(hammingEntry,9,1)
        self.entries.append(hammingEntry)

        grid.addWidget(QLabel("Shift:"),10,0,1,3)
        self.shiftEntry = QtGui.QLineEdit()
        self.shiftEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.shiftEntry.setText("0.00")
        self.shiftEntry.returnPressed.connect(self.apodPreview)
        grid.addWidget(self.shiftEntry,11,1)

        if self.father.current.data.dim > 1:
            grid.addWidget(QLabel("Shifting:"),12,0,1,3)
            self.shiftingEntry = QtGui.QLineEdit()
            self.shiftingEntry.setAlignment(QtCore.Qt.AlignHCenter)
            self.shiftingEntry.setText("0.00")
            self.shiftingEntry.returnPressed.connect(self.apodPreview)
            grid.addWidget(self.shiftingEntry,13,1)
            self.shiftingAxes = QtGui.QComboBox()
            self.shiftingValues = list(map(str,np.delete(range(1,self.father.current.data.dim+1),self.father.current.axes)))
            self.shiftingAxes.addItems(self.shiftingValues)
            self.shiftingAxes.currentIndexChanged.connect(self.apodPreview)
            grid.addWidget(self.shiftingAxes,14,1)
            
        self.singleSlice = QtGui.QCheckBox("Single slice")
        grid.addWidget(self.singleSlice,15,0,1,3)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,1,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,1,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)

    def checkEval(self,num):
        if self.ticks[num].isChecked():
            self.entries[num].setEnabled(True)
        else:
            self.entries[num].setEnabled(False)
        self.apodPreview()

    def setLor(self,value, *args):
        if self.available:
            self.entries[0].setText('%.2f' % (float(value)*self.maximum/self.RESOLUTION))
            if not self.ticks[0].isChecked():
                self.ticks[0].setChecked(1)
            self.apodPreview()

    def setGauss(self,value, *args):
        if self.available:
            self.entries[1].setText('%.2f' % (float(value)*self.maximum/self.RESOLUTION))
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
            self.entries[0].setText('%.2f' % lor)
            self.lorScale.setValue(round(lor*self.RESOLUTION/self.maximum))
        if self.ticks[1].isChecked():
            gauss = safeEval(self.entries[1].text())
            self.entries[1].setText('%.2f' % gauss)
            self.gaussScale.setValue(round(gauss*self.RESOLUTION/self.maximum))
        if self.ticks[2].isChecked():
            cos2 = safeEval(self.entries[2].text())
            self.entries[2].setText('%.2f' % cos2)
        if self.ticks[3].isChecked():
            hamming = safeEval(self.entries[3].text())
            self.entries[3].setText('%.2f' % hamming)
        shift = safeEval(self.shiftEntry.text())
        self.shiftEntry.setText('%.2f' % shift)
        if self.father.current.data.dim > 1:
            shifting = safeEval(self.shiftingEntry.text())
            self.shiftingEntry.setText('%.2f' % shifting)
            shiftingAxes = int(self.shiftingValues[self.shiftingAxes.currentIndex()])-1
        else:
            shiftingAxes = None
        self.available = True
        self.father.current.apodPreview(lor,gauss,cos2,hamming,shift,shifting,shiftingAxes)

    def stepLB(self,lorincr,gaussincr):
        self.entries[0].setText('%.2f' %(float(self.entries[0].text())+lorincr*self.lorstep))
        self.entries[1].setText('%.2f' %(float(self.entries[1].text())+gaussincr*self.gaussstep))
        if (lorincr!=0) and (not self.ticks[0].isChecked()):
            self.ticks[0].setChecked(1)
        if (gaussincr!=0) and (not self.ticks[1].isChecked()):
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
        hamming=None
        shifting = None
        shiftingAxes = 0
        if self.ticks[0].isChecked():
            lor = safeEval(self.entries[0].text())
        if self.ticks[1].isChecked():
            gauss = safeEval(self.entries[1].text())
        if self.ticks[2].isChecked():
            cos2 = safeEval(self.entries[2].text())
        if self.ticks[3].isChecked():
            hamming = safeEval(self.entries[3].text())
        shift = safeEval(self.shiftEntry.text())
        if self.father.current.data.dim > 1:
            shifting = safeEval(self.shiftingEntry.text())
            shiftingAxes = int(self.shiftingValues[self.shiftingAxes.currentIndex()])-1
        else:
            shiftingAxes = None
        self.father.redoList = []
        self.father.undoList.append(self.father.current.applyApod(lor,gauss,cos2,hamming,shift,shifting,shiftingAxes,(self.singleSlice.isChecked())))
        self.father.menuEnable()
        self.deleteLater()

#######################################################################################
class SizeWindow(QtGui.QWidget): 
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        parent.menuDisable()
        self.father = parent
        self.setWindowTitle("Set size")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Size:"),0,0)
        self.sizeVal = parent.current.data1D.shape[-1]
        self.sizeEntry = QtGui.QLineEdit()
        self.sizeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.sizeEntry.setText(str(self.sizeVal))
        self.sizeEntry.returnPressed.connect(self.sizePreview)
        grid.addWidget(self.sizeEntry,1,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,1,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,1,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)
 
    def sizePreview(self, *args):
        inp = safeEval(self.sizeEntry.text())
        if inp is not None:
            self.sizeVal = int(round(inp))
        if self.sizeVal < 1:
            self.sizeVal = 1
        self.sizeEntry.setText(str(self.sizeVal))
        self.father.current.setSizePreview(self.sizeVal)

    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.plotReset()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        inp = safeEval(self.sizeEntry.text())
        if inp is not None:
            self.sizeVal = int(round(inp))
        if self.sizeVal < 1:
            self.sizeVal = 1
        self.father.redoList = []
        self.father.undoList.append(self.father.current.applySize(self.sizeVal))
        self.father.sideframe.upd()
        self.father.menuEnable()
        self.deleteLater()
        
##########################################################################################
class SwapEchoWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Swap echo")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Echo position:"),0,0)
        self.posVal = int(round(0.5*len(parent.current.data1D)))
        self.posEntry = QtGui.QLineEdit()
        self.posEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.posEntry.setText(str(self.posVal))
        self.posEntry.returnPressed.connect(self.swapEchoPreview)
        grid.addWidget(self.posEntry,1,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,1,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,1,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.father.current.peakPick = True
        self.setGeometry(0,0,0,0)
 
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
        self.father.current.plotReset()
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
            print("not a valid index for swap echo")
        
    def picked(self,pos):
        self.father.current.setSwapEchoPreview(pos[0])
        self.posEntry.setText(str(pos[0]))
        self.father.current.peakPick = False

###########################################################################
class ShiftDataWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Shifting data")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Data points to shift:"),0,0,1,3)
        self.shiftVal = 0
        self.shiftEntry = QtGui.QLineEdit()
        self.shiftEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.shiftEntry.setText(str(self.shiftVal))
        self.shiftEntry.returnPressed.connect(self.shiftPreview)
        grid.addWidget(self.shiftEntry,1,1)
        leftShift = QtGui.QPushButton("<")
        leftShift.clicked.connect(self.stepDownShift)
        leftShift.setAutoRepeat(True)
        grid.addWidget(leftShift,1,0)
        rightShift = QtGui.QPushButton(">")
        rightShift.clicked.connect(self.stepUpShift)
        rightShift.setAutoRepeat(True)
        grid.addWidget(rightShift,1,2)
        self.singleSlice = QtGui.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice,1,0,1,2)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)

    def stepUpShift(self, *args):
        inp = safeEval(self.shiftEntry.text())
        if inp is not None:
            self.shiftVal = int(round(inp))
        self.shiftVal = self.shiftVal+1
        self.shiftEntry.setText(str(self.shiftVal))
        self.shiftPreview()

    def stepDownShift(self, *args):
        inp = safeEval(self.shiftEntry.text())
        if inp is not None:
            self.shiftVal = int(round(inp))
        self.shiftVal = self.shiftVal-1
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
        self.father.current.plotReset()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        inp = safeEval(self.shiftEntry.text())
        if inp is None:
            print("Not a valid input")
            return
        shift = int(round(inp))
        self.father.redoList = []
        self.father.undoList.append(self.father.current.applyShift(shift,(self.singleSlice.isChecked())))
        self.father.menuEnable()
        self.deleteLater()

#############################################################
class DCWindow(QtGui.QWidget): 
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Offset correction")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        self.startVal = int(round(0.8*parent.current.data1D.shape[-1]))
        self.endVal = parent.current.data1D.shape[-1]
        grid.addWidget(QLabel("Start point:"),0,0)
        self.startEntry = QtGui.QLineEdit()
        self.startEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.startEntry.setText(str(self.startVal))
        self.startEntry.returnPressed.connect(self.offsetPreview)
        grid.addWidget(self.startEntry,1,0)
        grid.addWidget(QLabel("End point:"),2,0)
        self.endEntry = QtGui.QLineEdit()
        self.endEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.endEntry.setText(str(self.endVal))
        self.endEntry.returnPressed.connect(self.offsetPreview)
        grid.addWidget(self.endEntry,3,0)
        grid.addWidget(QLabel("Offset:"),4,0)
        self.offsetEntry = QtGui.QLineEdit()
        self.offsetEntry.setAlignment(QtCore.Qt.AlignHCenter)
        val = parent.current.getdcOffset(int(round(0.8*parent.current.data1D.shape[-1])),parent.current.data1D.shape[-1])
        self.offsetEntry.setText('{:.2e}'.format(val))
        self.offsetEntry.returnPressed.connect(lambda arg: self.offsetPreview(arg,True))
        grid.addWidget(self.offsetEntry,5,0)
        self.singleSlice = QtGui.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice,1,0,1,2)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.father.current.peakPick = True
        self.setGeometry(0,0,0,0)

    def picked(self,pos,second=False): 
        dataLength = self.father.current.data1D.shape[-1]
        if second:
            inp = safeEval(self.startEntry.text())
            if inp is not None:
                self.startVal=int(round(inp))
            if self.startVal < 0:
                self.startVal = 0
            elif self.startVal > dataLength:
                self.startVal = dataLength
            self.startEntry.setText(str(self.startVal))
            self.endVal=pos[0]
            self.endEntry.setText(str(self.endVal))
            val = self.father.current.getdcOffset(self.startVal,self.endVal)
            self.offsetEntry.setText('{:.2e}'.format(val))
            self.father.current.dcOffset(val)
            self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
            self.father.current.peakPick = True
        else:
            self.startEntry.setText(str(pos[0]))
            inp = safeEval(self.endEntry.text())
            if inp is not None:
                self.endVal=int(round(inp))
            if self.endVal < 0:
                self.endVal = 0
            elif self.endVal > dataLength:
                self.endVal = dataLength
            self.startVal = pos[0]
            val = self.father.current.getdcOffset(self.startVal,self.endVal)
            self.offsetEntry.setText('{:.2e}'.format(val))
            self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos,True) 
            self.father.current.peakPick = True

    def offsetPreview(self, arg, inserted=False):
        if inserted:
            self.father.current.dcOffset(safeEval(self.offsetEntry.text()))
        else:
            dataLength = self.father.current.data1D.shape[-1]
            inp = safeEval(self.startEntry.text())
            if inp is not None:
                self.startVal=int(round(inp))
            if self.startVal < 0:
                self.startVal = 0
            elif self.startVal > dataLength:
                self.startVal = dataLength
            self.minVal.setText(str(self.startVal))
            inp = safeEval(self.endEntry.text())
            if inp is not None:
                self.endVal=int(round(inp))
            if self.endVal < 0:
                self.endVal = 0
            elif self.endVal > dataLength:
                self.endVal = dataLength
            self.maxEntry.setText(str(self.endVal))
            val = self.father.current.getdcOffset(self.startVal,self.endVal)
            self.offsetEntry.setText('{:.2e}'.format(val))
            self.father.current.dcOffset(val)

    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.current.upd()
        self.father.current.plotReset()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        inp = safeEval(self.offsetEntry.text())
        if inp is None:
            print("Not a valid value")
            return
        self.father.current.peakPickReset()
        self.father.redoList = []
        self.father.undoList.append(self.father.current.subtract(inp,self.singleSlice.isChecked()))
        self.father.menuEnable()
        self.deleteLater()

#############################################################
class BaselineWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Baseline correction")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Polynomial Degree:"),0,0,1,2)
        self.removeList = []
        self.degree = 3
        self.degreeEntry = QtGui.QLineEdit()
        self.degreeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.degreeEntry.setText(str(self.degree))
        self.degreeEntry.returnPressed.connect(self.setDegree)
        grid.addWidget(self.degreeEntry,1,0,1,2)
        resetButton = QtGui.QPushButton("&Reset")
        resetButton.clicked.connect(self.reset)
        grid.addWidget(resetButton,2,0)
        fitButton = QtGui.QPushButton("&Fit")
        fitButton.clicked.connect(self.preview)
        grid.addWidget(fitButton,2,1)
        self.singleSlice = QtGui.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice,1,0,1,2)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.father.current.peakPick = True
        self.setGeometry(0,0,0,0)

    def picked(self,pos):
        self.removeList.append(pos[0])
        self.father.current.previewRemoveList(self.removeList)
        self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
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
        self.father.current.previewBaseline(self.degree,self.removeList)
        self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
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
            print("Not a valid value")
            return
        self.father.current.peakPickReset()
        self.father.current.resetPreviewRemoveList()
        self.father.redoList = []
        self.father.undoList.append(self.father.current.applyBaseline(self.degree,self.removeList,self.singleSlice.isChecked()))
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

#############################################################
class regionWindow(QtGui.QWidget): 
    def __init__(self, parent,name):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle(name)
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        self.startVal = 0
        self.endVal = parent.current.data1D.shape[-1]
        grid.addWidget(QLabel("Start point:"),0,0)
        self.startEntry = QtGui.QLineEdit()
        self.startEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.startEntry.setText(str(self.startVal))
        self.startEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.startEntry,1,0)
        grid.addWidget(QLabel("End point:"),2,0)
        self.endEntry = QtGui.QLineEdit()
        self.endEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.endEntry.setText(str(self.endVal))
        self.endEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.endEntry,3,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.father.current.peakPick = True
        self.setGeometry(0,0,0,0)

    def picked(self,pos,second=False): 
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
            self.endVal=pos[0]
            self.endEntry.setText(str(self.endVal))
            self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
            self.father.current.peakPick = True
        else:
            self.startEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos,True) 
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
        self.maxEntry.setText(str(self.endVal))
        
    def apply(self,maximum,minimum):
        pass
    
    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.updAllFrames()
        self.father.menuEnable()
        self.deleteLater()
        
    def applyAndClose(self):
        self.father.current.peakPickReset()
        dataLength = self.father.current.data1D.shape[-1]
        inp = safeEval(self.startEntry.text())
        if inp is None:
            print("Not a valid value")
            return
        self.startVal = int(round(inp))
        if self.startVal < 0:
            self.startVal = 0
        elif self.startVal > dataLength:
            self.startVal = dataLength
        inp = safeEval(self.endEntry.text())
        if inp is None:
            print("Not a valid value")
            return
        self.endVal = int(round(inp))
        if self.endVal < 0:
            self.endVal = 0
        elif self.endVal > dataLength:
            self.endVal = dataLength
        self.apply(self.startVal,self.endVal)
        self.father.menuEnable()
        self.deleteLater()

############################################################
class integrateWindow(regionWindow): 
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Integrate')

    def apply(self,maximum,minimum):
        self.father.redoList = []
        self.father.undoList.append(self.father.current.integrate(minimum,maximum))
        self.father.updAllFrames()

############################################################
class sumWindow(regionWindow): 
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Sum')

    def apply(self,maximum,minimum):
        self.father.redoList = []
        self.father.undoList.append(self.father.current.sum(minimum,maximum))
        self.father.updAllFrames()
        
############################################################
class maxWindow(regionWindow): 
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Max')

    def apply(self,maximum,minimum):
        self.father.redoList = []
        self.father.undoList.append(self.father.current.maxMatrix(minimum,maximum))
        self.father.updAllFrames()

############################################################
class minWindow(regionWindow): 
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Min')

    def apply(self,maximum,minimum):
        self.father.redoList = []
        self.father.undoList.append(self.father.current.minMatrix(minimum,maximum))
        self.father.updAllFrames()
        
############################################################
class argmaxWindow(regionWindow):
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Max position')

    def apply(self,maximum,minimum):
        self.father.redoList = []
        self.father.undoList.append(self.father.current.argmaxMatrix(minimum,maximum))
        self.father.updAllFrames()

############################################################
class argminWindow(regionWindow):
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Min position')

    def apply(self,maximum,minimum):
        self.father.redoList = []
        self.father.undoList.append(self.father.current.argminMatrix(minimum,maximum))
        self.father.updAllFrames()
        
############################################################
class extractRegionWindow(regionWindow):
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Extract part')

    def apply(self,maximum,minimum):
        self.father.redoList = []
        self.father.undoList.append(self.father.current.getRegion(minimum,maximum))
        self.father.updAllFrames()

##############################################################
class DeleteWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Delete")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Indexes to delete:"),0,0)
        self.delEntry = QtGui.QLineEdit()
        self.delEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.delEntry.setText(str('0'))
        self.delEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.delEntry,1,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)
        
    def preview(self, *args):
        env = vars(np).copy()
        length = int(self.father.current.data1D.shape[-1])
        env['length']=length # so length can be used to in equations
        pos=np.array(eval(self.delEntry.text(),env)).flatten()                # find a better solution, also add catch for exceptions
        pos[pos<0]=pos[pos<0]+length
        if (pos > -1).all() and (pos < length).all():
            self.father.current.deletePreview(pos)
        else:
            print('Not all values are valid indexes to delete')
        
    def applyAndClose(self):
        env = vars(np).copy()
        length = int(self.father.current.data1D.shape[-1])
        env['length']=length # so length can be used to in equations
        pos=np.array(eval(self.delEntry.text(),env)).flatten()                # find a better solution, also add catch for exceptions
        pos[pos<0]=pos[pos<0]+length
        if (pos > -1).all() and (pos < length).all():
            self.father.redoList = []
            self.father.undoList.append(self.father.current.delete(pos))
            self.father.menuEnable()
            self.father.sideframe.upd()
            self.deleteLater()
        else:
            print('Not all values are valid indexes to delete')
        
    def closeEvent(self, *args):
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

##############################################################
class SplitWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Split")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Sections:"),0,0)
        self.splitEntry = QtGui.QLineEdit()
        self.splitEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.splitEntry.setText('1')
        self.splitEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.splitEntry,1,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)

    def preview(self, *args):
        val = safeEval(self.splitEntry.text())
        if val is not None:
            self.splitEntry.setText(str(int(round(val))))
        
    def applyAndClose(self):
        val = safeEval(self.splitEntry.text())
        if val is None:
            print("Not a valid value")
            return
        self.father.redoList = []
        self.father.undoList.append(self.father.current.split(int(round(val))))
        self.father.menuEnable()
        self.father.updAllFrames()
        self.deleteLater()
        
    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()
        
##############################################################
class ConcatenateWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Concatenate")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Concatenation axes:"),0,0)
        self.axesEntry = QtGui.QComboBox()
        self.axesEntry.addItems(np.array(np.arange(self.father.current.data.dim-1)+1,dtype=str))
        grid.addWidget(self.axesEntry,1,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)
        
    def applyAndClose(self):
        self.father.redoList = []
        self.father.undoList.append(self.father.current.concatenate(self.axesEntry.currentIndex()))
        self.father.menuEnable()
        self.father.updAllFrames()
        self.deleteLater()
        
    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

##############################################################
class InsertWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Insert")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Start insert at index:"),0,0)
        self.posEntry = QtGui.QLineEdit()
        self.posEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.posEntry.setText(str(self.father.current.data1D.shape[-1]))
        self.posEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.posEntry,1,0)
        grid.addWidget(QLabel("Workspace to insert:"),2,0)
        self.wsEntry = QtGui.QComboBox()
        self.wsEntry.addItems(self.father.mainProgram.workspaceNames)
        grid.addWidget(self.wsEntry,3,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)

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
            print("Not a valid value")
            return
        pos = int(round(pos))
        if pos > self.father.current.data1D.shape[-1]:
            pos = self.father.current.data1D.shape[-1]
        elif pos < 0:
            pos = 0
        ws = self.wsEntry.currentIndex()
        self.father.redoList = []
        self.father.undoList.append(self.father.current.insert(self.father.mainProgram.workspaces[ws].masterData.data,pos))
        self.father.menuEnable()
        self.father.sideframe.upd()
        self.deleteLater()
        
    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

##############################################################
class AddWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Add")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Workspace to add:"),0,0)
        self.wsEntry = QtGui.QComboBox()
        self.wsEntry.addItems(self.father.mainProgram.workspaceNames)
        grid.addWidget(self.wsEntry,1,0)
        self.singleSlice = QtGui.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice,1,0,1,2)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)

    def applyAndClose(self):
        ws = self.wsEntry.currentIndex()
        self.father.redoList = []
        self.father.undoList.append(self.father.current.add(self.father.mainProgram.workspaces[ws].masterData.data,self.singleSlice.isChecked()))
        self.father.menuEnable()
        self.father.sideframe.upd()
        self.deleteLater()
        
    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

##############################################################
class SubtractWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Subtract")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Workspace to subtract:"),0,0)
        self.wsEntry = QtGui.QComboBox()
        self.wsEntry.addItems(self.father.mainProgram.workspaceNames)
        grid.addWidget(self.wsEntry,1,0)
        self.singleSlice = QtGui.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice,1,0,1,2)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)

    def applyAndClose(self):
        ws = self.wsEntry.currentIndex()
        self.father.redoList = []
        self.father.undoList.append(self.father.current.subtract(self.father.mainProgram.workspaces[ws].masterData.data,self.singleSlice.isChecked()))
        self.father.menuEnable()
        self.father.sideframe.upd()
        self.deleteLater()
        
    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

##############################################################
class SNWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Signal to noise")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Start point noise:"),0,0)
        self.minNoiseEntry = QtGui.QLineEdit()
        self.minNoiseEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.minNoiseEntry.setText("0")
        self.minNoiseEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.minNoiseEntry,1,0)
        grid.addWidget(QLabel("End point noise:"),2,0)
        self.maxNoiseEntry = QtGui.QLineEdit()
        self.maxNoiseEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.maxNoiseEntry.setText(str(parent.current.data1D.shape[-1]))
        self.maxNoiseEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.maxNoiseEntry,3,0)
        grid.addWidget(QLabel("Start point signal:"),4,0)
        self.minEntry = QtGui.QLineEdit()
        self.minEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.minEntry.setText("0")
        self.minEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.minEntry,5,0)
        grid.addWidget(QLabel("End point signal:"),6,0)
        self.maxEntry = QtGui.QLineEdit()
        self.maxEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.maxEntry.setText(str(parent.current.data1D.shape[-1]))
        self.maxEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.maxEntry,7,0)
        grid.addWidget(QLabel("S/N:"),8,0)
        self.snEntry = QtGui.QLineEdit()
        self.snEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.snEntry.setText('0.0')
        grid.addWidget(self.snEntry,9,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Fit")
        okButton.clicked.connect(self.apply)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.father.current.peakPick = True
        self.setGeometry(0,0,0,0)
        
    def picked(self,pos,num=0): 
        if num == 0:
            self.minNoiseEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos,1) 
            self.father.current.peakPick = True
        elif num == 1:
            self.maxNoiseEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos,2) 
            self.father.current.peakPick = True
        elif num == 2:
            self.minEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos,3) 
            self.father.current.peakPick = True
        elif num == 3:
            self.maxEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos,0) 
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
            print("Not a valid value")
            return
        minimumNoise = int(round(inp))
        if minimumNoise < 0:
            minimumNoise = 0
        elif minimumNoise > dataLength:
            minimumNoise = dataLength
        self.minNoiseEntry.setText(str(minimumNoise))
        inp = safeEval(self.maxNoiseEntry.text())
        if inp is None:
            print("Not a valid value")
            return
        maximumNoise = int(round(inp))
        if maximumNoise < 0:
            maximumNoise = 0
        elif maximumNoise > dataLength:
            maximumNoise = dataLength
        self.maxNoiseEntry.setText(str(maximumNoise))
        inp = safeEval(self.minEntry.text())
        if inp is None:
            print("Not a valid value")
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            print("Not a valid value")
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.snEntry.setText(str(self.father.current.SN(minimumNoise,maximumNoise,minimum,maximum)))
        
    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.menuEnable()
        self.deleteLater()

##############################################################
class FWHMWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("FWHM")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Start point:"),0,0)
        self.minEntry = QtGui.QLineEdit()
        self.minEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.minEntry.setText("0")
        self.minEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.minEntry,1,0)
        grid.addWidget(QLabel("End point:"),2,0)
        self.maxEntry = QtGui.QLineEdit()
        self.maxEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.maxEntry.setText(str(parent.current.data1D.shape[-1]))
        self.maxEntry.returnPressed.connect(self.checkValues)
        grid.addWidget(self.maxEntry,3,0)
        if self.father.current.spec == 1:
            if self.father.current.axType == 0:
                grid.addWidget(QLabel("FWHM [Hz]:"),4,0)
            elif self.father.current.axType == 1:
                grid.addWidget(QLabel("FWHM [kHz]:"),4,0)
            elif self.father.current.axType == 2:
                grid.addWidget(QLabel("FWHM [MHz]:"),4,0)
            elif self.father.current.axType == 3:
                grid.addWidget(QLabel("FWHM [ppm]:"),4,0)
        else:
            if self.father.current.axType == 0:
                grid.addWidget(QLabel("FWHM [s]:"),4,0)
            elif self.father.current.axType == 1:
                grid.addWidget(QLabel("FWHM [ms]:"),4,0)
            elif self.father.current.axType == 2:
                grid.addWidget(QLabel(u"FWHM [\u03bcs]:"),4,0)
        self.fwhmEntry = QtGui.QLineEdit()
        self.fwhmEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.fwhmEntry.setText('0.0')
        grid.addWidget(self.fwhmEntry,5,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Fit")
        okButton.clicked.connect(self.apply)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.father.current.peakPick = True
        self.setGeometry(0,0,0,0)
        
    def picked(self,pos,num=0): 
        if num == 0:
            self.minEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos,1) 
            self.father.current.peakPick = True
        elif num == 1:
            self.maxEntry.setText(str(pos[0]))
            self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos,0) 
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
            print("Not a valid value")
            return
        minimum = int(round(inp))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minEntry.setText(str(minimum))
        inp = safeEval(self.maxEntry.text())
        if inp is None:
            print("Not a valid value")
            return
        maximum = int(round(inp))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxEntry.setText(str(maximum))
        self.fwhmEntry.setText(str(self.father.current.fwhm(minimum,maximum)))
        
    def closeEvent(self, *args):
        self.father.current.peakPickReset()
        self.father.menuEnable()
        self.deleteLater()        

################################################################
class ShearingWindow(QtGui.QWidget): 
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Shearing")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        options = list(map(str,range(1,self.father.masterData.dim+1)))
        grid.addWidget(QLabel("Shearing constant:"),0,0)
        self.shearEntry = QtGui.QLineEdit()
        self.shearEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.shearEntry.setText("0.0")
        self.shearEntry.returnPressed.connect(self.shearPreview)
        grid.addWidget(self.shearEntry,1,0)
        grid.addWidget(QLabel("Shearing direction:"),2,0)
        self.dirEntry = QtGui.QComboBox()
        self.dirEntry.addItems(options)
        self.dirEntry.setCurrentIndex(self.father.masterData.dim-2)
        grid.addWidget(self.dirEntry,3,0)
        grid.addWidget(QLabel("Shearing axis:"),4,0)
        self.axEntry = QtGui.QComboBox()
        self.axEntry.addItems(options)
        self.axEntry.setCurrentIndex(self.father.masterData.dim-1)
        grid.addWidget(self.axEntry,5,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)

    def shearPreview(self, *args):
        shear = safeEval(self.shearEntry.text())
        if inp is not None:
            self.shear.set(str(float(shear)))

    def closeEvent(self, *args):
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        shear = safeEval(self.shearEntry.text())
        if inp is None:
            print("Not a valid value")
            return
        axes = self.dirEntry.currentIndex()
        axes2 = self.axEntry.currentIndex()
        if axes == axes2:
            print("Axes can't be the same for shearing")
            return
        else:
            self.father.redoList = []
            self.father.undoList.append(self.father.current.shearing(float(shear),axes,axes2))
            self.father.menuEnable()
            self.deleteLater()

##########################################################################################
class MultiplyWindow(QtGui.QWidget): 
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Multiply")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Values:"),0,0)
        self.valEntry = QtGui.QLineEdit()
        self.valEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.valEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.valEntry,1,0)
        self.singleSlice = QtGui.QCheckBox("Single slice")
        layout.addWidget(self.singleSlice,1,0,1,2)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)

    def preview(self, *args):
        env = vars(np).copy()
        env['length']=int(self.father.current.data1D.shape[-1]) # so length can be used to in equations
        env['euro']=lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal,num)
        val=eval(self.valEntry.text(),env)                # find a better solution, also add catch for exceptions          
        self.father.current.multiplyPreview(np.array(val))

    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.plotReset()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        env = vars(np).copy()
        env['length']=int(self.father.current.data1D.shape[-1]) # so length can be used to in equations
        env['euro']=lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal,num)
        val=eval(self.valEntry.text(),env)                # find a better solution, also add catch for exceptions
        self.father.redoList = []
        self.father.undoList.append(self.father.current.multiply(np.array(val),self.singleSlice.isChecked()))
        self.father.menuEnable()
        self.deleteLater()

##########################################################################################
class XaxWindow(QtGui.QWidget): 
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("User defined x-axis")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        grid.addWidget(QLabel("Expression for x-axis values:"),0,0)
        self.valEntry = QtGui.QLineEdit()
        self.valEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.valEntry.returnPressed.connect(self.xaxPreview)
        grid.addWidget(self.valEntry,1,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(0,0,0,0)
        
    def xaxPreview(self, *args):
        env = vars(np).copy()
        env['length']=int(self.father.current.data1D.shape[-1]) # so length can be used to in equations
        env['euro']=lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal,num)
        val=eval(self.valEntry.text(),env)                # find a better solution, also add catch for exceptions          
        if not isinstance(val,(list,np.ndarray)):
            print("Input is not a list or array")
            return
        if len(val)!=self.father.current.data1D.shape[-1]:
            print("Length of input does not match length of data")
            return
        if not all(isinstance(x,(int,float)) for x in val):
            print("Array is not all of int or float type")
            return
        self.father.current.setXaxPreview(np.array(val))

    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.plotReset()
        self.father.current.showFid()
        self.father.menuEnable()
        self.deleteLater()

    def applyAndClose(self):
        env = vars(np).copy()
        env['length']=int(self.father.current.data1D.shape[-1]) # so length can be used to in equations
        env['euro']=lambda fVal, num=int(self.father.current.data1D.shape[-1]): euro(fVal,num)
        val=eval(self.valEntry.text(),env)                # find a better solution, also add catch for exceptions
        if not isinstance(val,(list,np.ndarray)):
            print("Input is not a list or array")
            return
        if len(val)!=self.father.current.data1D.shape[-1]:
            print("Length of input does not match length of data")
            return
        if not all(isinstance(x,(int,float)) for x in val):
            print("Array is not all of int or float type")
            return
        self.father.redoList = []
        self.father.undoList.append(self.father.current.setXax(np.array(val)))
        self.father.menuEnable()
        self.deleteLater()

##########################################################################################
class RefWindow(QtGui.QWidget): 
    def __init__(self, parent):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window| QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Reference")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        if parent.current.spec == 0:
            print('Setting ppm is only available for frequency data')
            self.deleteLater()
            return
        grid.addWidget(QLabel("Frequency [MHz]:"),0,0)
        self.freqEntry = QtGui.QLineEdit()
        self.freqEntry.setText("%.7f" % (self.father.current.ref*1e-6))
        self.freqEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.freqEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.freqEntry,1,0)
        grid.addWidget(QLabel("Reference [ppm]:"),2,0)
        self.refEntry = QtGui.QLineEdit()
        self.refEntry.setText("0.0")
        self.refEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.refEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.refEntry,3,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,2,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton,2,1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.father.current.peakPick = True
        self.setGeometry(0,0,0,0)

    def preview(self, *args): 
        freq = safeEval(self.freqEntry.text())
        ref = safeEval(self.refEntry.text())
        if freq is None or ref is None:
            return
        self.freqEntry.setText("%.7f" % (freq))
        self.refEntry.setText(str(ref))

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
            print("Not a valid value")
            return
        freq = freq*1e6
        self.father.redoList = []
        self.father.undoList.append(self.father.current.setRef(freq/(1.0+ref*1e-6)))
        self.father.menuEnable()
        self.deleteLater()
        
    def picked(self,pos): 
        self.freqEntry.setText("%.7f" % ((self.father.current.ref+self.father.current.xax[pos[0]])*1e-6))
        self.father.current.peakPickFunc = lambda pos,self=self: self.picked(pos)
        self.father.current.peakPick = True

root = QtGui.QApplication(sys.argv)
root.setWindowIcon(QtGui.QIcon(os.path.dirname(os.path.realpath(__file__))+'/logo.gif')) 
mainProgram = MainProgram(root)
mainProgram.setWindowTitle("ssNake")
mainProgram.show()
sys.exit(root.exec_())
