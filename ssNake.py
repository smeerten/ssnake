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

import numpy as np
import sys
if sys.version_info >= (3,0):
    from tkinter import *
    import tkinter as tk
    from tkinter.ttk import *
    from tkinter import messagebox as tkMessageBox
    from tkinter.filedialog import askopenfilename
    from tkinter.filedialog import asksaveasfile
    from tkinter.filedialog import asksaveasfilename
    from tkinter.simpledialog import askstring
else:
    from Tkinter import *
    import Tkinter as tk
    from ttk import *
    import tkMessageBox as tkMessageBox
    from tkFileDialog   import askopenfilename
    from tkFileDialog   import asksaveasfile
    from tkFileDialog   import asksaveasfilename
    from tkSimpleDialog import askstring
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib  import rcParams
rcParams['toolbar'] = 'None' 
import spectrum_classes as sc
import fitting as fit
import math
import copy
import os
from struct import unpack
import scipy.io
import json
import copy
import weakref

#------------
from safeEval import safeEval
from euro import euro

pi=math.pi

#one class to rule them all
class MainProgram:
    def __init__(self,root):
        self.root = root
        self.menubar = Menu(self.root)
        self.root.config(menu=self.menubar)
        self.workspaces = []
        self.workspaceNames = []
        self.workspaceNum = 0
        self.workspaceVar = StringVar()
        self.macros = {}
        self.figOrphanage = []
        #the file drop down menu
        self.filemenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        #the hotkeys for different commands
        self.root.protocol("WM_DELETE_WINDOW", self.kill)
        self.root.bind_all("<Control-q>", lambda extra: self.kill())
        self.root.bind_all("<Control-z>", self.undo)
        self.root.bind_all("<Control-y>", self.redo)
        self.root.bind_all("<Control-w>", self.destroyWorkspace)
        self.root.bind_all("<Control-d>", self.duplicateWorkspace)
        self.root.bind_all("<Control-Prior>", lambda args: self.stepWorkspace(-1))
        self.root.bind_all("<Control-Next>", lambda args: self.stepWorkspace(1))
        
        self.filemenu.add_command(label="Load", command=self.autoLoad)

        #the save drop down menu
        self.savemenu = Menu(self.filemenu, tearoff=0)
        self.filemenu.add_cascade(label="Save", menu=self.savemenu)
        self.savemenu.add_command(label="Save figure", command=self.saveFigure)
        self.savemenu.add_command(label="Save JSON", command=self.saveJSONFile)
        self.savemenu.add_command(label="Save MATLAB", command=self.saveMatlabFile)
        self.savemenu.add_command(label="Save as Simpson data", command=self.saveSimpsonFile)

        self.mainWindow = None
        self.workspacemenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Workspaces", menu=self.workspacemenu)
        self.workspacemenu.add_command(label="Duplicate", command=self.duplicateWorkspace)
        self.workspacemenu.add_command(label="Delete", command=self.destroyWorkspace)
        self.workspacemenu.add_command(label="Rename", command=self.renameWorkspace)
        self.activemenu = None

        #the macro menu
        self.macromenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Macros", menu=self.macromenu)
        self.macromenu.add_command(label="Start recording", command=self.macroCreate)
        self.macromenu.add_command(label="Stop recording", command=self.stopMacro)
        self.macrolistmenu = Menu(self.macromenu, tearoff=0)
        self.macromenu.add_cascade(label="Run", menu=self.macrolistmenu)
        self.macrodeletemenu = Menu(self.macromenu, tearoff=0)
        self.macromenu.add_cascade(label="Delete", menu=self.macrodeletemenu)
        self.macrosavemenu = Menu(self.macromenu, tearoff=0)
        self.macromenu.add_cascade(label="Save", menu=self.macrosavemenu)
        self.macromenu.add_command(label="Load", command=self.loadMacro)

        self.filemenu.add_command(label="Exit", command=self.kill)
        self.menuCheck()
        photo = PhotoImage(file=os.path.dirname(os.path.realpath(__file__))+'/logo.gif')
        self.logo = Label(self.root,image=photo)
        self.logo.image = photo # keep a reference!
        self.logo.pack()

    def kill(self):
        answer = tkMessageBox.askquestion("Question", "Are you sure you want to quit?")
        if answer == 'yes':
            self.root.quit()
            self.root.destroy()

    def donateFig(self,fig):
        self.figOrphanage.append(fig)

    def getFig(self):
        if len(self.figOrphanage)==0:
            return Figure()
        else:
            fig = self.figOrphanage.pop()
            fig.clf()
            return fig
            
    def menuCheck(self):
        if self.mainWindow is None:
            self.filemenu.entryconfig('Save',state='disabled')
            self.menubar.entryconfig('Workspaces',state='disabled')
            self.menubar.entryconfig('Macros',state='disabled')
        elif isinstance(self.mainWindow, Main1DWindow):
            self.menubar.entryconfig('Macros',state='normal')
            self.savemenu.entryconfig('Save figure',state='normal')
            if self.mainWindow.currentMacro is None:
                self.macromenu.entryconfig('Stop recording',state='disabled')
                self.macromenu.entryconfig('Start recording',state='normal')
            else:
                self.macromenu.entryconfig('Stop recording',state='normal')
                self.macromenu.entryconfig('Start recording',state='disabled')
            self.filemenu.entryconfig('Save',state='normal')
            self.menubar.entryconfig('Workspaces',state='normal')
        elif isinstance(self.mainWindow, fit.MainPlotWindow):
            self.filemenu.entryconfig('Save',state='normal')
            self.savemenu.entryconfig('Save figure',state='disabled')
            self.menubar.entryconfig('Workspaces',state='normal')
            self.menubar.entryconfig('Macros',state='disabled')
        else:
            self.filemenu.entryconfig('Save',state='normal')
            self.savemenu.entryconfig('Save figure',state='normal')
            self.menubar.entryconfig('Workspaces',state='normal')
            self.menubar.entryconfig('Macros',state='disabled')
            
    def askName(self):
        count = 0
        name = 'spectrum'+str(count)
        while name in self.workspaceNames:
            count += 1
            name = 'spectrum'+str(count)
        givenName = askstring('Spectrum name','Name:',initialvalue=name)
        while (givenName in self.workspaceNames) or givenName is '':
            print('Name exists')
            givenName = askstring('Spectrum name','Name:')
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
        givenName = askstring('Macro name','Name:',initialvalue=name)
        while (givenName in self.macros.keys()) or givenName is '':
            print('Name exists')
            givenName = askstring('Macro name','Name:')
        self.macros[givenName] = []
        self.mainWindow.redoMacro = []
        self.mainWindow.currentMacro = givenName
        self.macrolistmenu.add_command(label=givenName,command=lambda name=givenName: self.runMacro(name))
        self.macrosavemenu.add_command(label=givenName,command=lambda name=givenName: self.saveMacro(name))
        self.macrodeletemenu.add_command(label=givenName,command=lambda name=givenName: self.deleteMacro(name))
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
        fileName=asksaveasfilename(filetypes=(('JSON','.json'),))
        if not fileName:
            return
        with open(fileName,'w') as f:
            json.dump(self.macros[name],f)

    def deleteMacro(self,name):
        self.macrolistmenu.delete(name)
        self.macrosavemenu.delete(name)
        self.macrodeletemenu.delete(name)
        del self.macros[name]
        for i in self.workspaces:
            if i.currentMacro == name:
                i.currentMacro = None
        self.menuCheck()
            
    def loadMacro(self):
        filePath = askopenfilename()
        if len(filePath)==0:
            return
        count = 0
        name = 'macro'+str(count)
        while name in self.macros.keys():
            count += 1
            name = 'macro'+str(count)
        givenName = askstring('Macro name','Name:',initialvalue=name)
        while (givenName in self.macros.keys()) or givenName is '':
            print('Name exists')
            givenName = askstring('Macro name','Name:')
        with open(filePath,'r') as f:
            self.macros[givenName] = json.load(f)
        self.macrolistmenu.add_command(label=givenName,command=lambda name=givenName: self.runMacro(givenName))
        self.macrosavemenu.add_command(label=givenName,command=lambda name=givenName: self.saveMacro(givenName))
        self.menuCheck()
            
    def changeMainWindow(self, var):
        self.logo.pack_forget()
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
        self.workspaces.append(Main1DWindow(self.root,self,copy.deepcopy(self.mainWindow.get_masterData()),self.mainWindow.get_current(),name=name))
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
            self.logo.pack()
            self.updWorkspaceMenu(None)
            
    def updWorkspaceMenu(self,var):
        if self.activemenu is not None:
            self.workspacemenu.delete("Active")
        self.activemenu = Menu(self.workspacemenu, tearoff=0)
        self.workspacemenu.add_cascade(label="Active", menu=self.activemenu)
        self.workspaceVar.set(var)
        for i in self.workspaceNames:
            self.activemenu.add_radiobutton(label=i,variable=self.workspaceVar,value=i,command=lambda i=i: self.changeMainWindow(i))
        self.menuCheck()

    def autoLoad(self):
        filePath = askopenfilename()
        if len(filePath)==0:
            return
        direc = os.path.dirname(filePath)
        filename = os.path.basename(filePath)
        if filename.endswith('.fid') or filename.endswith('.spe'): 
            self.loading(4,filePath)
        elif filename.endswith('.json') or filename.endswith('.JSON'):
            self.loading(5,filePath)
        elif filename.endswith('.mat') or filename.endswith('.MAT'):
            self.loading(6,filePath)
        elif os.path.exists(direc+os.path.sep+'procpar') and os.path.exists(direc+os.path.sep+'fid'):
            self.loading(0,filePath)
        elif os.path.exists(direc+os.path.sep+'acqus') and (os.path.exists(direc+os.path.sep+'fid') or os.path.exists(direc+os.path.sep+'ser')):
            self.loading(1,filePath)
        elif os.path.exists(direc+os.path.sep+'acq') and os.path.exists(direc+os.path.sep+'data'):
            self.loading(2,filePath)
        elif os.path.exists(direc+os.path.sep+'acqu.par'):
            dirFiles = os.listdir(direc)
            files2D = [x for x in dirFiles if '.2d' in x]
            files1D = [x for x in dirFiles if '.1d' in x]
            if len(files2D)!=0 or len(files1D)!=0:
                self.loading(3,filePath)
        else:
            return
        
    def loading(self,num,filePath):
        name = self.askName()
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
            self.workspaces.append(Main1DWindow(self.root,self,masterData,name=name))
            self.workspaceNames.append(name)
            self.changeMainWindow(name)
            
    def LoadVarianFile(self,filePath):
        Dir = os.path.dirname(filePath) #convert path to file to path of folder
        #Extract Procpar data if it exist------------------------
        #Initilize standard values
        freq = 300e6
        sw   = 50e3
        sw1  = 50e3
        if os.path.exists(Dir+os.path.sep+'procpar'):
            with open(Dir+os.path.sep+'procpar', 'r') as f: #read entire procfile (data[0] gives first line)
                data = f.read().split('\n')
            for s in range(0,len(data)): #exctract info from procpar
                if data[s].startswith('sfrq '):
                    freq=float(data[s+1].split()[1])*1e6 #convert to MHz
                elif data[s].startswith('sw '):
                    sw=float(data[s+1].split()[1])
                elif data[s].startswith('sw1 '):
                    sw1=float(data[s+1].split()[1])
        else:
            print(Dir+os.path.sep+'procpar does not exits, used standard sw and freq')
        #Get fid data----------------------------- 
        if os.path.exists(Dir+os.path.sep+'fid'):    
            with open(Dir+os.path.sep+'fid', "rb") as f:
                raw = np.fromfile(f, np.int32,6) #read 6 steps, 32 bits
                nblocks = unpack('>l', raw[0])[0] #unpack bitstring using bigendian and as LONG interger
                ntraces = unpack('>l', raw[1])[0]
                npoints = unpack('>l', raw[2])[0]
                ebytes = unpack('>l', raw[3])[0]
                tbytes = unpack('>l', raw[4])[0]
                bbytes = unpack('>l', raw[5])[0]
                raw = np.fromfile(f, np.int16,2) #16bit, 2 steps
                vers_id = unpack('>h', raw[0])[0] #bigendian short
                status = unpack('>h', raw[1])[0]
                raw = np.fromfile(f, np.int32,1) 
                nbheaders = unpack('>l', raw[0])[0]
                SizeTD2 = npoints
                SizeTD1 = nblocks*ntraces
                a = []
                fid32 = int(bin(status)[-3]) #check if 32 bits, or float
                fidfloat = int(bin(status)[-4])
                for iter1 in range(0,nblocks): #now read all blocks
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
        if SizeTD1 is 1: #convert to 1D dat if the data is 1D (so no 1xnp data, but np)
            fid = fid[0][:]
            masterData=sc.Spectrum(fid,lambda self :self.LoadVarianFile(filePath),[freq],[sw])
        else: #For 2D data
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
        Dir = os.path.dirname(filePath) #convert path to file to path of folder
        if os.path.exists(Dir+os.path.sep+'acqus'):
            with open(Dir+os.path.sep+'acqus', 'r') as f: 
                data = f.read().split('\n')
            for s in range(0,len(data)): #exctract info from acqus
                if data[s].startswith('##$TD='):
                    sizeTD2 = int(data[s][6:])
                if data[s].startswith('##$SFO1='):
                    freq2 = float(data[s][8:])*1e6
                if data[s].startswith('##$SW_h='):
                    SW2 = float(data[s][8:])
                if data[s].startswith('##$BYTORDA='):
                    ByteOrder = int(data[s][11:]) #1 little endian, 0 big endian 
        sizeTD1=1 #Preset to one           
        if os.path.exists(Dir+os.path.sep+'acqu2s'): #read 2d pars if available
            with open(Dir+os.path.sep+'acqu2s', 'r') as f: 
                data2 = f.read().split('\n')
            for s in range(0,len(data2)): #exctract info from acqus
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
        if ByteOrder: #Bigendian if ByteOrder 1, otherwise smallendian
            RawInt=raw.newbyteorder('b')
        else:
            RawInt=raw.newbyteorder('l')
        ComplexData = np.array(RawInt[0:len(RawInt):2])+1j*np.array(RawInt[1:len(RawInt):2])
        spec = [False]
        if sizeTD1 is 1:
            #data = np.transpose(ComplexData)[0][:] #convert to 1D np.array
            masterData=sc.Spectrum(ComplexData,lambda self :self.LoadBrukerTopspin(filePath),[freq2],[SW2],spec)
        else:
            data = ComplexData.reshape(sizeTD1,sizeTD2/2)
            masterData=sc.Spectrum(data,lambda self :self.LoadBrukerTopspin(filePath),[freq1,freq2],[SW1,SW2],spec*2)
        return masterData
                
    def LoadChemFile(self,filePath):
        Dir = os.path.dirname(filePath)
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
        data = np.array(fid) #convert to numpy array
        spec = [False]                    
        if sizeTD1 is 1:
            data = data[0][:] #convert to 1D np.array
            masterData=sc.Spectrum(data,lambda self :self.LoadChemFile(filePath),[freq*1e6],[sw],spec)
        else:
            #data = np.transpose(data.reshape((sizeTD1,sizeTD2)))
            data = data.reshape((sizeTD1,sizeTD2))
            masterData=sc.Spectrum(data,lambda self :self.LoadChemFile(filePath),[freq*1e6]*2,[sw1,sw],spec*2)
        return masterData

    def LoadMagritek(self,filePath):
        #Magritek load script based on some Matlab files by Ole Brauckman
        Dir = os.path.dirname(filePath)
        DirFiles = os.listdir(Dir)
        Files2D = [x for x in DirFiles if '.2d' in x]
        Files1D = [x for x in DirFiles if '.1d' in x]
        H = dict(line.strip().split('=') for line in open(Dir+os.path.sep+'acqu.par','r'))
        sw = float(H['bandwidth                 '])*1000 #in kHz
        sizeTD2 = int(H['nrPnts                    '])
        freq = float(H['b1Freq                    '])
        if len(Files2D)==1:
            File=Files2D[0]
            sizeTD1 = int(H['nrSteps                   '])
            if 'bandwidth2                ' in H:
                sw1 = float(H['bandwidth2                ']) #Is already in kHz
            else:
                sw1 = 50e3 #set to default 50 kHz
            with open(Dir+os.path.sep+File,'rb') as f:
                raw = np.fromfile(f, np.float32)
            Data = raw[-2*sizeTD2*sizeTD1::] #Get last 2*sizeTD2*sizeTD1 points
            ComplexData = Data[0:Data.shape[0]:2]+1j*Data[1:Data.shape[0]:2]
            ComplexData = ComplexData.reshape((sizeTD1,sizeTD2))
            masterData=sc.Spectrum(ComplexData,lambda self :self.LoadMagritek(filePath),[freq*1e6]*2,[sw,sw1],[False]*2)
        elif len(Files1D)!=0:
            File = 'data.1d'
            with open(Dir+os.path.sep+File,'rb') as f:
                raw = np.fromfile(f, np.float32)
            Data = raw[-2*sizeTD2::] #Get last 2*sizeTD points
            ComplexData = Data[0:Data.shape[0]:2]+1j*Data[1:Data.shape[0]:2]
            masterData=sc.Spectrum(ComplexData,lambda self :self.LoadMagritek(filePath),[freq*1e6],[sw],[False])
        return masterData
            
    def LoadSimpsonFile(self,filePath):
        #Loads Simpson data (Fid or Spectrum) to the ssNake data format
        with open(filePath, 'r') as f: #read file
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
        if 'Normal' in FORMAT: #If normal format (e.g. not binary)
            data = []
            for iii in range(DataStart+1,DataEnd): #exctract data
                temp = Lines[iii].split()
                data.append(float(temp[0])+1j*float(temp[1]))
        elif 'BINARY' in FORMAT: #needs to be im-plemented
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
                        except: #if end of file, append zeroes
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
                        #Simpson to float
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
        data = np.array(data) #convert to numpy array
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

    def saveFigure(self):
        if self.mainWindow is not None:
            self.mainWindow.removeFromView()
        num = self.workspaces.index(self.mainWindow)
        self.mainWindow = fit.MainPlotWindow(self.root,self,self.mainWindow)
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
        
    def saveSimpsonFile(self):
        self.mainWindow.get_mainWindow().SaveSimpsonFile()
        
    def saveJSONFile(self):
        self.mainWindow.get_mainWindow().saveJSONFile()
        
    def saveMatlabFile(self):
        self.mainWindow.get_mainWindow().saveMatlabFile()
        
class Main1DWindow(Frame):
    def __init__(self,parent,mainProgram,masterData,duplicateCurrent=None,name=''):
        Frame.__init__(self,parent)
        self.name = name
        self.fig = mainProgram.getFig()
        self.canvas = FigureCanvasTkAgg(self.fig, master=weakref.proxy(self))
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="nswe")
        self.undoList = [] 
        self.redoList = []
        self.currentMacro = None
        self.redoMacro = []
        self.parent = parent
        self.mainProgram = mainProgram
        self.masterData = masterData
        if duplicateCurrent is not None:
            self.current = duplicateCurrent.copyCurrent(self,self.fig,self.canvas,masterData)
        else:
            self.current = sc.Current1D(self,self.fig,self.canvas,masterData)
        self.menubar = self.mainProgram.menubar
        self.sideframe=SideFrame(weakref.proxy(self))
        self.sideframe.grid(row=0,column=2,sticky='n')
        Separator(weakref.proxy(self),orient=VERTICAL).grid(row=0,column=1,rowspan=4,sticky='ns')
        self.bottomframe=BottomFrame(weakref.proxy(self))
        self.bottomframe.grid(row=1,column=0,sticky='w') 
        Separator(weakref.proxy(self),orient=HORIZONTAL).grid(row=2,sticky='ew')
        self.textframe=TextFrame(weakref.proxy(self))
        self.textframe.grid(row=3,column=0,sticky='s')  
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
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
        self.canvas._tkcanvas = None
        self.destroy()
        self.current.kill()
        del self.current
        del self.masterData
        del self.canvas
        self.mainProgram.donateFig(self.fig)

    def rescue(self):
        self.current.kill()
        self.current = sc.Current1D(self,self.current.fig,self.current.canvas,self.masterData)
        
    def removeFromView(self):
        self.menubar.delete("Edit")
        self.menubar.delete("Tools")
        self.menubar.delete("Matrix")
        self.menubar.delete("Fourier")
        self.menubar.delete("Fitting")
        self.menubar.delete("Combine")
        self.menubar.delete("Plot")
        self.editmenu.destroy()
        del self.editmenu
        self.matrixMenu.destroy()
        del self.matrixMenu
        self.plotMenu.destroy()
        del self.plotMenu
        self.toolMenu.destroy()
        del self.toolMenu
        self.fftMenu.destroy()
        del self.fftMenu
        self.fittingMenu.destroy()
        del self.fittingMenu
        self.combineMenu.destroy()
        self.pack_forget()

    def addToView(self):
	#the edit drop down menu
        self.editmenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Edit", menu=self.editmenu)
        self.editmenu.add_command(label="Undo", command=self.undo)
        self.editmenu.add_command(label="Redo", command=self.redo)
        self.editmenu.add_command(label="Reload", command=self.reloadLast)

	#the tool drop down menu
        self.toolMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Tools",menu=self.toolMenu)
        self.toolMenu.add_command(label="Real", command=self.real)
        self.toolMenu.add_command(label="Imag", command=self.imag)
        self.toolMenu.add_command(label="Abs", command=self.abs) 
        self.toolMenu.add_command(label="Apodize", command=self.createApodWindow)
        self.toolMenu.add_command(label="Phasing", command=self.createPhaseWindow)
        self.toolMenu.add_command(label="Sizing", command=self.createSizeWindow) 
        self.toolMenu.add_command(label="Swap Echo", command=self.createSwapEchoWindow)
        self.toolMenu.add_command(label="Shift Data", command=self.createShiftDataWindow)
        self.toolMenu.add_command(label="Offset correction", command=self.createDCWindow)
        self.toolMenu.add_command(label="Baseline correction", command=self.createBaselineWindow)
        self.toolMenu.add_command(label="States", command=self.states)
        self.toolMenu.add_command(label="States-TPPI", command=self.statesTPPI)
        self.toolMenu.add_command(label="Correct Bruker digital filter", command=self.BrukerDigital)
        #self.toolMenu.add_command(label="LPSVD", command=self.LPSVD)

        #the matrix drop down menu
        self.matrixMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Matrix",menu=self.matrixMenu)
        self.matrixMenu.add_command(label="Integrate", command=self.createIntegrateWindow)
        self.matrixMenu.add_command(label="Max", command=self.createMaxWindow)
        self.matrixMenu.add_command(label="Min", command=self.createMinWindow)
        self.matrixMenu.add_command(label="Max position", command=self.createArgMaxWindow)
        self.matrixMenu.add_command(label="Min position", command=self.createArgMinWindow)
        self.matrixMenu.add_command(label="Extract part", command=self.createRegionWindow)
        self.matrixMenu.add_command(label="Flip L/R", command=self.flipLR)
        self.matrixMenu.add_command(label="Delete", command=self.createDeleteWindow)
        self.matrixMenu.add_command(label="Split", command=self.createSplitWindow)
        self.matrixMenu.add_command(label="Multiply", command=self.createMultiplyWindow)
        self.matrixMenu.add_command(label="Concatenate", command=self.createConcatenateWindow)
        self.matrixMenu.add_command(label="Shearing", command=self.createShearingWindow)
        
        #the fft drop down menu
        self.fftMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Fourier",menu=self.fftMenu)
        self.fftMenu.add_command(label="Fourier transform", command=self.fourier)
        self.fftMenu.add_command(label="Real Fourier transform", command=self.realFourier)
        self.fftMenu.add_command(label="Fftshift", command=self.fftshift)
        self.fftMenu.add_command(label="Inv fftshift", command=self.invFftshift)
        self.fftMenu.add_command(label="Hilbert transform", command=self.hilbert)

	#the fitting drop down menu
        self.fittingMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Fitting",menu=self.fittingMenu)
        self.fittingMenu.add_command(label="S/N", command=self.createSNWindow)
        self.fittingMenu.add_command(label="FWHM", command=self.createFWHMWindow)
        self.fittingMenu.add_command(label="Relaxation Curve", command=self.createRelaxWindow)
        self.fittingMenu.add_command(label="Diffusion Curve", command=self.createDiffusionWindow)
        self.fittingMenu.add_command(label="Peak Deconvolution", command=self.createPeakDeconvWindow)
        self.fittingMenu.add_command(label="CSA tensor", command=self.createTensorDeconvWindow)
        self.fittingMenu.add_command(label="First order quadrupole", command=self.createQuad1DeconvWindow)
        self.fittingMenu.add_command(label="Second order quadrupole static", command=self.createQuad2StaticDeconvWindow)
        self.fittingMenu.add_command(label="Second order quadrupole MAS", command=self.createQuad2MASDeconvWindow)
        self.fittingMenu.add_command(label="Czjzek static", command=self.createQuad2StaticCzjzekWindow)
        self.fittingMenu.add_command(label="Czjzek MAS", command=self.createQuad2MASCzjzekWindow)
        
        #the combine drop down menu
        self.combineMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Combine",menu=self.combineMenu)
        self.combineMenu.add_command(label="Insert from workspace", command=self.createInsertWindow)
        self.combineMenu.add_command(label="Add", command=self.createAddWindow)
        self.combineMenu.add_command(label="Subtract", command=self.createSubtractWindow)

	#the plot drop down menu
        self.plotMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Plot",menu=self.plotMenu)
        self.plotMenu.add_command(label="1D plot", command=self.plot1D)
        self.plotMenu.add_command(label="Scatter plot", command=self.plotScatter)
        self.plotMenu.add_command(label="Stack plot", command=self.plotStack)
        self.plotMenu.add_command(label="Array plot", command=self.plotArray)
        self.plotMenu.add_command(label="Contour plot", command=self.plotContour)
        self.plotMenu.add_command(label="Skewed plot", command=self.plotSkewed)
        self.plotMenu.add_command(label="Set reference", command=self.createRefWindow)
        self.plotMenu.add_command(label="User x-axis", command=self.createXaxWindow)
        self.pack(fill=BOTH,expand=1)
        self.menuCheck()
        
    def menuEnable(self):
        for i in range(12):
            self.menubar.entryconfig(i,state='normal')
        self.sideframe.frameEnable()
        self.bottomframe.frameEnable()
        self.textframe.frameEnable()
        self.menuCheck()

    def menuCheck(self):
        if (len(self.masterData.data.shape) < 2):
            self.matrixMenu.entryconfig('Concatenate',state='disabled')
            self.matrixMenu.entryconfig('Shearing',state='disabled')
            self.plotMenu.entryconfig('Stack plot',state='disabled')
            self.plotMenu.entryconfig('Array plot',state='disabled')
            self.plotMenu.entryconfig('Contour plot',state='disabled')
            self.plotMenu.entryconfig('Skewed plot',state='disabled')
        else:
            self.matrixMenu.entryconfig('Concatenate',state='normal')
            self.matrixMenu.entryconfig('Shearing',state='normal')
            self.plotMenu.entryconfig('Stack plot',state='normal')
            self.plotMenu.entryconfig('Array plot',state='normal')
            self.plotMenu.entryconfig('Contour plot',state='normal')
            self.plotMenu.entryconfig('Skewed plot',state='normal')
        if not self.undoList:
            self.editmenu.entryconfig('Undo',state='disabled')
        else:
            self.editmenu.entryconfig('Undo',state='normal')
        if not self.redoList:
            self.editmenu.entryconfig('Redo',state='disabled')
        else:
            self.editmenu.entryconfig('Redo',state='normal')
        
    def menuDisable(self):
        for i in range(12):
            self.menubar.entryconfig(i,state='disabled')
        self.sideframe.frameDisable()
        self.bottomframe.frameDisable()
        self.textframe.frameDisable()

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
            elif iter1[0] == 'lpsvd':
                self.undoList.append(self.masterData.LPSVD(*iter1[1]))
            elif iter1[0] == 'states':
                self.undoList.append(self.masterData.states(*iter1[1]))
            elif iter1[0] == 'statesTPPI':
                self.undoList.append(self.masterData.statesTPPI(*iter1[1]))
            elif iter1[0] == 'baselineCorrection':
                self.undoList.append(self.masterData.baselineCorrection(*iter1[1]))
            elif iter1[0] == 'integrate':
                self.undoList.append(self.masterData.matrixManip(*iter1[1],which=0))
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
        name=asksaveasfilename(filetypes=(('JSON','.json'),))
        if not name:
            return
        struct = {}
        struct['dataReal'] = np.real(self.masterData.data).tolist()
        struct['dataImag'] = np.imag(self.masterData.data).tolist()
        struct['freq'] = self.masterData.freq.tolist()
        struct['sw'] = list(self.masterData.sw)
        struct['spec'] = list(self.masterData.spec)
        struct['wholeEcho'] = list(1*self.masterData.wholeEcho)
        struct['ref'] = np.array(self.masterData.ref,dtype=np.float).tolist()
        tmpXax = []
        for i in self.masterData.xaxArray:
            tmpXax.append(i.tolist())
        struct['xaxArray'] = tmpXax
        with open(name, 'w') as outfile:
            json.dump(struct, outfile)
        
    def saveMatlabFile(self):
        name=asksaveasfilename(filetypes=(('MATLAB file','.mat'),))
        if not name:
            return
        struct = {}
        struct['dim'] = self.masterData.dim #because matlab cannot remember the number of dimensions
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
        #TO DO:
        #Make sure that stat of second dimension (SPE/FID) is saved. This is not supported in original
        #Simpson format, but can it be included?
        try:
            if self.masterData.dim   > 2:
                print('Saving to Simpson format only allowed for 1D and 2D data!')
            else:
                if sum(self.masterData.spec)/len(self.masterData.spec)==1: #If all are true
                    f=asksaveasfile(mode='w',defaultextension=".spe")
                elif sum(self.masterData.spec) == 0: #If all are false (i.e all FID)
                    f=asksaveasfile(mode='w',defaultextension=".fid")
                if 'f' in locals(): #If no 'f', there is a mixed fid/spe format, which simpson does not support 
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
                        #for Line in self.masterData.data:
                        #    for SubLine in Line:
                        #        f.write(str(SubLine.real) +' '+str(SubLine.imag)+'\n')
                    f.write('END')
                    f.close()
        except:
             print('An error occured while saving to Simpson format')

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

    def createPhaseWindow(self):
        PhaseWindow(self)
        
    def createApodWindow(self):
        ApodWindow(self)

    def createSizeWindow(self):
        SizeWindow(self)

    def createSwapEchoWindow(self):
        SwapEchoWindow(self)

    def createShiftDataWindow(self):
        ShiftDataWindow(self)

    def createDCWindow(self):
        DCWindow(self)
        
    def createBaselineWindow(self):
        BaselineWindow(self)
        
    def createRefWindow(self):
        RefWindow(self)

    def createIntegrateWindow(self):
        integrateWindow(self)
        
    def createMultiplyWindow(self):
        MultiplyWindow(self)
        
    def createMaxWindow(self):
        maxWindow(self)
        
    def createMinWindow(self):
        minWindow(self)

    def createArgMaxWindow(self):
        argmaxWindow(self)
        
    def createArgMinWindow(self):
        argminWindow(self)
        
    def createRegionWindow(self):
        extractRegionWindow(self)

    def flipLR(self):
        self.redoList = []
        self.undoList.append(self.current.flipLR())
        self.menuCheck()
        
    def createDeleteWindow(self):
        DeleteWindow(self)
        
    def createSplitWindow(self):
        SplitWindow(self)
        
    def createConcatenateWindow(self):
        ConcatenateWindow(self)
        
    def createInsertWindow(self):
        InsertWindow(self)
        
    def createAddWindow(self):
        AddWindow(self)
        
    def createSubtractWindow(self):
        SubtractWindow(self)
        
    def createShearingWindow(self):
        if self.masterData.dim > 1:
            ShearingWindow(self)
        else:
            print('Data has too little dimensions for shearing transform')

    def BrukerDigital(self):
        FilePath = askopenfilename()
        if FilePath is not '': #if not canceled
            Dir = os.path.dirname(FilePath) #convert path to file to path of folder
            if os.path.exists(Dir+os.path.sep+'acqus'):
                with open(Dir+os.path.sep+'acqus', 'r') as f: 
                    data = f.read().split('\n')
                FilterCorrection = -1.0
                for s in range(0,len(data)): #exctract info from acqus
                    if data[s].startswith('##$GRPDLY='):
                        FilterCorrection = float(data[s][10:])
                    if data[s].startswith('##$DECIM='):
                        DECIM = int(data[s][9:])
                    if data[s].startswith('##$DSPFVS='):
                        DSPFVS = int(data[s][10:])
                if FilterCorrection == -1.0: #If the FilterCorrection has not been found in the acqus (old bruker format)
                    if DSPFVS == 10 or DSPFVS == 11 or DSPFVS == 12:#get from table
                        CorrectionList = [{'2':44.7500,'3':33.5000,'4':66.6250,'6':59.0833
                            ,'8':68.5625,'12':60.3750,'16':69.5313,'24':61.0208,'32':70.0156
                            ,'48':61.3438,'64':70.2578,'96':61.5052,'128':70.3789,'192':61.5859
                            ,'256':70.4395,'384':61.6263,'512':70.4697,'768':61.6465,'1024':70.4849,'1536':61.6566,'2048':70.4924},
                            {'2':46.0000,'3':36.5000,'4':48.0000,'6':50.1667,'8':53.2500,'12':69.5000,
                                        '16':72.2500,'24':70.1667,'32':72.7500,'48':70.5000,'64':73.0000,'96':70.6667,
                                        '128':72.5000,'192':71.3333,'256':72.2500,'384':71.6667,'512':72.1250,'768':71.8333,
                                        '1024':72.0625,'1536':71.9167,'2048':72.0313},{'2':46.311,'3':36.530,'4':47.870,'6':50.229,'8':53.289,'12':69.551,'16':71.600,
                                        '24':70.184,'32':72.138,'48':70.528,'64':72.348,'96':70.700,'128':72.524}]
                        #Take correction from database. Based on matNMR routine (Jacco van Beek), which is itself based 
                        #on a text by W. M. Westler and F. Abildgaard.
                        FilterCorrection = CorrectionList[10-DSPFVS][str(DECIM)]
                    else:
                        print('DSPFVS value not recognized (Bruker hardware version not known)')
                if FilterCorrection != -1.0: #If changed
                    self.redoList = []
                    self.undoList.append(self.current.applyPhase(0, FilterCorrection*2*np.pi))
                    self.menuCheck()
                    
    def LPSVD(self):
        self.redoList = []
        self.undoList.append(self.current.applyLPSVD())
        self.menuCheck()
                    
    def createSNWindow(self):
        SNWindow(self)
        
    def createFWHMWindow(self):
        FWHMWindow(self)
        
    def createRelaxWindow(self):
        self.mainProgram.createFitWindow(fit.RelaxWindow(self.parent,self.mainProgram,self.mainProgram.mainWindow))
        
    def createDiffusionWindow(self):
        self.mainProgram.createFitWindow(fit.DiffusionWindow(self.parent,self.mainProgram,self.mainProgram.mainWindow))
        
    def createPeakDeconvWindow(self):
        self.mainProgram.createFitWindow(fit.PeakDeconvWindow(self.parent,self.mainProgram,self.mainProgram.mainWindow))
        
    def createTensorDeconvWindow(self):
        self.mainProgram.createFitWindow(fit.TensorDeconvWindow(self.parent,self.mainProgram,self.mainProgram.mainWindow))

    def createQuad1DeconvWindow(self):
        self.mainProgram.createFitWindow(fit.Quad1DeconvWindow(self.parent,self.mainProgram,self.mainProgram.mainWindow))
        
    def createQuad2StaticDeconvWindow(self):
        self.mainProgram.createFitWindow(fit.Quad2DeconvWindow(self.parent,self.mainProgram,self.mainProgram.mainWindow))
        
    def createQuad2MASDeconvWindow(self):
        self.mainProgram.createFitWindow(fit.Quad2DeconvWindow(self.parent,self.mainProgram,self.mainProgram.mainWindow,True))

    def createQuad2StaticCzjzekWindow(self):
        self.mainProgram.createFitWindow(fit.Quad2CzjzekWindow(self.parent,self.mainProgram,self.mainProgram.mainWindow))
        
    def createQuad2MASCzjzekWindow(self):
        self.mainProgram.createFitWindow(fit.Quad2CzjzekWindow(self.parent,self.mainProgram,self.mainProgram.mainWindow,True))
        
    def plot1D(self):
        tmpcurrent = sc.Current1D(self,self.fig,self.canvas,self.masterData)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()

    def plotScatter(self):
        tmpcurrent = sc.CurrentScatter(self,self.fig,self.canvas,self.masterData)
        self.current.kill()
        del self.current
        self.current = tmpcurrent
        self.updAllFrames()

    def plotStack(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentStacked(self,self.fig,self.canvas,self.masterData)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")

    def plotArray(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentArrayed(self,self.fig,self.canvas,self.masterData)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")
            
    def plotContour(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentContour(self,self.fig,self.canvas,self.masterData)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")
            
    def plotSkewed(self):
        if len(self.masterData.data.shape) > 1:
            tmpcurrent = sc.CurrentSkewed(self,self.fig,self.canvas,self.masterData)
            self.current.kill()
            del self.current
            self.current = tmpcurrent
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")
            
    def createXaxWindow(self):
        XaxWindow(self)

    def updAllFrames(self):
        self.bottomframe.upd()
        self.sideframe.upd()

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
#the sideframe class which displays (if necessary) the position of the shown data relative to the full data matrix
class SideFrame(Frame):
    def __init__(self, parent):
        Frame.__init__(self,parent)
        self.parent = parent
        self.labels=[]
        self.entries=[]
        self.entryVars=[]
        self.buttons1=[]
        self.button1Var=IntVar()
        self.button1Var.set(0)
        #setting needed for 2d style plots
        self.buttons2=[]
        self.button2Var=IntVar()
        self.button2Var.set(1)
        self.plotIs2D = False
        self.spacing = StringVar()
        self.skew = StringVar()
        self.elev = StringVar()
        self.numLevels = StringVar()
        self.maxLevels = StringVar()
        self.minLevels = StringVar()
        self.from2D = StringVar()
        self.to2D = StringVar()
        self.step2D = StringVar()
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0)
        self.frame2 = Frame(self)
        self.sep = Separator(self,orient=HORIZONTAL).grid(row=1, sticky='ew')
        self.frame2.grid(row=2,column=0,sticky='nwe')
        self.frame2.grid_columnconfigure(0,weight=1)
        self.upd()

    def frameEnable(self):
        for child in self.frame1.winfo_children():
            child.configure(state='normal')
        for child in self.frame2.winfo_children():
            child.configure(state='normal')
            
    def frameDisable(self):
        for child in self.frame1.winfo_children():
            child.configure(state='disabled')
        for child in self.frame2.winfo_children():
            child.configure(state='disabled')
            
    def upd(self): #destroy the old widgets and create new ones
        current = self.parent.current
        self.shape = current.data.data.shape
        self.length = len(self.shape)
        if self.length < 2:
            if self.frame1 is not None:
                self.frame1.destroy()
                self.frame1 = Frame(self)
                self.frame1.grid(row=0,column=0)
            if self.frame2 is not None:
                self.frame2.destroy()
                self.frame2 = Frame(self)
                self.frame2.grid(row=2,column=0,sticky='nwe')
                self.frame2.grid_columnconfigure(0,weight=1)
            if self.sep is not None:
                self.sep.destroy()
                self.sep = Separator(self,orient=HORIZONTAL)
                self.sep.grid(row=1, sticky='ew')
        else:
            for widget in self.frame1.winfo_children():
                widget.destroy()
            for widget in self.frame2.winfo_children():
                widget.destroy()
        self.button1Var.set(current.axes)
        offset = 0
        self.plotIs2D = isinstance(current, (sc.CurrentStacked,sc.CurrentArrayed,sc.CurrentContour,sc.CurrentSkewed))
        if self.plotIs2D:
            offset = 1
            self.button2Var.set(current.axes2)
        self.labels = []
        self.entries=[]
        self.buttons1=[]
        self.buttons2=[]
        self.entryVars = []
        if self.length > 1:
            for num in range(self.length):
                self.buttons1.append(Radiobutton(self.frame1, variable=self.button1Var, value=num, command=lambda: self.setAxes(True)))
                self.buttons1[num].grid(row=num*2+1,column=0)
                if self.plotIs2D:
                    self.buttons2.append(Radiobutton(self.frame1, variable=self.button2Var, value=num, command=lambda: self.setAxes(False)))
                    self.buttons2[num].grid(row=num*2+1,column=1)
                self.labels.append(Label(self.frame1,text="D"+str(num+1))) 
                self.labels[num].grid(row=num*2,column=1+offset)
                self.entryVars.append(StringVar())
                if not self.plotIs2D:
                    if num < current.axes:
                        self.entryVars[num].set(str(current.locList[num]))
                    elif num == current.axes:
                        self.entryVars[num].set("0")
                    else:
                        self.entryVars[num].set(str(current.locList[num-1]))
                else:
                    if (num < current.axes) and (num < current.axes2):
                        self.entryVars[num].set(str(current.locList[num]))
                    elif (num == current.axes) or (num == current.axes2):
                        self.entryVars[num].set("0")
                    elif (num > current.axes) and (num > current.axes2):
                        self.entryVars[num].set(str(current.locList[num-2]))
                    else:
                        self.entryVars[num].set(str(current.locList[num-1]))
                self.entries.append(Spinbox(self.frame1,textvariable=self.entryVars[num],from_=0,to=self.shape[num]-1,justify="center",width=10,command=lambda event=None,num=num: self.getSlice(event,num)))
                self.entries[num].bind("<Return>", lambda event=None,num=num: self.getSlice(event,num)) 
                self.entries[num].bind("<KP_Enter>", lambda event=None,num=num: self.getSlice(event,num)) 
                self.entries[num].grid(row=num*2+1,column=1+offset)
            if isinstance(current, (sc.CurrentStacked,sc.CurrentArrayed,sc.CurrentSkewed)):
                if current.stackBegin is not None:
                    self.from2D.set(str(current.stackBegin))
                else:
                    self.from2D.set('0')
                if current.stackEnd is not None:
                    self.to2D.set(str(current.stackEnd))
                else:
                    self.to2D.set(str(self.shape[current.axes2]))
                if current.stackStep is not None:
                    self.step2D.set(str(current.stackStep))
                else:
                    self.step2D.set('1')
                Label(self.frame2,text="From").grid(row=1,column=0,sticky='n')
                self.fromSpin = Spinbox(self.frame2,textvariable=self.from2D,from_=0,to=int(self.to2D.get())-1,justify="center",width=10,command=self.setToFrom)
                self.fromSpin.bind("<Return>", self.setToFrom) 
                self.fromSpin.bind("<KP_Enter>", self.setToFrom)
                self.fromSpin.grid(row=2,column=0)
                Label(self.frame2,text="To").grid(row=3,column=0,sticky='n')
                self.toSpin = Spinbox(self.frame2,textvariable=self.to2D,from_=int(self.from2D.get())+1,to=self.shape[current.axes2],justify="center",width=10,command=self.setToFrom)
                self.toSpin.bind("<Return>", self.setToFrom) 
                self.toSpin.bind("<KP_Enter>", self.setToFrom) 
                self.toSpin.grid(row=4,column=0)
                Label(self.frame2,text="Step").grid(row=5,column=0,sticky='n')
                self.stepSpin = Spinbox(self.frame2,textvariable=self.step2D,from_=1,to=self.shape[current.axes2],justify="center",width=10,command=self.setToFrom)
                self.stepSpin.bind("<Return>", self.setToFrom) 
                self.stepSpin.bind("<KP_Enter>", self.setToFrom)
                self.stepSpin.grid(row=6,column=0)
                if isinstance(current, (sc.CurrentStacked,sc.CurrentArrayed)):
                    self.spacing.set('%.3g' % current.spacing)
                    Label(self.frame2,text="Spacing").grid(row=7,column=0,sticky='n')
                    self.spacingEntry = Entry(self.frame2,textvariable=self.spacing,justify="center",width=10)
                    self.spacingEntry.bind("<Return>", self.setSpacing) 
                    self.spacingEntry.bind("<KP_Enter>", self.setSpacing) 
                    self.spacingEntry.grid(row=8,column=0)

                if isinstance(current, (sc.CurrentSkewed)):
                    self.skew.set('%.2f' % current.skewed)
                    Label(self.frame2,text="Skew").grid(row=7,column=0,sticky='n')
                    self.skewEntry = Entry(self.frame2,textvariable=self.skew,justify="center",width=10)
                    self.skewEntry.bind("<Return>", self.setSkew) 
                    self.skewEntry.bind("<KP_Enter>", self.setSkew) 
                    self.skewEntry.grid(row=8,column=0)
                    self.elev.set('%.1f' % current.elevation)
                    Label(self.frame2,text="Elevation").grid(row=9,column=0,sticky='n')
                    self.elevEntry = Entry(self.frame2,textvariable=self.elev,justify="center",width=10)
                    self.elevEntry.bind("<Return>", self.setSkew) 
                    self.elevEntry.bind("<KP_Enter>", self.setSkew) 
                    self.elevEntry.grid(row=10,column=0)
            if isinstance(current, (sc.CurrentContour)):
                self.numLevels.set(str(current.numLevels))
                self.maxLevels.set(str(current.maxLevels*100.0))
                self.minLevels.set(str(current.minLevels*100.0))
                Label(self.frame2,text="Number of contours").grid(row=1,column=0,sticky='n')
                self.numLEntry = Entry(self.frame2,textvariable=self.numLevels,justify="center",width=10)
                self.numLEntry.bind("<Return>", self.setContour) 
                self.numLEntry.bind("<KP_Enter>", self.setContour) 
                self.numLEntry.grid(row=2,column=0)
                Label(self.frame2,text="Highest contour [%]").grid(row=3,column=0,sticky='n')
                self.maxLEntry = Entry(self.frame2,textvariable=self.maxLevels,justify="center",width=10)
                self.maxLEntry.bind("<Return>", self.setContour) 
                self.maxLEntry.bind("<KP_Enter>", self.setContour) 
                self.maxLEntry.grid(row=4,column=0)
                Label(self.frame2,text="Lowest contour [%]").grid(row=5,column=0,sticky='n')
                self.minLEntry = Entry(self.frame2,textvariable=self.minLevels,justify="center",width=10)
                self.minLEntry.bind("<Return>", self.setContour) 
                self.minLEntry.bind("<KP_Enter>", self.setContour) 
                self.minLEntry.grid(row=6,column=0)
                    
    def setToFrom(self, *args):
        current = self.parent.current
        if isinstance(current, (sc.CurrentStacked,sc.CurrentArrayed,sc.CurrentSkewed)):
            fromVar = int(safeEval(self.from2D.get()))
            toVar = int(safeEval(self.to2D.get()))
            stepVar = int(safeEval(self.step2D.get()))
            if fromVar < 0:
                fromVar = 0
            elif fromVar > self.shape[current.axes2]-1:
                fromVar = self.shape[current.axes2]-1
            if toVar > self.shape[current.axes2]:
                toVar = self.shape[current.axes2]
            elif toVar <= fromVar:
                toVar = fromVar + 1
            if stepVar < 1:
                stepVar = 1
            elif stepVar > self.shape[current.axes2]:
                stepVar = self.shape[current.axes2]
            current.stackSelect(fromVar,toVar,stepVar)
            self.from2D.set(str(fromVar))
            self.to2D.set(str(toVar))
            self.step2D.set(str(stepVar))
            self.fromSpin.config(to=toVar-1)
            self.toSpin.config(from_=fromVar+1)

    def setSpacing(self, *args):
        var =float(safeEval(self.spacing.get()))
        self.spacing.set('%.3g' % var)
        self.parent.current.setSpacing(var)

    def setSkew(self, *args):
        var =float(safeEval(self.skew.get()))
        self.skew.set('%.2f' % var)
        var2 =float(safeEval(self.elev.get()))
        self.elev.set('%.1f' % var2)
        self.parent.current.setSkewed(var,var2)
        
    def setContour(self, *args):
        var1 = int(round(safeEval(self.numLevels.get())))
        self.numLevels.set(str(var1))
        var2 =float(safeEval(self.maxLevels.get()))
        self.maxLevels.set('%.1f' % var2)
        var3 =float(safeEval(self.minLevels.get()))
        self.minLevels.set('%.1f' % var3)
        self.parent.current.setLevels(var1,var2/100.0,var3/100.0)
        
    def setAxes(self,first=True):
        if self.plotIs2D:
            axes= self.button1Var.get()
            axes2=self.button2Var.get()
            if axes==axes2:
                if first:
                    axes2 = self.parent.current.axes
                else:
                    axes = self.parent.current.axes2
            self.button2Var.set(axes2)
            self.getSlice(None, axes,True)
        else:
            self.getSlice(None, self.button1Var.get(),True)

    def getSlice(self, event, entryNum, button=False): #change the slice which is currently displayed
        if button:
            dimNum = entryNum
        elif not self.plotIs2D:
            if entryNum == self.parent.current.axes:
                if entryNum == self.length-1:
                    dimNum = self.length-2
                else:
                    dimNum = self.length-1
            else:
                dimNum = self.parent.current.axes
        else:
            dimNum = self.parent.current.axes

        locList=[]
        for num in range(self.length):
            appendLoc = True
            if self.plotIs2D and (num == self.button2Var.get()):
                appendLoc = False
            inp = safeEval(self.entryVars[num].get())
            if num == dimNum:
                pass
            else:
                if inp < -self.shape[num]:
                    val=int(round(-(self.shape[num])))
                elif inp >= self.shape[num]:
                    val=int(round(self.shape[num]-1))
                elif inp < 0:
                    val = int(round(self.shape[num] + inp)) 
                else:
                    val = int(round(inp))
                if appendLoc:
                    locList.append(val)
                self.entryVars[num].set(val)
        self.button1Var.set(dimNum)
        if self.plotIs2D:
            self.parent.current.setBlock(dimNum,self.button2Var.get(), locList)
        else:
            self.parent.current.setSlice(dimNum,locList)
        self.parent.bottomframe.upd()
        
################################################################################  
#the bottom frame holding the fourier button and stuff      
class BottomFrame(Frame):
    def __init__(self, parent):
        Frame.__init__(self,parent)
        self.parent = parent
        self.specVal = IntVar() #value for the time/freq radiobutton
        self.freqVal = StringVar() #value for frequency entybox
        self.swVal = StringVar() #value for sw entrybox
        self.plotOption = StringVar() #value for dropdown plot type box
        self.plotOption.set("Real")
        self.axisOption1 = StringVar()
        self.axisOption2 = StringVar()
        self.axis2Option1 = StringVar()
        self.axis2Option2 = StringVar()
        self.echoTick = IntVar()
        self.echoTick.set(0)
        Button(self, text="Fourier",command=self.parent.fourier).grid(row=0,column=0,rowspan=2)
        self.rb1 = Radiobutton(self,text="Time",variable=self.specVal,value=0,command=self.changeSpec)
        self.rb1.grid(row=0,column=1)
        self.rb2 = Radiobutton(self,text="Frequency",variable=self.specVal,value=1,command=self.changeSpec)
        self.rb2.grid(row=1,column=1)
        Checkbutton(self,text="Whole echo",variable=self.echoTick, command=self.setWholeEcho).grid(row=0,column=2,rowspan=2)
        Label(self,text="Freq [MHz]").grid(row=0,column=3)
        self.freqEntry = Entry(self,textvariable=self.freqVal,justify='center',width=10)
        self.freqEntry.bind("<Return>", self.changeFreq)
        self.freqEntry.bind("<KP_Enter>", self.changeFreq)
        self.freqEntry.grid(row=1,column=3)
        Label(self,text="Sweepwidth [kHz]").grid(row=0,column=4)
        self.swEntry = Entry(self,textvariable=self.swVal,justify='center',width=10)
        self.swEntry.bind("<Return>", self.changeFreq)
        self.swEntry.bind("<KP_Enter>", self.changeFreq)
        self.swEntry.grid(row=1,column=4)
        Label(self,text="Plot").grid(row=0,column=5)
        self.plotDrop = OptionMenu(self, self.plotOption,self.plotOption.get(),"Real", "Imag", "Both","Abs",command=self.changePlot)
        self.plotDrop.grid(row=1,column=5)
        Label(self,text="Axis").grid(row=0,column=6)
        self.axisDropTime = OptionMenu(self, self.axisOption1, self.axisOption1.get(), "s", "ms", u"\u03bcs",command=self.changeAxis)
        self.axisDropFreq = OptionMenu(self, self.axisOption2, self.axisOption2.get(), "Hz", "kHz", "MHz","ppm",command=self.changeAxis)
        self.axisDropTime.grid(row=1,column=6)
        self.axisDropFreq.grid(row=1,column=6)
        self.axisDropTime2 = OptionMenu(self, self.axis2Option1, self.axis2Option1.get(), "s", "ms", u"\u03bcs",command=self.changeAxis)
        self.axisDropFreq2 = OptionMenu(self, self.axis2Option2, self.axis2Option2.get(), "Hz", "kHz", "MHz","ppm",command=self.changeAxis)
        self.axisDropTime2.grid(row=2,column=6)
        self.axisDropFreq2.grid(row=2,column=6)
        self.swEntry
        self.upd()
        
    def frameEnable(self):
        for child in self.winfo_children():
            child.configure(state='normal')
            
    def frameDisable(self):
        for child in self.winfo_children():
            child.configure(state='disabled')
        
    def upd(self): #upd the values displayed in the bottom menu
        self.freqVal.set('%.6f' %(self.parent.current.freq/1000000.0)) #show in MHz
        self.swVal.set('%.6f' %(self.parent.current.sw/1000.0)) #show in kHz
        self.axisDropTime2.grid_forget()
        self.axisDropFreq2.grid_forget()
        if self.parent.current.spec==0:
            self.specVal.set(0)
            self.axisDropFreq.grid_forget()
            self.axisDropTime.grid(row=1,column=6)
            val = self.parent.current.axType
            if val == 0:
                self.axisOption1.set("s")
            elif val == 1:
                self.axisOption1.set("ms")
            elif val == 2:
                self.axisOption1.set( u"\u03bcs")
                
        elif self.parent.current.spec==1:
            self.specVal.set(1)
            self.axisDropTime.grid_forget()
            self.axisDropFreq.grid(row=1,column=6)
            val = self.parent.current.axType
            if val == 0:
                self.axisOption2.set("Hz")
            elif val == 1:
                self.axisOption2.set("kHz")
            elif val == 2:
                self.axisOption2.set("MHz")
            elif val == 3:
                self.axisOption2.set("ppm")
        if isinstance(self.parent.current,sc.CurrentContour):
            if self.parent.current.spec2==0:
                self.axisDropTime2.grid(row=2,column=6)
                val = self.parent.current.axType2
                if val == 0:
                    self.axis2Option1.set("s")
                elif val == 1:
                    self.axis2Option1.set("ms")
                elif val == 2:
                    self.axis2Option1.set( u"\u03bcs")
            elif self.parent.current.spec2==1:
                self.axisDropFreq2.grid(row=2,column=6)
                val = self.parent.current.axType2
                if val == 0:
                    self.axis2Option2.set("Hz")
                elif val == 1:
                    self.axis2Option2.set("kHz")
                elif val == 2:
                    self.axis2Option2.set("MHz")
                elif val == 3:
                    self.axis2Option2.set("ppm")
        if self.parent.current.wholeEcho:
            self.echoTick.set(1)
        else:
            self.echoTick.set(0)

    def setWholeEcho(self):
        self.parent.undoList.append(self.parent.current.setWholeEcho(self.echoTick.get()))
        self.parent.menuCheck()

    def changeSpec(self, *args): #change from time to spectral domain and vice versa
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.changeSpec(self.specVal.get()))
        self.upd()
        self.parent.menuCheck()

    def changeFreq(self, *args): #change the frequency and sw of the displayed axes
        freq = safeEval(self.freqVal.get())*1000000 #show in MHz
        sw = safeEval(self.swVal.get())*1000 #show in kHz
        if freq != 0 and sw != 0:
            self.parent.setFreq(freq,sw)
        self.upd()
    
    def changePlot(self, *args): #change the plot type
        pType = self.plotOption.get()
        if pType == "Real":
            self.parent.current.plotType=0
        elif pType == "Imag":
            self.parent.current.plotType=1
        elif pType == "Both":
            self.parent.current.plotType=2
        elif pType == "Abs":
            self.parent.current.plotType=3
        self.parent.current.showFid()

    def changeAxis(self, *args):
        if self.parent.current.spec == 0:
            pType = self.axisOption1.get()
            if pType == "s":
                self.parent.current.setAxType(0)
            elif pType == "ms":
                self.parent.current.setAxType(1)
            elif pType == u"\u03bcs":
                self.parent.current.setAxType(2)
        if self.parent.current.spec == 1:
            pType = self.axisOption2.get()
            if pType == "Hz":
                self.parent.current.setAxType(0)
            elif pType == "kHz":
                self.parent.current.setAxType(1)
            elif pType == "MHz":
                self.parent.current.setAxType(2)
            elif pType == "ppm":
                self.parent.current.setAxType('ppm')
        if isinstance(self.parent.current,sc.CurrentContour):
            if self.parent.current.spec2 == 0:
                pType = self.axis2Option1.get()
                if pType == "s":
                    self.parent.current.setAxType2(0)
                elif pType == "ms":
                    self.parent.current.setAxType2(1)
                elif pType == u"\u03bcs":
                    self.parent.current.setAxType2(2)
            if self.parent.current.spec2 == 1:
                pType = self.axis2Option2.get()
                if pType == "Hz":
                    self.parent.current.setAxType2(0)
                elif pType == "kHz":
                    self.parent.current.setAxType2(1)
                elif pType == "MHz":
                    self.parent.current.setAxType2(2)
                elif pType == "ppm":
                    self.parent.current.setAxType2('ppm')
        self.parent.current.showFid()

##################################################################
#the frame showing the get position data
class TextFrame(Frame):
    def __init__(self, parent):
        Frame.__init__(self,parent)
        self.parent = parent
        self.pos = StringVar()      #number of get_position data point
        self.pos.set(str(0))
        self.oldx = 0.0
        self.oldy = 0.0
        self.xpoint = StringVar()   #x value of the get_position data point
        self.xpoint.set(str(self.oldx))
        self.ypoint = StringVar()   #y value of the get_position data point
        self.ypoint.set(str(self.oldy))
        self.deltaxpoint = StringVar()   #x value of the get_position data point
        self.deltaxpoint.set(str(0.0))
        self.deltaypoint = StringVar()   #y value of the get_position data point
        self.deltaypoint.set(str(0.0))
        Button(self,text="get Position", command=self.getPosition).grid(row=0,column=0)
        Label(self,text="Position:").grid(row=0,column=1)
        Entry(self,textvariable=self.pos,justify='center',width=10).grid(row=0,column=2)
        Label(self,text="x-value:").grid(row=0,column=3)
        Entry(self,textvariable=self.xpoint,justify='center',width=10).grid(row=0,column=4)
        Label(self,text="y-value:").grid(row=0,column=5)
        Entry(self,textvariable=self.ypoint,justify='center',width=10).grid(row=0,column=6)
        Label(self,text=u"\u0394x:").grid(row=0,column=7)
        Entry(self,textvariable=self.deltaxpoint,justify='center',width=10).grid(row=0,column=8)
        Label(self,text=u"\u0394y:").grid(row=0,column=9)
        Entry(self,textvariable=self.deltaypoint,justify='center',width=10).grid(row=0,column=10)

    def frameEnable(self):
        for child in self.winfo_children():
            child.configure(state='normal')
            
    def frameDisable(self):
        for child in self.winfo_children():
            child.configure(state='disabled')
        
    def setLabels(self,position):
        self.deltaxpoint.set('%.3g' % np.abs(self.oldx-position[1]))
        self.deltaypoint.set('%.3g' % np.abs(self.oldy-position[2]))
        self.pos.set(str(position[0]))
        self.xpoint.set('%.3g' % position[1])
        self.ypoint.set('%.3g' % position[2])
        self.oldx = position[1]
        self.oldy = position[2]

    def getPosition(self, *args):
        self.parent.current.peakPickFunc = lambda pos,self=self: self.setLabels(pos) 
        self.parent.current.peakPick = True
        
#################################################################################   
class PhaseWindow(Toplevel): #a window for phasing the data
    def __init__(self, parent):
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Phasing")
        self.resizable(width=FALSE, height=FALSE)
        self.parent.menuDisable()
        #initialize variables for the widgets
        self.P1LIMIT=540.0 
        self.zeroVal = 0.0
        self.firstVal = 0.0
        self.refVal = 0.0
        self.zeroValue = StringVar()
        self.zeroValue.set("0.0")
        self.firstValue = StringVar()
        self.firstValue.set("0.0")
        self.refValue = StringVar()
        self.refValue.set("0.0")
        self.singleSlice = IntVar()
        self.singleSlice.set(0)
        #set stepsizes for the buttons
        self.phase0step = 1.0
        self.phase1step = 1.0
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0)
        Label(self.frame1,text="Zero order phasing").grid(row=0,column=0,columnspan=3)
        Button(self.frame1,text="Autophase 0th order",command=lambda: self.autophase(0)).grid(row=1,column=1)
        self.zeroEntry = Entry(self.frame1,textvariable=self.zeroValue,justify="center")
        self.zeroEntry.bind("<Return>", self.inputZeroOrder)
        self.zeroEntry.bind("<KP_Enter>", self.inputZeroOrder)
        self.zeroEntry.grid(row=2,column=1)
        tk.Button(self.frame1,text="<",repeatdelay=100, repeatinterval=1,command=lambda:self.stepPhase(-1,0)).grid(row=2,column=0)
        tk.Button(self.frame1,text=">",repeatdelay=100, repeatinterval=1,command=lambda:self.stepPhase(1,0)).grid(row=2,column=2)
        self.zeroScale=Scale(self.frame1, from_=-180, to=180,  orient="horizontal", command=self.setZeroOrder,length=300)
        self.zeroScale.grid(row=3,column=0,columnspan=3)
        Label(self.frame1,text="First order phasing").grid(row=4,column=0,columnspan=3)
        Button(self.frame1,text="Autophase 0th+1st order",command=lambda: self.autophase(1)).grid(row=5,column=1)
        self.firstEntry = Entry(self.frame1,textvariable=self.firstValue,justify="center")
        self.firstEntry.bind("<Return>", self.inputFirstOrder) 
        self.firstEntry.bind("<KP_Enter>", self.inputFirstOrder) 
        self.firstEntry.grid(row=6,column=1)
        tk.Button(self.frame1,text="<",repeatdelay=100, repeatinterval=1,command=lambda:self.stepPhase(0,-1)).grid(row=6,column=0)
        tk.Button(self.frame1,text=">",repeatdelay=100, repeatinterval=1,command=lambda:self.stepPhase(0,1)).grid(row=6,column=2)
        self.firstScale=Scale(self.frame1, from_=-self.P1LIMIT, to=self.P1LIMIT, orient="horizontal", command=self.setFirstOrder,length=300)
        self.firstScale.grid(row=7,column=0,columnspan=3)
        if self.parent.current.spec > 0:
            Label(self.frame1,text="Reference").grid(row=8,column=0,columnspan=3)
            self.refEntry = Entry(self.frame1,textvariable=self.refValue,justify="center")
            self.refEntry.bind("<Return>", self.inputRef) 
            self.refEntry.bind("<KP_Enter>", self.inputRef)
            self.refEntry.grid(row=9,column=1)
        Button(self.frame1, text="Pick reference", command=self.pickRef).grid(row=10,column=1)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1,column=0)
        Checkbutton(self.frame2, text="Single slice", variable=self.singleSlice).grid(row=0,column=0,columnspan=2)
        Button(self.frame2, text="Apply",command=self.applyPhaseAndClose).grid(row=1,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=1,column=1)      
        
    def setZeroOrder(self,value, *args): #function called by the zero order scale widget
        self.zeroVal = float(value)
        self.zeroValue.set('%.2f' % self.zeroVal)
        self.parent.current.setPhaseInter(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0)
        
    def inputZeroOrder(self, *args): #function called by the zero order entry widget
        inp = safeEval(self.zeroValue.get())
        self.zeroVal = np.mod(inp+180,360)-180
        self.zeroScale.set(self.zeroVal) #setting the scale to a value calls the previous function, so the phase of current doesn't need to be set here

    def setFirstOrder(self,value, *args): #function called by the first order scale widget
        newZero = (self.zeroVal-(float(value)-self.firstVal)*self.refVal/self.parent.current.sw) #calculate the new zero order phase depending on the reference
        self.zeroVal = np.mod(newZero+180,360)-180
        self.firstVal = float(value)
        self.firstValue.set('%.2f' % self.firstVal)
        self.zeroValue.set('%.2f' % self.zeroVal)
        self.zeroScale.set(self.zeroVal)

    def inputFirstOrder(self, *args): #function called by the first order entry widget
        tmpFirstVal = float(safeEval(self.firstValue.get()))
        if -self.P1LIMIT <= tmpFirstVal <= self.P1LIMIT:
            self.firstScale.set(tmpFirstVal) #setting the scale to a value calls the previous function, so the phase of current doesn't need to be set here
        elif -self.P1LIMIT > tmpFirstVal:
            if float(self.firstScale.get()) != -self.P1LIMIT:
                self.firstScale.set(-self.P1LIMIT)
            self.setFirstOrder(tmpFirstVal)
            self.parent.current.setPhaseInter(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0)
        elif self.P1LIMIT < tmpFirstVal:
            if float(self.firstScale.get()) != self.P1LIMIT:
                self.firstScale.set(self.P1LIMIT)
            self.setFirstOrder(tmpFirstVal)
            self.parent.current.setPhaseInter(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0)

    def autophase(self, num): #run the autophase for either zero order (0) or both orders (1)
        phases = self.parent.current.autoPhase(num)
        val = phases[0]/np.pi*180.0
        self.zeroVal=(np.mod(val+180,360)-180)
        self.zeroValue.set('%.2f' % self.zeroVal)
        self.zeroScale.set(self.zeroVal)
        if num == 1:
            val = phases[1]/np.pi*180.0
            self.firstVal = val
            self.firstValue.set('%.2f' % self.firstVal)
            self.inputFirstOrder()

    def stepPhase(self,phase0,phase1): #step phase from < and > keys
        inp = safeEval(self.zeroValue.get())+phase0*self.phase0step
        self.zeroVal = np.mod(inp+180,360)-180
        self.zeroScale.set(self.zeroVal)
        self.firstVal = safeEval(self.firstValue.get())+phase1*self.phase1step
        self.firstValue.set('%.2f' % self.firstVal)
        self.inputFirstOrder()

    def inputRef(self, *args): #set the reference from the entry widget
        self.refVal = safeEval(self.refValue.get())
        self.refValue.set('%.2f' % self.refVal)

    def setRef(self,value,*args):
        self.refVal = float(value)
        self.refValue.set('%.2f' % self.refVal)

    def pickRef(self, *args): #run the pick function to pick the reference value
        self.parent.current.peakPickFunc = lambda pos,self=self: self.setRef(self.parent.current.xax[pos[0]])
        self.parent.current.peakPick = True

    def cancelAndClose(self):
        self.parent.current.upd()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

    def applyPhaseAndClose(self):
        self.inputZeroOrder()
        self.inputFirstOrder()
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.applyPhase(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0,(self.singleSlice.get()==1)))
        self.parent.menuEnable()
        self.destroy()

################################################################
class ApodWindow(Toplevel): #a window for apodization
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Apodize")
        self.resizable(width=FALSE, height=FALSE)
        #initialize variables for the widgets
        self.lorTick = IntVar()
        self.lorVal = StringVar()
        self.lorVal.set("0.0")
        self.gaussTick = IntVar()
        self.gaussVal = StringVar()
        self.gaussVal.set("0.0")
        self.cos2Tick = IntVar()
        self.cos2Val = StringVar()
        self.cos2Val.set("1.0")
        self.hammingTick = IntVar()
        self.hammingVal = StringVar()
        self.hammingVal.set("1.0")
        self.shiftVal = StringVar()
        self.shiftVal.set("0.0")
        self.shiftingVal = StringVar()
        self.shiftingVal.set("0.0")
        if self.parent.current.data.dim > 1:
            options = list(map(str,np.delete(range(1,self.parent.current.data.dim+1),self.parent.current.axes)))
            self.shiftingAxes = StringVar()
            self.shiftingAxes.set(options[0])
        self.singleSlice = IntVar()
        self.singleSlice.set(0)
        #set stepsizes for the buttons
        self.lorstep = 1.0
        self.gaussstep = 1.0
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0)
        Label(self.frame1,text="Lorentzian").grid(row=0,column=0,columnspan=4)
        Checkbutton(self.frame1,variable=self.lorTick, command=lambda: self.checkEval(self.lorTick,self.lorEntry)).grid(row=1,column=1)
        self.lorEntry = Entry(self.frame1,textvariable=self.lorVal,justify="center", state='disabled')   
        tk.Button(self.frame1,text="<",repeatdelay=100, repeatinterval=1,command=lambda:self.stepLB(-0.5*self.parent.current.sw/(self.parent.current.data1D.shape[-1]),0)).grid(row=1,column=0)
        tk.Button(self.frame1,text=">",repeatdelay=100, repeatinterval=1,command=lambda:self.stepLB(0.5*self.parent.current.sw/(self.parent.current.data1D.shape[-1]),0)).grid(row=1,column=4)
        self.lorEntry.bind("<Return>", self.apodPreview)
        self.lorEntry.bind("<KP_Enter>", self.apodPreview)
        self.lorEntry.grid(row=1,column=2)
        self.lorScale=Scale(self.frame1, from_=0, to=100.0*self.parent.current.sw/(self.parent.current.data1D.shape[-1]),  orient="horizontal", command=self.setLor,length=200)
        self.lorScale.grid(row=2,column=1,columnspan=2)
        Label(self.frame1,text="Gaussian").grid(row=3,column=0,columnspan=4)
        Checkbutton(self.frame1,variable=self.gaussTick, command=lambda: self.checkEval(self.gaussTick,self.gaussEntry)).grid(row=4,column=1)
        self.gaussEntry = Entry(self.frame1,textvariable=self.gaussVal,justify="center", state='disabled')
        self.gaussEntry.bind("<Return>", self.apodPreview)
        self.gaussEntry.bind("<KP_Enter>", self.apodPreview)
        self.gaussEntry.grid(row=4,column=2)
        tk.Button(self.frame1,text="<",repeatdelay=100, repeatinterval=1,command=lambda:self.stepLB(0,-0.5*self.parent.current.sw/(self.parent.current.data1D.shape[-1]))).grid(row=4,column=0)
        tk.Button(self.frame1,text=">",repeatdelay=100, repeatinterval=1,command=lambda:self.stepLB(0,0.5*self.parent.current.sw/(self.parent.current.data1D.shape[-1]))).grid(row=4,column=4)
        self.gaussScale=Scale(self.frame1, from_=0, to=100.0*self.parent.current.sw/(self.parent.current.data1D.shape[-1]),  orient="horizontal", command=self.setGauss,length=200)
        self.gaussScale.grid(row=5,column=1,columnspan=2)
        Label(self.frame1,text="Cos^2").grid(row=6,column=0,columnspan=4)
        Checkbutton(self.frame1,variable=self.cos2Tick, command=lambda: self.checkEval(self.cos2Tick,self.cos2Entry)).grid(row=7,column=1)
        self.cos2Entry = Entry(self.frame1,textvariable=self.cos2Val,justify="center", state='disabled')
        self.cos2Entry.bind("<Return>", self.apodPreview)
        self.cos2Entry.bind("<KP_Enter>", self.apodPreview)
        self.cos2Entry.grid(row=7,column=2)
        Label(self.frame1,text="Hamming").grid(row=8,column=0,columnspan=4)
        Checkbutton(self.frame1,variable=self.hammingTick, command=lambda: self.checkEval(self.hammingTick,self.hammingEntry)).grid(row=9,column=1)
        self.hammingEntry = Entry(self.frame1,textvariable=self.hammingVal,justify="center", state='disabled')
        self.hammingEntry.bind("<Return>", self.apodPreview)
        self.hammingEntry.bind("<KP_Enter>", self.apodPreview)
        self.hammingEntry.grid(row=9,column=2)
        Label(self.frame1,text="Shift").grid(row=10,column=0,columnspan=4)
        self.shiftEntry = Entry(self.frame1,textvariable=self.shiftVal,justify="center")
        self.shiftEntry.bind("<Return>", self.apodPreview)
        self.shiftEntry.bind("<KP_Enter>", self.apodPreview)
        self.shiftEntry.grid(row=11,column=2)
        if self.parent.current.data.dim > 1:
            Label(self.frame1,text="Shifting").grid(row=12,column=0,columnspan=4)
            self.shiftingEntry = Entry(self.frame1,textvariable=self.shiftingVal,justify="center")
            self.shiftingEntry.bind("<Return>", self.apodPreview)
            self.shiftingEntry.bind("<KP_Enter>", self.apodPreview)
            self.shiftingEntry.grid(row=13,column=2)
            OptionMenu(self.frame1,self.shiftingAxes, self.shiftingAxes.get(),*options).grid(row=14,column=2)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1,column=0)
        Checkbutton(self.frame2, text="Single slice", variable=self.singleSlice).grid(row=0,column=0,columnspan=2)
        Button(self.frame2, text="Apply",command=self.applyApodAndClose).grid(row=1,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=1,column=1) 

    def checkEval(self,checkVar,entryVar): #change the state of the entry widget that the 
        if checkVar.get() == 0:
            entryVar.configure(state='disabled')
        else:
            entryVar.configure(state='normal')
        self.apodPreview()

    def setLor(self,value, *args): #set the lorentzian value from the scale widget
        if self.lorTick.get() == 0:
            self.lorTick.set(1)
            self.lorEntry.configure(state='normal')
        self.lorVal.set('%.2f' % float(value))
        self.apodPreview()

    def setGauss(self,value, *args): #set the gaussian value from the scale widget
        if self.gaussTick.get() == 0:
            self.gaussTick.set(1)
            self.gaussEntry.configure(state='normal')
        self.gaussVal.set('%.2f' % float(value))
        self.apodPreview()

    def apodPreview(self, *args): #display the apodization preview
        lor = None
        gauss = None
        cos2 = None
        hamming = None
        shifting = None
        shiftingAxes = 0
        if self.lorTick.get() == 1:
            lor = safeEval(self.lorVal.get())
            self.lorVal.set(lor)
        if self.gaussTick.get() == 1:
            gauss = safeEval(self.gaussVal.get())
            self.gaussVal.set(gauss)
        if self.cos2Tick.get() == 1:
            cos2 = safeEval(self.cos2Val.get())
            self.cos2Val.set(cos2)
        if self.hammingTick.get() == 1:
            hamming = safeEval(self.hammingVal.get())
            self.hammingVal.set(hamming)
        shift = safeEval(self.shiftVal.get())
        self.shiftVal.set(shift)
        if self.parent.current.data.dim > 1:
            shifting = safeEval(self.shiftingVal.get())
            self.shiftingVal.set(shifting)
            shiftingAxes = int(self.shiftingAxes.get())-1
        else:
            shiftingAxes = None
        self.parent.current.apodPreview(lor,gauss,cos2,hamming,shift,shifting,shiftingAxes)

    def stepLB(self,lorincr,gaussincr): #step linebroadening from < and > keys
        if lorincr!=0:
            self.lorScale.set(float(self.lorVal.get())+lorincr*self.lorstep)
        if gaussincr!=0:
            self.gaussScale.set(float(self.gaussVal.get())+gaussincr*self.gaussstep)

    def cancelAndClose(self):
        self.parent.current.upd()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

    def applyApodAndClose(self):
        lor = None
        gauss = None
        cos2 = None
        hamming=None
        shifting = None
        shiftingAxes = 0
        if self.lorTick.get() == 1:
            lor = safeEval(self.lorVal.get())
        if self.gaussTick.get() == 1:
            gauss = safeEval(self.gaussVal.get())
        if self.cos2Tick.get() == 1:
            cos2 = safeEval(self.cos2Val.get())
        if self.hammingTick.get() == 1:
            hamming = safeEval(self.hammingVal.get())
        shift = safeEval(self.shiftVal.get())
        if self.parent.current.data.dim > 1:
            shifting = safeEval(self.shiftingVal.get())
            self.shiftingVal.set(shifting)
            shiftingAxes = int(self.shiftingAxes.get())-1
        else:
            shiftingAxes = None
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.applyApod(lor,gauss,cos2,hamming,shift,shifting,shiftingAxes,(self.singleSlice.get()==1)))
        self.parent.menuEnable()
        self.destroy()

#######################################################################################
class SizeWindow(Toplevel): #a window for changing the size of the current dimension
    def __init__(self, parent):
        parent.menuDisable()
        #initialize variables for the widgets
        self.sizeVal = StringVar()
        self.sizeVal.set(str(parent.current.data1D.shape[-1]))
        #create a new window
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Set size")
        self.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Set size").grid(row=0,column=0,columnspan=2)
        self.sizeEntry = Entry(self.frame1,textvariable=self.sizeVal,justify="center")
        self.sizeEntry.bind("<Return>", self.sizePreview)
        self.sizeEntry.bind("<KP_Enter>", self.sizePreview)
        self.sizeEntry.grid(row=1,column=0,columnspan=2)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applySizeAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)
 
    def sizePreview(self, *args): #display the size preview from the entry widget value
        size = int(round(safeEval(self.sizeVal.get())))
        if size < 1:
            size = 1
        self.sizeVal.set(str(size))
        self.parent.current.setSizePreview(size)

    def cancelAndClose(self):
        self.parent.current.upd()
        self.parent.current.plotReset()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

    def applySizeAndClose(self):
        size = int(round(safeEval(self.sizeVal.get())))
        if size < 1:
            size = 1
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.applySize(size))
        self.parent.sideframe.upd()
        self.parent.menuEnable()
        self.destroy()

##########################################################################################
class SwapEchoWindow(Toplevel): #a window for changing the size of the current dimension
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        #initialize variables for the widgets
        self.posVal = StringVar()
        self.posVal.set(str(int(round(0.5*len(parent.current.data1D)))))
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Swap echo")
        self.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Position").grid(row=0,column=0,columnspan=2)
        self.posEntry = Entry(self.frame1,textvariable=self.posVal,justify="center")
        self.posEntry.bind("<Return>", self.swapEchoPreview)
        self.posEntry.bind("<KP_Enter>", self.swapEchoPreview)
        self.posEntry.grid(row=1,column=0,columnspan=2)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applySwapEchoAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)
        #activate the peak picking
        self.parent.current.peakPickFunc = lambda pos,self=self: self.pickedAndClose(pos) 
        self.parent.current.peakPick = True
 
    def swapEchoPreview(self, *args): #preview the swap echo result from the entry widget
        pos = int(round(safeEval(self.posVal.get())))
        if pos > 0 and pos < (self.parent.current.data1D.shape[-1]):
            self.parent.current.setSwapEchoPreview(pos)
            self.parent.current.peakPick = False
            
    def cancelAndClose(self):
        self.parent.current.peakPickReset()
        self.parent.current.upd()
        self.parent.current.plotReset()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

    def applySwapEchoAndClose(self):
        self.parent.current.peakPickReset()
        pos = int(round(safeEval(self.posVal.get())))
        if pos > 0 and pos < (self.parent.current.data1D.shape[-1]):
            self.parent.redoList = []
            self.parent.undoList.append(self.parent.current.applySwapEcho(pos))
            self.parent.bottomframe.upd()
            self.parent.menuEnable()
            self.destroy()
        else:
            print("not a valid index for swap echo")
        
    def pickedAndClose(self,pos): #apply directly if picked since another doesn't make pick doesn't make sense. find a good way to do both entry and picking in a proper way
        self.parent.current.setSwapEchoPreview(pos[0])
        self.posVal.set(str(pos[0]))
        self.parent.current.peakPick = False
        

###########################################################################
class ShiftDataWindow(Toplevel): #a window for shifting the data
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        #initialize variables for the widgets
        self.shiftVal = StringVar()
        self.shiftVal.set("0")
        self.singleSlice = IntVar()
        self.singleSlice.set(0)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Shift data")
        self.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="number of points to shift").grid(row=0,column=1)
        tk.Button(self.frame1,text="<",repeatdelay=100, repeatinterval=1,command=self.stepDownShift).grid(row=1,column=0)
        tk.Button(self.frame1,text=">",repeatdelay=100, repeatinterval=1,command=self.stepUpShift).grid(row=1,column=2)
        self.posEntry = Entry(self.frame1,textvariable=self.shiftVal,justify="center")
        self.posEntry.bind("<Return>", self.shiftPreview)
        self.posEntry.bind("<KP_Enter>", self.shiftPreview)
        self.posEntry.grid(row=1,column=1)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Checkbutton(self.frame2, text="Single slice", variable=self.singleSlice).grid(row=0,column=0,columnspan=2)
        Button(self.frame2, text="Apply",command=self.applyShiftAndClose).grid(row=1,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=1,column=1)

    def stepUpShift(self, *args):
        shift = int(round(safeEval(self.shiftVal.get())))
        shift = shift+1
        self.shiftVal.set(str(shift))
        self.shiftPreview()

    def stepDownShift(self, *args):
        shift = int(round(safeEval(self.shiftVal.get())))
        shift = shift-1
        self.shiftVal.set(str(shift))
        self.shiftPreview()

    def shiftPreview(self, *args): #preview a shifted spectrum from the entry widget
        shift = int(round(safeEval(self.shiftVal.get())))
        self.parent.current.setShiftPreview(shift)

    def cancelAndClose(self):
        self.parent.current.upd()
        self.parent.current.plotReset()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

    def applyShiftAndClose(self):
        shift = int(round(safeEval(self.shiftVal.get())))
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.applyShift(shift,(self.singleSlice.get()==1)))
        self.parent.menuEnable()
        self.destroy()

#############################################################
class DCWindow(Toplevel): #a window for changing the offset of the data
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        #initialize variables for the widgets
        self.minVal = StringVar()
        self.minVal.set(str(int(round(0.8*parent.current.data1D.shape[-1]))))
        self.maxVal = StringVar()
        self.maxVal.set(str(parent.current.data1D.shape[-1]))
        self.offsetVal = StringVar()
        self.offsetVal.set('{:.2e}'.format(parent.current.getdcOffset(int(round(0.8*parent.current.data1D.shape[-1])),parent.current.data1D.shape[-1])))
        self.singleSlice = IntVar()
        self.singleSlice.set(0)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Offset correction")
        self.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Start point").grid(row=0,column=0,columnspan=2)
        self.minEntry = Entry(self.frame1,textvariable=self.minVal,justify="center")
        self.minEntry.bind("<Return>", self.dcPreview)
        self.minEntry.bind("<KP_Enter>", self.dcPreview)
        self.minEntry.grid(row=1,column=0,columnspan=2)
        Label(self.frame1,text="End point").grid(row=2,column=0,columnspan=2)
        self.maxEntry = Entry(self.frame1,textvariable=self.maxVal,justify="center")
        self.maxEntry.bind("<Return>", self.dcPreview)
        self.maxEntry.bind("<KP_Enter>", self.dcPreview)
        self.maxEntry.grid(row=3,column=0,columnspan=2)
        Label(self.frame1,text="Offset").grid(row=4,column=0,columnspan=2)
        self.offsetEntry = Entry(self.frame1,textvariable=self.offsetVal,justify="center")
        self.offsetEntry.bind("<Return>", lambda arg: self.dcPreview(arg,True))
        self.offsetEntry.bind("<KP_Enter>", lambda arg: self.dcPreview(arg,True))
        self.offsetEntry.grid(row=5,column=0,columnspan=2)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Checkbutton(self.frame2, text="Single slice", variable=self.singleSlice).grid(row=0,column=0,columnspan=2)
        Button(self.frame2, text="Apply",command=self.applyDCAndClose).grid(row=1,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=1,column=1)
        #pick function
        self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.parent.current.peakPick = True

    def picked(self,pos,second=False): #pick a value alternating the first and second value determined by the second value.
        dataLength = self.parent.current.data1D.shape[-1]
        if second:
            minimum=int(round(safeEval(self.minVal.get())))
            if minimum < 0:
                minimum = 0
            elif minimum > dataLength:
                minimum = dataLength
            self.minVal.set(str(minimum))
            maximum=pos[0]
            self.maxVal.set(str(maximum))
            self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
            self.parent.current.peakPick = True
            val = self.parent.current.getdcOffset(minimum,maximum)
            self.offsetVal.set('{:.2e}'.format(val))
            self.parent.current.dcOffset(val)
        else:
            self.minVal.set(str(pos[0]))
            self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos,True) 
            self.parent.current.peakPick = True
            maximum=int(round(safeEval(self.maxVal.get())))
            if maximum < 0:
                maximum = 0
            elif maximum > dataLength:
                maximum = dataLength
            minimum = pos[0]
            val = self.parent.current.getdcOffset(minimum,maximum)
            self.offsetVal.set('{:.2e}'.format(val))

    def dcPreview(self, arg, inserted=False): #preview the dc offset correction
        if inserted:
            self.parent.current.dcOffset(safeEval(self.offsetVal.get()))
        else:
            dataLength = self.parent.current.data1D.shape[-1]
            minimum = int(round(safeEval(self.minVal.get())))
            if minimum < 0:
                minimum = 0
            elif minimum > dataLength:
                minimum = dataLength
            self.minVal.set(str(minimum))
            maximum = int(round(safeEval(self.maxVal.get())))
            if maximum < 0:
                maximum = 0
            elif maximum > dataLength:
                maximum = dataLength
            self.maxVal.set(str(maximum))
            val = self.parent.current.getdcOffset(minimum,maximum)
            self.offsetVal.set('{:.2e}'.format(val))
            self.parent.current.dcOffset(val)

    def cancelAndClose(self):
        self.parent.current.peakPickReset()
        self.parent.current.upd()
        self.parent.current.plotReset()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

    def applyDCAndClose(self):
        self.parent.current.peakPickReset()
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.subtract(safeEval(self.offsetVal.get()),(self.singleSlice.get()==1)))
        self.parent.menuEnable()
        self.destroy()

#############################################################
class BaselineWindow(Toplevel):
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.degreeVal = StringVar()
        self.degreeVal.set('3')
        self.singleSlice = IntVar()
        self.singleSlice.set(0)
        self.removeList = []
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Baseline correction")
        self.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Polynomial Degree:").grid(row=0,column=0,columnspan=2)
        self.degreeEntry = Entry(self.frame1,textvariable=self.degreeVal,justify="center")
        self.degreeEntry.bind("<Return>", self.setDegree)
        self.degreeEntry.bind("<KP_Enter>", self.setDegree)
        self.degreeEntry.grid(row=1,column=0,columnspan=2)
        Button(self.frame1, text="Fit",command=self.preview).grid(row=2,column=0)
        Button(self.frame1, text="Reset",command=self.reset).grid(row=2,column=1)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Checkbutton(self.frame2, text="Single slice", variable=self.singleSlice).grid(row=0,column=0,columnspan=2)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=1,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=1,column=1)
        self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.parent.current.peakPick = True

    def picked(self,pos):
        self.removeList.append(pos[0])
        self.parent.current.previewRemoveList(self.removeList)
        self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.parent.current.peakPick = True

    def setDegree(self):
        degree = safeEval(self.degreeVal.get())
        self.degreeVal.set(str(degree))
        
    def preview(self, *args): 
        degree = safeEval(self.degreeVal.get())
        self.degreeVal.set(str(degree))
        self.parent.current.previewBaseline(degree,self.removeList)
        self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.parent.current.peakPick = True

    def reset(self, *args):
        self.removeList = []
        self.parent.current.resetPreviewRemoveList()

    def cancelAndClose(self):
        self.parent.current.peakPickReset()
        self.parent.current.resetPreviewRemoveList()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

    def applyAndClose(self):
        degree = safeEval(self.degreeVal.get())
        self.parent.current.peakPickReset()
        self.parent.current.resetPreviewRemoveList()
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.applyBaseline(degree,self.removeList,(self.singleSlice.get()==1)))
        self.parent.current.upd()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()
        
#############################################################
class regionWindow(Toplevel): #A general region selection frame
    def __init__(self, parent,name):
        parent.menuDisable()
        Toplevel.__init__(self)
        #initialize variables for the widgets
        self.minVal = StringVar()
        self.minVal.set("0")
        self.maxVal = StringVar()
        self.maxVal.set(str(parent.current.data1D.shape[-1]))
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title(name)
        self.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Start point").grid(row=0,column=0,columnspan=2)
        self.minEntry = Entry(self.frame1,textvariable=self.minVal,justify="center")
        self.minEntry.bind("<Return>", self.checkValues)
        self.minEntry.bind("<KP_Enter>", self.checkValues)
        self.minEntry.grid(row=1,column=0,columnspan=2)
        Label(self.frame1,text="End point").grid(row=2,column=0,columnspan=2)
        self.maxEntry = Entry(self.frame1,textvariable=self.maxVal,justify="center")
        self.maxEntry.bind("<Return>", self.checkValues)
        self.maxEntry.bind("<KP_Enter>", self.checkValues)
        self.maxEntry.grid(row=3,column=0,columnspan=2)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)
        #pick function
        self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos)
        self.parent.current.peakPick = True

    def picked(self,pos,second=False): #pick a value alternating the first and second value determined by the second value.
        if second:
            dataLength = self.parent.current.data1D.shape[-1]
            minimum=int(round(safeEval(self.minVal.get())))
            if minimum < 0:
                minimum = 0
            elif minimum > dataLength:
                minimum = dataLength
            self.minVal.set(str(minimum))
            maximum=pos[0]
            self.maxVal.set(str(maximum))
            self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
            self.parent.current.peakPick = True
        else:
            self.minVal.set(str(pos[0]))
            self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos,True) 
            self.parent.current.peakPick = True

    def checkValues(self, *args): #not really preview but just to check the inputs
        dataLength = self.parent.current.data1D.shape[-1]
        minimum = int(round(safeEval(self.minVal.get())))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minVal.set(str(minimum))
        maximum = int(round(safeEval(self.maxVal.get())))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxVal.set(str(maximum))

    def cancelAndClose(self):
        self.parent.current.peakPickReset()
        self.parent.updAllFrames()
        self.parent.menuEnable()
        self.destroy()

    def apply(self,maximum,minimum):
        pass
        
    def applyAndClose(self):
        self.parent.current.peakPickReset()
        dataLength = self.parent.current.data1D.shape[-1]
        minimum = int(round(safeEval(self.minVal.get())))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        maximum = int(round(safeEval(self.maxVal.get())))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.apply(maximum,minimum)
        self.parent.menuEnable()
        self.destroy()
        
############################################################
class integrateWindow(regionWindow): #A window for obtaining the integral of a selected region
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Integrate')

    def apply(self,maximum,minimum):
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.integrate(minimum,maximum))
        self.parent.updAllFrames()
        
############################################################
class maxWindow(regionWindow): #A window for obtaining the max of a selected region
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Max')

    def apply(self,maximum,minimum):
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.maxMatrix(minimum,maximum))
        self.parent.updAllFrames()

############################################################
class minWindow(regionWindow): #A window for obtaining the min of a selected region
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Min')

    def apply(self,maximum,minimum):
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.minMatrix(minimum,maximum))
        self.parent.updAllFrames()
        
############################################################
class argmaxWindow(regionWindow): #A window for obtaining the max of a selected region
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Max position')

    def apply(self,maximum,minimum):
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.argmaxMatrix(minimum,maximum))
        self.parent.updAllFrames()

############################################################
class argminWindow(regionWindow): #A window for obtaining the min of a selected region
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Min position')

    def apply(self,maximum,minimum):
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.argminMatrix(minimum,maximum))
        self.parent.updAllFrames()
        
############################################################
class extractRegionWindow(regionWindow): #A window for obtaining a selected region
    def __init__(self, parent):
        regionWindow.__init__(self,parent,'Extract part')

    def apply(self,maximum,minimum):
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.getRegion(minimum,maximum))
        self.parent.updAllFrames()

##############################################################
class DeleteWindow(Toplevel):
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Delete")
        self.resizable(width=FALSE, height=FALSE)
        #initialize variables for the widgets
        self.pos = StringVar()
        self.pos.set('0')
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Indexes to delete").grid(row=0,column=0)
        self.posEntry = Entry(self.frame1,textvariable=self.pos,justify="center")
        self.posEntry.bind("<Return>", self.preview)
        self.posEntry.bind("<KP_Enter>", self.preview)
        self.posEntry.grid(row=1,column=0)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)

    def preview(self, *args):
        env = vars(np).copy()
        length = int(self.parent.current.data1D.shape[-1])
        env['length']=length # so length can be used to in equations
        pos=np.array(eval(self.pos.get(),env))                # find a better solution, also add catch for exceptions
        if (pos > -1).all() and (pos < length).all():
            self.parent.current.deletePreview(pos)
        else:
            print('Not all values are valid indexes to delete')
        
    def applyAndClose(self):
        env = vars(np).copy()
        length = int(self.parent.current.data1D.shape[-1])
        env['length']=length # so length can be used to in equations
        pos=np.array(eval(self.pos.get(),env))                # find a better solution, also add catch for exceptions
        if (pos > -1).all() and (pos < length).all():
            self.parent.redoList = []
            self.parent.undoList.append(self.parent.current.delete(pos))
            self.parent.menuEnable()
            self.parent.sideframe.upd()
            self.destroy()
        else:
            print('Not all values are valid indexes to delete')
        
    def cancelAndClose(self):
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

##############################################################
class SplitWindow(Toplevel):
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Split")
        self.resizable(width=FALSE, height=FALSE)
        #initialize variables for the widgets
        self.section = StringVar()
        self.section.set('1')
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Number of sections").grid(row=0,column=0)
        self.sectionEntry = Entry(self.frame1,textvariable=self.section,justify="center")
        self.sectionEntry.bind("<Return>", self.preview)
        self.sectionEntry.bind("<KP_Enter>", self.preview)
        self.sectionEntry.grid(row=1,column=0)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)

    def preview(self, *args):
        val = int(safeEval(self.section.get()))
        self.section.set(str(val))
        
    def applyAndClose(self):
        val = int(safeEval(self.section.get()))
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.split(val))
        self.parent.menuEnable()
        self.parent.updAllFrames()
        self.destroy()
        
    def cancelAndClose(self):
        self.parent.menuEnable()
        self.destroy()

##############################################################
class ConcatenateWindow(Toplevel):
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Concatenate")
        self.resizable(width=FALSE, height=FALSE)
        #initialize variables for the widgets
        self.axes = StringVar()
        self.axes.set('0')
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Concatenation axes").grid(row=0,column=0)
        self.axesEntry = Entry(self.frame1,textvariable=self.axes,justify="center")
        self.axesEntry.bind("<Return>", self.preview)
        self.axesEntry.bind("<KP_Enter>", self.preview)
        self.axesEntry.grid(row=1,column=0)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)

    def preview(self, *args):
        val = int(safeEval(self.axes.get()))
        if val < 0:
            val = 0
        elif val > (self.parent.current.data.dim-2):
            val = self.parent.current.data.dim-2
        self.axes.set(str(val))
        
    def applyAndClose(self):
        val = int(safeEval(self.axes.get()))
        if val < 0:
            print("concatenate axes is too small")
            return
        elif val > (self.parent.current.data.dim-2):
            print("concatenate axes is too large")
            return
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.concatenate(val))
        self.parent.menuEnable()
        self.parent.updAllFrames()
        self.destroy()
        
    def cancelAndClose(self):
        self.parent.menuEnable()
        self.destroy()
        
##############################################################
class InsertWindow(Toplevel):
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Insert")
        self.resizable(width=FALSE, height=FALSE)
        #initialize variables for the widgets
        self.pos = StringVar()
        self.pos.set(str(self.parent.current.data1D.shape[-1]))
        self.ws = StringVar()
        self.ws.set(self.parent.mainProgram.workspaceNames[0])
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Start insert at index").grid(row=0,column=0)
        self.posEntry = Entry(self.frame1,textvariable=self.pos,justify="center")
        self.posEntry.bind("<Return>", self.preview)
        self.posEntry.bind("<KP_Enter>", self.preview)
        self.posEntry.grid(row=1,column=0)
        Label(self.frame1,text="Workspace to insert").grid(row=2,column=0)
        OptionMenu(self.frame1,self.ws,self.parent.mainProgram.workspaceNames[0],*self.parent.mainProgram.workspaceNames).grid(row=3,column=0)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)

    def preview(self, *args):
        pos = int(round(safeEval(self.pos.get())))
        if pos > self.parent.current.data1D.shape[-1]:
            pos = self.parent.current.data1D.shape[-1]
        elif pos < 0:
            pos = 0
        self.pos.set(str(pos))
        
    def applyAndClose(self):
        pos = int(round(safeEval(self.pos.get())))
        if pos > self.parent.current.data1D.shape[-1]:
            pos = self.parent.current.data1D.shape[-1]
        elif pos < 0:
            pos = 0
        ws = self.parent.mainProgram.workspaceNames.index(self.ws.get())
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.insert(self.parent.mainProgram.workspaces[ws].masterData.data,pos))
        self.parent.menuEnable()
        self.parent.sideframe.upd()
        self.destroy()
        
    def cancelAndClose(self):
        self.parent.menuEnable()
        self.destroy()

##############################################################
class AddWindow(Toplevel):
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Add")
        self.resizable(width=FALSE, height=FALSE)
        #initialize variables for the widgets
        self.ws = StringVar()
        self.ws.set(self.parent.mainProgram.workspaceNames[0])
        self.singleSlice = IntVar()
        self.singleSlice.set(0)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Workspace to add").grid(row=0,column=0)
        OptionMenu(self.frame1,self.ws,self.parent.mainProgram.workspaceNames[0],*self.parent.mainProgram.workspaceNames).grid(row=1,column=0)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Checkbutton(self.frame2, text="Single slice", variable=self.singleSlice).grid(row=0,column=0,columnspan=2)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=1,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=1,column=1)
        
    def applyAndClose(self):
        ws = self.parent.mainProgram.workspaceNames.index(self.ws.get())
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.add(self.parent.mainProgram.workspaces[ws].masterData.data,(self.singleSlice.get()==1)))
        self.parent.menuEnable()
        self.parent.sideframe.upd()
        self.destroy()
        
    def cancelAndClose(self):
        self.parent.menuEnable()
        self.destroy()
        
##############################################################
class SubtractWindow(Toplevel):
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Subtract")
        self.resizable(width=FALSE, height=FALSE)
        #initialize variables for the widgets
        self.ws = StringVar()
        self.ws.set(self.parent.mainProgram.workspaceNames[0])
        self.singleSlice = IntVar()
        self.singleSlice.set(0)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Workspace to subtract").grid(row=0,column=0)
        OptionMenu(self.frame1,self.ws,self.parent.mainProgram.workspaceNames[0],*self.parent.mainProgram.workspaceNames).grid(row=1,column=0)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Checkbutton(self.frame2, text="Single slice", variable=self.singleSlice).grid(row=0,column=0,columnspan=2)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=1,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=1,column=1)
        
    def applyAndClose(self):
        ws = self.parent.mainProgram.workspaceNames.index(self.ws.get())
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.subtract(self.parent.mainProgram.workspaces[ws].masterData.data,(self.singleSlice.get()==1)))
        self.parent.menuEnable()
        self.parent.sideframe.upd()
        self.destroy()
        
    def cancelAndClose(self):
        self.parent.menuEnable()
        self.destroy()

##############################################################
class SNWindow(Toplevel):
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("S/N")
        self.resizable(width=FALSE, height=FALSE)
        #initialize variables for the widgets
        self.minNoiseVal = StringVar()
        self.minNoiseVal.set("0")
        self.maxNoiseVal = StringVar()
        self.maxNoiseVal.set(str(parent.current.data1D.shape[-1]))
        self.minVal = StringVar()
        self.minVal.set("0")
        self.maxVal = StringVar()
        self.maxVal.set(str(parent.current.data1D.shape[-1]))
        self.result = StringVar()
        self.result.set('0.0')
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Start point noise").grid(row=0,column=0,columnspan=2)
        self.minNoiseEntry = Entry(self.frame1,textvariable=self.minNoiseVal,justify="center")
        self.minNoiseEntry.bind("<Return>", self.checkValues)
        self.minNoiseEntry.bind("<KP_Enter>", self.checkValues)
        self.minNoiseEntry.grid(row=1,column=0,columnspan=2)
        Label(self.frame1,text="End point noise").grid(row=2,column=0,columnspan=2)
        self.maxNoiseEntry = Entry(self.frame1,textvariable=self.maxNoiseVal,justify="center")
        self.maxNoiseEntry.bind("<Return>", self.checkValues)
        self.maxNoiseEntry.bind("<KP_Enter>", self.checkValues)
        self.maxNoiseEntry.grid(row=3,column=0,columnspan=2)
        Label(self.frame1,text="Start point Signal").grid(row=4,column=0,columnspan=2)
        self.minEntry = Entry(self.frame1,textvariable=self.minVal,justify="center")
        self.minEntry.bind("<Return>", self.checkValues)
        self.minEntry.bind("<KP_Enter>", self.checkValues)
        self.minEntry.grid(row=5,column=0,columnspan=2)
        Label(self.frame1,text="End point Signal").grid(row=6,column=0,columnspan=2)
        self.maxEntry = Entry(self.frame1,textvariable=self.maxVal,justify="center")
        self.maxEntry.bind("<Return>", self.checkValues)
        self.maxEntry.bind("<KP_Enter>", self.checkValues)
        self.maxEntry.grid(row=7,column=0,columnspan=2)
        Label(self.frame1,text="S/N").grid(row=8,column=0,columnspan=2)
        Entry(self.frame1,textvariable=self.result,justify="center").grid(row=9,column=0,columnspan=2)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.apply).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)
        self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos)
        self.parent.current.peakPick = True
        
    def picked(self,pos,num=0): 
        if num == 0:
            self.minNoiseVal.set(str(pos[0]))
            self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos,1) 
            self.parent.current.peakPick = True
        elif num == 1:
            self.maxNoiseVal.set(str(pos[0]))
            self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos,2) 
            self.parent.current.peakPick = True
        elif num == 2:
            self.minVal.set(str(pos[0]))
            self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos,3) 
            self.parent.current.peakPick = True
        elif num == 3:
            self.maxVal.set(str(pos[0]))
            self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos,0) 
            self.parent.current.peakPick = True
            
    def checkValues(self, *args): 
        dataLength = self.parent.current.data1D.shape[-1]
        minimum = int(round(safeEval(self.minNoiseVal.get())))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minNoiseVal.set(str(minimum))
        maximum = int(round(safeEval(self.maxNoiseVal.get())))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxNoiseVal.set(str(maximum))
        minimum = int(round(safeEval(self.minVal.get())))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minVal.set(str(minimum))
        maximum = int(round(safeEval(self.maxVal.get())))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxVal.set(str(maximum))
        
    def apply(self):
        dataLength = self.parent.current.data1D.shape[-1]
        minimumNoise = int(round(safeEval(self.minNoiseVal.get())))
        if minimumNoise < 0:
            minimumNoise = 0
        elif minimumNoise > dataLength:
            minimumNoise = dataLength
        self.minNoiseVal.set(str(minimumNoise))
        maximumNoise = int(round(safeEval(self.maxNoiseVal.get())))
        if maximumNoise < 0:
            maximumNoise = 0
        elif maximumNoise > dataLength:
            maximumNoise = dataLength
        self.maxNoiseVal.set(str(maximumNoise))
        minimum = int(round(safeEval(self.minVal.get())))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minVal.set(str(minimum))
        maximum = int(round(safeEval(self.maxVal.get())))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxVal.set(str(maximum))
        self.result.set(str(self.parent.current.SN(minimumNoise,maximumNoise,minimum,maximum)))
        
    def cancelAndClose(self):
        self.parent.current.peakPickReset()
        self.parent.menuEnable()
        self.destroy()
        
##############################################################
class FWHMWindow(Toplevel):
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("FWHM")
        self.resizable(width=FALSE, height=FALSE)
        #initialize variables for the widgets
        self.minVal = StringVar()
        self.minVal.set("0")
        self.maxVal = StringVar()
        self.maxVal.set(str(parent.current.data1D.shape[-1]))
        self.result = StringVar()
        self.result.set('0.0')
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Start").grid(row=0,column=0,columnspan=2)
        self.minEntry = Entry(self.frame1,textvariable=self.minVal,justify="center")
        self.minEntry.bind("<Return>", self.checkValues)
        self.minEntry.bind("<KP_Enter>", self.checkValues)
        self.minEntry.grid(row=1,column=0,columnspan=2)
        Label(self.frame1,text="End").grid(row=2,column=0,columnspan=2)
        self.maxEntry = Entry(self.frame1,textvariable=self.maxVal,justify="center")
        self.maxEntry.bind("<Return>", self.checkValues)
        self.maxEntry.bind("<KP_Enter>", self.checkValues)
        self.maxEntry.grid(row=3,column=0,columnspan=2)
        Label(self.frame1,text="FWHM").grid(row=4,column=0,columnspan=2)
        Entry(self.frame1,textvariable=self.result,justify="center").grid(row=5,column=0,columnspan=2)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.apply).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)
        self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos)
        self.parent.current.peakPick = True
        
    def picked(self,pos,num=0): 
        if num == 0:
            self.minVal.set(str(pos[0]))
            self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos,1) 
            self.parent.current.peakPick = True
        elif num == 1:
            self.maxVal.set(str(pos[0]))
            self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos,0) 
            self.parent.current.peakPick = True
            self.apply()
            
    def checkValues(self, *args): 
        dataLength = self.parent.current.data1D.shape[-1]
        minimum = int(round(safeEval(self.minVal.get())))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minVal.set(str(minimum))
        maximum = int(round(safeEval(self.maxVal.get())))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxVal.set(str(maximum))
        
    def apply(self):
        dataLength = self.parent.current.data1D.shape[-1]
        minimum = int(round(safeEval(self.minVal.get())))
        if minimum < 0:
            minimum = 0
        elif minimum > dataLength:
            minimum = dataLength
        self.minVal.set(str(minimum))
        maximum = int(round(safeEval(self.maxVal.get())))
        if maximum < 0:
            maximum = 0
        elif maximum > dataLength:
            maximum = dataLength
        self.maxVal.set(str(maximum))
        self.result.set(str(self.parent.current.fwhm(minimum,maximum)))
        
    def cancelAndClose(self):
        self.parent.current.peakPickReset()
        self.parent.menuEnable()
        self.destroy()
        
################################################################
class ShearingWindow(Toplevel): 
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Shearing")
        self.resizable(width=FALSE, height=FALSE)
        #initialize variables for the widgets
        self.shear = StringVar()
        self.shear.set('0.0')
        options = list(map(str,range(1,self.parent.masterData.dim+1)))
        self.axes = StringVar()
        self.axes.set(str(self.parent.masterData.dim-1))
        self.axes2 = StringVar()
        self.axes2.set(str(self.parent.masterData.dim))
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Shearing constant").grid(row=0,column=0)
        self.shearEntry = Entry(self.frame1,textvariable=self.shear,justify="center")
        self.shearEntry.bind("<Return>", self.shearPreview)
        self.shearEntry.bind("<KP_Enter>", self.shearPreview)
        self.shearEntry.grid(row=1,column=0)
        Label(self.frame1,text="Shearing direction").grid(row=2,column=0)
        OptionMenu(self.frame1,self.axes, self.axes.get(),*options).grid(row=3,column=0)
        Label(self.frame1,text="Shearing axis").grid(row=4,column=0)
        OptionMenu(self.frame1,self.axes2,self.axes2.get(),*options).grid(row=5,column=0)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)

    def shearPreview(self, *args):
        shear = float(safeEval(self.shear.get()))
        self.shear.set(str(shear))

    def cancelAndClose(self):
        self.parent.menuEnable()
        self.destroy()

    def applyAndClose(self):
        shear = float(safeEval(self.shear.get()))
        axes = int(self.axes.get())-1
        axes2 = int(self.axes2.get())-1
        if axes == axes2:
            print("Axes can't be the same for shearing")
        else:
            self.parent.redoList = []
            self.parent.undoList.append(self.parent.current.shearing(shear,axes,axes2))
            self.parent.menuEnable()
            self.destroy()

##########################################################################################
class MultiplyWindow(Toplevel): 
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        #initialize variables for the widgets
        self.val = StringVar()
        self.singleSlice = IntVar()
        self.singleSlice.set(0)
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Multiply")
        self.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Values:").grid(row=0,column=0,columnspan=2)
        self.minEntry = Entry(self.frame1,textvariable=self.val,justify="center")
        self.minEntry.bind("<Return>", self.preview)
        self.minEntry.bind("<KP_Enter>", self.preview)
        self.minEntry.grid(row=1,column=0,columnspan=2)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Checkbutton(self.frame2, text="Single slice", variable=self.singleSlice).grid(row=0,column=0,columnspan=2)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=1,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=1,column=1)

    def preview(self, *args):
        env = vars(np).copy()
        env['length']=int(self.parent.current.data1D.shape[-1]) # so length can be used to in equations
        env['euro']=lambda fVal, num=int(self.parent.current.data1D.shape[-1]): euro(fVal,num)
        val=eval(self.val.get(),env)                # find a better solution, also add catch for exceptions          
        self.parent.current.multiplyPreview(np.array(val))

    def cancelAndClose(self):
        self.parent.current.upd()
        self.parent.current.plotReset()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

    def applyAndClose(self):
        env = vars(np).copy()
        env['length']=int(self.parent.current.data1D.shape[-1]) # so length can be used to in equations
        env['euro']=lambda fVal, num=int(self.parent.current.data1D.shape[-1]): euro(fVal,num)
        val=eval(self.val.get(),env)                # find a better solution, also add catch for exceptions
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.multiply(np.array(val),(self.singleSlice.get()==1)))
        self.parent.menuEnable()
        self.destroy()

            
##########################################################################################
class XaxWindow(Toplevel): #a window for setting the xax of the current data
    def __init__(self, parent):
        parent.menuDisable()
        Toplevel.__init__(self)
        #initialize variables for the widgets
        self.val = StringVar()
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("User defined x-axis")
        self.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Expression for x-axis values").grid(row=0,column=0,columnspan=2)
        self.minEntry = Entry(self.frame1,textvariable=self.val,justify="center")
        self.minEntry.bind("<Return>", self.xaxPreview)
        self.minEntry.bind("<KP_Enter>", self.xaxPreview)
        self.minEntry.grid(row=1,column=0,columnspan=2)

        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyXaxAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)

    def xaxPreview(self, *args):
        env = vars(np).copy()
        env['length']=int(self.parent.current.data1D.shape[-1]) # so length can be used to in equations
        env['euro']=lambda fVal, num=int(self.parent.current.data1D.shape[-1]): euro(fVal,num)
        val=eval(self.val.get(),env)                # find a better solution, also add catch for exceptions          
        if isinstance(val,(list,np.ndarray)):
            if len(val)==self.parent.current.data1D.shape[-1]:
                if all(isinstance(x,(int,float)) for x in val):
                    self.parent.current.setXaxPreview(np.array(val))
                else:
                    print("Array is not all of int or float type")
            else:
                print("Length of input does not match length of data")
        else:
            print("Input is not a list or array")

    def cancelAndClose(self):
        self.parent.current.upd()
        self.parent.current.plotReset()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

    def applyXaxAndClose(self):
        env = vars(np).copy()
        env['length']=int(self.parent.current.data1D.shape[-1]) # so length can be used to in equations
        env['euro']=lambda fVal, num=int(self.parent.current.data1D.shape[-1]): euro(fVal,num)
        val=eval(self.val.get(),env)                # find a better solution, also add catch for exceptions
        if isinstance(val,(list,np.ndarray)):
            if len(val)==self.parent.current.data1D.shape[-1]:
                if all(isinstance(x,(int,float)) for x in val):
                    self.parent.redoList = []
                    self.parent.undoList.append(self.parent.current.setXax(np.array(val)))
                    self.parent.menuEnable()
                    self.destroy()
                else:
                    print("Array is not all of int or float type")
            else:
                print("Length of input does not match length of data")
        else:
            print("Input is not a list or array")

##########################################################################################
class RefWindow(Toplevel): #a window for setting the ppm reference
    def __init__(self, parent):
        parent.menuDisable()
        if parent.current.spec == 0:
            print('Setting ppm is only available for frequency data')
        Toplevel.__init__(self)
        #initialize variables for the widgets
        self.freqVal = StringVar()
        self.freqVal.set(str(parent.current.freq))
        self.refVal = StringVar()
        self.refVal.set('0.0')
        self.parent = parent
        self.geometry('+0+0')
        self.transient(self.parent)
        self.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.title("Set ppm reference")
        self.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Frequency").grid(row=0,column=0)
        self.freqEntry = Entry(self.frame1,textvariable=self.freqVal,justify="center")
        self.freqEntry.bind("<Return>", self.preview)
        self.freqEntry.bind("<KP_Enter>", self.preview)
        self.freqEntry.grid(row=1,column=0)
        Label(self.frame1,text="ppm").grid(row=2,column=0)
        self.refEntry = Entry(self.frame1,textvariable=self.refVal,justify="center")
        self.refEntry.bind("<Return>", self.preview)
        self.refEntry.bind("<KP_Enter>", self.preview)
        self.refEntry.grid(row=3,column=0)
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)
        #activate the peak picking
        self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos)
        self.parent.current.peakPick = True
 
    def preview(self, *args): #fix the input values
        freq = safeEval(self.freqVal.get())
        ref = safeEval(self.refVal.get())
        self.freqVal.set(str(freq))
        self.refVal.set(str(ref))

    def cancelAndClose(self):
        self.parent.current.peakPickReset()
        self.parent.current.showFid()
        self.parent.menuEnable()
        self.destroy()

    def applyAndClose(self):
        self.parent.current.peakPickReset()
        freq = safeEval(self.freqVal.get())
        ref = safeEval(self.refVal.get())
        self.parent.redoList = []
        self.parent.undoList.append(self.parent.current.setRef(freq/(1.0+ref*1e-6)))
        self.parent.menuEnable()
        self.destroy()
        
    def picked(self,pos): 
        self.freqVal.set(str(self.parent.current.freq+self.parent.current.xax[pos[0]]))
        self.parent.current.peakPickFunc = lambda pos,self=self: self.picked(pos)
        self.parent.current.peakPick = True
            
#################################################################################    
#the main program
if __name__ == "__main__":
    root = Tk()
    img = tk.PhotoImage(file=os.path.dirname(os.path.realpath(__file__))+'/logo.gif')
    root.tk.call('wm', 'iconphoto', root._w, img)
    MainProgram(root)
    root.title("ssNake") 
    root.style = Style()
    root.style.theme_use("clam")
    #root.attributes('-zoomed', True)
    root.mainloop()
