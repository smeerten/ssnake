#!/usr/bin/env python

import numpy as np
import sys
if sys.version_info >= (3,0):
    from tkinter import *
    import tkinter as tk
    from tkinter.ttk import *
    from tkinter.filedialog import askopenfilename
    from tkinter.filedialog import asksaveasfile
    from tkinter.filedialog import asksaveasfilename
else:
    from Tkinter import *
    import Tkinter as tk
    from ttk import *
    from tkFileDialog   import askopenfilename
    from tkFileDialog   import asksaveasfile
    from tkFileDialog   import asksaveasfilename
import spectrum_classes as sc
import fitting as fit
import math
#For varianload
import os
from struct import unpack
#------------
from safeEval import safeEval

pi=math.pi

#one window to rule them all
class Main1DWindow(Frame):
    def __init__(self,parent):
        Frame.__init__(self,parent)
        self.undoList = [] #the list to hold all the undo lambda functions
        self.redoList = [] #the list to hold all the redo lambda functions
        self.parent = parent #remember your parents
        #create the menu
        self.menubar = Menu(self.parent)
        self.parent.config(menu=self.menubar)
        #the file drop down menu
        filemenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="New", command=self.NewFile)
        filemenu.add_command(label="Open...", command=self.OpenFile)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.quit)
        #the hotkeys for different commands
        self.bind_all("<Control-q>", lambda extra: self.quit())
        self.bind_all("<Control-z>", self.undo)
        self.bind_all("<Control-y>", self.redo)
	#the load drop down menu
        loadmenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Load", menu=loadmenu)
        loadmenu.add_command(label="Load Varian data", command=self.LoadVarianFile)
        loadmenu.add_command(label="Load Bruker Topspin/XWinNMR", command=self.LoadBrukerTopspin)
        loadmenu.add_command(label="Load Chemagnetics data", command=self.LoadChemFile)
        loadmenu.add_command(label="Load Simpson data", command=self.LoadSimpsonFile)
        
        #the save drop down menu
        savemenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Save", menu=savemenu)
        savemenu.add_command(label="Save figure", command=self.saveFigure)
        savemenu.add_command(label="Save as Simpson data", command=self.SaveSimpsonFile)
        
	#the edit drop down menu
        editmenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Edit", menu=editmenu)
        editmenu.add_command(label="Undo", command=self.undo)
        editmenu.add_command(label="Redo", command=self.redo)

	#the tool drop down menu
        toolMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Tools",menu=toolMenu)
        toolMenu.add_command(label="Real", command=self.real)
        toolMenu.add_command(label="Imag", command=self.imag)
        toolMenu.add_command(label="Abs", command=self.abs) 
        toolMenu.add_command(label="Apodize", command=self.createApodWindow)
        toolMenu.add_command(label="Phasing", command=self.createPhaseWindow)
        toolMenu.add_command(label="Sizing", command=self.createSizeWindow) 
        toolMenu.add_command(label="Swap Echo", command=self.createSwapEchoWindow)
        toolMenu.add_command(label="Shift Data", command=self.createShiftDataWindow)
        toolMenu.add_command(label="DC offset correction", command=self.createDCWindow)
        toolMenu.add_command(label="Correct Bruker digital filter", command=self.BrukerDigital)

        #the matrix drop down menu
        matrixMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Matrix",menu=matrixMenu)
        matrixMenu.add_command(label="Integrate", command=self.createIntegrateWindow)
        matrixMenu.add_command(label="Max", command=self.createMaxWindow)
        matrixMenu.add_command(label="Min", command=self.createMinWindow)
        matrixMenu.add_command(label="Extract part", command=self.createRegionWindow)
        matrixMenu.add_command(label="Flip L/R", command=self.flipLR)
        matrixMenu.add_command(label="Shearing", command=self.createShearingWindow)
        
        #the fft drop down menu
        fftMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Fourier",menu=fftMenu)
        fftMenu.add_command(label="Fourier transform", command=self.fourier)
        fftMenu.add_command(label="Fftshift", command=self.fftshift)
        fftMenu.add_command(label="Inv fftshift", command=self.invFftshift)
        fftMenu.add_command(label="Hilbert transform", command=self.hilbert)

	#the fitting drop down menu
        fittingMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Fitting",menu=fittingMenu)
        fittingMenu.add_command(label="Relaxation Curve", command=self.createRelaxWindow)
        fittingMenu.add_command(label="Peak Deconvolution", command=self.createPeakDeconvWindow)

	#the plot drop down menu
        plotMenu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Plot",menu=plotMenu)
        plotMenu.add_command(label="1D plot", command=self.plot1D)
        plotMenu.add_command(label="Stack plot", command=self.plotStack)
        plotMenu.add_command(label="Array plot", command=self.plotArray)
        plotMenu.add_command(label="Contour plot", command=self.plotContour)
        plotMenu.add_command(label="Skewed plot", command=self.plotSkewed)
        plotMenu.add_command(label="Set reference", command=self.createRefWindow)
        plotMenu.add_command(label="User x-axis", command=self.createXaxWindow)

        x=np.linspace(0,2*np.pi*10,1000)[:-1] #fake data
        x2=np.linspace(0,2*np.pi*10,200)[1:] #fake data
        test=np.exp(-1j*x)*np.exp(-1*x/10.0)#fake data
        test2=1-np.exp(-x2)
        #self.masterData=sc.Spectrum(np.array([np.array([test,test*2,test*3]),np.array([test*3,test*4,test*5])]),[600000000.0,500000000.0,400000000.0],[1000.0,2000.0,3000.0])#create a Spectrum instance with the fake data
        self.masterData=sc.Spectrum(np.outer(test2,test),[600000000.0,500000000.0],[1000.0,2000.0])
        self.current=sc.Current1D(self,self.masterData) 
        self.current.grid(row=0,column=0,sticky="nswe")
	#create the sideframe, bottomframe and textframe
        self.sideframe=SideFrame(self) 
        self.sideframe.grid(row=0,column=2,sticky='n')
        Separator(self,orient=VERTICAL).grid(row=0,column=1,rowspan=4,sticky='ns')
        self.bottomframe=BottomFrame(self)
        self.bottomframe.grid(row=1,column=0,sticky='w') 
        Separator(self,orient=HORIZONTAL).grid(row=2,sticky='ew')
        self.textframe=TextFrame(self)
        self.textframe.grid(row=3,column=0,sticky='s')  

        #all the functions that will be called from the menu and the extra frames

    def menuEnable(self):
        for i in range(12):
            self.menubar.entryconfig(i,state='normal')
        self.sideframe.frameEnable()
        self.bottomframe.frameEnable()
        self.textframe.frameEnable()

    def menuDisable(self):
        for i in range(12):
            self.menubar.entryconfig(i,state='disabled')
        self.sideframe.frameDisable()
        self.bottomframe.frameDisable()
        self.textframe.frameDisable()

    def NewFile(self):
        print("New File!") #to be added

    def OpenFile(self):
        name = askopenfilename() #to be added
        print(name)

    def LoadVarianFile(self):
        FilePath = askopenfilename()
        if FilePath is not '': #if not canceled
            Dir = os.path.dirname(FilePath) #convert path to file to path of folder
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
                try:
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
                        fid32 = bin(status)[-3] #check if 32 bits, or float
                        fidfloat = bin(status)[-4]
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
                    self.current.grid_remove()
                    self.current.destroy()
                    if SizeTD1 is 1: #convert to 1D dat if the data is 1D (so no 1xnp data, but np)
                        fid = fid[0][:]
                        self.masterData=sc.Spectrum(fid,[freq],[sw])
                        self.current=sc.Current1D(self,self.masterData) #create the Current1D instance from the Spectrum  
                    else: #For 2D data
                        self.masterData=sc.Spectrum(fid,[freq]*2,[sw]*2)
                        self.current=sc.Current1D(self,self.masterData) #create the Current1D instance from the Spectrum  
                    #add some check to see if current exists
                    self.current.grid(row=0,column=0,sticky="nswe")
                    self.updAllFrames()
                except:
                    print('Error loading Varian data from '+Dir+os.path.sep+'fid. No data loaded!')
            else: #If /fid does not exist
                print(Dir+os.path.sep+'fid does not exits, no Varian data loaded!')
                
                
                
    def LoadBrukerTopspin(self):
        FilePath = askopenfilename()
        if FilePath is not '': #if not canceled
            Dir = os.path.dirname(FilePath) #convert path to file to path of folder
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
            self.current.grid_remove()
            self.current.destroy()
            spec = [False]
            if sizeTD1 is 1:
                #data = np.transpose(ComplexData)[0][:] #convert to 1D np.array
                self.masterData=sc.Spectrum(ComplexData,[freq2],[SW2],spec)
                self.current=sc.Current1D(self,self.masterData) #create the Current1D instance from the Spectrum  
            else:
                data = ComplexData.reshape(sizeTD1,sizeTD2/2)
                self.masterData=sc.Spectrum(data,[freq1,freq2],[SW1,SW2],spec*2)
                self.current=sc.Current1D(self,self.masterData) #create the Current1D instance from the Spectrum  
                
            self.current.grid(row=0,column=0,sticky="nswe")
            self.updAllFrames()        
                
    def LoadChemFile(self):
        FileLocation = askopenfilename()
        Dir = os.path.dirname(FileLocation)
        if FileLocation is not '': #if not empty
           # try:
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
                self.current.grid_remove()
                self.current.destroy()
                spec = [False]
                    
                if sizeTD1 is 1:
                    data = data[0][:] #convert to 1D np.array
                    self.masterData=sc.Spectrum(data,[freq*1e6],[sw],spec)
                    self.current=sc.Current1D(self,self.masterData,0,[],0) #create the Current1D instance from the Spectrum  
                else:
                    #data = np.transpose(data.reshape((sizeTD1,sizeTD2)))
                    data = data.reshape((sizeTD1,sizeTD2))
                    self.masterData=sc.Spectrum(data,[freq*1e6]*2,[sw1,sw],spec*2)
                    self.current=sc.Current1D(self,self.masterData) #create the Current1D instance from the Spectrum  
                    
                self.current.grid(row=0,column=0,sticky="nswe")
                self.updAllFrames()
           # except:
           #     print('Error loading Chemagnetics data from '+Dir+os.path.sep+'data. No data loaded!')
        else:
            print(Dir+os.path.sep+'data does not exits, no Chemagnetics data loaded!')

    def LoadSimpsonFile(self):
        #Loads Simpson data (Fid or Spectrum) to the ssNake data format
        FileLocation = askopenfilename()
        if FileLocation is not '': #if not empty
            try:
                with open(FileLocation, 'r') as f: #read file
                    Lines = f.read().split('\n')
            
                NP, NI, SW, SW1, TYPE, FORMAT = 0,1,0,0,'','Normal'
                for s in range(0,len(Lines)):
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
                    elif Lines[s].startswith('DATA'):
                        DataStart = s
                    elif Lines[s].startswith('END'):
                        DataEnd = s
                
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
                self.current.grid_remove()
                self.current.destroy()
                if 'FID' in TYPE:
                    axis=0
                    spec = [False]
                elif 'SPE' in TYPE:
                    axis=1
                    spec = [True]
                    
                if NI is 1:
                    self.masterData=sc.Spectrum(data,[0],[SW],spec)
                    self.current=sc.Current1D(self,self.masterData,0,[],0) #create the Current1D instance from the Spectrum  
                else:
                    data = np.transpose(data.reshape((NP,NI)))
                    self.masterData=sc.Spectrum(data,[0,0],[SW,SW1],spec*2)
                    self.current=sc.Current1D(self,self.masterData) #create the Current1D instance from the Spectrum  
                    
                self.current.grid(row=0,column=0,sticky="nswe")
                self.updAllFrames()
            except:
                print('Error loading Simpson data from '+FileLocation+' . No data loaded!')

    def saveFigure(self):
        f=asksaveasfilename(filetypes=(('svg','.svg'),('png','.png'),('eps','.eps'),('jpg','.jpg'),('pdf','.pdf')))
        self.current.saveFigure(f)
                
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
                        for iii in range(0,Points[1]):
                            for jjj in range(0,Points[0]):
                                f.write(str(self.masterData.data[jjj][iii].real)+' '+ str(self.masterData.data[jjj][iii].imag)+'\n')
                        #for Line in self.masterData.data:
                        #    for SubLine in Line:
                        #        f.write(str(SubLine.real) +' '+str(SubLine.imag)+'\n')
                                
                    f.write('END')
                    f.close()
        except:
             print('An error occured while saving to Simpson format')
        
    def real(self):
        self.redoList = []
        self.undoList.append(self.masterData.real())
        self.current.upd()
        self.current.showFid()

    def imag(self):
        self.redoList = []
        self.undoList.append(self.masterData.imag())
        self.current.upd()
        self.current.showFid()

    def abs(self):
        self.redoList = []
        self.undoList.append(self.masterData.abs())
        self.current.upd()
        self.current.showFid()

    def fourier(self):
        self.redoList = []
        self.undoList.append(self.current.fourier())
        self.bottomframe.upd()

    def fftshift(self):
        self.redoList = []
        self.undoList.append(self.current.fftshift())
        self.updAllFrames()

    def invFftshift(self):
        self.redoList = []
        self.undoList.append(self.current.fftshift(inv=True))
        self.updAllFrames()
    
    def hilbert(self):
        self.redoList = []
        self.undoList.append(self.current.hilbert())

    def setFreq(self,freq,sw):
        self.redoList = []
        self.undoList.append(self.current.setFreq(freq,sw))

    def createPhaseWindow(self):
        PhaseWindow(self,self.current)
        
    def createApodWindow(self):
        ApodWindow(self,self.current)

    def createSizeWindow(self):
        SizeWindow(self,self.current)

    def createSwapEchoWindow(self):
        SwapEchoWindow(self,self.current)

    def createShiftDataWindow(self):
        ShiftDataWindow(self,self.current)

    def createDCWindow(self):
        DCWindow(self,self.current)

    def createRefWindow(self):
        RefWindow(self,self.current)

    def createIntegrateWindow(self):
        integrateWindow(self,self.current)
        
    def createMaxWindow(self):
        maxWindow(self,self.current)
        
    def createMinWindow(self):
        minWindow(self,self.current)
        
    def createRegionWindow(self):
        extractRegionWindow(self,self.current)

    def flipLR(self):
        self.redoList = []
        self.undoList.append(self.current.flipLR())

    def createShearingWindow(self):
        if self.masterData.dim > 1:
            ShearingWindow(self,self.current)
        else:
            print('Data has too little dimensions for shearing transform')
        
    def BrukerDigital(self):
        pass
#        FilePath = askopenfilename()
#        if FilePath is not '': #if not canceled
#            Dir = os.path.dirname(FilePath) #convert path to file to path of folder
#            if os.path.exists(Dir+os.path.sep+'acqus'):
#                with open(Dir+os.path.sep+'acqus', 'r') as f: 
#                    data = f.read().split('\n')
#                for s in range(0,len(data)): #exctract info from acqus
#                    if data[s].startswith('##$GRPDLY='):
#                        FilterCorrection = float(data[s][10:])
#                
#                self.redoList = []
#                self.undoList.append(self.current.applyPhase(0,-FilterCorrection*2*pi))

    def createRelaxWindow(self):
        root = fit.RelaxWindow(self.parent,self.current)
        root.title("Relaxation Curve") 
        root.grid_rowconfigure(0, weight=1)
        root.grid_columnconfigure(0, weight=1)

    def createPeakDeconvWindow(self):
        root = fit.PeakDeconvWindow(self.parent,self.current)
        root.title("Peak Deconvolution") 
        root.grid_rowconfigure(0, weight=1)
        root.grid_columnconfigure(0, weight=1)
       
    def plot1D(self):
        self.current.grid_remove()
        self.current.destroy()
        self.current = sc.Current1D(self,self.masterData)
        self.current.grid(row=0,column=0,sticky="nswe")
        self.updAllFrames()

    def plotStack(self):
        if len(self.masterData.data.shape) > 1:
            self.current.grid_remove()
            self.current.destroy()
            self.current = sc.CurrentStacked(self,self.masterData) 
            self.current.grid(row=0,column=0,sticky="nswe")
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")

    def plotArray(self):
        if len(self.masterData.data.shape) > 1:
            self.current.grid_remove()
            self.current.destroy()
            self.current = sc.CurrentArrayed(self,self.masterData) 
            self.current.grid(row=0,column=0,sticky="nswe")
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")
            
    def plotContour(self):
        if len(self.masterData.data.shape) > 1:
            self.current.grid_remove()
            self.current.destroy()
            self.current = sc.CurrentContour(self,self.masterData) 
            self.current.grid(row=0,column=0,sticky="nswe")
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")
            
    def plotSkewed(self):
        if len(self.masterData.data.shape) > 1:
            self.current.grid_remove()
            self.current.destroy()
            self.current = sc.CurrentSkewed(self,self.masterData) 
            self.current.grid(row=0,column=0,sticky="nswe")
            self.updAllFrames()
        else:
            print("Data does not have enough dimensions")
    def createXaxWindow(self):
        XaxWindow(self,self.current)

    def updAllFrames(self):
        self.bottomframe.upd()
        self.sideframe.upd()

    def undo(self, *args):
        if self.undoList:
            self.redoList.append(self.undoList.pop()(self.masterData))
            self.current.upd()
            self.current.plotReset()
            self.current.showFid()
            self.updAllFrames()
        else:
            print("no undo information")

    def redo(self, *args):
        if self.redoList:
            self.undoList.append(self.redoList.pop()(self.masterData))
            self.current.upd()
            self.current.plotReset()
            self.current.showFid()
            self.updAllFrames()
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
        Separator(self,orient=HORIZONTAL).grid(row=1, sticky='ew')
        self.frame2 = Frame(self)
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
        self.current = self.parent.current
        self.shape = self.current.data.data.shape
        self.length = len(self.shape)
        self.button1Var.set(self.current.axes)
        offset = 0
        self.plotIs2D = isinstance(self.current, (sc.CurrentStacked,sc.CurrentArrayed,sc.CurrentContour,sc.CurrentSkewed))
        if self.plotIs2D:
            offset = 1
            self.button2Var.set(self.current.axes2)
        for num in self.labels:
            num.destroy()
        self.labels = []
        for num in self.entries:
            num.destroy()
        self.entries=[]
        for num in self.buttons1:
            num.destroy()
        self.buttons1=[]
        for num in self.buttons2:
            num.destroy()
        self.buttons2=[]
        self.entryVars = []
        if self.length > 1:
            for num in range(self.length):
                self.buttons1.append(Radiobutton(self.frame1, variable=self.button1Var, value=num, command=lambda: self.setAxes(True)))
                self.buttons1[num].grid(row=num*2+1,column=0)
                if self.plotIs2D:
                    self.buttons2.append(Radiobutton(self.frame1, variable=self.button2Var, value=num, command=lambda: self.setAxes(False)))
                    self.buttons2[num].grid(row=num*2+1,column=1)
                self.labels.append(Label(self.frame1,text="TD"+str(num+1))) 
                self.labels[num].grid(row=num*2,column=1+offset)
                self.entryVars.append(StringVar())
                if not self.plotIs2D:
                    if num < self.current.axes:
                        self.entryVars[num].set(str(self.current.locList[num]))
                    elif num == self.current.axes:
                        self.entryVars[num].set("0")
                    else:
                        self.entryVars[num].set(str(self.current.locList[num-1]))
                else:
                    if (num < self.current.axes) and (num < self.current.axes2):
                        self.entryVars[num].set(str(self.current.locList[num]))
                    elif (num == self.current.axes) or (num == self.current.axes2):
                        self.entryVars[num].set("0")
                    elif (num > self.current.axes) or (num > self.current.axes2):
                        self.entryVars[num].set(str(self.current.locList[num-1]))
                    else:
                        self.entryVars[num].set(str(self.current.locList[num-2]))
                self.entries.append(Spinbox(self.frame1,textvariable=self.entryVars[num],from_=0,to=self.shape[num]-1,justify="center",command=lambda event=None,num=num: self.getSlice(event,num)))
                self.entries[num].bind("<Return>", lambda event=None,num=num: self.getSlice(event,num)) 
                self.entries[num].bind("<KP_Enter>", lambda event=None,num=num: self.getSlice(event,num)) 
                self.entries[num].grid(row=num*2+1,column=1+offset)
            for child in self.frame2.winfo_children():
                child.destroy()
            if isinstance(self.current, (sc.CurrentStacked,sc.CurrentArrayed,sc.CurrentSkewed)):
                if self.current.stackBegin is not None:
                    self.from2D.set(str(self.current.stackBegin))
                else:
                    self.from2D.set('0')
                if self.current.stackEnd is not None:
                    self.to2D.set(str(self.current.stackEnd))
                else:
                    self.to2D.set(str(self.shape[self.current.axes2]))
                if self.current.stackStep is not None:
                    self.step2D.set(str(self.current.stackStep))
                else:
                    self.step2D.set('1')
                Label(self.frame2,text="From").grid(row=1,column=0,sticky='n')
                self.fromSpin = Spinbox(self.frame2,textvariable=self.from2D,from_=0,to=int(self.to2D.get())-1,justify="center",command=self.setToFrom)
                self.fromSpin.bind("<Return>", self.setToFrom) 
                self.fromSpin.bind("<KP_Enter>", self.setToFrom)
                self.fromSpin.grid(row=2,column=0)
                Label(self.frame2,text="To").grid(row=3,column=0,sticky='n')
                self.toSpin = Spinbox(self.frame2,textvariable=self.to2D,from_=int(self.from2D.get())+1,to=self.shape[self.current.axes2],justify="center",command=self.setToFrom)
                self.toSpin.bind("<Return>", self.setToFrom) 
                self.toSpin.bind("<KP_Enter>", self.setToFrom) 
                self.toSpin.grid(row=4,column=0)
                Label(self.frame2,text="Step").grid(row=5,column=0,sticky='n')
                self.stepSpin = Spinbox(self.frame2,textvariable=self.step2D,from_=1,to=self.shape[self.current.axes2],justify="center",command=self.setToFrom)
                self.stepSpin.bind("<Return>", self.setToFrom) 
                self.stepSpin.bind("<KP_Enter>", self.setToFrom)
                self.stepSpin.grid(row=6,column=0)
                if isinstance(self.current, (sc.CurrentStacked,sc.CurrentArrayed)):
                    self.spacing.set('%.3e' % self.current.spacing)
                    Label(self.frame2,text="Spacing").grid(row=7,column=0,sticky='n')
                    self.spacingEntry = Entry(self.frame2,textvariable=self.spacing,justify="center")
                    self.spacingEntry.bind("<Return>", self.setSpacing) 
                    self.spacingEntry.bind("<KP_Enter>", self.setSpacing) 
                    self.spacingEntry.grid(row=8,column=0)

                if isinstance(self.current, (sc.CurrentSkewed)):
                    self.skew.set('%.2f' % self.current.skewed)
                    Label(self.frame2,text="Skew").grid(row=7,column=0,sticky='n')
                    self.skewEntry = Entry(self.frame2,textvariable=self.skew,justify="center")
                    self.skewEntry.bind("<Return>", self.setSkew) 
                    self.skewEntry.bind("<KP_Enter>", self.setSkew) 
                    self.skewEntry.grid(row=8,column=0)
                    self.elev.set('%.1f' % self.current.elevation)
                    Label(self.frame2,text="Elevation").grid(row=9,column=0,sticky='n')
                    self.elevEntry = Entry(self.frame2,textvariable=self.elev,justify="center")
                    self.elevEntry.bind("<Return>", self.setSkew) 
                    self.elevEntry.bind("<KP_Enter>", self.setSkew) 
                    self.elevEntry.grid(row=10,column=0)
            if isinstance(self.current, (sc.CurrentContour)):
                self.numLevels.set(str(self.current.numLevels))
                self.maxLevels.set(str(self.current.maxLevels*100.0))
                self.minLevels.set(str(self.current.minLevels*100.0))
                Label(self.frame2,text="Number of contours").grid(row=1,column=0,sticky='n')
                self.numLEntry = Entry(self.frame2,textvariable=self.numLevels,justify="center")
                self.numLEntry.bind("<Return>", self.setContour) 
                self.numLEntry.bind("<KP_Enter>", self.setContour) 
                self.numLEntry.grid(row=2,column=0)
                Label(self.frame2,text="Highest contour [%]").grid(row=3,column=0,sticky='n')
                self.maxLEntry = Entry(self.frame2,textvariable=self.maxLevels,justify="center")
                self.maxLEntry.bind("<Return>", self.setContour) 
                self.maxLEntry.bind("<KP_Enter>", self.setContour) 
                self.maxLEntry.grid(row=4,column=0)
                Label(self.frame2,text="Lowest contour [%]").grid(row=5,column=0,sticky='n')
                self.minLEntry = Entry(self.frame2,textvariable=self.minLevels,justify="center")
                self.minLEntry.bind("<Return>", self.setContour) 
                self.minLEntry.bind("<KP_Enter>", self.setContour) 
                self.minLEntry.grid(row=6,column=0)
                    
    def setToFrom(self, *args):
        if isinstance(self.current, (sc.CurrentStacked,sc.CurrentArrayed,sc.CurrentSkewed)):
            fromVar = int(safeEval(self.from2D.get()))
            toVar = int(safeEval(self.to2D.get()))
            stepVar = int(safeEval(self.step2D.get()))
            if fromVar < 0:
                fromVar = 0
            elif fromVar > self.shape[self.current.axes2]-1:
                fromVar = self.shape[self.current.axes2]-1
            if toVar > self.shape[self.current.axes2]:
                toVar = self.shape[self.current.axes2]
            elif toVar <= fromVar:
                toVar = fromVar + 1
            if stepVar < 1:
                stepVar = 1
            elif stepVar > self.shape[self.current.axes2]:
                stepVar = self.shape[self.current.axes2]
            self.current.stackSelect(fromVar,toVar,stepVar)
            self.from2D.set(str(fromVar))
            self.to2D.set(str(toVar))
            self.step2D.set(str(stepVar))
            self.fromSpin.config(to=toVar-1)
            self.toSpin.config(from_=fromVar+1)

    def setSpacing(self, *args):
        var =float(safeEval(self.spacing.get()))
        self.spacing.set('%.3e' % var)
        self.current.setSpacing(var)

    def setSkew(self, *args):
        var =float(safeEval(self.skew.get()))
        self.skew.set('%.2f' % var)
        var2 =float(safeEval(self.elev.get()))
        self.elev.set('%.1f' % var2)
        self.current.setSkewed(var,var2)
        
    def setContour(self, *args):
        var1 = int(round(safeEval(self.numLevels.get())))
        self.numLevels.set(str(var1))
        var2 =float(safeEval(self.maxLevels.get()))
        self.maxLevels.set('%.1f' % var2)
        var3 =float(safeEval(self.minLevels.get()))
        self.minLevels.set('%.1f' % var3)
        self.current.setLevels(var1,var2/100.0,var3/100.0)
        
    def setAxes(self,first=True):
        if self.plotIs2D:
            axes= self.button1Var.get()
            axes2=self.button2Var.get()
            if axes==axes2:
                if first:
                    axes2 = self.current.axes
                else:
                    axes = self.current.axes2
            self.button2Var.set(axes2)
            self.getSlice(None, axes,True)
        else:
            self.getSlice(None, self.button1Var.get(),True)

    def getSlice(self, event, entryNum, button=False): #change the slice which is currently displayed
        if button:
            dimNum = entryNum
        elif not self.plotIs2D:
            if entryNum == self.current.axes:
                if entryNum == self.length-1:
                    dimNum = self.length-2
                else:
                    dimNum = self.length-1
            else:
                dimNum = self.current.axes
        else:
            dimNum = self.current.axes

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
            self.current.setBlock(dimNum,self.button2Var.get(), locList)
        else:
            self.current.setSlice(dimNum,locList)
        self.parent.bottomframe.upd()

################################################################################  
#the bottom frame holding the fourier button and stuff      
class BottomFrame(Frame):
    def __init__(self, parent):
        Frame.__init__(self,parent)
        self.parent = parent
        self.current = parent.current
        self.specVal = IntVar() #value for the time/freq radiobutton
        self.freqVal = StringVar() #value for frequency entybox
        self.swVal = StringVar() #value for sw entrybox
        self.plotOption = StringVar() #value for dropdown plot type box
        self.plotOption.set("Real")
        self.axisOption1 = StringVar()
        self.axisOption2 = StringVar()
        self.echoTick = IntVar()
        self.echoTick.set(0)
        Button(self, text="Fourier",command=self.parent.fourier).grid(row=0,column=0,rowspan=2)
        self.rb1 = Radiobutton(self,text="Time",variable=self.specVal,value=0,command=self.changeSpec)
        self.rb1.grid(row=0,column=1)
        self.rb2 = Radiobutton(self,text="Frequency",variable=self.specVal,value=1,command=self.changeSpec)
        self.rb2.grid(row=1,column=1)
        Checkbutton(self,text="Whole echo",variable=self.echoTick, command=self.setWholeEcho).grid(row=0,column=2,rowspan=2)
        Label(self,text="Freq (MHz)").grid(row=0,column=3)
        self.freqEntry = Entry(self,textvariable=self.freqVal,justify='center')
        self.freqEntry.bind("<Return>", self.changeFreq)
        self.freqEntry.bind("<KP_Enter>", self.changeFreq)
        self.freqEntry.grid(row=1,column=3)
        Label(self,text="Sweepwidth (kHz)").grid(row=0,column=4)
        self.swEntry = Entry(self,textvariable=self.swVal,justify='center')
        self.swEntry.bind("<Return>", self.changeFreq)
        self.swEntry.bind("<KP_Enter>", self.changeFreq)
        self.swEntry.grid(row=1,column=4)
        Label(self,text="Plot").grid(row=0,column=5)
        self.plotDrop = OptionMenu(self, self.plotOption,"Real","Real", "Imag", "Both","Abs",command=self.changePlot)
        self.plotDrop.grid(row=1,column=5)
        Label(self,text="Axis").grid(row=0,column=6)
        self.axisDropTime = OptionMenu(self, self.axisOption1, "s", "s", "ms", u"\u03bcs",command=self.changeAxis)
        self.axisDropFreq = OptionMenu(self, self.axisOption2, "Hz", "Hz", "kHz", "MHz","ppm",command=self.changeAxis)
        self.axisDropTime.grid(row=1,column=6)
        self.axisDropFreq.grid(row=1,column=6)
        self.swEntry
        self.upd()

    def frameEnable(self):
        for child in self.winfo_children():
            child.configure(state='normal')
            
    def frameDisable(self):
        for child in self.winfo_children():
            child.configure(state='disabled')
        
    def upd(self): #upd the values displayed in the bottom menu
        self.current = self.parent.current
        self.freqVal.set(str(self.current.freq/1000000)) #show in MHz
        self.swVal.set(str(self.current.sw/1000)) #show in kHz
        if self.current.spec==0:
            self.specVal.set(0)
            self.axisDropFreq.grid_forget()
            self.axisDropTime.grid(row=1,column=6)
            val = self.current.axType
            if val == 0:
                self.axisOption1.set("s")
            elif val == 1:
                self.axisOption1.set("ms")
            elif val == 2:
                self.axisOption1.set( u"\u03bcs")
                
        elif self.current.spec==1:
            self.specVal.set(1)
            self.axisDropTime.grid_forget()
            self.axisDropFreq.grid(row=1,column=6)
            val = self.current.axType
            if val == 0:
                self.axisOption1.set("Hz")
            elif val == 1:
                self.axisOption1.set("kHz")
            elif val == 2:
                self.axisOption1.set("MHz")
            elif val == 3:
                self.axisOption1.set("ppm")
        if self.current.wholeEcho:
            self.echoTick.set(1)
        else:
            self.echoTick.set(0)

    def setWholeEcho(self):
        self.current.setWholeEcho(self.echoTick.get())

    def changeSpec(self, *args): #change from time to spectral domain and vice versa
        self.parent.redoList = []
        self.parent.undoList.append(self.current.changeSpec(self.specVal.get()))
        self.upd()

    def changeFreq(self, *args): #change the frequency and sw of the displayed axes
        freq = safeEval(self.freqVal.get())*1000000 #show in MHz
        sw = safeEval(self.swVal.get())*1000 #show in kHz
        if freq != 0 and sw != 0:
            self.parent.setFreq(freq,sw)
        self.upd()
    
    def changePlot(self, *args): #change the plot type
        pType = self.plotOption.get()
        if pType == "Real":
            self.current.plotType=0
        elif pType == "Imag":
            self.current.plotType=1
        elif pType == "Both":
            self.current.plotType=2
        elif pType == "Abs":
            self.current.plotType=3
        self.current.showFid()

    def changeAxis(self, *args):
        if self.current.spec == 0:
            pType = self.axisOption1.get()
            if pType == "s":
                self.current.setAxType(0)
            elif pType == "ms":
                self.current.setAxType(1)
            elif pType == u"\u03bcs":
                self.current.setAxType(2)
        if self.current.spec == 1:
            pType = self.axisOption2.get()
            if pType == "Hz":
                self.current.setAxType(0)
            elif pType == "kHz":
                self.current.setAxType(1)
            elif pType == "MHz":
                self.current.setAxType(2)
            elif pType == "ppm":
                self.current.setAxType('ppm')
        self.current.showFid()

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
        Entry(self,textvariable=self.pos,justify='center').grid(row=0,column=2)
        Label(self,text="x-value:").grid(row=0,column=3)
        Entry(self,textvariable=self.xpoint,justify='center').grid(row=0,column=4)
        Label(self,text="y-value:").grid(row=0,column=5)
        Entry(self,textvariable=self.ypoint,justify='center').grid(row=0,column=6)
        Label(self,text=u"\u0394x:").grid(row=0,column=7)
        Entry(self,textvariable=self.deltaxpoint,justify='center').grid(row=0,column=8)
        Label(self,text=u"\u0394y:").grid(row=0,column=9)
        Entry(self,textvariable=self.deltaypoint,justify='center').grid(row=0,column=10)

    def frameEnable(self):
        for child in self.winfo_children():
            child.configure(state='normal')
            
    def frameDisable(self):
        for child in self.winfo_children():
            child.configure(state='disabled')
        
    def setLabels(self,position):
        self.deltaxpoint.set('%.3e' % np.abs(self.oldx-position[1]))
        self.deltaypoint.set('%.3e' % np.abs(self.oldy-position[2]))
        self.pos.set(str(position[0]))
        self.xpoint.set('%.3e' % position[1])
        self.ypoint.set('%.3e' % position[2])
        self.oldx = position[1]
        self.oldy = position[2]

    def getPosition(self, *args):
        self.parent.current.peakPickFunc = lambda pos,self=self: self.setLabels(pos) 
        self.parent.current.peakPick = True
        
#################################################################################   
class PhaseWindow(Frame): #a window for phasing the data
    def __init__(self, parent,current):
        Frame.__init__(self, parent)
        parent.menuDisable()
        #initialize variables for the widgets
        self.zeroVal = 0.0
        self.firstVal = 0.0
        self.refVal = 0.0
        self.zeroValue = StringVar()
        self.zeroValue.set("0.00")
        self.firstValue = StringVar()
        self.firstValue.set("0.00")
        self.refValue = StringVar()
        self.refValue.set("0.0")
        #set stepsizes for the buttons
        self.phase0step = 1.0
        self.phase1step = 0.01
        #create a new window
        self.parent = parent
        self.current = current
        self.window = Toplevel(self)
        self.window.geometry('+0+0')
        self.window.transient(self.parent)
        self.window.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.window.title("Phasing")
        self.window.resizable(width=FALSE, height=FALSE)
        Label(self.window,text="Zero order phasing").grid(row=0,column=0,columnspan=3)
        Button(self.window,text="Autophase 0th order",command=lambda: self.autophase(0)).grid(row=1,column=1)
        self.zeroEntry = Entry(self.window,textvariable=self.zeroValue,justify="center")
        self.zeroEntry.bind("<Return>", self.inputZeroOrder)
        self.zeroEntry.bind("<KP_Enter>", self.inputZeroOrder)
        self.zeroEntry.grid(row=2,column=1)
        tk.Button(self.window,text="<",repeatdelay=100, repeatinterval=1,command=lambda:self.stepPhase(-1,0)).grid(row=2,column=0)
        tk.Button(self.window,text=">",repeatdelay=100, repeatinterval=1,command=lambda:self.stepPhase(1,0)).grid(row=2,column=2)
        self.zeroScale=Scale(self.window, from_=-180, to=180,  orient="horizontal", command=self.setZeroOrder,length=300)
        self.zeroScale.grid(row=3,column=0,columnspan=3)
        Label(self.window,text="First order phasing").grid(row=4,column=0,columnspan=3)
        Button(self.window,text="Autophase 0th+1st order",command=lambda: self.autophase(1)).grid(row=5,column=1)
        self.firstEntry = Entry(self.window,textvariable=self.firstValue,justify="center")
        self.firstEntry.bind("<Return>", self.inputFirstOrder) 
        self.firstEntry.bind("<KP_Enter>", self.inputFirstOrder) 
        self.firstEntry.grid(row=6,column=1)
        tk.Button(self.window,text="<",repeatdelay=100, repeatinterval=1,command=lambda:self.stepPhase(0,-1)).grid(row=6,column=0)
        tk.Button(self.window,text=">",repeatdelay=100, repeatinterval=1,command=lambda:self.stepPhase(0,1)).grid(row=6,column=2)
        self.firstScale=Scale(self.window, from_=-0.01*180*(self.current.data1D.shape[-1])/self.current.sw, to=0.01*180*self.current.data1D.shape[-1]/self.current.sw, orient="horizontal", command=self.setFirstOrder,length=300)
        self.firstScale.grid(row=7,column=0,columnspan=3)
        if self.current.spec > 0:
            Label(self.window,text="Reference").grid(row=8,column=0,columnspan=3)
            self.refEntry = Entry(self.window,textvariable=self.refValue,justify="center")
            self.refEntry.bind("<Return>", self.inputRef) 
            self.refEntry.bind("<KP_Enter>", self.inputRef)
            self.refEntry.grid(row=9,column=1)
            Button(self.window, text="Pick reference", command=self.pickRef).grid(row=10,column=1)
        Button(self.window, text="Apply",command=self.applyPhaseAndClose).grid(row=11,column=0)
        Button(self.window, text="Cancel",command=self.cancelAndClose).grid(row=11,column=2)      
        
    def setZeroOrder(self,value, *args): #function called by the zero order scale widget
        self.zeroVal = float(value)
        self.zeroValue.set('%.2f' % self.zeroVal)
        self.current.setPhaseInter(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0)
        
    def inputZeroOrder(self, *args): #function called by the zero order entry widget
        inp = safeEval(self.zeroValue.get())
        self.zeroVal = np.mod(inp+180,360)-180
        self.zeroScale.set(self.zeroVal) #setting the scale to a value calls the previous function, so the phase of current doesn't need to be set here

    def setFirstOrder(self,value, *args): #function called by the first order scale widget
        newZero = (self.zeroVal-(float(value)-self.firstVal)*self.refVal) #calculate the new zero order phase depending on the reference
        self.zeroVal = np.mod(newZero+180,360)-180
        self.firstVal = float(value)
        self.firstValue.set('%.2f' % self.firstVal)
        self.zeroValue.set('%.2f' % self.zeroVal)
        self.zeroScale.set(self.zeroVal)

    def inputFirstOrder(self, *args): #function called by the first order entry widget
        self.firstVal = safeEval(self.firstValue.get())
        self.firstScale.set(self.firstVal) #setting the scale to a value calls the previous function, so the phase of current doesn't need to be set here

    def autophase(self, num): #run the autophase for either zero order (0) or both orders (1)
        phases = self.current.autoPhase(num)
        if num == 0:
            self.zeroVal=(np.mod(phases[0]+180,360)-180)
            self.zeroValue.set('%.2f' % self.zeroVal)
            self.zeroScale.set(self.zeroVal)
        elif num == 1:
            self.zeroVal=(np.mod(phases[0]+180,360)-180)
            self.zeroValue.set('%.2f' % self.zeroVal)
            self.zeroScale.set(self.zeroVal)
            self.firstVal = phases[1]
            self.firstValue.set('%.2f' % self.firstVal)
            self.firstScale.set(self.firstVal)

    def stepPhase(self,phase0,phase1): #step phase from < and > keys
        inp = safeEval(self.zeroValue.get())+phase0*self.phase0step
        self.zeroVal = np.mod(inp+180,360)-180
        self.zeroScale.set(self.zeroVal)
        self.firstVal = safeEval(self.firstValue.get())+phase1*self.phase1step
        self.firstScale.set(self.firstVal)

    def inputRef(self, *args): #set the reference from the entry widget
        self.refVal = safeEval(self.refValue.get())
        self.refValue.set('%.2f' % self.refVal)

    def setRef(self,value,*args):
        self.refVal = float(value)
        self.refValue.set('%.2f' % self.refVal)

    def pickRef(self, *args): #run the pick function to pick the reference value
        self.current.peakPickFunc = lambda pos,self=self: self.setRef(self.current.xax[pos[0]])
        self.current.peakPick = True

    def cancelAndClose(self):
        self.current.upd()
        self.current.showFid()
        self.parent.menuEnable()
        self.window.destroy()

    def applyPhaseAndClose(self):
        self.parent.redoList = []
        self.parent.undoList.append(self.current.applyPhase(np.pi*self.zeroVal/180.0,np.pi*self.firstVal/180.0))
        self.parent.menuEnable()
        self.window.destroy()

################################################################
class ApodWindow(Frame): #a window for apodization
    def __init__(self, parent,current):
        parent.menuDisable()
        Frame.__init__(self, parent)
        self.parent = parent
        self.current = current
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
        if self.current.data.dim > 1:
            options = map(str,np.delete(range(self.current.data.dim),self.current.axes))
            self.shiftingAxes = StringVar()
            self.shiftingAxes.set(options[0])
        #set stepsizes for the buttons
        self.lorstep = 1.0
        self.gaussstep = 1.0
        self.window = Toplevel(self)
        self.window.geometry('+0+0')
        self.window.transient(self.parent)
        self.window.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.window.title("Apodize")
        self.window.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self.window)
        self.frame1.grid(row=0,column=0)
        Label(self.frame1,text="Lorentzian").grid(row=0,column=0,columnspan=4)
        Checkbutton(self.frame1,variable=self.lorTick, command=lambda: self.checkEval(self.lorTick,self.lorEntry)).grid(row=1,column=1)
        self.lorEntry = Entry(self.frame1,textvariable=self.lorVal,justify="center", state='disabled')   
        tk.Button(self.frame1,text="<",repeatdelay=100, repeatinterval=1,command=lambda:self.stepLB(-0.5*self.current.sw/(self.current.data1D.shape[-1]),0)).grid(row=1,column=0)
        tk.Button(self.frame1,text=">",repeatdelay=100, repeatinterval=1,command=lambda:self.stepLB(0.5*self.current.sw/(self.current.data1D.shape[-1]),0)).grid(row=1,column=4)
        self.lorEntry.bind("<Return>", self.apodPreview)
        self.lorEntry.bind("<KP_Enter>", self.apodPreview)
        self.lorEntry.grid(row=1,column=2)
        self.lorScale=Scale(self.frame1, from_=0, to=100.0*self.current.sw/(self.current.data1D.shape[-1]),  orient="horizontal", command=self.setLor,length=200)
        self.lorScale.grid(row=2,column=1,columnspan=2)
        Label(self.frame1,text="Gaussian").grid(row=3,column=0,columnspan=4)
        Checkbutton(self.frame1,variable=self.gaussTick, command=lambda: self.checkEval(self.gaussTick,self.gaussEntry)).grid(row=4,column=1)
        self.gaussEntry = Entry(self.frame1,textvariable=self.gaussVal,justify="center", state='disabled')
        self.gaussEntry.bind("<Return>", self.apodPreview)
        self.gaussEntry.bind("<KP_Enter>", self.apodPreview)
        self.gaussEntry.grid(row=4,column=2)
        tk.Button(self.frame1,text="<",repeatdelay=100, repeatinterval=1,command=lambda:self.stepLB(0,-0.5*self.current.sw/(self.current.data1D.shape[-1]))).grid(row=4,column=0)
        tk.Button(self.frame1,text=">",repeatdelay=100, repeatinterval=1,command=lambda:self.stepLB(0,0.5*self.current.sw/(self.current.data1D.shape[-1]))).grid(row=4,column=4)
        self.gaussScale=Scale(self.frame1, from_=0, to=100.0*self.current.sw/(self.current.data1D.shape[-1]),  orient="horizontal", command=self.setGauss,length=200)
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
        if self.current.data.dim > 1:
            Label(self.frame1,text="Shifting").grid(row=12,column=0,columnspan=4)
            self.shiftingEntry = Entry(self.frame1,textvariable=self.shiftingVal,justify="center")
            self.shiftingEntry.bind("<Return>", self.apodPreview)
            self.shiftingEntry.bind("<KP_Enter>", self.apodPreview)
            self.shiftingEntry.grid(row=13,column=2)
            OptionMenu(self.frame1,self.shiftingAxes, self.shiftingAxes.get(),*options).grid(row=14,column=2)
        self.frame2 = Frame(self.window)
        self.frame2.grid(row=1,column=0)
        Button(self.frame2, text="Apply",command=self.applyApodAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=2) 

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
        if self.current.data.dim > 1:
            shifting = safeEval(self.shiftingVal.get())
            self.shiftingVal.set(shifting)
            shiftingAxes = int(self.shiftingAxes.get())
        else:
            shiftingAxes = None
        self.current.apodPreview(lor,gauss,cos2,hamming,shift,shifting,shiftingAxes)

    def stepLB(self,lorincr,gaussincr): #step linebroadening from < and > keys
        if lorincr!=0:
            self.lorScale.set(float(self.lorVal.get())+lorincr*self.lorstep)
        if gaussincr!=0:
            self.gaussScale.set(float(self.gaussVal.get())+gaussincr*self.gaussstep)

    def cancelAndClose(self):
        self.current.upd()
        self.current.showFid()
        self.parent.menuEnable()
        self.window.destroy()

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
        if self.current.data.dim > 1:
            shifting = safeEval(self.shiftingVal.get())
            self.shiftingVal.set(shifting)
            shiftingAxes = int(self.shiftingAxes.get())
        else:
            shiftingAxes = None
        self.parent.redoList = []
        self.parent.undoList.append(self.current.applyApod(lor,gauss,cos2,hamming,shift,shifting,shiftingAxes))
        self.parent.menuEnable()
        self.window.destroy()

#######################################################################################
class SizeWindow(Frame): #a window for changing the size of the current dimension
    def __init__(self, parent,current):
        parent.menuDisable()
        #initialize variables for the widgets
        self.sizeVal = StringVar()
        self.sizeVal.set(str(current.data1D.shape[-1]))
        #create a new window
        Frame.__init__(self, parent)
        self.parent = parent
        self.current = current
        self.window = Toplevel(self)
        self.window.geometry('+0+0')
        self.window.transient(self.parent)
        self.window.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.window.title("Set size")
        self.window.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self.window)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Set size").grid(row=0,column=0,columnspan=2)
        self.sizeEntry = Entry(self.frame1,textvariable=self.sizeVal,justify="center")
        self.sizeEntry.bind("<Return>", self.sizePreview)
        self.sizeEntry.bind("<KP_Enter>", self.sizePreview)
        self.sizeEntry.grid(row=1,column=0,columnspan=2)
        self.frame2 = Frame(self.window)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applySizeAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)
 
    def sizePreview(self, *args): #display the size preview from the entry widget value
        size = int(round(safeEval(self.sizeVal.get())))
        if size < 1:
            size = 1
        self.sizeVal.set(str(size))
        self.current.setSizePreview(size)

    def cancelAndClose(self):
        self.current.upd()
        self.current.plotReset()
        self.current.showFid()
        self.parent.menuEnable()
        self.window.destroy()

    def applySizeAndClose(self):
        size = int(round(safeEval(self.sizeVal.get())))
        if size < 1:
            size = 1
        self.parent.redoList = []
        self.parent.undoList.append(self.current.applySize(size))
        self.parent.sideframe.upd()
        self.parent.menuEnable()
        self.window.destroy()

##########################################################################################
class SwapEchoWindow(Frame): #a window for changing the size of the current dimension
    def __init__(self, parent,current):
        parent.menuDisable()
        Frame.__init__(self, parent)
        #initialize variables for the widgets
        self.posVal = StringVar()
        self.posVal.set(str(int(round(0.5*len(current.data1D)))))
        self.parent = parent
        self.current = current
        self.window = Toplevel(self)
        self.window.geometry('+0+0')
        self.window.transient(self.parent)
        self.window.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.window.title("Swap echo")
        self.window.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self.window)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Position").grid(row=0,column=0,columnspan=2)
        self.posEntry = Entry(self.frame1,textvariable=self.posVal,justify="center")
        self.posEntry.bind("<Return>", self.swapEchoPreview)
        self.posEntry.bind("<KP_Enter>", self.swapEchoPreview)
        self.posEntry.grid(row=1,column=0,columnspan=2)
        self.frame2 = Frame(self.window)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applySwapEchoAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)
        #activate the peak picking
        self.current.peakPickFunc = lambda pos,self=self: self.pickedAndClose(pos) 
        self.current.peakPick = True
 
    def swapEchoPreview(self, *args): #preview the swap echo result from the entry widget
        pos = int(round(safeEval(self.posVal.get())))
        if pos > 0 and pos < (self.current.data1D.shape[-1]):
            self.current.setSwapEchoPreview(pos)

    def cancelAndClose(self):
        self.current.peakPickReset()
        self.current.upd()
        self.current.plotReset()
        self.current.showFid()
        self.parent.menuEnable()
        self.window.destroy()

    def applySwapEchoAndClose(self):
        self.current.peakPickReset()
        pos = int(round(safeEval(self.posVal.get())))
        if pos > 0 and pos < (self.current.data1D.shape[-1]):
            self.parent.redoList = []
            self.parent.undoList.append(self.current.applySwapEcho(pos))
            self.parent.bottomframe.upd()
            self.parent.menuEnable()
            self.window.destroy()
        else:
            print("not a valid index for swap echo")
        
    def pickedAndClose(self,pos): #apply directly if picked since another doesn't make pick doesn't make sense. find a good way to do both entry and picking in a proper way
        self.parent.redoList = []
        self.parent.undoList.append(self.current.applySwapEcho(pos[0]))
        self.parent.bottomframe.upd()
        self.parent.menuEnable()
        self.window.destroy()

###########################################################################
class ShiftDataWindow(Frame): #a window for shifting the data
    def __init__(self, parent,current):
        parent.menuDisable()
        Frame.__init__(self, parent)
        #initialize variables for the widgets
        self.shiftVal = StringVar()
        self.shiftVal.set("0")
        self.parent = parent
        self.current = current
        self.window = Toplevel(self)
        self.window.geometry('+0+0')
        self.window.transient(self.parent)
        self.window.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.window.title("Shift data")
        self.window.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self.window)
        self.frame1.grid(row=0)
        Label(self.frame1,text="number of points to shift").grid(row=0,column=1)
        tk.Button(self.frame1,text="<",repeatdelay=100, repeatinterval=1,command=self.stepDownShift).grid(row=1,column=0)
        tk.Button(self.frame1,text=">",repeatdelay=100, repeatinterval=1,command=self.stepUpShift).grid(row=1,column=2)
        self.posEntry = Entry(self.frame1,textvariable=self.shiftVal,justify="center")
        self.posEntry.bind("<Return>", self.shiftPreview)
        self.posEntry.bind("<KP_Enter>", self.shiftPreview)
        self.posEntry.grid(row=1,column=1)
        self.frame2 = Frame(self.window)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyShiftAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)

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
        self.current.setShiftPreview(shift)

    def cancelAndClose(self):
        self.current.upd()
        self.current.plotReset()
        self.current.showFid()
        self.parent.menuEnable()
        self.window.destroy()

    def applyShiftAndClose(self):
        shift = int(round(safeEval(self.shiftVal.get())))
        self.parent.redoList = []
        self.parent.undoList.append(self.current.applyShift(shift))
        self.parent.menuEnable()
        self.window.destroy()

#############################################################
class DCWindow(Frame): #a window for shifting the data
    def __init__(self, parent,current):
        parent.menuDisable()
        Frame.__init__(self, parent)
        #initialize variables for the widgets
        self.minVal = StringVar()
        self.minVal.set("0")
        self.maxVal = StringVar()
        self.maxVal.set(str(current.data1D.shape[-1]))
        self.parent = parent
        self.current = current
        self.window = Toplevel(self)
        self.window.geometry('+0+0')
        self.window.transient(self.parent)
        self.window.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.window.title("DC offset correction")
        self.window.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self.window)
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
        self.frame2 = Frame(self.window)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyDCAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)
        #pick function
        self.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
        self.current.peakPick = True

    def picked(self,pos,second=False): #pick a value alternating the first and second value determined by the second value.
        if second:
            dataLength = self.current.data1D.shape[-1]
            minimum=int(round(safeEval(self.minVal.get())))
            if minimum < 0:
                minimum = 0
            elif minimum > dataLength:
                minimum = dataLength
            self.minVal.set(str(minimum))
            maximum=pos[0]
            self.maxVal.set(str(maximum))
            self.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
            self.current.peakPick = True
            self.current.dcOffset(minimum,maximum)
        else:
            self.minVal.set(str(pos[0]))
            self.current.peakPickFunc = lambda pos,self=self: self.picked(pos,True) 
            self.current.peakPick = True

    def dcPreview(self, *args): #preview the dc offset correction
        dataLength = self.current.data1D.shape[-1]
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
        self.current.dcOffset(minimum,maximum)

    def cancelAndClose(self):
        self.current.peakPickReset()
        self.current.upd()
        self.current.plotReset()
        self.current.showFid()
        self.parent.menuEnable()
        self.window.destroy()

    def applyDCAndClose(self):
        self.current.peakPickReset()
        dataLength = self.current.data1D.shape[-1]
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
        self.parent.redoList = []
        self.parent.undoList.append(self.current.applydcOffset(minimum,maximum))
        self.parent.menuEnable()
        self.window.destroy()

#############################################################
class regionWindow(Toplevel): #A general region selection frame
    def __init__(self, parent,current,name):
        parent.menuDisable()
        Toplevel.__init__(self)
        #initialize variables for the widgets
        self.minVal = StringVar()
        self.minVal.set("0")
        self.maxVal = StringVar()
        self.maxVal.set(str(current.data1D.shape[-1]))
        self.parent = parent
        self.current = current
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
        self.current.peakPickFunc = lambda pos,self=self: self.picked(pos)
        self.current.peakPick = True

    def picked(self,pos,second=False): #pick a value alternating the first and second value determined by the second value.
        if second:
            dataLength = self.current.data1D.shape[-1]
            minimum=int(round(safeEval(self.minVal.get())))
            if minimum < 0:
                minimum = 0
            elif minimum > dataLength:
                minimum = dataLength
            self.minVal.set(str(minimum))
            maximum=pos[0]
            self.maxVal.set(str(maximum))
            self.current.peakPickFunc = lambda pos,self=self: self.picked(pos) 
            self.current.peakPick = True
        else:
            self.minVal.set(str(pos[0]))
            self.current.peakPickFunc = lambda pos,self=self: self.picked(pos,True) 
            self.current.peakPick = True

    def checkValues(self, *args): #not really preview but just to check the inputs
        dataLength = self.current.data1D.shape[-1]
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
        self.current.peakPickReset()
        self.parent.updAllFrames()
        self.parent.menuEnable()
        self.destroy()

    def apply(self,maximum,minimum):
        pass
        
    def applyAndClose(self):
        self.current.peakPickReset()
        dataLength = self.current.data1D.shape[-1]
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
        self.parent.menuEnable()
        self.apply(maximum,minimum)
        self.destroy()
        
############################################################
class integrateWindow(regionWindow): #A window for obtaining the integral of a selected region
    def __init__(self, parent,current):
        regionWindow.__init__(self,parent,current,'Integrate')

    def apply(self,maximum,minimum):
        self.parent.redoList = []
        self.parent.undoList.append(self.current.integrate(minimum,maximum))
        self.parent.current.grid_remove()
        self.parent.current.destroy()
        self.parent.current=sc.Current1D(self.parent,self.parent.masterData)
        self.parent.current.grid(row=0,column=0,sticky='nswe')
        self.parent.updAllFrames()
        
############################################################
class maxWindow(regionWindow): #A window for obtaining the max of a selected region
    def __init__(self, parent,current):
        regionWindow.__init__(self,parent,current,'Max')

    def apply(self,maximum,minimum):
        self.parent.redoList = []
        self.parent.undoList.append(self.current.maxMatrix(minimum,maximum))
        #self.parent.undoList.append(self.current.maxMatrix(minimum,maximum))
        self.parent.current.grid_remove()
        self.parent.current.destroy()
        self.parent.current=sc.Current1D(self.parent,self.parent.masterData)
        self.parent.current.grid(row=0,column=0,sticky='nswe')
        self.parent.updAllFrames()

############################################################
class minWindow(regionWindow): #A window for obtaining the min of a selected region
    def __init__(self, parent,current):
        regionWindow.__init__(self,parent,current,'Min')

    def apply(self,maximum,minimum):
        self.parent.redoList = []
        self.parent.undoList.append(self.current.minMatrix(minimum,maximum))
        self.parent.current.grid_remove()
        self.parent.current.destroy()
        self.parent.current=sc.Current1D(self.parent,self.parent.masterData)
        self.parent.current.grid(row=0,column=0,sticky='nswe')
        self.parent.updAllFrames()
        
############################################################
class extractRegionWindow(regionWindow): #A window for obtaining a selected region
    def __init__(self, parent,current):
        regionWindow.__init__(self,parent,current,'Extract part')

    def apply(self,maximum,minimum):
        self.parent.redoList = []
        self.parent.undoList.append(self.current.getRegion(minimum,maximum))
        self.parent.updAllFrames()

################################################################
class ShearingWindow(Frame): #a window for setting the xax of the current data
    def __init__(self, parent,current):
        parent.menuDisable()
        Frame.__init__(self, parent)
        self.parent = parent
        self.current = current
        #initialize variables for the widgets
        self.shear = StringVar()
        self.shear.set('0.0')
        options = map(str,range(self.current.data.dim))
        self.axes = StringVar()
        self.axes.set('0')
        self.axes2 = StringVar()
        self.axes2.set('1')
        self.parent = parent
        self.current = current
        self.window = Toplevel(self)
        self.window.geometry('+0+0')
        self.window.transient(self.parent)
        self.window.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.window.title("Shearing")
        self.window.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self.window)
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
        self.frame2 = Frame(self.window)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)

    def shearPreview(self, *args):
        shear = float(safeEval(self.shear.get()))
        self.shear.set(str(shear))

    def cancelAndClose(self):
        #self.current.upd()
        #self.current.plotReset()
        #self.current.showFid()
        self.parent.menuEnable()
        self.window.destroy()

    def applyAndClose(self):
        shear = float(safeEval(self.shear.get()))
        axes = int(self.axes.get())
        axes2 = int(self.axes2.get())
        if axes == axes2:
            print("Axes can't be the same for shearing")
        else:
            self.parent.redoList = []
            self.parent.undoList.append(self.current.shearing(shear,axes,axes2))
            self.parent.menuEnable()
            self.window.destroy()
        
##########################################################################################
class XaxWindow(Frame): #a window for setting the xax of the current data
    def __init__(self, parent,current):
        parent.menuDisable()
        Frame.__init__(self, parent)
        #initialize variables for the widgets
        self.val = StringVar()
        self.parent = parent
        self.current = current
        self.window = Toplevel(self)
        self.window.geometry('+0+0')
        self.window.transient(self.parent)
        self.window.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.window.title("User defined x-axis")
        self.window.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self.window)
        self.frame1.grid(row=0)
        Label(self.frame1,text="Expression for x-axis values").grid(row=0,column=0,columnspan=2)
        self.minEntry = Entry(self.frame1,textvariable=self.val,justify="center")
        self.minEntry.bind("<Return>", self.xaxPreview)
        self.minEntry.bind("<KP_Enter>", self.xaxPreview)
        self.minEntry.grid(row=1,column=0,columnspan=2)

        self.frame2 = Frame(self.window)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyXaxAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)

    def xaxPreview(self, *args):
        env = vars(np).copy()
        env['length']=int(self.current.data1D.shape[-1]) # so length can be used to in equations
        val=eval(self.val.get(),env)                # find a better solution, also add catch for exceptions          
        if isinstance(val,(list,np.ndarray)):
            if len(val)==self.current.data1D.shape[-1]:
                if all(isinstance(x,(int,float)) for x in val):
                    self.current.setXaxPreview(val)
                else:
                    print("Array is not all of int or float type")
            else:
                print("Length of input does not match length of data")
        else:
            print("Input is not a list or array")

    def cancelAndClose(self):
        self.current.upd()
        self.current.plotReset()
        self.current.showFid()
        self.parent.menuEnable()
        self.window.destroy()

    def applyXaxAndClose(self):
        env = vars(np).copy()
        env['length']=int(self.current.data1D.shape[-1]) # so length can be used to in equations
        val=eval(self.val.get(),env)                # find a better solution, also add catch for exceptions
        if isinstance(val,(list,np.ndarray)):
            if len(val)==self.current.data1D.shape[-1]:
                if all(isinstance(x,(int,float)) for x in val):
                    self.current.setXax(val)
                    self.parent.menuEnable()
                    self.window.destroy()
                else:
                    print("Array is not all of int or float type")
            else:
                print("Length of input does not match length of data")
        else:
            print("Input is not a list or array")

##########################################################################################
class RefWindow(Frame): #a window for setting the ppm reference
    def __init__(self, parent,current):
        parent.menuDisable()
        if current.spec == 0:
            print('Setting ppm is only available for frequency data')
        Frame.__init__(self, parent)
        #initialize variables for the widgets
        self.freqVal = StringVar()
        self.freqVal.set(str(current.freq))
        self.refVal = StringVar()
        self.refVal.set('0.0')
        self.parent = parent
        self.current = current
        self.window = Toplevel(self)
        self.window.geometry('+0+0')
        self.window.transient(self.parent)
        self.window.protocol("WM_DELETE_WINDOW", self.cancelAndClose)
        self.window.title("Set ppm reference")
        self.window.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self.window)
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
        self.frame2 = Frame(self.window)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Apply",command=self.applyAndClose).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.cancelAndClose).grid(row=0,column=1)
        #activate the peak picking
        self.current.peakPickFunc = lambda pos,self=self: self.picked(pos)
        self.current.peakPick = True
 
    def preview(self, *args): #fix the input values
        freq = safeEval(self.freqVal.get())
        ref = safeEval(self.refVal.get())
        self.freqVal.set(str(freq))
        self.refVal.set(str(ref))

    def cancelAndClose(self):
        self.current.peakPickReset()
        self.current.showFid()
        self.parent.menuEnable()
        self.window.destroy()

    def applyAndClose(self):
        self.current.peakPickReset()
        freq = safeEval(self.freqVal.get())
        ref = safeEval(self.refVal.get())
        self.parent.redoList = []
        self.parent.undoList.append(self.current.setRef(freq/(1.0+ref*1e-6)))
        self.parent.menuEnable()
        self.window.destroy()
        
    def picked(self,pos): 
        self.freqVal.set(str(self.current.freq+self.current.xax[pos[0]]))
        self.current.peakPickFunc = lambda pos,self=self: self.picked(pos)
        self.current.peakPick = True
            
#################################################################################    
#the main program
if __name__ == "__main__":
    root = Tk()
    mainWindow = Main1DWindow(root) #create an instance to control the main window
    mainWindow.pack(fill=BOTH,expand=1)
    mainWindow.rowconfigure(0, weight=1)
    mainWindow.grid_columnconfigure(0, weight=1)
    root.title("ssNake") 
    root.style = Style()
    root.style.theme_use("clam")
    root.attributes('-zoomed', True)
    root.mainloop()
