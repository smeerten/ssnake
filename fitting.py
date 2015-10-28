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
    from tkinter.filedialog import asksaveasfilename
else:
    from Tkinter import *
    import Tkinter as tk
    from ttk import *
    from tkFileDialog   import asksaveasfilename
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import weakref
import scipy.optimize
import scipy.ndimage
import math
import os.path
import copy
from safeEval import safeEval
from spectrumFrame import Plot1DFrame
from zcw import *


pi = math.pi

##############################################################################
class RelaxWindow(Frame): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        self.fig = self.mainProgram.getFig()
        self.canvas = FigureCanvasTkAgg(self.fig, master=weakref.proxy(self))
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="nswe")
        self.current = RelaxFrame(self,self.fig,self.canvas,oldMainWindow.current)
        self.paramframe = RelaxParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)

    def createNewData(self,data, axes):
        masterData = self.get_masterData()
        self.mainProgram.dataFromFit(data, copy.deepcopy(masterData.freq), copy.deepcopy(masterData.sw) , copy.deepcopy(masterData.spec), copy.deepcopy(masterData.wholeEcho), copy.deepcopy(masterData.ref), copy.deepcopy(masterData.xaxArray), axes)
        
    def rename(self,name):
        self.fig.suptitle(name)
        self.canvas.draw()
        self.oldMainWindow.rename(name)
        
    def buttonPress(self,event):
        self.current.buttonPress(event)

    def buttonRelease(self,event):
        self.current.buttonRelease(event)

    def pan(self,event):
        self.current.pan(event)

    def scroll(self,event):
        self.current.scroll(event)
        
    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()
    
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
        self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class RelaxFrame(Plot1DFrame): 
    def __init__(self, rootwindow,fig,canvas,current):
        axAdd=0
        self.ref = current.ref
        self.axType = current.axType
        self.freq = current.freq
        self.xax = current.xax
        self.data1D=current.getDisplayedData()
        self.plotType = 0
        self.logx = 0
        self.logy = 0
        Plot1DFrame.__init__(self,rootwindow,fig,canvas)
        self.current = current
        self.rootwindow = rootwindow
        self.plotReset()
        self.showPlot()

    def plotReset(self): 
        if self.plotType==0:
            miny = min(np.real(self.data1D))
            maxy = max(np.real(self.data1D))
        elif self.plotType==1:
            miny = min(np.imag(self.data1D))
            maxy = max(np.imag(self.data1D))
        elif self.plotType==2:
            miny = min(min(np.real(self.data1D)),min(np.imag(self.data1D)))
            maxy = max(max(np.real(self.data1D)),max(np.imag(self.data1D)))
        elif self.plotType==3:
            miny = min(np.abs(self.data1D))
            maxy = max(np.abs(self.data1D))
        else:
            miny=-1
            maxy=1
        differ = 0.05*(maxy-miny) 
        self.yminlim=miny-differ
        self.ymaxlim=maxy+differ
        axAdd = 0.0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.xminlim=min(self.xax*axMult+axAdd)
        self.xmaxlim=max(self.xax*axMult+axAdd)
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)

    def showPlot(self, tmpAx=None, tmpdata=None): 
        self.ax.cla()
        axAdd = 0.0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if tmpAx is not None:
            self.ax.plot(tmpAx*axMult+axAdd,tmpdata)
        self.ax.scatter(self.xax*axMult+axAdd,self.data1D)
        if self.logx==0:
            self.ax.set_xscale('linear')
        else:
           self.ax.set_xscale('log')
        if self.logy==0:
            self.ax.set_yscale('linear')
        else:
           self.ax.set_yscale('log') 
        if self.spec==0:
            if self.current.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec==1:
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
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        if self.logx==0:
            self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.logy==0:
            self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.canvas.draw()
        
    def scroll(self,event):
        if self.rightMouse:
            if self.logx == 0:
                middle = (self.xmaxlim+self.xminlim)/2.0
                width = self.xmaxlim-self.xminlim
                width = width*0.9**event.step
                self.xmaxlim = middle+width/2.0
                self.xminlim = middle-width/2.0
                if self.spec > 0 and not isinstance(self,spectrum_classes.CurrentArrayed):
                    self.ax.set_xlim(self.xmaxlim,self.xminlim)
                else:
                    self.ax.set_xlim(self.xminlim,self.xmaxlim)
            else:
                middle = (np.log(self.xmaxlim)+np.log(self.xminlim))/2.0
                width = np.log(self.xmaxlim)-np.log(self.xminlim)
                width = width*0.9**event.step
                self.xmaxlim = np.exp(middle+width/2.0)
                self.xminlim = np.exp(middle-width/2.0)
                if self.spec > 0 and not isinstance(self,spectrum_classes.CurrentArrayed):
                    self.ax.set_xlim(self.xmaxlim,self.xminlim)
                else:
                    self.ax.set_xlim(self.xminlim,self.xmaxlim)
        else:
            if self.logy == 0:
                middle = (self.ymaxlim+self.yminlim)/2.0
                width = self.ymaxlim-self.yminlim
                width = width*0.9**event.step
                self.ymaxlim = middle+width/2.0
                self.yminlim = middle-width/2.0
                if self.spec2 > 0 and isinstance(self,spectrum_classes.CurrentContour):
                    self.ax.set_ylim(self.ymaxlim,self.yminlim)
                else:
                    self.ax.set_ylim(self.yminlim,self.ymaxlim)
            else:
                middle = (np.log(self.ymaxlim)+np.log(self.yminlim))/2.0
                width = np.log(self.ymaxlim)-np.log(self.yminlim)
                width = width*0.9**event.step
                self.ymaxlim = np.exp(middle+width/2.0)
                self.yminlim = np.exp(middle-width/2.0)
                if self.spec2 > 0 and isinstance(self,spectrum_classes.CurrentContour):
                    self.ax.set_ylim(self.ymaxlim,self.yminlim)
                else:
                    self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()
        
    def buttonRelease(self,event):
        if event.button == 1:
            if self.peakPick:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                    self.rect[0]=None
                    self.peakPick = False
                    idx = np.argmin(np.abs(self.line_xdata-event.xdata))
                    if self.peakPickFunc is not None:
                        self.peakPickFunc((idx,self.line_xdata[idx],self.line_ydata[idx]))
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
                self.rect=[None,None,None,None]
                if self.zoomX2 is not None and self.zoomY2 is not None:
                    self.xminlim=min([self.zoomX1,self.zoomX2])
                    self.xmaxlim=max([self.zoomX1,self.zoomX2])
                    self.yminlim=min([self.zoomY1,self.zoomY2])
                    self.ymaxlim=max([self.zoomY1,self.zoomY2])
                    if self.spec > 0 and not isinstance(self,spectrum_classes.CurrentArrayed):
                        self.ax.set_xlim(self.xmaxlim,self.xminlim)
                    else:
                        self.ax.set_xlim(self.xminlim,self.xmaxlim)
                    if self.spec2 > 0 and isinstance(self,spectrum_classes.CurrentContour):
                        self.ax.set_ylim(self.ymaxlim,self.yminlim)
                    else:
                        self.ax.set_ylim(self.yminlim,self.ymaxlim)
                self.zoomX1=None
                self.zoomX2=None 
                self.zoomY1=None
                self.zoomY2=None 
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw()

    def pan(self,event):
        if self.rightMouse and self.panX is not None and self.panY is not None:
            if self.logx==0 and self.logy==0:
                inv = self.ax.transData.inverted()
                point = inv.transform((event.x,event.y))
                x=point[0]
                y=point[1]
            else:
                x=event.xdata
                y=event.ydata
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
            if self.spec > 0 and not isinstance(self,spectrum_classes.CurrentArrayed):
                self.ax.set_xlim(self.xmaxlim,self.xminlim)
            else:
                self.ax.set_xlim(self.xminlim,self.xmaxlim)
            if self.spec2 > 0 and isinstance(self,spectrum_classes.CurrentContour):
                self.ax.set_ylim(self.ymaxlim,self.yminlim)
            else:
                self.ax.set_ylim(self.yminlim,self.ymaxlim)
            self.canvas.draw()
        elif self.peakPick:
            if self.rect[0] is not None:
                self.rect[0].remove()
                self.rect[0]=None
            if event.xdata is not None:
                self.rect[0]=self.ax.axvline(event.xdata,c='k',linestyle='--')
            self.canvas.draw()
        elif self.leftMouse and (self.zoomX1 is not None) and (self.zoomY1 is not None):
            if self.logx==0 and self.logy==0:
                inv = self.ax.transData.inverted()
                point = inv.transform((event.x,event.y))
                self.zoomX2=point[0]
                self.zoomY2=point[1]
            else:
                self.zoomX2=event.xdata
                self.zoomY2=event.ydata
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
                self.rect=[None,None,None,None]
            self.rect[0],=self.ax.plot([self.zoomX1,self.zoomX2],[self.zoomY2,self.zoomY2],'k',clip_on=False)
            self.rect[1],=self.ax.plot([self.zoomX1,self.zoomX2],[self.zoomY1,self.zoomY1],'k',clip_on=False)
            self.rect[2],=self.ax.plot([self.zoomX1,self.zoomX1],[self.zoomY1,self.zoomY2],'k',clip_on=False)
            self.rect[3],=self.ax.plot([self.zoomX2,self.zoomX2],[self.zoomY1,self.zoomY2],'k',clip_on=False)
            self.canvas.draw()
            
    def setLog(self,logx,logy):
        self.logx = logx
        self.logy = logy
        if self.logx==0:
            self.ax.set_xscale('linear')
        else:
           self.ax.set_xscale('log')
        if self.logy==0:
            self.ax.set_yscale('linear')
        else:
           self.ax.set_yscale('log') 
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()

#################################################################################
class RelaxParamFrame(Frame): 
    def __init__(self, parent, rootwindow): 
        self.parent = parent
        self.rootwindow = rootwindow
        self.ampVal = StringVar()
        self.ampVal.set("%.3g" % np.amax(self.parent.data1D))
        self.ampTick = IntVar()
        self.constVal = StringVar()
        self.constVal.set("1.0")
        self.constTick = IntVar()
        self.numExp = StringVar()
        self.xlog=IntVar()
        self.ylog=IntVar()
        Frame.__init__(self, rootwindow)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0,sticky='n')
        self.optframe = Frame(self)
        self.optframe.grid(row=0,column=1,sticky='n')
        self.frame2 = Frame(self)
        self.frame2.grid(row=0,column=2,sticky='n')
        self.frame3 = Frame(self)
        self.frame3.grid(row=0,column=3,sticky='n')
        Button(self.frame1, text="Sim",command=self.sim).grid(row=0,column=0)
        Button(self.frame1, text="Fit",command=self.fit).grid(row=1,column=0)
        Button(self.frame1, text="Fit all",command=self.fitAll).grid(row=2,column=0)
        Button(self.frame1, text="Cancel",command=rootwindow.cancel).grid(row=3,column=0)
        Label(self.frame2,text="Amplitude").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.ampTick).grid(row=1,column=0)
        Entry(self.frame2,textvariable=self.ampVal,justify="center",width=10).grid(row=1,column=1)
        Label(self.frame2,text="Constant").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.constTick).grid(row=3,column=0)
        Entry(self.frame2,textvariable=self.constVal,justify="center",width=10).grid(row=3,column=1)
        OptionMenu(self.frame3, self.numExp, "1","1", "2", "3","4",command=self.changeNum).grid(row=0,column=0,columnspan=4)
        Label(self.frame3,text="Coefficient").grid(row=1,column=0,columnspan=2)
        Label(self.frame3,text="T [s]").grid(row=1,column=2,columnspan=2)
        Checkbutton(self.optframe,variable=self.xlog,text='x-log',command=self.setLog).grid(row=0,column=0)
        Checkbutton(self.optframe,variable=self.ylog,text='y-log',command=self.setLog).grid(row=1,column=0)
        
        self.coeffVal = []
        self.coeffTick = []
        self.T1Val = []
        self.T1Tick = []
        self.coeffCheck = []
        self.coeffEntries = []
        self.T1Check = []
        self.T1Entries = []
        for i in range(4):
            self.coeffVal.append(StringVar())
            self.coeffVal[i].set("-1.0")
            self.coeffTick.append(IntVar())
            self.T1Val.append(StringVar())
            self.T1Val[i].set("1.0")
            self.T1Tick.append(IntVar())
            self.coeffCheck.append(Checkbutton(self.frame3,variable=self.coeffTick[i]))
            self.coeffCheck[i].grid(row=i+2,column=0)
            self.coeffEntries.append(Entry(self.frame3,textvariable=self.coeffVal[i],justify="center",width=10))
            self.coeffEntries[i].grid(row=i+2,column=1)
            self.T1Check.append(Checkbutton(self.frame3,variable=self.T1Tick[i]))
            self.T1Check[i].grid(row=i+2,column=2)
            self.T1Entries.append(Entry(self.frame3,textvariable=self.T1Val[i],justify="center",width=10))
            self.T1Entries[i].grid(row=i+2,column=3)
            if i > 0:
                self.coeffCheck[i].grid_remove()
                self.coeffEntries[i].grid_remove()
                self.T1Check[i].grid_remove()
                self.T1Entries[i].grid_remove()

    def setLog(self, *args):
        self.parent.setLog(self.xlog.get(),self.ylog.get())
                
    def changeNum(self,*args):
        val = int(self.numExp.get())
        for i in range(4):
            if i < val:
                self.coeffCheck[i].grid()
                self.coeffEntries[i].grid()
                self.T1Check[i].grid()
                self.T1Entries[i].grid()
            else:
                self.coeffCheck[i].grid_remove()
                self.coeffEntries[i].grid_remove()
                self.T1Check[i].grid_remove()
                self.T1Entries[i].grid_remove()

    def fitFunc(self, x, *param):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        testFunc = np.zeros(len(x))
        if struc[0]:
            amplitude = param[0]
            param=np.delete(param,[0])
        else:
            amplitude = argu[0]
            argu=np.delete(argu,[0])
        if struc[1]:
            constant = param[0]
            param=np.delete(param,[0])
        else:
            constant = argu[0]
            argu=np.delete(argu,[0])
        for i in range(1,numExp+1):
            if struc[2*i]:
                coeff = param[0]
                param=np.delete(param,[0])
            else:
                coeff= argu[0]
                argu=np.delete(argu,[0])
            if struc[2*i+1]:
                T1 = param[0]
                param=np.delete(param,[0])
            else:
                T1= argu[0]
                argu=np.delete(argu,[0])
            testFunc += coeff*np.exp(-1.0*x/T1) 
        return amplitude*(constant+testFunc)

    def fit(self,*args):
        struc = []
        guess = []
        argu = []
        numExp = int(self.numExp.get())
        outCoeff = np.zeros(numExp)
        outT1 = np.zeros(numExp)
        if self.ampTick.get() == 0:
            guess.append(safeEval(self.ampVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.ampVal.get())
            argu.append(inp)
            outAmp = inp
            self.ampVal.set('%.3g' % inp)
            struc.append(False)
        if self.constTick.get() == 0:
            guess.append(safeEval(self.constVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.constVal.get())
            argu.append(inp)
            outConst = inp
            self.constVal.set('%.3g' % inp)
            struc.append(False)
        for i in range(numExp):
            if self.coeffTick[i].get() == 0:
                guess.append(safeEval(self.coeffVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.coeffVal[i].get())
                argu.append(inp)
                outCoeff[i] = inp
                self.coeffVal[i].set('%.3g' % inp)
                struc.append(False)
            if self.T1Tick[i].get() == 0:
                guess.append(safeEval(self.T1Val[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.T1Val[i].get())
                argu.append(inp)
                outT1[i] = inp
                self.T1Val[i].set('%.3g' % inp)
                struc.append(False)
        self.args=(numExp,struc,argu)
        fitVal = scipy.optimize.curve_fit(self.fitFunc,self.parent.xax, self.parent.data1D,guess)
        counter = 0
        if struc[0]:
            self.ampVal.set('%.3g' % fitVal[0][counter])
            outAmp = fitVal[0][counter]
            counter +=1
        if struc[1]:
            self.constVal.set('%.3g' % fitVal[0][counter])
            outConst = fitVal[0][counter]
            counter +=1
        for i in range(1,numExp+1):
            if struc[2*i]:
                self.coeffVal[i-1].set('%.3g' % fitVal[0][counter])
                outCoeff[i-1] = fitVal[0][counter]
                counter += 1
            if struc[2*i+1]:
                self.T1Val[i-1].set('%.3g' % fitVal[0][counter])
                outT1[i-1] = fitVal[0][counter]
                counter += 1
        self.disp(outAmp,outConst,outCoeff,outT1)

    def fitAll(self, *args):
        FitAllSelectionWindow(self,["Amplitude","Constant","Coefficient","T"])
        
    def fitAllFunc(self,outputs):
        struc = []
        guess = []
        argu = []
        numExp = int(self.numExp.get())
        outCoeff = np.zeros(numExp)
        outT1 = np.zeros(numExp)
        if self.ampTick.get() == 0:
            guess.append(safeEval(self.ampVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.ampVal.get())
            argu.append(inp)
            outAmp = inp
            self.ampVal.set('%.3g' % inp)
            struc.append(False)
        if self.constTick.get() == 0:
            guess.append(safeEval(self.constVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.constVal.get())
            argu.append(inp)
            outConst = inp
            self.constVal.set('%.3g' % inp)
            struc.append(False)
        for i in range(numExp):
            if self.coeffTick[i].get() == 0:
                guess.append(safeEval(self.coeffVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.coeffVal[i].get())
                argu.append(inp)
                outCoeff[i] = inp
                self.coeffVal[i].set('%.3g' % inp)
                struc.append(False)
            if self.T1Tick[i].get() == 0:
                guess.append(safeEval(self.T1Val[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.T1Val[i].get())
                argu.append(inp)
                outT1[i] = inp
                self.T1Val[i].set('%.3g' % inp)
                struc.append(False)
        self.args=(numExp,struc,argu)
        fullData = self.parent.current.data.data
        axes = self.parent.current.axes
        dataShape = fullData.shape
        dataShape2 = np.delete(dataShape,axes)
        rolledData = np.rollaxis(fullData,axes)
        intOutputs = np.array(outputs,dtype=int)
        numOutputs = np.sum(intOutputs[:2]) + numExp*np.sum(intOutputs[2:])
        outputData = np.zeros((np.product(dataShape2),numOutputs),dtype=complex)
        counter2 = 0
        for j in rolledData.reshape(dataShape[axes],np.product(dataShape2)).T:
            try:
                fitVal = scipy.optimize.curve_fit(self.fitFunc,self.parent.xax, np.real(j),guess)
            except:
                fitVal = [[0]*10]
            counter = 0
            if struc[0]:
                outAmp = fitVal[0][counter]
                counter +=1
            if struc[1]:
                outConst = fitVal[0][counter]
                counter +=1
            for i in range(1,numExp+1):
                if struc[2*i]:
                    outCoeff[i-1] = fitVal[0][counter]
                    counter += 1
                if struc[2*i+1]:
                    outT1[i-1] = fitVal[0][counter]
                    counter += 1
            outputArray = []
            if outputs[0]:
                outputArray = np.concatenate((outputArray,[outAmp]))
            if outputs[1]:
                outputArray = np.concatenate((outputArray,[outConst]))
            if outputs[2]:
                outputArray = np.concatenate((outputArray,outCoeff))
            if outputs[3]:
                outputArray = np.concatenate((outputArray,outT1))
            outputData[counter2] = outputArray
            counter2 += 1
        newShape = np.concatenate((np.array(dataShape2),[numOutputs]))
        self.rootwindow.createNewData(np.rollaxis(outputData.reshape(newShape),-1,axes), axes)
        
    def sim(self):
        numExp = int(self.numExp.get())
        outAmp = safeEval(self.ampVal.get())
        outConst = safeEval(self.constVal.get())
        outCoeff = []
        outT1 = []
        for i in range(numExp):
            outCoeff.append(safeEval(self.coeffVal[i].get()))
            outT1.append(safeEval(self.T1Val[i].get()))
        self.disp(outAmp,outConst,outCoeff,outT1)
        
    def disp(self, outAmp,outConst,outCoeff, outT1):
        numCurve = 100 #number of points in output curve
        outCurve = np.zeros(numCurve)
        if self.xlog.get() == 1:
            x = np.logspace(np.log(min(self.parent.xax)),np.log(max(self.parent.xax)),numCurve)
        else:
            x = np.linspace(min(self.parent.xax),max(self.parent.xax),numCurve)
        for i in range(len(outCoeff)):
            outCurve += outCoeff[i]*np.exp(-x/outT1[i])
        self.parent.showPlot(x, outAmp*(outConst+outCurve))

##############################################################################
class DiffusionWindow(Frame): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        self.fig = self.mainProgram.getFig()
        self.canvas = FigureCanvasTkAgg(self.fig, master=weakref.proxy(self))
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="nswe")
        self.current = DiffusionFrame(self,self.fig,self.canvas,oldMainWindow.current)
        self.paramframe = DiffusionParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)

    def rename(self,name):
        self.fig.suptitle(name)
        self.canvas.draw()
        self.oldMainWindow.rename(name)
        
    def buttonPress(self,event):
        self.current.buttonPress(event)

    def buttonRelease(self,event):
        self.current.buttonRelease(event)

    def pan(self,event):
        self.current.pan(event)

    def scroll(self,event):
        self.current.scroll(event)
        
    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()
    
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
        self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class DiffusionFrame(Plot1DFrame): 
    def __init__(self, rootwindow,fig,canvas,current):
        axAdd=0
        self.ref = current.ref
        self.axType = current.axType
        self.freq = current.freq
        self.xax = current.xax
        self.data1D=current.getDisplayedData()
        self.plotType = 0
        self.logx = 0
        self.logy = 0
        Plot1DFrame.__init__(self,rootwindow,fig,canvas)
        self.current = current
        self.rootwindow = rootwindow
        self.plotReset()
        self.showPlot()

    def plotReset(self): 
        if self.plotType==0:
            miny = min(np.real(self.data1D))
            maxy = max(np.real(self.data1D))
        elif self.plotType==1:
            miny = min(np.imag(self.data1D))
            maxy = max(np.imag(self.data1D))
        elif self.plotType==2:
            miny = min(min(np.real(self.data1D)),min(np.imag(self.data1D)))
            maxy = max(max(np.real(self.data1D)),max(np.imag(self.data1D)))
        elif self.plotType==3:
            miny = min(np.abs(self.data1D))
            maxy = max(np.abs(self.data1D))
        else:
            miny=-1
            maxy=1
        differ = 0.05*(maxy-miny) 
        self.yminlim=miny-differ
        self.ymaxlim=maxy+differ
        axAdd = 0.0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.xminlim=min(self.xax*axMult+axAdd)
        self.xmaxlim=max(self.xax*axMult+axAdd)
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)

    def showPlot(self, tmpAx=None, tmpdata=None): 
        self.ax.cla()
        axAdd = 0.0
        if self.spec == 1:
            if self.ppm:
                axAdd = (self.freq-self.ref)/self.ref*1e6
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if tmpAx is not None:
            self.ax.plot(tmpAx*axMult+axAdd,tmpdata)
        self.ax.scatter(self.xax*axMult+axAdd,self.data1D)
        if self.logx==0:
            self.ax.set_xscale('linear')
        else:
           self.ax.set_xscale('log')
        if self.logy==0:
            self.ax.set_yscale('linear')
        else:
           self.ax.set_yscale('log') 
        if self.spec==0:
            if self.current.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec==1:
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
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        if self.logx==0:
            self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.logy==0:
            self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.canvas.draw()
        
    def scroll(self,event):
        if self.rightMouse:
            if self.logx == 0:
                middle = (self.xmaxlim+self.xminlim)/2.0
                width = self.xmaxlim-self.xminlim
                width = width*0.9**event.step
                self.xmaxlim = middle+width/2.0
                self.xminlim = middle-width/2.0
                if self.spec > 0 and not isinstance(self,spectrum_classes.CurrentArrayed):
                    self.ax.set_xlim(self.xmaxlim,self.xminlim)
                else:
                    self.ax.set_xlim(self.xminlim,self.xmaxlim)
            else:
                middle = (np.log(self.xmaxlim)+np.log(self.xminlim))/2.0
                width = np.log(self.xmaxlim)-np.log(self.xminlim)
                width = width*0.9**event.step
                self.xmaxlim = np.exp(middle+width/2.0)
                self.xminlim = np.exp(middle-width/2.0)
                if self.spec > 0 and not isinstance(self,spectrum_classes.CurrentArrayed):
                    self.ax.set_xlim(self.xmaxlim,self.xminlim)
                else:
                    self.ax.set_xlim(self.xminlim,self.xmaxlim)
        else:
            if self.logy == 0:
                middle = (self.ymaxlim+self.yminlim)/2.0
                width = self.ymaxlim-self.yminlim
                width = width*0.9**event.step
                self.ymaxlim = middle+width/2.0
                self.yminlim = middle-width/2.0
                if self.spec2 > 0 and isinstance(self,spectrum_classes.CurrentContour):
                    self.ax.set_ylim(self.ymaxlim,self.yminlim)
                else:
                    self.ax.set_ylim(self.yminlim,self.ymaxlim)
            else:
                middle = (np.log(self.ymaxlim)+np.log(self.yminlim))/2.0
                width = np.log(self.ymaxlim)-np.log(self.yminlim)
                width = width*0.9**event.step
                self.ymaxlim = np.exp(middle+width/2.0)
                self.yminlim = np.exp(middle-width/2.0)
                if self.spec2 > 0 and isinstance(self,spectrum_classes.CurrentContour):
                    self.ax.set_ylim(self.ymaxlim,self.yminlim)
                else:
                    self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()
        
    def buttonRelease(self,event):
        if event.button == 1:
            if self.peakPick:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                    self.rect[0]=None
                    self.peakPick = False
                    idx = np.argmin(np.abs(self.line_xdata-event.xdata))
                    if self.peakPickFunc is not None:
                        self.peakPickFunc((idx,self.line_xdata[idx],self.line_ydata[idx]))
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
                self.rect=[None,None,None,None]
                if self.zoomX2 is not None and self.zoomY2 is not None:
                    self.xminlim=min([self.zoomX1,self.zoomX2])
                    self.xmaxlim=max([self.zoomX1,self.zoomX2])
                    self.yminlim=min([self.zoomY1,self.zoomY2])
                    self.ymaxlim=max([self.zoomY1,self.zoomY2])
                    if self.spec > 0 and not isinstance(self,spectrum_classes.CurrentArrayed):
                        self.ax.set_xlim(self.xmaxlim,self.xminlim)
                    else:
                        self.ax.set_xlim(self.xminlim,self.xmaxlim)
                    if self.spec2 > 0 and isinstance(self,spectrum_classes.CurrentContour):
                        self.ax.set_ylim(self.ymaxlim,self.yminlim)
                    else:
                        self.ax.set_ylim(self.yminlim,self.ymaxlim)
                self.zoomX1=None
                self.zoomX2=None 
                self.zoomY1=None
                self.zoomY2=None 
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw()

    def pan(self,event):
        if self.rightMouse and self.panX is not None and self.panY is not None:
            if self.logx==0 and self.logy==0:
                inv = self.ax.transData.inverted()
                point = inv.transform((event.x,event.y))
                x=point[0]
                y=point[1]
            else:
                x=event.xdata
                y=event.ydata
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
            if self.spec > 0 and not isinstance(self,spectrum_classes.CurrentArrayed):
                self.ax.set_xlim(self.xmaxlim,self.xminlim)
            else:
                self.ax.set_xlim(self.xminlim,self.xmaxlim)
            if self.spec2 > 0 and isinstance(self,spectrum_classes.CurrentContour):
                self.ax.set_ylim(self.ymaxlim,self.yminlim)
            else:
                self.ax.set_ylim(self.yminlim,self.ymaxlim)
            self.canvas.draw()
        elif self.peakPick:
            if self.rect[0] is not None:
                self.rect[0].remove()
                self.rect[0]=None
            if event.xdata is not None:
                self.rect[0]=self.ax.axvline(event.xdata,c='k',linestyle='--')
            self.canvas.draw()
        elif self.leftMouse and (self.zoomX1 is not None) and (self.zoomY1 is not None):
            if self.logx==0 and self.logy==0:
                inv = self.ax.transData.inverted()
                point = inv.transform((event.x,event.y))
                self.zoomX2=point[0]
                self.zoomY2=point[1]
            else:
                self.zoomX2=event.xdata
                self.zoomY2=event.ydata
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
                self.rect=[None,None,None,None]
            self.rect[0],=self.ax.plot([self.zoomX1,self.zoomX2],[self.zoomY2,self.zoomY2],'k',clip_on=False)
            self.rect[1],=self.ax.plot([self.zoomX1,self.zoomX2],[self.zoomY1,self.zoomY1],'k',clip_on=False)
            self.rect[2],=self.ax.plot([self.zoomX1,self.zoomX1],[self.zoomY1,self.zoomY2],'k',clip_on=False)
            self.rect[3],=self.ax.plot([self.zoomX2,self.zoomX2],[self.zoomY1,self.zoomY2],'k',clip_on=False)
            self.canvas.draw()
            
    def setLog(self,logx,logy):
        self.logx = logx
        self.logy = logy
        if self.logx==0:
            self.ax.set_xscale('linear')
        else:
           self.ax.set_xscale('log')
        if self.logy==0:
            self.ax.set_yscale('linear')
        else:
           self.ax.set_yscale('log') 
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()

#################################################################################
class DiffusionParamFrame(Frame): 
    def __init__(self, parent, rootwindow): 
        self.parent = parent
        self.gammaVal = StringVar()
        self.gammaVal.set("42.576")
        self.deltaVal = StringVar()
        self.deltaVal.set("1.0")
        self.triangleVal = StringVar()
        self.triangleVal.set("1.0")
        self.ampVal = StringVar()
        self.ampVal.set("%.3g" % np.amax(self.parent.data1D))
        self.ampTick = IntVar()
        self.constVal = StringVar()
        self.constVal.set("0.0")
        self.constTick = IntVar()
        self.constTick.set(1)
        self.numExp = StringVar()
        self.xlog=IntVar()
        self.ylog=IntVar()
        Frame.__init__(self, rootwindow)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0,sticky='n')
        self.optframe = Frame(self)
        self.optframe.grid(row=0,column=1,sticky='n')
        self.frame2 = Frame(self)
        self.frame2.grid(row=0,column=2,sticky='n')
        self.frame3 = Frame(self)
        self.frame3.grid(row=0,column=3,sticky='n')
        self.frame4 = Frame(self)
        self.frame4.grid(row=0,column=4,sticky='n')
        Button(self.frame1, text="Sim",command=self.sim).grid(row=0,column=0)
        Button(self.frame1, text="Fit",command=self.fit).grid(row=1,column=0)
        Button(self.frame1, text="Cancel",command=rootwindow.cancel).grid(row=2,column=0)
        Label(self.frame2,text=u"\u03b3 [MHz/T]").grid(row=0,column=0)
        Entry(self.frame2,textvariable=self.gammaVal,justify="center",width=10).grid(row=1,column=0)
        Label(self.frame2,text=u"\u03b4 [s]").grid(row=2,column=0)
        Entry(self.frame2,textvariable=self.deltaVal,justify="center",width=10).grid(row=3,column=0)
        Label(self.frame2,text=u"\u0394 [s]").grid(row=4,column=0)
        Entry(self.frame2,textvariable=self.triangleVal,justify="center",width=10).grid(row=5,column=0)
        Label(self.frame3,text="Amplitude").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame3,variable=self.ampTick).grid(row=1,column=0)
        Entry(self.frame3,textvariable=self.ampVal,justify="center",width=10).grid(row=1,column=1)
        Label(self.frame3,text="Constant").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame3,variable=self.constTick).grid(row=3,column=0)
        Entry(self.frame3,textvariable=self.constVal,justify="center",width=10).grid(row=3,column=1)
        OptionMenu(self.frame4, self.numExp, "1","1", "2", "3","4",command=self.changeNum).grid(row=0,column=0,columnspan=4)
        Label(self.frame4,text="Coefficient").grid(row=1,column=0,columnspan=2)
        Label(self.frame4,text="D").grid(row=1,column=2,columnspan=2)
        Checkbutton(self.optframe,variable=self.xlog,text='x-log',command=self.setLog).grid(row=0,column=0)
        Checkbutton(self.optframe,variable=self.ylog,text='y-log',command=self.setLog).grid(row=1,column=0)
        
        self.coeffVal = []
        self.coeffTick = []
        self.DVal = []
        self.DTick = []
        self.coeffCheck = []
        self.coeffEntries = []
        self.DCheck = []
        self.DEntries = []
        for i in range(4):
            self.coeffVal.append(StringVar())
            self.coeffVal[i].set("1.0")
            self.coeffTick.append(IntVar())
            self.DVal.append(StringVar())
            self.DVal[i].set("1.0e-9")
            self.DTick.append(IntVar())
            self.coeffCheck.append(Checkbutton(self.frame4,variable=self.coeffTick[i]))
            self.coeffCheck[i].grid(row=i+2,column=0)
            self.coeffEntries.append(Entry(self.frame4,textvariable=self.coeffVal[i],justify="center",width=10))
            self.coeffEntries[i].grid(row=i+2,column=1)
            self.DCheck.append(Checkbutton(self.frame4,variable=self.DTick[i]))
            self.DCheck[i].grid(row=i+2,column=2)
            self.DEntries.append(Entry(self.frame4,textvariable=self.DVal[i],justify="center",width=10))
            self.DEntries[i].grid(row=i+2,column=3)
            if i > 0:
                self.coeffCheck[i].grid_remove()
                self.coeffEntries[i].grid_remove()
                self.DCheck[i].grid_remove()
                self.DEntries[i].grid_remove()

    def setLog(self, *args):
        self.parent.setLog(self.xlog.get(),self.ylog.get())
                
    def changeNum(self,*args):
        val = int(self.numExp.get())
        for i in range(4):
            if i < val:
                self.coeffCheck[i].grid()
                self.coeffEntries[i].grid()
                self.DCheck[i].grid()
                self.DEntries[i].grid()
            else:
                self.coeffCheck[i].grid_remove()
                self.coeffEntries[i].grid_remove()
                self.DCheck[i].grid_remove()
                self.DEntries[i].grid_remove()

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
            param=np.delete(param,[0])
        else:
            amplitude = argu[0]
            argu=np.delete(argu,[0])
        if struc[1]:
            constant = param[0]
            param=np.delete(param,[0])
        else:
            constant = argu[0]
            argu=np.delete(argu,[0])
        for i in range(1,numExp+1):
            if struc[2*i]:
                coeff = param[0]
                param=np.delete(param,[0])
            else:
                coeff= argu[0]
                argu=np.delete(argu,[0])
            if struc[2*i+1]:
                D = param[0]
                param=np.delete(param,[0])
            else:
                D = argu[0]
                argu=np.delete(argu,[0])
            testFunc += coeff*np.exp(-(gamma*delta*x)**2*D*(triangle-delta/3.0)) 
        return amplitude*(constant+testFunc)

    def fit(self,*args):
        struc = []
        guess = []
        argu = []
        numExp = int(self.numExp.get())
        outCoeff = np.zeros(numExp)
        outD = np.zeros(numExp)
        gamma = safeEval(self.gammaVal.get())
        self.gammaVal.set('%.3g' % gamma)
        delta = safeEval(self.deltaVal.get())
        self.deltaVal.set('%.3g' % delta)
        triangle = safeEval(self.triangleVal.get())
        self.triangleVal.set('%.3g' % triangle)
        if self.ampTick.get() == 0:
            guess.append(safeEval(self.ampVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.ampVal.get())
            argu.append(inp)
            outAmp = inp
            self.ampVal.set('%.3g' % inp)
            struc.append(False)
        if self.constTick.get() == 0:
            guess.append(safeEval(self.constVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.constVal.get())
            argu.append(inp)
            outConst = inp
            self.constVal.set('%.3g' % inp)
            struc.append(False)
        for i in range(numExp):
            if self.coeffTick[i].get() == 0:
                guess.append(safeEval(self.coeffVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.coeffVal[i].get())
                argu.append(inp)
                outCoeff[i] = inp
                self.coeffVal[i].set('%.3g' % inp)
                struc.append(False)
            if self.DTick[i].get() == 0:
                guess.append(safeEval(self.DVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.DVal[i].get())
                argu.append(inp)
                outD[i] = inp
                self.DVal[i].set('%.3g' % inp)
                struc.append(False)
        self.args=(numExp, struc, argu, gamma, delta, triangle)
        fitVal = scipy.optimize.curve_fit(self.fitFunc, self.parent.xax, self.parent.data1D, guess)
        counter = 0
        if struc[0]:
            self.ampVal.set('%.3g' % fitVal[0][counter])
            outAmp = fitVal[0][counter]
            counter +=1
        if struc[1]:
            self.constVal.set('%.3g' % fitVal[0][counter])
            outConst = fitVal[0][counter]
            counter +=1
        for i in range(1,numExp+1):
            if struc[2*i]:
                self.coeffVal[i-1].set('%.3g' % fitVal[0][counter])
                outCoeff[i-1] = fitVal[0][counter]
                counter += 1
            if struc[2*i+1]:
                self.DVal[i-1].set('%.3g' % fitVal[0][counter])
                outD[i-1] = fitVal[0][counter]
                counter += 1
        self.disp(outAmp, outConst, outCoeff, outD, gamma, delta, triangle)

    def sim(self):
        numExp = int(self.numExp.get())
        outAmp = safeEval(self.ampVal.get())
        outConst = safeEval(self.constVal.get())
        gamma = safeEval(self.gammaVal.get())
        delta = safeEval(self.deltaVal.get())
        triangle = safeEval(self.triangleVal.get())
        outCoeff = []
        outD = []
        for i in range(numExp):
            outCoeff.append(safeEval(self.coeffVal[i].get()))
            outD.append(safeEval(self.DVal[i].get()))
        self.disp(outAmp,outConst,outCoeff,outD,gamma,delta,triangle)
        
    def disp(self, outAmp, outConst, outCoeff, outD, gamma, delta, triangle):
        numCurve = 100 
        outCurve = np.zeros(numCurve)
        if self.xlog.get() == 1:
            x = np.logspace(np.log(min(self.parent.xax)),np.log(max(self.parent.xax)),numCurve)
        else:
            x = np.linspace(min(self.parent.xax),max(self.parent.xax),numCurve)
        for i in range(len(outCoeff)):
            outCurve += outCoeff[i]*np.exp(-(gamma*delta*x)**2*outD[i]*(triangle-delta/3.0))
        self.parent.showPlot(x, outAmp*(outConst+outCurve))
        
##############################################################################
class PeakDeconvWindow(Frame): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        self.fig = self.mainProgram.getFig()
        self.canvas = FigureCanvasTkAgg(self.fig, master=weakref.proxy(self))
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="nswe")
        self.current = PeakDeconvFrame(self,self.fig,self.canvas,oldMainWindow.current)
        self.paramframe = PeakDeconvParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        
    def rename(self,name):
        self.fig.suptitle(name)
        self.canvas.draw()
        self.oldMainWindow.rename(name)
        
    def buttonPress(self,event):
        self.current.buttonPress(event)

    def buttonRelease(self,event):
        self.current.buttonRelease(event)

    def pan(self,event):
        self.current.pan(event)

    def scroll(self,event):
        self.current.scroll(event)
        
    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class PeakDeconvFrame(Plot1DFrame): 
    def __init__(self, rootwindow,fig,canvas,current):
        Plot1DFrame.__init__(self,rootwindow,fig,canvas)
        self.data1D = current.getDisplayedData()
        self.current = current
        self.spec = self.current.spec
        self.xax = self.current.xax
        self.plotType=0
        self.rootwindow = rootwindow
        self.pickNum = 0
        self.pickWidth = False
        self.plotReset()
        self.showPlot()

    def plotReset(self): 
        a=self.fig.gca()
        if self.plotType==0:
            miny = min(np.real(self.data1D))
            maxy = max(np.real(self.data1D))
        elif self.plotType==1:
            miny = min(np.imag(self.data1D))
            maxy = max(np.imag(self.data1D))
        elif self.plotType==2:
            miny = min(min(np.real(self.data1D)),min(np.imag(self.data1D)))
            maxy = max(max(np.real(self.data1D)),max(np.imag(self.data1D)))
        elif self.plotType==3:
            miny = min(np.abs(self.data1D))
            maxy = max(np.abs(self.data1D))
        else:
            miny=-1
            maxy=1
        differ = 0.05*(maxy-miny) 
        self.yminlim=miny-differ
        self.ymaxlim=maxy+differ
        axAdd = 0
        if self.spec == 1:
            if self.current.ppm:
                axAdd = (self.current.freq-self.current.ref)/self.current.ref*1e6
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xminlim=min(self.xax*axMult+axAdd)
        self.xmaxlim=max(self.xax*axMult+axAdd)
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim,self.xminlim)
        else:
            a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        a=self.fig.gca()
        a.cla()
        axAdd = 0
        if self.spec == 1:
            if self.current.ppm:
                axAdd = (self.current.freq-self.current.ref)/self.current.ref*1e6
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.line_xdata = self.xax*axMult+axAdd
        self.line_ydata = self.data1D
        a.plot(self.xax*axMult+axAdd,self.data1D)
        if tmpAx is not None:
            a.plot(tmpAx*axMult+axAdd,tmpdata)
        for i in range(len(tmpAx2)):
            a.plot(tmpAx2[i]*axMult+axAdd,tmpdata2[i])
        if self.spec==0:
            if self.current.axType == 0:
                a.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                a.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                a.set_xlabel(r'Time [$\mu$s]')
            else:
                a.set_xlabel('User defined')
        elif self.spec==1:
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
            a.set_xlim(self.xmaxlim,self.xminlim)
        else:
            a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()

    def togglePick(self,var):
        if var==1:
            self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
            self.peakPick = True
        else:
            self.peakPickFunc = None
            self.peakPick = False
        
    def pickDeconv(self, pos):
        if self.pickWidth:
            axAdd = 0
            if self.current.spec == 1:
                if self.current.ppm:
                    axAdd = (self.current.freq-self.current.ref)/self.current.ref*1e6
                    axMult = 1e6/self.current.ref
                else:
                    axMult = 1.0/(1000.0**self.current.axType)
            elif self.current.spec == 0:
                axMult = 1000.0**self.current.axType
            width = (2*abs(float(self.rootwindow.paramframe.posVal[self.pickNum].get())-pos[1])-axAdd)/axMult
            self.rootwindow.paramframe.ampVal[self.pickNum].set("%.3g" %(float(self.rootwindow.paramframe.ampVal[self.pickNum].get())*width))
            self.rootwindow.paramframe.widthVal[self.pickNum].set("%.3g" % abs(width))
            self.pickNum += 1
            self.pickWidth = False
        else:
            self.rootwindow.paramframe.posVal[self.pickNum].set("%.3g" %pos[1])
            left = pos[0] - 10 
            if left < 0:
                left = 0
            right = pos[0] + 10
            if right >= len(self.data1D) :
                right = len(self.data1D)-1
            self.rootwindow.paramframe.ampVal[self.pickNum].set("%.3g" %(pos[2]*0.5*np.pi))
            if self.pickNum < 10:
                self.rootwindow.paramframe.numExp.set(str(self.pickNum+1))
                self.rootwindow.paramframe.changeNum()
            self.pickWidth = True
        if self.pickNum < 10:
            self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
            self.peakPick = True 

#################################################################################
class PeakDeconvParamFrame(Frame): #a frame for the relaxtion parameters
    def __init__(self, parent, rootwindow): 
        self.parent = parent
        self.axAdd = 0
        if self.parent.current.spec == 1:
            if self.parent.current.ppm:
                self.axAdd = (self.parent.current.freq-self.parent.current.ref)/self.parent.current.ref*1e6
                self.axMult = 1e6/self.parent.current.ref
            else:
                self.axMult = 1.0/(1000.0**self.parent.current.axType)
        elif self.parent.current.spec == 0:
            self.axMult = 1000.0**self.parent.current.axType
        self.bgrndVal = StringVar()
        self.bgrndTick = IntVar()
        self.slopeVal = StringVar()
        self.slopeTick = IntVar()
        self.numExp = StringVar()
        self.pickTick = IntVar()
        Frame.__init__(self, rootwindow)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0,sticky='n')
        self.frame2 = Frame(self)
        self.frame2.grid(row=0,column=1,sticky='n')
        self.frame3 = Frame(self)
        self.frame3.grid(row=0,column=2,sticky='n')
        Button(self.frame1, text="Sim",command=self.sim).grid(row=0,column=0)
        Button(self.frame1, text="Fit",command=self.fit).grid(row=1,column=0)
        Button(self.frame1, text="Cancel",command=rootwindow.cancel).grid(row=2,column=0)
        Button(self.frame1, text="Reset",command=self.reset).grid(row=0,column=1)
        Checkbutton(self.frame1, variable=self.pickTick,text="Picking",command=self.togglePick).grid(row=1,column=1)
        Label(self.frame2,text="Bgrnd").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.bgrndTick).grid(row=1,column=0)
        Entry(self.frame2,textvariable=self.bgrndVal,justify="center",width=10).grid(row=1,column=1)
        Label(self.frame2,text="Slope").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.slopeTick).grid(row=3,column=0)
        Entry(self.frame2,textvariable=self.slopeVal,justify="center",width=10).grid(row=3,column=1)
        OptionMenu(self.frame3, self.numExp, "1","1", "2", "3","4","5","6","7","8","9","10",command=self.changeNum).grid(row=0,column=0,columnspan=6)
        Label(self.frame3,text="Position").grid(row=1,column=0,columnspan=2)
        Label(self.frame3,text="Amplitude").grid(row=1,column=2,columnspan=2)
        Label(self.frame3,text="Lorentz [Hz]").grid(row=1,column=4,columnspan=2)
        Label(self.frame3,text="Gauss [Hz]").grid(row=1,column=6,columnspan=2)
        self.posVal = []
        self.posTick = []
        self.ampVal = []
        self.ampTick = []
        self.widthVal = []
        self.widthTick = []
        self.gaussVal = []
        self.gaussTick = []
        self.posCheck = []
        self.posEntries = []
        self.ampCheck = []
        self.ampEntries = []
        self.widthCheck = []
        self.widthEntries = []
        self.gaussCheck = []
        self.gaussEntries = []

        for i in range(10):
            self.posVal.append(StringVar())
            self.posTick.append(IntVar())
            self.ampVal.append(StringVar())
            self.ampTick.append(IntVar())
            self.widthVal.append(StringVar())
            self.widthTick.append(IntVar())
            self.gaussVal.append(StringVar())
            self.gaussTick.append(IntVar())
            self.posCheck.append(Checkbutton(self.frame3,variable=self.posTick[i]))
            self.posCheck[i].grid(row=i+2,column=0)
            self.posEntries.append(Entry(self.frame3,textvariable=self.posVal[i],justify="center",width=10))
            self.posEntries[i].grid(row=i+2,column=1)
            self.ampCheck.append(Checkbutton(self.frame3,variable=self.ampTick[i]))
            self.ampCheck[i].grid(row=i+2,column=2)
            self.ampEntries.append(Entry(self.frame3,textvariable=self.ampVal[i],justify="center",width=10))
            self.ampEntries[i].grid(row=i+2,column=3)
            self.widthCheck.append(Checkbutton(self.frame3,variable=self.widthTick[i]))
            self.widthCheck[i].grid(row=i+2,column=4)
            self.widthEntries.append(Entry(self.frame3,textvariable=self.widthVal[i],justify="center",width=10))
            self.widthEntries[i].grid(row=i+2,column=5)
            self.gaussCheck.append(Checkbutton(self.frame3,variable=self.gaussTick[i]))
            self.gaussCheck[i].grid(row=i+2,column=6)
            self.gaussEntries.append(Entry(self.frame3,textvariable=self.gaussVal[i],justify="center",width=10))
            self.gaussEntries[i].grid(row=i+2,column=7)
        self.reset()

    def reset(self):
        self.parent.pickNum=0
        self.bgrndVal.set("0.0")
        self.bgrndTick.set(1)
        self.slopeVal.set("0.0")
        self.slopeTick.set(1)
        self.numExp.set('1')
        self.pickTick.set(1)
        for i in range(10):
            self.posVal[i].set("0.0")
            self.posTick[i].set(0)
            self.ampVal[i].set("1.0")
            self.ampTick[i].set(0)
            self.widthVal[i].set("1.0")
            self.widthTick[i].set(0)
            self.gaussVal[i].set("0.0")
            self.gaussTick[i].set(1)
            if i > 0:
                self.posCheck[i].grid_remove()
                self.posEntries[i].grid_remove()
                self.ampCheck[i].grid_remove()
                self.ampEntries[i].grid_remove()
                self.widthCheck[i].grid_remove()
                self.widthEntries[i].grid_remove()
                self.gaussCheck[i].grid_remove()
                self.gaussEntries[i].grid_remove()
        self.togglePick()
        self.parent.pickWidth=False
        self.parent.showPlot()

    def changeNum(self,*args):
        val = int(self.numExp.get())
        for i in range(10):
            if i < val:
                self.posCheck[i].grid()
                self.posEntries[i].grid()
                self.ampCheck[i].grid()
                self.ampEntries[i].grid()
                self.widthCheck[i].grid()
                self.widthEntries[i].grid()
                self.gaussCheck[i].grid()
                self.gaussEntries[i].grid()
            else:
                self.posCheck[i].grid_remove()
                self.posEntries[i].grid_remove()
                self.ampCheck[i].grid_remove()
                self.ampEntries[i].grid_remove()
                self.widthCheck[i].grid_remove()
                self.widthEntries[i].grid_remove()
                self.gaussCheck[i].grid_remove()
                self.gaussEntries[i].grid_remove()

    def togglePick(self):
        self.parent.togglePick(self.pickTick.get())
                
    def fitFunc(self, x, *param):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        testFunc = np.zeros(len(self.parent.data1D))
        if struc[0]:
            bgrnd = param[0]
            param=np.delete(param,[0])
        else:
            bgrnd = argu[0]
            argu=np.delete(argu,[0])
        if struc[1]:
            slope = param[0]
            param=np.delete(param,[0])
        else:
            slope = argu[0]
            argu=np.delete(argu,[0])
        for i in range(numExp):
            if struc[4*i+2]:
                pos = param[0]
                param=np.delete(param,[0])
            else:
                pos= argu[0]
                argu=np.delete(argu,[0])
            if struc[4*i+3]:
                amp = param[0]
                param=np.delete(param,[0])
            else:
                amp= argu[0]
                argu=np.delete(argu,[0])
            if struc[4*i+4]:
                width = abs(param[0])
                param=np.delete(param,[0])
            else:
                width = argu[0]
                argu=np.delete(argu,[0])
            if struc[4*i+5]:
                gauss = abs(param[0])
                param=np.delete(param,[0])
            else:
                gauss = argu[0]
                argu=np.delete(argu,[0])
            t=np.arange(len(x))/self.parent.current.sw
            timeSignal = np.exp(1j*2*np.pi*t*((pos-self.axAdd)/self.axMult))/len(x)*np.exp(-np.pi*width*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
            testFunc += amp*np.real(np.fft.fftshift(np.fft.fft(timeSignal)))
            #testFunc += amp/np.pi*0.5*width/((x-(pos-self.axAdd)/self.axMult)**2+(0.5*width)**2)
        testFunc += bgrnd+slope*x
        return testFunc

    def fit(self,*args):
        struc = []
        guess = []
        argu = []
        numExp = int(self.numExp.get())
        outPos = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if self.bgrndTick.get() == 0:
            guess.append(safeEval(self.bgrndVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.bgrndVal.get())
            argu.append(inp)
            outBgrnd = inp
            self.bgrndVal.set('%.3g' % inp)
            struc.append(False)
        if self.slopeTick.get() == 0:
            guess.append(safeEval(self.slopeVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.slopeVal.get())
            argu.append(inp)
            outSlope = inp
            self.slopeVal.set('%.3g' % inp)
            struc.append(False)
        for i in range(numExp):
            if self.posTick[i].get() == 0:
                guess.append(safeEval(self.posVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.posVal[i].get())
                argu.append(inp)
                outPos[i] = inp
                self.posVal[i].set('%.3g' % inp)
                struc.append(False)
            if self.ampTick[i].get() == 0:
                guess.append(safeEval(self.ampVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.ampVal[i].get())
                argu.append(inp)
                outAmp[i] = inp
                self.ampVal[i].set('%.3g' % inp)
                struc.append(False)
            if self.widthTick[i].get() == 0:
                guess.append(abs(safeEval(self.widthVal[i].get())))
                struc.append(True)
            else:
                inp = abs(safeEval(self.widthVal[i].get()))
                argu.append(inp)
                outWidth[i] = inp
                self.widthVal[i].set('%.3g' % inp)
                struc.append(False)
            if self.gaussTick[i].get() == 0:
                guess.append(abs(safeEval(self.gaussVal[i].get())))
                struc.append(True)
            else:
                inp = abs(safeEval(self.gaussVal[i].get()))
                argu.append(inp)
                outGauss[i] = inp
                self.gaussVal[i].set('%.3g' % inp)
                struc.append(False)
        self.args = (numExp,struc,argu)
        fitVal = scipy.optimize.curve_fit(self.fitFunc,self.parent.xax,self.parent.data1D,p0=guess)
        counter = 0
        if struc[0]:
            self.bgrndVal.set('%.3g' % fitVal[0][counter])
            outBgrnd = fitVal[0][counter]
            counter +=1
        if struc[1]:
            self.slopeVal.set('%.3g' % fitVal[0][counter])
            outSlope = fitVal[0][counter]
            counter +=1
        for i in range(numExp):
            if struc[4*i+2]:
                self.posVal[i].set('%.3g' % fitVal[0][counter])
                outPos[i] = fitVal[0][counter]
                counter += 1
            if struc[4*i+3]:
                self.ampVal[i].set('%.3g' % fitVal[0][counter])
                outAmp[i] = fitVal[0][counter]
                counter += 1
            if struc[4*i+4]:
                self.widthVal[i].set('%.3g' % abs(fitVal[0][counter]))
                outWidth[i] = abs(fitVal[0][counter])
                counter += 1
            if struc[4*i+5]:
                self.gaussVal[i].set('%.3g' % abs(fitVal[0][counter]))
                outGauss[i] = abs(fitVal[0][counter])
                counter += 1
        self.disp(outBgrnd, outSlope, outAmp, outPos, outWidth, outGauss)

    def sim(self):
        numExp = int(self.numExp.get())
        outPos = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        outBgrnd = safeEval(self.bgrndVal.get())
        outSlope = safeEval(self.slopeVal.get())
        for i in range(numExp):
            outPos[i] = safeEval(self.posVal[i].get())
            outAmp[i] = safeEval(self.ampVal[i].get())
            outWidth[i] = abs(safeEval(self.widthVal[i].get()))
            self.widthVal[i].set('%.3g' % outWidth[i])
            outGauss[i] = abs(safeEval(self.gaussVal[i].get()))
            self.gaussVal[i].set('%.3g' % outGauss[i])
        self.disp(outBgrnd, outSlope, outAmp, outPos, outWidth, outGauss)

    def disp(self, outBgrnd, outSlope, outAmp, outPos, outWidth, outGauss):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x=[]
        t=np.arange(len(tmpx))/self.parent.current.sw
        for i in range(len(outAmp)):
            x.append(tmpx)
            timeSignal = np.exp(1j*2*np.pi*t*((outPos[i]-self.axAdd)/self.axMult))/len(tmpx)*np.exp(-np.pi*outWidth[i]*t)*np.exp(-((np.pi*outGauss[i]*t)**2)/(4*np.log(2)))
            y = outAmp[i]*np.real(np.fft.fftshift(np.fft.fft(timeSignal)))
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        self.parent.showPlot(tmpx, outCurve, x, outCurvePart)

##############################################################################
class TensorDeconvWindow(Frame): #a window for fitting relaxation data
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        self.fig = self.mainProgram.getFig()
        self.canvas = FigureCanvasTkAgg(self.fig, master=weakref.proxy(self))
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="nswe")
        self.current = TensorDeconvFrame(self,self.fig,self.canvas,oldMainWindow.current)
        self.paramframe = TensorDeconvParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        
    def rename(self,name):
        self.fig.suptitle(name)
        self.canvas.draw()
        self.oldMainWindow.rename(name)
        
    def buttonPress(self,event):
        self.current.buttonPress(event)

    def buttonRelease(self,event):
        self.current.buttonRelease(event)

    def pan(self,event):
        self.current.pan(event)

    def scroll(self,event):
        self.current.scroll(event)
        
    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class TensorDeconvFrame(Plot1DFrame): #a window for fitting relaxation data
    def __init__(self, rootwindow,fig,canvas,current):
        Plot1DFrame.__init__(self,rootwindow,fig,canvas)
        self.data1D = current.getDisplayedData()
        self.current = current
        self.spec = self.current.spec
        axAdd = 0
        if self.spec == 1:
            if self.current.ppm:
                axAdd = (self.current.freq-self.current.ref)/self.current.ref*1e6
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xax = self.current.xax*axMult+axAdd
        self.plotType=0
        self.rootwindow = rootwindow
        self.pickNum = 0
        self.pickNum2 = 0
        self.plotReset()
        self.showPlot()

    def plotReset(self): #set the plot limits to min and max values
        if self.plotType==0:
            miny = min(np.real(self.data1D))
            maxy = max(np.real(self.data1D))
        elif self.plotType==1:
            miny = min(np.imag(self.data1D))
            maxy = max(np.imag(self.data1D))
        elif self.plotType==2:
            miny = min(min(np.real(self.data1D)),min(np.imag(self.data1D)))
            maxy = max(max(np.real(self.data1D)),max(np.imag(self.data1D)))
        elif self.plotType==3:
            miny = min(np.abs(self.data1D))
            maxy = max(np.abs(self.data1D))
        else:
            miny=-1
            maxy=1
        differ = 0.05*(maxy-miny) #amount to add to show all datapoints (10%)
        self.yminlim=miny-differ
        self.ymaxlim=maxy+differ
        self.xminlim=min(self.xax)
        self.xmaxlim=max(self.xax)
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        self.ax.cla()
        self.line_xdata = self.xax
        self.line_ydata = self.data1D
        self.ax.plot(self.xax,self.data1D)
        self.ax.plot(tmpAx,tmpdata)
        for i in range(len(tmpAx2)):
            self.ax.plot(tmpAx2[i],tmpdata2[i])
        if self.spec==0:
            if self.current.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec==1:
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
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()

    def togglePick(self,var):
        if var==1:
            self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
            self.peakPick = True
        else:
            self.peakPickFunc = None
            self.peakPick = False
        
    def pickDeconv(self, pos):
        if self.pickNum2 == 0:
            if self.pickNum < 10:
                self.rootwindow.paramframe.numExp.set(str(self.pickNum+1))
                self.rootwindow.paramframe.changeNum()
            self.rootwindow.paramframe.t11Val[self.pickNum].set("%.2g" %self.current.xax[pos[0]])
            self.pickNum2 = 1
        elif self.pickNum2 == 1:
            self.rootwindow.paramframe.t22Val[self.pickNum].set("%.2g" %self.current.xax[pos[0]])
            self.pickNum2 = 2
        elif self.pickNum2 == 2:
            self.rootwindow.paramframe.t33Val[self.pickNum].set("%.2g" %self.current.xax[pos[0]])
            self.pickNum2 = 0
            self.pickNum += 1
        if self.pickNum < 10:
            self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
            self.peakPick = True 

#################################################################################
class TensorDeconvParamFrame(Frame): #a frame for the relaxtion parameters
    def __init__(self, parent, rootwindow): 
        self.parent = parent
        self.bgrndVal = StringVar()
        self.bgrndTick = IntVar()
        self.slopeVal = StringVar()
        self.slopeTick = IntVar()
        self.numExp = StringVar()
        self.chengVal = StringVar()
        self.pickTick = IntVar()
        self.pickTick.set(1)
        Frame.__init__(self, rootwindow)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0,sticky='n')
        self.optframe = Frame(self)
        self.optframe.grid(row=0,column=1,sticky='n')
        self.frame2 = Frame(self)
        self.frame2.grid(row=0,column=2,sticky='n')
        self.frame3 = Frame(self)
        self.frame3.grid(row=0,column=3,sticky='n')
        Button(self.frame1, text="Sim",command=self.sim).grid(row=0,column=0)
        Button(self.frame1, text="Fit",command=self.fit).grid(row=1,column=0)
        Button(self.frame1, text="Cancel",command=rootwindow.cancel).grid(row=2,column=0)
        Button(self.frame1, text="Reset",command=self.reset).grid(row=0,column=1)
        Checkbutton(self.frame1, variable=self.pickTick,text="Picking",command=self.togglePick).grid(row=1,column=1)
        Label(self.optframe,text="Cheng").grid(row=0,column=0)
        self.chengEntry = Entry(self.optframe,textvariable=self.chengVal,justify="center",width=10)
        self.chengEntry.grid(row=1,column=0)
        self.chengEntry.bind("<Return>", self.setCheng) 
        self.chengEntry.bind("<KP_Enter>", self.setCheng)
        Label(self.frame2,text="Bgrnd").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.bgrndTick).grid(row=1,column=0)
        Entry(self.frame2,textvariable=self.bgrndVal,justify="center",width=10).grid(row=1,column=1)
        Label(self.frame2,text="Slope").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.slopeTick).grid(row=3,column=0)
        Entry(self.frame2,textvariable=self.slopeVal,justify="center",width=10).grid(row=3,column=1)
        OptionMenu(self.frame3, self.numExp, "1","1", "2", "3","4","5","6","7","8","9","10",command=self.changeNum).grid(row=0,column=0,columnspan=6)
        Label(self.frame3,text="T11").grid(row=1,column=0,columnspan=2)
        Label(self.frame3,text="T22").grid(row=1,column=2,columnspan=2)
        Label(self.frame3,text="T33").grid(row=1,column=4,columnspan=2)
        Label(self.frame3,text="Amplitude").grid(row=1,column=6,columnspan=2)
        Label(self.frame3,text="Lorentz [Hz]").grid(row=1,column=8,columnspan=2)
        Label(self.frame3,text="Gauss [Hz]").grid(row=1,column=10,columnspan=2)
        self.t11Val = []
        self.t11Tick = []
        self.t22Val = []
        self.t22Tick = []
        self.t33Val = []
        self.t33Tick = []
        self.ampVal = []
        self.ampTick = []
        self.widthVal = []
        self.widthTick = []
        self.gaussVal = []
        self.gaussTick = []
        self.t11Check = []
        self.t11Entries = []
        self.t22Check = []
        self.t22Entries = []
        self.t33Check = []
        self.t33Entries = []
        self.ampCheck = []
        self.ampEntries = []
        self.widthCheck = []
        self.widthEntries = []
        self.gaussCheck = []
        self.gaussEntries = []

        for i in range(10):
            self.t11Val.append(StringVar())
            self.t11Tick.append(IntVar())
            self.t22Val.append(StringVar())
            self.t22Tick.append(IntVar())
            self.t33Val.append(StringVar())
            self.t33Tick.append(IntVar())
            self.ampVal.append(StringVar())
            self.ampTick.append(IntVar())
            self.widthVal.append(StringVar())
            self.widthTick.append(IntVar())
            self.gaussVal.append(StringVar())
            self.gaussTick.append(IntVar())
            self.t11Check.append(Checkbutton(self.frame3,variable=self.t11Tick[i]))
            self.t11Check[i].grid(row=i+2,column=0)
            self.t11Entries.append(Entry(self.frame3,textvariable=self.t11Val[i],justify="center",width=10))
            self.t11Entries[i].grid(row=i+2,column=1)
            self.t22Check.append(Checkbutton(self.frame3,variable=self.t22Tick[i]))
            self.t22Check[i].grid(row=i+2,column=2)
            self.t22Entries.append(Entry(self.frame3,textvariable=self.t22Val[i],justify="center",width=10))
            self.t22Entries[i].grid(row=i+2,column=3)
            self.t33Check.append(Checkbutton(self.frame3,variable=self.t33Tick[i]))
            self.t33Check[i].grid(row=i+2,column=4)
            self.t33Entries.append(Entry(self.frame3,textvariable=self.t33Val[i],justify="center",width=10))
            self.t33Entries[i].grid(row=i+2,column=5)
            self.ampCheck.append(Checkbutton(self.frame3,variable=self.ampTick[i]))
            self.ampCheck[i].grid(row=i+2,column=6)
            self.ampEntries.append(Entry(self.frame3,textvariable=self.ampVal[i],justify="center",width=10))
            self.ampEntries[i].grid(row=i+2,column=7)
            self.widthCheck.append(Checkbutton(self.frame3,variable=self.widthTick[i]))
            self.widthCheck[i].grid(row=i+2,column=8)
            self.widthEntries.append(Entry(self.frame3,textvariable=self.widthVal[i],justify="center",width=10))
            self.widthEntries[i].grid(row=i+2,column=9)
            self.gaussCheck.append(Checkbutton(self.frame3,variable=self.gaussTick[i]))
            self.gaussCheck[i].grid(row=i+2,column=10)
            self.gaussEntries.append(Entry(self.frame3,textvariable=self.gaussVal[i],justify="center",width=10))
            self.gaussEntries[i].grid(row=i+2,column=11)
        self.reset()

    def reset(self):
        self.parent.pickNum=0
        self.parent.pickNum2=0
        self.bgrndVal.set("0.0")
        self.bgrndTick.set(1)
        self.slopeVal.set("0.0")
        self.slopeTick.set(1)
        self.cheng = 15
        self.chengVal.set(str(self.cheng))
        self.pickTick.set(1)
        for i in range(10):
            self.t11Val[i].set("0.0")
            self.t11Tick[i].set(0)
            self.t22Val[i].set("0.0")
            self.t22Tick[i].set(0)
            self.t33Val[i].set("0.0")
            self.t33Tick[i].set(0)
            self.ampVal[i].set("1.0")
            self.ampTick[i].set(0)
            self.widthVal[i].set("10.0")
            self.widthTick[i].set(1)
            self.gaussVal[i].set("0.0")
            self.gaussTick[i].set(1)
            if i > 0:
                self.t11Check[i].grid_remove()
                self.t11Entries[i].grid_remove()
                self.t22Check[i].grid_remove()
                self.t22Entries[i].grid_remove()
                self.t33Check[i].grid_remove()
                self.t33Entries[i].grid_remove()
                self.ampCheck[i].grid_remove()
                self.ampEntries[i].grid_remove()
                self.widthCheck[i].grid_remove()
                self.widthEntries[i].grid_remove()
                self.gaussCheck[i].grid_remove()
                self.gaussEntries[i].grid_remove()
        self.togglePick()
        self.parent.showPlot()

    def togglePick(self):
        self.parent.togglePick(self.pickTick.get())
                
    def setCheng(self,*args):
        self.cheng = int(safeEval(self.chengVal.get()))
        self.chengVal.set(str(self.cheng))
                
    def changeNum(self,*args):
        val = int(self.numExp.get())
        for i in range(10):
            if i < val:
                self.t11Check[i].grid()
                self.t11Entries[i].grid()
                self.t22Check[i].grid()
                self.t22Entries[i].grid()
                self.t33Check[i].grid()
                self.t33Entries[i].grid()
                self.ampCheck[i].grid()
                self.ampEntries[i].grid()
                self.widthCheck[i].grid()
                self.widthEntries[i].grid()
                self.gaussCheck[i].grid()
                self.gaussEntries[i].grid()
            else:
                self.t11Check[i].grid_remove()
                self.t11Entries[i].grid_remove()
                self.t22Check[i].grid_remove()
                self.t22Entries[i].grid_remove()
                self.t33Check[i].grid_remove()
                self.t33Entries[i].grid_remove()
                self.ampCheck[i].grid_remove()
                self.ampEntries[i].grid_remove()
                self.widthCheck[i].grid_remove()
                self.widthEntries[i].grid_remove()
                self.gaussCheck[i].grid_remove()
                self.gaussEntries[i].grid_remove()

    def tensorFunc(self, x, t11, t22, t33, lor, gauss):
        t11=t11*self.multt11
        t22=t22*self.multt22
        t33=t33*self.multt33
        v=t11+t22+t33
        length =len(x)
        t=np.arange(length)/self.parent.current.sw
        final = np.zeros(length)
        mult=v/(self.parent.current.sw)*length
        x1=np.array(np.round(mult)+np.floor(length/2.0),dtype=int)
        weight = self.weight[np.logical_and(x1>=0,x1<length)]
        x1 = x1[np.logical_and(x1>=0,x1<length)]
        final = np.bincount(x1,weight,length)
        apod = np.exp(-np.pi*lor*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
        apod[-1:-(len(apod)/2+1):-1]=apod[:len(apod)/2]
        I=np.real(np.fft.fft(np.fft.ifft(final)*apod))
        return I
                
    def fitFunc(self, param, x, y):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        testFunc = np.zeros(len(self.parent.data1D))
        if struc[0]:
            bgrnd = param[0]
            param=np.delete(param,[0])
        else:
            bgrnd = argu[0]
            argu=np.delete(argu,[0])
        if struc[1]:
            slope = param[0]
            param=np.delete(param,[0])
        else:
            slope = argu[0]
            argu=np.delete(argu,[0])
        for i in range(numExp):
            if struc[6*i+2]:
                t11 = param[0]
                param=np.delete(param,[0])
            else:
                t11= argu[0]
                argu=np.delete(argu,[0])
            if struc[6*i+3]:
                t22 = param[0]
                param=np.delete(param,[0])
            else:
                t22= argu[0]
                argu=np.delete(argu,[0])
            if struc[6*i+4]:
                t33 = param[0]
                param=np.delete(param,[0])
            else:
                t33= argu[0]
                argu=np.delete(argu,[0])
            if struc[6*i+5]:
                amp = param[0]
                param=np.delete(param,[0])
            else:
                amp= argu[0]
                argu=np.delete(argu,[0])
            if struc[6*i+6]:
                width = abs(param[0])
                param=np.delete(param,[0])
            else:
                width = argu[0]
                argu=np.delete(argu,[0])
            if struc[6*i+7]:
                gauss = abs(param[0])
                param=np.delete(param,[0])
            else:
                gauss = argu[0]
                argu=np.delete(argu,[0])
            testFunc += amp*self.tensorFunc(x,t11,t22,t33,width,gauss)
        testFunc += bgrnd+slope*x
        return np.sum((np.real(testFunc)-y)**2)

    def disp(self,outBgrnd,outSlope,outt11,outt22,outt33,outAmp,outWidth,outGauss):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x=[]
        for i in range(len(outt11)):
            x.append(tmpx)
            y =  outAmp[i]*self.tensorFunc(tmpx,outt11[i],outt22[i],outt33[i],outWidth[i],outGauss[i])
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        self.parent.showPlot(tmpx, outCurve, x, outCurvePart)
    
    def fit(self,*args):
        self.setCheng()
        struc = []
        guess = []
        argu = []
        numExp = int(self.numExp.get())
        outt11 = np.zeros(numExp)
        outt22 = np.zeros(numExp)
        outt33 = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if self.bgrndTick.get() == 0:
            guess.append(safeEval(self.bgrndVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.bgrndVal.get())
            argu.append(inp)
            outBgrnd = inp
            self.bgrndVal.set('%.2g' % inp)
            struc.append(False)
        if self.slopeTick.get() == 0:
            guess.append(safeEval(self.slopeVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.slopeVal.get())
            argu.append(inp)
            outSlope = inp
            self.slopeVal.set('%.2g' % inp)
            struc.append(False)
        for i in range(numExp):
            if self.t11Tick[i].get() == 0:
                guess.append(safeEval(self.t11Val[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.t11Val[i].get())
                argu.append(inp)
                outt11[i] = inp
                self.t11Val[i].set('%.2g' % inp)
                struc.append(False)
            if self.t22Tick[i].get() == 0:
                guess.append(safeEval(self.t22Val[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.t22Val[i].get())
                argu.append(inp)
                outt22[i] = inp
                self.t22Val[i].set('%.2g' % inp)
                struc.append(False)
            if self.t33Tick[i].get() == 0:
                guess.append(safeEval(self.t33Val[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.t33Val[i].get())
                argu.append(inp)
                outt33[i] = inp
                self.t33Val[i].set('%.2g' % inp)
                struc.append(False)
            if self.ampTick[i].get() == 0:
                guess.append(safeEval(self.ampVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.ampVal[i].get())
                argu.append(inp)
                outAmp[i] = inp
                self.ampVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.widthTick[i].get() == 0:
                guess.append(abs(safeEval(self.widthVal[i].get())))
                struc.append(True)
            else:
                inp = abs(safeEval(self.widthVal[i].get()))
                argu.append(inp)
                outWidth[i] = inp
                self.widthVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.gaussTick[i].get() == 0:
                guess.append(abs(safeEval(self.gaussVal[i].get())))
                struc.append(True)
            else:
                inp = abs(safeEval(self.gaussVal[i].get()))
                argu.append(inp)
                outGauss[i] = inp
                self.gaussVal[i].set('%.2g' % inp)
                struc.append(False)
        self.args = (numExp,struc,argu)
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.multt11=np.sin(theta)**2*np.cos(phi)**2
        self.multt22=np.sin(theta)**2*np.sin(phi)**2
        self.multt33=np.cos(theta)**2
        fitVal = scipy.optimize.fmin(self.fitFunc,guess, args=(self.parent.xax,np.real(self.parent.data1D)))
        counter = 0
        if struc[0]:
            self.bgrndVal.set('%.2g' % fitVal[counter])
            outBgrnd = fitVal[counter]
            counter +=1
        if struc[1]:
            self.slopeVal.set('%.2g' % fitVal[counter])
            outSlope = fitVal[counter]
            counter +=1
        for i in range(numExp):
            if struc[6*i+2]:
                self.t11Val[i].set('%.2g' % fitVal[counter])
                outt11[i] = fitVal[counter]
                counter += 1
            if struc[6*i+3]:
                self.t22Val[i].set('%.2g' % fitVal[counter])
                outt22[i] = fitVal[counter]
                counter += 1
            if struc[6*i+4]:
                self.t33Val[i].set('%.2g' % fitVal[counter])
                outt33[i] = fitVal[counter]
                counter += 1
            if struc[6*i+5]:
                self.ampVal[i].set('%.2g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[6*i+6]:
                self.widthVal[i].set('%.2g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
            if struc[6*i+7]:
                self.gaussVal[i].set('%.2g' % abs(fitVal[counter]))
                outGauss[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd,outSlope,outt11,outt22,outt33,outAmp,outWidth,outGauss)

    def sim(self):
        self.setCheng()
        numExp = int(self.numExp.get())
        bgrnd = safeEval(self.bgrndVal.get())
        slope = safeEval(self.slopeVal.get())
        t11 = np.zeros(numExp)
        t22 = np.zeros(numExp)
        t33 = np.zeros(numExp)
        amp = np.zeros(numExp)
        width = np.zeros(numExp)
        gauss = np.zeros(numExp)
        for i in range(numExp):
            t11[i] = safeEval(self.t11Val[i].get())
            t22[i] = safeEval(self.t22Val[i].get())
            t33[i] = safeEval(self.t33Val[i].get())
            amp[i] = safeEval(self.ampVal[i].get())
            width[i] = safeEval(self.widthVal[i].get())
            gauss[i] = safeEval(self.gaussVal[i].get())
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.multt11=np.sin(theta)**2*np.cos(phi)**2
        self.multt22=np.sin(theta)**2*np.sin(phi)**2
        self.multt33=np.cos(theta)**2
        self.disp(bgrnd,slope,t11,t22,t33,amp,width,gauss)
        
##############################################################################
class Quad1DeconvWindow(Frame): #a window for fitting relaxation data
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        self.fig = self.mainProgram.getFig()
        self.canvas = FigureCanvasTkAgg(self.fig, master=weakref.proxy(self))
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="nswe")
        self.current = Quad1DeconvFrame(self,self.fig,self.canvas,oldMainWindow.current)
        self.paramframe = Quad1DeconvParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)

    def rename(self,name):
        self.fig.suptitle(name)
        self.canvas.draw()
        self.oldMainWindow.rename(name)

    def buttonPress(self,event):
        self.current.buttonPress(event)

    def buttonRelease(self,event):
        self.current.buttonRelease(event)

    def pan(self,event):
        self.current.pan(event)

    def scroll(self,event):
        self.current.scroll(event)
        
    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class Quad1DeconvFrame(Plot1DFrame): #a window for fitting relaxation data
    def __init__(self, rootwindow,fig,canvas,current):
        Plot1DFrame.__init__(self,rootwindow,fig,canvas)
        self.data1D = current.getDisplayedData()
        self.current = current
        self.spec = self.current.spec
        self.xax = self.current.xax
        self.plotType=0
        self.rootwindow = rootwindow
        self.pickNum = 0
        self.pickNum2 = 0
        self.plotReset()
        self.showPlot()

    def plotReset(self): #set the plot limits to min and max values
        if self.plotType==0:
            miny = min(np.real(self.data1D))
            maxy = max(np.real(self.data1D))
        elif self.plotType==1:
            miny = min(np.imag(self.data1D))
            maxy = max(np.imag(self.data1D))
        elif self.plotType==2:
            miny = min(min(np.real(self.data1D)),min(np.imag(self.data1D)))
            maxy = max(max(np.real(self.data1D)),max(np.imag(self.data1D)))
        elif self.plotType==3:
            miny = min(np.abs(self.data1D))
            maxy = max(np.abs(self.data1D))
        else:
            miny=-1
            maxy=1
        differ = 0.05*(maxy-miny) #amount to add to show all datapoints (10%)
        self.yminlim=miny-differ
        self.ymaxlim=maxy+differ
        axAdd = 0
        if self.spec == 1:
            if self.current.ppm:
                axAdd = (self.current.freq-self.current.ref)/self.current.ref*1e6
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xminlim=min(self.xax*axMult+axAdd)
        self.xmaxlim=max(self.xax*axMult+axAdd)
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        self.ax.cla()
        axAdd = 0
        if self.spec == 1:
            if self.current.ppm:
                axAdd = (self.current.freq-self.current.ref)/self.current.ref*1e6
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.line_xdata = self.xax*axMult+axAdd
        self.line_ydata = self.data1D
        self.ax.plot(self.xax*axMult+axAdd,self.data1D)
        if tmpAx is not None:
            self.ax.plot(tmpAx*axMult+axAdd,tmpdata)
        for i in range(len(tmpAx2)):
            self.ax.plot(tmpAx2[i]*axMult+axAdd,tmpdata2[i])
        if self.spec==0:
            if self.current.axType == 0:
                self.ax.set_xlabel('Time [s]')
            elif self.current.axType == 1:
                self.ax.set_xlabel('Time [ms]')
            elif self.current.axType == 2:
                self.ax.set_xlabel(r'Time [$\mu$s]')
            else:
                self.ax.set_xlabel('User defined')
        elif self.spec==1:
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
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()

#################################################################################
class Quad1DeconvParamFrame(Frame): #a frame for the quadrupole parameters
    
    Ioptions = ['1','3/2','2','5/2','3','7/2','4','9/2']
    
    def __init__(self, parent, rootwindow):
        self.parent = parent
        self.bgrndVal = StringVar()
        self.bgrndVal.set("0.0")
        self.bgrndTick = IntVar()
        self.bgrndTick.set(1)
        self.slopeVal = StringVar()
        self.slopeVal.set("0.0")
        self.slopeTick = IntVar()
        self.slopeTick.set(1)
        self.numExp = StringVar()
        self.cheng = 15
        self.chengVal = StringVar()
        self.chengVal.set(str(self.cheng))
        self.IVal = StringVar()
        self.IVal.set("3/2")
        Frame.__init__(self, rootwindow)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0,sticky='n')
        self.optframe = Frame(self)
        self.optframe.grid(row=0,column=1,sticky='n')
        self.frame2 = Frame(self)
        self.frame2.grid(row=0,column=2,sticky='n')
        self.frame3 = Frame(self)
        self.frame3.grid(row=0,column=3,sticky='n')
        Button(self.frame1, text="Sim",command=self.sim).grid(row=0,column=0)
        Button(self.frame1, text="Fit",command=self.fit).grid(row=1,column=0)
        Button(self.frame1, text="Cancel",command=rootwindow.cancel).grid(row=2,column=0)
        Label(self.optframe,text="Cheng").grid(row=0,column=0)
        self.chengEntry = Entry(self.optframe,textvariable=self.chengVal,justify="center",width=10)
        self.chengEntry.grid(row=1,column=0)
        self.chengEntry.bind("<Return>", self.setCheng) 
        self.chengEntry.bind("<KP_Enter>", self.setCheng)
        Label(self.optframe,text="I").grid(row=0,column=1)
        OptionMenu(self.optframe,self.IVal,self.IVal.get(),*self.Ioptions).grid(row=1,column=1)
        Label(self.frame2,text="Bgrnd").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.bgrndTick).grid(row=1,column=0)
        Entry(self.frame2,textvariable=self.bgrndVal,justify="center",width=10).grid(row=1,column=1)
        Label(self.frame2,text="Slope").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.slopeTick).grid(row=3,column=0)
        Entry(self.frame2,textvariable=self.slopeVal,justify="center",width=10).grid(row=3,column=1)
        OptionMenu(self.frame3, self.numExp, "1","1", "2", "3","4","5","6","7","8","9","10",command=self.changeNum).grid(row=0,column=0,columnspan=6)
        Label(self.frame3,text="Pos").grid(row=1,column=0,columnspan=2)
        Label(self.frame3,text="Cq [MHz]").grid(row=1,column=2,columnspan=2)
        Label(self.frame3,text="Eta").grid(row=1,column=4,columnspan=2)
        Label(self.frame3,text="Amplitude").grid(row=1,column=6,columnspan=2)
        Label(self.frame3,text="Lorentz [Hz]").grid(row=1,column=8,columnspan=2)
        Label(self.frame3,text="Gauss [Hz]").grid(row=1,column=10,columnspan=2)
        self.posVal = []
        self.posTick = []
        self.cqVal = []
        self.cqTick = []
        self.etaVal = []
        self.etaTick = []
        self.ampVal = []
        self.ampTick = []
        self.widthVal = []
        self.widthTick = []
        self.gaussVal = []
        self.gaussTick = []
        self.posCheck = []
        self.posEntries = []
        self.cqCheck = []
        self.cqEntries = []
        self.etaCheck = []
        self.etaEntries = []
        self.ampCheck = []
        self.ampEntries = []
        self.widthCheck = []
        self.widthEntries = []
        self.gaussCheck = []
        self.gaussEntries = []
        for i in range(10):
            self.posVal.append(StringVar())
            self.posVal[i].set("0.0")
            self.posTick.append(IntVar())
            self.cqVal.append(StringVar())
            self.cqVal[i].set("0.0")
            self.cqTick.append(IntVar())
            self.etaVal.append(StringVar())
            self.etaVal[i].set("0.0")
            self.etaTick.append(IntVar())
            self.ampVal.append(StringVar())
            self.ampVal[i].set("1.0")
            self.ampTick.append(IntVar())
            self.widthVal.append(StringVar())
            self.widthVal[i].set("10.0")
            self.widthTick.append(IntVar())
            self.widthTick[i].set(1)
            self.gaussVal.append(StringVar())
            self.gaussVal[i].set("0.0")
            self.gaussTick.append(IntVar())
            self.gaussTick[i].set(1)
            self.posCheck.append(Checkbutton(self.frame3,variable=self.posTick[i]))
            self.posCheck[i].grid(row=i+2,column=0)
            self.posEntries.append(Entry(self.frame3,textvariable=self.posVal[i],justify="center",width=10))
            self.posEntries[i].grid(row=i+2,column=1)
            self.cqCheck.append(Checkbutton(self.frame3,variable=self.cqTick[i]))
            self.cqCheck[i].grid(row=i+2,column=2)
            self.cqEntries.append(Entry(self.frame3,textvariable=self.cqVal[i],justify="center",width=10))
            self.cqEntries[i].grid(row=i+2,column=3)
            self.etaCheck.append(Checkbutton(self.frame3,variable=self.etaTick[i]))
            self.etaCheck[i].grid(row=i+2,column=4)
            self.etaEntries.append(Entry(self.frame3,textvariable=self.etaVal[i],justify="center",width=10))
            self.etaEntries[i].grid(row=i+2,column=5)
            self.ampCheck.append(Checkbutton(self.frame3,variable=self.ampTick[i]))
            self.ampCheck[i].grid(row=i+2,column=6)
            self.ampEntries.append(Entry(self.frame3,textvariable=self.ampVal[i],justify="center",width=10))
            self.ampEntries[i].grid(row=i+2,column=7)
            self.widthCheck.append(Checkbutton(self.frame3,variable=self.widthTick[i]))
            self.widthCheck[i].grid(row=i+2,column=8)
            self.widthEntries.append(Entry(self.frame3,textvariable=self.widthVal[i],justify="center",width=10))
            self.widthEntries[i].grid(row=i+2,column=9)
            self.gaussCheck.append(Checkbutton(self.frame3,variable=self.gaussTick[i]))
            self.gaussCheck[i].grid(row=i+2,column=10)
            self.gaussEntries.append(Entry(self.frame3,textvariable=self.gaussVal[i],justify="center",width=10))
            self.gaussEntries[i].grid(row=i+2,column=11)
            if i > 0:
                self.posCheck[i].grid_remove()
                self.posEntries[i].grid_remove()
                self.cqCheck[i].grid_remove()
                self.cqEntries[i].grid_remove()
                self.etaCheck[i].grid_remove()
                self.etaEntries[i].grid_remove()
                self.ampCheck[i].grid_remove()
                self.ampEntries[i].grid_remove()
                self.widthCheck[i].grid_remove()
                self.widthEntries[i].grid_remove()
                self.gaussCheck[i].grid_remove()
                self.gaussEntries[i].grid_remove()

    def checkI(self,I):
        if I == '1':
            return 1.0
        elif I == '3/2':
            return 1.5
        elif I == '2':
            return 2
        elif I == '5/2':
            return 2.5
        elif I == '3':
            return 3
        elif I == '7/2':
            return 3.5
        elif I == '4':
            return 4
        elif I == '9/2':
            return 4.5
                
    def setCheng(self,*args):
        self.cheng = int(safeEval(self.chengVal.get()))
        self.chengVal.set(str(self.cheng))
                
    def changeNum(self,*args):
        val = int(self.numExp.get())
        for i in range(10):
            if i < val:
                self.posCheck[i].grid()
                self.posEntries[i].grid()
                self.cqCheck[i].grid()
                self.cqEntries[i].grid()
                self.etaCheck[i].grid()
                self.etaEntries[i].grid()
                self.ampCheck[i].grid()
                self.ampEntries[i].grid()
                self.widthCheck[i].grid()
                self.widthEntries[i].grid()
                self.gaussCheck[i].grid()
                self.gaussEntries[i].grid()
            else:
                self.posCheck[i].grid_remove()
                self.posEntries[i].grid_remove()
                self.cqCheck[i].grid_remove()
                self.cqEntries[i].grid_remove()
                self.etaCheck[i].grid_remove()
                self.etaEntries[i].grid_remove()
                self.ampCheck[i].grid_remove()
                self.ampEntries[i].grid_remove()
                self.widthCheck[i].grid_remove()
                self.widthEntries[i].grid_remove()
                self.gaussCheck[i].grid_remove()
                self.gaussEntries[i].grid_remove()

    def tensorFunc(self, x, I, pos, cq, eta, width, gauss):
        m=np.arange(-I,I)
        v=[]
        weights=[]
        for i in m:
            tmp = (cq/(4*I*(2*I-1))*(I*(I+1)-3*(i+1)**2))-(cq/(4*I*(2*I-1))*(I*(I+1)-3*(i)**2))
            v=np.append(v,tmp*(self.angleStuff1-eta*self.angleStuff2)+pos)
            weights=np.append(weights,self.weight)
        length =len(x)
        t=np.arange(length)/self.parent.current.sw
        final = np.zeros(length)
        mult=v/(self.parent.current.sw)*length
        x1=np.array(np.round(mult)+np.floor(length/2),dtype=int)
        weights = weights[np.logical_and(x1>=0,x1<length)]
        x1 = x1[np.logical_and(x1>=0,x1<length)]
        final = np.bincount(x1,weights,length)
        apod = np.exp(-np.pi*width*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
        apod[-1:-(len(apod)/2+1):-1]=apod[:len(apod)/2]
        inten=np.real(np.fft.fft(np.fft.ifft(final)*apod))
        return inten
                
    def fitFunc(self, param, x, y):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        I = self.args[3]
        testFunc = np.zeros(len(self.parent.data1D))
        if struc[0]:
            bgrnd = param[0]
            param=np.delete(param,[0])
        else:
            bgrnd = argu[0]
            argu=np.delete(argu,[0])
        if struc[1]:
            slope = param[0]
            param=np.delete(param,[0])
        else:
            slope = argu[0]
            argu=np.delete(argu,[0])
        for i in range(numExp):
            if struc[6*i+2]:
                pos = param[0]
                param=np.delete(param,[0])
            else:
                pos= argu[0]
                argu=np.delete(argu,[0])
            if struc[6*i+3]:
                cq = param[0]
                param=np.delete(param,[0])
            else:
                cq= argu[0]
                argu=np.delete(argu,[0])
            if struc[6*i+4]:
                eta = param[0]
                param=np.delete(param,[0])
            else:
                eta= argu[0]
                argu=np.delete(argu,[0])
            if struc[6*i+5]:
                amp = param[0]
                param=np.delete(param,[0])
            else:
                amp= argu[0]
                argu=np.delete(argu,[0])
            if struc[6*i+6]:
                width = abs(param[0])
                param=np.delete(param,[0])
            else:
                width = argu[0]
                argu=np.delete(argu,[0])
            if struc[6*i+7]:
                gauss = abs(param[0])
                param=np.delete(param,[0])
            else:
                gauss = argu[0]
                argu=np.delete(argu,[0])
            testFunc += amp*self.tensorFunc(x,I,pos,cq,eta,width,gauss)
        testFunc += bgrnd+slope*x
        return np.sum((np.real(testFunc)-y)**2)

    def disp(self,outBgrnd,outSlope,outI,outPos,outCq,outEta,outAmp,outWidth,outGauss):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x=[]
        for i in range(len(outPos)):
            x.append(tmpx)
            y =  outAmp[i]*self.tensorFunc(tmpx,outI,outPos[i],outCq[i],outEta[i],outWidth[i],outGauss[i])
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        self.parent.showPlot(tmpx, outCurve, x, outCurvePart)

    def setAngleStuff(self):
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.angleStuff1 = 0.5*(3*np.cos(theta)**2-1)
        self.angleStuff2 = np.cos(2*phi)*(np.cos(theta)**2-1)
        
    def fit(self,*args):
        self.setCheng()
        struc = []
        guess = []
        argu = []
        numExp = int(self.numExp.get())
        outPos = np.zeros(numExp)
        outCq = np.zeros(numExp)
        outEta = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        I = self.checkI(self.IVal.get())
        if self.bgrndTick.get() == 0:
            guess.append(safeEval(self.bgrndVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.bgrndVal.get())
            argu.append(inp)
            outBgrnd = inp
            self.bgrndVal.set('%.2g' % inp)
            struc.append(False)
        if self.slopeTick.get() == 0:
            guess.append(safeEval(self.slopeVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.slopeVal.get())
            argu.append(inp)
            outSlope = inp
            self.slopeVal.set('%.2g' % inp)
            struc.append(False)
        for i in range(numExp):
            if self.posTick[i].get() == 0:
                guess.append(safeEval(self.posVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.posVal[i].get())
                argu.append(inp)
                outPos[i] = inp
                self.posVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.cqTick[i].get() == 0:
                guess.append(safeEval(self.cqVal[i].get())*1e6)
                struc.append(True)
            else:
                inp = safeEval(self.cqVal[i].get())
                argu.append(inp*1e6)
                outCq[i] = inp*1e6
                self.cqVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.etaTick[i].get() == 0:
                guess.append(safeEval(self.etaVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.etaVal[i].get())
                argu.append(inp)
                outEta[i] = inp
                self.etaVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.ampTick[i].get() == 0:
                guess.append(safeEval(self.ampVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.ampVal[i].get())
                argu.append(inp)
                outAmp[i] = inp
                self.ampVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.widthTick[i].get() == 0:
                guess.append(abs(safeEval(self.widthVal[i].get())))
                struc.append(True)
            else:
                inp = abs(safeEval(self.widthVal[i].get()))
                argu.append(inp)
                outWidth[i] = inp
                self.widthVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.gaussTick[i].get() == 0:
                guess.append(abs(safeEval(self.gaussVal[i].get())))
                struc.append(True)
            else:
                inp = abs(safeEval(self.gaussVal[i].get()))
                argu.append(inp)
                outGauss[i] = inp
                self.gaussVal[i].set('%.2g' % inp)
                struc.append(False)
        self.args = (numExp,struc,argu,I)
        self.setAngleStuff()
        fitVal = scipy.optimize.fmin(self.fitFunc,guess, args=(self.parent.xax,np.real(self.parent.data1D)))
        counter = 0
        if struc[0]:
            self.bgrndVal.set('%.2g' % fitVal[counter])
            outBgrnd = fitVal[counter]
            counter +=1
        if struc[1]:
            self.slopeVal.set('%.2g' % fitVal[counter])
            outSlope = fitVal[counter]
            counter +=1
        for i in range(numExp):
            if struc[6*i+2]:
                self.posVal[i].set('%.2g' % fitVal[counter])
                outPos[i] = fitVal[counter]
                counter += 1
            if struc[6*i+3]:
                self.cqVal[i].set('%.2g' % (fitVal[counter]*1e-6))
                outCq[i] = fitVal[counter]
                counter += 1
            if struc[6*i+4]:
                self.etaVal[i].set('%.2g' % fitVal[counter])
                outEta[i] = fitVal[counter]
                counter += 1
            if struc[6*i+5]:
                self.ampVal[i].set('%.2g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[6*i+6]:
                self.widthVal[i].set('%.2g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
            if struc[6*i+7]:
                self.gaussVal[i].set('%.2g' % abs(fitVal[counter]))
                outGauss[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd,outSlope,I,outPos,outCq,outEta,outAmp,outWidth,outGauss)

    def sim(self):
        self.setCheng()
        numExp = int(self.numExp.get())
        bgrnd = safeEval(self.bgrndVal.get())
        slope = safeEval(self.slopeVal.get())
        pos = np.zeros(numExp)
        cq = np.zeros(numExp)
        eta = np.zeros(numExp)
        amp = np.zeros(numExp)
        width = np.zeros(numExp)
        gauss = np.zeros(numExp)
        I = self.checkI(self.IVal.get())
        for i in range(numExp):
            pos[i] = safeEval(self.posVal[i].get())
            cq[i] = safeEval(self.cqVal[i].get())*1e6
            eta[i] = safeEval(self.etaVal[i].get())
            amp[i] = safeEval(self.ampVal[i].get())
            width[i] = safeEval(self.widthVal[i].get())
            gauss[i] = safeEval(self.gaussVal[i].get())
        self.setAngleStuff()
        self.disp(bgrnd,slope,I,pos,cq,eta,amp,width,gauss)

##############################################################################
class Quad2DeconvWindow(Frame): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow,mas=False):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        self.fig = self.mainProgram.getFig()
        self.canvas = FigureCanvasTkAgg(self.fig, master=weakref.proxy(self))
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="nswe")
        self.current = Quad1DeconvFrame(self,self.fig,self.canvas,oldMainWindow.current)
        if mas:
            self.paramframe = Quad2MASDeconvParamFrame(self.current,self)
        else:
            self.paramframe = Quad2StaticDeconvParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        
    def rename(self,name):
        self.fig.suptitle(name)
        self.canvas.draw()
        self.oldMainWindow.rename(name)
        
    def buttonPress(self,event):
        self.current.buttonPress(event)

    def buttonRelease(self,event):
        self.current.buttonRelease(event)

    def pan(self,event):
        self.current.pan(event)

    def scroll(self,event):
        self.current.scroll(event)
        
    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)
         
#################################################################################
class Quad2StaticDeconvParamFrame(Quad1DeconvParamFrame): 

    Ioptions = ['3/2','5/2','7/2','9/2']
    
    def __init__(self, parent, rootwindow):
        Quad1DeconvParamFrame.__init__(self,parent,rootwindow)
        
    def setAngleStuff(self):
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.angleStuff1 = -27/8.0*np.cos(theta)**4+15/4.0*np.cos(theta)**2-3/8.0
        self.angleStuff2 = (-9/4.0*np.cos(theta)**4+2*np.cos(theta)**2+1/4.0)*np.cos(2*phi)
        self.angleStuff3 = -1/2.0*np.cos(theta)**2+1/3.0+(-3/8.0*np.cos(theta)**4+3/4.0*np.cos(theta)**2-3/8.0)*np.cos(2*phi)**2

    def tensorFunc(self, x, I, pos, cq, eta, width, gauss):
        v = -cq**2/(6.0*self.parent.current.freq)*(I*(I+1)-3/4.0)*(self.angleStuff1+self.angleStuff2*eta+self.angleStuff3*eta**2)+pos
        length =len(x)
        t=np.arange(length)/self.parent.current.sw
        final = np.zeros(length)
        mult=v/(self.parent.current.sw)*length
        x1=np.array(np.round(mult)+np.floor(length/2),dtype=int)
        weights = self.weight[np.logical_and(x1>=0,x1<length)]
        x1 = x1[np.logical_and(x1>=0,x1<length)]
        final = np.bincount(x1,weights,length)
        apod = np.exp(-np.pi*width*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
        apod[-1:-(len(apod)/2+1):-1]=apod[:len(apod)/2]
        inten=np.real(np.fft.fft(np.fft.ifft(final)*apod))
        return inten
    
#################################################################################
class Quad2MASDeconvParamFrame(Quad2StaticDeconvParamFrame): 
    Ioptions = ['3/2','5/2','7/2','9/2']
    
    def __init__(self, parent, rootwindow):
        Quad2StaticDeconvParamFrame.__init__(self,parent,rootwindow)
        
    def setAngleStuff(self):
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.angleStuff1 = 21/16.0*np.cos(theta)**4-9/8.0*np.cos(theta)**2+5/16.0
        self.angleStuff2 = (-7/8.0*np.cos(theta)**4+np.cos(theta)**2-1/8.0)*np.cos(2*phi)
        self.angleStuff3 = 1/12.0*np.cos(theta)**2+(+7/48.0*np.cos(theta)**4-7/24.0*np.cos(theta)**2+7/48.0)*np.cos(2*phi)**2

##############################################################################
class Quad2CzjzekWindow(Frame): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow,mas=False):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        self.fig = self.mainProgram.getFig()
        self.canvas = FigureCanvasTkAgg(self.fig, master=weakref.proxy(self))
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="nswe")
        self.current = Quad1DeconvFrame(self,self.fig,self.canvas,oldMainWindow.current)
        if mas:
            self.paramframe = Quad2MASCzjzekParamFrame(self.current,self)
        else:
            self.paramframe = Quad2StaticCzjzekParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        
    def rename(self,name):
        self.fig.suptitle(name)
        self.canvas.draw()
        self.oldMainWindow.rename(name)
        
    def buttonPress(self,event):
        self.current.buttonPress(event)

    def buttonRelease(self,event):
        self.current.buttonRelease(event)

    def pan(self,event):
        self.current.pan(event)

    def scroll(self,event):
        self.current.scroll(event)
        
    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)

#################################################################################
class Quad2StaticCzjzekParamFrame(Frame): 

    Ioptions = ['3/2','5/2','7/2','9/2']

    def __init__(self, parent, rootwindow):
        self.parent = parent
        self.bgrndVal = StringVar()
        self.bgrndVal.set("0.0")
        self.bgrndTick = IntVar()
        self.bgrndTick.set(1)
        self.slopeVal = StringVar()
        self.slopeVal.set("0.0")
        self.slopeTick = IntVar()
        self.slopeTick.set(1)
        self.numExp = StringVar()
        self.cheng = 15
        self.chengVal = StringVar()
        self.chengVal.set(str(self.cheng))
        self.cqGridVal = StringVar()
        self.cqGridVal.set('50')
        self.etaGridVal = StringVar()
        self.etaGridVal.set('10')
        self.IVal = StringVar()
        self.IVal.set('3/2')
        Frame.__init__(self, rootwindow)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0,sticky='n')
        self.optframe = Frame(self)
        self.optframe.grid(row=0,column=1,sticky='n')
        self.frame2 = Frame(self)
        self.frame2.grid(row=0,column=2,sticky='n')
        self.frame3 = Frame(self)
        self.frame3.grid(row=0,column=3,sticky='n')
        Button(self.frame1, text="Sim",command=self.sim).grid(row=0,column=0)
        Button(self.frame1, text="Fit",command=self.fit).grid(row=1,column=0)
        Button(self.frame1, text="Cancel",command=rootwindow.cancel).grid(row=2,column=0)
        Label(self.optframe,text="Cheng").grid(row=0,column=0)
        self.chengEntry = Entry(self.optframe,textvariable=self.chengVal,justify="center",width=10)
        self.chengEntry.grid(row=1,column=0)
        self.chengEntry.bind("<Return>", self.setCheng) 
        self.chengEntry.bind("<KP_Enter>", self.setCheng)
        Label(self.optframe,text="Cq grid size").grid(row=2,column=0)
        self.cqGridEntry = Entry(self.optframe,textvariable=self.cqGridVal,justify="center",width=10)
        self.cqGridEntry.grid(row=3,column=0)
        self.cqGridEntry.bind("<Return>", self.setGrid) 
        self.cqGridEntry.bind("<KP_Enter>", self.setGrid)
        Label(self.optframe,text=u"\u03b7 grid size").grid(row=4,column=0)
        self.etaGridEntry = Entry(self.optframe,textvariable=self.etaGridVal,justify="center",width=10)
        self.etaGridEntry.grid(row=5,column=0)
        self.etaGridEntry.bind("<Return>", self.setGrid) 
        self.etaGridEntry.bind("<KP_Enter>", self.setGrid)
        Label(self.optframe,text="I").grid(row=0,column=1)
        OptionMenu(self.optframe,self.IVal,self.IVal.get(),*self.Ioptions).grid(row=1,column=1)
        Label(self.frame2,text="Bgrnd").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.bgrndTick).grid(row=1,column=0)
        Entry(self.frame2,textvariable=self.bgrndVal,justify="center",width=10).grid(row=1,column=1)
        Label(self.frame2,text="Slope").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.slopeTick).grid(row=3,column=0)
        Entry(self.frame2,textvariable=self.slopeVal,justify="center",width=10).grid(row=3,column=1)
        OptionMenu(self.frame3, self.numExp, "1","1", "2", "3","4","5","6","7","8","9","10",command=self.changeNum).grid(row=0,column=0,columnspan=6)
        Label(self.frame3,text="d").grid(row=1,column=0)
        Label(self.frame3,text="Pos").grid(row=1,column=1,columnspan=2)
        Label(self.frame3,text=u"\u03c3 [MHz]").grid(row=1,column=3,columnspan=2)
        Label(self.frame3,text="Amplitude").grid(row=1,column=5,columnspan=2)
        Label(self.frame3,text="Lorentz [Hz]").grid(row=1,column=7,columnspan=2)
        Label(self.frame3,text="Gauss [Hz]").grid(row=1,column=9,columnspan=2)
        self.posVal = []
        self.posTick = []
        self.sigmaVal = []
        self.sigmaTick = []
        self.dVal = []
        self.ampVal = []
        self.ampTick = []
        self.widthVal = []
        self.widthTick = []
        self.gaussVal = []
        self.gaussTick = []
        self.posCheck = []
        self.posEntries = []
        self.sigmaCheck = []
        self.sigmaEntries = []
        self.dEntries = []
        self.ampCheck = []
        self.ampEntries = []
        self.widthCheck = []
        self.widthEntries = []
        self.gaussCheck = []
        self.gaussEntries = []
        for i in range(10):
            self.posVal.append(StringVar())
            self.posVal[i].set("0.0")
            self.posTick.append(IntVar())
            self.sigmaVal.append(StringVar())
            self.sigmaVal[i].set("1.0")
            self.sigmaTick.append(IntVar())
            self.dVal.append(StringVar())
            self.dVal[i].set("5")
            self.ampVal.append(StringVar())
            self.ampVal[i].set("1.0")
            self.ampTick.append(IntVar())
            self.widthVal.append(StringVar())
            self.widthVal[i].set("10.0")
            self.widthTick.append(IntVar())
            self.widthTick[i].set(1)
            self.gaussVal.append(StringVar())
            self.gaussVal[i].set("0.0")
            self.gaussTick.append(IntVar())
            self.gaussTick[i].set(1)
            self.dEntries.append(Entry(self.frame3,textvariable=self.dVal[i],justify="center",width=10))
            self.dEntries[i].grid(row=i+2,column=0)
            self.posCheck.append(Checkbutton(self.frame3,variable=self.posTick[i]))
            self.posCheck[i].grid(row=i+2,column=1)
            self.posEntries.append(Entry(self.frame3,textvariable=self.posVal[i],justify="center",width=10))
            self.posEntries[i].grid(row=i+2,column=2)
            self.sigmaCheck.append(Checkbutton(self.frame3,variable=self.sigmaTick[i]))
            self.sigmaCheck[i].grid(row=i+2,column=3)
            self.sigmaEntries.append(Entry(self.frame3,textvariable=self.sigmaVal[i],justify="center",width=10))
            self.sigmaEntries[i].grid(row=i+2,column=4)
            self.ampCheck.append(Checkbutton(self.frame3,variable=self.ampTick[i]))
            self.ampCheck[i].grid(row=i+2,column=5)
            self.ampEntries.append(Entry(self.frame3,textvariable=self.ampVal[i],justify="center",width=10))
            self.ampEntries[i].grid(row=i+2,column=6)
            self.widthCheck.append(Checkbutton(self.frame3,variable=self.widthTick[i]))
            self.widthCheck[i].grid(row=i+2,column=7)
            self.widthEntries.append(Entry(self.frame3,textvariable=self.widthVal[i],justify="center",width=10))
            self.widthEntries[i].grid(row=i+2,column=8)
            self.gaussCheck.append(Checkbutton(self.frame3,variable=self.gaussTick[i]))
            self.gaussCheck[i].grid(row=i+2,column=9)
            self.gaussEntries.append(Entry(self.frame3,textvariable=self.gaussVal[i],justify="center",width=10))
            self.gaussEntries[i].grid(row=i+2,column=10)
            if i > 0:
                self.posCheck[i].grid_remove()
                self.posEntries[i].grid_remove()
                self.sigmaCheck[i].grid_remove()
                self.sigmaEntries[i].grid_remove()
                self.dEntries[i].grid_remove()
                self.ampCheck[i].grid_remove()
                self.ampEntries[i].grid_remove()
                self.widthCheck[i].grid_remove()
                self.widthEntries[i].grid_remove()
                self.gaussCheck[i].grid_remove()
                self.gaussEntries[i].grid_remove()

    def checkI(self,I):
        if I == '1':
            return 1.0
        elif I == '3/2':
            return 1.5
        elif I == '2':
            return 2
        elif I == '5/2':
            return 2.5
        elif I == '3':
            return 3
        elif I == '7/2':
            return 3.5
        elif I == '4':
            return 4
        elif I == '9/2':
            return 4.5
                
    def setCheng(self,*args):
        self.cheng = int(safeEval(self.chengVal.get()))
        self.chengVal.set(str(self.cheng))

    def setGrid(self, *args):
        self.cqGridVal(str(int(safeEval(self.cqGridVal.get()))))
        self.etaGridVal(str(int(safeEval(self.etaGridVal.get()))))
        
    def changeNum(self,*args):
        val = int(self.numExp.get())
        for i in range(10):
            if i < val:
                self.posCheck[i].grid()
                self.posEntries[i].grid()
                self.sigmaCheck[i].grid()
                self.sigmaEntries[i].grid()
                self.dEntries[i].grid()
                self.ampCheck[i].grid()
                self.ampEntries[i].grid()
                self.widthCheck[i].grid()
                self.widthEntries[i].grid()
                self.gaussCheck[i].grid()
                self.gaussEntries[i].grid()
            else:
                self.posCheck[i].grid_remove()
                self.posEntries[i].grid_remove()
                self.sigmaCheck[i].grid_remove()
                self.sigmaEntries[i].grid_remove()
                self.dEntries[i].grid_remove()
                self.ampCheck[i].grid_remove()
                self.ampEntries[i].grid_remove()
                self.widthCheck[i].grid_remove()
                self.widthEntries[i].grid_remove()
                self.gaussCheck[i].grid_remove()
                self.gaussEntries[i].grid_remove()

    def bincounting(self, x1, weight, length):
        weights = weight[np.logical_and(x1>=0,x1<length)]
        x1 = x1[np.logical_and(x1>=0,x1<length)]
        return np.fft.ifft(np.bincount(x1,weights,length))
    
    def genLib(self, length, I, maxCq, numCq, numEta):
        self.cq, self.eta = np.meshgrid(np.linspace(0,maxCq,numCq),np.linspace(0,1,numEta))
        cq=self.cq[...,None]
        eta = self.eta[...,None]
        v = -cq**2/(6.0*self.parent.current.freq)*(I*(I+1)-3/4.0)*(self.angleStuff1+self.angleStuff2*eta+self.angleStuff3*eta**2)
        mult=v/(self.parent.current.sw)*length
        x1=np.array(np.round(mult)+np.floor(length/2),dtype=int)
        self.lib = np.apply_along_axis(self.bincounting,2,x1,self.weight,length)

    def tensorFunc(self,sigma,d,pos,width,gauss):
        cq = self.cq
        eta = self.eta
        czjzek = cq**(d-1)*eta/(np.sqrt(2*np.pi)*sigma)*(1-eta**2/9.0)*np.exp(-cq**2/(2.0*sigma**2)*(1+eta**2/3.0))
        czjzek = czjzek/np.sum(czjzek)
        fid = np.sum(self.lib*czjzek[...,None],axis=(0,1))
        t=np.arange(len(fid))/self.parent.current.sw
        apod = np.exp(-np.pi*width*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
        apod[-1:-(len(apod)/2+1):-1]=apod[:len(apod)/2]
        return scipy.ndimage.interpolation.shift(np.real(np.fft.fft(fid*apod)),len(fid)*pos/self.parent.current.sw)
                
    def fitFunc(self, param, x, y):
        numExp = self.args[0]
        struc = self.args[1]
        argu = self.args[2]
        testFunc = np.zeros(len(self.parent.data1D))
        if struc[0]:
            bgrnd = param[0]
            param=np.delete(param,[0])
        else:
            bgrnd = argu[0]
            argu=np.delete(argu,[0])
        if struc[1]:
            slope = param[0]
            param=np.delete(param,[0])
        else:
            slope = argu[0]
            argu=np.delete(argu,[0])
        for i in range(numExp):
            d = argu[0]
            argu = np.delete(argu,[0])
            if struc[5*i+2]:
                pos = param[0]
                param = np.delete(param,[0])
            else:
                pos= argu[0]
                argu = np.delete(argu,[0])
            if struc[5*i+3]:
                sigma = param[0]
                param = np.delete(param,[0])
            else:
                sigma = argu[0]
                argu = np.delete(argu,[0])
            if struc[5*i+4]:
                amp = param[0]
                param = np.delete(param,[0])
            else:
                amp= argu[0]
                argu = np.delete(argu,[0])
            if struc[5*i+5]:
                width = abs(param[0])
                param = np.delete(param,[0])
            else:
                width = argu[0]
                argu = np.delete(argu,[0])
            if struc[5*i+6]:
                gauss = abs(param[0])
                param = np.delete(param,[0])
            else:
                gauss = argu[0]
                argu = np.delete(argu,[0])
            testFunc += amp*self.tensorFunc(sigma,d,pos,width,gauss)
        testFunc += bgrnd+slope*x
        return np.sum((np.real(testFunc)-y)**2)

    def disp(self,outBgrnd,outSlope,outPos,outSigma,outD,outAmp,outWidth,outGauss):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x=[]
        for i in range(len(outPos)):
            x.append(tmpx)
            y =  outAmp[i]*self.tensorFunc(outSigma[i],outD[i],outPos[i],outWidth[i],outGauss[i])
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        self.parent.showPlot(tmpx, outCurve, x, outCurvePart)

    def setAngleStuff(self):
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.angleStuff1 = -27/8.0*np.cos(theta)**4+15/4.0*np.cos(theta)**2-3/8.0
        self.angleStuff2 = (-9/4.0*np.cos(theta)**4+2*np.cos(theta)**2+1/4.0)*np.cos(2*phi)
        self.angleStuff3 = -1/2.0*np.cos(theta)**2+1/3.0+(-3/8.0*np.cos(theta)**4+3/4.0*np.cos(theta)**2-3/8.0)*np.cos(2*phi)**2
        
    def fit(self,*args):
        self.setCheng()
        struc = []
        guess = []
        argu = []
        maxCq = 0.0
        I = self.checkI(self.IVal.get())
        numExp = int(self.numExp.get())
        outPos = np.zeros(numExp)
        outSigma = np.zeros(numExp)
        outD = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if self.bgrndTick.get() == 0:
            guess.append(safeEval(self.bgrndVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.bgrndVal.get())
            argu.append(inp)
            outBgrnd = inp
            self.bgrndVal.set('%.2g' % inp)
            struc.append(False)
        if self.slopeTick.get() == 0:
            guess.append(safeEval(self.slopeVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.slopeVal.get())
            argu.append(inp)
            outSlope = inp
            self.slopeVal.set('%.2g' % inp)
            struc.append(False)
        for i in range(numExp):
            inp = int(safeEval(self.dVal[i].get()))
            if inp < 1:
                inp = 1
            elif inp > 5:
                inp = 5
            argu.append(inp)
            outD[i] = inp
            self.dVal[i].set('%.2g' % inp)
            if self.posTick[i].get() == 0:
                guess.append(safeEval(self.posVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.posVal[i].get())
                argu.append(inp)
                outPos[i] = inp
                self.posVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.sigmaTick[i].get() == 0:
                inp = safeEval(self.sigmaVal[i].get())
                maxCq = max(maxCq,inp*1e6)
                guess.append(inp*1e6)
                struc.append(True)
            else:
                inp = safeEval(self.sigmaVal[i].get())
                maxCq = max(maxCq,inp*1e6)
                argu.append(inp*1e6)
                outSigma[i] = inp
                self.sigmaVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.ampTick[i].get() == 0:
                guess.append(safeEval(self.ampVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.ampVal[i].get())
                argu.append(inp)
                outAmp[i] = inp
                self.ampVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.widthTick[i].get() == 0:
                guess.append(abs(safeEval(self.widthVal[i].get())))
                struc.append(True)
            else:
                inp = abs(safeEval(self.widthVal[i].get()))
                argu.append(inp)
                outWidth[i] = inp
                self.widthVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.gaussTick[i].get() == 0:
                guess.append(abs(safeEval(self.gaussVal[i].get())))
                struc.append(True)
            else:
                inp = abs(safeEval(self.gaussVal[i].get()))
                argu.append(inp)
                outGauss[i] = inp
                self.gaussVal[i].set('%.2g' % inp)
                struc.append(False)
        self.args = (numExp,struc,argu)
        self.setAngleStuff()
        numCq = int(safeEval(self.cqGridVal.get()))
        numEta = int(safeEval(self.etaGridVal.get()))
        self.genLib(len(self.parent.xax), I, maxCq*4.0, numCq, numEta)
        fitVal = scipy.optimize.fmin(self.fitFunc,guess, args=(self.parent.xax,np.real(self.parent.data1D)))
        counter = 0
        if struc[0]:
            self.bgrndVal.set('%.2g' % fitVal[counter])
            outBgrnd = fitVal[counter]
            counter +=1
        if struc[1]:
            self.slopeVal.set('%.2g' % fitVal[counter])
            outSlope = fitVal[counter]
            counter +=1
        for i in range(numExp):
            if struc[5*i+2]:
                self.posVal[i].set('%.2g' % fitVal[counter])
                outPos[i] = fitVal[counter]
                counter += 1
            if struc[5*i+3]:
                self.sigmaVal[i].set('%.2g' % (fitVal[counter]*1e-6))
                outCq[i] = fitVal[counter]
                counter += 1
            if struc[5*i+4]:
                self.ampVal[i].set('%.2g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[5*i+5]:
                self.widthVal[i].set('%.2g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
            if struc[5*i+6]:
                self.gaussVal[i].set('%.2g' % abs(fitVal[counter]))
                outGauss[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd,outSlope,outPos,outSigma,outD,outAmp,outWidth,outGauss)

    def sim(self):
        self.setCheng()
        numExp = int(self.numExp.get())
        bgrnd = safeEval(self.bgrndVal.get())
        slope = safeEval(self.slopeVal.get())
        pos = np.zeros(numExp)
        sigma = np.zeros(numExp)
        d = np.zeros(numExp)
        amp = np.zeros(numExp)
        width = np.zeros(numExp)
        gauss = np.zeros(numExp)
        I = self.checkI(self.IVal.get())
        for i in range(numExp):
            pos[i] = safeEval(self.posVal[i].get())
            sigma[i] = safeEval(self.sigmaVal[i].get())*1e6
            d[i] = safeEval(self.dVal[i].get())
            amp[i] = safeEval(self.ampVal[i].get())
            width[i] = safeEval(self.widthVal[i].get())
            gauss[i] = safeEval(self.gaussVal[i].get())
        self.setAngleStuff()
        numCq = int(safeEval(self.cqGridVal.get()))
        numEta = int(safeEval(self.etaGridVal.get()))
        self.genLib(len(self.parent.xax), I, max(sigma)*4.0, numCq, numEta)
        self.disp(bgrnd,slope,pos,sigma,d,amp,width,gauss)
    
#################################################################################
class Quad2MASCzjzekParamFrame(Quad2StaticCzjzekParamFrame): 

    Ioptions = ['3/2','5/2','7/2','9/2']
    
    def __init__(self, parent, rootwindow):
        Quad2StaticCzjzekParamFrame.__init__(self,parent,rootwindow)
        
    def setAngleStuff(self):
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.angleStuff1 = 21/16.0*np.cos(theta)**4-9/8.0*np.cos(theta)**2+5/16.0
        self.angleStuff2 = (-7/8.0*np.cos(theta)**4+np.cos(theta)**2-1/8.0)*np.cos(2*phi)
        self.angleStuff3 = 1/12.0*np.cos(theta)**2+(+7/48.0*np.cos(theta)**4-7/24.0*np.cos(theta)**2+7/48.0)*np.cos(2*phi)**2
         
#####################################################################################
class MainPlotWindow(Frame):
    def __init__(self,parent,mainProgram,oldMainWindow):
        Frame.__init__(self,parent)
        self.parent = parent 
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.fig = oldMainWindow.current.fig
        self.canvas = oldMainWindow.current.canvas
        self.ax = oldMainWindow.current.ax
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0,sticky="nw")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame1)
        self.canvas.get_tk_widget().pack(fill=BOTH,expand=1)
        self.frame2 = Frame(self)
        self.frame2.grid(row=0,column=1,sticky="ne")
        Label(self.frame2,text='Title').grid(row=0,column=0)
        self.title = StringVar()
        self.titleBackup = oldMainWindow.name
        self.title.set(self.titleBackup)
        self.titleEntry = Entry(self.frame2,textvariable=self.title,justify="center")
        self.titleEntry.bind("<Return>", self.updatePlot) 
        self.titleEntry.bind("<KP_Enter>", self.updatePlot) 
        self.titleEntry.grid(row=1,column=0)
        Label(self.frame2,text='x-label').grid(row=2,column=0)
        self.xlabel = StringVar()
        self.xlabelBackup = self.ax.get_xlabel()
        self.xlabel.set(self.xlabelBackup)
        self.xlabelEntry = Entry(self.frame2,textvariable=self.xlabel,justify="center")
        self.xlabelEntry.bind("<Return>", self.updatePlot) 
        self.xlabelEntry.bind("<KP_Enter>", self.updatePlot) 
        self.xlabelEntry.grid(row=3,column=0)
        Label(self.frame2,text='y-label').grid(row=4,column=0)
        self.ylabel = StringVar()
        self.ylabelBackup = self.ax.get_ylabel()
        self.ylabel.set(self.ylabelBackup)
        self.ylabelEntry = Entry(self.frame2,textvariable=self.ylabel,justify="center")
        self.ylabelEntry.bind("<Return>", self.updatePlot) 
        self.ylabelEntry.bind("<KP_Enter>", self.updatePlot) 
        self.ylabelEntry.grid(row=5,column=0)
        self.xlimBackup = self.ax.get_xlim()
        Label(self.frame2,text='x-limit left').grid(row=6,column=0)
        self.xlimLeft = StringVar()
        self.xlimLeft.set(self.xlimBackup[0])
        self.xlimLeftEntry = Entry(self.frame2,textvariable=self.xlimLeft,justify="center")
        self.xlimLeftEntry.bind("<Return>", self.updatePlot) 
        self.xlimLeftEntry.bind("<KP_Enter>", self.updatePlot) 
        self.xlimLeftEntry.grid(row=7,column=0)
        Label(self.frame2,text='x-limit right').grid(row=8,column=0)
        self.xlimRight = StringVar()
        self.xlimRight.set(self.xlimBackup[1])
        self.xlimRightEntry = Entry(self.frame2,textvariable=self.xlimRight,justify="center")
        self.xlimRightEntry.bind("<Return>", self.updatePlot) 
        self.xlimRightEntry.bind("<KP_Enter>", self.updatePlot) 
        self.xlimRightEntry.grid(row=9,column=0)
        self.ylimBackup = self.ax.get_ylim()
        Label(self.frame2,text='y-limit down').grid(row=10,column=0)
        self.ylimLeft = StringVar()
        self.ylimLeft.set(self.ylimBackup[0])
        self.ylimLeftEntry = Entry(self.frame2,textvariable=self.ylimLeft,justify="center")
        self.ylimLeftEntry.bind("<Return>", self.updatePlot) 
        self.ylimLeftEntry.bind("<KP_Enter>", self.updatePlot) 
        self.ylimLeftEntry.grid(row=11,column=0)
        Label(self.frame2,text='y-limit up').grid(row=12,column=0)
        self.ylimRight = StringVar()
        self.ylimRight.set(self.ylimBackup[1])
        self.ylimRightEntry = Entry(self.frame2,textvariable=self.ylimRight,justify="center")
        self.ylimRightEntry.bind("<Return>", self.updatePlot) 
        self.ylimRightEntry.bind("<KP_Enter>", self.updatePlot) 
        self.ylimRightEntry.grid(row=13,column=0)
        
        Label(self.frame2,text='Width [cm]').grid(row=26,column=0)
        self.widthBackup, self.heightBackup = self.fig.get_size_inches()
        self.widthBackup = self.widthBackup*2.54
        self.heightBackup = self.heightBackup*2.54
        self.width = StringVar()
        self.width.set(self.widthBackup)
        self.widthEntry = Entry(self.frame2,textvariable=self.width,justify="center")
        self.widthEntry.bind("<Return>", self.updatePlot) 
        self.widthEntry.bind("<KP_Enter>", self.updatePlot) 
        self.widthEntry.grid(row=27,column=0)
        Label(self.frame2,text='height [cm]').grid(row=28,column=0)
        self.height = StringVar()
        self.height.set(self.heightBackup)
        self.heightEntry = Entry(self.frame2,textvariable=self.height,justify="center")
        self.heightEntry.bind("<Return>", self.updatePlot) 
        self.heightEntry.bind("<KP_Enter>", self.updatePlot) 
        self.heightEntry.grid(row=29,column=0)
        Label(self.frame2,text='file type').grid(row=30,column=0)
        self.fileType = StringVar()
        self.fileType.set('svg')
        OptionMenu(self.frame2, self.fileType, self.fileType.get(), 'svg', 'png', 'eps', 'jpg', 'pdf').grid(row=31,column=0)
        self.inFrame = Frame(self.frame2)
        self.inFrame.grid(row=32,column=0)
        Button(self.inFrame,text='Save',command=self.save).grid(row=0,column=0)
        Button(self.inFrame,text='Cancel',command=self.cancel).grid(row=0,column=1)
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
    def rename(self,name):
        self.oldMainWindow.rename(name)
        
    def updatePlot(self, *args):
        self.fig.suptitle(self.titleEntry.get())
        self.ax.set_xlabel(self.xlabelEntry.get())
        self.ax.set_ylabel(self.ylabelEntry.get())
        self.ax.set_xlim((safeEval(self.xlimLeft.get()),safeEval(self.xlimRight.get())))
        self.ax.set_ylim((safeEval(self.ylimLeft.get()),safeEval(self.ylimRight.get())))
        self.fig.set_size_inches((int(safeEval(self.width.get()))/2.54,int(safeEval(self.height.get()))/2.54))
        self.fig.canvas.draw()

    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def kill(self):
        self.oldMainWindow.kill()
        self.destroy()
        
    def save(self):
        self.updatePlot()
        f=asksaveasfilename()
        if f:
            f=os.path.splitext(f)[0]+'.'+self.fileType.get()
            self.fig.savefig(f)
        self.cancel()

    def cancel(self):
        self.fig.suptitle(self.titleBackup)
        self.ax.set_xlabel(self.xlabelBackup)
        self.ax.set_ylabel(self.ylabelBackup)
        self.ax.set_xlim((self.xlimBackup[0],self.xlimBackup[1]))
        self.ax.set_ylim((self.ylimBackup[0],self.ylimBackup[1]))
        self.fig.set_size_inches((self.widthBackup/2.54,self.heightBackup/2.54))
        self.mainProgram.closeSaveFigure(self.oldMainWindow)

######################################################################

class FitAllSelectionWindow(Toplevel): #a window to select wich data fields should be saved after fitting
    def __init__(self, parent, fitNames):
        Toplevel.__init__(self)
        self.parent = parent
        self.geometry('+0+0')
        self.transient()
        self.title("Select output")
        self.resizable(width=FALSE, height=FALSE)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0)
        self.values = []
        for i in range(len(fitNames)):
            self.values.append(IntVar())
            self.values[-1].set(0)
            Checkbutton(self.frame1, text=fitNames[i], variable=self.values[i]).grid(row=i,column=0,sticky='w')
        self.frame2 = Frame(self)
        self.frame2.grid(row=1)
        Button(self.frame2, text="Fit",command=self.fit).grid(row=0,column=0)
        Button(self.frame2, text="Cancel",command=self.destroy).grid(row=0,column=1)
        
    def fit(self):
        returnVals = []
        for i in self.values:
            returnVals.append(i.get())
        self.parent.fitAllFunc(np.array(returnVals,dtype=bool))
        self.destroy()
