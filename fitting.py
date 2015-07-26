#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
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
import scipy.optimize
import math
from safeEval import safeEval
from spectrumFrame import Plot1DFrame
from zcw import *

pi = math.pi

##############################################################################
class RelaxWindow(Frame): #a window for fitting relaxation data
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.current = RelaxFrame(self,oldMainWindow.current)
        self.current.grid(row=0,column=0,sticky='nswe')
        self.current.rowconfigure(0, weight=1)
        self.current.columnconfigure(0, weight=1)
        self.current.fig.set_size_inches(oldMainWindow.current.fig.get_size_inches())
        self.paramframe = RelaxParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def cancel(self):
        self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class RelaxFrame(Plot1DFrame): #a window for fitting relaxation data
    def __init__(self, rootwindow,current):
        axAdd=0
        if current.spec == 1:
            if current.ppm:
                axAdd = (current.freq-current.ref)/current.ref*1e6
                axMult = 1e6/current.ref
            else:
                axMult = 1.0/(1000.0**current.axType)
        elif current.spec == 0:
            axMult = 1000.0**current.axType
        self.xax = current.xax*axMult+axAdd
        self.data1D=current.getDisplayedData()
        self.plotType = 0
        self.logx = 0
        self.logy = 0
        Plot1DFrame.__init__(self,rootwindow)
        self.current = current
        self.rootwindow = rootwindow
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
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)

    def showPlot(self, tmpAx=None, tmpdata=None): 
        self.ax.cla()
        self.ax.plot(tmpAx,tmpdata)
        self.ax.scatter(self.xax,self.data1D)
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
class RelaxParamFrame(Frame): #a frame for the relaxtion parameters
    def __init__(self, parent, rootwindow): 
        self.parent = parent
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
        Button(self.frame1, text="Cancel",command=rootwindow.cancel).grid(row=2,column=0)
        Label(self.frame2,text="Amplitude").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.ampTick).grid(row=1,column=0)
        Entry(self.frame2,textvariable=self.ampVal,justify="center").grid(row=1,column=1)
        Label(self.frame2,text="Constant").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.constTick).grid(row=3,column=0)
        Entry(self.frame2,textvariable=self.constVal,justify="center").grid(row=3,column=1)
        OptionMenu(self.frame3, self.numExp, "1","1", "2", "3","4",command=self.changeNum).grid(row=0,column=0,columnspan=4)
        Label(self.frame3,text="Coefficient").grid(row=1,column=0,columnspan=2)
        Label(self.frame3,text="T").grid(row=1,column=2,columnspan=2)
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
            self.coeffEntries.append(Entry(self.frame3,textvariable=self.coeffVal[i],justify="center"))
            self.coeffEntries[i].grid(row=i+2,column=1)
            self.T1Check.append(Checkbutton(self.frame3,variable=self.T1Tick[i]))
            self.T1Check[i].grid(row=i+2,column=2)
            self.T1Entries.append(Entry(self.frame3,textvariable=self.T1Val[i],justify="center"))
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
                T1 = param[0]
                param=np.delete(param,[0])
            else:
                T1= argu[0]
                argu=np.delete(argu,[0])
            testFunc += coeff*np.exp(-1.0*x/T1) 
        return amplitude*(constant+testFunc)

    def fit(self,*args):
        #structure of the fitting arguments is : [amp,cost, coeff1, t1, coeff2, t2, coeff3, t3, coeff4, t4]
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
        x = np.linspace(min(self.parent.xax),max(self.parent.xax),numCurve)
        for i in range(len(outCoeff)):
            outCurve += outCoeff[i]*np.exp(-x/outT1[i])
        self.parent.showPlot(x, outAmp*(outConst+outCurve))
        
##############################################################################
class PeakDeconvWindow(Frame): #a window for fitting relaxation data
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.current = PeakDeconvFrame(self,oldMainWindow.current)
        self.current.grid(row=0,column=0,sticky='nswe')
        self.current.rowconfigure(0, weight=1)
        self.current.columnconfigure(0, weight=1)
        self.paramframe = PeakDeconvParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class PeakDeconvFrame(Plot1DFrame): #a window for fitting relaxation data
    def __init__(self, rootwindow,current):
        Plot1DFrame.__init__(self,rootwindow)
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
        self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
        self.peakPick = True
        self.pickNum = 0
        self.plotReset()
        self.showPlot()

    def plotReset(self): #set the plot limits to min and max values
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
        differ = 0.05*(maxy-miny) #amount to add to show all datapoints (10%)
        self.yminlim=miny-differ
        self.ymaxlim=maxy+differ
        self.xminlim=min(self.xax)
        self.xmaxlim=max(self.xax)
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim,self.xminlim)
        else:
            a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        a=self.fig.gca()
        a.cla()
        self.line = a.plot(self.xax,self.data1D)
        a.plot(tmpAx,tmpdata)
        for i in range(len(tmpAx2)):
            a.plot(tmpAx2[i],tmpdata2[i])
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

    def pickDeconv(self, pos):
        self.rootwindow.paramframe.posVal[self.pickNum].set("%.2g" %pos[1])
        left = pos[0] - 10 #number of points to find the maximum in
        if left < 0:
            left = 0
        right = pos[0] + 10 #number of points to find the maximum in
        if right >= len(self.data1D) :
            right = len(self.data1D)-1
        width = self.current.fwhm(left,right)
        self.rootwindow.paramframe.ampVal[self.pickNum].set("%.2g" %(pos[2]*width*0.5*np.pi))
        self.rootwindow.paramframe.widthVal[self.pickNum].set("%.2g" % width)
        if self.pickNum < 10:
            self.rootwindow.paramframe.numExp.set(str(self.pickNum+1))
            self.rootwindow.paramframe.changeNum()
        self.pickNum += 1
        if self.pickNum < 10:
            self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
            self.peakPick = True 

#################################################################################
class PeakDeconvParamFrame(Frame): #a frame for the relaxtion parameters
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
        Label(self.frame2,text="Bgrnd").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.bgrndTick).grid(row=1,column=0)
        Entry(self.frame2,textvariable=self.bgrndVal,justify="center").grid(row=1,column=1)
        Label(self.frame2,text="Slope").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.slopeTick).grid(row=3,column=0)
        Entry(self.frame2,textvariable=self.slopeVal,justify="center").grid(row=3,column=1)
        OptionMenu(self.frame3, self.numExp, "1","1", "2", "3","4","5","6","7","8","9","10",command=self.changeNum).grid(row=0,column=0,columnspan=6)
        Label(self.frame3,text="Position").grid(row=1,column=0,columnspan=2)
        Label(self.frame3,text="Amplitude").grid(row=1,column=2,columnspan=2)
        Label(self.frame3,text="Width").grid(row=1,column=4,columnspan=2)
        self.posVal = []
        self.posTick = []
        self.ampVal = []
        self.ampTick = []
        self.widthVal = []
        self.widthTick = []
        self.posCheck = []
        self.posEntries = []
        self.ampCheck = []
        self.ampEntries = []
        self.widthCheck = []
        self.widthEntries = []

        for i in range(10):
            self.posVal.append(StringVar())
            self.posVal[i].set("0.0")
            self.posTick.append(IntVar())
            self.ampVal.append(StringVar())
            self.ampVal[i].set("1.0")
            self.ampTick.append(IntVar())
            self.widthVal.append(StringVar())
            self.widthVal[i].set("1.0")
            self.widthTick.append(IntVar())
            self.posCheck.append(Checkbutton(self.frame3,variable=self.posTick[i]))
            self.posCheck[i].grid(row=i+2,column=0)
            self.posEntries.append(Entry(self.frame3,textvariable=self.posVal[i],justify="center"))
            self.posEntries[i].grid(row=i+2,column=1)
            self.ampCheck.append(Checkbutton(self.frame3,variable=self.ampTick[i]))
            self.ampCheck[i].grid(row=i+2,column=2)
            self.ampEntries.append(Entry(self.frame3,textvariable=self.ampVal[i],justify="center"))
            self.ampEntries[i].grid(row=i+2,column=3)
            self.widthCheck.append(Checkbutton(self.frame3,variable=self.widthTick[i]))
            self.widthCheck[i].grid(row=i+2,column=4)
            self.widthEntries.append(Entry(self.frame3,textvariable=self.widthVal[i],justify="center"))
            self.widthEntries[i].grid(row=i+2,column=5)
            if i > 0:
                self.posCheck[i].grid_remove()
                self.posEntries[i].grid_remove()
                self.ampCheck[i].grid_remove()
                self.ampEntries[i].grid_remove()
                self.widthCheck[i].grid_remove()
                self.widthEntries[i].grid_remove()

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
            else:
                self.posCheck[i].grid_remove()
                self.posEntries[i].grid_remove()
                self.ampCheck[i].grid_remove()
                self.ampEntries[i].grid_remove()
                self.widthCheck[i].grid_remove()
                self.widthEntries[i].grid_remove()

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
        for i in range(1,numExp+1):
            if struc[3*i-1]:
                pos = param[0]
                param=np.delete(param,[0])
            else:
                pos= argu[0]
                argu=np.delete(argu,[0])
            if struc[3*i]:
                amp = param[0]
                param=np.delete(param,[0])
            else:
                amp= argu[0]
                argu=np.delete(argu,[0])
            if struc[3*i+1]:
                width = abs(param[0])
                param=np.delete(param,[0])
            else:
                width = argu[0]
                argu=np.delete(argu,[0])
            testFunc += amp/np.pi*0.5*width/((x-pos)**2+(0.5*width)**2)
        testFunc += bgrnd+slope*x
        return testFunc

    def fit(self,*args):
        #structure of the fitting arguments is : [amp,cost, coeff1, t1, coeff2, t2, coeff3, t3, coeff4, t4]
        struc = []
        guess = []
        argu = []
        numExp = int(self.numExp.get())
        outPos = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
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
        for i in range(1,numExp+1):
            if struc[3*i-1]:
                self.posVal[i-1].set('%.3g' % fitVal[0][counter])
                outPos[i-1] = fitVal[0][counter]
                counter += 1
            if struc[3*i]:
                self.ampVal[i-1].set('%.3g' % fitVal[0][counter])
                outAmp[i-1] = fitVal[0][counter]
                counter += 1
            if struc[3*i+1]:
                self.widthVal[i-1].set('%.3g' % abs(fitVal[0][counter]))
                outWidth[i-1] = abs(fitVal[0][counter])
                counter += 1
        self.disp(outBgrnd, outSlope, outAmp, outPos, outWidth)

    def sim(self):
        numExp = int(self.numExp.get())
        outPos = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outBgrnd = safeEval(self.bgrndVal.get())
        outSlope = safeEval(self.slopeVal.get())
        for i in range(numExp):
            outPos[i] = safeEval(self.posVal[i].get())
            outAmp[i] = safeEval(self.ampVal[i].get())
            outWidth[i] = abs(safeEval(self.widthVal[i].get()))
        self.disp(outBgrnd, outSlope, outAmp, outPos, outWidth)

    def disp(self, outBgrnd, outSlope, outAmp, outPos, outWidth):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x=[]
        for i in range(len(outAmp)):
            x.append(tmpx)
            y =  outAmp[i]/np.pi*0.5*outWidth[i]/((tmpx-outPos[i])**2+(0.5*outWidth[i])**2)
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        self.parent.showPlot(tmpx, outCurve, x, outCurvePart)

##############################################################################
class TensorDeconvWindow(Frame): #a window for fitting relaxation data
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.current = TensorDeconvFrame(self,oldMainWindow.current)
        self.current.grid(row=0,column=0,sticky='nswe')
        self.current.rowconfigure(0, weight=1)
        self.current.columnconfigure(0, weight=1)
        self.paramframe = TensorDeconvParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class TensorDeconvFrame(Plot1DFrame): #a window for fitting relaxation data
    def __init__(self, rootwindow,current):
        Plot1DFrame.__init__(self,rootwindow)
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
        self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
        self.peakPick = True
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
        self.line = self.ax.plot(self.xax,self.data1D)
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
        self.chengEntry = Entry(self.optframe,textvariable=self.chengVal,justify="center")
        self.chengEntry.grid(row=1,column=0)
        self.chengEntry.bind("<Return>", self.setCheng) 
        self.chengEntry.bind("<KP_Enter>", self.setCheng)
        Label(self.frame2,text="Bgrnd").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.bgrndTick).grid(row=1,column=0)
        Entry(self.frame2,textvariable=self.bgrndVal,justify="center").grid(row=1,column=1)
        Label(self.frame2,text="Slope").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.slopeTick).grid(row=3,column=0)
        Entry(self.frame2,textvariable=self.slopeVal,justify="center").grid(row=3,column=1)
        OptionMenu(self.frame3, self.numExp, "1","1", "2", "3","4","5","6","7","8","9","10",command=self.changeNum).grid(row=0,column=0,columnspan=6)
        Label(self.frame3,text="T11").grid(row=1,column=0,columnspan=2)
        Label(self.frame3,text="T22").grid(row=1,column=2,columnspan=2)
        Label(self.frame3,text="T33").grid(row=1,column=4,columnspan=2)
        Label(self.frame3,text="Amplitude").grid(row=1,column=6,columnspan=2)
        Label(self.frame3,text="Width").grid(row=1,column=8,columnspan=2)
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

        for i in range(10):
            self.t11Val.append(StringVar())
            self.t11Val[i].set("0.0")
            self.t11Tick.append(IntVar())
            self.t22Val.append(StringVar())
            self.t22Val[i].set("0.0")
            self.t22Tick.append(IntVar())
            self.t33Val.append(StringVar())
            self.t33Val[i].set("0.0")
            self.t33Tick.append(IntVar())
            self.ampVal.append(StringVar())
            self.ampVal[i].set("1.0")
            self.ampTick.append(IntVar())
            self.widthVal.append(StringVar())
            self.widthVal[i].set("10.0")
            self.widthTick.append(IntVar())
            self.widthTick[i].set(1)
            self.t11Check.append(Checkbutton(self.frame3,variable=self.t11Tick[i]))
            self.t11Check[i].grid(row=i+2,column=0)
            self.t11Entries.append(Entry(self.frame3,textvariable=self.t11Val[i],justify="center"))
            self.t11Entries[i].grid(row=i+2,column=1)
            self.t22Check.append(Checkbutton(self.frame3,variable=self.t22Tick[i]))
            self.t22Check[i].grid(row=i+2,column=2)
            self.t22Entries.append(Entry(self.frame3,textvariable=self.t22Val[i],justify="center"))
            self.t22Entries[i].grid(row=i+2,column=3)
            self.t33Check.append(Checkbutton(self.frame3,variable=self.t33Tick[i]))
            self.t33Check[i].grid(row=i+2,column=4)
            self.t33Entries.append(Entry(self.frame3,textvariable=self.t33Val[i],justify="center"))
            self.t33Entries[i].grid(row=i+2,column=5)
            self.ampCheck.append(Checkbutton(self.frame3,variable=self.ampTick[i]))
            self.ampCheck[i].grid(row=i+2,column=6)
            self.ampEntries.append(Entry(self.frame3,textvariable=self.ampVal[i],justify="center"))
            self.ampEntries[i].grid(row=i+2,column=7)
            self.widthCheck.append(Checkbutton(self.frame3,variable=self.widthTick[i]))
            self.widthCheck[i].grid(row=i+2,column=8)
            self.widthEntries.append(Entry(self.frame3,textvariable=self.widthVal[i],justify="center"))
            self.widthEntries[i].grid(row=i+2,column=9)
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

    def tensorFunc(self, x, t11, t22, t33, width):
        t11=t11*self.multt11
        t22=t22*self.multt22
        t33=t33*self.multt33
        v=t11+t22+t33
        length =len(x)
        t=np.arange(length)/self.parent.current.sw
        final = np.zeros(length)
        mult=v/(self.parent.current.sw)*length
        x1=np.round(mult)
        weight = self.weight[np.logical_and(x1>-length,x1<length)]
        x1 = x1[np.logical_and(x1>-length,x1<length)]
        for i in range(len(x1)):
            final[x1[i]] += self.weight[i]
        apod = np.exp(-width*t)
        apod[-1:-(len(apod)/2+1):-1]=apod[:len(apod)/2]
        I=np.real(np.fft.fftshift(np.fft.fft(np.fft.ifft(final)*apod)))
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
            if struc[5*i+2]:
                t11 = param[0]
                param=np.delete(param,[0])
            else:
                t11= argu[0]
                argu=np.delete(argu,[0])
            if struc[5*i+3]:
                t22 = param[0]
                param=np.delete(param,[0])
            else:
                t22= argu[0]
                argu=np.delete(argu,[0])
            if struc[5*i+4]:
                t33 = param[0]
                param=np.delete(param,[0])
            else:
                t33= argu[0]
                argu=np.delete(argu,[0])
            if struc[5*i+5]:
                amp = param[0]
                param=np.delete(param,[0])
            else:
                amp= argu[0]
                argu=np.delete(argu,[0])
            if struc[5*i+6]:
                width = abs(param[0])
                param=np.delete(param,[0])
            else:
                width = argu[0]
                argu=np.delete(argu,[0])
            testFunc += amp*self.tensorFunc(x,t11,t22,t33,width)
        testFunc += bgrnd+slope*x
        return np.sum((np.real(testFunc)-y)**2)

    def disp(self,outBgrnd,outSlope,outt11,outt22,outt33,outAmp,outWidth):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x=[]
        for i in range(len(outt11)):
            x.append(tmpx)
            y =  outAmp[i]*self.tensorFunc(tmpx,outt11[i],outt22[i],outt33[i],outWidth[i])
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
            if struc[5*i+2]:
                self.t11Val[i].set('%.2g' % fitVal[counter])
                outt11[i] = fitVal[counter]
                counter += 1
            if struc[5*i+3]:
                self.t22Val[i].set('%.2g' % fitVal[counter])
                outt22[i] = fitVal[counter]
                counter += 1
            if struc[5*i+4]:
                self.t33Val[i].set('%.2g' % fitVal[counter])
                outt33[i] = fitVal[counter]
                counter += 1
            if struc[5*i+5]:
                self.ampVal[i].set('%.2g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[5*i+6]:
                self.widthVal[i].set('%.2g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd,outSlope,outt11,outt22,outt33,outAmp,outWidth)

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
        for i in range(numExp):
            t11[i] = safeEval(self.t11Val[i].get())
            t22[i] = safeEval(self.t22Val[i].get())
            t33[i] = safeEval(self.t33Val[i].get())
            amp[i] = safeEval(self.ampVal[i].get())
            width[i] = safeEval(self.widthVal[i].get())
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.multt11=np.sin(theta)**2*np.cos(phi)**2
        self.multt22=np.sin(theta)**2*np.sin(phi)**2
        self.multt33=np.cos(theta)**2
        self.disp(bgrnd,slope,t11,t22,t33,amp,width)
        
##############################################################################
class Quad1DeconvWindow(Frame): #a window for fitting relaxation data
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.current = Quad1DeconvFrame(self,oldMainWindow.current)
        self.current.grid(row=0,column=0,sticky='nswe')
        self.current.rowconfigure(0, weight=1)
        self.current.columnconfigure(0, weight=1)
        self.paramframe = Quad1DeconvParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class Quad1DeconvFrame(Plot1DFrame): #a window for fitting relaxation data
    def __init__(self, rootwindow,current):
        Plot1DFrame.__init__(self,rootwindow)
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
        self.line = self.ax.plot(self.xax,self.data1D)
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
        self.chengEntry = Entry(self.optframe,textvariable=self.chengVal,justify="center")
        self.chengEntry.grid(row=1,column=0)
        self.chengEntry.bind("<Return>", self.setCheng) 
        self.chengEntry.bind("<KP_Enter>", self.setCheng)
        Label(self.frame2,text="Bgrnd").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.bgrndTick).grid(row=1,column=0)
        Entry(self.frame2,textvariable=self.bgrndVal,justify="center").grid(row=1,column=1)
        Label(self.frame2,text="Slope").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.slopeTick).grid(row=3,column=0)
        Entry(self.frame2,textvariable=self.slopeVal,justify="center").grid(row=3,column=1)
        OptionMenu(self.frame3, self.numExp, "1","1", "2", "3","4","5","6","7","8","9","10",command=self.changeNum).grid(row=0,column=0,columnspan=6)
        Label(self.frame3,text="I").grid(row=1,column=0,columnspan=2)
        Label(self.frame3,text="Pos").grid(row=1,column=1,columnspan=2)
        Label(self.frame3,text="Cq").grid(row=1,column=3,columnspan=2)
        Label(self.frame3,text="Eta").grid(row=1,column=5,columnspan=2)
        Label(self.frame3,text="Amplitude").grid(row=1,column=7,columnspan=2)
        Label(self.frame3,text="Width").grid(row=1,column=9,columnspan=2)
        self.IVal = []
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
        self.IDrop = []
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

        for i in range(10):
            self.IVal.append(StringVar())
            self.IVal[i].set("3/2")
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
            self.IDrop.append(OptionMenu(self.frame3,self.IVal[i],self.IVal[i].get(),*self.Ioptions))
            self.IDrop[i].grid(row=i+2,column=0)
            self.posCheck.append(Checkbutton(self.frame3,variable=self.posTick[i]))
            self.posCheck[i].grid(row=i+2,column=1)
            self.posEntries.append(Entry(self.frame3,textvariable=self.posVal[i],justify="center"))
            self.posEntries[i].grid(row=i+2,column=2)
            self.cqCheck.append(Checkbutton(self.frame3,variable=self.cqTick[i]))
            self.cqCheck[i].grid(row=i+2,column=3)
            self.cqEntries.append(Entry(self.frame3,textvariable=self.cqVal[i],justify="center"))
            self.cqEntries[i].grid(row=i+2,column=4)
            self.etaCheck.append(Checkbutton(self.frame3,variable=self.etaTick[i]))
            self.etaCheck[i].grid(row=i+2,column=5)
            self.etaEntries.append(Entry(self.frame3,textvariable=self.etaVal[i],justify="center"))
            self.etaEntries[i].grid(row=i+2,column=6)
            self.ampCheck.append(Checkbutton(self.frame3,variable=self.ampTick[i]))
            self.ampCheck[i].grid(row=i+2,column=7)
            self.ampEntries.append(Entry(self.frame3,textvariable=self.ampVal[i],justify="center"))
            self.ampEntries[i].grid(row=i+2,column=8)
            self.widthCheck.append(Checkbutton(self.frame3,variable=self.widthTick[i]))
            self.widthCheck[i].grid(row=i+2,column=9)
            self.widthEntries.append(Entry(self.frame3,textvariable=self.widthVal[i],justify="center"))
            self.widthEntries[i].grid(row=i+2,column=10)
            if i > 0:
                self.IDrop[i].grid_remove()
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
                self.IDrop[i].grid()
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
            else:
                self.IDrop[i].grid_remove()
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

    def tensorFunc(self, x, I, pos, cq, eta, width):
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
        x1=np.round(mult)
        weights = weights[np.logical_and(x1>-length,x1<length)]
        x1 = x1[np.logical_and(x1>-length,x1<length)]
        for i in range(len(x1)):
            final[x1[i]] += weights[i]
        apod = np.exp(-width*t)
        apod[-1:-(len(apod)/2+1):-1]=apod[:len(apod)/2]
        Inten=np.real(np.fft.fftshift(np.fft.fft(np.fft.ifft(final)*apod)))
        return Inten
                
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
            if struc[5*i+2]:
                pos = param[0]
                param=np.delete(param,[0])
            else:
                pos= argu[0]
                argu=np.delete(argu,[0])
            if struc[5*i+3]:
                cq = param[0]
                param=np.delete(param,[0])
            else:
                cq= argu[0]
                argu=np.delete(argu,[0])
            if struc[5*i+4]:
                eta = param[0]
                param=np.delete(param,[0])
            else:
                eta= argu[0]
                argu=np.delete(argu,[0])
            if struc[5*i+5]:
                amp = param[0]
                param=np.delete(param,[0])
            else:
                amp= argu[0]
                argu=np.delete(argu,[0])
            if struc[5*i+6]:
                width = abs(param[0])
                param=np.delete(param,[0])
            else:
                width = argu[0]
                argu=np.delete(argu,[0])
            testFunc += amp*self.tensorFunc(x,I[i],pos,cq,eta,width)
        testFunc += bgrnd+slope*x
        return np.sum((np.real(testFunc)-y)**2)

    def disp(self,outBgrnd,outSlope,outI,outPos,outCq,outEta,outAmp,outWidth):
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x=[]
        for i in range(len(outPos)):
            x.append(tmpx)
            y =  outAmp[i]*self.tensorFunc(tmpx,outI[i],outPos[i],outCq[i],outEta[i],outWidth[i])
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
        I = []
        numExp = int(self.numExp.get())
        outPos = np.zeros(numExp)
        outCq = np.zeros(numExp)
        outEta = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
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
                outpos[i] = inp
                self.posVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.cqTick[i].get() == 0:
                guess.append(safeEval(self.cqVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.cqVal[i].get())
                argu.append(inp)
                outcq[i] = inp
                self.cqVal[i].set('%.2g' % inp)
                struc.append(False)
            if self.etaTick[i].get() == 0:
                guess.append(safeEval(self.etaVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.etaVal[i].get())
                argu.append(inp)
                outeta[i] = inp
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
            I.append(self.checkI(self.IVal[i].get()))
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
            if struc[5*i+2]:
                self.posVal[i].set('%.2g' % fitVal[counter])
                outpos[i] = fitVal[counter]
                counter += 1
            if struc[5*i+3]:
                self.cqVal[i].set('%.2g' % fitVal[counter])
                outcq[i] = fitVal[counter]
                counter += 1
            if struc[5*i+4]:
                self.etaVal[i].set('%.2g' % fitVal[counter])
                outeta[i] = fitVal[counter]
                counter += 1
            if struc[5*i+5]:
                self.ampVal[i].set('%.2g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[5*i+6]:
                self.widthVal[i].set('%.2g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd,outSlope,outpos,outcq,outeta,outAmp,outWidth)

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
        I = np.zeros(numExp)
        for i in range(numExp):
            pos[i] = safeEval(self.posVal[i].get())
            cq[i] = safeEval(self.cqVal[i].get())
            eta[i] = safeEval(self.etaVal[i].get())
            amp[i] = safeEval(self.ampVal[i].get())
            width[i] = safeEval(self.widthVal[i].get())
            I[i] = self.checkI(self.IVal[i].get())
        self.setAngleStuff()
        self.disp(bgrnd,slope,I,pos,cq,eta,amp,width)

##############################################################################
class Quad2DeconvWindow(Frame): #a window for fitting second order quadrupole lineshapes data
    def __init__(self, rootwindow,mainProgram,oldMainWindow,mas=False):
        Frame.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.current = Quad1DeconvFrame(self,oldMainWindow.current)
        self.current.grid(row=0,column=0,sticky='nswe')
        self.current.rowconfigure(0, weight=1)
        self.current.columnconfigure(0, weight=1)
        if mas:
            self.paramframe = Quad2MASDeconvParamFrame(self.current,self)
        else:
            self.paramframe = Quad2StaticDeconvParamFrame(self.current,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')
        self.rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)
         
#################################################################################
class Quad2StaticDeconvParamFrame(Quad1DeconvParamFrame): #a frame for the quadrupole parameters

    Ioptions = ['3/2','5/2','7/2','9/2']
    
    def __init__(self, parent, rootwindow):
        Quad1DeconvParamFrame.__init__(self,parent,rootwindow)
        
    def setAngleStuff(self):
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.angleStuff1 = -27/8.0*np.cos(theta)**4+15/4.0*np.cos(theta)**2-3/8.0
        self.angleStuff2 = (-9/4.0*np.cos(theta)**4+2*np.cos(theta)**2+1/4.0)*np.cos(2*phi)
        self.angleStuff3 = -1/2.0*np.cos(theta)**2+1/3.0+(-3/8.0*np.cos(theta)**4+3/4.0*np.cos(theta)**2-3/8.0)*np.cos(2*phi)**2

    def tensorFunc(self, x, I, pos, cq, eta, width):
        v = -cq**2/(6.0*self.parent.current.freq)*(I*(I+1)-3/4.0)*(self.angleStuff1+self.angleStuff2*eta+self.angleStuff3*eta**2)+pos
        length =len(x)
        t=np.arange(length)/self.parent.current.sw
        final = np.zeros(length)
        mult=v/(self.parent.current.sw)*length
        x1=np.round(mult)
        weights = self.weight[np.logical_and(x1>-length,x1<length)]
        x1 = x1[np.logical_and(x1>-length,x1<length)]
        for i in range(len(x1)):
            final[x1[i]] += weights[i]
        apod = np.exp(-width*t)
        apod[-1:-(len(apod)/2+1):-1]=apod[:len(apod)/2]
        Inten=np.real(np.fft.fftshift(np.fft.fft(np.fft.ifft(final)*apod)))
        return Inten
    
#################################################################################
class Quad2MASDeconvParamFrame(Quad2StaticDeconvParamFrame): #a frame for the quadrupole parameters
    Ioptions = ['3/2','5/2','7/2','9/2']
    
    def __init__(self, parent, rootwindow):
        Quad2StaticDeconvParamFrame.__init__(self,parent,rootwindow)
        
    def setAngleStuff(self):
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.angleStuff1 = 21/16.0*np.cos(theta)**4-9/8.0*np.cos(theta)**2+5/16.0
        self.angleStuff2 = (-7/8.0*np.cos(theta)**4+np.cos(theta)**2-1/8.0)*np.cos(2*phi)
        self.angleStuff3 = 1/12.0*np.cos(theta)**2+(+7/48.0*np.cos(theta)**4-7/24.0*np.cos(theta)**2+7/48.0)*np.cos(2*phi)**2
    
#####################################################################################
class MainPlotWindow(Frame):
    def __init__(self,parent,mainProgram,oldMainWindow):
        Frame.__init__(self,parent)
        self.parent = parent #remember your parents
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
        self.titleBackup = self.ax.get_title()
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
        
        Label(self.frame2,text='Width [inches]').grid(row=26,column=0)
        self.widthBackup, self.heightBackup = self.fig.get_size_inches()
        self.width = StringVar()
        self.width.set(self.widthBackup)
        self.widthEntry = Entry(self.frame2,textvariable=self.width,justify="center")
        self.widthEntry.bind("<Return>", self.updatePlot) 
        self.widthEntry.bind("<KP_Enter>", self.updatePlot) 
        self.widthEntry.grid(row=27,column=0)
        Label(self.frame2,text='height [inches]').grid(row=28,column=0)
        self.height = StringVar()
        self.height.set(self.heightBackup)
        self.heightEntry = Entry(self.frame2,textvariable=self.height,justify="center")
        self.heightEntry.bind("<Return>", self.updatePlot) 
        self.heightEntry.bind("<KP_Enter>", self.updatePlot) 
        self.heightEntry.grid(row=29,column=0)
        self.inFrame = Frame(self.frame2)
        self.inFrame.grid(row=30,column=0)
        Button(self.inFrame,text='Save',command=self.save).grid(row=0,column=0)
        Button(self.inFrame,text='Cancel',command=self.cancel).grid(row=0,column=1)

    def updatePlot(self, *args):
        self.ax.set_title(self.titleEntry.get())
        self.ax.set_xlabel(self.xlabelEntry.get())
        self.ax.set_ylabel(self.ylabelEntry.get())
        self.ax.set_xlim((safeEval(self.xlimLeft.get()),safeEval(self.xlimRight.get())))
        self.ax.set_ylim((safeEval(self.ylimLeft.get()),safeEval(self.ylimRight.get())))
        self.fig.set_size_inches((int(safeEval(self.width.get())),int(safeEval(self.height.get()))))
        self.fig.canvas.draw()
        
    def addToView(self):
        self.pack(fill=BOTH,expand=1)

    def removeFromView(self):
        self.pack_forget()

    def save(self):
        self.updatePlot()
        f=asksaveasfilename(filetypes=(('svg','.svg'),('png','.png'),('eps','.eps'),('jpg','.jpg'),('pdf','.pdf')))
        if f:
            self.fig.savefig(f)
        self.cancel()

    def cancel(self):
        self.ax.set_title(self.titleBackup)
        self.ax.set_xlabel(self.xlabelBackup)
        self.ax.set_ylabel(self.ylabelBackup)
        self.fig.set_size_inches((self.widthBackup,self.heightBackup))
        self.mainProgram.closeSaveFigure(self.oldMainWindow)
