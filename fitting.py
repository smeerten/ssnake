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
else:
    from Tkinter import *
    import Tkinter as tk
    from ttk import *
import scipy.optimize
import math
from safeEval import safeEval
from spectrumFrame import Plot1DFrame

pi = math.pi

##############################################################################
class RelaxWindow(Toplevel): #a window for fitting relaxation data
    def __init__(self, rootwindow,current):
        Toplevel.__init__(self,rootwindow)
        window = RelaxFrame(self,current)
        window.grid(row=0,column=0,sticky='nswe')
        window.rowconfigure(0, weight=1)
        window.columnconfigure(0, weight=1)
        self.paramframe = RelaxParamFrame(window,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')

#################################################################################   
class RelaxFrame(Plot1DFrame): #a window for fitting relaxation data
    def __init__(self, rootwindow,current):
        self.xax=current.xax
        self.data1D=current.getDisplayedData()
        self.plotType = 0
        Plot1DFrame.__init__(self,rootwindow)
        self.current = current
        self.rootwindow = rootwindow
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
        a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)

    def showPlot(self, tmpAx=None, tmpdata=None): 
        a=self.fig.gca()
        a.cla()
        a.plot(tmpAx,tmpdata)
        a.scatter(self.xax,self.data1D)
        a.set_xlabel('Time [s]')
        a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        a.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        a.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.canvas.draw()

#################################################################################
class RelaxParamFrame(Frame): #a frame for the relaxtion parameters
    def __init__(self, parent, rootwindow): 
        self.parent = parent
        self.ampVal = StringVar()
        self.ampVal.set("%.3e" % np.amax(self.parent.data1D))
        self.ampTick = IntVar()
        self.constVal = StringVar()
        self.constVal.set("1.0")
        self.constTick = IntVar()
        self.numExp = StringVar()
        Frame.__init__(self, rootwindow)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0,sticky='n')
        self.frame2 = Frame(self)
        self.frame2.grid(row=0,column=1,sticky='n')
        self.frame3 = Frame(self)
        self.frame3.grid(row=0,column=2,sticky='n')
        Button(self.frame1, text="Fit",command=self.fit).grid(row=0,column=0)
        Button(self.frame1, text="Cancel",command=rootwindow.destroy).grid(row=1,column=0)
        Label(self.frame2,text="Amplitude").grid(row=0,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.ampTick).grid(row=1,column=0)
        Entry(self.frame2,textvariable=self.ampVal,justify="center").grid(row=1,column=1)
        Label(self.frame2,text="Constant").grid(row=2,column=0,columnspan=2)
        Checkbutton(self.frame2,variable=self.constTick).grid(row=3,column=0)
        Entry(self.frame2,textvariable=self.constVal,justify="center").grid(row=3,column=1)
        OptionMenu(self.frame3, self.numExp, "1","1", "2", "3","4",command=self.changeNum).grid(row=0,column=0,columnspan=4)
        Label(self.frame3,text="Coefficient").grid(row=1,column=0,columnspan=2)
        Label(self.frame3,text="T").grid(row=1,column=2,columnspan=2)
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
        numCurve = 100 #number of points in output curve
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
            self.ampVal.set('%.3e' % inp)
            struc.append(False)
        if self.constTick.get() == 0:
            guess.append(safeEval(self.constVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.constVal.get())
            argu.append(inp)
            outConst = inp
            self.constVal.set('%.2f' % inp)
            struc.append(False)
        for i in range(numExp):
            if self.coeffTick[i].get() == 0:
                guess.append(safeEval(self.coeffVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.coeffVal[i].get())
                argu.append(inp)
                outCoeff[i] = inp
                self.coeffVal[i].set('%.2f' % inp)
                struc.append(False)
            if self.T1Tick[i].get() == 0:
                guess.append(safeEval(self.T1Val[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.T1Val[i].get())
                argu.append(inp)
                outT1[i] = inp
                self.T1Val[i].set('%.3e' % inp)
                struc.append(False)
        #fitVal = scipy.optimize.minimize(self.fitFunc,guess,(numExp,struc,argu),'Nelder-Mead',options={'xtol':1e-5,'ftol':1e-5})
        self.args=(numExp,struc,argu)
        fitVal = scipy.optimize.curve_fit(self.fitFunc,self.parent.xax, self.parent.data1D,guess)
        counter = 0
        if struc[0]:
            self.ampVal.set('%.3e' % fitVal[0][counter])
            outAmp = fitVal[0][counter]
            counter +=1
        if struc[1]:
            self.constVal.set('%.2f' % fitVal[0][counter])
            outConst = fitVal[0][counter]
            counter +=1
        for i in range(1,numExp+1):
            if struc[2*i]:
                self.coeffVal[i-1].set('%.2f' % fitVal[0][counter])
                outCoeff[i-1] = fitVal[0][counter]
                counter += 1
            if struc[2*i+1]:
                self.T1Val[i-1].set('%.3e' % fitVal[0][counter])
                outT1[i-1] = fitVal[0][counter]
                counter += 1
        outCurve = np.zeros(numCurve)
        x = np.linspace(min(self.parent.xax),max(self.parent.xax),numCurve)
        for i in range(numExp):
            outCurve += outCoeff[i]*np.exp(-x/outT1[i])
        self.parent.showPlot(x, outAmp*(outConst+outCurve))

##############################################################################
class PeakDeconvWindow(Toplevel): #a window for fitting relaxation data
    def __init__(self, rootwindow,current):
        Toplevel.__init__(self,rootwindow)
        window = PeakDeconvFrame(self,current)
        window.grid(row=0,column=0,sticky='nswe')
        window.rowconfigure(0, weight=1)
        window.columnconfigure(0, weight=1)
        self.paramframe = PeakDeconvParamFrame(window,self)
        self.paramframe.grid(row=1,column=0,sticky='sw')

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
        self.rootwindow.paramframe.posVal[self.pickNum].set("%.2e" %pos[1])
        self.rootwindow.paramframe.ampVal[self.pickNum].set("%.2e" %pos[2])
        left = pos[0] - 10 #number of points to find the maximum in
        if left < 0:
            left = 0
        right = pos[0] + 10 #number of points to find the maximum in
        if right >= len(self.data1D) :
            right = len(self.data1D)-1
        self.rootwindow.paramframe.widthVal[self.pickNum].set("%.2e" % self.current.fwhm(left,right))
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
        self.slopeVal = StringVar()
        self.slopeVal.set("0.0")
        self.slopeTick = IntVar()
        self.numExp = StringVar()
        Frame.__init__(self, rootwindow)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0,sticky='n')
        self.frame2 = Frame(self)
        self.frame2.grid(row=0,column=1,sticky='n')
        self.frame3 = Frame(self)
        self.frame3.grid(row=0,column=2,sticky='n')
        Button(self.frame1, text="Fit",command=self.fit).grid(row=0,column=0)
        Button(self.frame1, text="Cancel",command=rootwindow.destroy).grid(row=1,column=0)
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
            testFunc += amp/(1.0+((x-pos)/(0.5*width))**2)
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
            self.bgrndVal.set('%.2e' % inp)
            struc.append(False)
        if self.slopeTick.get() == 0:
            guess.append(safeEval(self.slopeVal.get()))
            struc.append(True)
        else:
            inp = safeEval(self.slopeVal.get())
            argu.append(inp)
            outSlope = inp
            self.slopeVal.set('%.2e' % inp)
            struc.append(False)
        for i in range(numExp):
            if self.posTick[i].get() == 0:
                guess.append(safeEval(self.posVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.posVal[i].get())
                argu.append(inp)
                outPos[i] = inp
                self.posVal[i].set('%.2e' % inp)
                struc.append(False)
            if self.ampTick[i].get() == 0:
                guess.append(safeEval(self.ampVal[i].get()))
                struc.append(True)
            else:
                inp = safeEval(self.ampVal[i].get())
                argu.append(inp)
                outAmp[i] = inp
                self.ampVal[i].set('%.2e' % inp)
                struc.append(False)
            if self.widthTick[i].get() == 0:
                guess.append(abs(safeEval(self.widthVal[i].get())))
                struc.append(True)
            else:
                inp = abs(safeEval(self.widthVal[i].get()))
                argu.append(inp)
                outWidth[i] = inp
                self.widthVal[i].set('%.2e' % inp)
                struc.append(False)
        self.args = (numExp,struc,argu)
        fitVal = scipy.optimize.curve_fit(self.fitFunc,self.parent.xax,self.parent.data1D,p0=guess)
        counter = 0
        if struc[0]:
            self.bgrndVal.set('%.2e' % fitVal[0][counter])
            outBgrnd = fitVal[0][counter]
            counter +=1
        if struc[1]:
            self.slopeVal.set('%.2e' % fitVal[0][counter])
            outSlope = fitVal[0][counter]
            counter +=1
        for i in range(1,numExp+1):
            if struc[3*i-1]:
                self.posVal[i-1].set('%.2e' % fitVal[0][counter])
                outPos[i-1] = fitVal[0][counter]
                counter += 1
            if struc[3*i]:
                self.ampVal[i-1].set('%.2e' % fitVal[0][counter])
                outAmp[i-1] = fitVal[0][counter]
                counter += 1
            if struc[3*i+1]:
                self.widthVal[i-1].set('%.2e' % abs(fitVal[0][counter]))
                outWidth[i-1] = abs(fitVal[0][counter])
                counter += 1 
        tmpx = self.parent.xax
        outCurveBase = outBgrnd + tmpx*outSlope
        outCurve = outCurveBase.copy()
        outCurvePart = []
        x=[]
        for i in range(numExp):
            x.append(tmpx)
            y =  outAmp[i]/(1.0+((tmpx-outPos[i])/(0.5*outWidth[i]))**2)
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        self.parent.showPlot(tmpx, outCurve, x, outCurvePart)
