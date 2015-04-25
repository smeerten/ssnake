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

pi = math.pi

#################################################################################   
class RelaxWindow(Frame): #a window for shifting the data
    def __init__(self, parent,current):
        self.parent = parent
        self.current = current
        self.ax=np.linspace(0,2*np.pi*10,10)[:-1] #fake data
        self.data=1-1.2*np.exp(-1.0*self.ax/10.0)#fake data
        self.leftMouse = False        #is the left mouse button currently pressed
        self.panX = None              #start position of dragging the spectrum
        self.panY = None              #start position of dragging the spectrum 
        self.zoomX1 = None            #first corner of the zoombox 
        self.zoomY1 = None            #first corner of the zoombox
        self.zoomX2 = None            #second corner of the zoombox
        self.zoomY2 = None            #second corner of the zoombox
        self.rect=[None,None,None,None]      #lines for zooming or peak picking
        self.rightMouse = False              #is the right mouse button currently pressed
        self.root = Toplevel(self.parent.parent) #create the main window
        self.root.title("Relaxation Curve") 
        self.root.style = Style()
        self.root.style.theme_use("clam")
        self.root.columnconfigure(0, weight=1) #make the main window sizable
        self.root.rowconfigure(0, weight=1)
        self.root.attributes('-zoomed', True)
        Label(self.root, text="Amplitude * (Constant + Coefficient*exp(-time/T) + ...)").grid(row=1,column=0,sticky='n')
        self.paramframe = RelaxParamFrame(self)
        self.paramframe.grid(row=2,column=0,sticky='sw') 
        self.fig = Figure()
        self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(row=0,column=0,sticky="nswe")
        #connect click events to the canvas
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        self.plotReset()
        self.showPlot()

    def plotReset(self): #set the plot limits to min and max values
        a=self.fig.gca()
        miny=min(self.data)
        maxy=max(self.data)
        differ = 0.05*(maxy-miny) #amount to add to show all datapoints (10%)
        self.yminlim=miny-differ
        self.ymaxlim=maxy+differ
        self.xminlim=min(self.ax)
        self.xmaxlim=max(self.ax)
        a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)

    def showPlot(self, tmpAx=None, tmpdata=None): 
        a=self.fig.gca()
        a.cla()
        a.plot(tmpAx,tmpdata)
        a.scatter(self.ax,self.data)
        a.set_title("Relaxation Curve")
        a.set_xlabel('X axis label')
        a.set_ylabel('Y label')
        a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        a.get_xaxis().get_major_formatter().set_powerlimits((-2, 2))
        a.get_yaxis().get_major_formatter().set_powerlimits((-2, 2))
        self.canvas.draw()

    ################
    # mouse events #
    ################

    def scroll(self,event):
        a=self.fig.gca()
        if self.rightMouse:
            middle = (self.xmaxlim+self.xminlim)/2.0
            width = self.xmaxlim-self.xminlim
            width = width*0.9**event.step
            self.xmaxlim = middle+width/2.0
            self.xminlim = middle-width/2.0
            a.set_xlim(self.xminlim,self.xmaxlim)
        else:
            middle = (self.ymaxlim+self.yminlim)/2.0
            width = self.ymaxlim-self.yminlim
            width = width*0.9**event.step
            self.ymaxlim = middle+width/2.0
            self.yminlim = middle-width/2.0
            a.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()

    def buttonPress(self,event):
        if event.button == 1:
            self.leftMouse = True
            self.zoomX1 = event.xdata
            self.zoomY1 = event.ydata
        elif (event.button == 3) and event.dblclick:
            self.plotReset()
        elif event.button == 3:
            self.rightMouse = True
            self.panX = event.xdata
            self.panY = event.ydata

    def buttonRelease(self,event):
        a=self.fig.gca()
        if event.button == 1:
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
                a.set_xlim(self.xminlim,self.xmaxlim)
                a.set_ylim(self.yminlim,self.ymaxlim)
            self.zoomX1=None
            self.zoomX2=None #WF: should also be cleared, memory of old zoom
            self.zoomY1=None
            self.zoomY2=None #WF: should also be cleared, memory of old zoom
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw()

    def pan(self,event):
        if self.rightMouse and self.panX is not None and self.panY is not None:
            a=self.fig.gca()
            inv = a.transData.inverted()
            point = inv.transform((event.x,event.y))
            diffx = point[0]-self.panX
            diffy = point[1]-self.panY
            self.xmaxlim = self.xmaxlim-diffx
            self.xminlim = self.xminlim-diffx
            self.ymaxlim = self.ymaxlim-diffy
            self.yminlim = self.yminlim-diffy
            a.set_xlim(self.xminlim,self.xmaxlim)
            a.set_ylim(self.yminlim,self.ymaxlim)
            self.canvas.draw()
        elif self.leftMouse and (self.zoomX1 is not None) and (self.zoomY1 is not None):
            a=self.fig.gca()
            inv = a.transData.inverted()
            point = inv.transform((event.x,event.y))
            self.zoomX2 =  point[0]
            self.zoomY2 = point[1]
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
            self.rect[0],=a.plot([self.zoomX1,self.zoomX2],[self.zoomY2,self.zoomY2],'k',clip_on=False)
            self.rect[1],=a.plot([self.zoomX1,self.zoomX2],[self.zoomY1,self.zoomY1],'k',clip_on=False)
            self.rect[2],=a.plot([self.zoomX1,self.zoomX1],[self.zoomY1,self.zoomY2],'k',clip_on=False)
            self.rect[3],=a.plot([self.zoomX2,self.zoomX2],[self.zoomY1,self.zoomY2],'k',clip_on=False)
            self.canvas.draw()

#################################################################################
class RelaxParamFrame(Frame): #a window for shifting the data
    def __init__(self, parent): 
        self.parent = parent
        self.ampVal = StringVar()
        self.ampVal.set(str(np.amax(self.parent.data)))
        self.ampTick = IntVar()
        self.constVal = StringVar()
        self.constVal.set("1.0")
        self.constTick = IntVar()
        self.numExp = StringVar()
        Frame.__init__(self, self.parent.root)
        self.frame1 = Frame(self)
        self.frame1.grid(row=0,column=0,sticky='n')
        self.frame2 = Frame(self)
        self.frame2.grid(row=0,column=1,sticky='n')
        self.frame3 = Frame(self)
        self.frame3.grid(row=0,column=2,sticky='n')
        Button(self.frame1, text="Fit",command=self.fit).grid(row=0,column=0)
        Button(self.frame1, text="Cancel",command=self.parent.root.destroy).grid(row=1,column=0)
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

    def fitFunc(self, param, numExp, struc, argu):
        testFunc = np.zeros(len(self.parent.data))
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
            testFunc += coeff*np.exp(-1.0*self.parent.ax/T1) 
        return sum((self.parent.data - amplitude*(constant+testFunc))**2)


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
            self.ampVal.set('%.2f' % inp)
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
                self.T1Val[i].set('%.2f' % inp)
                struc.append(False)
        fitVal = scipy.optimize.minimize(self.fitFunc,guess,(numExp,struc,argu),'Nelder-Mead')
        counter = 0
        if struc[0]:
            self.ampVal.set('%.2f' % fitVal['x'][counter])
            outAmp = fitVal['x'][counter]
            counter +=1
        if struc[1]:
            self.constVal.set('%.2f' % fitVal['x'][counter])
            outConst = fitVal['x'][counter]
            counter +=1
        for i in range(1,numExp+1):
            if struc[2*i]:
                self.coeffVal[i-1].set('%.2f' % fitVal['x'][counter])
                outCoeff[i-1] = fitVal['x'][counter]
                counter += 1
            if struc[2*i+1]:
                self.T1Val[i-1].set('%.2f' % fitVal['x'][counter])
                outT1[i-1] = fitVal['x'][counter]
                counter += 1
        outCurve = np.zeros(numCurve)
        x = np.linspace(min(self.parent.ax),max(self.parent.ax),numCurve)
        for i in range(numExp):
            outCurve += outCoeff[i]*np.exp(-x/outT1[i])
        self.parent.showPlot(x, outAmp*(outConst+outCurve))
