import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec

import sys
if sys.version_info >= (3,0):
    from tkinter import *
else:
    from Tkinter import *
import spectrum_classes


#########################################################################################################
#the class from which the 1d data is displayed, the operations which only edit the content of this class are for previewing
class Plot1DFrame(Frame):
    def __init__(self, root):
        Frame.__init__(self,root)
        self.fig = Figure()           #figure
        if isinstance(self,spectrum_classes.CurrentContour):
            gs = gridspec.GridSpec(2, 2,width_ratios=[3,1],height_ratios=[1,3])
            self.ax = self.fig.add_subplot(gs[2])
            self.x_ax = self.fig.add_subplot(gs[0],sharex=self.ax)
            self.y_ax = self.fig.add_subplot(gs[3],sharey=self.ax) 
        else:
            self.ax = self.fig.add_subplot(111) 
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(fill=BOTH,expand=1)
        self.leftMouse = False        #is the left mouse button currently pressed
        self.panX = None              #start position of dragging the spectrum
        self.panY = None              #start position of dragging the spectrum 
        self.zoomX1 = None            #first corner of the zoombox 
        self.zoomY1 = None            #first corner of the zoombox
        self.zoomX2 = None            #second corner of the zoombox
        self.zoomY2 = None            #second corner of the zoombox
        self.rect=[None,None,None,None]      #lines for zooming or peak picking
        self.rightMouse = False              #is the right mouse button currently pressed
        self.peakPick = False         #currently peakPicking
        self.peakPickFunc = None      #the function that needs to be called after peakPicking
        #connect click events to the canvas
        self.canvas.mpl_connect('button_press_event', self.buttonPress)      
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        #variables to be initialized
        self.spec = 0

    def plotReset(self): #this function needs to be overriden by the classes who inherit from Plot1DFrame
        pass

    ################
    # mouse events #
    ################

    def peakPickReset(self):
        self.rect=[None,None,None,None]
        self.peakPick=False
        self.peakPickFunc = None

    def scroll(self,event):
        if self.rightMouse:
            middle = (self.xmaxlim+self.xminlim)/2.0
            width = self.xmaxlim-self.xminlim
            width = width*0.9**event.step
            self.xmaxlim = middle+width/2.0
            self.xminlim = middle-width/2.0
            if self.spec > 0:
                self.ax.set_xlim(self.xmaxlim,self.xminlim)
            else:
                self.ax.set_xlim(self.xminlim,self.xmaxlim)
        else:
            middle = (self.ymaxlim+self.yminlim)/2.0
            width = self.ymaxlim-self.yminlim
            width = width*0.9**event.step
            self.ymaxlim = middle+width/2.0
            self.yminlim = middle-width/2.0
            self.ax.set_ylim(self.yminlim,self.ymaxlim)
        self.canvas.draw()

    def buttonPress(self,event):
        if event.button == 1 and not self.peakPick:
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
        if event.button == 1:
            if self.peakPick:
                if self.rect[0] is not None:
                    self.rect[0].remove()
                    self.rect[0]=None
                    self.peakPick = False
                    minim = np.min(np.abs(self.line[0].get_xdata()-event.xdata))
                    minPos = 0
                    for i in range(1,len(self.line)):
                        minimNew = np.min(np.abs(self.line[i].get_xdata()-event.xdata))
                        if minimNew < minim:
                            minim = minimNew
                            minPos = i
                    xdata = self.line[minPos].get_xdata()
                    ydata = self.line[minPos].get_ydata()
                    idx = np.argmin(np.abs(xdata-event.xdata))
                    if self.peakPickFunc is not None:
                        self.peakPickFunc((idx,xdata[idx],ydata[idx]))
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
                    if self.spec > 0:
                        self.ax.set_xlim(self.xmaxlim,self.xminlim)
                    else:
                        self.ax.set_xlim(self.xminlim,self.xmaxlim)
                    self.ax.set_ylim(self.yminlim,self.ymaxlim)
                self.zoomX1=None
                self.zoomX2=None #WF: should also be cleared, memory of old zoom
                self.zoomY1=None
                self.zoomY2=None #WF: should also be cleared, memory of old zoom
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw()

    def pan(self,event):
        if self.rightMouse and self.panX is not None and self.panY is not None:
            inv = self.ax.transData.inverted()
            point = inv.transform((event.x,event.y))
            diffx = point[0]-self.panX
            diffy = point[1]-self.panY
            self.xmaxlim = self.xmaxlim-diffx
            self.xminlim = self.xminlim-diffx
            self.ymaxlim = self.ymaxlim-diffy
            self.yminlim = self.yminlim-diffy
            if self.spec > 0:
                self.ax.set_xlim(self.xmaxlim,self.xminlim)
            else:
                self.ax.set_xlim(self.xminlim,self.xmaxlim)
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
            inv = self.ax.transData.inverted()
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
            self.rect[0],=self.ax.plot([self.zoomX1,self.zoomX2],[self.zoomY2,self.zoomY2],'k',clip_on=False)
            self.rect[1],=self.ax.plot([self.zoomX1,self.zoomX2],[self.zoomY1,self.zoomY1],'k',clip_on=False)
            self.rect[2],=self.ax.plot([self.zoomX1,self.zoomX1],[self.zoomY1,self.zoomY2],'k',clip_on=False)
            self.rect[3],=self.ax.plot([self.zoomX2,self.zoomX2],[self.zoomY1,self.zoomY2],'k',clip_on=False)
            self.canvas.draw()
