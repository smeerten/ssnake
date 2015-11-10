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
from PyQt4 import QtGui, QtCore
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import scipy.optimize
import scipy.ndimage
import math
import os.path
import copy
from safeEval import *
from spectrumFrame import Plot1DFrame
from zcw import *

pi = math.pi

##############################################################################
class RelaxWindow(QtGui.QWidget): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        QtGui.QWidget.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        self.fig = Figure()
        self.canvas =  FigureCanvas(self.fig)
        grid = QtGui.QGridLayout(self)
        grid.addWidget(self.canvas,0,0)
        self.current = RelaxFrame(self,self.fig,self.canvas,oldMainWindow.current)
        self.paramframe = RelaxParamFrame(self.current,self)
        grid.addWidget(self.paramframe,1,0)
        grid.setColumnStretch(0,1)
        grid.setRowStretch(0,1)
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
        self.show()

    def removeFromView(self):
        self.hide()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
        self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class RelaxFrame(Plot1DFrame): 
    def __init__(self, rootwindow,fig,canvas,current):
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
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.xminlim=min(self.xax*axMult)
        self.xmaxlim=max(self.xax*axMult)
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)

    def showPlot(self, tmpAx=None, tmpdata=None): 
        self.ax.cla()
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if tmpAx is not None:
            self.ax.plot(tmpAx*axMult,tmpdata)
        self.ax.scatter(self.xax*axMult,self.data1D)
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
class RelaxParamFrame(QtGui.QWidget): 
    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1,0,0)
        grid.addLayout(self.optframe,0,1)
        grid.addLayout(self.frame2,0,2)
        grid.addLayout(self.frame3,0,3)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton,0,0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton,1,0)
        fitAllButton = QtGui.QPushButton("Fit all")
        fitAllButton.clicked.connect(self.fitAll)
        self.frame1.addWidget(fitAllButton,2,0)
        cancelButton = QtGui.QPushButton("Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton,3,0)
        self.frame1.setColumnStretch(10,1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel("Amplitude:"),0,0,1,2)
        self.ampTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.ampTick,1,0)
        self.ampEntry = QtGui.QLineEdit()
        self.ampEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ampEntry.setText("%.3g" % np.amax(self.parent.data1D))
        self.frame2.addWidget(self.ampEntry,1,1)
        self.frame2.addWidget(QLabel("Constant:"),2,0,1,2)
        self.constTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.constTick,3,0)
        self.constEntry = QtGui.QLineEdit()
        self.constEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.constEntry.setText("1.0")
        self.frame2.addWidget(self.constEntry,3,1)
        self.frame2.setColumnStretch(10,1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1','2','3','4'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp,0,0,1,2)
        self.frame3.addWidget(QLabel("Coefficient:"),1,0,1,2)
        self.frame3.addWidget(QLabel("T [s]:"),1,2,1,2)
        self.frame3.setColumnStretch(10,1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.xlog = QtGui.QCheckBox('x-log')
        self.xlog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.xlog,0,0,QtCore.Qt.AlignTop)
        self.ylog = QtGui.QCheckBox('y-log')
        self.ylog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.ylog,1,0,QtCore.Qt.AlignTop)
        self.optframe.setColumnStretch(10,1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        grid.setColumnStretch(10,1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        
        self.coeffTicks = []
        self.coeffEntries = []
        self.t1Ticks = []
        self.t1Entries = []
        for i in range(4):
            self.coeffTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.coeffTicks[i],i+2,0)
            self.coeffEntries.append(QtGui.QLineEdit())
            self.coeffEntries[i].setText("-1.0")
            self.coeffEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.coeffEntries[i],i+2,1)
            self.t1Ticks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.t1Ticks[i],i+2,2)
            self.t1Entries.append(QtGui.QLineEdit())
            self.t1Entries[i].setText("1.0")
            self.t1Entries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.t1Entries[i],i+2,3)
            if i > 0:
                self.coeffTicks[i].hide()
                self.coeffEntries[i].hide()
                self.t1Ticks[i].hide()
                self.t1Entries[i].hide()

    def setLog(self, *args):
        self.parent.setLog(self.xlog.isChecked(),self.ylog.isChecked())
                
    def changeNum(self,*args):
        val = self.numExp.currentIndex()+1
        for i in range(4):
            if i < val:
                self.coeffTicks[i].show()
                self.coeffEntries[i].show()
                self.t1Ticks[i].show()
                self.t1Entries[i].show()
            else:
                self.coeffTicks[i].hide()
                self.coeffEntries[i].hide()
                self.t1Ticks[i].hide()
                self.t1Entries[i].hide()

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

    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.ampEntry.text())
        if inp is None:
            return False
        self.ampEntry.setText('%.3g' % inp)
        inp = safeEval(self.constEntry.text())
        if inp is None:
            return False
        self.constEntry.setText('%.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.coeffEntries[i].text())
            if inp is None:
                return False
            self.coeffEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.t1Entries[i].text())
            if inp is None:
                return False
            self.t1Entries[i].setText('%.3g' % inp)
        return True
    
    def fit(self,*args):
        if not self.checkInputs():
            print("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outCoeff = np.zeros(numExp)
        outT1 = np.zeros(numExp)
        if not self.ampTick.isChecked():
            guess.append(float(self.ampEntry.text()))
            struc.append(True)
        else:
            outAmp = float(self.ampEntry.text())
            argu.append(outAmp)
            struc.append(False)
        if not self.constTick.isChecked():
            guess.append(float(self.constEntry.text()))
            struc.append(True)
        else:
            outConst = float(self.constEntry.text())
            argu.append(outConst)
            struc.append(False)
        for i in range(numExp):
            if not self.coeffTicks[i].isChecked():
                guess.append(float(self.coeffEntries[i].text()))
                struc.append(True)
            else:
                outCoeff[i] = float(self.coeffEntries[i].text())
                argu.append(outCoeff[i])
                struc.append(False)
            if not self.t1Ticks[i].isChecked():
                guess.append(float(self.t1Entries[i].text()))
                struc.append(True)
            else:
                outT1[i] = safeEval(self.t1Entries[i].text())
                argu.append(outT1[i])
                struc.append(False)
        self.args=(numExp,struc,argu)
        fitVal = scipy.optimize.curve_fit(self.fitFunc,self.parent.xax, self.parent.data1D,guess)
        counter = 0
        if struc[0]:
            self.ampEntry.setText('%.3g' % fitVal[0][counter])
            outAmp = fitVal[0][counter]
            counter +=1
        if struc[1]:
            self.constEntry.setText('%.3g' % fitVal[0][counter])
            outConst = fitVal[0][counter]
            counter +=1
        for i in range(1,numExp+1):
            if struc[2*i]:
                self.coeffEntries[i-1].setText('%.3g' % fitVal[0][counter])
                outCoeff[i-1] = fitVal[0][counter]
                counter += 1
            if struc[2*i+1]:
                self.t1Entries[i-1].setText('%.3g' % fitVal[0][counter])
                outT1[i-1] = fitVal[0][counter]
                counter += 1
        self.disp(outAmp,outConst,outCoeff,outT1)

    def fitAll(self, *args):
        FitAllSelectionWindow(self,["Amplitude","Constant","Coefficient","T"])
        
    def fitAllFunc(self,outputs):
        if not self.checkInputs():
            print("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outCoeff = np.zeros(numExp)
        outT1 = np.zeros(numExp)
        if not self.ampTick.isChecked():
            guess.append(float(self.ampEntry.text()))
            struc.append(True)
        else:
            outAmp = float(self.ampEntry.text())
            argu.append(outAmp)
            struc.append(False)
        if not self.constTick.isChecked():
            guess.append(float(self.constEntry.text()))
            struc.append(True)
        else:
            outConst = float(self.constEntry.text())
            argu.append(outConst)
            struc.append(False)
        for i in range(numExp):
            if not self.coeffTicks[i].isChecked():
                guess.append(float(self.coeffEntries[i].text()))
                struc.append(True)
            else:
                outCoeff[i] = float(self.coeffEntries[i].text())
                argu.append(outCoeff[i])
                struc.append(False)
            if not self.t1Ticks[i].isChecked():
                guess.append(float(self.t1Entries[i].text()))
                struc.append(True)
            else:
                outT1[i] = safeEval(self.t1Entries[i].text())
                argu.append(outT1[i])
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
        numExp = self.numExp.currentIndex()+1
        outAmp = safeEval(self.ampEntry.text())
        outConst = safeEval(self.constEntry.text())
        if outAmp is None or outConst is None:
            print("One of the inputs is not valid")
            return
        outCoeff = []
        outT1 = []
        for i in range(numExp):
            outCoeff.append(safeEval(self.coeffEntries[i].text()))
            outT1.append(safeEval(self.t1Entries[i].text()))
            if outCoeff[i] is None or outT1[i] is None:
                print("One of the inputs is not valid")
                return
        self.disp(outAmp,outConst,outCoeff,outT1)
        
    def disp(self, outAmp,outConst,outCoeff, outT1):
        numCurve = 100 #number of points in output curve
        outCurve = np.zeros(numCurve)
        if self.xlog.isChecked():
            x = np.logspace(np.log(min(self.parent.xax)),np.log(max(self.parent.xax)),numCurve)
        else:
            x = np.linspace(min(self.parent.xax),max(self.parent.xax),numCurve)
        for i in range(len(outCoeff)):
            outCurve += outCoeff[i]*np.exp(-x/outT1[i])
        self.parent.showPlot(x, outAmp*(outConst+outCurve))

##############################################################################
class DiffusionWindow(QtGui.QWidget): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        QtGui.QWidget.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        #self.fig = self.mainProgram.getFig()
        self.fig = Figure()
        self.canvas =  FigureCanvas(self.fig)
        grid = QtGui.QGridLayout(self)
        grid.addWidget(self.canvas,0,0)
        self.current = DiffusionFrame(self,self.fig,self.canvas,oldMainWindow.current)
        self.paramframe = DiffusionParamFrame(self.current,self)
        grid.addWidget(self.paramframe,1,0)
        grid.setColumnStretch(0,1)
        grid.setRowStretch(0,1)
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
        self.show()

    def removeFromView(self):
        self.hide()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
        self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class DiffusionFrame(Plot1DFrame): 
    def __init__(self, rootwindow,fig,canvas,current):
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
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        self.xminlim=min(self.xax*axMult)
        self.xmaxlim=max(self.xax*axMult)
        self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)

    def showPlot(self, tmpAx=None, tmpdata=None): 
        self.ax.cla()
        if self.spec == 1:
            if self.ppm:
                axMult = 1e6/self.ref
            else:
                axMult = 1.0/(1000.0**self.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.axType
        if tmpAx is not None:
            self.ax.plot(tmpAx*axMult,tmpdata)
        self.ax.scatter(self.xax*axMult,self.data1D)
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
class DiffusionParamFrame(QtGui.QWidget): 
    def __init__(self, parent, rootwindow): 
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        self.frame4 = QtGui.QGridLayout()
        grid.addLayout(self.frame1,0,0)
        grid.addLayout(self.optframe,0,1)
        grid.addLayout(self.frame2,0,2)
        grid.addLayout(self.frame3,0,3)
        grid.addLayout(self.frame4,0,4)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton,0,0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton,1,0)
        #fitAllButton = QtGui.QPushButton("Fit all")
        #fitAllButton.clicked.connect(self.fitAll)
        #self.frame1.addWidget(fitAllButton,2,0)
        cancelButton = QtGui.QPushButton("Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton,3,0)
        self.frame1.setColumnStretch(10,1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel(u"\u03b3 [MHz/T]:"),0,0)
        self.gammaEntry = QtGui.QLineEdit()
        self.gammaEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.gammaEntry.setText("42.576")
        self.frame2.addWidget(self.gammaEntry,1,0)
        self.frame2.addWidget(QLabel(u"\u03b4 [s]:"),2,0)
        self.deltaEntry = QtGui.QLineEdit()
        self.deltaEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.deltaEntry.setText("1.0")
        self.frame2.addWidget(self.deltaEntry,3,0)
        self.frame2.addWidget(QLabel(u"\u0394 [s]:"),4,0)
        self.triangleEntry = QtGui.QLineEdit()
        self.triangleEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.triangleEntry.setText("1.0")
        self.frame2.addWidget(self.triangleEntry,5,0)
        self.frame2.setColumnStretch(10,1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.frame3.addWidget(QLabel("Amplitude:"),0,0,1,2)
        self.ampTick = QtGui.QCheckBox('')
        self.frame3.addWidget(self.ampTick,1,0)
        self.ampEntry = QtGui.QLineEdit()
        self.ampEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ampEntry.setText("%.3g" % np.amax(self.parent.data1D))
        self.frame3.addWidget(self.ampEntry,1,1)
        self.frame3.addWidget(QLabel("Constant:"),2,0,1,2)
        self.constTick = QtGui.QCheckBox('')
        self.constTick.setChecked(True)
        self.frame3.addWidget(self.constTick,3,0)
        self.constEntry = QtGui.QLineEdit()
        self.constEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.constEntry.setText("0.0")
        self.frame3.addWidget(self.constEntry,3,1)
        self.frame3.setColumnStretch(10,1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1','2','3','4'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame4.addWidget(self.numExp,0,0,1,2)
        self.frame4.addWidget(QLabel("Coefficient:"),1,0,1,2)
        self.frame4.addWidget(QLabel("D [m^2/s]:"),1,2,1,2)
        self.frame4.setColumnStretch(20,1)
        self.frame4.setAlignment(QtCore.Qt.AlignTop)
        self.xlog = QtGui.QCheckBox('x-log')
        self.xlog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.xlog,0,0)
        self.ylog = QtGui.QCheckBox('y-log')
        self.ylog.stateChanged.connect(self.setLog)
        self.optframe.addWidget(self.ylog,1,0)
        self.optframe.setColumnStretch(10,1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        grid.setColumnStretch(10,1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        self.coeffEntries = []
        self.coeffTicks = []
        self.dEntries = []
        self.dTicks = []
        for i in range(4):
            self.coeffTicks.append(QtGui.QCheckBox(''))
            self.frame4.addWidget(self.coeffTicks[i],i+2,0)
            self.coeffEntries.append(QtGui.QLineEdit())
            self.coeffEntries[i].setText("1.0")
            self.coeffEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame4.addWidget(self.coeffEntries[i],i+2,1)
            self.dTicks.append(QtGui.QCheckBox(''))
            self.frame4.addWidget(self.dTicks[i],i+2,2)
            self.dEntries.append(QtGui.QLineEdit())
            self.dEntries[i].setText("1.0e-9")
            self.dEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame4.addWidget(self.dEntries[i],i+2,3)
            if i > 0:
                self.coeffTicks[i].hide()
                self.coeffEntries[i].hide()
                self.dTicks[i].hide()
                self.dEntries[i].hide()
                
    def setLog(self, *args):
        self.parent.setLog(self.xlog.isChecked(),self.ylog.isChecked())
                
    def changeNum(self,*args):
        val = self.numExp.currentIndex()+1
        for i in range(4):
            if i < val:
                self.coeffTicks[i].show()
                self.coeffEntries[i].show()
                self.dTicks[i].show()
                self.dEntries[i].show()
            else:
                self.coeffTicks[i].hide()
                self.coeffEntries[i].hide()
                self.dTicks[i].hide()
                self.dEntries[i].hide()

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
    
    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.gammaEntry.text())
        if inp is None:
            return False
        self.gammaEntry.setText('%.3g' % inp)
        inp = safeEval(self.deltaEntry.text())
        if inp is None:
            return False
        self.deltaEntry.setText('%.3g' % inp)
        inp = safeEval(self.triangleEntry.text())
        if inp is None:
            return False
        self.triangleEntry.setText('%.3g' % inp)
        inp = safeEval(self.ampEntry.text())
        if inp is None:
            return False
        self.ampEntry.setText('%.3g' % inp)
        inp = safeEval(self.constEntry.text())
        if inp is None:
            return False
        self.constEntry.setText('%.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.coeffEntries[i].text())
            if inp is None:
                return False
            self.coeffEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.dEntries[i].text())
            if inp is None:
                return False
            self.dEntries[i].setText('%.3g' % inp)
        return True
    
    def fit(self,*args):
        if not self.checkInputs():
            print("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outCoeff = np.zeros(numExp)
        outD = np.zeros(numExp)
        gamma = float(self.gammaEntry.text())
        delta = float(self.deltaEntry.text())
        triangle = float(self.triangleEntry.text())
        if not self.ampTick.isChecked():
            guess.append(float(self.ampEntry.text()))
            struc.append(True)
        else:
            outAmp = float(self.ampEntry.text())
            argu.append(outAmp)
            struc.append(False)
        if not self.constTick.isChecked():
            guess.append(float(self.constEntry.text()))
            struc.append(True)
        else:
            outConst = float(self.constEntry.text())
            argu.append(outConst)
            struc.append(False)
        for i in range(numExp):
            if not self.coeffTicks[i].isChecked():
                guess.append(float(self.coeffEntries[i].text()))
                struc.append(True)
            else:
                outCoeff[i] = float(self.coeffEntries[i].text())
                argu.append(outCoeff[i])
                struc.append(False)
            if not self.dTicks[i].isChecked():
                guess.append(float(self.dEntries[i].text()))
                struc.append(True)
            else:
                outD[i] = float(self.dEntries[i].text())
                argu.append(outD[i])
                struc.append(False)
        self.args=(numExp, struc, argu, gamma, delta, triangle)
        fitVal = scipy.optimize.curve_fit(self.fitFunc, self.parent.xax, self.parent.data1D, guess)
        counter = 0
        if struc[0]:
            self.ampEntry.setText('%.3g' % fitVal[0][counter])
            outAmp = fitVal[0][counter]
            counter +=1
        if struc[1]:
            self.constEntry.setText('%.3g' % fitVal[0][counter])
            outConst = fitVal[0][counter]
            counter +=1
        for i in range(1,numExp+1):
            if struc[2*i]:
                self.coeffEntries[i-1].setText('%.3g' % fitVal[0][counter])
                outCoeff[i-1] = fitVal[0][counter]
                counter += 1
            if struc[2*i+1]:
                self.dEntries[i-1].set('%.3g' % fitVal[0][counter])
                outD[i-1] = fitVal[0][counter]
                counter += 1
        self.disp(outAmp, outConst, outCoeff, outD, gamma, delta, triangle)

    def sim(self):
        numExp = self.numExp.currentIndex()+1
        outAmp = safeEval(self.ampEntry.text())
        outConst = safeEval(self.constEntry.text())
        gamma = safeEval(self.gammaEntry.text())
        delta = safeEval(self.deltaEntry.text())
        triangle = safeEval(self.triangleEntry.text())
        if not np.isfinite([outAmp, outConst, gamma, delta, triangle]).all():
            print("One of the inputs is not valid")
            return
        outCoeff = []
        outD = []
        for i in range(numExp):
            outCoeff.append(safeEval(self.coeffEntries[i].text()))
            outD.append(safeEval(self.dEntries[i].text()))
            if outCoeff[i] is None or outD[i] is None:
                print("One of the inputs is not valid")
                return
        self.disp(outAmp,outConst,outCoeff,outD,gamma,delta,triangle)
        
    def disp(self, outAmp, outConst, outCoeff, outD, gamma, delta, triangle):
        numCurve = 100 
        outCurve = np.zeros(numCurve)
        if self.xlog.isChecked():
            x = np.logspace(np.log(min(self.parent.xax)),np.log(max(self.parent.xax)),numCurve)
        else:
            x = np.linspace(min(self.parent.xax),max(self.parent.xax),numCurve)
        for i in range(len(outCoeff)):
            outCurve += outCoeff[i]*np.exp(-(gamma*delta*x)**2*outD[i]*(triangle-delta/3.0))
        self.parent.showPlot(x, outAmp*(outConst+outCurve))
        
##############################################################################
class PeakDeconvWindow(QtGui.QWidget): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        QtGui.QWidget.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        #self.fig = self.mainProgram.getFig()
        self.fig = Figure()
        self.canvas =  FigureCanvas(self.fig)
        grid = QtGui.QGridLayout(self)
        grid.addWidget(self.canvas,0,0)
        self.current = PeakDeconvFrame(self,self.fig,self.canvas,oldMainWindow.current)
        self.paramframe = PeakDeconvParamFrame(self.current,self)
        grid.addWidget(self.paramframe,1,0)
        grid.setColumnStretch(0,1)
        grid.setRowStretch(0,1)
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
        self.show()

    def removeFromView(self):
        self.hide()

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
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xminlim=min(self.xax*axMult)
        self.xmaxlim=max(self.xax*axMult)
        if self.spec > 0 :
            a.set_xlim(self.xmaxlim,self.xminlim)
        else:
            a.set_xlim(self.xminlim,self.xmaxlim)
        a.set_ylim(self.yminlim,self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        a=self.fig.gca()
        a.cla()
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.line_xdata = self.xax*axMult
        self.line_ydata = self.data1D
        a.plot(self.xax*axMult,self.data1D)
        if tmpAx is not None:
            a.plot(tmpAx*axMult,tmpdata)
        for i in range(len(tmpAx2)):
            a.plot(tmpAx2[i]*axMult,tmpdata2[i])
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
        self.peakPickReset()
        if var==1:
            self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
            self.peakPick = True
        else:
            self.peakPickFunc = None
            self.peakPick = False
        
    def pickDeconv(self, pos):
        if self.pickWidth:
            if self.current.spec == 1:
                if self.current.ppm:
                    axMult = 1e6/self.current.ref
                else:
                    axMult = 1.0/(1000.0**self.current.axType)
            elif self.current.spec == 0:
                axMult = 1000.0**self.current.axType
            width = (2*abs(float(self.rootwindow.paramframe.posEntries[self.pickNum].text())-pos[1]))/axMult
            self.rootwindow.paramframe.ampEntries[self.pickNum].setText("%.3g" %(float(self.rootwindow.paramframe.ampEntries[self.pickNum].text())*width))
            self.rootwindow.paramframe.lorEntries[self.pickNum].setText("%.3g" % abs(width))
            self.pickNum += 1
            self.pickWidth = False
        else:
            self.rootwindow.paramframe.posEntries[self.pickNum].setText("%.3g" %pos[1])
            left = pos[0] - 10 
            if left < 0:
                left = 0
            right = pos[0] + 10
            if right >= len(self.data1D) :
                right = len(self.data1D)-1
            self.rootwindow.paramframe.ampEntries[self.pickNum].setText("%.3g" %(pos[2]*0.5*np.pi))
            if self.pickNum < 10:
                self.rootwindow.paramframe.numExp.setCurrentIndex(self.pickNum)
                self.rootwindow.paramframe.changeNum()
            self.pickWidth = True
        if self.pickNum < 10:
            self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
            self.peakPick = True 

#################################################################################
class PeakDeconvParamFrame(QtGui.QWidget): 
    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        if self.parent.current.spec == 1:
            self.axAdd = self.parent.current.freq-self.parent.current.ref
            if self.parent.current.ppm:
                self.axMult = 1e6/self.parent.current.ref
            else:
                self.axMult = 1.0/(1000.0**self.parent.current.axType)
        elif self.parent.current.spec == 0:
            self.axMult = 1000.0**self.parent.current.axType
            self.axAdd = 0
        self.frame1 = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1,0,0)
        grid.addLayout(self.frame2,0,1)
        grid.addLayout(self.frame3,0,2)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton,0,0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton,1,0)
        cancelButton = QtGui.QPushButton("Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton,3,0)
        resetButton = QtGui.QPushButton("Reset")
        resetButton.clicked.connect(self.reset)
        self.frame1.addWidget(resetButton,0,1)
        self.pickTick = QtGui.QCheckBox("Pick")
        self.pickTick.stateChanged.connect(self.togglePick)
        self.frame1.addWidget(self.pickTick,1,1)
        self.frame1.setColumnStretch(10,1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel("Bgrnd:"),0,0,1,2)
        self.bgrndTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.bgrndTick,1,0)
        self.bgrndEntry = QtGui.QLineEdit()
        self.bgrndEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.bgrndEntry.setText("0.0")
        self.frame2.addWidget(self.bgrndEntry,1,1)
        self.frame2.addWidget(QLabel("Slope:"),2,0,1,2)
        self.slopeTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.slopeTick,3,0)
        self.slopeEntry = QtGui.QLineEdit()
        self.slopeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.slopeEntry.setText("0.0")
        self.frame2.addWidget(self.slopeEntry,3,1)
        self.frame2.setColumnStretch(10,1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1','2','3','4','5','6','7','8','9','10'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp,0,0,1,2)
        self.frame3.addWidget(QLabel("Position:"),1,0,1,2)
        self.frame3.addWidget(QLabel("Amplitude:"),1,2,1,2)
        self.frame3.addWidget(QLabel("Lorentz [Hz]:"),1,4,1,2)
        self.frame3.addWidget(QLabel("Gauss [Hz]:"),1,6,1,2)
        self.frame3.setColumnStretch(20,1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.posTicks = []
        self.posEntries = []
        self.ampTicks = []
        self.ampEntries = []
        self.lorTicks = []
        self.lorEntries = []
        self.gaussTicks = []
        self.gaussEntries = []
        for i in range(10):
            self.posTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.posTicks[i],i+2,0)
            self.posEntries.append(QtGui.QLineEdit())
            self.posEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.posEntries[i],i+2,1)
            self.ampTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.ampTicks[i],i+2,2)
            self.ampEntries.append(QtGui.QLineEdit())
            self.ampEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.ampEntries[i],i+2,3)
            self.lorTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.lorTicks[i],i+2,4)
            self.lorEntries.append(QtGui.QLineEdit())
            self.lorEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.lorEntries[i],i+2,5)
            self.gaussTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.gaussTicks[i],i+2,6)
            self.gaussEntries.append(QtGui.QLineEdit())
            self.gaussEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.gaussEntries[i],i+2,7)
        self.reset()
        grid.setColumnStretch(10,1)
        grid.setAlignment(QtCore.Qt.AlignLeft)

    def reset(self):
        self.parent.pickNum=0
        self.bgrndEntry.setText("0.0")
        self.bgrndTick.setChecked(True)
        self.slopeEntry.setText("0.0")
        self.slopeTick.setChecked(True)
        self.numExp.setCurrentIndex(0)
        self.pickTick.setChecked(True)
        for i in range(10):
            self.posEntries[i].setText("0.0")
            self.posTicks[i].setChecked(False)
            self.ampEntries[i].setText("1.0")
            self.ampTicks[i].setChecked(False)
            self.lorEntries[i].setText("1.0")
            self.lorTicks[i].setChecked(False)
            self.gaussEntries[i].setText("0.0")
            self.gaussTicks[i].setChecked(True)
            if i > 0:
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()
        self.togglePick()
        self.parent.pickWidth=False
        self.parent.showPlot()

    def changeNum(self,*args):
        val = self.numExp.currentIndex()+1
        for i in range(10):
            if i < val:
                self.posTicks[i].show()
                self.posEntries[i].show()
                self.ampTicks[i].show()
                self.ampEntries[i].show()
                self.lorTicks[i].show()
                self.lorEntries[i].show()
                self.gaussTicks[i].show()
                self.gaussEntries[i].show()
            else:
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()

    def togglePick(self):
        self.parent.togglePick(self.pickTick.isChecked())
                
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
            timeSignal = np.exp(1j*2*np.pi*t*(pos/self.axMult-self.axAdd))/len(x)*np.exp(-np.pi*width*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
            testFunc += amp*np.real(np.fft.fftshift(np.fft.fft(timeSignal)))
        testFunc += bgrnd+slope*x
        return testFunc

    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.bgrndEntry.text())
        if inp is None:
            return False
        self.bgrndEntry.setText('%.3g' % inp)
        inp = safeEval(self.slopeEntry.text())
        if inp is None:
            return False
        self.slopeEntry.setText('%.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.posEntries[i].text())
            if inp is None:
                return False
            self.posEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.ampEntries[i].text())
            if inp is None:
                return False
            self.ampEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.lorEntries[i].text())
            if inp is None:
                return False
            self.lorEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.gaussEntries[i].text())
            if inp is None:
                return False
            self.gaussEntries[i].setText('%.3g' % inp)
        return True
    
    def fit(self,*args):
        if not self.checkInputs():
            print("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = float(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            if not self.posTicks[i].isChecked():
                guess.append(float(self.posEntries[i].text()))
                struc.append(True)
            else:
                outPos[i] = float(self.posEntries[i].text())
                argu.append(outPos[i])
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].isChecked():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp,struc,argu)
        fitVal = scipy.optimize.curve_fit(self.fitFunc,self.parent.xax,self.parent.data1D,p0=guess)
        counter = 0
        if struc[0]:
            self.bgrndEntry.setText('%.3g' % fitVal[0][counter])
            outBgrnd = fitVal[0][counter]
            counter +=1
        if struc[1]:
            self.slopeEntry.setText('%.3g' % fitVal[0][counter])
            outSlope = fitVal[0][counter]
            counter +=1
        for i in range(numExp):
            if struc[4*i+2]:
                self.posEntries[i].setText('%.3g' % fitVal[0][counter])
                outPos[i] = fitVal[0][counter]
                counter += 1
            if struc[4*i+3]:
                self.ampEntries[i].setText('%.3g' % fitVal[0][counter])
                outAmp[i] = fitVal[0][counter]
                counter += 1
            if struc[4*i+4]:
                self.lorEntries[i].setText('%.3g' % abs(fitVal[0][counter]))
                outWidth[i] = abs(fitVal[0][counter])
                counter += 1
            if struc[4*i+5]:
                self.gaussEntries[i].setText('%.3g' % abs(fitVal[0][counter]))
                outGauss[i] = abs(fitVal[0][counter])
                counter += 1
        self.disp(outBgrnd, outSlope, outAmp, outPos, outWidth, outGauss)

    def sim(self):
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        outBgrnd = safeEval(self.bgrndEntry.text())
        outSlope = safeEval(self.slopeEntry.text())
        if outBgrnd is None or outSlope is None:
            print("One of the inputs is not valid")
            return
        for i in range(numExp):
            outPos[i] = safeEval(self.posEntries[i].text())
            outAmp[i] = safeEval(self.ampEntries[i].text())
            outWidth[i] = abs(safeEval(self.lorEntries[i].text()))
            outGauss[i] = abs(safeEval(self.gaussEntries[i].text()))
            if not np.isfinite([outPos[i],outAmp[i],outWidth[i],outGauss[i]]).all():
                print("One of the inputs is not valid")
                return
            self.lorEntries[i].setText('%.3g' % outWidth[i])
            self.gaussEntries[i].setText('%.3g' % outGauss[i])
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
            timeSignal = np.exp(1j*2*np.pi*t*(outPos[i]/self.axMult-self.axAdd))/len(tmpx)*np.exp(-np.pi*outWidth[i]*t)*np.exp(-((np.pi*outGauss[i]*t)**2)/(4*np.log(2)))
            y = outAmp[i]*np.real(np.fft.fftshift(np.fft.fft(timeSignal)))
            outCurvePart.append(outCurveBase + y)
            outCurve += y
        self.parent.showPlot(tmpx, outCurve, x, outCurvePart)

##############################################################################
class TensorDeconvWindow(QtGui.QWidget): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        QtGui.QWidget.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        #self.fig = self.mainProgram.getFig()
        self.fig = Figure()
        self.canvas =  FigureCanvas(self.fig)
        grid = QtGui.QGridLayout(self)
        grid.addWidget(self.canvas,0,0)
        self.current = TensorDeconvFrame(self,self.fig,self.canvas,oldMainWindow.current)
        self.paramframe = TensorDeconvParamFrame(self.current,self)
        grid.addWidget(self.paramframe,1,0)
        grid.setColumnStretch(0,1)
        grid.setRowStretch(0,1)
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
        self.show()

    def removeFromView(self):
        self.hide()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class TensorDeconvFrame(Plot1DFrame): 
    def __init__(self, rootwindow,fig,canvas,current):
        Plot1DFrame.__init__(self,rootwindow,fig,canvas)
        self.data1D = current.getDisplayedData()
        self.current = current
        self.spec = self.current.spec
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xax = self.current.xax*axMult
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
        self.peakPickReset()
        if var==1:
            self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
            self.peakPick = True
        else:
            self.peakPickFunc = None
            self.peakPick = False
        
    def pickDeconv(self, pos):
        if self.pickNum2 == 0:
            if self.pickNum < 10:
                self.rootwindow.paramframe.numExp.setCurrentIndex(self.pickNum)
                self.rootwindow.paramframe.changeNum()
            self.rootwindow.paramframe.t11Entries[self.pickNum].setText("%.2g" %self.current.xax[pos[0]])
            self.pickNum2 = 1
        elif self.pickNum2 == 1:
            self.rootwindow.paramframe.t22Entries[self.pickNum].setText("%.2g" %self.current.xax[pos[0]])
            self.pickNum2 = 2
        elif self.pickNum2 == 2:
            self.rootwindow.paramframe.t33Entries[self.pickNum].setText("%.2g" %self.current.xax[pos[0]])
            self.pickNum2 = 0
            self.pickNum += 1
        if self.pickNum < 10:
            self.peakPickFunc = lambda pos,self=self: self.pickDeconv(pos) 
            self.peakPick = True 

#################################################################################
class TensorDeconvParamFrame(QtGui.QWidget): 
    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        self.cheng = 15
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        if self.parent.current.spec == 1:
            self.axAdd = self.parent.current.freq-self.parent.current.ref
        elif self.parent.current.spec == 0:
            self.axAdd = 0
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1,0,0)
        grid.addLayout(self.optframe,0,1)
        grid.addLayout(self.frame2,0,2)
        grid.addLayout(self.frame3,0,3)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton,0,0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton,1,0)
        cancelButton = QtGui.QPushButton("Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton,3,0)
        resetButton = QtGui.QPushButton("Reset")
        resetButton.clicked.connect(self.reset)
        self.frame1.addWidget(resetButton,0,1)
        self.pickTick = QtGui.QCheckBox("Pick")
        self.pickTick.stateChanged.connect(self.togglePick)
        self.frame1.addWidget(self.pickTick,1,1)
        self.frame1.setColumnStretch(10,1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.optframe.addWidget(QLabel("Cheng"),0,0)
        self.chengEntry = QtGui.QLineEdit()
        self.chengEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.chengEntry.setText(str(self.cheng))
        self.optframe.addWidget(self.chengEntry,1,0)
        self.optframe.setColumnStretch(10,1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel("Bgrnd"),0,0,1,2)
        self.bgrndTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.bgrndTick,1,0)
        self.bgrndEntry = QtGui.QLineEdit()
        self.bgrndEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.bgrndEntry.setText("0.0")
        self.frame2.addWidget(self.bgrndEntry,1,1)
        self.frame2.addWidget(QLabel("Slope"),2,0,1,2)
        self.slopeTick = QtGui.QCheckBox('')
        self.frame2.addWidget(self.slopeTick,3,0)
        self.slopeEntry = QtGui.QLineEdit()
        self.slopeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.slopeEntry.setText("0.0")
        self.frame2.addWidget(self.slopeEntry,3,1)
        self.frame2.setColumnStretch(10,1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1','2','3','4','5','6','7','8','9','10'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp,0,0,1,2)
        self.frame3.addWidget(QLabel("T11"),1,0,1,2)
        self.frame3.addWidget(QLabel("T22"),1,2,1,2)
        self.frame3.addWidget(QLabel("T33"),1,4,1,2)
        self.frame3.addWidget(QLabel("Amplitude"),1,6,1,2)
        self.frame3.addWidget(QLabel("Lorentz [Hz]"),1,8,1,2)
        self.frame3.addWidget(QLabel("Gauss [Hz]"),1,10,1,2)
        self.frame3.setColumnStretch(20,1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.t11Entries = []
        self.t11Ticks = []
        self.t22Entries = []
        self.t22Ticks = []
        self.t33Entries = []
        self.t33Ticks = []
        self.ampEntries = []
        self.ampTicks = []
        self.lorEntries = []
        self.lorTicks = []
        self.gaussEntries = []
        self.gaussTicks = []
        for i in range(10):
            self.t11Ticks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.t11Ticks[i],i+2,0)
            self.t11Entries.append(QtGui.QLineEdit())
            self.t11Entries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.t11Entries[i],i+2,1)
            self.t22Ticks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.t22Ticks[i],i+2,2)
            self.t22Entries.append(QtGui.QLineEdit())
            self.t22Entries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.t22Entries[i],i+2,3)
            self.t33Ticks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.t33Ticks[i],i+2,4)
            self.t33Entries.append(QtGui.QLineEdit())
            self.t33Entries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.t33Entries[i],i+2,5)
            self.ampTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.ampTicks[i],i+2,6)
            self.ampEntries.append(QtGui.QLineEdit())
            self.ampEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.ampEntries[i],i+2,7)
            self.lorTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.lorTicks[i],i+2,8)
            self.lorEntries.append(QtGui.QLineEdit())
            self.lorEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.lorEntries[i],i+2,9)
            self.gaussTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.gaussTicks[i],i+2,10)
            self.gaussEntries.append(QtGui.QLineEdit())
            self.gaussEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.frame3.addWidget(self.gaussEntries[i],i+2,11)
        self.reset()
        grid.setColumnStretch(10,1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        
    def reset(self):
        self.parent.pickNum=0
        self.parent.pickNum2=0
        self.bgrndEntry.setText("0.0")
        self.bgrndTick.setChecked(True)
        self.slopeEntry.setText("0.0")
        self.slopeTick.setChecked(True)
        self.cheng = 15
        self.chengEntry.setText(str(self.cheng))
        self.pickTick.setChecked(True)
        for i in range(10):
            self.t11Entries[i].setText("0.0")
            self.t11Ticks[i].setChecked(False)
            self.t22Entries[i].setText("0.0")
            self.t22Ticks[i].setChecked(False)
            self.t33Entries[i].setText("0.0")
            self.t33Ticks[i].setChecked(False)
            self.ampEntries[i].setText("1.0")
            self.ampTicks[i].setChecked(False)
            self.lorEntries[i].setText("10.0")
            self.lorTicks[i].setChecked(True)
            self.gaussEntries[i].setText("0.0")
            self.gaussTicks[i].setChecked(True)
            if i > 0:
                self.t11Ticks[i].hide()
                self.t11Entries[i].hide()
                self.t22Ticks[i].hide()
                self.t22Entries[i].hide()
                self.t33Ticks[i].hide()
                self.t33Entries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()
        self.togglePick()
        self.parent.showPlot()

    def togglePick(self):
        self.parent.togglePick(self.pickTick.isChecked())
                
    def setCheng(self,*args):
        inp = safeEval(self.chengEntry.text())
        if inp is None:
            self.cheng = 15
        else:
            self.cheng = int(inp)
        self.chengEntry.setText(str(self.cheng))
                
    def changeNum(self,*args):
        val = self.numExp.currentIndex()+1
        for i in range(10):
            if i < val:
                self.t11Ticks[i].show()
                self.t11Entries[i].show()
                self.t22Ticks[i].show()
                self.t22Entries[i].show()
                self.t33Ticks[i].show()
                self.t33Entries[i].show()
                self.ampTicks[i].show()
                self.ampEntries[i].show()
                self.lorTicks[i].show()
                self.lorEntries[i].show()
                self.gaussTicks[i].show()
                self.gaussEntries[i].show()
            else:
                self.t11Ticks[i].hide()
                self.t11Entries[i].hide()
                self.t22Ticks[i].hide()
                self.t22Entries[i].hide()
                self.t33Ticks[i].hide()
                self.t33Entries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()

    def tensorFunc(self, x, t11, t22, t33, lor, gauss):
        t11=t11*self.multt11
        t22=t22*self.multt22
        t33=t33*self.multt33
        v=t11+t22+t33-self.axAdd
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

    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.bgrndEntry.text())
        if inp is None:
            return False
        self.bgrndEntry.setText('%.3g' % inp)
        inp = safeEval(self.slopeEntry.text())
        if inp is None:
            return False
        self.slopeEntry.setText('%.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.t11Entries[i].text())
            if inp is None:
                return False
            self.t11Entries[i].setText('%.3g' % inp)
            inp = safeEval(self.t22Entries[i].text())
            if inp is None:
                return False
            self.t22Entries[i].setText('%.3g' % inp)
            inp = safeEval(self.t33Entries[i].text())
            if inp is None:
                return False
            self.t33Entries[i].setText('%.3g' % inp)
            inp = safeEval(self.ampEntries[i].text())
            if inp is None:
                return False
            self.ampEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.lorEntries[i].text())
            if inp is None:
                return False
            self.lorEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.gaussEntries[i].text())
            if inp is None:
                return False
            self.gaussEntries[i].setText('%.3g' % inp)
        return True
        
    def fit(self,*args):
        self.setCheng()
        if not self.checkInputs():
            print("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outt11 = np.zeros(numExp)
        outt22 = np.zeros(numExp)
        outt33 = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = float(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            if not self.t11Ticks[i].isChecked():
                guess.append(float(self.t11Entries[i].text()))
                struc.append(True)
            else:
                outt11[i] = float(self.t11Entries[i].text())
                argu.append(outt11[i])
                struc.append(False)
            if not self.t22Ticks[i].isChecked():
                guess.append(float(self.t22Entries[i].text()))
                struc.append(True)
            else:
                outt22[i] = float(self.t22Entries[i].text())
                argu.append(outt22[i])
                struc.append(False)
            if not self.t33Ticks[i].isChecked():
                guess.append(float(self.t33Entries[i].text()))
                struc.append(True)
            else:
                outt33[i] = float(self.t33Entries[i].text())
                argu.append(outt33[i])
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].isChecked():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp,struc,argu)
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.multt11=np.sin(theta)**2*np.cos(phi)**2
        self.multt22=np.sin(theta)**2*np.sin(phi)**2
        self.multt33=np.cos(theta)**2
        fitVal = scipy.optimize.fmin(self.fitFunc,guess, args=(self.parent.xax,np.real(self.parent.data1D)))
        counter = 0
        if struc[0]:
            self.bgrndEntry.setText('%.2g' % fitVal[counter])
            outBgrnd = fitVal[counter]
            counter +=1
        if struc[1]:
            self.slopeEntry.setText('%.2g' % fitVal[counter])
            outSlope = fitVal[counter]
            counter +=1
        for i in range(numExp):
            if struc[6*i+2]:
                self.t11Entries[i].setText('%.2g' % fitVal[counter])
                outt11[i] = fitVal[counter]
                counter += 1
            if struc[6*i+3]:
                self.t22Entries[i].setText('%.2g' % fitVal[counter])
                outt22[i] = fitVal[counter]
                counter += 1
            if struc[6*i+4]:
                self.t33Entries[i].setText('%.2g' % fitVal[counter])
                outt33[i] = fitVal[counter]
                counter += 1
            if struc[6*i+5]:
                self.ampEntries[i].setText('%.2g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[6*i+6]:
                self.lorEntries[i].setText('%.2g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
            if struc[6*i+7]:
                self.gaussEntries[i].setText('%.2g' % abs(fitVal[counter]))
                outGauss[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd,outSlope,outt11,outt22,outt33,outAmp,outWidth,outGauss)

    def sim(self):
        self.setCheng()
        numExp = self.numExp.currentIndex()+1
        bgrnd = safeEval(self.bgrndEntry.text())
        slope = safeEval(self.slopeEntry.text())
        if bgrnd is None or slope is None:
            print("One of the inputs is not valid")
            return
        t11 = np.zeros(numExp)
        t22 = np.zeros(numExp)
        t33 = np.zeros(numExp)
        amp = np.zeros(numExp)
        width = np.zeros(numExp)
        gauss = np.zeros(numExp)
        for i in range(numExp):
            t11[i] = safeEval(self.t11Entries[i].text())
            t22[i] = safeEval(self.t22Entries[i].text())
            t33[i] = safeEval(self.t33Entries[i].text())
            amp[i] = safeEval(self.ampEntries[i].text())
            width[i] = safeEval(self.lorEntries[i].text())
            gauss[i] = safeEval(self.gaussEntries[i].text())
            if not np.isfinite([t11[i], t22[i], t33[i], amp[i], width[i], gauss[i]]).all():
                print("One of the inputs is not valid")
                return
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.multt11=np.sin(theta)**2*np.cos(phi)**2
        self.multt22=np.sin(theta)**2*np.sin(phi)**2
        self.multt33=np.cos(theta)**2
        self.disp(bgrnd,slope,t11,t22,t33,amp,width,gauss)
        
##############################################################################
class Quad1DeconvWindow(QtGui.QWidget): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow):
        QtGui.QWidget.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        #self.fig = self.mainProgram.getFig()
        self.fig = Figure()
        self.canvas =  FigureCanvas(self.fig)
        grid = QtGui.QGridLayout(self)
        grid.addWidget(self.canvas,0,0)
        self.current = Quad1DeconvFrame(self,self.fig,self.canvas,oldMainWindow.current)
        self.paramframe = Quad1DeconvParamFrame(self.current,self)
        grid.addWidget(self.paramframe,1,0)
        grid.setColumnStretch(0,1)
        grid.setRowStretch(0,1)
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
        self.show()

    def removeFromView(self):
        self.hide()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)
        
#################################################################################   
class Quad1DeconvFrame(Plot1DFrame):
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
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.xminlim=min(self.xax*axMult)
        self.xmaxlim=max(self.xax*axMult)
        if self.spec > 0 :
            self.ax.set_xlim(self.xmaxlim,self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim,self.xmaxlim)
        self.ax.set_ylim(self.yminlim,self.ymaxlim)
        
    def showPlot(self, tmpAx=None, tmpdata=None, tmpAx2=[], tmpdata2=[]): 
        self.ax.cla()
        if self.spec == 1:
            if self.current.ppm:
                axMult = 1e6/self.current.ref
            else:
                axMult = 1.0/(1000.0**self.current.axType)
        elif self.spec == 0:
            axMult = 1000.0**self.current.axType
        self.line_xdata = self.xax*axMult
        self.line_ydata = self.data1D
        self.ax.plot(self.xax*axMult,self.data1D)
        if tmpAx is not None:
            self.ax.plot(tmpAx*axMult,tmpdata)
        for i in range(len(tmpAx2)):
            self.ax.plot(tmpAx2[i]*axMult,tmpdata2[i])
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
class Quad1DeconvParamFrame(QtGui.QWidget):
    
    Ioptions = ['1','3/2','2','5/2','3','7/2','4','9/2']
    
    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        self.cheng = 15
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        if self.parent.current.spec == 1:
            self.axAdd = self.parent.current.freq-self.parent.current.ref
        elif self.parent.current.spec == 0:
            self.axAdd = 0
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1,0,0)
        grid.addLayout(self.optframe,0,1)
        grid.addLayout(self.frame2,0,2)
        grid.addLayout(self.frame3,0,3)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton,0,0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton,1,0)
        cancelButton = QtGui.QPushButton("Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton,3,0)
        self.frame1.setColumnStretch(10,1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.optframe.addWidget(QLabel("Cheng:"),0,0)
        self.chengEntry = QtGui.QLineEdit()
        self.chengEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.chengEntry.setText(str(self.cheng))
        self.optframe.addWidget(self.chengEntry,1,0)
        self.optframe.addWidget(QLabel("I:"),0,1)
        self.IEntry = QtGui.QComboBox()
        self.IEntry.addItems(self.Ioptions)
        self.IEntry.setCurrentIndex(1)
        self.optframe.addWidget(self.IEntry,1,1)
        self.optframe.setColumnStretch(10,1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel("Bgrnd:"),0,0,1,2)
        self.bgrndTick = QtGui.QCheckBox('')
        self.bgrndTick.setChecked(True)
        self.frame2.addWidget(self.bgrndTick,1,0)
        self.bgrndEntry = QtGui.QLineEdit()
        self.bgrndEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.bgrndEntry.setText("0.0")
        self.frame2.addWidget(self.bgrndEntry,1,1)
        self.frame2.addWidget(QLabel("Slope:"),2,0,1,2)
        self.slopeTick = QtGui.QCheckBox('')
        self.slopeTick.setChecked(True)
        self.frame2.addWidget(self.slopeTick,3,0)
        self.slopeEntry = QtGui.QLineEdit()
        self.slopeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.slopeEntry.setText("0.0")
        self.frame2.addWidget(self.slopeEntry,3,1)
        self.frame2.setColumnStretch(10,1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1','2','3','4','5','6','7','8','9','10'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp,0,0,1,2)
        self.frame3.addWidget(QLabel("Pos:"),1,0,1,2)
        self.frame3.addWidget(QLabel("Cq [MHz]:"),1,2,1,2)
        self.frame3.addWidget(QLabel(u"\u03b7:"),1,4,1,2)
        self.frame3.addWidget(QLabel("Amplitude:"),1,6,1,2)
        self.frame3.addWidget(QLabel("Lorentz [Hz]:"),1,8,1,2)
        self.frame3.addWidget(QLabel("Gauss [Hz]:"),1,10,1,2)
        self.frame3.setColumnStretch(20,1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.posEntries = []
        self.posTicks = []
        self.cqEntries = []
        self.cqTicks = []
        self.etaEntries = []
        self.etaTicks = []
        self.ampEntries = []
        self.ampTicks = []
        self.lorEntries = []
        self.lorTicks = []
        self.gaussEntries = []
        self.gaussTicks = []
        for i in range(10):
            self.posTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.posTicks[i],i+2,0)
            self.posEntries.append(QtGui.QLineEdit())
            self.posEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.posEntries[i].setText("0.0")
            self.frame3.addWidget(self.posEntries[i],i+2,1)
            self.cqTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.cqTicks[i],i+2,2)
            self.cqEntries.append(QtGui.QLineEdit())
            self.cqEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.cqEntries[i].setText("0.0")
            self.frame3.addWidget(self.cqEntries[i],i+2,3)
            self.etaTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.etaTicks[i],i+2,4)
            self.etaEntries.append(QtGui.QLineEdit())
            self.etaEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.etaEntries[i].setText("0.0")
            self.frame3.addWidget(self.etaEntries[i],i+2,5)
            self.ampTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.ampTicks[i],i+2,6)
            self.ampEntries.append(QtGui.QLineEdit())
            self.ampEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.ampEntries[i].setText("1.0")
            self.frame3.addWidget(self.ampEntries[i],i+2,7)
            self.lorTicks.append(QtGui.QCheckBox(''))
            self.lorTicks[i].setChecked(True)
            self.frame3.addWidget(self.lorTicks[i],i+2,8)
            self.lorEntries.append(QtGui.QLineEdit())
            self.lorEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.lorEntries[i].setText("10.0")
            self.frame3.addWidget(self.lorEntries[i],i+2,9)
            self.gaussTicks.append(QtGui.QCheckBox(''))
            self.gaussTicks[i].setChecked(True)
            self.frame3.addWidget(self.gaussTicks[i],i+2,10)
            self.gaussEntries.append(QtGui.QLineEdit())
            self.gaussEntries[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.gaussEntries[i].setText("0.0")
            self.frame3.addWidget(self.gaussEntries[i],i+2,11)

            if i > 0:
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.cqTicks[i].hide()
                self.cqEntries[i].hide()
                self.etaTicks[i].hide()
                self.etaEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()
        grid.setColumnStretch(10,1)
        grid.setAlignment(QtCore.Qt.AlignLeft)
        
    def checkI(self,I):
        return I*0.5+1
                
    def setCheng(self,*args):
        inp = safeEval(self.chengEntry.text())
        if inp is None:
            print("One of the inputs is not valid")
            return
        self.cheng = int(inp)
        self.chengEntry.setText(str(self.cheng))
                
    def changeNum(self,*args):
        val = self.numExp.currentIndex()+1
        for i in range(10):
            if i < val:
                self.posTicks[i].show()
                self.posEntries[i].show()
                self.cqTicks[i].show()
                self.cqEntries[i].show()
                self.etaTicks[i].show()
                self.etaEntries[i].show()
                self.ampTicks[i].show()
                self.ampEntries[i].show()
                self.lorTicks[i].show()
                self.lorEntries[i].show()
                self.gaussTicks[i].show()
                self.gaussEntries[i].show()
            else:
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.cqTicks[i].hide()
                self.cqEntries[i].hide()
                self.etaTicks[i].hide()
                self.etaEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()

    def tensorFunc(self, x, I, pos, cq, eta, width, gauss):
        m=np.arange(-I,I)
        v=[]
        weights=[]
        pos = pos - self.axAdd
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
        
    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.bgrndEntry.text())
        if inp is None:
            return False
        self.bgrndEntry.setText('%.3g' % inp)
        inp = safeEval(self.slopeEntry.text())
        if inp is None:
            return False
        self.slopeEntry.setText('%.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.posEntries[i].text())
            if inp is None:
                return False
            self.posEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.cqEntries[i].text())
            if inp is None:
                return False
            self.cqEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.etaEntries[i].text())
            if inp is None:
                return False
            self.etaEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.ampEntries[i].text())
            if inp is None:
                return False
            self.ampEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.lorEntries[i].text())
            if inp is None:
                return False
            self.lorEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.gaussEntries[i].text())
            if inp is None:
                return False
            self.gaussEntries[i].setText('%.3g' % inp)
        return True
    
    def fit(self,*args):
        self.setCheng()
        if not self.checkInputs():
            print("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outCq = np.zeros(numExp)
        outEta = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        I = self.checkI(self.IEntry.currentIndex())
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = float(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            if not self.posTicks[i].isChecked():
                guess.append(float(self.posEntries[i].text()))
                struc.append(True)
            else:
                outPos[i] = float(self.posEntries[i].text())
                argu.append(outPos[i])
                struc.append(False)
            if not self.cqTicks[i].isChecked():
                guess.append(float(self.cqEntries[i].text())*1e6)
                struc.append(True)
            else:
                outCq[i] = float(self.cqEntries[i].text())*1e6
                argu.append(outCq[i])
                struc.append(False)
            if not self.etaTicks[i].isChecked():
                guess.append(float(self.etaEntries[i].text()))
                struc.append(True)
            else:
                outEta[i] = float(self.etaEntries[i].text())
                argu.append(outEta[i])
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].isChecked():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp,struc,argu,I)
        self.setAngleStuff()
        fitVal = scipy.optimize.fmin(self.fitFunc,guess, args=(self.parent.xax,np.real(self.parent.data1D)))
        counter = 0
        if struc[0]:
            self.bgrndEntry.setText('%.2g' % fitVal[counter])
            outBgrnd = fitVal[counter]
            counter +=1
        if struc[1]:
            self.slopeEntry.setText('%.2g' % fitVal[counter])
            outSlope = fitVal[counter]
            counter +=1
        for i in range(numExp):
            if struc[6*i+2]:
                self.posEntries[i].setText('%.2g' % fitVal[counter])
                outPos[i] = fitVal[counter]
                counter += 1
            if struc[6*i+3]:
                self.cqEntries[i].setText('%.2g' % (fitVal[counter]*1e-6))
                outCq[i] = fitVal[counter]
                counter += 1
            if struc[6*i+4]:
                self.etaEntries[i].setText('%.2g' % abs(fitVal[counter]))
                outEta[i] = fitVal[counter]
                counter += 1
            if struc[6*i+5]:
                self.ampEntries[i].setText('%.2g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[6*i+6]:
                self.lorEntries[i].setText('%.2g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
            if struc[6*i+7]:
                self.gaussEntries[i].setText('%.2g' % abs(fitVal[counter]))
                outGauss[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd,outSlope,I,outPos,outCq,outEta,outAmp,outWidth,outGauss)

    def sim(self):
        self.setCheng()
        if not self.checkInputs():
            print("One of the inputs is not valid")
            return
        numExp = self.numExp.currentIndex()
        bgrnd = float(self.bgrndEntry.text())
        slope = float(self.slopeEntry.text())
        pos = np.zeros(numExp)
        cq = np.zeros(numExp)
        eta = np.zeros(numExp)
        amp = np.zeros(numExp)
        width = np.zeros(numExp)
        gauss = np.zeros(numExp)
        I = self.checkI(self.IEntry.currentIndex())
        for i in range(numExp):
            pos[i] = float(self.posEntries[i].text())
            cq[i] = float(self.cqEntries[i].text())*1e6
            eta[i] = float(self.etaEntries[i].text())
            amp[i] = float(self.ampEntries[i].text())
            width[i] = float(self.lorEntries[i].text())
            gauss[i] = float(self.gaussEntries[i].text())
        self.setAngleStuff()
        self.disp(bgrnd,slope,I,pos,cq,eta,amp,width,gauss)

##############################################################################
class Quad2DeconvWindow(QtGui.QWidget): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow,mas=False):
        QtGui.QWidget.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        #self.fig = self.mainProgram.getFig()
        self.fig = Figure()
        self.canvas =  FigureCanvas(self.fig)
        grid = QtGui.QGridLayout(self)
        grid.addWidget(self.canvas,0,0)
        self.current = Quad1DeconvFrame(self,self.fig,self.canvas,oldMainWindow.current)
        if mas:
            self.paramframe = Quad2MASDeconvParamFrame(self.current,self)
        else:
            self.paramframe = Quad2StaticDeconvParamFrame(self.current,self)
        grid.addWidget(self.paramframe,1,0)
        grid.setColumnStretch(0,1)
        grid.setRowStretch(0,1)
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
        self.show()

    def removeFromView(self):
        self.hide()

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
        
    def checkI(self,I):
        return I*1.0+1.5
    
    def setAngleStuff(self):
        phi,theta,self.weight = zcw_angles(self.cheng,symm=2)
        self.angleStuff1 = -27/8.0*np.cos(theta)**4+15/4.0*np.cos(theta)**2-3/8.0
        self.angleStuff2 = (-9/4.0*np.cos(theta)**4+2*np.cos(theta)**2+1/4.0)*np.cos(2*phi)
        self.angleStuff3 = -1/2.0*np.cos(theta)**2+1/3.0+(-3/8.0*np.cos(theta)**4+3/4.0*np.cos(theta)**2-3/8.0)*np.cos(2*phi)**2

    def tensorFunc(self, x, I, pos, cq, eta, width, gauss):
        pos = pos - self.axAdd
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
class Quad2CzjzekWindow(QtGui.QWidget): 
    def __init__(self, rootwindow,mainProgram,oldMainWindow,mas=False):
        QtGui.QWidget.__init__(self,rootwindow)
        self.mainProgram = mainProgram
        self.oldMainWindow = oldMainWindow
        self.name = self.oldMainWindow.name
        self.fig = Figure()
        self.canvas =  FigureCanvas(self.fig)
        grid = QtGui.QGridLayout(self)
        grid.addWidget(self.canvas,0,0)
        self.current = Quad1DeconvFrame(self,self.fig,self.canvas,oldMainWindow.current)
        if mas:
            self.paramframe = Quad2MASCzjzekParamFrame(self.current,self)
        else:
            self.paramframe = Quad2StaticCzjzekParamFrame(self.current,self)
        grid.addWidget(self.paramframe,1,0)
        grid.setColumnStretch(0,1)
        grid.setRowStretch(0,1)
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
        self.show()

    def removeFromView(self):
        self.hide()

    def kill(self):
        self.current.kill()
        self.oldMainWindow.kill()
        del self.current
        
    def cancel(self):
         self.mainProgram.closeFitWindow(self.oldMainWindow)

#################################################################################
class Quad2StaticCzjzekParamFrame(QtGui.QWidget): 

    Ioptions = ['3/2','5/2','7/2','9/2']

    def __init__(self, parent, rootwindow):
        QtGui.QWidget.__init__(self, rootwindow)
        self.parent = parent
        self.rootwindow = rootwindow
        self.cheng = 15
        grid = QtGui.QGridLayout(self)
        self.setLayout(grid)
        if self.parent.current.spec == 1:
            self.axAdd = self.parent.current.freq-self.parent.current.ref
        elif self.parent.current.spec == 0:
            self.axAdd = 0
        self.frame1 = QtGui.QGridLayout()
        self.optframe = QtGui.QGridLayout()
        self.frame2 = QtGui.QGridLayout()
        self.frame3 = QtGui.QGridLayout()
        grid.addLayout(self.frame1,0,0)
        grid.addLayout(self.optframe,0,1)
        grid.addLayout(self.frame2,0,2)
        grid.addLayout(self.frame3,0,3)
        simButton = QtGui.QPushButton("Sim")
        simButton.clicked.connect(self.sim)
        self.frame1.addWidget(simButton,0,0)
        fitButton = QtGui.QPushButton("Fit")
        fitButton.clicked.connect(self.fit)
        self.frame1.addWidget(fitButton,1,0)
        cancelButton = QtGui.QPushButton("Cancel")
        cancelButton.clicked.connect(rootwindow.cancel)
        self.frame1.addWidget(cancelButton,3,0)
        self.frame1.setColumnStretch(10,1)
        self.frame1.setAlignment(QtCore.Qt.AlignTop)
        self.optframe.addWidget(QLabel("Cheng:"),0,0)
        self.chengEntry = QtGui.QLineEdit()
        self.chengEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.chengEntry.setText(str(self.cheng))
        self.optframe.addWidget(self.chengEntry,1,0)
        self.optframe.addWidget(QLabel("I"),0,1)
        self.IEntry = QtGui.QComboBox()
        self.IEntry.addItems(self.Ioptions)
        self.IEntry.setCurrentIndex(1)
        self.optframe.addWidget(self.IEntry,1,1)
        self.optframe.addWidget(QLabel("Cq grid size:"),2,0)
        self.cqGridEntry = QtGui.QLineEdit()
        self.cqGridEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.cqGridEntry.setText("50")
        self.cqGridEntry.returnPressed.connect(self.setGrid)
        self.optframe.addWidget(self.cqGridEntry,3,0)
        self.optframe.addWidget(QLabel(u"\u03b7 grid size:"),4,0)
        self.etaGridEntry = QtGui.QLineEdit()
        self.etaGridEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.etaGridEntry.setText("10")
        self.etaGridEntry.returnPressed.connect(self.setGrid)
        self.optframe.addWidget(self.etaGridEntry,5,0)
        self.optframe.setColumnStretch(10,1)
        self.optframe.setAlignment(QtCore.Qt.AlignTop)
        self.frame2.addWidget(QLabel("Bgrnd:"),0,0,1,2)
        self.bgrndTick = QtGui.QCheckBox('')
        self.bgrndTick.setChecked(True)
        self.frame2.addWidget(self.bgrndTick,1,0)
        self.bgrndEntry = QtGui.QLineEdit()
        self.bgrndEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.bgrndEntry.setText("0.0")
        self.frame2.addWidget(self.bgrndEntry,1,1)
        self.frame2.addWidget(QLabel("Slope:"),2,0,1,2)
        self.slopeTick = QtGui.QCheckBox('')
        self.slopeTick.setChecked(True)
        self.frame2.addWidget(self.slopeTick,3,0)
        self.slopeEntry = QtGui.QLineEdit()
        self.slopeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.slopeEntry.setText("0.0")
        self.frame2.addWidget(self.slopeEntry,3,1)
        self.frame2.setColumnStretch(10,1)
        self.frame2.setAlignment(QtCore.Qt.AlignTop)
        self.numExp = QtGui.QComboBox()
        self.numExp.addItems(['1','2','3','4','5','6','7','8','9','10'])
        self.numExp.currentIndexChanged.connect(self.changeNum)
        self.frame3.addWidget(self.numExp,0,0)
        self.frame3.addWidget(QLabel("d:"),1,0)
        self.frame3.addWidget(QLabel("Pos:"),1,1,1,2)
        self.frame3.addWidget(QLabel(u"\u03c3 [MHz]:"),1,3,1,2)
        self.frame3.addWidget(QLabel("Amplitude:"),1,5,1,2)
        self.frame3.addWidget(QLabel("Lorentz [Hz]:"),1,7,1,2)
        self.frame3.addWidget(QLabel("Gauss [Hz]:"),1,9,1,2)
        self.frame3.setColumnStretch(20,1)
        self.frame3.setAlignment(QtCore.Qt.AlignTop)
        self.posEntries = []
        self.posTicks = []
        self.sigmaEntries = []
        self.sigmaTicks = []
        self.dEntries = []
        self.ampEntries = []
        self.ampTicks = []
        self.lorEntries = []
        self.lorTicks = []
        self.gaussEntries = []
        self.gaussTicks = []
        for i in range(10):
            self.dEntries.append(QtGui.QLineEdit())
            self.dEntries[i].setText("5")
            self.frame3.addWidget(self.dEntries[i],i+2,0)
            self.posTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.posTicks[i],i+2,1)
            self.posEntries.append(QtGui.QLineEdit())
            self.posEntries[i].setText("0.0")
            self.frame3.addWidget(self.posEntries[i],i+2,2)
            self.sigmaTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.sigmaTicks[i],i+2,3)
            self.sigmaEntries.append(QtGui.QLineEdit())
            self.sigmaEntries[i].setText("1.0")
            self.frame3.addWidget(self.sigmaEntries[i],i+2,4)
            self.ampTicks.append(QtGui.QCheckBox(''))
            self.frame3.addWidget(self.ampTicks[i],i+2,5)
            self.ampEntries.append(QtGui.QLineEdit())
            self.ampEntries[i].setText("1.0")
            self.frame3.addWidget(self.ampEntries[i],i+2,6)
            self.lorTicks.append(QtGui.QCheckBox(''))
            self.lorTicks[i].setChecked(True)
            self.frame3.addWidget(self.lorTicks[i],i+2,7)
            self.lorEntries.append(QtGui.QLineEdit())
            self.lorEntries[i].setText("10.0")
            self.frame3.addWidget(self.lorEntries[i],i+2,8)
            self.gaussTicks.append(QtGui.QCheckBox(''))
            self.gaussTicks[i].setChecked(True)
            self.frame3.addWidget(self.gaussTicks[i],i+2,9)
            self.gaussEntries.append(QtGui.QLineEdit())
            self.gaussEntries[i].setText("0.0")
            self.frame3.addWidget(self.gaussEntries[i],i+2,10)

            if i > 0:
                self.dEntries[i].hide()
                self.posTicks[i].hide()
                self.posEntries[i].hide()
                self.sigmaTicks[i].hide()
                self.sigmaEntries[i].hide()
                self.ampTicks[i].hide()
                self.ampEntries[i].hide()
                self.lorTicks[i].hide()
                self.lorEntries[i].hide()
                self.gaussTicks[i].hide()
                self.gaussEntries[i].hide()
        grid.setColumnStretch(10,1)
        grid.setAlignment(QtCore.Qt.AlignLeft)

    def checkI(self,I):
        return I*1.0+1.5
                
    def setCheng(self,*args):
        inp = safeEval(self.chengEntry.text())
        if inp is None:
            self.cheng = 15
        else:
            self.cheng = int(inp)
        self.chengEntry.setText(str(self.cheng))

    def setGrid(self, *args):
        inp = safeEval(self.cqGridEntry.text())
        if inp is not None:
            self.cqGridEntry(str(int(inp)))
        inp = safeEval(self.etaGridEntry.text())
        if inp is not None:
            self.etaGridEntry(str(int(inp)))
        
    def changeNum(self,*args):
        val = self.numExp.currentIndex()+1
        for i in range(10):
            if i < val:
                self.posCheck[i].show()
                self.posEntries[i].show()
                self.sigmaCheck[i].show()
                self.sigmaEntries[i].show()
                self.dEntries[i].show()
                self.ampCheck[i].show()
                self.ampEntries[i].show()
                self.lorCheck[i].show()
                self.lorEntries[i].show()
                self.gaussCheck[i].show()
                self.gaussEntries[i].show()
            else:
                self.posCheck[i].hide()
                self.posEntries[i].hide()
                self.sigmaCheck[i].hide()
                self.sigmaEntries[i].hide()
                self.dEntries[i].hide()
                self.ampCheck[i].hide()
                self.ampEntries[i].hide()
                self.lorCheck[i].hide()
                self.lorEntries[i].hide()
                self.gaussCheck[i].hide()
                self.gaussEntries[i].hide()

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
        pos = pos - self.axAdd
        cq = self.cq
        eta = self.eta
        czjzek = cq**(d-1)*eta/(np.sqrt(2*np.pi)*sigma)*(1-eta**2/9.0)*np.exp(-cq**2/(2.0*sigma**2)*(1+eta**2/3.0))
        czjzek = czjzek/np.sum(czjzek)
        fid = np.sum(self.lib*czjzek[...,None],axis=(0,1))
        t=np.arange(len(fid))/self.parent.current.sw
        apod = np.exp(-np.pi*width*t)*np.exp(-((np.pi*gauss*t)**2)/(4*np.log(2)))
        apod[-1:-(len(apod)/2+1):-1]=apod[:len(apod)/2]
        return scipy.ndimage.interpolation.shift(np.real(np.fft.fft(fid*apod)),len(fid)*pos/self.parent.current.sw)
    
    def checkInputs(self):
        numExp = self.numExp.currentIndex()+1
        inp = safeEval(self.bgrndEntry.text())
        if inp is None:
            return False
        self.bgrndEntry.setText('%.3g' % inp)
        inp = safeEval(self.slopeEntry.text())
        if inp is None:
            return False
        self.slopeEntry.setText('%.3g' % inp)
        for i in range(numExp):
            inp = safeEval(self.dEntries[i].text())
            if inp is None:
                return False
            if inp < 1:
                inp = 1
            elif inp > 5:
                inp = 5
            self.dEntries[i].setText(str(int(inp)))
            inp = safeEval(self.posEntries[i].text())
            if inp is None:
                return False
            self.posEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.sigmaEntries[i].text())
            if inp is None:
                return False
            self.sigmaEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.ampEntries[i].text())
            if inp is None:
                return False
            self.ampEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.lorEntries[i].text())
            if inp is None:
                return False
            self.lorEntries[i].setText('%.3g' % inp)
            inp = safeEval(self.gaussEntries[i].text())
            if inp is None:
                return False
            self.gaussEntries[i].setText('%.3g' % inp)
        return True
    
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
        if not self.setGrid():
            print("One of the inputs is not valid")
            return
        if not self.checkInputs():
            print("One of the inputs is not valid")
            return
        struc = []
        guess = []
        argu = []
        maxCq = 0.0
        I = self.checkI(self.IEntry.currentIndex())
        numExp = self.numExp.currentIndex()+1
        outPos = np.zeros(numExp)
        outSigma = np.zeros(numExp)
        outD = np.zeros(numExp)
        outAmp = np.zeros(numExp)
        outWidth = np.zeros(numExp)
        outGauss = np.zeros(numExp)
        if not self.bgrndTick.isChecked():
            guess.append(float(self.bgrndEntry.text()))
            struc.append(True)
        else:
            outBgrnd = float(self.bgrndEntry.text())
            argu.append(outBgrnd)
            struc.append(False)
        if not self.slopeTick.isChecked():
            guess.append(float(self.slopeEntry.text()))
            struc.append(True)
        else:
            outSlope = safeEval(self.slopeEntry.text())
            argu.append(outSlope)
            struc.append(False)
        for i in range(numExp):
            inp = int(self.dEntries[i].text())
            if inp < 1:
                inp = 1
            elif inp > 5:
                inp = 5
            argu.append(inp)
            outD[i] = inp
            self.dEntries[i].setText('%.2g' % inp)
            if not self.posTicks[i].isChecked():
                guess.append(float(self.posEntries[i].text()))
                struc.append(True)
            else:
                outPos[i] = float(self.posEntries[i].text())
                argu.append(outPos[i])
                struc.append(False)
            if not self.sigmaTicks[i].isChecked():
                inp = float(self.sigmaEntries[i].text())
                maxCq = max(maxCq,inp*1e6)
                guess.append(inp*1e6)
                struc.append(True)
            else:
                inp = float(self.sigmaEntries[i].text())
                maxCq = max(maxCq,inp*1e6)
                argu.append(inp*1e6)
                outSigma[i] = inp
                struc.append(False)
            if not self.ampTicks[i].isChecked():
                guess.append(float(self.ampEntries[i].text()))
                struc.append(True)
            else:
                outAmp[i] = float(self.ampEntries[i].text())
                argu.append(outAmp[i])
                struc.append(False)
            if not self.lorTicks[i].isChecked():
                guess.append(abs(float(self.lorEntries[i].text())))
                struc.append(True)
            else:
                outWidth[i] = abs(float(self.lorEntries[i].text()))
                argu.append(outWidth[i])
                struc.append(False)
            if not self.gaussTicks[i].text():
                guess.append(abs(float(self.gaussEntries[i].text())))
                struc.append(True)
            else:
                outGauss[i] = abs(float(self.gaussEntries[i].text()))
                argu.append(outGauss[i])
                struc.append(False)
        self.args = (numExp,struc,argu)
        self.setAngleStuff()
        numCq = int(self.cqGridEntry.text())
        numEta = int(self.etaGridEntry.text())
        self.genLib(len(self.parent.xax), I, maxCq*4.0, numCq, numEta)
        fitVal = scipy.optimize.fmin(self.fitFunc,guess, args=(self.parent.xax,np.real(self.parent.data1D)))
        counter = 0
        if struc[0]:
            self.bgrndEntry.setText('%.2g' % fitVal[counter])
            outBgrnd = fitVal[counter]
            counter +=1
        if struc[1]:
            self.slopeEntry.setText('%.2g' % fitVal[counter])
            outSlope = fitVal[counter]
            counter +=1
        for i in range(numExp):
            if struc[5*i+2]:
                self.posEntries[i].setText('%.2g' % fitVal[counter])
                outPos[i] = fitVal[counter]
                counter += 1
            if struc[5*i+3]:
                self.sigmaEntries[i].setText('%.2g' % (fitVal[counter]*1e-6))
                outCq[i] = fitVal[counter]
                counter += 1
            if struc[5*i+4]:
                self.ampEntries[i].setText('%.2g' % fitVal[counter])
                outAmp[i] = fitVal[counter]
                counter += 1
            if struc[5*i+5]:
                self.lorEntries[i].setText('%.2g' % abs(fitVal[counter]))
                outWidth[i] = abs(fitVal[counter])
                counter += 1
            if struc[5*i+6]:
                self.gaussEntries[i].setText('%.2g' % abs(fitVal[counter]))
                outGauss[i] = abs(fitVal[counter])
                counter += 1
        self.disp(outBgrnd,outSlope,outPos,outSigma,outD,outAmp,outWidth,outGauss)

    def sim(self):
        self.setCheng()
        if not self.setGrid():
            print("One of the inputs is not valid")
            return
        if not self.checkInputs():
            print("One of the inputs is not valid")
            return
        numExp = self.numExp.currentIndex()+1
        bgrnd = float(self.bgrndEntry.text())
        slope = float(self.slopeEntry.text())
        pos = np.zeros(numExp)
        sigma = np.zeros(numExp)
        d = np.zeros(numExp)
        amp = np.zeros(numExp)
        width = np.zeros(numExp)
        gauss = np.zeros(numExp)
        I = self.checkI(self.IEntry.currentIndex())
        for i in range(numExp):
            pos[i] = safeEval(self.posEntries[i].text())
            sigma[i] = safeEval(self.sigmaEntries[i].text())*1e6
            d[i] = safeEval(self.dEntries[i].text())
            amp[i] = safeEval(self.ampEntries[i].text())
            width[i] = safeEval(self.lorEntries[i].text())
            gauss[i] = safeEval(self.gaussEntries[i].text())
        self.setAngleStuff()
        numCq = int(self.cqGridEntry.text())
        numEta = int(self.etaGridEntry.text())
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
class MainPlotWindow(QtGui.QWidget):
    def __init__(self,father,oldMainWindow):
        QtGui.QWidget.__init__(self,father)
        self.father = father 
        self.oldMainWindow = oldMainWindow
        self.fig = oldMainWindow.current.fig
        self.canvas = FigureCanvas(self.fig)
        #self.canvas = oldMainWindow.current.canvas
        self.ax = oldMainWindow.current.ax
        grid = QtGui.QGridLayout(self)
        grid.addWidget(self.canvas,0,0)
        self.frame1 = QtGui.QGridLayout()
        grid.addLayout(self.frame1,0,1)
        self.frame1.addWidget(QLabel("Title:"),0,0)
        self.titleEntry = QtGui.QLineEdit()
        self.titleEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.titleBackup = oldMainWindow.name
        self.titleEntry.setText(self.titleBackup)
        self.titleEntry.returnPressed.connect(self.updatePlot)
        self.frame1.addWidget(self.titleEntry,1,0)
        self.frame1.addWidget(QLabel("x-label:"),2,0)
        self.xlabelEntry = QtGui.QLineEdit()
        self.xlabelEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.xlabelBackup = self.ax.get_xlabel()
        self.xlabelEntry.setText(self.xlabelBackup)
        self.xlabelEntry.returnPressed.connect(self.updatePlot)
        self.frame1.addWidget(self.xlabelEntry,3,0)
        self.frame1.addWidget(QLabel("y-label:"),4,0)
        self.ylabelEntry = QtGui.QLineEdit()
        self.ylabelEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ylabelBackup = self.ax.get_ylabel()
        self.ylabelEntry.setText(self.ylabelBackup)
        self.ylabelEntry.returnPressed.connect(self.updatePlot)
        self.frame1.addWidget(self.ylabelEntry,5,0)
        self.xlimBackup = self.ax.get_xlim()
        self.frame1.addWidget(QLabel("x-limit left:"),6,0)
        self.xlimLeftEntry = QtGui.QLineEdit()
        self.xlimLeftEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.xlimLeftEntry.setText(str(self.xlimBackup[0]))
        self.xlimLeftEntry.returnPressed.connect(self.updatePlot)
        self.frame1.addWidget(self.xlimLeftEntry,7,0)
        self.frame1.addWidget(QLabel("x-limit right:"),8,0)
        self.xlimRightEntry = QtGui.QLineEdit()
        self.xlimRightEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.xlimRightEntry.setText(str(self.xlimBackup[1]))
        self.xlimRightEntry.returnPressed.connect(self.updatePlot)
        self.frame1.addWidget(self.xlimRightEntry,9,0)
        self.ylimBackup = self.ax.get_ylim()
        self.frame1.addWidget(QLabel("y-limit left:"),10,0)
        self.ylimLeftEntry = QtGui.QLineEdit()
        self.ylimLeftEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ylimLeftEntry.setText(str(self.ylimBackup[0]))
        self.ylimLeftEntry.returnPressed.connect(self.updatePlot)
        self.frame1.addWidget(self.ylimLeftEntry,11,0)
        self.frame1.addWidget(QLabel("y-limit right:"),12,0)
        self.ylimRightEntry = QtGui.QLineEdit()
        self.ylimRightEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ylimRightEntry.setText(str(self.ylimBackup[1]))
        self.ylimRightEntry.returnPressed.connect(self.updatePlot)
        self.frame1.addWidget(self.ylimRightEntry,13,0)

        self.widthBackup, self.heightBackup = self.fig.get_size_inches()
        self.widthBackup = self.widthBackup*2.54
        self.heightBackup = self.heightBackup*2.54
        self.frame1.addWidget(QLabel("Width [cm]:"),26,0)
        self.widthEntry = QtGui.QLineEdit()
        self.widthEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.widthEntry.setText(str(self.widthBackup))
        self.widthEntry.returnPressed.connect(self.updatePlot)
        self.frame1.addWidget(self.widthEntry,27,0)
        self.frame1.addWidget(QLabel("Height [cm]:"),28,0)
        self.heightEntry = QtGui.QLineEdit()
        self.heightEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.heightEntry.setText(str(self.heightBackup))
        self.heightEntry.returnPressed.connect(self.updatePlot)
        self.frame1.addWidget(self.heightEntry,29,0)
        self.frame1.addWidget(QLabel("File type:"),30,0)
        self.filetypeEntry = QtGui.QComboBox()
        self.fileOptions = ['svg', 'png', 'eps', 'jpg', 'pdf']
        self.filetypeEntry.addItems(self.fileOptions)
        self.frame1.addWidget(self.filetypeEntry,31,0)
        self.inFrame = QtGui.QGridLayout()
        self.frame1.addLayout(self.inFrame,32,0)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.cancel)
        self.inFrame.addWidget(cancelButton,0,0)
        okButton = QtGui.QPushButton("&Save")
        okButton.clicked.connect(self.save)
        self.inFrame.addWidget(okButton,0,1)
        grid.setColumnStretch(0,1)
        grid.setRowStretch(0,1)
        
    def rename(self,name):
        self.oldMainWindow.rename(name)
        
    def updatePlot(self, *args):
        self.fig.suptitle(self.titleEntry.text())
        self.ax.set_xlabel(self.xlabelEntry.text())
        self.ax.set_ylabel(self.ylabelEntry.text())
        self.ax.set_xlim((safeEval(self.xlimLeftEntry.text()),safeEval(self.xlimRightEntry.text())))
        self.ax.set_ylim((safeEval(self.ylimLeftEntry.text()),safeEval(self.ylimRightEntry.text())))
        #self.fig.set_size_inches((int(safeEval(self.widthEntry.text()))/2.54,int(safeEval(self.heightEntry.text()))/2.54))
        self.fig.canvas.draw()

    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()
        
    def addToView(self):
        self.show()

    def removeFromView(self):
        self.hide()

    def kill(self):
        self.oldMainWindow.kill()
        self.destroy()
        
    def save(self):
        self.updatePlot()
        self.fig.set_size_inches((int(safeEval(self.widthEntry.text()))/2.54,int(safeEval(self.heightEntry.text()))/2.54))
        f = QtGui.QFileDialog.getSaveFileName(self, 'Save File')
        if f:
            f=os.path.splitext(f)[0]+'.'+self.fileOptions[self.filetypeEntry.currentIndex()]
            self.fig.savefig(f)
        self.cancel()

    def cancel(self):
        self.fig.suptitle(self.titleBackup)
        self.ax.set_xlabel(self.xlabelBackup)
        self.ax.set_ylabel(self.ylabelBackup)
        self.ax.set_xlim((self.xlimBackup[0],self.xlimBackup[1]))
        self.ax.set_ylim((self.ylimBackup[0],self.ylimBackup[1]))
        self.fig.set_size_inches((self.widthBackup/2.54,self.heightBackup/2.54))
        self.father.closeSaveFigure(self.oldMainWindow)

######################################################################

class FitAllSelectionWindow(QtGui.QWidget): 
    def __init__(self, parent, fitNames):
        QtGui.QWidget.__init__(self,parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Select output")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid,0,0,1,2)
        
        self.ticks = []
        for i in range(len(fitNames)):
            self.ticks.append(QtGui.QCheckBox(fitNames[i]))
            grid.addWidget(self.ticks[i],i,0)
        
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton,1,0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.fit)
        layout.addWidget(okButton,1,1)
        self.show()
        self.setFixedSize(self.size())
        
    def fit(self):
        returnVals = []
        for i in self.ticks:
            returnVals.append(i.isChecked())
        self.deleteLater()
        self.father.fitAllFunc(np.array(returnVals,dtype=bool))
