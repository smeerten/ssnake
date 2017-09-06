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

import numpy as np
import matplotlib.gridspec as gridspec
try:
    from PyQt4 import QtGui, QtCore
    from PyQt4 import QtGui as QtWidgets    
except ImportError:
    from PyQt5 import QtGui, QtCore, QtWidgets
import spectrum_classes
import warnings
#########################################################################################################
# the class from which the 1d data is displayed, the operations which only edit the content of this class are for previewing


class Plot1DFrame(object):

    def __init__(self, root, fig, canvas):
        self.root = root
        self.fig = fig
        self.canvas = canvas
        self.fig.clf()
        if isinstance(self, spectrum_classes.CurrentContour):
            self.gs = gridspec.GridSpec(2, 2, width_ratios=[self.root.father.defaultWidthRatio, 1], height_ratios=[1, self.root.father.defaultHeightRatio])
            self.ax = self.fig.add_subplot(self.gs[2])
#            self.x_ax = self.fig.add_subplot(gs[0], sharex=self.ax,axisbg='#bfbfbf')
#            self.y_ax = self.fig.add_subplot(gs[3], sharey=self.ax,axisbg='#bfbfbf')
            self.x_ax = self.fig.add_subplot(self.gs[0], sharex=self.ax,axisbg='none',frameon=False)
            self.y_ax = self.fig.add_subplot(self.gs[3], sharey=self.ax,axisbg='none',frameon=False)
            self.fig.subplots_adjust(hspace=0)
            self.fig.subplots_adjust(wspace=0)
            self.x_ax.axes.get_xaxis().set_visible(False)
            self.x_ax.axes.get_yaxis().set_visible(False)
            self.y_ax.axes.get_xaxis().set_visible(False)
            self.y_ax.axes.get_yaxis().set_visible(False)            
        #elif isinstance(self, spectrum_classes.CurrentSkewed):
        #    self.ax = self.fig.add_subplot(1, 1, 1, projection='3d')
        #    self.ax.disable_mouse_rotation()
        else:
            self.ax = self.fig.add_subplot(111)
        self.leftMouse = False  # is the left mouse button currently pressed
        self.panX = None  # start position of dragging the spectrum
        self.panY = None  # start position of dragging the spectrum
        self.zoomX1 = None  # first corner of the zoombox
        self.zoomY1 = None  # first corner of the zoombox
        self.zoomX2 = None  # second corner of the zoombox
        self.zoomY2 = None  # second corner of the zoombox
        self.rect = [None, None, None, None]  # lines for zooming or peak picking
        self.rightMouse = False  # is the right mouse button currently pressed
        self.peakPick = False  # currently peakPicking (if 2 display cross)
        self.peakPickFunc = None  # the function that needs to be called after peakPicking
        # variables to be initialized
        self.spec = 0
        self.spec2 = 0

    def kill(self):
        pass

    def plotReset(self):  # this function needs to be overriden by the classes who inherit from Plot1DFrame
        pass

    ################
    # mouse events #
    ################

    def peakPickReset(self):
        if self.rect[0] is not None:
            try:
                self.rect[0].remove()
            except:
                pass
            self.canvas.draw_idle()
        if self.rect[1] is not None:
            try:
                self.rect[1].remove()
            except:
                pass
            self.canvas.draw_idle()
        self.rect = [None, None, None, None]
        self.peakPick = False
        self.peakPickFunc = None

    def scroll(self, event):
        modifiers = QtWidgets.QApplication.keyboardModifiers()
        if modifiers == QtCore.Qt.ShiftModifier:
            self.altScroll(event)
        else:
            if self.rightMouse:
                middle = (self.xmaxlim + self.xminlim) / 2.0
                width = self.xmaxlim - self.xminlim
                if modifiers == QtCore.Qt.ControlModifier:
                    width = width * 0.6**event.step
                else:
                    width = width * 0.9**event.step
                self.xmaxlim = middle + width / 2.0
                self.xminlim = middle - width / 2.0
                if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                    self.ax.set_xlim(self.xmaxlim, self.xminlim)
                else:
                    self.ax.set_xlim(self.xminlim, self.xmaxlim)
            else:
                if getattr(self, '__module__', None) == spectrum_classes.__name__:
                    noZeroScroll = isinstance(self, (spectrum_classes.CurrentContour, spectrum_classes.CurrentStacked)) or not self.root.father.defaultZeroScroll
                else:
                    noZeroScroll = not self.rootwindow.mainProgram.defaultZeroScroll
                if noZeroScroll:
                    middle = (self.ymaxlim + self.yminlim) / 2.0
                    width = self.ymaxlim - self.yminlim
                    if modifiers == QtCore.Qt.ControlModifier:
                        width = width * 0.6**event.step
                    else:
                        width = width * 0.9**event.step
                    self.ymaxlim = middle + width / 2.0
                    self.yminlim = middle - width / 2.0
                else:
                    if modifiers == QtCore.Qt.ControlModifier:
                        self.ymaxlim *= 0.6**event.step
                        self.yminlim *= 0.6**event.step
                    else:
                        self.ymaxlim *= 0.9**event.step
                        self.yminlim *= 0.9**event.step
                if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
                    self.ax.set_ylim(self.ymaxlim, self.yminlim)
                else:
                    self.ax.set_ylim(self.yminlim, self.ymaxlim)
            self.canvas.update()
            self.canvas.draw_idle()

    def altScroll(self, event):
        pass

    def altReset(self):
        pass

    def buttonPress(self, event):
        if event.button == 1 and not self.peakPick:
            self.leftMouse = True
            #if isinstance(self, spectrum_classes.CurrentSkewed):
            #    a = self.ax.format_coord(event.xdata, event.ydata)
            #    self.zoomX1 = float(a[2:14])
            #    self.zoomY1 = float(a[18:30])
            #else:
            self.zoomX1 = event.xdata
            self.zoomY1 = event.ydata
        elif (event.button == 3) and event.dblclick:
            modifiers = QtWidgets.QApplication.keyboardModifiers()
            if modifiers == QtCore.Qt.ShiftModifier:
                self.altReset()
            else:
                self.plotReset()
        elif event.button == 3:
            self.rightMouse = True
            #if isinstance(self, spectrum_classes.CurrentSkewed):
            #    a = self.ax.format_coord(event.xdata, event.ydata)
            #    self.panX = float(a[2:14])
            #    self.panY = float(a[18:30])
            #else:
            self.panX = event.xdata
            self.panY = event.ydata

    def buttonRelease(self, event):
        if event.button == 1:
            if self.peakPick:
                if self.rect[0] is not None:
                    try:
                        self.rect[0].remove()
                    finally:
                        self.rect[0] = None
                    if self.rect[1] is not None:
                        try:
                            self.rect[1].remove()
                        finally:
                            self.rect[1] = None
                    self.peakPick = False
                    idx = np.argmin(np.abs(self.line_xdata - event.xdata))
                    if self.peakPickFunc is not None:
                        self.peakPickFunc((idx, self.line_xdata[idx], self.line_ydata[idx]))
                    if not self.peakPick:  # check if peakpicking is still required
                        self.peakPickFunc = None
            else:
                self.leftMouse = False
                try:
                    if self.rect[0] is not None:
                        self.rect[0].remove()
                    if self.rect[1] is not None:
                        self.rect[1].remove()
                    if self.rect[2] is not None:
                        self.rect[2].remove()
                    if self.rect[3] is not None:
                        self.rect[3].remove()
                finally:
                    self.rect = [None, None, None, None]
                if self.zoomX2 is not None and self.zoomY2 is not None:
                    self.xminlim = min([self.zoomX1, self.zoomX2])
                    self.xmaxlim = max([self.zoomX1, self.zoomX2])
                    self.yminlim = min([self.zoomY1, self.zoomY2])
                    self.ymaxlim = max([self.zoomY1, self.zoomY2])
                    if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                        self.ax.set_xlim(self.xmaxlim, self.xminlim)
                    else:
                        self.ax.set_xlim(self.xminlim, self.xmaxlim)
                    if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
                        self.ax.set_ylim(self.ymaxlim, self.yminlim)
                    else:
                        self.ax.set_ylim(self.yminlim, self.ymaxlim)
                self.zoomX1 = None
                self.zoomX2 = None
                self.zoomY1 = None
                self.zoomY2 = None
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw_idle()

    def pan(self, event):
        modifiers = QtWidgets.QApplication.keyboardModifiers()
        if self.rightMouse and self.panX is not None and self.panY is not None:
            #if isinstance(self, spectrum_classes.CurrentSkewed):
            #    a = self.ax.format_coord(event.xdata, event.ydata)
            #    diffx = float(a[2:14]) - self.panX
            #    diffy = float(a[18:30]) - self.panY
            #else:
            inv = self.ax.transData.inverted()
            point = inv.transform((event.x, event.y))
            diffx = point[0] - self.panX
            diffy = point[1] - self.panY
            
            if modifiers == QtCore.Qt.ControlModifier:
                self.xmaxlim = self.xmaxlim - diffx
                self.xminlim = self.xminlim - diffx
            elif modifiers == QtCore.Qt.ShiftModifier:
                self.ymaxlim = self.ymaxlim - diffy
                self.yminlim = self.yminlim - diffy
            else:
                self.xmaxlim = self.xmaxlim - diffx
                self.xminlim = self.xminlim - diffx
                self.ymaxlim = self.ymaxlim - diffy
                self.yminlim = self.yminlim - diffy
            if self.spec > 0 and not isinstance(self, spectrum_classes.CurrentArrayed):
                self.ax.set_xlim(self.xmaxlim, self.xminlim)
            else:
                self.ax.set_xlim(self.xminlim, self.xmaxlim)
            if self.spec2 > 0 and isinstance(self, spectrum_classes.CurrentContour):
                self.ax.set_ylim(self.ymaxlim, self.yminlim)
            else:
                self.ax.set_ylim(self.yminlim, self.ymaxlim)
            self.canvas.draw_idle()
        elif self.peakPick:
            if self.rect[0] is not None:
                try:
                    self.rect[0].remove()
                except:
                    pass
                self.rect[0] = None
            if self.rect[1] is not None:
                try:
                    self.rect[1].remove()
                except:
                    pass
                self.rect[1] = None
            if event.xdata is not None:
                self.rect[0] = self.ax.axvline(event.xdata, c='k', linestyle='--')
            if self.peakPick == 2:
                if event.ydata is not None:
                    self.rect[1] = self.ax.axhline(event.ydata, c='k', linestyle='--')
            self.canvas.draw_idle()
        elif self.leftMouse and (self.zoomX1 is not None) and (self.zoomY1 is not None):
            #if isinstance(self, spectrum_classes.CurrentSkewed):
            #    a = self.ax.format_coord(event.xdata, event.ydata)
            #    self.zoomX2 = float(a[2:14])
            #    self.zoomY2 = float(a[18:30])
            #else:
            inv = self.ax.transData.inverted()
            point = inv.transform((event.x, event.y))
            self.zoomX2 = point[0]
            self.zoomY2 = point[1]
            if self.rect[0] is not None:
                try:
                    if self.rect[0] is not None:
                        self.rect[0].remove()
                    if self.rect[1] is not None:
                        self.rect[1].remove()
                    if self.rect[2] is not None:
                        self.rect[2].remove()
                    if self.rect[3] is not None:
                        self.rect[3].remove()
                finally:
                    self.rect = [None, None, None, None]
            self.rect[0], = self.ax.plot([self.zoomX1, self.zoomX2], [self.zoomY2, self.zoomY2], 'k', clip_on=False)
            self.rect[1], = self.ax.plot([self.zoomX1, self.zoomX2], [self.zoomY1, self.zoomY1], 'k', clip_on=False)
            self.rect[2], = self.ax.plot([self.zoomX1, self.zoomX1], [self.zoomY1, self.zoomY2], 'k', clip_on=False)
            self.rect[3], = self.ax.plot([self.zoomX2, self.zoomX2], [self.zoomY1, self.zoomY2], 'k', clip_on=False)
            self.canvas.draw_idle()
