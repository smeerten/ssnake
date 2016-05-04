#!/usr/bin/env python

# Copyright 2016 Bas van Meerten and Wouter Franssen

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

from PyQt4 import QtGui, QtCore
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from widgetClasses import *
import os
import copy

#####################################################################################
class MainPlotWindow(QtGui.QWidget):
    def __init__(self, father, oldMainWindow):
        QtGui.QWidget.__init__(self, father)
        self.father = father 
        self.oldMainWindow = oldMainWindow
        self.fig = oldMainWindow.current.fig
        self.canvas = FigureCanvas(self.fig)
        #self.canvas = oldMainWindow.current.canvas
        self.ax = oldMainWindow.current.ax
        grid = QtGui.QGridLayout(self)
        scroll2 = QtGui.QScrollArea()
        grid.addWidget(scroll2, 0, 0)
        scroll2.setWidget(self.canvas)
        #grid.addWidget(self.canvas, 0, 0)
        self.frame1 = QtGui.QGridLayout()
        grid.addLayout(self.frame1, 0, 1)
        scroll = QtGui.QScrollArea()
        self.frame1.addWidget(scroll, 0, 0)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        content = QtGui.QWidget()
        self.optionFrame = QtGui.QGridLayout(content)
        self.optionFrame.addWidget(QLabel("Title:"), 0, 0)
        self.titleEntry = QtGui.QLineEdit()
        self.titleEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.titleBackup = oldMainWindow.masterData.name
        self.titleEntry.setText(self.titleBackup)
        self.titleEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.titleEntry, 1, 0)
        self.optionFrame.addWidget(QLabel("x-label:"), 2, 0)
        self.xlabelEntry = QtGui.QLineEdit()
        self.xlabelEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.xlabelBackup = self.ax.get_xlabel()
        self.xlabelEntry.setText(self.xlabelBackup)
        self.xlabelEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.xlabelEntry, 3, 0)
        self.optionFrame.addWidget(QLabel("y-label:"), 4, 0)
        self.ylabelEntry = QtGui.QLineEdit()
        self.ylabelEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ylabelBackup = self.ax.get_ylabel()
        self.ylabelEntry.setText(self.ylabelBackup)
        self.ylabelEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.ylabelEntry, 5, 0)
        self.xlimBackup = self.ax.get_xlim()
        self.optionFrame.addWidget(QLabel("x-limit left:"), 6, 0)
        self.xlimLeftEntry = QtGui.QLineEdit()
        self.xlimLeftEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.xlimLeftEntry.setText(str(self.xlimBackup[0]))
        self.xlimLeftEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.xlimLeftEntry, 7, 0)
        self.optionFrame.addWidget(QLabel("x-limit right:"), 8, 0)
        self.xlimRightEntry = QtGui.QLineEdit()
        self.xlimRightEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.xlimRightEntry.setText(str(self.xlimBackup[1]))
        self.xlimRightEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.xlimRightEntry, 9, 0)
        self.ylimBackup = self.ax.get_ylim()
        self.optionFrame.addWidget(QLabel("y-limit down:"), 10, 0)
        self.ylimLeftEntry = QtGui.QLineEdit()
        self.ylimLeftEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ylimLeftEntry.setText(str(self.ylimBackup[0]))
        self.ylimLeftEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.ylimLeftEntry, 11, 0)
        self.optionFrame.addWidget(QLabel("y-limit up:"), 12, 0)
        self.ylimRightEntry = QtGui.QLineEdit()
        self.ylimRightEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ylimRightEntry.setText(str(self.ylimBackup[1]))
        self.ylimRightEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.ylimRightEntry, 13, 0)
        self.widthBackup, self.heightBackup = self.fig.get_size_inches()
        self.widthBackup = self.widthBackup*2.54
        self.heightBackup = self.heightBackup*2.54
        self.optionFrame.addWidget(QLabel("Width [cm]:"), 26, 0)
        self.widthEntry = QtGui.QLineEdit()
        self.widthEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.widthEntry.setText(str(self.widthBackup))
        self.widthEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.widthEntry, 27, 0)
        self.optionFrame.addWidget(QLabel("Height [cm]:"), 28, 0)
        self.heightEntry = QtGui.QLineEdit()
        self.heightEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.heightEntry.setText(str(self.heightBackup))
        self.heightEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.heightEntry, 29, 0)
        self.titleFontSizeBackup = 12
        self.optionFrame.addWidget(QLabel("Title font size:"), 30, 0)
        self.titleFontSizeEntry = QtGui.QLineEdit()
        self.titleFontSizeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.titleFontSizeEntry.setText(str(self.titleFontSizeBackup))
        self.titleFontSizeEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.titleFontSizeEntry, 31, 0)
        self.xlabelFontSizeBackup = 12
        self.optionFrame.addWidget(QLabel("X-label font size:"), 32, 0)
        self.xlabelFontSizeEntry = QtGui.QLineEdit()
        self.xlabelFontSizeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.xlabelFontSizeEntry.setText(str(self.xlabelFontSizeBackup))
        self.xlabelFontSizeEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.xlabelFontSizeEntry, 33, 0)
        self.ylabelFontSizeBackup = 12
        self.optionFrame.addWidget(QLabel("Y-label font size:"), 34, 0)
        self.ylabelFontSizeEntry = QtGui.QLineEdit()
        self.ylabelFontSizeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ylabelFontSizeEntry.setText(str(self.ylabelFontSizeBackup))
        self.ylabelFontSizeEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.ylabelFontSizeEntry, 35, 0)
        self.xtickFontSizeBackup = 12
        self.optionFrame.addWidget(QLabel("X-ticks font size:"), 36, 0)
        self.xtickFontSizeEntry = QtGui.QLineEdit()
        self.xtickFontSizeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.xtickFontSizeEntry.setText(str(self.xtickFontSizeBackup))
        self.xtickFontSizeEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.xtickFontSizeEntry, 37, 0)
        self.ytickFontSizeBackup = 12
        self.optionFrame.addWidget(QLabel("Y-ticks font size:"), 38, 0)
        self.ytickFontSizeEntry = QtGui.QLineEdit()
        self.ytickFontSizeEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.ytickFontSizeEntry.setText(str(self.ytickFontSizeBackup))
        self.ytickFontSizeEntry.returnPressed.connect(self.updatePlot)
        self.optionFrame.addWidget(self.ytickFontSizeEntry, 39, 0)
        self.legend = self.ax.legend()
        self.legend.draggable(True)
        self.legendPos = 'best'
        self.legendTextList = []
        for line in self.legend.get_texts():
            self.legendTextList.append(line.get_text())
        self.legend.set_visible(False)
        self.legendCheck = QtGui.QCheckBox('Legend')
        self.legendCheck.stateChanged.connect(self.updateLegend)
        self.optionFrame.addWidget(self.legendCheck, 40, 0)
        legendButton = QtGui.QPushButton('Legend settings')
        legendButton.clicked.connect(lambda : LegendWindow(self))
        self.optionFrame.addWidget(legendButton, 41, 0)
        
        self.inFrame = QtGui.QGridLayout()
        self.frame1.addLayout(self.inFrame, 1, 0)
        self.inFrame.addWidget(QLabel("File type:"), 0, 0, 1, 2)
        self.filetypeEntry = QtGui.QComboBox()
        self.fileOptions = ['svg', 'png', 'eps', 'jpg', 'pdf']
        self.filetypeEntry.addItems(self.fileOptions)
        self.inFrame.addWidget(self.filetypeEntry, 1, 0, 1, 2)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.cancel)
        self.inFrame.addWidget(cancelButton, 2, 0)
        okButton = QtGui.QPushButton("&Save")
        okButton.clicked.connect(self.save)
        self.inFrame.addWidget(okButton, 2, 1)
        grid.setColumnStretch(0, 1)
        grid.setRowStretch(0, 1)
        self.optionFrame.setAlignment(QtCore.Qt.AlignTop)
        self.grid = grid
        scroll.setWidget(content)
        self.updatePlot()
        
    def rename(self, name):
        self.oldMainWindow.rename(name)

    def updateLegend(self, *args):
        if self.legendCheck.isChecked():
            self.legend = self.ax.legend(self.legendTextList, loc=self.legendPos)
            self.legend.draggable(True)
        else:
            if self.legend is not None:
                self.legend.set_visible(False)
        self.updatePlot()
        
    def updatePlot(self, *args):
        self.fig.suptitle(self.titleEntry.text(), fontsize=safeEval(self.titleFontSizeEntry.text()))
        self.ax.set_xlabel(self.xlabelEntry.text(), fontsize=safeEval(self.xlabelFontSizeEntry.text()))
        self.ax.set_ylabel(self.ylabelEntry.text(), fontsize=safeEval(self.ylabelFontSizeEntry.text()))
        self.ax.set_xlim((safeEval(self.xlimLeftEntry.text()), safeEval(self.xlimRightEntry.text())))
        self.ax.set_ylim((safeEval(self.ylimLeftEntry.text()), safeEval(self.ylimRightEntry.text())))
        self.ax.tick_params(axis='x', labelsize=safeEval(self.xtickFontSizeEntry.text()))
        self.ax.tick_params(axis='y', labelsize=safeEval(self.ytickFontSizeEntry.text()))
        self.fig.set_size_inches((int(safeEval(self.widthEntry.text()))/2.54, int(safeEval(self.heightEntry.text()))/2.54))
        self.canvas.draw()
        self.canvas.adjustSize()

    def get_mainWindow(self):
        return self.oldMainWindow
        
    def get_masterData(self):
        return self.oldMainWindow.get_masterData()
    
    def get_current(self):
        return self.oldMainWindow.get_current()

    def kill(self):
        for i in reversed(range(self.grid.count())): 
            self.grid.itemAt(i).widget().deleteLater()
        self.grid.deleteLater()
        self.oldMainWindow.kill()
        del self.fig
        del self.canvas
        self.deleteLater()
        
    def save(self):
        self.updatePlot()
        self.fig.set_size_inches((int(safeEval(self.widthEntry.text()))/2.54, int(safeEval(self.heightEntry.text()))/2.54))
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum] #Set name of file to be saved to workspace name to start
        f = QtGui.QFileDialog.getSaveFileName(self, 'Save File', self.father.LastLocation+os.path.sep+WorkspaceName+'.'+self.fileOptions[self.filetypeEntry.currentIndex()])
        if f:
            self.father.LastLocation = os.path.dirname(f)
            self.fig.savefig(f)
            self.cancel()

    def cancel(self):
        self.fig.suptitle(self.titleBackup, fontsize=self.titleFontSizeBackup)
        self.ax.set_xlabel(self.xlabelBackup, fontsize=self.xlabelFontSizeBackup)
        self.ax.set_ylabel(self.ylabelBackup, fontsize=self.ylabelFontSizeBackup)
        self.ax.set_xlim((self.xlimBackup[0], self.xlimBackup[1]))
        self.ax.set_ylim((self.ylimBackup[0], self.ylimBackup[1]))
        self.ax.tick_params(axis='x', labelsize=self.xtickFontSizeBackup)
        self.ax.tick_params(axis='y', labelsize=self.ytickFontSizeBackup)
        if self.legend is not None:
            self.legend.set_visible(False)
        self.fig.set_size_inches((self.widthBackup/2.54, self.heightBackup/2.54))
        self.grid.deleteLater()
        del self.canvas
        del self.fig
        self.father.closeSaveFigure(self.oldMainWindow)
        self.deleteLater()

#####################################################################################################################
class LegendWindow(QtGui.QWidget):
    def __init__(self, parent):
        QtGui.QWidget.__init__(self, parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Legend")
        layout = QtGui.QGridLayout(self)
        grid = QtGui.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(QLabel("Legend position:"), 0, 0)
        self.posVal = self.father.legendPos
        self.posEntry = QtGui.QLineEdit()
        self.posEntry.setAlignment(QtCore.Qt.AlignHCenter)
        self.posEntry.setText(self.posVal)
        self.posEntry.returnPressed.connect(self.preview)
        grid.addWidget(self.posEntry, 1, 0)
        grid.addWidget(QLabel("Legend:"), 2, 0)
        self.father.legendCheck.setChecked(True)
        self.father.legendTextList = []
        self.spinBox = QtGui.QSpinBox()
        self.spinBox.setMaximum(len(self.father.legendTextList)-1)
        self.spinBox.valueChanged.connect(self.changeEdit)
        grid.addWidget(self.spinBox, 3, 0)
        self.legendEditList = []
        for i in range(len(self.father.legendTextList)):
            self.legendEditList.append(QtGui.QLineEdit())
            self.legendEditList[i].setAlignment(QtCore.Qt.AlignHCenter)
            self.legendEditList[i].setText(self.father.legendTextList[i])
            self.legendEditList[i].returnPressed.connect(self.preview)
            grid.addWidget(self.legendEditList[i], 4, 0)
            self.legendEditList[i].setVisible(False)
        self.legendEditList[0].setVisible(True)
        cancelButton = QtGui.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtGui.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 1, 1)
        self.show()
        self.setFixedSize(self.size())
        self.setGeometry(self.frameSize().width()-self.geometry().width(), self.frameSize().height()-self.geometry().height(), 0, 0)

    def changeEdit(self, num):
        for i in range(len(self.legendEditList)):
            self.legendEditList[i].setVisible(False)
        self.legendEditList[num].setVisible(True)
        
    def preview(self, *args):
        tmp = copy.deepcopy(self.father.legendTextList)
        for i in range(len(self.legendEditList)):
            tmp[i] = self.legendEditList[i].text()
        env = vars(np).copy()
        try:
            inp = eval(self.posEntry.text(), env)                
        except:
            inp = self.posEntry.text()
        self.father.ax.legend(tmp, loc=inp)
        self.father.legend.draggable(True)
        self.father.canvas.draw()
        
    def closeEvent(self, *args):
        self.deleteLater()
        self.father.updateLegend()
        
    def applyAndClose(self):
        for i in range(len(self.legendEditList)):
            self.father.legendTextList[i] = self.legendEditList[i].text()
        env = vars(np).copy()
        try:
            inp = eval(self.posEntry.text(), env)                
        except:
            inp = self.posEntry.text()
        self.father.legendPos = inp
        self.closeEvent()
        
 
