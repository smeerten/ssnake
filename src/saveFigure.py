#!/usr/bin/env python3

# Copyright 2016 - 2024 Bas van Meerten and Wouter Franssen

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
import os
import copy
import matplotlib
matplotlib.rc('svg', fonttype='none')
from matplotlib.colors import colorConverter
import widgetClasses as wc
from safeEval import safeEval
from ssNake import QtGui, QtCore, QtWidgets, FigureCanvas

#####################################################################################


class SaveFigureWindow(QtWidgets.QWidget):
    """
    The window for saving a figure.
    This window replaces the original figure.
    """

    def __init__(self, father, oldMainWindow):
        """
        Initializes the SaveFigureWindow.

        Parameters
        ----------
        father : MainProgram
            The main program of ssNake.
        oldMainWindow : Main1DWindow
            The figure window that this window replaces.
        """
        super(SaveFigureWindow, self).__init__(father)
        self.father = father
        self.oldMainWindow = oldMainWindow
        self.rename = self.oldMainWindow.rename                  # Forward function
        self.get_masterData = self.oldMainWindow.get_masterData  # Forward function
        self.get_current = self.oldMainWindow.get_current        # Forward function
        self.fig = oldMainWindow.current.fig
        self.canvas = FigureCanvas(self.fig)
        self.canvas.mpl_connect('pick_event', self.pickHandler)
        self.ax = oldMainWindow.current.ax
        grid = QtWidgets.QGridLayout(self)
        scroll2 = QtWidgets.QScrollArea()
        grid.addWidget(scroll2, 0, 0)
        scroll2.setWidget(self.canvas)
        self.frame1 = QtWidgets.QGridLayout()
        grid.addLayout(self.frame1, 0, 1)
        scroll = QtWidgets.QScrollArea()
        self.frame1.addWidget(scroll, 0, 0)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        content = QtWidgets.QWidget()
        self.optionFrame = QtWidgets.QGridLayout(content)
        self.optionFrame.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        self.labelGroup = QtWidgets.QGroupBox('Labels:')
        self.labelFrame = QtWidgets.QGridLayout()
        self.labelFrame.addWidget(wc.QLeftLabel('Title:'), 0, 0)
        self.titleBackup = oldMainWindow.get_masterData().name
        self.titleEntry = wc.QLineEdit(self.titleBackup, self.updatePlot)
        self.labelFrame.addWidget(self.titleEntry, 0, 1)
        self.labelFrame.addWidget(wc.QLeftLabel("x:"), 1, 0)
        self.xlabelBackup = self.ax.get_xlabel()
        self.xlabelEntry = wc.QLineEdit(self.xlabelBackup, self.updatePlot)
        self.labelFrame.addWidget(self.xlabelEntry, 1, 1)
        self.labelFrame.addWidget(wc.QLeftLabel("y:"), 2, 0)
        self.ylabelBackup = self.ax.get_ylabel()
        self.ylabelEntry = wc.QLineEdit(self.ylabelBackup, self.updatePlot)
        self.labelFrame.addWidget(self.ylabelEntry, 2, 1)
        self.labelGroup.setLayout(self.labelFrame)
        self.optionFrame.addWidget(self.labelGroup, 0, 0)
        # limits
        self.limitsGroup = QtWidgets.QGroupBox('Limits:')
        self.limitsFrame = QtWidgets.QGridLayout()
        self.xlimBackup = self.ax.get_xlim()
        self.limitsFrame.addWidget(wc.QLeftLabel("x left:"), 0, 0)
        self.xlimLeftEntry = wc.QLineEdit(self.xlimBackup[0], self.updatePlot)
        self.limitsFrame.addWidget(self.xlimLeftEntry, 0, 1)
        self.limitsFrame.addWidget(wc.QLeftLabel("x right:"), 1, 0)
        self.xlimRightEntry = wc.QLineEdit(self.xlimBackup[1], self.updatePlot)
        self.limitsFrame.addWidget(self.xlimRightEntry, 1, 1)
        self.ylimBackup = self.ax.get_ylim()
        self.limitsFrame.addWidget(wc.QLeftLabel("y down:"), 2, 0)
        self.ylimLeftEntry = wc.QLineEdit(self.ylimBackup[0], self.updatePlot)
        self.limitsFrame.addWidget(self.ylimLeftEntry, 2, 1)
        self.limitsFrame.addWidget(wc.QLeftLabel("y up:"), 3, 0)
        self.ylimRightEntry = wc.QLineEdit(self.ylimBackup[1], self.updatePlot)
        self.limitsFrame.addWidget(self.ylimRightEntry, 3, 1)
        self.limitsGroup.setLayout(self.limitsFrame)
        self.optionFrame.addWidget(self.limitsGroup, 1, 0)
        # Dimensions
        self.dimensionsGroup = QtWidgets.QGroupBox('Dimensions:')
        self.dimensionsFrame = QtWidgets.QGridLayout()
        self.widthBackup, self.heightBackup = self.fig.get_size_inches()
        self.widthBackup = self.widthBackup * 2.54
        self.heightBackup = self.heightBackup * 2.54
        self.dimensionsFrame.addWidget(wc.QLeftLabel("Width [cm]:"), 0, 0)
        self.widthEntry = QtWidgets.QDoubleSpinBox()
        self.widthEntry.setSingleStep(0.1)
        self.widthEntry.setMinimum(0)
        self.widthEntry.setMaximum(1e6)
        self.widthEntry.setValue(self.widthBackup)
        self.widthEntry.valueChanged.connect(self.updatePlot)
        self.dimensionsFrame.addWidget(self.widthEntry, 0, 1)
        self.dimensionsFrame.addWidget(wc.QLeftLabel("Height [cm]:"), 1, 0)
        self.heightEntry = QtWidgets.QDoubleSpinBox()
        self.heightEntry.setSingleStep(0.1)
        self.heightEntry.setMinimum(0)
        self.heightEntry.setMaximum(1e6)
        self.heightEntry.setValue(self.heightBackup)
        self.heightEntry.valueChanged.connect(self.updatePlot)
        self.dimensionsFrame.addWidget(self.heightEntry, 1, 1)
        self.dimensionsFrame.addWidget(wc.QLeftLabel("dpi:"), 2, 0)
        self.dpiEntry = QtWidgets.QSpinBox()
        self.dpiEntry.setSingleStep(1)
        self.dpiEntry.setMinimum(0)
        self.dpiEntry.setMaximum(int(1e6))
        self.dpiEntry.setValue(int(self.fig.dpi))
        self.dpiEntry.valueChanged.connect(self.updatePlot)
        self.dimensionsFrame.addWidget(self.dpiEntry, 2, 1)
        self.dimensionsGroup.setLayout(self.dimensionsFrame)
        self.optionFrame.addWidget(self.dimensionsGroup, 2, 0)
        # Font size
        self.fontGroup = QtWidgets.QGroupBox('Font sizes:')
        self.fontFrame = QtWidgets.QGridLayout()
        self.mainFontLabel = wc.QLeftLabel("Main:")
        self.fontFrame.addWidget(self.mainFontLabel, 0, 0)
        self.mainFontSizeEntry = QtWidgets.QDoubleSpinBox()
        self.mainFontSizeEntry.setSingleStep(0.1)
        self.mainFontSizeEntry.setMinimum(0)
        self.mainFontSizeEntry.setValue(self.ax.xaxis.get_label().get_fontsize())
        self.mainFontSizeEntry.valueChanged.connect(self.updatePlot)
        self.fontFrame.addWidget(self.mainFontSizeEntry, 0, 1)
        self.fontDetailsCheck = QtWidgets.QCheckBox('Details')
        self.fontDetailsCheck.stateChanged.connect(self.fontCheckChanged)
        self.fontFrame.addWidget(self.fontDetailsCheck, 1, 0)
        self.titleFontSizeBackup = self.ax.title.get_fontsize()
        self.titleFontLabel = wc.QLeftLabel("Title:")
        self.titleFontLabel.hide()
        self.fontFrame.addWidget(self.titleFontLabel, 2, 0)
        self.titleFontSizeEntry = QtWidgets.QDoubleSpinBox()
        self.titleFontSizeEntry.setSingleStep(0.1)
        self.titleFontSizeEntry.setMinimum(0)
        self.titleFontSizeEntry.setValue(self.titleFontSizeBackup)
        self.titleFontSizeEntry.valueChanged.connect(self.updatePlot)
        self.titleFontSizeEntry.hide()
        self.fontFrame.addWidget(self.titleFontSizeEntry, 2, 1)
        self.xlabelFontSizeBackup = self.ax.xaxis.get_label().get_fontsize()
        self.xlabelFontLabel = wc.QLeftLabel("X-label:")
        self.xlabelFontLabel.hide()
        self.fontFrame.addWidget(self.xlabelFontLabel, 3, 0)
        self.xlabelFontSizeEntry = QtWidgets.QDoubleSpinBox()
        self.xlabelFontSizeEntry.setSingleStep(0.1)
        self.xlabelFontSizeEntry.setMinimum(0)
        self.xlabelFontSizeEntry.setValue(self.xlabelFontSizeBackup)
        self.xlabelFontSizeEntry.valueChanged.connect(self.updatePlot)
        self.xlabelFontSizeEntry.hide()
        self.fontFrame.addWidget(self.xlabelFontSizeEntry, 3, 1)
        self.ylabelFontSizeBackup = self.ax.yaxis.get_label().get_fontsize()
        self.ylabelFontLabel = wc.QLeftLabel("Y-label:")
        self.ylabelFontLabel.hide()
        self.fontFrame.addWidget(self.ylabelFontLabel, 4, 0)
        self.ylabelFontSizeEntry = QtWidgets.QDoubleSpinBox()
        self.ylabelFontSizeEntry.setSingleStep(0.1)
        self.ylabelFontSizeEntry.setMinimum(0)
        self.ylabelFontSizeEntry.setValue(self.ylabelFontSizeBackup)
        self.ylabelFontSizeEntry.valueChanged.connect(self.updatePlot)
        self.ylabelFontSizeEntry.hide()
        self.fontFrame.addWidget(self.ylabelFontSizeEntry, 4, 1)
        self.xtickFontSizeBackup = self.ax.xaxis.get_ticklabels()[0].get_fontsize()
        self.xtickFontLabel = wc.QLeftLabel("X-ticks:")
        self.xtickFontLabel.hide()
        self.fontFrame.addWidget(self.xtickFontLabel, 5, 0)
        self.xtickFontSizeEntry = QtWidgets.QDoubleSpinBox()
        self.xtickFontSizeEntry.setSingleStep(0.1)
        self.xtickFontSizeEntry.setMinimum(0)
        self.xtickFontSizeEntry.setValue(self.xtickFontSizeBackup)
        self.xtickFontSizeEntry.valueChanged.connect(self.updatePlot)
        self.xtickFontSizeEntry.hide()
        self.fontFrame.addWidget(self.xtickFontSizeEntry, 5, 1)
        self.ytickFontSizeBackup = self.ax.yaxis.get_ticklabels()[0].get_fontsize()
        self.ytickFontLabel = wc.QLeftLabel("Y-ticks:")
        self.ytickFontLabel.hide()
        self.fontFrame.addWidget(self.ytickFontLabel, 6, 0)
        self.ytickFontSizeEntry = QtWidgets.QDoubleSpinBox()
        self.ytickFontSizeEntry.setSingleStep(0.1)
        self.ytickFontSizeEntry.setMinimum(0)
        self.ytickFontSizeEntry.setValue(self.ytickFontSizeBackup)
        self.ytickFontSizeEntry.valueChanged.connect(self.updatePlot)
        self.ytickFontSizeEntry.hide()
        self.fontFrame.addWidget(self.ytickFontSizeEntry, 6, 1)
        self.legend = self.ax.legend()
        if self.legend is not None:              # Fix for matplotlib 2.0, were for contour self.legend becomes None
            if not self.legend.get_texts():
                self.legend.set_visible(False)
                self.legend = None
            else:
                self.legendFontSizeBackup = self.legend.get_texts()[0].get_fontsize()
                self.legendFontLabel = wc.QLeftLabel("Legend:")
                self.legendFontLabel.hide()
                self.fontFrame.addWidget(self.legendFontLabel, 7, 0)
                self.legendFontSizeEntry = QtWidgets.QDoubleSpinBox()
                self.legendFontSizeEntry.setSingleStep(0.1)
                self.legendFontSizeEntry.setMinimum(0)
                self.legendFontSizeEntry.setValue(self.legendFontSizeBackup)
                self.legendFontSizeEntry.valueChanged.connect(self.updatePlot)
                self.legendFontSizeEntry.hide()
                self.fontFrame.addWidget(self.legendFontSizeEntry, 7, 1)
        self.fontGroup.setLayout(self.fontFrame)
        self.optionFrame.addWidget(self.fontGroup, 3, 0)
        # Legend
        if self.legend is not None:
            if self.oldMainWindow.current.__class__.__name__ == 'CurrentMulti':  # If from multiplot
                order = list(self.oldMainWindow.current.viewSettings['extraOffset'])
                order.append(0)
                self.legendOrder = list(np.argsort(order))[::-1]
            elif self.oldMainWindow.current.__class__.__name__ == 'CurrentStacked':
                self.legendOrder = list(np.arange(0, len(self.legend.get_texts())))[::-1]
            else:
                self.legendOrder = list(np.arange(0, len(self.legend.get_texts())))
            try:
                self.legend.set_draggable(True)
            except AttributeError:
                self.legend.draggable(True) # For older Matplotlib versions
            self.legendPos = 'best'
            self.legendTextList = []
            for line in self.legend.get_texts():
                self.legendTextList.append(line.get_text())
            self.legend.set_visible(False)
            self.legendGroup = QtWidgets.QGroupBox('Legend:')
            self.legendGroup.setCheckable(True)
            self.legendGroup.setChecked(False)
            self.legendGroup.toggled.connect(self.updatePlot)
            self.legendFrame = QtWidgets.QGridLayout()
            legendButton = QtWidgets.QPushButton('Legend settings')
            legendButton.clicked.connect(lambda: LegendWindow(self))
            self.legendFrame.addWidget(legendButton, 0, 0)
            self.legendGroup.setLayout(self.legendFrame)
            self.optionFrame.addWidget(self.legendGroup, 4, 0)
        self.xticksToggle = QtWidgets.QCheckBox('X-ticks')
        self.xticksToggle.setChecked(True)
        self.xticksToggle.stateChanged.connect(self.updatePlot)
        self.optionFrame.addWidget(self.xticksToggle, 5, 0)
        self.yticksToggle = QtWidgets.QCheckBox('Y-ticks')
        self.yticksToggle.setChecked(True)
        self.yticksToggle.stateChanged.connect(self.updatePlot)
        self.optionFrame.addWidget(self.yticksToggle, 6, 0)
        self.frameToggle = QtWidgets.QCheckBox('Box')
        self.frameToggle.setChecked(True)
        self.frameToggle.stateChanged.connect(self.updatePlot)
        self.optionFrame.addWidget(self.frameToggle, 7, 0)
        self.inFrame = QtWidgets.QGridLayout()
        self.frame1.addLayout(self.inFrame, 1, 0)
        self.inFrame.addWidget(wc.QLabel("File type:"), 0, 0, 1, 2)
        self.filetypeEntry = QtWidgets.QComboBox()
        self.fileOptions = ['svg', 'png', 'eps', 'jpg', 'pdf']
        self.filetypeEntry.addItems(self.fileOptions)
        self.inFrame.addWidget(self.filetypeEntry, 1, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Close")
        cancelButton.clicked.connect(self.cancel)
        okButton = QtWidgets.QPushButton("&Save")
        okButton.clicked.connect(self.save)
        box = QtWidgets.QDialogButtonBox()
        box.setOrientation(QtCore.Qt.Horizontal)
        box.addButton(cancelButton, QtWidgets.QDialogButtonBox.RejectRole)
        box.addButton(okButton, QtWidgets.QDialogButtonBox.AcceptRole)
        self.inFrame.addWidget(box, 2, 0)
        grid.setColumnStretch(0, 1)
        grid.setRowStretch(0, 1)
        self.optionFrame.setAlignment(QtCore.Qt.AlignTop)
        self.grid = grid
        scroll.setWidget(content)
        scroll.setMinimumWidth(content.sizeHint().width() + scroll.verticalScrollBar().sizeHint().width())
        self.updatePlot()

    def fontCheckChanged(self, val):
        """
        Shows or hides the details of the font settings.

        Parameters
        ----------
        val : bool
            When True the font details are shown, otherwise they are hidden.
        """
        if val:  # If active
            self.mainFontLabel.setEnabled(False)
            self.mainFontSizeEntry.setEnabled(False)
            self.titleFontLabel.show()
            self.titleFontSizeEntry.show()
            self.xlabelFontLabel.show()
            self.xlabelFontSizeEntry.show()
            self.ylabelFontLabel.show()
            self.ylabelFontSizeEntry.show()
            self.xtickFontLabel.show()
            self.xtickFontSizeEntry.show()
            self.ytickFontLabel.show()
            self.ytickFontSizeEntry.show()
            if self.legend is not None:
                self.legendFontLabel.show()
                self.legendFontSizeEntry.show()
        else:
            self.mainFontLabel.setEnabled(True)
            self.mainFontSizeEntry.setEnabled(True)
            self.titleFontLabel.hide()
            self.titleFontSizeEntry.hide()
            self.xlabelFontLabel.hide()
            self.xlabelFontSizeEntry.hide()
            self.ylabelFontLabel.hide()
            self.ylabelFontSizeEntry.hide()
            self.xtickFontLabel.hide()
            self.xtickFontSizeEntry.hide()
            self.ytickFontLabel.hide()
            self.ytickFontSizeEntry.hide()
            if self.legend is not None:
                self.legendFontLabel.hide()
                self.legendFontSizeEntry.hide()
        self.updatePlot()

    def updateLegend(self, *args):
        """
        Updates the figure legend.
        """
        if self.legend is None:
            return
        if self.legendGroup.isChecked():
            orderedLines = [self.ax.lines[x] for x in self.legendOrder]
            orderedLegendText = [self.legendTextList[x] for x in self.legendOrder]
            if self.fontDetailsCheck.checkState():  # If details checked
                size = self.legendFontSizeEntry.value()
            else:
                size = self.mainFontSizeEntry.value()
            self.legend = self.ax.legend(orderedLines, orderedLegendText, framealpha=1.0, loc=self.legendPos, prop={'size': size})
            try:
                self.legend.set_draggable(True)
            except AttributeError:
                self.legend.draggable(True) # For older Matplotlib versions
        else:
            self.legend.set_visible(False)

    def updatePlot(self, *args):
        """
        Updates the plot.
        """
        if not self.xlimLeftEntry.text():
            self.xlimLeftEntry.setText(str(self.xlimBackup[0]))
        if not self.xlimRightEntry.text():
            self.xlimRightEntry.setText(str(self.xlimBackup[1]))
        if not self.ylimLeftEntry.text():
            self.ylimLeftEntry.setText(str(self.ylimBackup[0]))
        if not self.ylimRightEntry.text():
            self.ylimRightEntry.setText(str(self.ylimBackup[1]))
        if self.fontDetailsCheck.checkState():  # If details checked
            self.fig.suptitle(self.titleEntry.text(), fontsize=self.titleFontSizeEntry.value())
            self.ax.set_xlabel(self.xlabelEntry.text(), fontsize=self.xlabelFontSizeEntry.value())
            self.ax.set_ylabel(self.ylabelEntry.text(), fontsize=self.ylabelFontSizeEntry.value())
            self.ax.set_xlim((safeEval(self.xlimLeftEntry.text(), Type='FI'), safeEval(self.xlimRightEntry.text(), Type='FI')))
            self.ax.set_ylim((safeEval(self.ylimLeftEntry.text(), Type='FI'), safeEval(self.ylimRightEntry.text(), Type='FI')))
            self.ax.tick_params(axis='x', labelsize=self.xtickFontSizeEntry.value())
            self.ax.xaxis.get_offset_text().set_fontsize(self.xtickFontSizeEntry.value())
            self.ax.tick_params(axis='y', labelsize=self.ytickFontSizeEntry.value())
            self.ax.yaxis.get_offset_text().set_fontsize(self.ytickFontSizeEntry.value())
        else:
            self.fig.suptitle(self.titleEntry.text(), fontsize=self.mainFontSizeEntry.value())
            self.ax.set_xlabel(self.xlabelEntry.text(), fontsize=self.mainFontSizeEntry.value())
            self.ax.set_ylabel(self.ylabelEntry.text(), fontsize=self.mainFontSizeEntry.value())
            self.ax.set_xlim((safeEval(self.xlimLeftEntry.text(), Type='FI'), safeEval(self.xlimRightEntry.text(), Type='FI')))
            self.ax.set_ylim((safeEval(self.ylimLeftEntry.text(), Type='FI'), safeEval(self.ylimRightEntry.text(), Type='FI')))
            self.ax.tick_params(axis='x', labelsize=self.mainFontSizeEntry.value())
            self.ax.xaxis.get_offset_text().set_fontsize(self.mainFontSizeEntry.value())
            self.ax.tick_params(axis='y', labelsize=self.mainFontSizeEntry.value())
            self.ax.yaxis.get_offset_text().set_fontsize(self.mainFontSizeEntry.value())
            if self.legend is not None:
                self.legend.prop = {'size': self.mainFontSizeEntry.value()}
        xticksbool = self.xticksToggle.isChecked()
        yticksbool = self.yticksToggle.isChecked()
        if self.frameToggle.isChecked():
            self.ax.spines['top'].set_visible(True)
            self.ax.spines['right'].set_visible(True)
            self.ax.spines['bottom'].set_visible(True)
            self.ax.spines['left'].set_visible(True)
        else:
            self.ax.spines['top'].set_visible(False)
            self.ax.spines['right'].set_visible(False)
            if xticksbool:
                self.ax.spines['bottom'].set_visible(True)
            else:
                self.ax.spines['bottom'].set_visible(False)
            if yticksbool:
                self.ax.spines['left'].set_visible(True)
            else:
                self.ax.spines['left'].set_visible(False)
        self.ax.tick_params(bottom=xticksbool, labelbottom=xticksbool, left=yticksbool, labelleft=yticksbool)
        if not xticksbool:
            self.ax.xaxis.get_offset_text().set_visible(False)
        else:
            self.ax.xaxis.get_offset_text().set_visible(True)
        if not yticksbool:
            self.ax.yaxis.get_offset_text().set_visible(False)
        else:
            self.ax.yaxis.get_offset_text().set_visible(True)
        self.updateLegend()
        self.fig.set_size_inches(self.widthEntry.value() / 2.54, self.heightEntry.value() / 2.54)
        self.canvas.draw()
        self.canvas.adjustSize()

    def pickHandler(self, pickEvent):
        """
        The handler for the peak picking in the save figure window.

        Parameters
        ----------
        pickEvent : QEvent
            The peak picking event.
        """
        if pickEvent.mouseevent.dblclick and (pickEvent.mouseevent.button == 1):
            if isinstance(pickEvent.artist, matplotlib.lines.Line2D):
                EditLineWindow(self, pickEvent.artist)

    def get_mainWindow(self):
        """Returns the mainwindow that has been replaced by the save figure window."""
        return self.oldMainWindow

    def kill(self):
        """
        Closes the save figure window.
        """
        for i in reversed(range(self.grid.count())):
            item = self.grid.itemAt(i).widget()
            if item is not None:
                item.deleteLater()
        self.grid.deleteLater()
        self.oldMainWindow.kill()
        del self.fig
        del self.canvas
        self.deleteLater()

    def save(self):
        """
        Asks the user for a filename and saves the plot in the format set in the window.
        """
        self.updatePlot()
        self.fig.set_size_inches(self.widthEntry.value() / 2.54, self.heightEntry.value() / 2.54)
        WorkspaceName = self.father.workspaceNames[self.father.workspaceNum]  # Set name of file to be saved to workspace name to start
        f = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', self.father.lastLocation + os.path.sep + WorkspaceName + '.' + self.fileOptions[self.filetypeEntry.currentIndex()], filter='(*.' + self.fileOptions[self.filetypeEntry.currentIndex()] + ')')
        if isinstance(f, tuple):
            f = f[0]
        if f:
            self.father.lastLocation = os.path.dirname(f)
            dpi = self.dpiEntry.value()
            if dpi is None:
                dpi = self.fig.dpi
            self.fig.savefig(f, format=self.fileOptions[self.filetypeEntry.currentIndex()], dpi=dpi)
            if self.fileOptions[self.filetypeEntry.currentIndex()] == 'svg':
                with open(f) as fd:  # workarround for stroke miter limit
                    s = fd.read()
                with open(f, 'w') as fd:
                    fd.write(s.replace('stroke-miterlimit:100000;', ''))

    def cancel(self):
        """
        Closes the save figure by the cancel button.
        """
        if self.oldMainWindow.current.viewSettings['showTitle']:
            self.fig.suptitle(self.titleBackup, fontsize=self.titleFontSizeBackup)
        else:
            self.fig.suptitle('', fontsize=self.titleFontSizeBackup)
        self.ax.set_xlabel(self.xlabelBackup, fontsize=self.xlabelFontSizeBackup)
        self.ax.set_ylabel(self.ylabelBackup, fontsize=self.ylabelFontSizeBackup)
        self.ax.set_xlim((self.xlimBackup[0], self.xlimBackup[1]))
        self.ax.set_ylim((self.ylimBackup[0], self.ylimBackup[1]))
        self.ax.tick_params(axis='x', labelsize=self.xtickFontSizeBackup)
        self.ax.xaxis.get_offset_text().set_fontsize(self.xtickFontSizeBackup)
        self.ax.tick_params(axis='y', labelsize=self.ytickFontSizeBackup)
        self.ax.yaxis.get_offset_text().set_fontsize(self.ytickFontSizeBackup)
        if self.legend is not None:
            self.legend.set_visible(False)
        self.ax.spines['top'].set_visible(True)
        self.ax.spines['right'].set_visible(True)
        self.ax.spines['bottom'].set_visible(True)
        self.ax.spines['left'].set_visible(True)
        self.ax.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True)
        self.ax.xaxis.get_offset_text().set_visible(True)
        self.ax.yaxis.get_offset_text().set_visible(True)
        self.fig.set_size_inches((self.widthBackup / 2.54, self.heightBackup / 2.54))
        self.grid.deleteLater()
        self.canvas.draw()
        del self.canvas
        del self.fig
        self.father.closeSaveFigure(self.oldMainWindow)
        self.deleteLater()

#####################################################################################################################


class LegendWindow(QtWidgets.QWidget):
    """
    The window with the settings of the plot legend.
    """

    def __init__(self, parent):
        """
        Initializes the legend window.

        Parameters
        ----------
        parent : SaveFigureWindow
            The save figure window from which the legend window is called.
        """
        super(LegendWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle("Legend")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Legend position:"), 0, 0)
        self.posVal = self.father.legendPos
        self.posEntry = wc.QLineEdit(self.posVal, self.preview)
        grid.addWidget(self.posEntry, 1, 0)
        grid.addWidget(wc.QLabel("Legend order:"), 2, 0)
        self.orderVal = self.father.legendOrder
        self.orderEntry = wc.QLineEdit(self.orderVal, self.preview)
        grid.addWidget(self.orderEntry, 3, 0)
        grid.addWidget(wc.QLabel("Legend:"), 4, 0)
        self.father.legendGroup.setChecked(True)
        self.spinBox = QtWidgets.QSpinBox()
        self.spinBox.setMaximum(len(self.father.legendTextList) - 1)
        self.spinBox.valueChanged.connect(self.changeEdit)
        grid.addWidget(self.spinBox, 5, 0)
        self.legendEditList = []
        for i in range(len(self.father.legendTextList)):
            self.legendEditList.append(wc.QLineEdit(self.father.legendTextList[i], self.preview))
            grid.addWidget(self.legendEditList[i], 6, 0)
            self.legendEditList[i].setVisible(False)
        self.legendEditList[0].setVisible(True)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        layout.addWidget(okButton, 1, 1)
        self.show()
        self.setFixedSize(self.size())
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def changeEdit(self, num):
        """
        Changes the line for which the legend label is edited.

        Parameters
        ----------
        num : int
            The line number.
        """
        for i in range(len(self.legendEditList)):
            self.legendEditList[i].setVisible(False)
        self.legendEditList[num].setVisible(True)

    def preview(self, *args):
        """
        Preview the changes made to the legend.
        """
        tmp = copy.deepcopy(self.father.legendTextList)
        order = eval(self.orderEntry.text())
        for i in range(len(self.legendEditList)):
            tmp[i] = self.legendEditList[i].text()
        env = vars(np).copy()
        try:
            inp = eval(self.posEntry.text(), env)
        except Exception:
            inp = self.posEntry.text()
        orderedLines = [self.father.ax.lines[x] for x in order]
        orderedLegendText = [tmp[x] for x in order]
        self.father.ax.legend(orderedLines, orderedLegendText, loc=inp)
        try:
            self.father.legend.set_draggable(True)
        except AttributeError:
            self.father.legend.draggable(True) # For older Matplotlib versions
        self.father.canvas.draw()

    def closeEvent(self, *args):
        """
        Closes the legend window.
        """
        self.deleteLater()
        self.father.updatePlot()

    def applyAndClose(self):
        """
        Applies the changes made to the legend and closes the window.
        """
        for i in range(len(self.legendEditList)):
            self.father.legendTextList[i] = self.legendEditList[i].text()
        self.father.legendOrder = eval(self.orderEntry.text())
        env = vars(np).copy()
        try:
            inp = eval(self.posEntry.text(), env)
        except Exception:
            inp = self.posEntry.text()
        self.father.legendPos = inp
        self.closeEvent()


#####################################################################################################################


class EditLineWindow(QtWidgets.QWidget):
    """
    The window for editing the line styles.
    """

    def __init__(self, parent, line=None):
        """
        Initializes the line edit window.

        Parameters
        ----------
        parent : SaveFigureWindow
            The save figure window from which the line edit window is called.
        line : Line2D, optional
            The line for which to show the line style.
        """
        super(EditLineWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.canvas = self.father.canvas
        self.ax = self.father.ax
        self.lineList = self.ax.lines
        if line is None:
            self.line = line[0]
        else:
            self.line = line
        self.setWindowTitle("Line")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        grid.addWidget(wc.QLabel("Index:"), 0, 0)
        self.indexSpinBox = QtWidgets.QSpinBox()
        self.indexSpinBox.setMaximum(len(self.lineList) - 1)
        self.indexSpinBox.setValue(self.lineList.index(self.line))
        self.indexSpinBox.valueChanged.connect(self.setIndex)
        grid.addWidget(self.indexSpinBox, 1, 0)
        grid.addWidget(wc.QLabel("Line:"), 2, 0)
        colorbutton = QtWidgets.QPushButton("Color", self)
        colorbutton.clicked.connect(self.setColor)
        grid.addWidget(colorbutton, 3, 0)
        grid.addWidget(wc.QLabel("Linewidth:"), 4, 0)
        self.lwSpinBox = QtWidgets.QDoubleSpinBox()
        self.lwSpinBox.setSingleStep(0.1)
        self.lwSpinBox.valueChanged.connect(self.setLineWidth)
        grid.addWidget(self.lwSpinBox, 5, 0)
        grid.addWidget(wc.QLabel("Linestyle:"), 6, 0)
        self.LINESTYLES = ['-', '--', '-.', ':', 'None']
        self.LINENAMES = ['solid', 'dashed', 'dashdot', 'dotted', 'none']
        self.lineDrop = QtWidgets.QComboBox()
        self.lineDrop.addItems(self.LINENAMES)
        self.lineDrop.activated.connect(self.setLineStyle)
        grid.addWidget(self.lineDrop, 7, 0)
        grid.addWidget(wc.QLabel("Marker:"), 8, 0)
        self.MARKERSTYLES = ['o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', '+', 'x', 'd', 'None']
        self.MARKERNAMES = ['circle', 'triangle down', 'triangle up', 'triangle left', 'triangle right', 'tri up', 'tri down', 'tri left', 'tri right', 'octagon', 'square', 'pentagon', 'star', 'hexagon', 'plus', 'x', 'diamond', 'none']
        self.markerDrop = QtWidgets.QComboBox()
        self.markerDrop.addItems(self.MARKERNAMES)
        self.markerDrop.activated.connect(self.setMarker)
        grid.addWidget(self.markerDrop, 9, 0)
        colorbutton = QtWidgets.QPushButton("Facecolor", self)
        colorbutton.clicked.connect(self.setFaceColor)
        grid.addWidget(colorbutton, 10, 0)
        colorbutton = QtWidgets.QPushButton("Edgecolor", self)
        colorbutton.clicked.connect(self.setEdgeColor)
        grid.addWidget(colorbutton, 11, 0)
        grid.addWidget(wc.QLabel("Markersize:"), 12, 0)
        self.msSpinBox = QtWidgets.QDoubleSpinBox()
        self.msSpinBox.setSingleStep(0.1)
        self.msSpinBox.valueChanged.connect(self.setMarkerSize)
        grid.addWidget(self.msSpinBox, 13, 0)
        cancelButton = QtWidgets.QPushButton("&Close")
        cancelButton.clicked.connect(self.cancelAndClose)
        layout.addWidget(cancelButton, 1, 0)
        okButton = QtWidgets.QPushButton("&Apply")
        okButton.clicked.connect(self.apply)
        layout.addWidget(okButton, 1, 1)
        self.setup()
        self.show()

    def setup(self):
        """
        Stores backup information about the styles to restore the figure.
        """
        self.backupColor = colorConverter.to_rgba(self.line.get_color())
        self.backupLineWidth = self.line.get_linewidth()
        self.lwSpinBox.setValue(self.backupLineWidth)
        self.backupLineStyle = self.line.get_linestyle()
        self.lineDrop.setCurrentIndex(self.LINESTYLES.index(self.backupLineStyle))
        self.backupMarker = self.line.get_marker()
        if self.backupMarker == '':
            self.backupMarker = 'None'
        self.markerDrop.setCurrentIndex(self.MARKERSTYLES.index(self.backupMarker))
        self.backupFaceColor = colorConverter.to_rgba(self.line.get_markerfacecolor())
        self.backupEdgeColor = colorConverter.to_rgba(self.line.get_markeredgecolor())
        self.backupMarkerSize = self.line.get_markersize()
        self.msSpinBox.setValue(self.backupMarkerSize)

    def setIndex(self, val):
        """
        Changes the line of which the information is shown.

        Parameters
        ----------
        val : int
            The line number of which to show information.
        """
        self.reset()
        self.line = self.lineList[val]
        self.setup()

    def setColor(self, *args):
        """
        Sets the color of the selected line.
        """
        color = QtGui.QColor()
        color.setRgbF(*self.backupColor)
        color = QtWidgets.QColorDialog.getColor(color, self, 'Color', QtWidgets.QColorDialog.ColorDialogOption(1))
        if not color.isValid():
            return
        self.line.set_color(color.getRgbF())
        self.canvas.draw()

    def setEdgeColor(self, *args):
        """
        Sets the edgecolor of the selected line.
        """
        color = QtGui.QColor()
        color.setRgbF(*self.backupEdgeColor)
        color = QtWidgets.QColorDialog.getColor(color, self, 'Edgecolor', QtWidgets.QColorDialog.ColorDialogOption(1))
        if not color.isValid():
            return
        self.line.set_markeredgecolor(color.getRgbF())
        self.canvas.draw()

    def setFaceColor(self, *args):
        """
        Sets the facecolor of the selected line.
        """
        color = QtGui.QColor()
        color.setRgbF(*self.backupFaceColor)
        color = QtWidgets.QColorDialog.getColor(color, self, 'Facecolor', QtWidgets.QColorDialog.ColorDialogOption(1))
        if not color.isValid():
            return
        self.line.set_markerfacecolor(color.getRgbF())
        self.canvas.draw()

    def setLineWidth(self, val):
        """
        Sets the linewidth of the selected line.

        Parameters
        ----------
        val : float
            Linewidth in points.
        """
        self.line.set_linewidth(val)
        self.canvas.draw()

    def setLineStyle(self, val):
        """
        Sets the line style of the selected line.

        Parameters
        ----------
        val : {'-', '--', '-.', ':', '', (offset, on-off-seq), ...}
            The line style.
        """
        self.line.set_linestyle(self.LINESTYLES[val])
        self.canvas.draw()

    def setMarker(self, val):
        """
        Sets the marker of the selected line.

        Parameters
        ----------
        val : str
            The marker style.
        """
        self.line.set_marker(self.MARKERSTYLES[val])
        self.canvas.draw()

    def setMarkerSize(self, val):
        """
        Sets the marker size of the selected line.

        Parameters
        ----------
        val : float
            The marker size in points.
        """
        self.line.set_markersize(val)
        self.canvas.draw()

    def reset(self):
        """
        Resets the line to the backup values.
        """
        self.line.set_color(self.backupColor)
        self.line.set_linewidth(self.backupLineWidth)
        self.line.set_linestyle(self.backupLineStyle)
        self.line.set_marker(self.backupMarker)
        self.line.set_markeredgecolor(self.backupEdgeColor)
        self.line.set_markerfacecolor(self.backupFaceColor)
        self.line.set_markersize(self.backupMarkerSize)
        self.canvas.draw()

    def closeEvent(self, *args):
        """
        Closes the line edit window.
        """
        self.deleteLater()
        self.father.updatePlot()

    def apply(self):
        """
        Applies the changes to the line.
        """
        self.setup()

    def cancelAndClose(self):
        """
        Cancels the changes to the line and closes the line edit window.
        """
        self.reset()
        self.closeEvent()
