#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

import sys
from PyQt5 import QtGui, QtCore, QtWidgets
QT = 5
import os
import math
import loadIsotopes

def safeEval(inp, length=None):
    try:
        return eval(str(inp))
    except Exception:
        return None

SPINNAMES = ['0', '1/2', '1', '3/2', '2',
             '5/2', '3', '7/2', '4', '9/2',
             '5', '11/2', '6', '13/2', '7',
             '15/2', '8', '17/2', '9']
SPINCOLORS = ['white', 'blue', 'orange', 'green', 'yellow',
              'red', 'lime', 'olive', 'lightBlue', 'maroon',
              'gray', 'pink', 'magenta', 'black', 'cyan',
              'crimson', 'navy', 'beige']
GAMMASCALE = 100 / 42.576
ELECTRONSCALE = GAMMASCALE * 28.02495266


isoPath = os.path.dirname(os.path.realpath(__file__)) + os.path.sep + "IsotopeProperties"
isotopes = loadIsotopes.getIsotopeInfo(isoPath)

# Create a list of structures containing the isotope information
ATOMNUM = max(isotopes['atomNum'])
MASTERISOTOPELIST = []
LONGEST = 0
for i in range(ATOMNUM):
    isotopeEntry = {'name': [],
                      'fullName': [],
                      'mass': [],
                      'spin': [],
                      'abundance': [],
                      'gamma': [],
                      'q': [],
                      'freqRatio': [],
                      'refSample': [],
                      'sampleCondition': [],
                      'linewidthFactor': [],
                      'lifetime': [],
                      'sensitivity': []}
    Vals = [j for j,val in enumerate(isotopes['atomNum']) if val==i+1] #Position of all elements with val == i+1
    for elem in Vals:
        isotopeEntry['name'].append(isotopes['name'][elem])
        isotopeEntry['fullName'].append(isotopes['fullName'][elem])
        isotopeEntry['mass'].append(isotopes['atomMass'][elem])
        isotopeEntry['spin'].append(isotopes['spin'][elem])
        isotopeEntry['abundance'].append(isotopes['abundance'][elem])
        isotopeEntry['gamma'].append(isotopes['gamma'][elem])
        isotopeEntry['q'].append(isotopes['q'][elem])
        isotopeEntry['freqRatio'].append(isotopes['freqRatio'][elem])
        isotopeEntry['refSample'].append(isotopes['refSample'][elem])
        isotopeEntry['sampleCondition'].append(isotopes['sampleCondition'][elem])
        isotopeEntry['linewidthFactor'].append(isotopes['linewidthFactor'][elem])
        isotopeEntry['lifetime'].append(isotopes['lifetime'][elem])
        isotopeEntry['sensitivity'].append(isotopes['sensitivity'][elem])

    LONGEST = max(LONGEST, len(Vals))
    if len(Vals) > 0:
        if isotopeEntry['mass'].count(None) == len(Vals): #If all 'None'
            isotopeEntry['mass'] = None
        MASTERISOTOPELIST.append(isotopeEntry)
    else:
        MASTERISOTOPELIST.append(None)
nameList = sorted(set(isotopes['name']))

class PeriodicTable(QtWidgets.QWidget):

    def __init__(self):
        super(PeriodicTable, self).__init__()
        self.freqConst = 6
        self.windowList = []
        self.resetIso()
        self.initUI()
        self.upd()

    def resetIso(self):
        self.isoSelect = [0] * ATOMNUM
        for i in range(ATOMNUM):
            if MASTERISOTOPELIST[i] is not None:
                if MASTERISOTOPELIST[i]['sensitivity'].count(None) == len(MASTERISOTOPELIST[i]['sensitivity']): #If all 'None'
                    self.isoSelect[i] = 0
                else:
                    tmp = [float(x) for x in MASTERISOTOPELIST[i]['sensitivity'] if x is not None]
                    max_value = max(tmp)
                    self.isoSelect[i] = tmp.index(max_value)
            else:
                self.isoSelect[i] = None

    def initUI(self):
        grid = QtWidgets.QGridLayout()
        self.setLayout(grid)
        count1 = 0
        count2 = 0
        groupList = []
        self.labelList = []
        self.freqEditList = []
        self.legendEntries = []
        grid.addWidget(PtQLabel('B<sub>0</sub>[T]:'), 0, 2)
        self.b0Entry = PtQLineEdit()
        self.b0Entry.returnPressed.connect(self.setB0)
        grid.addWidget(self.b0Entry, 0, 3)
        grid.addWidget(PtQLabel('e[GHz]:'), 1, 2)
        self.electronEntry = PtQLineEdit()
        self.electronEntry.returnPressed.connect(self.setElectron)
        grid.addWidget(self.electronEntry, 1, 3)
        self.electronEntry.returnPressed.connect(self.setElectron)
        grid.addWidget(self.electronEntry, 1, 3)
        self.detailsPush = QtWidgets.QPushButton('Details')
        self.detailsPush.clicked.connect(lambda : self.openWindow(None, 0))
        grid.addWidget(self.detailsPush, 1, 4,1,2)
        self.listPush = QtWidgets.QPushButton('List')
        self.listPush.clicked.connect(lambda : self.openList())
        grid.addWidget(self.listPush, 1, 6,1,2)
        grid.addWidget(PtQLabel('Spin:'), 0, 4)
        for i in range(len(SPINNAMES)):
            tmpWidget = QtWidgets.QWidget()
            self.legendEntries.append(PtQLineEdit(SPINNAMES[i]))
            self.legendEntries[-1].setReadOnly(True)
            grid.addWidget(self.legendEntries[-1], 0, i + 5)
            self.legendEntries[-1].hide()
        for i in range(ATOMNUM):
            groupList.append(QtWidgets.QWidget())
            groupList[-1].mouseDoubleClickEvent = lambda arg, i=i: self.openWindow(arg, i)
            grid.addWidget(groupList[-1], count1 + 2, count2)
            count2 += 1
            if count1 == 0 and count2 == 1:
                count2 = 17
            elif (count1 == 1 or count1 == 2) and count2 == 2:
                count2 = 12
            elif (count1 == 5 or count1 == 6) and count2 == 2:
                count1 += 2
            elif (count1 == 7 or count1 == 8) and count2 == 16:
                count2 = 2
                count1 -= 2
            elif count2 == 18:
                count2 = 0
                count1 += 1
            grid2 = QtWidgets.QGridLayout()
            grid2.setSpacing(0)
            grid2.setContentsMargins(0, 0, 0, 0)
            groupList[i].setLayout(grid2)
            cl = self.palette().color(QtGui.QPalette.Base)
            groupList[i].setStyleSheet('background-color: rgb(' + str(cl.red()) +',' + str(cl.green()) + ',' + str(cl.blue()) + ');')
            self.labelList.append(PtQLabel())
            grid2.addWidget(self.labelList[-1], 0, 0)
            self.freqEditList.append(PtQLineEdit())
            self.freqEditList[-1].returnPressed.connect(lambda i=i: self.setFreq(i))
            grid2.addWidget(self.freqEditList[-1], 1, 0)
        self.setWindowTitle('NMR table')
        self.show()

    def upd(self):
        self.updWindows()
        self.b0Entry.setText('%0.2f' % (self.freqConst * GAMMASCALE))
        self.electronEntry.setText('%0.2f' % (self.freqConst * ELECTRONSCALE))
        self.spinSet = set([])
        for i in range(ATOMNUM):
            if MASTERISOTOPELIST[i] is not None:
                if MASTERISOTOPELIST[i]['mass'] is not None:
                    self.labelList[i].setText(str(i + 1) + ': <sup>' + str(int(MASTERISOTOPELIST[i]['mass'][int(self.isoSelect[i])])) + '</sup>' + MASTERISOTOPELIST[i]['name'][int(self.isoSelect[i])])
                    self.freqEditList[i].setText('%0.2f' % (self.freqConst * MASTERISOTOPELIST[i]['freqRatio'][int(self.isoSelect[i])]))
                    color = QtGui.QColor(SPINCOLORS[int(2 * MASTERISOTOPELIST[i]['spin'][int(self.isoSelect[i])])])
                    colorA = QtGui.QColor()
                    colorA.setHsl(color.hslHue(), color.hslSaturation(), 245)
                    self.freqEditList[i].setStyleSheet('border-style: solid;color: rgb(0, 0, 0); border-width: 2px; border-color: rgb' + repr(color.getRgb()[:-1]) + '; background-color: rgb' + repr(colorA.getRgb()[:-1]) + ';')
                    self.spinSet.add(MASTERISOTOPELIST[i]['spin'][int(self.isoSelect[i])])
                else:
                    self.labelList[i].setText(str(i + 1) + ': ' + MASTERISOTOPELIST[i]['name'][int(self.isoSelect[i])])
                    self.freqEditList[i].setStyleSheet('border-style: solid;color: rgb(0, 0, 0); border-width: 2px;' )
            else:
                self.labelList[i].setText(str(i + 1) + ': ')
                self.freqEditList[i].setText('')
                self.freqEditList[i].setStyleSheet('border-style: solid;color: rgb(0, 0, 0); border-width: 2px;' )
        self.updLegend()

    def updLegend(self):
        sortSpinList = [x * 2 for x in sorted(list(self.spinSet))]
        for legendEntry in self.legendEntries:
            legendEntry.hide()
        for i in range(len(sortSpinList)):
            index = int(sortSpinList[i])
            color = QtGui.QColor(SPINCOLORS[index])
            colorA = QtGui.QColor()
            colorA.setHsl(color.hslHue(), color.hslSaturation(), 245)
            self.legendEntries[i].setStyleSheet('color: rgb(0, 0, 0);border-style: solid; border-width: 2px; border-color: rgb' + repr(color.getRgb()[:-1]) + '; background-color: rgb' + repr(colorA.getRgb()[:-1]) + ';')
            self.legendEntries[i].setText(SPINNAMES[index])
            self.legendEntries[i].show()

    def updWindows(self):
        for win in self.windowList:
            win.upd()

    def openWindow(self, event, n):
        self.windowList.append(DetailWindow(self, n))

    def openList(self):
        self.windowList.append(ListWindow(self))

    def removeWindow(self, win):
        self.windowList.remove(win)

    def setFreq(self, n):
        if MASTERISOTOPELIST[n] is not None and MASTERISOTOPELIST[n]['mass'] is not None:
            val = safeEval(self.freqEditList[n].text())
            if val is not None:
                self.freqConst = val / MASTERISOTOPELIST[n]['freqRatio'][int(self.isoSelect[n])]
                self.upd()

    def setB0(self):
        val = safeEval(self.b0Entry.text())
        if val is not None:
            self.freqConst = val / GAMMASCALE
            self.upd()

    def setElectron(self):
        val = safeEval(self.electronEntry.text())
        if val is not None:
            self.freqConst = val / ELECTRONSCALE
            self.upd()

    def closeEvent(self, event):
        for win in self.windowList:
            win.close()
        super(PeriodicTable, self).closeEvent(event)


class DetailWindow(QtWidgets.QWidget):

    def __init__(self, parent, n=0):
        super(DetailWindow, self).__init__()
        self.father = parent
        self.n = n
        self.refAtom = 0
        self.refIso = 0
        self.initUI()
        self.atomSelect(n + 1)
        self.show()

    def initUI(self):
        self.setWindowTitle('Details')
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
        grid2 = QtWidgets.QGridLayout()
        self.grid.addLayout(grid2, 0, 0, 1, 6)
        grid2.addWidget(PtQLabel('N:'), 0, 0)
        self.nSpinBox = QtWidgets.QSpinBox()
        self.nSpinBox.setMinimum(1)
        self.nSpinBox.setMaximum(ATOMNUM)
        self.nSpinBox.setValue(self.n)
        self.nSpinBox.valueChanged.connect(self.atomSelect)
        grid2.addWidget(self.nSpinBox, 0, 1)
        grid2.addWidget(PtQLabel('Element:'), 0, 2)
        self.nameLabel = QtWidgets.QComboBox()
        self.nameLabel.addItems(list(nameList))
        self.nameLabel.currentIndexChanged[str].connect(self.atomSelectName)
        grid2.addWidget(self.nameLabel, 0, 3)
        grid2.addWidget(PtQLabel('Name:'), 0, 4)
        self.fullNameLabel = PtQLabel()
        grid2.addWidget(self.fullNameLabel, 0, 5)
        grid2.setColumnStretch(6, 1)
        self.refLabel = QtWidgets.QComboBox()
        self.refLabel.addItems(list(isotopes['formatName']))
        self.refLabel.currentIndexChanged[str].connect(self.refSelect)
        self.grid.addWidget(self.refLabel, 0, 7)
        self.grid.addWidget(PtQLabel('Mass:'), 1, 1)
        self.grid.addWidget(PtQLabel('Spin:'), 1, 2)
        self.grid.addWidget(PtQLabel('Abundance [%]:'), 1, 3)
        self.grid.addWidget(PtQLabel(u'γ [10<sup>7</sup> rad s<sup>-1</sup> T<sup>-1</sup>]:'), 1, 4)
        self.grid.addWidget(PtQLabel('Q [fm<sup>2</sup>]:'), 1, 5)
        self.grid.addWidget(PtQLabel('Frequency [MHz]:'), 1, 6)
        self.grid.addWidget(PtQLabel('Sensitivity:'), 1, 7)
        self.grid.addWidget(PtQLabel('Sample:'), 1, 8)
        self.grid.addWidget(PtQLabel('Condition:'), 1, 9)
        self.grid.addWidget(PtQLabel('Linewidth:'), 1, 10)
        self.buttongroup = QtWidgets.QButtonGroup()
        self.buttongroup.buttonClicked.connect(self.changeSelect)
        self.radiobuttons = []
        self.massLabels = []
        self.spinLabels = []
        self.abundanceLabels = []
        self.gammaLabels = []
        self.qLabels = []
        self.freqEntries = []
        self.sensLabels = []
        self.sampleLabels = []
        self.conditionLabels = []
        self.linewidthLabels = []
        for i in range(LONGEST):
            self.radiobuttons.append(QtWidgets.QRadioButton())
            self.buttongroup.addButton(self.radiobuttons[-1], i)
            self.grid.addWidget(self.radiobuttons[-1], i + 2, 0)
            self.massLabels.append(PtQLabel())
            self.massLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.massLabels[-1], i + 2, 1)
            self.spinLabels.append(PtQLabel())
            self.spinLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.spinLabels[-1], i + 2, 2)
            self.abundanceLabels.append(PtQLabel())
            self.abundanceLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.abundanceLabels[-1], i + 2, 3)
            self.gammaLabels.append(PtQLabel())
            self.gammaLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.gammaLabels[-1], i + 2, 4)
            self.qLabels.append(PtQLabel())
            self.qLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.qLabels[-1], i + 2, 5)
            self.freqEntries.append(PtQLineEdit())
            self.freqEntries[-1].setFixedWidth(self.freqEntries[-1].sizeHint().width())
            self.freqEntries[-1].returnPressed.connect(lambda i=i: self.setFreq(i))
            self.grid.addWidget(self.freqEntries[-1], i + 2, 6)
            self.sensLabels.append(PtQLabel())
            self.sensLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.sensLabels[-1], i + 2, 7)
            self.sampleLabels.append(PtQLabel())
            self.sampleLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.sampleLabels[-1], i + 2, 8)
            self.conditionLabels.append(PtQLabel())
            self.conditionLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.conditionLabels[-1], i + 2, 9)
            self.linewidthLabels.append(PtQLabel())
            self.linewidthLabels[-1].setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            self.grid.addWidget(self.linewidthLabels[-1], i + 2, 10)
        self.grid.setRowStretch(LONGEST + 2, 1)
        self.grid.setColumnStretch(12, 1)

    def upd(self):
        self.atomSelect(self.n + 1)

    def atomSelectName(self, name):
        for i in range(len(MASTERISOTOPELIST)):
            var = MASTERISOTOPELIST[i]
            if var is not None:
                if name == var['name'][0]:
                    self.nSpinBox.setValue(i + 1)
                    return

    def atomSelect(self, n):
        self.n = n - 1
        atomProp = MASTERISOTOPELIST[self.n]
        if atomProp is None:
            self.display(0)
            return
        index = self.nameLabel.findText(atomProp['name'][0])
        self.nameLabel.setCurrentIndex(index)
        self.fullNameLabel.setText('<b>' + atomProp['fullName'][0] + '</b>')
        if atomProp['mass'] is None:
            self.display(0)
            return
        self.radiobuttons[int(self.father.isoSelect[self.n])].setChecked(True)
        num = len(atomProp['mass'])
        for i in range(num):
            self.massLabels[i].setText(str(int(atomProp['mass'][i])))
            self.spinLabels[i].setText(SPINNAMES[int(2 * atomProp['spin'][i])])
            if atomProp['abundance'][i] == None:
                if atomProp['lifetime'][i] == '-':
                    self.abundanceLabels[i].setText('-')
                else:
                    self.abundanceLabels[i].setText('- [' + atomProp['lifetime'][i] + ']')
            else:
                if atomProp['lifetime'][i] == '-':
                    self.abundanceLabels[i].setText(str(atomProp['abundance'][i]))
                else:
                    self.abundanceLabels[i].setText(str(atomProp['abundance'][i]) + ' [' + atomProp['lifetime'][i] + ']')
            self.gammaLabels[i].setText(str(atomProp['gamma'][i]))
            if atomProp['q'][i] is None:
                self.qLabels[i].setText('-')
            else:
                self.qLabels[i].setText(str(atomProp['q'][i]))
            self.freqEntries[i].setText(str(self.father.freqConst * atomProp['freqRatio'][i]))
            if sys.version_info < (3,):  # check version for possible unicode tricks
                self.sampleLabels[i].setText(atomProp['refSample'][i].decode('utf-8'))
                self.conditionLabels[i].setText(atomProp['sampleCondition'][i].decode('utf-8'))
            else:
                self.sampleLabels[i].setText(atomProp['refSample'][i])
                self.conditionLabels[i].setText(atomProp['sampleCondition'][i])
            if atomProp['linewidthFactor'][i] is None:
                self.linewidthLabels[i].setText('-')
            else:
                self.linewidthLabels[i].setText('%#2.2f' % atomProp['linewidthFactor'][i])
            spin1 = atomProp['spin'][i]
            spin2 = MASTERISOTOPELIST[self.refAtom]['spin'][self.refIso]
            if atomProp['sensitivity'][i] is None or MASTERISOTOPELIST[self.refAtom]['sensitivity'][self.refIso] is None :
                self.sensLabels[i].setText('-')
            else:
                sens = atomProp['sensitivity'][i] / MASTERISOTOPELIST[self.refAtom]['sensitivity'][self.refIso]
                self.sensLabels[i].setText('%0.4g' % sens)
        self.display(num)

    def refSelect(self, name):
        n = isotopes['formatName'].index(name)
        self.refAtom = int(isotopes['atomNum'][n] - 1)
        if MASTERISOTOPELIST[self.refAtom]['mass'] is not None:
            self.refIso = MASTERISOTOPELIST[self.refAtom]['mass'].index(isotopes['atomMass'][n])
        self.upd()

    def changeSelect(self, n):
        self.father.isoSelect[self.n] = self.buttongroup.checkedId()
        self.father.upd()

    def setFreq(self, i):
        val = safeEval(self.freqEntries[i].text())
        if val is not None:
            self.father.freqConst = val / MASTERISOTOPELIST[self.n]['freqRatio'][i]
            self.father.upd()
            self.atomSelect(self.n + 1)

    def display(self, num):
        for i in range(num):
            self.radiobuttons[i].show()
            self.massLabels[i].show()
            self.spinLabels[i].show()
            self.abundanceLabels[i].show()
            self.gammaLabels[i].show()
            self.qLabels[i].show()
            self.freqEntries[i].show()
            self.sampleLabels[i].show()
            self.conditionLabels[i].show()
            self.linewidthLabels[i].show()
            self.sensLabels[i].show()
        for i in range(num, LONGEST):
            self.radiobuttons[i].hide()
            self.massLabels[i].hide()
            self.spinLabels[i].hide()
            self.abundanceLabels[i].hide()
            self.gammaLabels[i].hide()
            self.qLabels[i].hide()
            self.freqEntries[i].hide()
            self.sampleLabels[i].hide()
            self.conditionLabels[i].hide()
            self.linewidthLabels[i].hide()
            self.sensLabels[i].hide()

    def closeEvent(self, event):
        super(DetailWindow, self).closeEvent(event)
        self.father.removeWindow(self)

class tableItem(QtWidgets.QTableWidgetItem):

    def __init__(self, parent, *args, **kwargs):
        super(tableItem, self).__init__(parent, *args, **kwargs)
        self.setFlags(QtCore.Qt.ItemIsEnabled)

class ListWindow(QtWidgets.QWidget):

    def __init__(self, parent):
        super(ListWindow, self).__init__()
        self.setWindowTitle('List')
        self.father = parent

        grid = QtWidgets.QGridLayout(self)


        grid.addWidget(PtQLabel('Order:'), 0, 0)
        self.orderType = QtWidgets.QComboBox()
        self.orderType.addItems(['Element number','Frequency (Descending)','Frequency (Ascending)',  'Spin (Descending)', 'Spin (Ascending)','Q (Descending)', 'Q (Ascending)',
            'Abundance (Descending)', 'Abundance (Ascending)','Sensitivity (Descending)', 'Sensitivity (Ascending)'])
        self.orderType.currentIndexChanged.connect(self.upd)
        grid.addWidget(self.orderType,0,1)

        self.table = QtWidgets.QTableWidget(1,8)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.table.verticalHeader().setVisible(False)
        self.table.setHorizontalHeaderLabels(['Nucleus','Name','Spin','Abundance [%]','Q [fm^2]','Frequency ratio [%]',
            'Frequency [MHz]','Sensitivity [1H]'])

        grid.addWidget(self.table, 1, 0,1,6)
        #grid.setRowStretch(3, 2)
        grid.setColumnStretch(5, 1)
        self.spinName = ['1/2','1','3/2','2','5/2','3','7/2','4','9/2','5','6','7']
        self.spinVals = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7]
        self.upd()
        self.resize(900, 800)
        self.show()

    def upd(self):
        count = 0
        isotopes = []
        for elem in MASTERISOTOPELIST:
            if elem['mass'] is not None:
                for i in range(len(elem['mass'])):
                    isotopes.append( {'mass':elem['mass'][i],  'name': elem['name'][i], 'fullName':elem['fullName'][i],
                        'spin': elem['spin'][i], 'abundance': elem['abundance'][i],'q': elem['q'][i],'freqRatio': elem['freqRatio'][i],
                        'gamma': elem['gamma'][i],})

                    if elem['sensitivity'][i] is None or MASTERISOTOPELIST[0]['sensitivity'][0] is None:
                        sens = {'sens': None}
                    else:
                        sensTmp = elem['sensitivity'][i] / MASTERISOTOPELIST[0]['sensitivity'][0]
                        sens = {'sens': sensTmp}
                    isotopes[-1].update(sens)

        orderType = self.orderType.currentIndex()
        if orderType != 0:
            orderings = [['freqRatio', True],['freqRatio', False],['spin', True],['spin', False],['q', True],['q', False],['abundance', True],['abundance', False],
                    ['sens', True],['sens', False]]
            actions = orderings[orderType - 1]

            tmp = []
            for elem in isotopes:
                tmp.append(elem[actions[0]])
            tmp = [-1e32 if x is None else x for x in tmp] #Make None to be the lowest value for the sort
            order = sorted(range(len(tmp)), key=tmp.__getitem__) #argsort without numpy
            if actions[1]:
                order = order[::-1]
            isotopes = [isotopes[x] for x in order]

        for elem in isotopes:
            self.table.setRowCount(count + 1)
            self.table.setItem(count,0,tableItem(str(int(elem['mass'])) + elem['name']))
            self.table.setItem(count,1,tableItem(elem['fullName']))
            Spin = elem['spin']
            Spin = self.spinName[self.spinVals.index(Spin)]
            self.table.setItem(count,2,tableItem(Spin))
            Abun = elem['abundance']
            if Abun is None:
                Abun = '-'
            else:
                Abun = str(Abun)
            self.table.setItem(count,3,tableItem(Abun))
            Q = elem['q']
            if Q is None:
                Q = '-'
            else:
                Q = str(Q)
            self.table.setItem(count,4,tableItem(Q))
            self.table.setItem(count,5,tableItem(str(elem['freqRatio'])))
            self.table.setItem(count,6,tableItem(str(elem['freqRatio'] * self.father.freqConst)))
            sens = elem['sens']
            if sens is None:
                sens = '-'
            else:
                sens = str(sens)
            self.table.setItem(count,7,tableItem(sens))
            count += 1

    def closeEvent(self, event):
        super(ListWindow, self).closeEvent(event)
        self.father.removeWindow(self)

class PtQLabel(QtWidgets.QLabel):
    def __init__(self, parent=None):
        super(PtQLabel, self).__init__(parent)
        self.setAlignment(QtCore.Qt.AlignCenter)


class PtQLineEdit(QtWidgets.QLineEdit):
    def __init__(self, parent=None):
        super(PtQLineEdit, self).__init__(parent)
        self.setAlignment(QtCore.Qt.AlignCenter)


def main():
    app = QtWidgets.QApplication(sys.argv)
    ex = PeriodicTable()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
