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

try:
    from PyQt4 import QtGui, QtCore
    from PyQt4 import QtGui as QtWidgets
except ImportError:
    from PyQt5 import QtGui, QtCore, QtWidgets
from safeEval import safeEval


class SsnakeTabs(QtWidgets.QTabWidget):
    # A tab widget were tabs can be closed with the middle mouse button

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.MidButton:
            index = self.tabBar().tabAt(event.pos())
            if index >= 0:
                self.tabCloseRequested.emit(index)


class MyEventFilter(QtCore.QObject):

    def __init__(self, root, *args):
        super(MyEventFilter, self).__init__(*args)
        self.root = root

    def eventFilter(self, receiver, event):
        if event.type() == QtCore.QEvent.KeyPress:
            if event.key() == QtCore.Qt.Key_Z:
                if (event.modifiers() & QtCore.Qt.ControlModifier) and (event.modifiers() & QtCore.Qt.ShiftModifier):
                    self.root.redo()
                    return True
                elif event.modifiers() == (QtCore.Qt.ControlModifier):
                    self.root.undo()
                    return True
        return False


class ToolWindows(QtWidgets.QWidget):

    NAME = ""
    PICK = False
    SINGLESLICE = False

    def __init__(self, parent):
        super(ToolWindows, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle(self.NAME)
        self.layout = QtWidgets.QGridLayout(self)
        self.grid = QtWidgets.QGridLayout()
        self.layout.addLayout(self.grid, 0, 0, 1, 2)
        if self.SINGLESLICE:
            self.singleSlice = QtWidgets.QCheckBox("Single slice")
            self.layout.addWidget(self.singleSlice, 1, 0, 1, 2)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        self.layout.addWidget(cancelButton, 2, 0)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        okButton.setFocus()
        self.layout.addWidget(okButton, 2, 1)
        self.show()
        self.setFixedSize(self.size())
        self.father.menuDisable()
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height() - self.geometry().height(), 0, 0)

    def applyFunc(self):
        pass
    
    def applyAndClose(self):
        self.applyFunc()
        self.closeEvent()
        
    def closeEvent(self, *args):
        self.father.current.upd()
        self.father.current.showFid()
        self.father.menuEnable()
        self.father.updAllFrames()
        self.deleteLater()

class SliceValidator(QtGui.QValidator):

    def validate(self, string, position):
        string = str(string)
        try:
            int(safeEval(string))
            return (QtGui.QValidator.Acceptable, string, position)
        except:
            return (QtGui.QValidator.Intermediate, string, position)


class SliceSpinBox(QtWidgets.QSpinBox):

    def __init__(self, parent, minimum, maximum, *args, **kwargs):
        self.validator = SliceValidator()
        super(SliceSpinBox, self).__init__(parent, *args, **kwargs)
        self.setMinimum(minimum)
        self.setMaximum(maximum)
        self.setKeyboardTracking(False)

    def validate(self, text, position):
        return self.validator.validate(text, position)

    def fixup(self, text):
        return self.validator.fixup(text)

    def valueFromText(self, text):
        inp = int(safeEval(str(text)))
        if inp < 0:
            inp = inp + self.maximum() + 1
        return inp

    def textFromValue(self, value):
        inp = int(value)
        if inp < 0:
            inp = inp + self.maximum() + 1
        return str(inp)


class QLabel(QtWidgets.QLabel):

    def __init__(self, parent, *args, **kwargs):
        super(QLabel, self).__init__(parent, *args, **kwargs)
        self.setAlignment(QtCore.Qt.AlignCenter)


class QLeftLabel(QtWidgets.QLabel):

    def __init__(self, parent, *args, **kwargs):
        super(QLeftLabel, self).__init__(parent, *args, **kwargs)
        self.setAlignment(QtCore.Qt.AlignLeft)
        self.setAlignment(QtCore.Qt.AlignVCenter)


class specialProgressBar(QtWidgets.QProgressBar):
    
    def __init__(self):
        super(specialProgressBar, self).__init__()
        self.setAlignment(QtCore.Qt.AlignCenter)
        self._text = 'Inactive'

    def setText(self, text):
        self._text = text

    def text(self):
        return self._text
