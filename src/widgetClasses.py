#!/usr/bin/env python

# Copyright 2016 Bas van Meerten and Wouter Franssen

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
        QtCore.QObject.__init__(self, *args)
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
        QtWidgets.QDoubleSpinBox.__init__(self, parent, *args, **kwargs)
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
        QtWidgets.QLabel.__init__(self, parent, *args, **kwargs)
        self.setAlignment(QtCore.Qt.AlignCenter)


class specialProgressBar(QtWidgets.QProgressBar):
    def __init__(self):
        QtWidgets.QProgressBar.__init__(self)
        self.setAlignment(QtCore.Qt.AlignCenter)
        self._text = 'Inactive'

    def setText(self, text):
        self._text = text

    def text(self):
        return self._text