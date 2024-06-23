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

import os
import sys
from safeEval import safeEval
from ssNake import QtGui, QtCore, QtWidgets

class SsnakeTabs(QtWidgets.QTabWidget):
    """
    A reimplementation of the PyQt QTabWidget.
    A tab widget were tabs can be closed with the middle mouse button.
    """

    def mousePressEvent(self, event):
        """
        Reimplementation from QTabWidget.
        Middle mousebutton closes the tab.

        Parameters
        ----------
        event : QMouseEvent
            The mouse event.
        """
        if event.button() == QtCore.Qt.MidButton:
            index = self.tabBar().tabAt(event.pos())
            if index >= 0:
                self.tabCloseRequested.emit(index)

class SsnakeTreeWidget(QtWidgets.QTreeView):
    """
    A reimplementation of the PyQt QTreeView.
    Allows the loading of files and directories in ssNake.
    """

    def __init__(self, parent):
        """
        Initializes the SsnakeTreeWidget.

        Parameters
        ----------
        parent : QWidget
            Parent of the treeview.
        """
        super(SsnakeTreeWidget, self).__init__(parent)
        self.father = parent
        self.dirmodel = QtWidgets.QFileSystemModel()
        self.dirmodel.setRootPath('')
        # Don't show files, just folders
        self.dirmodel.setFilter(QtCore.QDir.NoDotAndDotDot | QtCore.QDir.AllDirs| QtCore.QDir.Files| QtCore.QDir.Drives)
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.openMenu)
        self.setModel(self.dirmodel)
        self.setRootIndex(self.dirmodel.index(''))
        self.header().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
        self.header().setStretchLastSection(False)
        self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        #self.setToolTip(index.model()->data(index,Qt:isplayRole).toString());
        # Don't show columns for size, file type, and last modified
        self.setHeaderHidden(True)
        self.hideColumn(1)
        self.hideColumn(2)
        self.hideColumn(3)
        self.expand_all(self.dirmodel.index(self.father.lastLocation))

    def mouseDoubleClickEvent(self, event):
        """
        Reimplementation from QTreeView.
        Middle mousebutton loads the data.
        Left mousebutton loads the data from a file.

        Parameters
        ----------
        event : QMouseEvent
            The mouse event.
        """
        index = self.indexAt(event.pos())
        path = self.dirmodel.filePath(index)
        if event.button() == QtCore.Qt.MidButton:
            self.father.loadData([path])
        elif event.button() == QtCore.Qt.LeftButton and not self.dirmodel.isDir(index):
            self.father.loadData([path])
        super(SsnakeTreeWidget, self).mouseDoubleClickEvent(event)

    def mousePressEvent(self, event):
        """
        Reimplementation from QTreeView.
        Middle mousebutton loads the data.

        Parameters
        ----------
        event : QMouseEvent
            The mouse event.
        """
        if event.button() == QtCore.Qt.MidButton:
            index = self.indexAt(event.pos())
            path = self.dirmodel.filePath(index)
            self.father.loadData([path])
        else: #If not, let the QTreeView handle the event
            super(SsnakeTreeWidget, self).mousePressEvent(event)

    def expand_all(self, index):
        """
        Expand all folders up to a certain index.

        Parameters
        ----------
        event : QModelIndex
            The index to which the folder should be expanded in the treeview.
        """
        path = self.dirmodel.filePath(index)
        pathOld = '-1'
        while pathOld != path:
            self.setExpanded(self.dirmodel.index(path), True)
            pathOld = path
            path = os.path.dirname(path)

    def openMenu(self, position):
        """
        Creates the contextmenu of the treeview at a given position.
        The contents of the contextmenu depend on whether the selected item is a file or a directory.

        Parameters
        ----------
        position : QPoint
            The position to place the contextmenu.
        """
        index = self.selectedIndexes()
        path = [self.dirmodel.filePath(x) for x in index]
        menu = QtWidgets.QMenu()
        if len(path) == 1:
            if self.dirmodel.isDir(index[0]):
                menu.addAction("Load Directory", lambda: self.father.loadData(path))
            else:
                menu.addAction("Load File", lambda: self.father.loadData(path))
                menu.addAction("Open File Externally", lambda: self.openExtAct(path))
            menu.addAction("Open in File Browser", lambda: self.openBrowser(path))
        else:
            menu.addAction("Load Selection", lambda: self.father.loadData(path))
        menu.exec_(self.viewport().mapToGlobal(position))

    def openExtAct(self, fileNames):
        """
        Opens a file externally.

        Parameters
        ----------
        fileNames : list of str
            The first string from the list is used as a path to open externally.
        """
        fileName = fileNames[0]
        if sys.platform.startswith('linux'):
            os.system("xdg-open " + '"' + fileName + '"')
        elif sys.platform.startswith('darwin'):
            os.system("open " + '"' + fileName + '"')
        elif sys.platform.startswith('win'):
            os.startfile(fileName)

    def openBrowser(self, path):
        """
        Opens a directory using the filebrowser.
        When the a file is given the containing directory is used.

        Parameters
        ----------
        path : list of str
            The first string from the list is used as the path to open in a filebrowser.
        """
        path = os.path.dirname(path[0])
        if sys.platform.startswith('linux'):
            os.system("xdg-open " + '"' + path + '"')
        elif sys.platform.startswith('darwin'):
            os.system("open " + '"' + path + '"')
        elif sys.platform.startswith('win'):
            os.startfile(path)


class SsnakeSlider(QtWidgets.QSlider):
    """
    A reimplementation of the PyQt QSlider.
    The behaviour is modified when the Control or Shift buttons are pressed.
    """

    def wheelEvent(self, event):
        """
        Reimplementation from QSlider.
        When the Control button is held the stepsize is multiplied by 10.
        When the Shift button is held the stepsize is multiplied by 100.
        When both Control and Shift buttons are held the stepsize is multiplied by 1000.

        Parameters
        ----------
        event : QMouseEvent
            The mouse event.
        """
        delta = event.angleDelta().y()
        step = self.singleStep() * 3
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier and QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            step *= 1000
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            step *= 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            step *= 100
        if delta > 0:
            self.setValue(self.value() + step)
        else:
            self.setValue(self.value() - step)


class SplitterEventFilter(QtCore.QObject):
    """
    The event filter for the splitter between the treeview and the canvas.
    """

    def __init__(self, root, *args):
        """
        Initializes the splitter event filter.

        Parameters
        ----------
        root : QSplitter
            The splitter to which the event filter is installed.
        *args
            Additional arguments are passed to QObject.
        """
        super(SplitterEventFilter, self).__init__(*args)
        self.root = root
        self.sizeBak = 0

    def eventFilter(self, receiver, event):
        """
        The event filter function.
        A double mouseclick or a single middle mouseclick minimizes or restores the treeview.

        Parameters
        ----------
        receiver : QObject
            The receiver.
            Not used.
        event : QEvent
            The event.

        Returns
        -------
        bool
            True if the event matched the resize event, False otherwise.
        """
        select = False
        if event.type() == QtCore.QEvent.MouseButtonDblClick:
            select = True
        # If single click with middle mouse button
        if event.type() == QtCore.QEvent.MouseButtonPress:
            if event.button() == QtCore.Qt.MidButton:
                select = True
        if select:
            sizes = self.root.sizes()
            if sizes[0] == 0:
                self.root.setSizes([self.sizeBak, 1])
            else:
                self.sizeBak = sizes[0]
                self.root.setSizes([0, 1])
            return True
        return False

class SsnakeEventFilter(QtCore.QObject):
    """
    The event filter of ssNake for undo and redo.
    """

    def __init__(self, root, *args):
        """
        Initializes the ssNake event filter.

        Parameters
        ----------
        root : MainProgram
            The main program to which the event filter is installed.
        *args
            Additional arguments are passed to QObject.
        """
        super(SsnakeEventFilter, self).__init__(*args)
        self.root = root

    def eventFilter(self, receiver, event):
        """
        The event filter function.
        Control + Z triggers the undo function.
        Control + Shift + Z triggers the redo function.

        Parameters
        ----------
        receiver : QObject
            The receiver.
            Not used.
        event : QEvent
            The event.

        Returns
        -------
        bool
            True if the event matched undo or redo, False otherwise.
        """
        if event.type() == QtCore.QEvent.KeyPress:
            if event.key() == QtCore.Qt.Key_Z:
                if (event.modifiers() & QtCore.Qt.ControlModifier) and (event.modifiers() & QtCore.Qt.ShiftModifier):
                    self.root.redo()
                    return True
                if event.modifiers() == (QtCore.Qt.ControlModifier):
                    self.root.undo()
                    return True
        return False


class ToolWindow(QtWidgets.QWidget):
    """
    The base class of the toolwindows.
    Implements the basic features shared by all toolwindows.
    Toolwindows inherit this class and alter its behaviour by the constants or by reimplementing functions.
    """

    NAME = ""              # The name displayed in the title of the window
    PICK = False           # Does the window use peak picking
    SINGLESLICE = False    # Should the single slice button be displayed
    BROWSE = False         # Should the window have a browse button
    RESIZABLE = False      # should the window be resizable
    MENUDISABLE = True     # Should the window disable the menu of the main window
    APPLYANDCLOSE = True   # Should the window close after the ok button is pressed
    CANCELNAME = "&Cancel" # The name on the cancel button
    OKNAME = "&Ok"         # The name on the ok button

    def __init__(self, parent):
        """
        Initializes the ToolWindow.

        Parameters
        ----------
        parent : Main1DWindow or AbstractParamFrame
            Parent of the toolwindow.
        """
        super(ToolWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        self.setWindowTitle(self.NAME)
        self.layout = QtWidgets.QGridLayout(self)
        self.grid = QtWidgets.QGridLayout()
        self.box = QtWidgets.QDialogButtonBox()
        if self.BROWSE:
            self.layout.addLayout(self.grid, 0, 0, 1, 3)
        else:
            self.layout.addLayout(self.grid, 0, 0, 1, 2)
        if self.SINGLESLICE:
            self.singleSlice = QtWidgets.QCheckBox("Single slice")
            self.layout.addWidget(self.singleSlice, 1, 0, 1, 2)
        if self.BROWSE:
            self.browseButton = QtWidgets.QPushButton("&Browse")
            self.browseButton.clicked.connect(self.browse)
            self.box.addButton(self.browseButton, QtWidgets.QDialogButtonBox.ActionRole)
        self.cancelButton = QtWidgets.QPushButton(self.CANCELNAME)
        self.cancelButton.clicked.connect(self.closeEvent)
        self.okButton = QtWidgets.QPushButton(self.OKNAME)
        self.okButton.clicked.connect(self.applyAndClose)
        self.okButton.setFocus()
        self.box.addButton(self.cancelButton, QtWidgets.QDialogButtonBox.RejectRole)
        self.box.addButton(self.okButton, QtWidgets.QDialogButtonBox.AcceptRole)
        self.show()
        self.layout.addWidget(self.box, 3, 0)
        if not self.RESIZABLE:
            self.setFixedSize(self.size())
        if self.MENUDISABLE:
            self.father.menuEnable(False)
        self.setGeometry(self.frameSize().width() - self.geometry().width(), self.frameSize().height(), 0, 0)

    def browse(self):
        """
        Dummy function for the Browse button.
        Should be reimplemented by the toolwindows using BROWSE=True.
        """
        pass

    def applyFunc(self):
        """
        Dummy function for the apply function.
        Should be reimplemented by the toolwindows.
        """
        pass

    def applyAndClose(self):
        """
        Runs the apply function and when APPLYANDCLOSE is set, closes the window.
        """
        self.applyFunc()
        if self.APPLYANDCLOSE:
            self.closeEvent()

    def closeEvent(self, *args):
        """
        Updates the view and closes the toolwindow.

        Parameters
        ----------
        *args
            Any arguments are ignored.
        """
        self.father.current.upd()
        self.father.current.showFid()
        if self.MENUDISABLE:
            self.father.menuEnable(True)
        self.father.updAllFrames()
        self.deleteLater()


class SliceValidator(QtGui.QValidator):
    """
    A reimplementation of the QValidator.
    It uses the safeEval to validate a string for slice selection.
    """

    def validate(self, string, position):
        """
        Validates a given string using safeEval.

        Parameters
        ----------
        string : str
            String to be validated.
        position : int
            Position of the cursor.

        Returns
        -------
        State
            Acceptable if the string is parsable by safeEval, Intermediate otherwise.
        string : str
            The input string.
        position : int
            The input position of the cursor.
        """
        string = str(string)
        try:
            int(safeEval(string))
            return (QtGui.QValidator.Acceptable, string, position)
        except Exception:
            return (QtGui.QValidator.Intermediate, string, position)


class SliceSpinBox(QtWidgets.QSpinBox):
    """
    A reimplementation of the QSpinBox.
    This spinbox is designed for the slice selection of spectra.
    """

    def __init__(self, parent, minimum, maximum, *args, **kwargs):
        """
        Initializes the SliceSpinBox.

        Parameters
        ----------
        parent : SideFrame
            The sideframe which contains the SliceSpinBox.
        minimum : int
            The minimum value of the spinbox.
        maximum : int
            The maximum value of the spinbox.
        *args
            Additional arguments are passed to QSpinBox.
        **kwargs
            Keyword arguments are passed to QSpinBox
        """
        self.validator = SliceValidator()
        self.validate = self.validator.validate
        self.fixup = self.validator.fixup
        super(SliceSpinBox, self).__init__(parent, *args, **kwargs)
        self.setMinimum(minimum)
        self.setMaximum(maximum)
        self.setKeyboardTracking(False)

    def valueFromText(self, text):
        """
        Parses a string to a slice number.

        Parameters
        ----------
        text : str
            String to be parsed.

        Returns
        -------
        int
            The slice number.
        """
        inp = int(safeEval(str(text)))
        if inp < 0:
            inp = inp + self.maximum() + 1
        return inp

    def textFromValue(self, value):
        """
        Parses a value to a slice number.

        Parameters
        ----------
        text : int or float
            Value to be parsed.

        Returns
        -------
        int
            The slice number.
        """
        inp = int(value)
        if inp < 0:
            inp = inp + self.maximum() + 1
        return str(inp)

    def stepBy(self, steps):
        """
        Increases or decreases the spinbox by a given value.
        When the Control button is held the stepsize is multiplied 10.
        When the Shift button is held the stepsize is multiplied 100.
        When both the Control and Shift buttons are held the stepsize is multiplied 1000.

        Parameters
        ----------
        steps : int
            The number of steps to increase the spinbox
        """
        mod = 1
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier and QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            mod = 1000
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            mod = 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            mod = 100
        self.setValue(self.value() + mod * steps)

    def wheelEvent(self, event):
        """
        The function for handling the scroll event on the spinbox.

        Parameters
        ----------
        event : QWheelEvent
            The event on the spinbox.
        """
        delta = event.angleDelta().y() + event.angleDelta().x()
        step = 1
        if delta > 0:
            self.stepBy(step)
        else:
            self.stepBy(-step)
        event.accept()


class SsnakeDoubleSpinBox(QtWidgets.QDoubleSpinBox):
    """
    A reimplementation of the QDoubleSpinBox.
    """

    def stepBy(self, steps):
        """
        Increases or decreases the spinbox by a given value.
        When the Control button is held the stepsize is multiplied 10.
        When the Shift button is held the stepsize is multiplied 100.
        When both the Control and Shift buttons are held the stepsize is multiplied 1000.
        When the Alt button is held the inverse

        Parameters
        ----------
        steps : int or float
            The value to increase the spinbox
        """
        mod = 1
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier and QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            mod = 1000
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            mod = 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            mod = 100
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.AltModifier:
            self.setValue(self.value() + self.singleStep() * steps / mod)
        else:
            self.setValue(self.value() + self.singleStep() * mod * steps)

    def wheelEvent(self, event):
        """
        The function for handling the scroll event on the spinbox.

        Parameters
        ----------
        event : QWheelEvent
            The event on the spinbox.
        """
        delta = event.angleDelta().y() + event.angleDelta().x()
        step = 1
        if delta > 0:
            self.stepBy(step)
        else:
            self.stepBy(-step)
        event.accept()


class SsnakeSpinBox(QtWidgets.QSpinBox):
    """
    A reimplementation of the QSpinBox.
    """

    def stepBy(self, steps):
        """
        Increases or decreases the spinbox by a given value.
        When the Control button is held the stepsize is multiplied 10.
        When the Shift button is held the stepsize is multiplied 100.
        When both the Control and Shift buttons are held the stepsize is multiplied 1000.

        Parameters
        ----------
        steps : int
            The value to increase the spinbox
        """
        mod = 1
        if QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier and QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            mod = 1000
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ControlModifier:
            mod = 10
        elif QtWidgets.qApp.keyboardModifiers() & QtCore.Qt.ShiftModifier:
            mod = 100
        self.setValue(self.value() + self.singleStep() * mod * steps)

    def wheelEvent(self, event):
        """
        The function for handling the scroll event on the spinbox.

        Parameters
        ----------
        event : QWheelEvent
            The event on the spinbox.
        """
        delta = event.angleDelta().y() + event.angleDelta().x()
        step = 1
        if delta > 0:
            self.stepBy(step)
        else:
            self.stepBy(-step)
        event.accept()


class QLabel(QtWidgets.QLabel):
    """
    A reimplementation of the QLabel.
    The text is center aligned by default.
    """

    def __init__(self, parent, *args, **kwargs):
        super(QLabel, self).__init__(parent, *args, **kwargs)
        self.setAlignment(QtCore.Qt.AlignCenter)


class QSelectLabel(QtWidgets.QLabel):
    """
    A reimplementation of the QLabel.
    The text is selectable by default.
    """

    def __init__(self, parent, *args, **kwargs):
        super(QSelectLabel, self).__init__(parent, *args, **kwargs)
        self.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)


class QLeftLabel(QtWidgets.QLabel):
    """
    A reimplementation of the QLabel.
    The text is left aligned by default.
    """

    def __init__(self, parent, *args, **kwargs):
        super(QLeftLabel, self).__init__(parent, *args, **kwargs)
        self.setAlignment(QtCore.Qt.AlignLeft)
        self.setAlignment(QtCore.Qt.AlignVCenter)


class QLineEdit(QtWidgets.QLineEdit):
    """
    A reimplementation of the QLineEdit.
    The text is center aligned by default.
    """

    def __init__(self, text='', func=None, parent=None):
        """
        Initializes QLineEdit.

        Parameters
        ----------
        text : str, optional
            The initial text, by default this is an empty string.
        func : function, optional
            The function that is called when return is pressed.
            By default no function is called.
        parent : QWidget, optional
            The parent widget of the lineedit.
            This object is passed to QLineEdit.
            By default None is used.
        """
        super(QLineEdit, self).__init__(parent)
        self.setText(str(text))
        self.setAlignment(QtCore.Qt.AlignCenter)
        if func is not None:
            self.returnPressed.connect(func)

    def setText(self,*args, **kwargs):
        super(QLineEdit, self).setText(*args,**kwargs)
        self.setCursorPosition(0)


class FitQLineEdit(QLineEdit):
    """
    A reimplementation of the QLineEdit designed for the fitting parameter frame.
    Allows connecting parameters from the contextmenu.
    """

    def __init__(self, fitParent, paramName, *args):
        """
        Initializes FitQLineEdit.

        Parameters
        ----------
        parent : AbstractParamFrame
            The fitting parameter frame of the lineedit.
        paramName : str
            The name of this lineedit.
        *args
            Additional arguments are passed to QLineEdit.
        """
        super(FitQLineEdit, self).__init__(*args)
        self.fitParent = fitParent
        self.paramName = paramName

    def contextMenuEvent(self, event):
        """
        Creates the context menu, with the Connect Parameter option.

        Parameters
        ----------
        event : QEvent
            The event.
        """
        menu = self.createStandardContextMenu()
        menu.addAction('Connect Parameter', self.connectParams)
        menu.exec_(event.globalPos())

    def connectParams(self, *args):
        """
        Opens the ConnectParamsWindow.

        Parameters
        ----------
        *args
            Additional arguments are ignored.
        """
        parametertxtlist = self.fitParent.rootwindow.getParamTextList()
        ConnectParamsWindow(self, parametertxtlist, self.paramName, self.fitParent.rootwindow.getTabNames(), self.fitParent.rootwindow.getCurrentTabName(), self.setConnect)

    def setConnect(self, inpTuple):
        self.setText(str(inpTuple))


class ConnectParamsWindow(QtWidgets.QWidget):
    """
    The window for connecting fitting parameters.
    """

    def __init__(self, parent, paramTextList, paramName, spectrumNames, currentSpectrum, returnFunc):
        """
        Initializes the ConnectParamsWindow.

        Parameters
        ----------
        parent : QWidget
            The parent of the ConnectParamsWindow. This value is passed to QWidget.
        paramTextList : list of str
            List with the names of all parameters.
        paramName : str
            Name of the lineedit from which this window was opened.
        spectrumNames : list of str
            A list of the names of the spectra being fit.
        currentSpectrum : str
            The name of the spectrum that is currently open.
        returnFunc : function
            The function that should be called when the ok button is pressed.
        """
        super(ConnectParamsWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.setWindowTitle("Connect Parameter")
        self.paramTextList = paramTextList
        self.paramText = paramName
        self.spectrumNames = spectrumNames
        self.currentSpectrum = self.spectrumNames.index(currentSpectrum)
        self.returnFunc = returnFunc
        self.layout = QtWidgets.QGridLayout(self)
        self.grid = QtWidgets.QGridLayout()
        self.layout.addLayout(self.grid, 0, 0, 1, 2)
        self.grid.addWidget(QtWidgets.QLabel("Parameter:"), 0, 0)
        self.paramNameEntry = QtWidgets.QComboBox()
        self.paramNameEntry.addItems(self.paramTextList)
        self.paramNameEntry.setCurrentIndex(self.paramTextList.index(self.paramText))
        self.grid.addWidget(self.paramNameEntry, 0, 1)
        self.grid.addWidget(QtWidgets.QLabel("Data:"), 1, 0)
        self.spectrumNameEntry = QtWidgets.QComboBox()
        self.spectrumNameEntry.addItems(self.spectrumNames)
        self.spectrumNameEntry.setCurrentIndex(self.currentSpectrum)
        self.grid.addWidget(self.spectrumNameEntry, 1, 1)
        self.grid.addWidget(QtWidgets.QLabel("Line:"), 2, 0)
        self.lineEntry = QtWidgets.QSpinBox()
        self.grid.addWidget(self.lineEntry, 2, 1)
        self.grid.addWidget(QtWidgets.QLabel("Multiplier:"), 3, 0)
        self.multEntry = QLineEdit("1.0")
        self.grid.addWidget(self.multEntry, 3, 1)
        self.grid.addWidget(QtWidgets.QLabel("Offset:"), 4, 0)
        self.addEntry = QLineEdit("0.0")
        self.grid.addWidget(self.addEntry, 4, 1)
        self.cancelButton = QtWidgets.QPushButton("&Cancel")
        self.cancelButton.clicked.connect(self.closeEvent)
        self.okButton = QtWidgets.QPushButton("&Ok")
        self.okButton.clicked.connect(self.applyAndClose)
        self.okButton.setFocus()
        self.box = QtWidgets.QDialogButtonBox()
        self.box.addButton(self.cancelButton, QtWidgets.QDialogButtonBox.RejectRole)
        self.box.addButton(self.okButton, QtWidgets.QDialogButtonBox.AcceptRole)
        self.layout.addWidget(self.box, 2, 0)
        self.show()

    def applyAndClose(self):
        """
        Runs the returnFunc and closes the window.
        """
        paramName = self.paramTextList[self.paramNameEntry.currentIndex()]
        returnTuple = (paramName, self.lineEntry.value(), safeEval(self.multEntry.text()), safeEval(self.addEntry.text()), self.spectrumNameEntry.currentIndex())
        self.closeEvent()
        self.returnFunc(returnTuple)

    def closeEvent(self, *args):
        self.deleteLater()
