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
import json
import zipfile
import tempfile
import shutil
if sys.version_info >= (3, 0):
    from urllib.request import urlopen, urlretrieve
else:
    from urllib import urlopen, urlretrieve
import ssNake as sc
from ssNake import QtCore, QtWidgets


class UpdateWindow(QtWidgets.QWidget):
    """
    The window for updating ssNake.
    """

    def __init__(self, parent):
        """
        Initializes the update window.

        Parameters
        ----------
        parent : MainProgram
            The mainprogram object of ssNake.
        """
        super(UpdateWindow, self).__init__(parent)
        self.setWindowFlags(QtCore.Qt.Window | QtCore.Qt.Tool)
        self.father = parent
        try:
            req = urlopen('https://api.github.com/repos/smeerten/ssnake/tags')
            if sys.version_info >= (3, 0):
                info = json.loads(str(req.read().decode('utf-8')))
            else:
                info = json.loads(str(req.read()))
            req.close()
            self.nameList = [u'develop']
            self.urlList = [u'https://api.github.com/repos/smeerten/ssnake/zipball/develop']
            for i, _ in enumerate(info):
                self.nameList.append(info[i]['name'])
                self.urlList.append(info[i]['zipball_url'])
        except Exception:
            raise sc.SsnakeException("Could not connect to the server")
        self.setWindowTitle("Update ssNake")
        layout = QtWidgets.QGridLayout(self)
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid, 0, 0, 1, 2)
        message_label = QtWidgets.QLabel("<b>Note: Updating to versions after v1.5 needs to be done manually via the <a href='https://gitlab.science.ru.nl/mrrc/nmrzoo/ssnake'>new repository</a>.</b>")
        message_label.setWordWrap(True)
        message_label.setOpenExternalLinks(True)
        grid.addWidget(message_label, 0, 0)
        grid.addWidget(QtWidgets.QLabel("Update to version:"), 1, 0)
        self.versionDrop = QtWidgets.QComboBox(parent=self)
        self.versionDrop.addItems(self.nameList)
        self.versionDrop.setCurrentIndex(1)
        grid.addWidget(self.versionDrop, 2, 0)
        grid.addWidget(QtWidgets.QLabel("Current version: " + self.father.VERSION), 3, 0)
        cancelButton = QtWidgets.QPushButton("&Cancel")
        cancelButton.clicked.connect(self.closeEvent)
        okButton = QtWidgets.QPushButton("&Ok")
        okButton.clicked.connect(self.applyAndClose)
        box = QtWidgets.QDialogButtonBox()
        box.addButton(cancelButton, QtWidgets.QDialogButtonBox.RejectRole)
        box.addButton(okButton, QtWidgets.QDialogButtonBox.AcceptRole)
        layout.addWidget(box, 2, 0)
        layout.setColumnStretch(1, 1)
        self.show()

    def closeEvent(self, *args):
        """
        Closes the update window.
        """
        self.deleteLater()

    def applyAndClose(self, *args):
        """
        Asks the user to update and closes the window.
        """
        ssnake_location = os.path.dirname(os.path.dirname(__file__))
        try:
            os.mkdir(ssnake_location + os.path.sep + 'test')
            os.rmdir(ssnake_location + os.path.sep + 'test')
        except Exception:
            QtWidgets.QMessageBox.critical(QtWidgets.QWidget(), 'Update failed', 'You do not have write permission in the ssNake directory.', QtWidgets.QMessageBox.Ok)
            return
        message = "Update is going to replace all files in " + str(ssnake_location) + "\n ssNake needs to be restarted for the changes to take effect.\n Are you sure?"
        reply = QtWidgets.QMessageBox.question(self, 'Update', message, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
        if reply == QtWidgets.QMessageBox.Yes:
            self.deleteLater()
            progress = QtWidgets.QProgressDialog("Downloading...", "Cancel", 0, 3)
            progress.show()
            tempDir = tempfile.mkdtemp()
            filehandle, _ = urlretrieve(self.urlList[self.versionDrop.currentIndex()])
            progress.setValue(1)
            progress.setLabelText("Extracting...")
            zip_file = zipfile.ZipFile(filehandle, 'r')
            tmpPath = os.path.join(tempDir, zip_file.namelist()[0])
            zip_file.extractall(tempDir)
            progress.setValue(2)
            progress.setLabelText("Copying to destination...")
            shutil.rmtree(ssnake_location)
            shutil.move(tmpPath, ssnake_location)
            progress.setValue(2)
            progress.setLabelText("Cleaning up...")
            shutil.rmtree(tempDir)
            progress.close()
