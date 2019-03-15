![Alt text](src/Icons/logo.gif?raw=true)

ssNake
======

ssNake is a versatile tool for processing and analysing NMR spectra, focused on solid-state NMR.

Requirements
------------

ssNake requires:
- [python](http://python.org/download/) >= 2.7 or [python](http://python.org/download/) >= 3.4

And the following python packages are required[1]:
- [numpy](http://sourceforge.net/projects/numpy/files/NumPy/) >= 1.11.0
- [matplotlib](http://matplotlib.org/) >= 1.4.2
- [scipy](http://sourceforge.net/projects/scipy/files/scipy/) >= 0.14.1
- [PyQt4](http://www.riverbankcomputing.com/software/pyqt/download) >= 4.11.4
- [h5py](http://www.h5py.org/) >= 2.5.0 (for loading Matlab data)

On Ubuntu and Debian these packages can be installed using the package manager:
```
sudo apt-get install python python-numpy python-matplotlib python-scipy python-qt4 python-h5py
```

On Windows (and macOS) these packages can easily be installed by downloading [Anaconda](https://www.anaconda.com/distribution/).
During installation, you should enable the 'Add Anaconda to the system PATH environment variable' box.
This makes sure ssNake can find the python executable installed by Anaconda.

[1]: The program might work on older versions, but they have not been tested.

Installation
------------

###Linux and macOS###

To install ssNake, copy the ssNake directory to your favourite location (/usr/local/, for example).
ssNake can then be run by executing 'python /InstallPath/src/ssNake.py'.
Aliases or symlinks can be used to create a shortcut to start the program.
When multiple versions of Python are installed make sure that the correct one is used to start ssNake.
This can be done either by modifying the PATH environment variable or by starting ssNake with the full path to the Python binary, for example by using '/anaconda3/bin/python /InstallPath/src/ssNake.py'

###Windows###

On Windows systems, ssNake can be installed using the 'Windows installer' that we supply. This holds
a compiled version of both ssNake, and the relevant python libraries, and can be run without any
other requirements. For the installer, please visit our [website](https://www.ru.nl/science/solidstatenmr/software/ssnake/).

Users that have installed python via the Anaconda progrom descibed above can do the following:
To install ssNake, copy the ssNake directory to your favourite location (C:\Program Files\, for example).
ssNake can then be run by double clicking on the 'WindowsRun.bat' file.
Alternatively, you can execute the 'WindowsInstall.vbs' file from within the ssNake directory.
This creates shortcuts on your desktop and in the startmenu, which you can use to run ssNake.

Contributing
------------

1. Fork it
2. Create your feature branch
3. Commit your changes
4. Push to the branch
5. Submit a pull request

Creators
--------

**Bas van Meerten**

**Wouter Franssen**

Contact
-------
For question and suggestions mail to: ssnake@science.ru.nl

License
-------

This project is licensed under the GNU General Public License v3.0 - See [LICENSE.md](LICENSE.md) for details.
