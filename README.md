![Alt text](src/Icons/logo.gif?raw=true)

New repository location
=======================
From the release of version 1.5, official ssNake development will continue [here](https://gitlab.science.ru.nl/mrrc/nmrzoo/ssnake).

ssNake
======

ssNake is a versatile tool for processing and analysing NMR spectra, focused on solid-state NMR.

Requirements
------------

ssNake requires:
- [python](http://python.org/download/) >= 3.4

And the following python packages are required[1]:
- [numpy](http://sourceforge.net/projects/numpy/files/NumPy/) >= 1.11.0
- [matplotlib](http://matplotlib.org/) >= 1.5.0
- [scipy](http://sourceforge.net/projects/scipy/files/scipy/) >= 0.14.1
- [PyQt5](http://www.riverbankcomputing.com/software/pyqt/download) >= 5.15.0
- [h5py](http://www.h5py.org/) >= 2.5.0 (for loading Matlab data)

On Ubuntu and Debian these packages can be installed using the package manager:
```
sudo apt-get install python3 python3-numpy python3-matplotlib python3-scipy python3-pyqt5 python3-h5py
```

On Windows (and macOS) these packages can easily be installed by downloading [Anaconda](https://www.anaconda.com/distribution/).
During installation, you should enable the 'Add Anaconda to the system PATH environment variable' box.
This makes sure ssNake can find the python executable installed by Anaconda.

[1]: The program might work on older versions, but they have not been tested.

Installation
------------

### Linux and macOS ###

To install ssNake, copy the ssNake directory to your favourite location (/usr/local/, for example).
ssNake can then be run by executing 'python3 /InstallPath/src/ssNake.py'.
Aliases or symlinks can be used to create a shortcut to start the program.
When multiple versions of Python are installed make sure that the correct one is used to start ssNake.
This can be done either by modifying the PATH environment variable or by starting ssNake with the full path to the Python binary, for example by using '/anaconda3/bin/python /InstallPath/src/ssNake.py'

### Windows ###

On Windows systems, ssNake can be installed using the 'Windows installer' that we supply. This holds
a compiled version of both ssNake, and the relevant python libraries, and can be run without any
other requirements. The installer can be found on GitHub under 'Releases'.

Users that have installed python via the Anaconda progrom descibed above can do the following:
To install ssNake, copy the ssNake directory to your favourite location (C:\Program Files\, for example).
ssNake can then be run by double clicking on the 'WindowsRun.bat' file.
Alternatively, you can execute the 'WindowsInstall.vbs' file from within the ssNake directory.
This creates shortcuts on your desktop and in the startmenu, which you can use to run ssNake.

Contributing
------------

1. Fork it
2. Create your feature branch (preferably based on the develop branch)
3. Commit your changes
4. Push to the branch
5. Submit a pull request

Citations
---------

If you use the ssNake software please cite:

**ssNake: A cross-platform open-source NMR data processing and fitting application**  
S.G.J.van Meerten, W.M.J.Franssen, A.P.M.Kentgens  
Journal of Magnetic Resonance, Volume 301, 56-66 (2019)  
https://doi.org/10.1016/j.jmr.2019.02.006  

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
