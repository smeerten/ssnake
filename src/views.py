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

import copy
import numpy as np
from matplotlib.pyplot import get_cmap, colormaps
import matplotlib
import matplotlib.ticker as ticker
import spectrum as sc
from spectrumFrame import PlotFrame
import reimplement as reim

COLORMAPLIST = ['seismic','gray', 'BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'rainbow', 'jet']


COLORMAPLIST = [i for i in COLORMAPLIST if i in colormaps()]
COLORRANGELIST = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'coolwarm', 'Spectral', 'rainbow', 'jet']
COLORRANGELIST = [i for i in COLORRANGELIST if i in colormaps()]
COLORRANGELIST = ['none'] + COLORRANGELIST

COLORCYCLE = list(matplotlib.rcParams['axes.prop_cycle'])
COLORCONVERTER = matplotlib.colors.ColorConverter()


##################################################################################################
# the class from which the 1d data is displayed, the operations which only edit the content of this class are for previewing


class Current1D(PlotFrame):

    X_RESIZE = False
    Y_RESIZE = False
    NDIM_PLOT = 1  # Number of dimensions in the plot
    MARKER = ''
    LINESTYLE = '-'

    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        super(Current1D, self).__init__(root, fig, canvas)
        self.data = data  # the actual spectrum instance
        self.data1D = None  # the data1D
        if duplicateCurrent is None:
            self.axes = np.array([len(self.data.shape()) - 1], dtype=int)
            self.resetLocList()
            self.viewSettings = {"plotType": 0,
                                 "showTitle": self.root.father.defaultShowTitle,
                                 "axType": np.array([self.root.father.defaultUnits] * self.NDIM_PLOT, dtype=int),
                                 "ppm": np.array([self.root.father.defaultPPM] * self.NDIM_PLOT, dtype=bool),   # display frequency as ppm
                                 "color": self.root.father.defaultColor,
                                 "linewidth": self.root.father.defaultLinewidth,
                                 "minXTicks": self.root.father.defaultMinXTicks,
                                 "minYTicks": self.root.father.defaultMinYTicks,
                                 "grids": self.root.father.defaultGrids,
                                 "colorRange": self.root.father.defaultColorRange,
                                 "colorMap": self.root.father.defaultColorMap,
                                 "contourConst": self.root.father.defaultContourConst,
                                 "contourColors": [self.root.father.defaultPosColor, self.root.father.defaultNegColor],
                                 "pColorMap": self.root.father.defaultPColorMap,
                                 "diagonalBool": self.root.father.defaultDiagonalBool,
                                 "diagonalMult": self.root.father.defaultDiagonalMult,
                                 "contourSign": 0,
                                 "contourType": 0,
                                 "numLevels": 20,
                                 "minLevels": 0.1,
                                 "maxLevels": 1,
                                 "limitType": 0,
                                 "multiValue": 1.5,
                                 "projTop": 0,
                                 "projRight": 0,
                                 "projLimitsBool": False,
                                 "projLimits": [None, None, None, None],
                                 "projPos": [0, 0],
                                 # CurrentMulti variables
                                 "extraData": [],
                                 "extraLoc": [],
                                 "extraColor": [],
                                 "extraName": [],
                                 "extraAxes": [],
                                 "extraScale": [],
                                 "extraOffset": [],
                                 "extraShift": [],
                                 "extraShift2": [],
                                 # CurrentStacked variables
                                 "stackBegin": None,
                                 "stackEnd": None,
                                 "stackStep": None,
                                 "spacing": 0}
            self.upd()  # get the first slice of data
            if self.viewSettings['showTitle']:
                self.fig.suptitle(self.data.name)
            else:
                self.fig.suptitle('')
            self.resetSpacing(False)
            self.startUp()
        else:
            self.axes = self.fixAxes(copy.deepcopy(duplicateCurrent.axes))
            self.locList = copy.deepcopy(duplicateCurrent.locList)
            self.viewSettings = copy.deepcopy(duplicateCurrent.viewSettings)
            self.viewSettings.update({"extraData": [],
                                      "extraLoc": [],
                                      "extraColor": [],
                                      "extraName": [],
                                      "extraAxes": [],
                                      "extraScale": [],
                                      "extraOffset": [],
                                      "extraShift": [],
                                      "extraShift2": []})
            if len(self.viewSettings['axType']) != self.NDIM_PLOT:
                diff = self.NDIM_PLOT - len(self.viewSettings['axType'])
                self.viewSettings['axType'] = np.append(np.array([self.root.father.defaultUnits] * diff, dtype=int), self.viewSettings['axType'][-self.NDIM_PLOT:])
                self.viewSettings['ppm'] = np.append(np.array([self.root.father.defaultPPM] * diff, dtype=bool), self.viewSettings['ppm'][-self.NDIM_PLOT:])
            if type(self) not in (CurrentStacked, CurrentArrayed):
                self.viewSettings.update({"stackBegin": None, "stackEnd": None, "stackStep": None})
            self.xminlim = duplicateCurrent.xminlim
            self.xmaxlim = duplicateCurrent.xmaxlim
            self.yminlim = duplicateCurrent.yminlim
            self.ymaxlim = duplicateCurrent.ymaxlim
            xReset = self.X_RESIZE or duplicateCurrent.X_RESIZE
            yReset = self.Y_RESIZE or duplicateCurrent.Y_RESIZE
            self.upd()  # get the first slice of data
            if (self.viewSettings["stackStep"] is None) and (type(self) in (CurrentStacked, CurrentArrayed)):
                if self.data1D.ndim() > 1 and self.data1D.shape()[0] > 100:
                    self.viewSettings["stackStep"] = self.data1D.shape()[0] // 100 + 1
                    self.upd()
            if self.viewSettings['showTitle']:
                self.fig.suptitle(self.data.name)
            else:
                self.fig.suptitle('')
            if type(self) != type(duplicateCurrent):
                self.resetSpacing(False)
            self.startUp(xReset, yReset)

    def resetSpacing(self, *args):
        # Dummy function
        pass

    def shape(self):
        """
        Returns the shape of the data

        Returns
        -------
        tuple:
            The shape along each dimension
        """
        return self.data1D.shape()

    def len(self, dim=-1):
        """
        Returns the length of a specific dimension of the data. By default
        the last dimension (-1) is selected.

        Parameters
        ----------
        dim (optional = -1): int
            The dimension of which the length should be returned.

        Returns
        -------
        int:
            The length
        """
        return self.shape()[dim]

    def ndim(self):
        """
        Returns the number of dimensions of the data.

        Returns
        -------
        int:
            The number of dimension
        """
        return self.data1D.ndim()

    def freq(self, dim=-1):
        """
        Returns the spectrum frequency of a specified dimension
        (i.e. centre frequency) in hertz.
        By default, the last dimension (-1) is selected.

        Parameters
        ----------
        dim (optional = -1): int
            The dimension of which the frequency should be returned.

        Returns
        -------
        float:
            The frequency
        """
        return self.data1D.freq[dim]

    def ref(self, dim=-1):
        """
        Returns the reference frequency in Hz for a selected dimension.
        By default, the last dimension (-1) is selected.
        If no reference is defined for the dimension, the centre frequency
        is retuned instead.

        Parameters
        ----------
        dim (optional = -1): int
            The dimension of which the reference frequency should be returned.

        Returns
        -------
        float:
            The reference frequency of the dimension
        """
        if self.data1D.ref[dim] is not None:
            return self.data1D.ref[dim]
        return self.data1D.freq[dim]

    def sw(self, dim=-1):
        """
        Returns the spectral width  (sweep width) in Hz for a selected dimension.
        By default, the last dimension (-1) is selected.

        Parameters
        ----------
        dim (optional = -1): int
            The dimension of which the sw should be returned.

        Returns
        -------
        float:
            The spectral width of the dimension
        """
        return self.data1D.sw[dim]

    def spec(self, dim=-1):
        """
        Returns the state of a specified dimension (spectrum or FID).
        By default, the last dimension (-1) is selected.

        Parameters
        ----------
        dim (optional = -1): int
            The dimension of which the spec/fid state should be returned.

        Returns
        -------
        bool:
            True for spectrum, False for FID.
        """
        return self.data1D.spec[dim]

    def xax(self, dim=-1):
        """
        Returns the x-axis for the selected dimension.
        By default, the last dimension (-1) is selected.

        Parameters
        ----------
        dim (optional = -1): int
            The dimension of which the x-axis state should be returned.

        Returns
        -------
        ndarray:
            1D array with the x-axis.
        """
        return self.data1D.xaxArray[dim]

    def wholeEcho(self, dim=-1):
        """
        Returns the whole Echo boolean for the selected dimension.
        By default, the last dimension (-1) is selected.

        Parameters
        ----------
        dim (optional = -1): int
            The dimension of which the whole echo state should be returned.

        Returns
        -------
        bool:
            The whole echo state of the dimension.
        """
        return self.data1D.wholeEcho[dim]

    def getCurrentAxMult(self, axis=-1):
        """
        Returns the axis multiplier (i.e. unit) for the selected dimension.
        By default, the last dimension (-1) is selected.

        Parameters
        ----------
        axis (optional = -1): int
            The dimension of which the multiplier should be returned.

        Returns
        -------
        float:
            The multiplier for the dimension.
        """
        return self.getAxMult(self.spec(axis), self.getAxType(), self.getppm(), self.freq(axis), self.ref(axis))

    def getAxType(self, num=-1):
        """
        Returns the axis type (i.e. unit) for the selected dimension.
        By default, the last dimension (-1) is selected.

        The type can be 0,1,2 or 3.
        For a spectrum axis:
            0: Hz
            1: kHz
            2: MHz
            3: ppm

        For an FID axis:
            0: s
            1: ms
            2: us

        Parameters
        ----------
        num (optional = -1): int
            The dimension of which the type should be returned.

        Returns
        -------
        int:
            The type index for the dimension.
        """
        return self.viewSettings["axType"][num]

    def getppm(self, num=-1):
        """
        Get the ppm boolean of the selected dimension.
        By default, the last dimension (-1) is selected.
        Is always False for FID axes.

        Parameters
        ----------
        num (optional = -1): int
            The dimension of which the ppm boolean should be returned.

        Returns
        -------
        bool:
            True if the axis is in ppm, False if not.
        """
        if self.data1D.ref[num] == 0.0 or self.data1D.freq[num] == 0.0:
            return False
        return self.viewSettings["ppm"][num]

    def getDataType(self, data):
        """
        Returns the input data in the current representation
        (real, imag, both, abs).

        Parameters
        ----------
        data: ndarray
            The data that should be converted to the representation

        Returns
        -------
        ndarray:
            The changed data
        """

        typeList = [np.real, np.imag, np.array, np.abs]
        return typeList[self.viewSettings["plotType"]](data)

    def startUp(self, xReset=True, yReset=True):
        """
        Plot the data, with optional resets of the x an y axis limits.

        Parameters
        ----------
        xReset (optional = True): bool
            Whether to reset the x-axis limits
        yReset (optional = True): bool
            Whether to reset the y-axis limits
        """
        self.showFid()  # plot the data
        self.plotReset(xReset, yReset)  # reset the axes limits

    def fixAxes(self, axes):
        """
        Fixes the axes.

        Parameters
        ----------
        axes: list
            List of the axes indexes that should be changed.

        Returns
        -------
        list:
            The new axis list
        """
        if len(axes) < self.NDIM_PLOT:
            fullAxes = np.arange(self.data.ndim())
            fullAxes = np.delete(fullAxes, axes)
            diff = self.NDIM_PLOT - len(axes)
            return np.append(fullAxes[-diff:], axes[-self.NDIM_PLOT:])
        elif len(axes) > self.NDIM_PLOT:
            return axes[-self.NDIM_PLOT:]
        return axes

    def rename(self, name):
        """
        Rename the data set.

        Parameters
        ----------
        name: str
            The new name
        """
        self.data.rename(name)
        if self.viewSettings['showTitle']:
            self.fig.suptitle(self.data.name)
        else:
            self.fig.suptitle('')
        self.canvas.draw()

    def copyCurrent(self, root, fig, canvas, data):
        """
        Make a copy of the current data structure and the
        associated canvas and axis information.

        Parameters
        ----------
        root: main1d window class
        fig: matplotlib figure
        canvas: matplotlib figure canvas
        data: spectrum class instance
        """
        return self.__class__(root, fig, canvas, data, self)

    def upd(self):
        """
        Get new data from the data instance.

        Returns
        -------
        bool:
            True if update was OK
        """
        if self.data.ndim() <= self.axes[-1]:
            self.axes = np.array([len(self.data.shape()) - 1])
        if len(self.locList) != self.data.ndim():
            self.resetLocList()
        try:
            self.data1D = self.data.getSlice(self.axes, self.locList)
            if self.data1D is None:
                self.root.rescue()
        except Exception:
            self.resetLocList()
            self.data1D = self.data.getSlice(self.axes, self.locList)
        if self.ref() == 0.0 or self.freq() == 0.0: #Check if ppm is allowed
            self.viewSettings["ppm"][-1] = False
        return True

    def setSlice(self, axes, locList):
        """
        Change the slice

        Parameters
        ----------
        axes: ndarray
            List of the axes along which the data is plot
            (contains a single entry for this 1D plot)
        locList: list
            List with the location (i.e. data positions) of each dimension
            of the data.
        """
        axesSame = True
        if not np.array_equal(self.axes, axes):
            axesSame = False
            self.axes = axes
        self.locList = locList
        self.upd()
        self.showFid()
        if not axesSame:
            self.plotReset()

    def resetLocList(self):
        """
        Resets the locations (i.e. data positions) of all dimensions.
        Default values are [0] for all dimensions.
        """
        self.locList = [0] * self.data.ndim()

    def getSelect(self):
        """
        Get thew slice operator to construct the 1D data from the total
        data.

        Returns
        -------
        ndarray:
            Object array with the slices for each dimension.
        """
        tmp = np.array(self.locList, dtype=object)
        tmp[self.axes] = slice(None)
        return tmp

    def setGrids(self, grids):
        """
        Turn the plot grid along x and y axis on or off.

        Parameters
        ----------
        grids: list of booleans
            [xgrid,ygrid]
        """
        self.viewSettings["grids"] = grids

    def setDiagonal(self, diagonalBool=None, diagonalMult=None):
        """
        Change the plotting of a diagonal line in a 2D plot.

        Parameters
        ----------
        diagonalBool (optional = None): bool
            True if diagonal should be plot
        diagonalMult (optional = None): float
            Multiplier of the diagonal.
        """
        if diagonalBool is not None:
            self.viewSettings["diagonalBool"] = diagonalBool
        if diagonalMult is not None:
            self.viewSettings["diagonalMult"] = diagonalMult
        self.showFid()

    def reload(self, *args):
        """
        Reloads the data, updates and resets the plot.
        """
        self.data.reload()
        self.upd()
        self.showFid()
        self.plotReset()
        self.root.addMacro(['reload'])

    def isComplex(self, *args):
        """
        Checks if an axis is complex. The last axis (direct dimension) is always complex.
        """
        return self.data.isComplex(*args)

    def setNoUndo(self, val):
        """
        Sets the 'no undo' flag. When this flag is enabled
        no undo information is saved.

        Parameters
        ----------
        val: bool
           The no undo flag
        """
        self.root.addMacro(['setNoUndo', (val,)])
        self.data.setNoUndo(val)

    def real(self, *args):
        """
        Takes the real value along the current axis.
        """
        self.root.addMacro(['real', (self.axes[-1] - self.data.ndim(), )])
        self.data.real(self.axes[-1])
        self.upd()
        self.showFid()

    def imag(self, *args):
        """
        Takes the imaginary value along the current axis.
        """
        self.root.addMacro(['imag', (self.axes[-1] - self.data.ndim(), )])
        self.data.imag(self.axes[-1])
        self.upd()
        self.showFid()

    def abs(self, *args):
        """
        Takes the absolute value along the current axis.
        """
        self.root.addMacro(['abs', (self.axes[-1] - self.data.ndim(), )])
        self.data.abs(self.axes[-1])
        self.upd()
        self.showFid()

    def conj(self, *args):
        """
        Takes the complex conjugate value along the current axis.
        """
        self.root.addMacro(['conj', (self.axes[-1] - self.data.ndim(), )])
        self.data.conj(self.axes[-1])
        self.upd()
        self.showFid()

    def setPhaseInter(self, phase0in, phase1in, phase2in):
        """
        Interactive changing the phase without editing the actual data.

        Parameters
        ----------
        phase0in: float
            The 0th order phase
        phase1in: float
            The 1st order phase
        """
        phase0 = float(phase0in)
        phase1 = float(phase1in)
        phase2 = float(phase2in)
        self.data1D.phase(phase0, phase1, phase2, -1)
        self.showFid()
        self.upd()

    def applyPhase(self, phase0, phase1, phase2=0, select=False):
        """
        Phase the data.

        Parameters
        ----------
        phase0in: float
            The 0th order phase
        phase1in: float
            The 1st order phase
        select (optional = False): bool
            If true, apply only for current slice
        """
        phase0 = float(phase0)
        phase1 = float(phase1)
        phase2 = float(phase2)
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['phase', (phase0, phase1, phase2, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.phase(phase0, phase1, phase2, self.axes[-1], selectSlice)
        self.upd()
        self.showFid()

    def correctDFilter(self):
        """
        Corrects the digital filter via first order phasing.
        The filter value (if any) is set upon loading the data.
        """
        self.root.addMacro(['correctDFilter', (self.axes[-1] - self.data.ndim(),)])
        self.data.correctDFilter(self.axes[-1])
        self.upd()
        self.showFid()

    def complexFourier(self):
        """
        Complex Fourier transform along the current axis.
        Redraw the plot.
        """
        self.root.addMacro(['complexFourier', (self.axes[-1] - self.data.ndim(), )])
        self.data.complexFourier(self.axes[-1])
        self.upd()
        if isinstance(self, (CurrentStacked, CurrentArrayed)):
            self.resetSpacing()
        self.showFid()
        self.plotReset()

    def realFourier(self):
        """
        Real Fourier transform along the current axis.
        Redraw the plot.
        """
        self.root.addMacro(['realFourier', (self.axes[-1] - self.data.ndim(), )])
        self.data.realFourier(self.axes[-1])
        self.upd()
        if isinstance(self, (CurrentStacked, CurrentArrayed)):
            self.resetSpacing()
        self.showFid()
        self.plotReset()

    def fftshift(self, inv=False):  # fftshift the actual data and replot
        """
        fftshift along the current axis.
        Redraw the plot.

        Parameters
        ----------
        inv: bool
            True if inverse shift, False if normal.
        """
        self.root.addMacro(['fftshift', (self.axes[-1] - self.data.ndim(), inv)])
        self.data.fftshift(self.axes[-1], inv)
        self.upd()
        self.showFid()

    def apodPreview(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxis=None):  
        """
        Preview the effect of apodization without touching the data structure.

        Parameters
        ----------
        lor (optional = None): float
            Lorentzian broadening in Hz
        gauss (optional = None): float
            Gaussian broadening in Hz
        cos2 (optional = [None, None]): list
            List of floats. The first value is the frequency (two times the number of periods in the time domain).
            The second value is the phase shift in degrees.
        hamming (optional = None): float
            Hamming apodization frequency
        shift (optional = 0.0): float
            The amount of time shift to the function (positive is to the right).
        shifting (optional = 0.0): float
            A shift in time of the function as a function of the x-axis values along shiftingAxis.
        shiftingAxis (optional = None): int
            The dimension for the shifting.
        """
        if cos2 is None:
            cos2 = [None, None]
        y = self.data1D.data.copy()
        preview = True
        if shiftingAxis is None:
            curve = self.data1D.apodize(lor, gauss, cos2, hamming, shift, shifting, None, -1, preview=preview)
        else:
            if (self.ndim() > 1) and (shiftingAxis == self.axes[-2]):
                curve = self.data1D.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxis, -1, preview=preview)
            else:
                if shiftingAxis == self.axes[-1]:
                    raise sc.SpectrumException('shiftingAxis cannot be equal to axis')
                shift += shifting * self.locList[shiftingAxis] / self.data.sw[shiftingAxis]
                curve = self.data1D.apodize(lor, gauss, cos2, hamming, shift, 0.0, None, -1, preview=preview)
        if self.spec() == 0 and not isinstance(self, CurrentContour):
            tmp = self.getDataType(y.getHyperData(0))
            scale = np.max([np.real(tmp), np.imag(tmp)])
            self.showFid(y, curve[0], scale*np.array(curve[1]), extraColor=['g'])
        else:
            self.showFid(y)
        self.upd()
    
    def applyApod(self, lor=None, gauss=None, cos2=None, hamming=None, shift=0.0, shifting=0.0, shiftingAxis=0, select=False):  # apply the apodization to the actual data
        """
        Apply apodization effects.

        Parameters
        ----------
        lor (optional = None): float
            Lorentzian broadening in Hz
        gauss (optional = None): float
            Gaussian broadening in Hz
        cos2 (optional = [None, None]): list
            List of floats. The first value is the frequency (two times the number of periods in the time domain).
            The second value is the phase shift in degrees.
        hamming (optional = None): float
            Hamming apodization frequency
        shift (optional = 0.0): float
            The amount of time shift to the function (positive is to the right).
        shifting (optional = 0.0): float
            A shift in time of the function as a function of the x-axis values along shiftingAxis.
        shiftingAxis (optional = None): int
            The dimension for the shifting.
        select (optional = False): boolean
            If True, apply only to the current slice.
        """
        if cos2 is None:
            cos2 = [None, None]
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['apodize', (lor, gauss, cos2, hamming, shift, shifting, shiftingAxis, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.apodize(lor, gauss, cos2, hamming, shift, shifting, shiftingAxis, self.axes[-1], selectSlice)
        self.upd()
        self.showFid()

    def setFreq(self, freq, sw):  # set the frequency of the actual data
        """
        Set the center frequency and spectral width (sweep width).

        Parameters
        ----------
        freq: float
            Center frequency in Hz
        sw: float
            Spectral width in Hz

        """
        self.root.addMacro(['setFreq', (freq, sw, self.axes[-1] - self.data.ndim())])
        self.data.setFreq(freq, sw, self.axes[-1])
        self.upd()
        self.showFid()

    def scaleSw(self, scale):
        """
        Scale the spectral width (sweep width).

        Parameters
        ----------
        scale: float
            The multiplier
        """
        self.root.addMacro(['scaleSw', (scale, self.axes[-1] - self.data.ndim())])
        self.data.scaleSw(scale, self.axes[-1])
        self.upd()
        self.showFid()

    def setRef(self, ref):
        """
        Set the reference frequency (0 Hz position).

        Parameters
        ----------
        ref: float
            The reference frequency
        """
        oldref = self.ref()
        self.root.addMacro(['setRef', (ref, self.axes[-1] - self.data.ndim())])
        self.data.setRef(ref, self.axes[-1])
        if ref is None:
            ref = self.freq()
        self.upd()
        val = self.getAxType()
        if self.spec() == 1:
            if self.getppm():
                self.xminlim = (self.xminlim * oldref * 10**-6 + oldref - ref) / (ref * 10**-6)
                self.xmaxlim = (self.xmaxlim * oldref * 10**-6 + oldref - ref) / (ref * 10**-6)
            else:
                self.xminlim = self.xminlim + (oldref - ref) / 10**(val * 3)  # set new limits, and scale for axis type
                self.xmaxlim = self.xmaxlim + (oldref - ref) / 10**(val * 3)
        self.showFid()

    def regrid(self, limits, numPoints):
        """
        Regrind along the current dimension. This creates a new x-axis, and interpolates to
        construct the new y-values.

        Parameters
        ----------
        limits: list
            List with the minimum and maximum position of the new x-axis (in Hz).
        numPoints: int
            Number of points in the new x-axis
        """
        self.root.addMacro(['regrid', (limits, numPoints, self.axes[-1] - self.data.ndim())])
        self.data.regrid(limits, numPoints, self.axes[-1])
        self.upd()
        self.showFid()
        self.plotReset()

    def SN(self, minNoise, maxNoise, minPeak, maxPeak):
        """
        Get the signal-to-noise ratio of a specific region.

        Parameters
        ----------
        minNoise: int
            Minimum position of the noise region.
        maxNoise: int
            Maximum position of the noise region.
        minPeak: int
            Minimum position of the region with the peak inside.
        maxPeak: int
            Maximum position of the region with the peak inside.

        Returns
        -------
        float:
            The SNR
        """
        minN = min(minNoise, maxNoise)
        maxN = max(minNoise, maxNoise)
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpData = np.real(self.getDataType(tmpData))
        return np.max(tmpData[minP:maxP]) / (np.std(tmpData[minN:maxN]))

    def fwhm(self, minPeak, maxPeak, level, unitType=None):
        """
        Get the Full Width at Half Maximum (FWHM) for a peak in
        a specific region. The level can be pout in by the user,
        so a level of 0.5 will give the FWHM, while other values
        report the width at other heights.

        Parameters
        ----------
        minPeak: int
            Minimum position of the region with the peak inside.
        maxPeak: int
            Maximum position of the region with the peak inside.
        level: float
            The level for which the width is to be calculated. For
            the half width, this should be 0.5.
        unitType (optional = None): int
            The unit type, only needs to be supplied when the type should
            be different than in the current plot.
            0: Hz
            1: kHz
            2: MHz
            3: ppm

        Returns
        -------
        float:
            The width
        """

        from scipy.interpolate import UnivariateSpline
        if unitType is None:
            axType = self.getAxType()
            ppm = self.getppm()
        else:
            axType = unitType
            if unitType == 3:  # ppm
                ppm = 1
            else:
                ppm = 0
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpData = np.real(self.getDataType(tmpData))
        maxPos = np.argmax(tmpData[minP:maxP])
        x = self.xax() * self.getAxMult(self.spec(), axType, ppm, self.freq(), self.ref())
        maxX = x[minP:maxP][maxPos]
        spline = UnivariateSpline(x, tmpData - level * tmpData[minP:maxP][maxPos], s=0)
        zeroPos = spline.roots()
        left = zeroPos[zeroPos > maxX]
        right = zeroPos[zeroPos < maxX]
        if right.size > 0 and left.size > 0:
            return abs(left[0] - right[-1])
        return 0.0

    def COM(self, minPeak, maxPeak, unitType=None):  # Centre of Mass
        """
        Get the centre of mass (COM) for a peak in
        a specific region.

        Parameters
        ----------
        minPeak: int
            Minimum position of the region with the peak inside.
        maxPeak: int
            Maximum position of the region with the peak inside.
        unitType (optional = None): int
            The unit type, only needs to be supplied when the type should
            be different than in the current plot.
            0: Hz
            1: kHz
            2: MHz
            3: ppm

        Returns
        -------
        float:
            The centre of mass
        """
        dim = len(minPeak)
        ppm = []
        axType = []
        if unitType is None:
            for i in range(dim):
                axType.append(self.getAxType(-(i + 1)))
                ppm.append(self.getppm(-(i+1)))
        else:
            axType = unitType
            for i in range(dim):
                if unitType[i] == 3:
                    ppm.append(1)
                else:
                    ppm.append(0)
        for i in range(dim):
            if minPeak[i] > maxPeak[i]:
                minPeak[i], maxPeak[i] = maxPeak[i], minPeak[i]
        tmpData = self.data1D.getHyperData(0)
        tmpData = np.real(self.getDataType(tmpData))
        #Slice data
        slc = ()
        for i in range(dim):
            slc = slc + (slice(minPeak[-(i+1)], maxPeak[-(i+1)], None),)
        tmpData = tmpData[slc]
        COM = []
        axes = list(range(dim))
        for i in range(dim):
            sumAxes = list(axes)
            del sumAxes[-(i+1)]
            tmpAxis = self.xax(-(i+1))
            tmpAxis = tmpAxis[minPeak[i]:maxPeak[i]] * self.getAxMult(self.spec(-(i+1)), axType[i], ppm[i], self.freq(-(i+1)), self.ref(-(i+1)))
            tmpDataNew = np.sum(tmpData, tuple(sumAxes))
            # COM = 1/M *sum(m_i * r_i)
            COM.append(1.0 / np.sum(tmpDataNew) * np.sum(tmpDataNew * tmpAxis))
        return COM

    def Integrals(self, minPeak, maxPeak):
        """
        Get the integral for a region.

        Parameters
        ----------
        minPeak: int
            Minimum position of the region.
        maxPeak: int
            Maximum position of the region.

        Returns
        -------
        float:
            Integral
        ndarray:
            The x-axis of the selected range
        ndarray:
            The cumulative sum of the y-values of the range
        float:
            The absolute maximum of y-values
        """
        dim = len(minPeak)
        for i in range(dim): #Check max > min
            if minPeak[i] > maxPeak[i]:
                minPeak[i], maxPeak[i] = maxPeak[i], minPeak[i]
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-dim) + (slice(None), ) * dim]
        tmpAxis = self.xax()
        totShape = tmpData.shape
        tmpData = np.real(self.getDataType(tmpData))
        maxim = np.max(np.abs(tmpData))
        tmpAxis = tmpAxis[minPeak[0]:maxPeak[0]+1]
        slc = tuple()
        for i in reversed(range(dim)): #Make slice operator along all dimensions
            slc = slc + (slice(minPeak[i],maxPeak[i] + 1),)
        tmpData = tmpData[slc] #slice data
        if self.spec() == 0 and dim == 1:
            intSum = np.cumsum(tmpData)
        elif self.spec() == 1 and dim == 1:
            intSum = np.cumsum(tmpData[-1::-1])[-1::-1]
        else:
            intSum = None
        inte = np.sum(tmpData)
        for i in range(dim): #Scale sum for each integrated dimension
            if self.spec(-i - 1) == 0:
                inte /= self.sw(-i - 1)
            else:
                inte *= self.sw(-i - 1) / (1.0 * totShape[-i - 1])
        return inte, tmpAxis, intSum, maxim

    def MaxMin(self, minPeak, maxPeak, type='max'):
        """
        Get the maximum or minimum value of a region.

        Parameters
        ----------
        minPeak: int
            Minimum position of the region.
        maxPeak: int
            Maximum position of the region.
        type (optional = 'max'): str
            'max' or 'min'

        Returns
        -------
        float:
            The min/max value of the region
        """
        minP = min(minPeak, maxPeak)
        maxP = max(minPeak, maxPeak)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpData = np.real(self.getDataType(tmpData))
        if type == 'max':
            return np.max(tmpData[minP:maxP])
        if type == 'min':
            return np.min(tmpData[minP:maxP])

    def integralsPreview(self, x, y, maxim):
        """
        Plot preview lines for the integrals tool.

        Parameters
        ----------
        x: list of ndarrays
            x-values of the extra lines
        y: list of ndarrays
            y-values of the extra lines
        maxim: float
            Maximum values the each integrated regions

        """
        xNew = []
        yNew = []
        scale = 0
        for num, _ in enumerate(x):
            if x[num] is not None and y[num] is not None:
                xNew.append(x[num])
                yNew.append(y[num])
                scale = np.max([scale, abs(yNew[-1][0]), abs(yNew[-1][-1])])
        for num, _ in enumerate(yNew):
            yNew[num] = yNew[num] / scale * maxim
        self.showFid(extraX=xNew, extraY=yNew, extraColor=['g'])

    def resizePreview(self, size, pos):
        """
        Preview a resizing (zero filling) operation on the data.

        Parameters
        ----------
        size: int
            The new number of data points
        pos: int
            The position where any data points should be inserted
            or removed.
        """
        self.data1D.resize(size, pos, -1)
        self.showFid()
        if not self.spec():
            self.plotReset(True, False)
        self.upd()

    def resize(self, size, pos):  # set size to the actual data
        """
        Apply a resizing (zero filling) operation on the data.

        Parameters
        ----------
        size: int
            The new number of data points
        pos: int
            The position where any data points should be inserted
            or removed.
        """
        self.root.addMacro(['resize', (size, pos, self.axes[-1] - self.data.ndim())])
        self.data.resize(size, pos, self.axes[-1])
        self.upd()
        self.showFid()
        if not self.spec():
            self.plotReset(True, False)

    def lpsvd(self, nPredict, maxFreq, forward, numPoints):
        """
        Apply linear prediction on the data. Both forward and backwards predictions
        are supported.

        Parameters
        ----------
        nPredict : int
            The number of datapoints to predict.
        maxFreq : int
            The maximum number of frequencies to take from the SVD.
        forward : bool
            If True, a forward prediction is performed, otherwise a backward prediction is performed.
        numPoints : int, optional
            The number of points to use for SVD.
        """
        self.root.addMacro(['lpsvd', (nPredict, maxFreq, forward, numPoints, self.axes[-1] - self.data.ndim())])
        self.data.lpsvd(nPredict, maxFreq, forward, numPoints, self.axes[-1])
        self.upd()
        self.showFid()
        if not self.spec():
            self.plotReset(True, False)

    def setSpec(self, val):
        """
        Change the time/spectrum domain type of the current dimension.

        Parameters
        ----------
        val: bool
            True: spectrum, False: FID
        """
        self.root.addMacro(['setSpec', (val, self.axes[-1] - self.data.ndim())])
        self.data.setSpec(val, self.axes[-1])
        self.upd()
        if isinstance(self, CurrentArrayed):
            self.resetSpacing()
        self.showFid()
        self.plotReset()

    def swapEcho(self, idx):
        """
        Apply a swap echo operation.

        Parameters
        ----------
        idx: int
            Position where the swap should occur.
        """
        self.root.addMacro(['swapEcho', (idx, self.axes[-1] - self.data.ndim())])
        self.data.swapEcho(idx, self.axes[-1])
        self.upd()
        self.showFid()

    def swapEchoPreview(self, idx):
        """
        Preview a swap echo operation.

        Parameters
        ----------
        idx: int
            Position where the swap should occur.
        """
        self.data1D.swapEcho(idx, -1)
        self.showFid()
        self.upd()

    def setWholeEcho(self, value):
        """
        Set the Whole Echo toggle for the current dimension.

        Parameters
        ----------
        value: bool
            The new Whole Echo value
        """
        valBool = value != 0
        self.root.addMacro(['setWholeEcho', (valBool, self.axes[-1] - self.data.ndim())])
        self.data.setWholeEcho(valBool, self.axes[-1])
        self.upd()

    def shift(self, shift, select=False):
        """
        Shifts the data along the current dimension.
        Shifting is always done in the time domain (i.e. is the current data is a spectrum
        it is Fourier transformed back and forward to do the shift).

        Parameters
        ----------
        shift: int
            The amount of data points to shift (negative is left shift, positive right shift)
        select (optional = False): boolean
            If True, apply only to the current slice.
        """
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['shift', (shift, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.shift(shift, self.axes[-1], selectSlice)
        self.upd()
        self.showFid()

    def shiftPreview(self, shift):
        """
        Shows a preview of the shift data operation.

        Parameters
        ----------
        shift: int
            The amount of data points to shift (negative is left shift, positive right shift)
        """
        self.data1D.shift(shift, -1)
        self.showFid()
        self.upd()

    def roll(self, shift, select=False, shift_axis=True):
        """
        Circularly rolls the data along the current dimension. Non-integer shift values are allowed.

        Parameters
        ----------
        shift: float
            The amount of data points to roll (negative is left roll, positive right roll)
        select (optional = False): boolean
            If True, apply only to the current slice.
        shift_axis: boolean 
            If True then shifts the axis when rolling data. Applies on frequency domain only and if select is False only.
        """
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['roll', (shift, self.axes[-1] - self.data.ndim(), selectSlice, shift_axis)])
        self.data.roll(shift, self.axes[-1], selectSlice, shift_axis)
        self.upd()
        self.showFid()

    def rollPreview(self, shift, shift_axis=False):
        """
        Shows a preview of the roll data operation.

        Parameters
        ----------
        shift: float
            The amount of data points to roll (negative is left roll, positive right roll)
        shift_axis: boolean 
            If True then shifts the axis when rolling data. Applies on frequency domain only and if select is False only.
        """
        self.data1D.roll(shift, -1, shift_axis=shift_axis)
        self.showFid()
        self.upd()

    def align(self, pos1, pos2):
        """
        Aligns the maximum of each slice along this dimension, within the pos1-pos2 region.


        Parameters
        ----------
        pos1: int
            First data position.
        pos2: int
            Second data position.
        """
        self.root.addMacro(['align', (pos1, pos2, self.axes[-1] - self.data.ndim())])
        self.data.align(pos1, pos2, self.axes[-1])
        self.upd()
        self.showFid()

    def getdcOffset(self, pos1, pos2):
        """
        Gets the average (i.e. DC offset) of a selected region.

        Parameters
        ----------
        pos1: int
            First data position.
        pos2: int
            Second data position.

        Returns
        -------
        complex value:
            The offset value
        """
        minPos = int(min(pos1, pos2))
        maxPos = int(max(pos1, pos2))
        if minPos != maxPos:
            tmpData = self.data1D.data[(len(self.shape()) - 1) * (slice(None), ) + (slice(minPos, maxPos), )]
            return np.mean(tmpData.getHyperData(0))
        return 0

    def dcOffset(self, offset):
        """
        Corrects the DC offset (i.e. subtracts the offset value from all data points).

        Parameters
        ----------
        offset: complex value
            The amount of offset that is to be subtracted.
        """
        self.data1D.subtract([offset])
        self.showFid()
        self.upd()

    def baselineFunctionFit(self, x, data, bArray, degree, type):
        """
        Fit a function through the selected data.

        Parameters
        ----------
        x: ndarray
            The x-axis
        data: ndarray
            The data along this slice
        bArray: ndarray, boolian
            The points that should be used
        degree: int
            Number of polynomial orders
        type: str
            Either 'poly' for polynomial fit, or 'sin/cos' for sine/cosine

        Returns
        -------
        ndarray:
            The fitted polynomial
        """
        if type == 'poly':
            import numpy.polynomial.polynomial as poly
            polyCoeff = poly.polyfit(x[bArray], data[bArray], degree)
            return poly.polyval(x, polyCoeff)
        elif type == 'sin/cos':
            import numpy.linalg as ln
            fit = np.ones([len(x),degree *2 + 1])
            x_fake = np.linspace(0, 2 * np.pi, len(x)) #cos/sine always in 0-->2pi regime
            for order in range(degree):
                fit[:,order * 2 + 1] = np.cos((order + 1) * x_fake)
                fit[:,order * 2 + 2] = np.sin((order + 1) * x_fake)
            values = ln.lstsq(fit[bArray], data[bArray] ,rcond = None)[0] #Perform leats squares fit
            return np.sum(fit * values, axis = 1)

    def baselineCorrectionAll(self, degree, removeList, type, invert=False):
        """
        Correct baseline of a series of data

        Parameters
        ----------
        degree: int
            Polynomial degree
        removeList: list
            Indexes of points not include in function fit
        type: str
            Either 'poly' for polynomial fit, or 'sin/cos' for sine/cosine
        invert (optional = False): boolean
            If True, the removeList is treated as an include list (i.e. inverting the selection)           
        """
        tmpAx = np.arange(self.len())
        bArray = np.array([True] * self.len())
        for i in range(int(np.floor(len(removeList) / 2.0))):
            minVal = min(removeList[2 * i], removeList[2 * i + 1])
            maxVal = max(removeList[2 * i], removeList[2 * i + 1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        if invert:
            bArray = np.logical_not(bArray)
        y = np.apply_along_axis(lambda data: self.baselineFunctionFit(self.xax(), data, bArray, degree, type), self.axes[-1], self.data.getHyperData(0))
        y = np.real(self.getDataType(y))
        self.root.addMacro(['subtract', (y,)])
        self.data.subtract(y)

    def baselineCorrection(self, degree, removeList, type, select=False, invert=False):
        """
        Correct baseline of spectrum/fid.

        Parameters
        ----------
        degree: int
            Polynomial degree
        removeList: list
            Indexes of points not include in function fit
        type: str
            Either 'poly' for polynomial fit, or 'sin/cos' for sine/cosine
        select (optional = False): boolean
            If True, apply only to the current slice.
        invert (optional = False): boolean
            If True, the removeList is treated as an include list (i.e. inverting the selection)           
        """
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpAx = np.arange(self.len())
        bArray = np.array([True] * self.len())
        for i in range(int(np.floor(len(removeList) / 2.0))):
            minVal = min(removeList[2 * i], removeList[2 * i + 1])
            maxVal = max(removeList[2 * i], removeList[2 * i + 1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        if invert:
            bArray = np.logical_not(bArray)
        y = self.baselineFunctionFit(self.xax(), tmpData, bArray, degree, type)
        y = np.real(self.getDataType(y))
        self.root.addMacro(['baselineCorrection', (y, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.baselineCorrection(y, self.axes[-1], select=selectSlice, degree=degree, type=type)

    def previewBaselineCorrection(self, degree, removeList, type, invert=False):
        """
        Preview the baseline correction of a spectrum/fid.

        Parameters
        ----------
        degree: int
            Polynomial degree
        removeList: list
            Indexes of points not include in function fit
        type: str
            Either 'poly' for polynomial fit, or 'sin/cos' for sine/cosine
        invert (optional = False): boolean
            If True, the removeList is treated as an include list (i.e. inverting the selection)           
        """
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        tmpAx = np.arange(self.len())
        bArray = np.array([True] * self.len())
        for i in range(int(np.floor(len(removeList) / 2.0))):
            minVal = min(removeList[2 * i], removeList[2 * i + 1])
            maxVal = max(removeList[2 * i], removeList[2 * i + 1])
            bArray = np.logical_and(bArray, np.logical_or((tmpAx < minVal), (tmpAx > maxVal)))
        if invert:
            bArray = np.logical_not(bArray)
        y = self.baselineFunctionFit(self.xax(), tmpData, bArray, degree, type)
        y = np.real(self.getDataType(y))
        self.resetPreviewRemoveList()
        if self.NDIM_PLOT > 1:
            if isinstance(self, CurrentContour):
                self.showFid()
            else:
                self.showFid(extraX=[self.xax()], extraY=[y]*self.len(-2), extraColor=['g'])
        else:
            self.showFid(extraX=[self.xax()], extraY=[y], extraColor=['g'])
        self.previewRemoveList(removeList, invert)
        self.upd()

    def previewRemoveList(self, removeList, invert=False):
        """
        Preview the removelist of a baseline correction.

        Parameters
        ----------
        removeList: list
            Indexes of points not include in polyfit
        invert (optional = False): boolean
            If True, the removeList is treated as an include list (i.e. inverting the selection)           
        """
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        self.resetPreviewRemoveList()
        lineColor = 'r'
        if invert:
            lineColor = 'w'
            self.removeListLines.append(self.ax.axvspan(self.xax()[0] * axMult, self.xax()[-1] * axMult, color='r'))
        for i in range(int(np.floor(len(removeList) / 2.0))):
            self.removeListLines.append(self.ax.axvspan(self.xax()[removeList[2 * i]] * axMult, self.xax()[removeList[2 * i + 1]] * axMult, color=lineColor))
        if len(removeList) % 2:
            self.removeListLines.append(self.ax.axvline(self.xax()[removeList[-1]] * axMult, c=lineColor, linestyle='--'))
        self.canvas.draw()

    def resetPreviewRemoveList(self):
        """
        Resets the preview remove list of the baseline correction.
        """
        if hasattr(self, 'removeListLines'):
            for i in self.removeListLines:
                i.remove()
        self.removeListLines = []

    def states(self):
        """
        Performs a States data conversion along the current dimension.
        """
        self.root.addMacro(['states', (self.axes[-1] - self.data.ndim(), )])
        self.data.states(self.axes[-1])
        self.upd()
        self.showFid()

    def statesTPPI(self):
        """
        Performs a States-TPPI data conversion along the current dimension.
        """
        self.root.addMacro(['statesTPPI', (self.axes[-1] - self.data.ndim(), )])
        self.data.statesTPPI(self.axes[-1])
        self.upd()
        self.showFid()

    def echoAntiEcho(self):
        """
        Performs an Echo-Antiecho data conversion along the current dimension.
        """
        self.root.addMacro(['echoAntiEcho', (self.axes[-1] - self.data.ndim(), )])
        self.data.echoAntiEcho(self.axes[-1])
        self.upd()
        self.showFid()

    def matrixFuncs(self, func, name, pos1, pos2, newSpec=False):
        """
        Function that handles multiple matrix operations.

        name: str
            Name of operation: integrate, sum, max, min, argmax, argmin, averge
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        newSpec (optional = False): boolean
            If True, a new spectrum class is returned, holding the output of the matrix function.

        Returns
        -------
        Spectrum object:
            If newSpec is True, a new spectrum object is returned.
        """
        if newSpec:
            tmpData = copy.deepcopy(self.data)
            func(tmpData, pos1, pos2, self.axes[-1])
            return tmpData
        self.root.addMacro([name, (pos1, pos2, self.axes[-1] - self.data.ndim(), )])
        func(self.data, pos1, pos2, self.axes[-1])
        if self.upd():
            self.showFid()
            self.plotReset()

    def integrate(self, pos1, pos2, newSpec=False):
        """
        Integrate all slices over the selected region.

        Parameters
        ----------
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        newSpec (optional = False): boolean
            If True, a new spectrum class is returned, holding the output of the matrix function.

        Returns
        -------
        Spectrum object:
            If newSpec is True, a new spectrum object is returned. Else, 'None' is returned.
        """
        return self.matrixFuncs(lambda obj, *args: obj.integrate(*args), 'integrate', pos1, pos2, newSpec)

    def sum(self, pos1, pos2, newSpec=False):
        return self.matrixFuncs(lambda obj, *args: obj.sum(*args), 'sum', pos1, pos2, newSpec)

    def max(self, pos1, pos2, newSpec=False):
        """
        Get the maximum of all slices over the selected region.

        Parameters
        ----------
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        newSpec (optional = False): boolean
            If True, a new spectrum class is returned, holding the output of the matrix function.

        Returns
        -------
        Spectrum object:
            If newSpec is True, a new spectrum object is returned. Else, 'None' is returned.
        """
        return self.matrixFuncs(lambda obj, *args: obj.max(*args), 'max', pos1, pos2, newSpec)

    def min(self, pos1, pos2, newSpec=False):
        """
        Get the minimum of all slices over the selected region.

        Parameters
        ----------
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        newSpec (optional = False): boolean
            If True, a new spectrum class is returned, holding the output of the matrix function.

        Returns
        -------
        Spectrum object:
            If newSpec is True, a new spectrum object is returned. Else, 'None' is returned.
        """
        return self.matrixFuncs(lambda obj, *args: obj.min(*args), 'min', pos1, pos2, newSpec)

    def argmax(self, pos1, pos2, newSpec=False):
        """
        Get the arg maxmimum of all slices over the selected region.

        Parameters
        ----------
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        newSpec (optional = False): boolean
            If True, a new spectrum class is returned, holding the output of the matrix function.

        Returns
        -------
        Spectrum object:
            If newSpec is True, a new spectrum object is returned. Else, 'None' is returned.
        """
        return self.matrixFuncs(lambda obj, *args: obj.argmax(*args), 'argmax', pos1, pos2, newSpec)

    def argmin(self, pos1, pos2, newSpec=False):
        """
        Get the arg minimum of all slices over the selected region.

        Parameters
        ----------
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        newSpec (optional = False): boolean
            If True, a new spectrum class is returned, holding the output of the matrix function.

        Returns
        -------
        Spectrum object:
            If newSpec is True, a new spectrum object is returned. Else, 'None' is returned.
        """
        return self.matrixFuncs(lambda obj, *args: obj.argmin(*args), 'argmin', pos1, pos2, newSpec)

    def average(self, pos1, pos2, newSpec=False):
        """
        Get the average of all slices over the selected region.

        Parameters
        ----------
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        newSpec (optional = False): boolean
            If True, a new spectrum class is returned, holding the output of the matrix function.

        Returns
        -------
        Spectrum object:
            If newSpec is True, a new spectrum object is returned. Else, 'None' is returned.
        """
        return self.matrixFuncs(lambda obj, *args: obj.average(*args), 'average', pos1, pos2, newSpec)

    def flipLR(self):
        """
        Flip (i.e. mirror) the data over the x-axis.
        """
        self.data.flipLR(self.axes[-1])
        self.upd()
        self.showFid()
        self.root.addMacro(['flipLR', (self.axes[-1] - self.data.ndim(), )])

    def concatenate(self, axes):
        """
        Concatenate the data along an input axis.

        Parameters
        ----------
        axes: int
            The concatenation axis

        """
        self.data.concatenate(axes)
        self.upd()
        self.showFid()
        self.plotReset()
        self.root.addMacro(['concatenate', (axes - self.data.ndim() - 1, )])

    def split(self, sections):
        """
        Split the data long the current dimension in part

        Parameters
        ----------
        sections: int
            Split in this number of sections.
        """
        self.data.split(sections, self.axes[-1])
        self.upd()
        self.showFid()
        self.plotReset()
        self.root.addMacro(['split', (sections, self.axes[-1] - self.data.ndim() + 1)])

    def diff(self):
        """
        Get the difference between data points along the current dimension.
        This reduces the size of the data along this dimension by 1.

        Diff is taken from left to right in the fid, but from right to left in the spectrum (i.e. the spectrum is
        displayed in a mirrored way).
        """
        self.data.diff(self.axes[-1])
        self.upd()
        self.showFid()
        self.root.addMacro(['diff', (self.axes[-1] - self.data.ndim(), )])

    def cumsum(self):
        """
        Get the cumulative sum of the data points along the current dimension.

        Cumsum is taken from left to right in the fid, but from right to left in the spectrum (i.e. the spectrum is
        displayed in a mirrored way).
        """
        self.data.cumsum(self.axes[-1])
        self.upd()
        self.showFid()
        self.root.addMacro(['cumsum', (self.axes[-1] - self.data.ndim(), )])

    def insert(self, data, pos):
        """
        Insert data at a position along the current dimension.

        Parameters
        ----------
        data: hypercomplex data class
            The data to be inserted
        pos: int
            Position where the data should be inserted            
        """

        self.root.addMacro(['insert', (data, pos, self.axes[-1] - self.data.ndim())])
        self.data.insert(data, pos, self.axes[-1])
        self.upd()
        self.showFid()
        self.plotReset()

    def delete(self, pos):
        """
        Delete the data from specified positions along the current dimension.

        Parameters
        ----------
        pos: int or array_like
            The indices to remove.
        """
        self.data.delete(pos, self.axes[-1])
        self.upd()
        self.showFid()
        self.root.addMacro(['delete', (pos, self.axes[-1] - self.data.ndim())])

    def deletePreview(self, pos):
        """
        Preview the delete operation.

        Parameters
        ----------
        pos: int or array_like
            The indices to remove.
        """
        self.data1D.delete(pos, -1)
        self.showFid()
        self.upd()

    def add(self, data, select=False):
        """
        Adds (sums) data to the current data.

        Parameters
        ----------
        data: hypercomplex data class
            The data to be added.
        select (optional = False): boolean
            If True, apply only to the current slice.
        """
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['add', (data, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.add(data, self.axes[-1], select=selectSlice)
        self.upd()
        self.showFid()

    def subtract(self, data, select=False):
        """
        Subtract data from to the current data.

        Parameters
        ----------
        data: hypercomplex data class
            The data to be subtracted.
        select (optional = False): boolean
            If True, apply only to the current slice.
        """
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['subtract', (data, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.subtract(data, self.axes[-1], select=selectSlice)
        self.upd()
        self.showFid()

    def multiply(self, data, select=False):
        """
        Multiply the current data with extra data. Note that a complex data
        multiplication is used.

        Parameters
        ----------
        data: hypercomplex data class
            The data to be multiplied with.
        select (optional = False): boolean
            If True, apply only to the current slice.
        """
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['multiply', (data, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.multiply(data, self.axes[-1], select=selectSlice)
        self.upd()
        self.showFid()

    def multiplyPreview(self, data):
        """
        Preview the multiplication of data to the current data.

        Parameters
        ----------
        data: hypercomplex data class
            The data to be multiplied with.
        """
        self.data1D.multiply(data, -1)
        self.showFid()
        self.upd()

    def divide(self, data, select=False):
        """
        Divide the current data with extra data. Note that a complex data
        division is used.

        Parameters
        ----------
        data: hypercomplex data class
            The data to be divided with.
        select (optional = False): boolean
            If True, apply only to the current slice.
        """
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['divide', (data, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.divide(data, self.axes[-1], select=selectSlice)
        self.upd()
        self.showFid()

    def normalize(self, value, scale, type, select=False):
        """
        Normalize the data, relative to a selected region.
        Different types are supported: maximum, minimum, and integral.

        Parameters
        ----------
        value: float
            Value necessary to scale the selected region to 1.
        scale: float
            Extra scaling to get to this value
        type: int
            0:, integral. 1: max. 2: min
        select (optional = False): boolean
            If True, apply only to the current slice.
        """
        if select:
            selectSlice = self.getSelect()
        else:
            selectSlice = slice(None)
        self.root.addMacro(['normalize', (value, scale, type, self.axes[-1] - self.data.ndim(), selectSlice)])
        self.data.normalize(value, scale, type, self.axes[-1], select=selectSlice)
        self.upd()
        self.showFid()

    def subtractAvg(self, pos1, pos2):
        """
        Subtract the average of a region off all slices from that slice. This can be used
        to correct an offset, for example.

        Parameters
        ----------
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        """
        self.root.addMacro(['subtractAvg', (pos1, pos2, self.axes[-1] - self.data.ndim())])
        self.data.subtractAvg(pos1, pos2, self.axes[-1])
        self.upd()
        self.showFid()

    def subtractAvgPreview(self, pos1, pos2):
        """
        Preview of the effect of the "subtractAvg" function.

        Parameters
        ----------
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        """
        self.data1D.subtractAvg(pos1, pos2, -1)
        self.showFid()
        self.upd()

    def extract(self, pos1, pos2, newSpec=False, step=1):
        """
        Extract a region from the data (i.e. remove data outside the region.
        Can return a new data class if newSpec==True.

        Parameters
        ----------
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        newSpec (optional = False): boolean
            Return a new data class object
        step : int, optional
            Steps between indices to use

        Returns
        -------
        Spectrum object:
            If newSpec is True, a new spectrum object is returned. Else, 'None' is returned.
        """
        if newSpec:
            tmpData = copy.deepcopy(self.data)
            tmpData.extract(pos1, pos2, self.axes[-1], step)
            return tmpData
        self.root.addMacro(['extract', (pos1, pos2, self.axes[-1] - self.data.ndim(), step)])
        self.data.extract(pos1, pos2, self.axes[-1], step)
        self.upd()
        self.showFid()
        self.plotReset()

    def fiddle(self, pos1, pos2, lb):
        """
        Do a reference deconvolution using the fiddle algorithm.

        Parameters
        ----------
        pos1: int
            First data point limit of the selected region.
        pos2: int
            Second data point limit of the selected region.
        lb: float
            Added line broadening (Lorentzian) in Hz. This is needed to avoid artifacts. 
        """
        minPos = min(pos1, pos2)
        maxPos = max(pos1, pos2)
        tmpData = self.data1D.getHyperData(0)
        tmpData = tmpData[(0,)*(self.ndim()-1) + (slice(None), )]
        refSpec = np.zeros(self.len())
        refSpec[minPos:maxPos] = np.real(tmpData[minPos:maxPos])
        self.root.addMacro(['fiddle', (refSpec.tolist(), lb, self.axes[-1] - self.data.ndim())])
        self.data.fiddle(refSpec, lb, self.axes[-1])
        self.upd()
        self.showFid()

    def shearing(self, shear, axes, axes2, toRef=False):
        """
        Apply a shearing transform to the data.

        Parameters
        ----------
        shear: float
            Shearing constant
        axes: int
            The axis number over which the amount of shearing differs.
        axes2: int
            The axis over which the data must be rolled.
        toRef (optional = False): boolean
            Whether shearing should be relative to the reference (otherwise to the centre of the spectrum)
        """
        self.root.addMacro(['shear', (shear, axes - self.data.ndim(), axes2 - self.data.ndim()), toRef])
        self.data.shear(shear, axes, axes2, toRef)
        self.upd()
        self.showFid()

    def reorder(self, pos, newLength):
        """
        Reorder the current data to the new positions. Missing points are set to zero, and 
        zeroes are appended to reach 'newLength'.

        Parameters
        ----------
        pos: list of ints
            List with the new positions of each data point
        newLength: int
            The new length of the data
        """
        self.root.addMacro(['reorder', (pos, newLength, self.axes[-1] - self.data.ndim())])
        self.data.reorder(pos, newLength, self.axes[-1])
        self.upd()
        self.showFid()

    def ffm(self, posList, typeVal):
        """
        Apply the Fast Forward Maximum Entropy reconstruction method (for NUS data).

        Parameters
        ----------
        pos: list of ints
            List with the measured (non-zero) positions
        typeVal : {0, 1, 2}
            The type of data to be reconstructed.
            0=complex, 1=States or States-TPPI, 2=TPPI.
        """
        self.root.addMacro(['ffm', (posList, typeVal, self.axes[-1] - self.data.ndim())])
        self.data.ffm(posList, typeVal, self.axes[-1])
        self.upd()
        self.showFid()

    def clean(self, posList, typeVal, gamma, threshold, maxIter):
        """
        Apply the CLEAN reconstruction method (for NUS data).

        Parameters
        ----------
        posList : array_like
            A list of indices that are recorded datapoints.
            All other datapoints will be reconstructed.
        typeVal : {0, 1, 2}
            The type of data to be reconstructed.
            0=complex, 1=States or States-TPPI, 2=TPPI.
        gamma : float
            Gamma value of the CLEAN calculation.
        threshold : float
            Stopping limit (0 < x < 1) (stop if residual intensity below this point).
        maxIter : int
            Maximum number of iterations.
        """
        self.root.addMacro(['clean', (posList, typeVal, self.axes[-1] - self.data.ndim(), gamma, threshold, maxIter)])
        self.data.clean(posList, typeVal, self.axes[-1], gamma, threshold, maxIter)
        self.upd()
        self.showFid()

    def ist(self, posList, typeVal, threshold, maxIter, tracelimit):
        """
        Apply the IST (Iterative Soft Thresholding) reconstruction method (for NUS data).

        Parameters
        ----------
        posList : array_like
            A list of indices that are recorded datapoints.
            All other datapoints will be reconstructed.
        typeVal : {0, 1, 2}
            The type of data to be reconstructed.
            0=complex, 1=States or States-TPPI, 2=TPPI.
        threshold : float
            threshold. The level (0 < x < 1) at which the data is cut every iteration.
        maxIter : int
            Maximum number of iterations.
        tracelimit : float
            Stopping limit (0 < x < 1) (stop if residual intensity below this point).
        """
        self.root.addMacro(['ist', (posList, typeVal, self.axes[-1] - self.data.ndim(), threshold, maxIter, tracelimit)])
        self.data.ist(posList, typeVal, self.axes[-1], threshold, maxIter, tracelimit)
        self.upd()
        self.showFid()

    def autoPhase(self, phaseNum):
        """
        Automatically phase the data along the current dimension.
        This function returns the phasing answers, but does not execute the phasing yet.

        Parameters
        ----------
        phaseNum: int
            Order up to which to perform the autophasing.
            For 0 only zero order phasing is performed, for 1 both zero and first order phasing is performed.

        Returns
        -------
        list:
            List with 0th and 1st order phase.
        """
        phases = self.data1D.autoPhase(phaseNum, -1, [0]*self.ndim(), returnPhases=True)
        self.upd()
        return phases

    def directAutoPhase(self, phaseNum):
        """
        Automatically phase the data along the current dimension.
        This function applies the phasing directly to the data.

        Parameters
        ----------
        phaseNum: int
            Order up to which to perform the autophasing.
            For 0 only zero order phasing is performed, for 1 both zero and first order phasing is performed.
        """
        tmpLocList = copy.copy(self.locList)
        if self.ndim() > 1:
            if self.viewSettings["stackBegin"] is None:
                tmpLocList[self.axes[-1]] = 0
            else:
                tmpLocList[self.axes[-1]] = self.viewSettings["stackBegin"]
        self.root.addMacro(['autoPhase', (phaseNum, self.axes[-1] - self.data.ndim(), tmpLocList)])
        self.data.autoPhase(phaseNum, self.axes[-1], tmpLocList)
        self.upd()
        self.showFid()

    def autoPhaseAll(self, phaseNum):
        """
        Automatically phase the data along the current dimension, for
        each trace individually.

        Parameters
        ----------
        phaseNum: int
            Order up to which to perform the autophasing.
            For 0 only zero order phasing is performed, for 1 both zero and first order phasing is performed.
        """
        self.root.addMacro(['autoPhaseAll', (phaseNum, self.axes[-1] - self.data.ndim())])
        self.data.autoPhaseAll(phaseNum, self.axes[-1])
        self.upd()
        self.showFid()

    def setXaxPreview(self, xax):
        """
        Preview the plot with a new x-axis.

        Parameters
        ----------
        xax : array_like
            The x-axis.
            It should have the same length as the size of the data along dimension axis.
        """
        self.data1D.setXax(xax, -1)
        self.showFid()
        self.plotReset()
        self.upd()

    def setXax(self, xax):
        """
        Change the x-axis of the data.

        Parameters
        ----------
        xax : array_like
            The x-axis.
            It should have the same length as the size of the data along dimension axis.
        """
        self.root.addMacro(['setXax', (xax, self.axes[-1] - self.data.ndim())])
        self.data.setXax(xax, self.axes[-1])
        self.upd()
        self.showFid()
        self.plotReset()

    def setAxType(self, val, update=True, num=-1):
        """
        Change the axis type if the x-axis.

        The type can be 0,1,2 or 3.
        For a spectrum axis:
            0: Hz
            1: kHz
            2: MHz
            3: ppm

        For an FID axis:
            0: s
            1: ms
            2: us

        Parameters
        ----------
        val: int
            The new axis type
        update (optional = True): boolean
            If True, update the displays with the new axis.
        num (optional = -1): int
            Which axis to change (default -1 is the x-axis, -2 would be the y axis, etc.)
        """
        oldAxMult = self.getAxMult(self.spec(num), self.getAxType(num), self.getppm(num), self.freq(num), self.ref(num))
        if val == 3:
            self.viewSettings["ppm"][num] = True
        else:
            self.viewSettings["ppm"][num] = False
            self.viewSettings["axType"][num] = val
        newAxMult = self.getAxMult(self.spec(num), self.getAxType(num), self.getppm(num), self.freq(num), self.ref(num))
        if num in (-1, self.ndim()-1):
            self.xminlim = self.xminlim * newAxMult / oldAxMult
            self.xmaxlim = self.xmaxlim * newAxMult / oldAxMult
        elif num in (-2, self.ndim()-2):
            self.yminlim = self.yminlim * newAxMult / oldAxMult
            self.ymaxlim = self.ymaxlim * newAxMult / oldAxMult
        if update:
            self.showFid()

    def hilbert(self):
        """
        Apply a Hilbert transform along the current dimension.
        This reconstructs the imaginary part based on the real part of the data.
        """
        self.root.addMacro(['hilbert', (self.axes[-1] - self.data.ndim(), )])
        self.data.hilbert(self.axes[-1])
        self.upd()
        self.showFid()

    def getColorMap(self):
        """
        Get the current color map

        Returns
        -------
        colormap object
        """
        return COLORMAPLIST.index(self.viewSettings["colorMap"])

    def setColorMap(self, num):
        """
        Set the color map to an input number.

        Parameters
        ----------
        num: int
            The number of the color map.
        """
        self.viewSettings["colorMap"] = COLORMAPLIST[num]

    def getPColorMap(self):
        """
        Get the current color map

        Returns
        -------
        colormap object
        """
        return COLORMAPLIST.index(self.viewSettings["pColorMap"])

    def setPColorMap(self, num):
        """
        Set the color map to an input number.

        Parameters
        ----------
        num: int
            The number of the color map.
        """
        self.viewSettings["pColorMap"] = COLORMAPLIST[num]


    def getColorRange(self):
        """
        Returns the name of the current color range.
        """
        return COLORRANGELIST.index(self.viewSettings["colorRange"])

    def setColorRange(self, num):
        """
        Set the color range to an input number.

        Parameters
        ----------
        num: int
            The number of the color range.
        """
        self.viewSettings["colorRange"] = COLORRANGELIST[num]

    def setColor(self, color):
        """
        Set line color.

        Parameters
        ----------
        color: string
            Color string like '#1F77B4'
        """
        self.viewSettings["color"] = color

    def setLw(self, lw):
        """
        Set the line width of the plot line.

        Parameters
        lw: float:
            The new line width
        """
        self.viewSettings["linewidth"] = lw

    def setTickNum(self, x, y):
        """
        Set suggested minimum number of ticks for the x and y axis.

        Parameters
        ----------
        x: int
            Number of x ticks
        y: int
            Number of y ticks
        """
        self.viewSettings["minXTicks"] = x
        self.viewSettings["minYTicks"] = y

    def setContourColors(self, colors):
        """
        Sets the positive/negative colors for the contour plot.

        Parameters
        ----------
        colors: list of color strings
            The positive/negative color string
        """
        self.viewSettings["contourColors"] = colors

    def setContourConst(self, constant):
        """
        Parameters
        ----------
        constant: bool
            If True, use constant colors for negative/positive contours. 
            If False, use the color gradient
        """
        self.viewSettings["contourConst"] = constant

    def getOOM(self):
        """
        Get order of magnitude for the intensity of the current data.

        Returns
        -------
        int:
            Order of magnitude
        """
        absVal = np.max(np.abs(self.data.getHyperData(0)))
        if absVal == 0.0:
            return 1
        return int(np.floor(np.log10(absVal)))

    def showFid(self, oldData=None, extraX=None, extraY=None, extraColor=None):
        """
        Display the data
        
        Parameters
        ----------
        oldData (optional = None): hypercomplex data type
            The old data, to display under the current (i.e. during apodization).
        extraX (optional = None): list of ndarrays
            List of extra x-axes for 1D data curves
        extray (optional = None): list of ndarrays
            List of extra intensity data for 1D data curves
        extraColor (optional = None): list of color strings
            List of color strings for the extra data
        """
        self.peakPickReset()
        tmpdata = self.data1D.getHyperData(0)
        self.ax.cla()
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        self.line_xdata = [self.xax() * axMult]
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        if self.len() == 1:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = self.MARKER
            linestyle = self.LINESTYLE
        if oldData is not None:
            tmp = np.real(self.getDataType(oldData.getHyperData(0)))
            self.line_xdata_extra.append(self.xax() * axMult)
            self.line_ydata_extra.append(tmp)
            self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, c=(0,0,0,0.2), linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if extraX is not None:
            for num, _ in enumerate(extraX):
                self.line_xdata_extra.append(extraX[num] * axMult)
                self.line_ydata_extra.append(extraY[num])
                if extraColor is None:
                    self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker='', linestyle='-', linewidth=self.viewSettings["linewidth"], picker=True)
                else:
                    if len(extraColor) < len(extraY):
                        color = extraColor[0]
                    else:
                        color = extraColor[num]
                    self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker='', linestyle='-', c=color, linewidth=self.viewSettings["linewidth"], picker=True)
        tmpdata = self.getDataType(tmpdata)
        if self.viewSettings["plotType"] == 2:
            self.line_xdata.append(self.line_xdata[-1])
            self.line_ydata = [np.imag(tmpdata), np.real(tmpdata)]
            self.ax.plot(self.line_xdata[-2], self.line_ydata[-2], marker=marker, linestyle=linestyle, c='#FF7F0E', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
        else:
            self.line_ydata = [np.real(tmpdata)]
        self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
        self.ax.set_xlabel(self.getLabel(self.spec(), self.axes[-1], self.getAxType(), self.getppm()))
        if self.logx:
            self.ax.set_xscale('log')
        else:
            self.ax.set_xscale('linear')
            self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        if self.logy:
            self.ax.set_yscale('log')
        else:
            self.ax.set_yscale('linear')
            self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        if self.spec() > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.setTicks()
        self.canvas.draw()

    def setTicks(self, Xset=True, Yset=True):
        """
        Set ticks for the current plot.

        Parameters
        ----------
        Xset (optional = True): boolean
            If False, do not draw x-ticks
        Yset (optional = True): boolean
            If False, do not draw y-ticks
        """
        if  matplotlib.__version__[0] > '1':
            if Xset:
                self.ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=1, steps=[1, 2, 2.5, 5, 10], min_n_ticks=self.viewSettings["minXTicks"]))
            if Yset:
                self.ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=1, steps=[1, 2, 2.5, 5, 10], min_n_ticks=self.viewSettings["minYTicks"]))

    def plotReset(self, xReset=True, yReset=True):  
        """
        Reset plot limits.

        Parameters
        ----------
        xReset (optional = True): boolean
            If True, reset the x-axis limits
        yReset (optional = True): boolean
            If True, reset the y-axis limits
        """
        miny = np.min(self.line_ydata)
        maxy = np.max(self.line_ydata)
        for line in self.line_ydata_extra:
            miny = min(miny, min(line))
            maxy = max(maxy, max(line))
        if miny == maxy: # Prevents setting the limits equal
            miny -= 0.01
            maxy += 0.01
        if isinstance(self, CurrentContour):
            self.plotReset_x_ax()
            self.plotReset_y_ax()            
            differ = 0
        else:
            differ = 0.05 * (maxy - miny)  # amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny - differ
            self.ymaxlim = maxy + differ
        if self.ref() == 0.0:
            self.viewSettings["ppm"][-1] = False
        minx = np.min(self.line_xdata)
        maxx = np.max(self.line_xdata)
        for line in self.line_xdata_extra:
            minx = min(minx, min(line))
            maxx = max(maxx, max(line))
        if minx == maxx: # Prevents setting the limits equal
            minx -= 0.01
            maxx += 0.01
        if xReset:
            self.xminlim = minx
            self.xmaxlim = maxx
        if self.spec() > 0 and type(self) is not CurrentArrayed:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        if isinstance(self, CurrentContour): #If contour: reverse y-axis
            if  self.spec(-2) > 0:
                self.ax.set_ylim(self.ymaxlim, self.yminlim)
        self.canvas.draw()


#########################################################################################################
class CurrentScatter(Current1D):

    X_RESIZE = False
    Y_RESIZE = False
    MARKER = 'o'
    LINESTYLE = 'none'


#########################################################################################################
# the class from which the data of multiple spectra is displayed, the operations which only edit the content of this class are for previewing


class CurrentMulti(Current1D):

    X_RESIZE = False
    Y_RESIZE = True

    def setExtraSlice(self, extraNum, axes, locList):
        """
        Change the slice of one of the extra data sets

        Parameters
        ----------
        extraNum: int
            Index of the extra data
        axes: 1darray
            The new axis
        locList: list
            New slice information for this data
        """
        self.viewSettings["extraAxes"][extraNum] = axes
        self.viewSettings["extraLoc"][extraNum] = locList

    def addExtraData(self, data, name):
        """
        Add extra data to the multiview

        Parameters
        ----------
        data: spectrum class data
            The extra data
        name: str
            The name of the extra data
        """
        self.viewSettings["extraName"].append(name)
        self.viewSettings["extraData"].append(data)
        self.viewSettings["extraLoc"].append([0] * (len(self.viewSettings["extraData"][-1].shape())))
        self.viewSettings["extraColor"].append(COLORCONVERTER.to_rgb(COLORCYCLE[np.mod(len(self.viewSettings["extraData"]), len(COLORCYCLE))]['color']))  # find a good color system
        self.viewSettings["extraAxes"].append([len(data.shape()) - 1])
        self.viewSettings["extraScale"].append(1.0)
        self.viewSettings["extraOffset"].append(0.0)
        self.viewSettings["extraShift"].append(0.0)
        self.showFid()

    def delExtraData(self, num):
        """
        Delete extra data

        Parameters
        ----------
        num: int
            Index of the data to be removed.
        """
        del self.viewSettings["extraData"][num]
        del self.viewSettings["extraLoc"][num]
        del self.viewSettings["extraColor"][num]
        del self.viewSettings["extraName"][num]
        del self.viewSettings["extraAxes"][num]
        del self.viewSettings["extraScale"][num]
        del self.viewSettings["extraOffset"][num]
        del self.viewSettings["extraShift"][num]
        self.showFid()

    def setExtraColor(self, num, color):
        """
        Set the color of a specified extra data set

        Parameters
        ----------
        num: int
            Index of the extra data
        color: tuple
            Color tuple (R,G,B,Alpha) of the new color
        """
        self.viewSettings["extraColor"][num] = color
        self.showFid()

    def getExtraColor(self, num):
        """
        Returns the colour tuple for a specified extra data

        Parameters
        ----------
        num: int
            Index of the extra data

        Returns
        -------
        tuple:
            The colour tuple
        """
        return tuple(np.array(255 * np.array(self.viewSettings["extraColor"][num]), dtype=int))

    def resetLocList(self):
        """
        Resets the location list (slices) of all data.
        """
        super(CurrentMulti, self).resetLocList()
        self.resetExtraLocList()

    def setExtraScale(self, num, scale):
        """
        Set the vertical scaling of additional plotted data.

        Parameters
        ----------
        num: int
            Index of the extra data
        scale: float
            The new scaling factor
        """
        self.viewSettings["extraScale"][num] = scale
        self.showFid()

    def setExtraOffset(self, num, offset):
        """
        Set the vertical offset of additional plotted data.

        Parameters
        ----------
        num: int
            Index of the extra data
        offset: float
            The new offset
        """
        self.viewSettings["extraOffset"][num] = offset
        self.showFid()

    def setExtraShift(self, num, shift):
        """
        Set the x offset of additional plotted data.

        Parameters
        ----------
        num: int
            Index of the extra data
        shift: float
            The new shift in units of the current axis
        """
        self.viewSettings["extraShift"][num] = shift
        self.showFid()
        
    def resetExtraLocList(self, num=None):
        """
        Resets the location list (active slice) of all extra data or
        a specific data set.

        Parameters
        ----------
        num (optional = None): int
            Index of the data to be adjusted. If None, reset all.
        """
        if num is None:
            for i in range(len(self.viewSettings["extraLoc"])):
                self.viewSettings["extraLoc"][i] = [0] * (len(self.viewSettings["extraData"][i].shape()))
        else:
            self.viewSettings["extraLoc"][num] = [0] * (len(self.viewSettings["extraData"][num].shape()))

    def showFid(self, oldData=None, extraX=None, extraY=None, extraColor=None):
        """
        Plot all data.

        Parameters
        ----------
        oldData (optional = None): hypercomplex data type
            The old data, to display under the current (i.e. during apodization).
        extraX (optional = None): list of ndarrays
            List of extra x-axes for 1D data curves
        extray (optional = None): list of ndarrays
            List of extra intensity data for 1D data curves
        extraColor (optional = None): list of color strings
            List of color strings for the extra data
        """
        self.peakPickReset()
        tmpdata = self.data1D.getHyperData(0)
        self.ax.cla()
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        self.line_xdata = [self.xax() * axMult]
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        for i in range(len(self.viewSettings["extraData"])):
            data = self.viewSettings["extraData"][i]
            try:
                if self.viewSettings["extraData"][i].ndim() <= self.viewSettings["extraAxes"][i][-1]:
                    self.viewSettings["extraAxes"][i][-1] = len(self.viewSettings["extraData"][i].shape()) - 1
                    self.resetExtraLocList(i)
                extraData1D = data.getSlice(self.viewSettings["extraAxes"][i], self.viewSettings["extraLoc"][i])
            except Exception:
                self.resetExtraLocList(i)
                extraData1D = data.getSlice(self.viewSettings["extraAxes"][i], self.viewSettings["extraLoc"][i])
            xax = extraData1D.xaxArray[-1]
            spec = extraData1D.spec[-1]
            freq = extraData1D.freq[-1]
            if extraData1D.ref[-1] is not None:
                ref = extraData1D.ref[-1]
            else:
                ref = extraData1D.freq[-1]
            self.line_xdata_extra.append(self.viewSettings["extraShift"][i] + xax * self.getAxMult(spec, self.getAxType(), self.getppm(), freq, ref))
            if extraData1D.getHyperData(0).shape[-1] == 1:
                marker = 'o'
                linestyle = 'none'
            else:
                marker = self.MARKER
                linestyle = self.LINESTYLE
            extraData = np.real(self.getDataType(extraData1D.getHyperData(0)))
            self.line_ydata_extra.append(extraData * self.viewSettings["extraScale"][i] + self.viewSettings["extraOffset"][i])
            self.ax.plot(self.line_xdata_extra[-1],
                         self.line_ydata_extra[-1],
                         marker=marker, linestyle=linestyle,
                         c=self.viewSettings["extraColor"][i],
                         linewidth=self.viewSettings["linewidth"],
                         label=data.name,
                         picker=True)
        if self.len() == 1:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = ''
            linestyle = '-'
        if oldData is not None:
            tmp = np.real(self.getDataType(oldData.getHyperData(0)))
            self.line_xdata_extra.append(self.xax() * axMult)
            self.line_ydata_extra.append(tmp)
            self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, c=(0,0,0,0.2), linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if extraX is not None:
            for num, _ in enumerate(extraX):
                self.line_xdata_extra.append(extraX[num] * axMult)
                self.line_ydata_extra.append(extraY[num])
                if extraColor is None:
                    self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], picker=True)
                else:
                    if len(extraColor) < len(extraY):
                        color = extraColor[0]
                    else:
                        color = extraColor[num]
                    self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], c=color, picker=True)
        tmpdata = self.getDataType(tmpdata)
        if self.viewSettings["plotType"] == 2:
            self.line_xdata.append(self.line_xdata[-1])
            self.line_ydata = [np.imag(tmpdata), np.real(tmpdata)]
            self.ax.plot(self.line_xdata[-2], self.line_ydata[-2], marker=marker, linestyle=linestyle, c='#FF7F0E', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
        else:
            self.line_ydata = [np.real(tmpdata)]
        self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c=self.viewSettings["color"], linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
        self.ax.set_xlabel(self.getLabel(self.spec(), self.axes[-1], self.getAxType(), self.getppm()))
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.setTicks()
        if self.spec() > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()

#########################################################################################################
# the class from which the stacked data is displayed, the operations which only edit the content of this class are for previewing


class CurrentStacked(Current1D):

    X_RESIZE = False
    Y_RESIZE = True
    ZERO_SCROLL_ALLOWED = False
    NDIM_PLOT = 2

    def startUp(self, xReset=True, yReset=True):
        """
        Run when starting this plot.

        Parameters
        ----------
        xReset (optional = True): boolean
            Reset the x-axis if True

        yReset (optional = True): boolean
            Reset the y-axis if True
        """
        self.showFid()
        self.plotReset(xReset, yReset)

    def upd(self):  
        """
        Get new data from the data instance
        """
        if self.data.ndim() < 2:
            self.root.rescue()
            return False
        if self.data.ndim() <= self.axes[-1] or self.data.ndim() <= self.axes[-2] or self.axes[-1] == self.axes[-2]:
            self.axes = np.array([len(self.data.shape()) - 2, len(self.data.shape()) - 1])
        if len(self.locList) != self.data.ndim():
            self.resetLocList()
        stack = [reim.floatSlice(self.viewSettings["stackBegin"], self.viewSettings["stackEnd"], self.viewSettings["stackStep"]), slice(None)]
        try:
            self.data1D = self.data.getSlice(self.axes, self.locList, stack)
            if self.data1D is None:
                self.root.rescue()
        except Exception:
            self.resetLocList()
            self.data1D = self.data.getSlice(self.axes, self.locList, stack)
        if self.ref(-1) == 0.0 or self.freq(-1) == 0.0: #Check if ppm is allowed
            self.viewSettings["ppm"][-1] = False
        if self.ref(-2) == 0.0 or self.freq(-2) == 0.0: #Check if ppm is allowed
            self.viewSettings["ppm"][-2] = False
        return True

    def stackSelect(self, stackBegin, stackEnd, stackStep):
        """
        Select which data to plot in the stack plot.
        The data indexes go from stackBegin to stackEnd with stackStep as
        step size.

        Parameters
        ----------
        stackBegin: int
            Bgin value of the series
        stackEnd: int
            End value of the series. Note that in python, this value is not inluded (0:2 gives 0,1)
        stackStep: int
            Step size
        """
        self.viewSettings["stackBegin"] = stackBegin
        self.viewSettings["stackEnd"] = stackEnd
        self.viewSettings["stackStep"] = stackStep
        self.upd()
        self.showFid()
        self.plotReset(self.X_RESIZE, self.Y_RESIZE)

    def setSpacing(self, spacing):
        """
        Sets the vertical spacing between the different traces.

        Parameters
        ----------
        spacing: float
            The new spacing
        """
        self.viewSettings["spacing"] = spacing
        self.showFid()
        self.plotReset(self.X_RESIZE, self.Y_RESIZE)

    def resetSpacing(self, zlims=True):
        """
        Reset plot spacing

        Parameters
        ----------
        zlims:
            Not used
        """
        difference = np.diff(self.data1D.getHyperData(0), axis=0)
        if difference.size == 0:
            self.viewSettings["spacing"] = 0
        else:
            difference = self.getDataType(difference)
            tmpData = self.getDataType(self.data1D.getHyperData(0))
            difference = np.min((np.real(difference), np.imag(difference)))
            amp = np.max((np.real(tmpData), np.imag(tmpData))) - np.min((np.real(tmpData), np.imag(tmpData)))
            self.viewSettings["spacing"] = np.abs(difference) + 0.1 * amp

    def altScroll(self, event):
        """
        Scroll spacing

        Parameters
        ----------
        event: mouse event            
        """
        self.viewSettings["spacing"] = self.viewSettings["spacing"] * 1.1**event.step
        self.root.sideframe.scrollSpacing(self.viewSettings["spacing"])
        self.showFid()

    def altReset(self):
        """
        Reset the spacing.
        """
        self.resetSpacing()
        self.root.sideframe.scrollSpacing(self.viewSettings["spacing"])
        self.showFid()

    def showFid(self, oldData=None, extraX=None, extraY=None, extraColor=None):
        """
        Plot the data.

        Parameters
        ----------
        oldData (optional = None): hypercomplex data type
            The old data, to display under the current (i.e. during apodization).
        extraX (optional = None): list of ndarrays
            List of extra x-axes for 1D data curves
        extray (optional = None): list of ndarrays
            List of extra intensity data for 1D data curves
        extraColor (optional = None): list of color strings
            List of color strings for the extra data
        """
        self.peakPickReset()
        tmpdata = self.data1D.getHyperData(0)
        self.ax.cla()
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        tmp_line_xdata = self.xax() * axMult
        self.line_xdata = []
        self.line_ydata = []
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        if self.len() == 1:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = self.MARKER
            linestyle = self.LINESTYLE
        if oldData is not None:
            for num in range(len(oldData.getHyperData(0))):
                tmp = np.real(self.getDataType(oldData.getHyperData(0)[num]))
                self.line_xdata_extra.append(tmp_line_xdata)
                self.line_ydata_extra.append(num * self.viewSettings["spacing"] + tmp)
                self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, c=(0,0,0,0.2), linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if extraX is not None:
            tmpx = extraX[0] * axMult
            for num, _ in enumerate(extraY):
                self.line_xdata_extra.append(tmpx)
                self.line_ydata_extra.append(num * self.viewSettings["spacing"] + extraY[num])
                if extraColor is None:
                    self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], picker=True)
                else:
                    if len(extraColor) < len(extraY):
                        color = extraColor[0]
                    else:
                        color = extraColor[num]
                    self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], c=color, picker=True)
        tmpdata = self.getDataType(tmpdata)
        if self.viewSettings["colorRange"] == 'none':
            colorRange = None
        else:
            colorRange = get_cmap(self.viewSettings["colorRange"])
        for num, _ in enumerate(tmpdata):
            if self.viewSettings["plotType"] == 2:
                self.line_xdata.append(tmp_line_xdata)
                self.line_ydata.append(num * self.viewSettings["spacing"] + np.imag(tmpdata[num]))
                self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c='#FF7F0E', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
            self.line_xdata.append(tmp_line_xdata)
            self.line_ydata.append(num * self.viewSettings["spacing"] + np.real(tmpdata[num]))
            if colorRange is None:
                color = self.viewSettings["color"]
            else:
                color = colorRange(num/float(len(tmpdata)))
            self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c=color, linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
        self.ax.set_xlabel(self.getLabel(self.spec(), self.axes[-1], self.getAxType(), self.getppm()))
        if self.spec() > 0:
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.setTicks()
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.canvas.draw()

#########################################################################################################
# the class from which the arrayed data is displayed, the operations which only edit the content of this class are for previewing

class CurrentArrayed(CurrentStacked):

    X_RESIZE = True
    Y_RESIZE = False
    INVERT_X = False
    ZERO_SCROLL_ALLOWED = True

    def __init__(self, root, fig, canvas, data, duplicateCurrent=None):
        if duplicateCurrent is not None:
            if isinstance(duplicateCurrent, CurrentArrayed):
                self.zminlim = duplicateCurrent.zminlim
                self.zmaxlim = duplicateCurrent.zmaxlim
            else:
                # The z-axes limits are in xax units unlike the x-axes and y-axes limits
                axMult = self.getAxMult(duplicateCurrent.spec(),
                                        duplicateCurrent.getAxType(),
                                        duplicateCurrent.getppm(),
                                        duplicateCurrent.freq(),
                                        duplicateCurrent.ref())
                self.zminlim = (duplicateCurrent.xminlim) / axMult
                self.zmaxlim = (duplicateCurrent.xmaxlim) / axMult
        super(CurrentArrayed, self).__init__(root, fig, canvas, data, duplicateCurrent)

    def startUp(self, xReset=True, yReset=True):
        """
        Run when starting this plot.

        Parameters
        ----------
        xReset (optional = True): boolean
            Reset the x-axis if True

        yReset (optional = True): boolean
            Reset the y-axis if True
        """
        self.showFid()
        self.plotReset(xReset, yReset)

    def setAxType(self, val, update=True, num=-1):
        #Reimplement of base function. Prevent change of yaxis limits
        """
        Change the axis type if the x-axis.

        The type can be 0,1,2 or 3.
        For a spectrum axis:
            0: Hz
            1: kHz
            2: MHz
            3: ppm

        For an FID axis:
            0: s
            1: ms
            2: us

        Parameters
        ----------
        val: int
            The new axis type
        update (optional = True): boolean
            If True, update the displays with the new axis.
        num (optional = -1): int
            Which axis to change (default -1 is the x-axis, -2 would be the y axis, etc.)
        """
        yminlimBack = self.yminlim
        ymaxlimBack = self.ymaxlim
        super(CurrentArrayed, self).setAxType(val, False, num)
        self.yminlim = yminlimBack
        self.ymaxlim = ymaxlimBack
        if update:
            self.showFid()

    def resetSpacing(self, zlims=True):
        """
        Reset spacing

        Parameters
        ----------
        zlims (optional = True): boolean
            If True, reset the limits of the individual x-axes too.
        """
        if zlims:
            self.zminlim = min(self.xax())
            self.zmaxlim = max(self.xax())
        xaxZlims = (self.xax() > self.zminlim) & (self.xax() < self.zmaxlim)
        self.viewSettings["spacing"] = (self.xax()[xaxZlims][-1] - self.xax()[xaxZlims][0]) * 1.1

    def showFid(self, oldData=None, extraX=None, extraY=None, extraColor=None):
        """
        Plot the data.

        Parameters
        ----------
        oldData (optional = None): hypercomplex data type
            The old data, to display under the current (i.e. during apodization).
        extraX (optional = None): list of ndarrays
            List of extra x-axes for 1D data curves
        extray (optional = None): list of ndarrays
            List of extra intensity data for 1D data curves
        extraColor (optional = None): list of color strings
            List of color strings for the extra data
        """
        self.peakPickReset()
        tmpdata = self.data1D.getHyperData(0)
        self.ax.cla()
        if self.spec() > 0:
            direc = slice(None, None, -1)
        else:
            direc = slice(None, None, 1)
        xaxZlims = (self.xax() > self.zminlim) & (self.xax() < self.zmaxlim)
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        axMult2 = self.getAxMult(self.spec(-2), self.getAxType(-2), self.getppm(-2), self.freq(-2), self.ref(-2))
        self.line_xdata = []
        self.line_ydata = []
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        if self.len() == 1:
            marker = 'o'
            linestyle = 'none'
        else:
            marker = self.MARKER
            linestyle = self.LINESTYLE
        if oldData is not None:
            for num in range(len(oldData.getHyperData(0))):
                tmp = np.real(self.getDataType(oldData.getHyperData(0)[num]))
                self.line_xdata_extra.append((num * self.viewSettings["spacing"] + self.xax()[xaxZlims]) * axMult)
                self.line_ydata_extra.append(tmp[xaxZlims][direc])
                self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, c=(0,0,0,0.2), linewidth=self.viewSettings["linewidth"], label=self.data.name + '_old', picker=True)
        if extraX is not None:
            extraZlims = (extraX[0] > self.zminlim) & (extraX[0] < self.zmaxlim)
            for num, _ in enumerate(extraY):
                self.line_xdata_extra.append((num * self.viewSettings["spacing"] + extraX[0][extraZlims]) * axMult)
                self.line_ydata_extra.append(extraY[num][extraZlims][direc])
                if extraColor is None:
                    self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], picker=True)
                else:
                    if len(extraColor) < len(extraY):
                        color = extraColor[0]
                    else:
                        color = extraColor[num]
                    self.ax.plot(self.line_xdata_extra[-1], self.line_ydata_extra[-1], marker=marker, linestyle=linestyle, linewidth=self.viewSettings["linewidth"], c=color, picker=True)
        tmpdata = self.getDataType(tmpdata)
        ticksPos = []
        if self.viewSettings["colorRange"] == 'none':
            colorRange = None
        else:
            colorRange = get_cmap(self.viewSettings["colorRange"])
        for num, _ in enumerate(tmpdata):
            if self.viewSettings["plotType"] == 2:
                self.line_xdata.append((num * self.viewSettings["spacing"] + self.xax()[xaxZlims]) * axMult)
                self.line_ydata.append(np.imag(tmpdata[num][xaxZlims])[direc])
                self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c='#FF7F0E', linewidth=self.viewSettings["linewidth"], label=self.data.name + '_imag', picker=True)
            self.line_xdata.append((num * self.viewSettings["spacing"] + self.xax()[xaxZlims]) * axMult)
            self.line_ydata.append(np.real(tmpdata[num][xaxZlims])[direc])
            if colorRange is None:
                color = self.viewSettings["color"]
            else:
                color = colorRange(num/float(len(tmpdata)))
            self.ax.plot(self.line_xdata[-1], self.line_ydata[-1], marker=marker, linestyle=linestyle, c=color, linewidth=self.viewSettings["linewidth"], label=self.data.name, picker=True)
            pos = (num * self.viewSettings["spacing"] + 0.5 * (self.xax()[xaxZlims][-1] + self.xax()[xaxZlims][0])) * axMult
            ticksPos.append(pos)
        self.ax.set_xticks(ticksPos)
        self.ax.set_xticklabels([('%#.3g') % x for x in self.xax(-2) * axMult2])
        self.ax.set_xlim(self.xminlim, self.xmaxlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.setTicks(Xset=False)
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.canvas.draw()

######################################################################################################


def add_diagonal(axes, mult, *line_args, **line_kwargs):
    """
    Add a diagonal to the plot.

    Parameters
    ----------
    axes: matplotlib axes
    mult: diagonal multiplier (1 is true diagonal).
    """
    identity, = axes.plot([], [], *line_args, **line_kwargs)

    def callback(axes):
        low_x, high_x = axes.get_xlim()
        low_y, high_y = axes.get_ylim()
        if mult == 0:
            low = low_y
            high = high_y
        else:
            low_y = low_y / float(mult)
            high_y = high_y / float(mult)
            low = max(low_x, low_y)
            high = min(high_x, high_y)
        identity.set_data([low, high], [low * mult, high * mult])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes

#########################################################################################################
# the class from which the contour data is displayed, the operations which only edit the content of this class are for previewing


class CurrentContour(CurrentStacked):

    X_RESIZE = False
    Y_RESIZE = True
    GRID_PLOT = True
    INVERT_Y = True
    ZERO_SCROLL_ALLOWED = False

    def startUp(self, xReset=True, yReset=True):
        """
        Run when starting this plot.

        Parameters
        ----------
        xReset (optional = True): boolean
            Reset the x-axis if True

        yReset (optional = True): boolean
            Reset the y-axis if True
        """
        self.showFid()
        self.plotReset(xReset, yReset)

    def altScroll(self, event):
        """
        Scroll contour level limit.

        Parameters
        ----------
        event: mouse event            
        """
        minLevels = self.viewSettings["minLevels"] / 1.1**event.step
        if minLevels > 1:
            minLevels = 1
        if self.viewSettings["maxLevels"] > 1:
            self.viewSettings["maxLevels"] = 1
        self.viewSettings["minLevels"] = minLevels
        self.root.sideframe.minLEntry.setText(format(self.viewSettings["minLevels"] * 100, '.7g'))
        self.root.sideframe.maxLEntry.setText(str(self.viewSettings["maxLevels"] * 100))
        self.showFid()

    def plotReset_x_ax(self):
        if not self.line_xProjData:
            return
        minz = min([np.nanmin(i) for i in self.line_xProjData])
        maxz = max([np.nanmax(i) for i in self.line_xProjData])
        if np.isnan(minz): # If there are no acceptable values minz and maxz will be NaN
            minz = -0.01
            maxz = 0.01
        if minz == maxz: # Prevents setting the limits equal
            minz -= 0.01
            maxz += 0.01
        differ = 0.05 * (maxz - minz)  # amount to add to show all datapoints (10%)
        self.zminlim_x_ax = minz - differ
        self.zmaxlim_x_ax = maxz + differ
        self.x_ax.set_ylim(self.zminlim_x_ax, self.zmaxlim_x_ax)
        self.canvas.draw()

    def plotReset_y_ax(self):
        if not self.line_yProjData:
            return
        minz = min([np.nanmin(i) for i in self.line_yProjData])
        maxz = max([np.nanmax(i) for i in self.line_yProjData])
        if np.isnan(minz): # If there are no acceptable values minz and maxz will be NaN
            minz = -0.01
            maxz = 0.01
        if minz == maxz: # Prevents setting the limits equal
            minz -= 0.01
            maxz += 0.01
        differ = 0.05 * (maxz - minz)  # amount to add to show all datapoints (10%)
        self.zminlim_y_ax = minz - differ
        self.zmaxlim_y_ax = maxz + differ
        self.y_ax.set_xlim(self.zminlim_y_ax, self.zmaxlim_y_ax)
        self.canvas.draw()

    def setLevels(self, numLevels, maxLevels, minLevels, limitType, contourSign, contourType, multiValue):
        """
        Sets the contour settings

        Parameters
        ----------
        numLevels: int
            Number of contours
        maxLevels: float
            Maximum value (1 is max) of the contours
        minLevels: float
            Minimum level of the contours
        limitType: int
            0: relative to current 2D slice 1: relative to full data
        contourSign: int
            0: both, 1: + only 2: - only
        contourType: int
            0: linear 1: multiplier
        multiValue: float
            Value of the multiplier
        """
        if multiValue < 1:
            raise sc.SpectrumException("Contour level multiplier cannot be below 1")
        self.viewSettings["numLevels"] = numLevels
        self.viewSettings["maxLevels"] = maxLevels
        self.viewSettings["minLevels"] = minLevels
        self.viewSettings["limitType"] = limitType
        self.viewSettings["contourSign"] = contourSign
        self.viewSettings["contourType"] = contourType
        self.viewSettings["multiValue"] = multiValue
        self.showFid()

    def setProjLimits(self, ProjBool, Limits):
        """
        Set projection limits (i.e. ranges).

        Parameters
        ----------
        ProjBool: boolean
            If True, projection ranges are taken into account.
        Limits: list of 4 ints
            Slice positions that limit the projections
        """
        self.viewSettings["projLimits"] = Limits
        self.viewSettings["projLimitsBool"] = ProjBool
        self.clearProj()
        self.showAllProj()

    def setProjPos(self, pos):
        """
        Sets the projection slice, if a specific slice is plot as the projection.

        Parameters
        ----------
        pos: list of ints
            The slices to be taken
        """
        self.viewSettings["projPos"] = pos
        self.clearProj()
        self.showAllProj()

    def setProjType(self, val, direc):
        """
        Set the type of projection

        Parameters
        ----------
        val: int
            The type. 0: sum 1: max 2: min 3: off 4: slice 5: diagonal
        direct: int
            1: top 2: right
        """
        if direc == 1:
            self.viewSettings["projTop"] = val
            self.clearProj()
            self.showAllProj()
            self.plotReset_x_ax()
        if direc == 2:
            self.viewSettings["projRight"] = val
            self.clearProj()
            self.showAllProj()
            self.plotReset_y_ax()

    def setProjTraces(self, val, direc):
        """
        Set a specific trace for a projection.

        Parameters
        ----------
        val: int
            The trace index
        direct: int
            1: top 2: right
        """
        self.viewSettings["projPos"][direc] = val
        self.clearProj()
        self.showAllProj()

    def integralsPreview(self, xMin, xMax, yMin, yMax):
        """
        Draw different rectanglur patches, for a preview of 
        the intergral selection tool.

        Parameters
        ----------
        xMin: list of int
            Minimum x positions of the rectangles
        xMax: list of int
            Maximum x positions of the rectangles
        yMin: list of int
            Minimum y positions of the rectangles
        yMax: list of int
            Maximum y positions of the rectangles
        """
        nPatches = min(len(xMin), len(xMax), len(yMin), len(yMax))
        self.resetPreviewRemoveList()
        xax = self.xax()
        yax = self.xax(-2)
        xaxMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        yaxMult = self.getAxMult(self.spec(-2), self.getAxType(-2), self.getppm(-2), self.freq(-2), self.ref(-2))
        for i in range(nPatches):
            color = 'C'+str((i+2)%10)
            xminTmp = xax[xMin[i]] * xaxMult
            xmaxTmp = xax[xMax[i]] * xaxMult
            yminTmp = yax[yMin[i]] * yaxMult
            ymaxTmp = yax[yMax[i]] * yaxMult
            self.removeListLines.append(self.ax.fill([xminTmp, xminTmp, xmaxTmp, xmaxTmp], [yminTmp, ymaxTmp, ymaxTmp, yminTmp], color=color, fill=False, linestyle='--')[0])
        self.canvas.draw()

    def showFid(self, oldData=None, extraX=None, extraY=None, extraZ=None, extraColor=None):
        """
        Plot the data.

        Parameters
        ----------
        oldData (optional = None): hypercomplex data type
            The old data, to display under the current (i.e. during apodization).
        extraX (optional = None): list of ndarrays
            List of extra x-axes for data curves
        extraY (optional = None): list of ndarrays
            List of extra y-axes for data curves
        extraZ (optional = None): list of ndarrays
            List of extra intensity data for 1D data curves
        extraColor (optional = None): list of color strings
            List of color strings for the extra data
        """
        # The oldData and extra plots are not displayed in the contourplot for now
        self.line_xdata_extra = []
        self.line_ydata_extra = []
        self.line_zdata_extra = []
        self.line_color_extra = []
        self.differ = None
        self.peakPickReset()
        tmpdata = np.real(self.getDataType(self.data1D.getHyperData(0)))
        if self.viewSettings["limitType"] == 0:
            self.differ = np.max(np.abs(tmpdata))
        else:
            self.differ = np.max(np.abs(np.ravel(self.data.getHyperData(0))))
        self.ax.cla()
        self.clearProj()
        for i in range(len(self.viewSettings["extraData"])):
            data = self.viewSettings["extraData"][i]
            try:
                if self.viewSettings["extraData"][i].ndim() <= self.viewSettings["extraAxes"][i][-1]:
                    self.viewSettings["extraAxes"][i][-1] = len(self.viewSettings["extraData"][i].shape()) - 1
                    self.resetExtraLocList(i)
                extraData1D = data.getSlice(self.viewSettings["extraAxes"][i], self.viewSettings["extraLoc"][i])
            except Exception:
                self.resetExtraLocList(i)
                extraData1D = data.getSlice(self.viewSettings["extraAxes"][i], self.viewSettings["extraLoc"][i])
            xax = extraData1D.xaxArray[-1]
            spec = extraData1D.spec[-1]
            freq = extraData1D.freq[-1]
            if extraData1D.ref[-1] is not None:
                ref = extraData1D.ref[-1]
            else:
                ref = extraData1D.freq[-1]
            if extraData1D.ref[-2] is not None:
                ref2 = extraData1D.ref[-2]
            else:
                ref2 = extraData1D.freq[-2]
            axMult = self.getAxMult(extraData1D.spec[-1], self.getAxType(), self.getppm(), extraData1D.freq[-1], ref)
            axMult2 = self.getAxMult(extraData1D.spec[-2], self.getAxType(-2), self.getppm(-2), extraData1D.freq[-2], ref2)
            self.line_xdata_extra.append(self.viewSettings["extraShift"][i] + extraData1D.xaxArray[-1] * axMult)
            self.line_ydata_extra.append(self.viewSettings["extraShift2"][i] + extraData1D.xaxArray[-2] * axMult2)
            extraData = np.real(self.getDataType(extraData1D.getHyperData(0))) * self.viewSettings["extraScale"][i]
            self.line_zdata_extra.append(extraData)
            self.line_color_extra.append(self.viewSettings["extraColor"][i])
            self.plotContour(self.line_xdata_extra[-1], self.line_ydata_extra[-1], self.line_zdata_extra[-1], color=[self.viewSettings["extraColor"][i],tuple(j+(1-j)*0.5 for j  in self.viewSettings["extraColor"][i])])
        axMult = self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
        axMult2 = self.getAxMult(self.spec(-2), self.getAxType(-2), self.getppm(-2), self.freq(-2), self.ref(-2))
        if self.viewSettings["diagonalBool"]:
            add_diagonal(self.ax, self.viewSettings["diagonalMult"], c='k', ls='--')
        if oldData is not None:
            tmp = np.real(self.getDataType(oldData.getHyperData(0)))
            self.line_xdata_extra.append(self.xax() * axMult)
            self.line_ydata_extra.append(self.xax(-2) * axMult2)
            self.line_zdata_extra.append(tmp)
            self.line_color_extra.append('k')
            self.plotContour(self.line_xdata_extra[-1], self.line_ydata_extra[-1], self.line_zdata_extra[-1], color=[(0,0,0,0.2), (0,0,0,0.2)])
        if extraX is not None:
            for num, _ in enumerate(extraX):
                self.line_xdata_extra.append(extraX[num] * axMult)
                self.line_ydata_extra.append(extraY[num] * axMult2)
                self.line_zdata_extra.append(extraZ[num])
                if extraColor is None:
                    self.line_color_extra.append(extraColor)
                    self.plotContour(self.line_xdata_extra[-1], self.line_ydata_extra[-1], self.line_zdata_extra[-1])
                else:
                    if len(extraColor) < len(extraY):
                        color = extraColor[0]
                    else:
                        color = extraColor[num]
                    self.line_color_extra.append(color)
                    self.plotContour(self.line_xdata_extra[-1], self.line_ydata_extra[-1], self.line_zdata_extra[-1], color=[color, tuple(j+(1-j)*0.5 for j  in color)])
        self.line_xdata = [self.xax() * axMult]
        self.line_ydata = [self.xax(-2) * axMult2]
        self.line_zdata = [tmpdata]
        if isinstance(self, CurrentMultiContour):
            tmpColor = COLORCONVERTER.to_rgb(self.viewSettings["contourColors"][0])
            self.plotContour(self.line_xdata[-1], self.line_ydata[-1], self.line_zdata[-1], color=[tmpColor, tuple(j+(1-j)*0.5 for j  in tmpColor)])
        else:
            self.plotContour(self.line_xdata[-1], self.line_ydata[-1], self.line_zdata[-1])
        self.showAllProj()
        self.ax.set_xlabel(self.getLabel(self.spec(), self.axes[-1], self.getAxType(), self.getppm()))
        self.ax.set_ylabel(self.getLabel(self.spec(-2), self.axes[-2], self.getAxType(-2), self.getppm(-2)))
        if self.spec():
            self.ax.set_xlim(self.xmaxlim, self.xminlim)
        else:
            self.ax.set_xlim(self.xminlim, self.xmaxlim)
        if self.spec(-2):
            self.ax.set_ylim(self.ymaxlim, self.yminlim)
        else:
            self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.ax.get_xaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.ax.get_yaxis().get_major_formatter().set_powerlimits((-4, 4))
        self.x_ax.get_yaxis().get_major_formatter().set_powerlimits((-2, 2))
        self.y_ax.get_xaxis().get_major_formatter().set_powerlimits((-2, 2))
        self.ax.xaxis.grid(self.viewSettings["grids"][0])
        self.ax.yaxis.grid(self.viewSettings["grids"][1])
        self.setTicks()
        self.canvas.draw()
    
    def plotContour(self, line_xdata, line_ydata, line_zdata, color=None, updateOnly=False):  
        """
        Make the contour plot

        Parameters
        ----------
        line_xdata: 1darray
            xaxis
        line_ydata: 1darray
            yaxis
        line_zdata: 2darray
            Intensity (z) data
        color (optional = None): list of colors
            If not None, positive and negative contour colors should be in here
        updateOnly (optional = False): booleans
            If True, update only the contour plot
        """
        if color is None and self.viewSettings["contourConst"]:
            color = self.viewSettings["contourColors"]
        X, Y = np.meshgrid(line_xdata, line_ydata)
        if updateOnly:  # Set some extra stuff if only the contour plot needs updating
            del self.ax.collections[:]  # Clear all plot collections
        if self.viewSettings["contourType"] == 0:  # if linear
            contourLevels = np.linspace(self.viewSettings["minLevels"] * self.differ, self.viewSettings["maxLevels"] * self.differ, self.viewSettings["numLevels"])
        elif self.viewSettings["contourType"] == 1:  # if Multiplier
            contourLevels = [self.viewSettings["minLevels"] * self.differ]
            while contourLevels[-1] < self.viewSettings["maxLevels"] * self.differ and len(contourLevels) < self.viewSettings["numLevels"]:
                contourLevels.append(contourLevels[-1] * self.viewSettings["multiValue"])
            contourLevels = np.array(contourLevels)
        # Trim matrix of unused rows/columns for more efficient contour plotting
        PlotPositive = False
        if self.viewSettings["contourSign"] == 0 or self.viewSettings["contourSign"] == 1:
            if line_zdata.shape[0] > 2:  # if size 2 or lower, convolve fails, just take whole data then
                YposMax = np.where(np.convolve(np.max(line_zdata, 1) > contourLevels[0], [True, True, True], 'same'))[0]
            else:
                YposMax = np.arange(line_zdata.shape[0])
            if YposMax.size > 0:  # if number of positive contours is non-zero
                if line_zdata.shape[1] > 2:
                    XposMax = np.where(np.convolve(np.max(line_zdata, 0) > contourLevels[0], [True, True, True], 'same'))[0]
                else:
                    XposMax = np.arange(line_zdata.shape[1])
                PlotPositive = True
        PlotNegative = False
        if self.viewSettings["contourSign"] == 0 or self.viewSettings["contourSign"] == 2:
            if self.viewSettings["plotType"] != 3:  # for Absolute plot no negative
                if line_zdata.shape[0] > 2:
                    YposMin = np.where(np.convolve(np.min(line_zdata, 1) < -contourLevels[0], [True, True, True], 'same'))[0]
                else:
                    YposMin = np.arange(line_zdata.shape[0])
                if YposMin.size > 0:  # if number of negative contours is non-zero
                    if line_zdata.shape[1] > 2:
                        XposMin = np.where(np.convolve(np.min(line_zdata, 0) < -contourLevels[0], [True, True, True], 'same'))[0]
                    else:
                        XposMin = np.arange(line_zdata.shape[1])
                    PlotNegative = True
        vmax = max(np.abs(self.viewSettings["minLevels"] * self.differ), np.abs(self.viewSettings["maxLevels"] * self.differ))
        vmin = -vmax
        if len(line_xdata) > 1 and len(line_ydata) > 1: # Do not plot if too few points
            if color is not None:
                if PlotPositive:
                    self.ax.contour(X[YposMax[:, None], XposMax], Y[YposMax[:, None], XposMax], line_zdata[YposMax[:, None], XposMax], colors=[color[0]], levels=contourLevels, vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], linestyles='solid')
                if PlotNegative:
                    self.ax.contour(X[YposMin[:, None], XposMin], Y[YposMin[:, None], XposMin], line_zdata[YposMin[:, None], XposMin], colors=[color[1]], levels=-contourLevels[::-1], vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], linestyles='solid')
            else:
                if PlotPositive:
                    self.ax.contour(X[YposMax[:, None], XposMax], Y[YposMax[:, None], XposMax], line_zdata[YposMax[:, None], XposMax], cmap=get_cmap(self.viewSettings["colorMap"]), levels=contourLevels, vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], linestyles='solid')
                if PlotNegative:
                    self.ax.contour(X[YposMin[:, None], XposMin], Y[YposMin[:, None], XposMin], line_zdata[YposMin[:, None], XposMin], cmap=get_cmap(self.viewSettings["colorMap"]), levels=-contourLevels[::-1], vmax=vmax, vmin=vmin, linewidths=self.viewSettings["linewidth"], linestyles='solid')
        self.setTicks()
        if updateOnly:
            self.canvas.draw()

    def clearProj(self):
        """
        Clear the projections.
        """
        self.x_ax.cla()
        self.y_ax.cla()

    def showAllProj(self):
        self.line_xProjData = []
        self.line_yProjData = []
        for i, _ in enumerate(self.line_xdata_extra):
            self.showProj(self.line_xdata_extra[i], self.line_ydata_extra[i], self.line_zdata_extra[i], self.line_color_extra[i])
        self.showProj()
        
    def showProj(self, line_xdata=None, line_ydata=None, line_zdata=None, color=None):
        """
        Show the projections

        Parameters
        ----------
        line_xdata (optional = None): 1darray
            xaxis
        line_ydata (optional = None): 1darray
            yaxis
        line_zdata (optional = None): 2darray
            Intensity (z) data
        color  (optional = None): color string
            Colour of the projection lines
        """
        xLimOld = self.x_ax.get_xlim()
        if line_xdata is None:
            x = self.line_xdata[-1]
            y = self.line_ydata[-1]
            tmpdata = self.line_zdata[-1]
        else:
            x = line_xdata
            y = line_ydata
            tmpdata = line_zdata
        yLimOld = self.y_ax.get_ylim()
        if color is None:
            color = self.viewSettings["color"]
        Limits = self.viewSettings["projLimits"]
        topSlice = (slice(None), slice(None))
        rightSlice = (slice(None), slice(None))
        if self.viewSettings["projLimitsBool"] is True:
            if Limits[0] < Limits[1]:
                topSlice = (slice(Limits[0], Limits[1] + 1), slice(None))
            elif Limits[0] > Limits[1]:
                topSlice = (slice(Limits[1], Limits[0] + 1), slice(None))
            else:
                topSlice = (slice(Limits[0], Limits[0] + 1), slice(None))
            if Limits[2] < Limits[3]:
                rightSlice = (slice(None), slice(Limits[2], Limits[3] + 1))
            elif Limits[2] > Limits[3]:
                rightSlice = (slice(None), slice(Limits[3], Limits[2] + 1))
            else:
                rightSlice = (slice(None), slice(Limits[2], Limits[2] + 1))
        if self.viewSettings["projTop"] == 0:
            xprojdata = np.sum(tmpdata[topSlice], axis=0)
        elif self.viewSettings["projTop"] == 1:
            xprojdata = np.max(tmpdata[topSlice], axis=0)
        elif self.viewSettings["projTop"] == 2:
            xprojdata = np.min(tmpdata[topSlice], axis=0)
        elif self.viewSettings["projTop"] == 4:
            if self.viewSettings["projPos"][0] >= tmpdata.shape[self.axes[-2]]:
                self.viewSettings["projPos"][0] = tmpdata.shape[self.axes[-2]] - 1
            elif self.viewSettings["projPos"][0] < 0:
                self.viewSettings["projPos"][0] = 0
            xprojdata = tmpdata[self.viewSettings["projPos"][0]]
        elif self.viewSettings["projTop"] == 5:
            indices2 = np.searchsorted(y, x)
            indices1 = indices2 -1
            indices2[indices2 >= len(y)] = len(y)-1
            indices1[indices1 < 0] = 0
            dist1 = x - y[indices1]
            dist2 = y[indices2] - x
            distSum = (dist1 + dist2)
            xprojdata = dist2 * tmpdata[indices1, np.arange(len(x))] + dist1 * tmpdata[indices2, np.arange(len(x))]
            xprojdata[distSum != 0.0] /= distSum[distSum != 0.0]
            xprojdata[x > np.max(y)] = np.nan
            xprojdata[x < np.min(y)] = np.nan
        if self.viewSettings["projTop"] != 3:
            self.line_xProjData.append(xprojdata)
            self.x_ax.plot(x, self.line_xProjData[-1], color=color, linewidth=self.viewSettings["linewidth"], picker=True)
            #xmin, xmax = np.nanmin(xprojdata), np.nanmax(xprojdata)
            #self.x_ax.set_ylim([xmin - 0.15 * (xmax - xmin), xmax + 0.05 * (xmax - xmin)])  # Set projection limits, and force 15% whitespace below plot
            self.x_ax.set_ylim(self.zminlim_x_ax, self.zmaxlim_x_ax)
            self.x_ax.set_xlim(xLimOld)
        if self.viewSettings["projRight"] == 0:
            yprojdata = np.sum(tmpdata[rightSlice], axis=1)
        elif self.viewSettings["projRight"] == 1:
            yprojdata = np.max(tmpdata[rightSlice], axis=1)
        elif self.viewSettings["projRight"] == 2:
            yprojdata = np.min(tmpdata[rightSlice], axis=1)
        elif self.viewSettings["projRight"] == 4:
            if self.viewSettings["projPos"][1] >= tmpdata.shape[self.axes[-1]]:
                self.viewSettings["projPos"][1] = tmpdata.shape[self.axes[-1]] - 1
            elif self.viewSettings["projPos"][1] < 0:
                self.viewSettings["projPos"][1] = 0
            yprojdata = tmpdata[:, self.viewSettings["projPos"][1]]
        elif self.viewSettings["projTop"] == 5:
            indices2 = np.searchsorted(x, y)
            indices1 = indices2 -1
            indices2[indices2 >= len(x)] = len(x)-1
            indices1[indices1 < 0] = 0
            dist1 = y - x[indices1]
            dist2 = x[indices2] - y
            distSum = (dist1 + dist2)
            yprojdata = dist2 * tmpdata[np.arange(len(y)), indices1] + dist1 * tmpdata[np.arange(len(y)), indices2]
            yprojdata[distSum != 0.0] /= distSum[distSum != 0.0]
            yprojdata[y > np.max(x)] = np.nan
            yprojdata[y < np.min(x)] = np.nan
        if self.viewSettings["projRight"] != 3:
            self.line_yProjData.append(yprojdata)
            self.y_ax.plot(self.line_yProjData[-1], y, color=color, linewidth=self.viewSettings["linewidth"], picker=True)
            #ymin, ymax = np.nanmin(yprojdata), np.nanmax(yprojdata)
            #self.y_ax.set_xlim([ymin - 0.15 * (ymax - ymin), ymax + 0.05 * (ymax - ymin)])  # Set projection limits, and force 15% whitespace below plot
            self.y_ax.set_xlim(self.zminlim_y_ax, self.zmaxlim_y_ax)
            self.y_ax.set_ylim(yLimOld)
        self.setTicks()
        self.canvas.draw()

    # The peakpicking function needs to be changed for contour plots
    def buttonRelease(self, event):
        if event.button == 1:
            if self.peakPick:
                if self.rect[1] is not None:
                    self.rect[1].remove()
                    self.rect[1] = None
                if self.rect[0] is not None:
                    self.rect[0].remove()
                    self.rect[0] = None
                    self.peakPick = False
                    xdata = self.xax() * self.getAxMult(self.spec(), self.getAxType(), self.getppm(), self.freq(), self.ref())
                    ydata = self.xax(-2) * self.getAxMult(self.spec(-2), self.getAxType(-2), self.getppm(-2), self.freq(-2), self.ref(-2))
                    idx = np.argmin(np.abs(xdata - event.xdata))
                    idy = np.argmin(np.abs(ydata - event.ydata))
                    if self.peakPickFunc is not None:
                        tmpdata = np.real(self.getDataType(self.data1D.getHyperData(0)[idy, idx]))
                        self.peakPickFunc((idx, xdata[idx], tmpdata, idy, ydata[idy]))
                    if not self.peakPick:  # check if peakpicking is still required
                        self.peakPickFunc = None
            else:
                self.leftMouse = False
                if self.rect[0] is not None:
                    self.rect[0].remove()
                if self.rect[1] is not None:
                    self.rect[1].remove()
                if self.rect[2] is not None:
                    self.rect[2].remove()
                if self.rect[3] is not None:
                    self.rect[3].remove()
                self.rect = [None, None, None, None]
                if self.zoomX2 is not None and self.zoomY2 is not None:
                    if self.zoomAx == self.y_ax:
                        self.zminlim_y_ax = min([self.zoomX1, self.zoomX2])
                        self.zmaxlim_y_ax = max([self.zoomX1, self.zoomX2])
                        self.y_ax.set_xlim(self.zminlim_y_ax, self.zmaxlim_y_ax)
                    else:
                        self.xminlim = min([self.zoomX1, self.zoomX2])
                        self.xmaxlim = max([self.zoomX1, self.zoomX2])
                    if self.zoomAx == self.x_ax:
                        self.zminlim_x_ax = min([self.zoomY1, self.zoomY2])
                        self.zmaxlim_x_ax = max([self.zoomY1, self.zoomY2])
                        self.x_ax.set_ylim(self.zminlim_x_ax, self.zmaxlim_x_ax)
                    else:
                        self.yminlim = min([self.zoomY1, self.zoomY2])
                        self.ymaxlim = max([self.zoomY1, self.zoomY2])
                    if self.spec() > 0:
                        self.ax.set_xlim(self.xmaxlim, self.xminlim)
                    else:
                        self.ax.set_xlim(self.xminlim, self.xmaxlim)
                    if self.spec(-2) > 0:
                        self.ax.set_ylim(self.ymaxlim, self.yminlim)
                    else:
                        self.ax.set_ylim(self.yminlim, self.ymaxlim)
                self.zoomX1 = None
                self.zoomX2 = None
                self.zoomY1 = None
                self.zoomY2 = None
        elif event.button == 3:
            self.rightMouse = False
        self.canvas.draw()


class CurrentColour2D(CurrentContour):
    """
    2D colour plot class. Currently a child of CurrentContour. Probably a lower level 2D class
    should be made.
    """

    X_RESIZE = False
    Y_RESIZE = True
    GRID_PLOT = True
    INVERT_Y = True
    ZERO_SCROLL_ALLOWED = False
    
    def plotContour(self, line_xdata, line_ydata, line_zdata, color=None, updateOnly=False):  
        """
        Make the contour plot

        Parameters
        ----------
        line_xdata: 1darray
            xaxis
        line_ydata: 1darray
            yaxis
        line_zdata: 2darray
            Intensity (z) data
        color (optional = None): list of colors
            If not None, positive and negative contour colors should be in here
        updateOnly (optional = False): booleans
            If True, update only the contour plot
        """
        if updateOnly:  # Set some extra stuff if only the contour plot needs updating
            del self.ax.collections[:]  # Clear all plot collections

        vmax = np.max(np.abs(line_zdata))
        vmin = -vmax

        self.ax.imshow(np.flipud(line_zdata), extent=[line_xdata[0],line_xdata[-1],line_ydata[0],line_ydata[-1]],
                aspect='auto',cmap=get_cmap(self.viewSettings["pColorMap"]),vmax=vmax,vmin=vmin,interpolation='hanning')
        self.setTicks()
        if updateOnly:
            self.canvas.draw()

            
#########################################################################################################


class CurrentMultiContour(CurrentContour):

    def setExtraSlice(self, extraNum, axes, locList):
        """
        Change the slice of one of the extra data sets

        Parameters
        ----------
        extraNum: int
            Index of the extra data
        axes: 1darray
            The new axis
        locList: list
            New slice information for this data
        """
        self.viewSettings["extraAxes"][extraNum] = axes
        self.viewSettings["extraLoc"][extraNum] = locList

    def addExtraData(self, data, name):
        """
        Add extra data to the multiview

        Parameters
        ----------
        data: spectrum class data
            The extra data
        name: str
            The name of the extra data
        """
        if data.ndim() < 2:
            raise sc.SpectrumException("Data requires at least 2 dimensions to be added")
        self.viewSettings["extraName"].append(name)
        self.viewSettings["extraData"].append(data)
        self.viewSettings["extraLoc"].append([0] * (len(self.viewSettings["extraData"][-1].shape())))
        self.viewSettings["extraColor"].append(COLORCONVERTER.to_rgb(COLORCYCLE[np.mod(len(self.viewSettings["extraData"]), len(COLORCYCLE))]['color']))  # find a good color system
        self.viewSettings["extraAxes"].append([len(data.shape()) - 2, len(data.shape()) - 1])
        self.viewSettings["extraScale"].append(1.0)
        self.viewSettings["extraShift"].append(0.0)
        self.viewSettings["extraShift2"].append(0.0)
        self.showFid()

    def delExtraData(self, num):
        """
        Delete extra data

        Parameters
        ----------
        num: int
            Index of the data to be removed.
        """
        del self.viewSettings["extraData"][num]
        del self.viewSettings["extraLoc"][num]
        del self.viewSettings["extraColor"][num]
        del self.viewSettings["extraName"][num]
        del self.viewSettings["extraAxes"][num]
        del self.viewSettings["extraScale"][num]
        del self.viewSettings["extraShift"][num]
        del self.viewSettings["extraShift2"][num]
        self.showFid()

    def setExtraColor(self, num, color):
        """
        Set the color of a specified extra data set

        Parameters
        ----------
        num: int
            Index of the extra data
        color: tuple
            Color tuple (R,G,B,Alpha) of the new color
        """
        self.viewSettings["extraColor"][num] = color
        self.showFid()

    def getExtraColor(self, num):
        """
        Returns the colour tuple for a specified extra data

        Parameters
        ----------
        num: int
            Index of the extra data

        Returns
        -------
        tuple:
            The colour tuple
        """
        return tuple(np.array(255 * np.array(self.viewSettings["extraColor"][num]), dtype=int))

    def setExtraScale(self, num, scale):
        """
        Set the vertical scaling of additional plotted data.

        Parameters
        ----------
        num: int
            Index of the extra data
        scale: float
            The new scaling factor
        """
        self.viewSettings["extraScale"][num] = scale
        self.showFid()

    def setExtraShift(self, num, shift):
        """
        Set the x offset of additional plotted data.

        Parameters
        ----------
        num: int
            Index of the extra data
        shift: float
            The new shift in units of the current axis
        """
        self.viewSettings["extraShift"][num] = shift
        self.showFid()

    def setExtraShift2(self, num, shift):
        """
        Set the y offset of additional plotted data.

        Parameters
        ----------
        num: int
            Index of the extra data
        shift: float
            The new shift in units of the current axis
        """
        self.viewSettings["extraShift2"][num] = shift
        self.showFid()
    
    def resetLocList(self):
        """
        Resets the location list (slices) of all data.
        """
        super(CurrentMultiContour, self).resetLocList()
        self.resetExtraLocList()

    def resetExtraLocList(self, num=None):
        """
        Resets the location list (active slice) of all extra data or
        a specific data set.

        Parameters
        ----------
        num (optional = None): int
            Index of the data to be adjusted. If None, reset all.
        """
        if num is None:
            for i in range(len(self.viewSettings["extraLoc"])):
                self.viewSettings["extraLoc"][i] = [0] * (len(self.viewSettings["extraData"][i].shape()))
        else:
            self.viewSettings["extraLoc"][num] = [0] * (len(self.viewSettings["extraData"][num].shape()))


