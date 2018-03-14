#!/usr/bin/env python

# Copyright 2016 - 2018 Bas van Meerten and Wouter Franssen

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

#########################################################################
# the hyper complex data class

class HComplexData(object):

    def __init__(self, data=None, hyper=None):
        if data is None:
            self.data = np.array([], dtype=complex)
            self.hyper = np.array([])
        else:
            if hyper is None:
                # Data is not hypercomplex
                self.data = np.array([data], dtype=complex)
                self.hyper = np.array([0]) # The data is only complex in the last dimension
            else:
                if len(hyper) != len(data):
                    raise RuntimeError('')
                self.data = np.array(data, dtype=complex)
                self.hyper = np.array(hyper)
                
    def __len__(self):
        if self.data:
            return self.data[0]

    def __add__(self, x):
        if isinstance(x, HComplexData):
            tmpData = np.setdiff1d(x.hyper, self.hyper)
            tmpData[np.isin(self.hyper, x.hyper, assume_unique=True)] += x.data[np.isin(x.hyper, self.hyper, assume_unique=True)]
            tmpData
        else:
            return HComplexData(self.data+x, self.hyper)

    def __radd__(self, x):
        return self.__add__(x)

    def ndim(self):
        return self.data.ndim - 1 # One extra dimension to contain the hypercomplex information

    def isComplex(self, axis):
        if axis == (self.ndim()-1):
            return True
        return bool(np.max(self.hyper) & (2**axis))
        
