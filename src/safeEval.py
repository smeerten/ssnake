#!/usr/bin/env python

# Copyright 2016 - 2019 Bas van Meerten and Wouter Franssen

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
import re
import hypercomplex as hc
import scipy.special

def safeEval(inp, length=None, keywords=[], type='All', x=None):
    env = vars(np).copy()
    env.update(vars(hc).copy())
    env.update(vars(scipy.special).copy())
    env.update(vars(scipy.integrate).copy())
    env["locals"] = None
    env["globals"] = None
    env["__name__"] = None
    env["__file__"] = None
    env["__builtins__"] = {'None': None, 'False': False, 'True':True} # None
    env["slice"] = slice
    if length is not None:
        env["length"] = length
    if x is not None:
        env["x"] = x
    inp = re.sub('([0-9]+)[kK]', '\g<1>*1024', str(inp))
    for i in keywords:
        if i in inp:
            return inp
    try:
        val = eval(inp, env)
        if isinstance(val, str):
            return None
        if type == 'All':
            return val
        elif type == 'FI':  #single float/int type
            if isinstance(val, (float, int)) and not np.isnan(val) and not np.isinf(val):
                return val
            return None
        elif type == 'C': #single complex number
            if isinstance(val, (float, int, complex)) and not np.isnan(val) and not np.isinf(val):
                return val
            return None
    except Exception:
        return None
