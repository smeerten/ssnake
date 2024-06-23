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

import re
import numpy as np
import scipy.special
import hypercomplex as hc

def safeEval(inp, length=None, Type='All', x=None):
    """
    Creates a more restricted eval environment.
    Note that this method is still not acceptable to process strings from untrusted sources.

    Parameters
    ----------
    inp : str
        String to evaluate.
    length : int or float, optional
        The variable length will be set to this value.
        By default the variable length is not set.
    Type : {'All', 'FI', 'C'}, optional
        Type of expected output. 'All' will return all types, 'FI' will return a float or int, and 'C' will return a complex number.
        By default Type is set to 'All'
    x : array_like, optional
        The variable x is set to this variable,
        By default the variable x is not used.

    Returns
    -------
    Object
        The result of the evaluated string.
    """
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
    try:
        val = eval(inp, env)
        if isinstance(val, str):
            return None
        if Type == 'All':
            return val
        if Type == 'FI':  #single float/int type
            if isinstance(val, (float, int)) and not np.isnan(val) and not np.isinf(val):
                return val
            return None
        if Type == 'C': #single complex number
            if isinstance(val, (float, int, complex)) and not np.isnan(val) and not np.isinf(val):
                return val
            return None
    except Exception:
        return None
