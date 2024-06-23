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


def floatSlice(*args):
    """
    Function to create Slice objects with float input values.

    Parameters
    ----------
    *args
        All arguments that aren't None are converted to int.
        The arguments are then used to create the Slice object.

    Returns
    -------
    Slice
        The output Slice object.
    """
    tmp = ()
    for arg in args:
        if arg is None:
            tmp += (arg,)
        else:
            tmp += (int(arg),)
    return slice(*tmp)
