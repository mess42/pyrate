#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and    Uwe Lippmann  uwe.lippmann@web.de

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

class OptimizableVariable(object):
    """
    Class that contains a float. Used to get a pointer on a variable.
    """
    def __init__(self,name, value=0.0, status=False):
        self.name =  name
        self.val = value
        self.status = status

class ClassWithOptimizableVariables(object):
    """
    Virtual class for all classes that contain an optimizable variable.
    """
    def __init__(self):
        self.listOfOptimizableVariables = []

    def createOptimizableVariable(self, name, value = 0.0, status=False):
        """
        Creates a new float that can be used with optimization algorithms.

        :param name: name of the variable (str)
        :param value: initial value of the variable (float)
        :param status: this variable may be changed during optimization (bool)

        :return newvar: new float variable container (OptimizableVariable)
        """
        # to do: improve status so that it can also handle solves

        newvar = OptimizableVariable(name, value, status)

        self.listOfOptimizableVariables.append(newvar)

        return newvar
        
    def getAllOptimizableVariables(self):
        return self.listOfOptimizableVariables




















