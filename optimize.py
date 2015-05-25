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

import numpy as np


class OptimizableVariable(object):
    """
    Class that contains a float. Used to get a pointer on a variable.
    """
    def __init__(self, name, value=0.0, status=False):
        self.name = name
        self.val = value
        self.status = status


class ClassWithOptimizableVariables(object):
    """
    Virtual class for all classes that contain an optimizable variable.
    """
    def __init__(self):
        self.listOfOptimizableVariables = []

    def createOptimizableVariable(self, name, value=0.0, status=False):
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

    def getAllOptimizableValues(self):
        return [a.val for a in self.getAllOptimizableVariables()]

    def getAllOptimiziableStates(self):
        return [a.status for a in self.getAllOptimizableVariables()]

    def getActiveOptimizableVariables(self):
        ind = np.array(self.getAllOptimiziableStates())
        return np.array(self.getAllOptimizableVariables())[ind]

    def getActiveOptimizableValues(self):
        ind = np.array(self.getAllOptimiziableStates())
        return np.array(self.getAllOptimizableValues())[ind]

    def setStatus(self, name, status=True):
        N = len(self.listOfOptimizableVariables)

        names = []
        for i in np.arange(N):
            names.append(self.listOfOptimizableVariables[i].name)

        if name in names:
            i = names.index(name)
            self.listOfOptimizableVariables[i].status = status
        else:
            print "not found"


def optimizeNewton1D(s, meritfunction, iterations=1, dxFactor=1.00001):
    """
    1 dimensional Newton optimizer for an optical system.

    :param s: initial OpticalSystem object, will be overwritten
    :param meritfunction: pointer on the merit function
    :param iterations: number of iterations (int)
    :param dxFactor: factor determining an infinitessimal change (float)

    :return s: optimized OpticalSystem object
    """

    if dxFactor <= 1:
        print "Warning: dxFactor must be larger than 1. Setting value to 1.00001"
        dxFactor = 1.00001

    optVars = s.getActiveOptimizableVariables()

    for i in np.arange(iterations):
        for v in np.arange(len(optVars)):
            merit0 = meritfunction(s)
            var = optVars[v]
            varvalue0 = var.val
            varvalue1 = var.val * dxFactor
            var.val = varvalue1
            merit1 = meritfunction(s)

            m = (merit1 - merit0) / (varvalue1 - varvalue0)
            n = merit0 - m * varvalue0
            varvalue2 = - n / m  # Newton method for next iteration value

            var.val = varvalue2
            merit2 = meritfunction(s)

            guard = 0  # guard element to prevent freezing
            while merit2 > merit0 and varvalue2 / varvalue0 > dxFactor and guard < 1000:
                varvalue2 = 0.5 * (varvalue2 + varvalue0)
                var.val = varvalue2
                merit2 = meritfunction(s)
                guard += 1
    return s
