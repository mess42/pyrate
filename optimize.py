#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and Johannes Hartung j.hartung@gmx.net
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
from scipy.optimize import minimize
import warnings



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

        #print "create opt var: ",name, " val ", value, " status ", status

        newvar = OptimizableVariable(name, value, status)

        self.listOfOptimizableVariables.append(newvar)

        return newvar

    def getAllOptimizableVariables(self):
        #print "getAllOptVariables: ", [i.name for i in self.listOfOptimizableVariables]
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

        #print "setstatus listofvars: ", [i.name + " " + str(i.val) + " " + str(i.status)  for i in self.listOfOptimizableVariables]

    def copyOptimizableVariables(self, otherClassWithOptVars):
        """ Helper function due to non optimal communication of optvars between surface, material and shape
            should be removed later """

        try:
            varsToRemove = otherClassWithOptVars.getAllOptimizableVariables()
            #print "other listofoptvars: ", [i.name for i in varsToRemove]

            for v in varsToRemove:
                self.listOfOptimizableVariables.remove(v)
        except:
            pass

        #print "orig listofoptvars: ", [i.name for i in self.listOfOptimizableVariables]


        # add optimizable variables other object
        self.listOfOptimizableVariables += otherClassWithOptVars.getAllOptimizableVariables()

        #print "new listofoptvars: ", [i.name for i in self.listOfOptimizableVariables]

def optimizeSciPyInterface(s, meritfunction, **kwargs):
    # TODO: needs heavy testing
    """
    Scipy interface to minimize

    :param s: initial optical system
    :param meritfunction: pointer to merit function
    :param kwargs: keyword arguments to minimize from scipy

    :return s: optimized optical system object
    """

    optVars = s.getActiveOptimizableVariables()
    x0 = np.array([ov.val for ov in optVars])

    print "Starting values: ", x0
    print "Options to optimizer: ", kwargs['options']

    def MeritFunctionWrapper(x, *args):
        s2 = args[0]
        ovars = s2.getActiveOptimizableVariables()
        # ask optical reference system for his optimizable variables during every evaluation
        # TODO: wrapper on optical system level?

        for index, ov in enumerate(ovars):
            ov.val = x[index]
        return meritfunction(s2)

    res = minimize(MeritFunctionWrapper, x0, args=(s,), method=kwargs['method'], options=kwargs['options'])

    print "result: ", res

    for index, ov in enumerate(optVars):
        ov.val = res.x[index]


    return s

def optimizeNewton1D(s, meritfunction, iterations=1, dx=1e-6):
    """
    1 dimensional Newton optimizer for an optical system.

    :param s: initial OpticalSystem object, will be overwritten
    :param meritfunction: pointer on the merit function
    :param iterations: number of iterations (int)
    :param dx: infinitessimal change (float) for derivative

    :return s: optimized OpticalSystem object
    """

    warnings.filterwarnings("error", category=RuntimeWarning)
    # is needed such that we can intercept a RuntimeWarning during optimization procedure as an exception

    if dx <= 0:
        dx = 1e-6

    optVars = s.getActiveOptimizableVariables()

    print [i.name + ": " + str(i.val) + "; " + str(i.status) for i in optVars]

    for i in np.arange(iterations):

        print "Iteration: ", i, " merit: ", meritfunction(s)
        #print [i.name + ": " + str(i.val) + "; " + str(i.status) for i in optVars]


        for v in np.arange(len(optVars)):

            retry = True
            retrycount = 0

            while retry:

                try: # try block for intercepting RuntimeWarnings during evaluation of merit function
                    retrycount += 1
                    merit0 = meritfunction(s)
                    var = optVars[v]
                    varvalue0 = var.val
                    varvalue1 = var.val + dx*(1.0 - 2.0*np.random.random()) # also allow values in the negative direction
                    var.val = varvalue1

                    merit1 = meritfunction(s)



                    if (abs(merit1 - merit0) > 1e-16):
                        m = (merit1 - merit0) / (varvalue1 - varvalue0)
                        n = merit0 - m * varvalue0
                        varvalue2 = - n / m  # Newton method for next iteration value
                        print "m: ", m, " n: ", n, " v0: ", varvalue0, " v1: ", varvalue1, " newv: ", varvalue2
                    else:
                        varvalue2 = varvalue1

                    var.val = varvalue2
                    merit2 = meritfunction(s)

                    guard = 0  # guard element to prevent freezing
                    while merit2 > merit0 and varvalue2 - varvalue0 > dx and guard < 1000:
                        varvalue2 = 0.5 * (varvalue2 + varvalue0)
                        var.val = varvalue2
                        merit2 = meritfunction(s)
                        guard += 1

                    retry = False

                except RuntimeWarning:
                    print "Runtime waring occured retrying ..."
                    var.val = varvalue0 # reseting system variable
                    if retrycount < 10: # if it does not work after 10 trys than set variable to standard value and end loop
                        retry = True
                    else:
                        retry = False




    warnings.filterwarnings("default") # reset warning behavior to default value


    return s
