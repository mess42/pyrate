#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de
               and    Thomas Heinze t.heinze@uni-jena.de
               and    others

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

from ..core.log import BaseLogger
from ..core.iterators import OptimizableVariableActiveCollector


class Optimizer(BaseLogger):
    '''
    Easy optimization interface. All variables are public such that a quick
    attachment of other meritfunctions or other update functions with other
    parameters is possible.
    '''
    def __init__(self, classwithoptvariables,
                 meritfunction, backend,
                 name="", kind="optimizer", updatefunction=None):

        def noupdate(cl):
            pass

        super(Optimizer, self).__init__(name=name, kind=kind)
        self.collector = OptimizableVariableActiveCollector(
                classwithoptvariables)
        self.meritfunction = meritfunction  # function to minimize
        if updatefunction is None:
            updatefunction = noupdate
        self.updatefunction = updatefunction  # function to be called to update classwithoptvariables
        self.setBackend(backend)  # eats vector performs optimization, returns vector
        # scipy Nelder-Mead, scipy ..., evolutionary, genetic, ...

        self.meritparameters = {}
        self.updateparameters = {}
        self.number_of_calls = 0 # how often is the merit function called during one run?

    def setBackend(self, backend):
        self.__backend = backend
        self.__backend.init(self.MeritFunctionWrapper)

    backend = property(fget=None, fset=setBackend)

    def MeritFunctionWrapper(self, x):
        """
        Merit function wrapper for backend.
        Notice that x and length of active values must have the same size

        :param x (np.array): active variable values
        :param meritfunction (function): meritfunction depending on s

        :return value of the merit function
        """
        self.number_of_calls += 1
        self.collector.fromNumpyArrayTransformed(x)
        self.updatefunction(self.collector.class_instance,
                            **self.updateparameters)
        res = self.meritfunction(self.collector.class_instance,
                                 **self.meritparameters)
        self.debug("call number " + str(self.number_of_calls) + " meritfunction: " + str(res))
        return res

    def run(self):
        '''
        Funtion to perform a certain number of optimization steps.
        '''
        self.info("optimizer run start")
        x0 = self.collector.toNumpyArrayTransformed()

        self.info("initial x: " + str(x0))
        self.info("initial merit: " + str(self.MeritFunctionWrapper(x0)))
        self.debug("calling backend run")
        xfinal = self.__backend.run(x0)
        self.debug("finished backend run")
        self.info("final x: " + str(xfinal))
        self.info("final merit: " + str(self.MeritFunctionWrapper(xfinal)))
        self.collector.fromNumpyArrayTransformed(xfinal)
        self.info("called merit function " + str(self.number_of_calls) +
                  " times.")
        self.number_of_calls = 0
        self.info("optimizer run finished")
        return self.collector.class_instance
