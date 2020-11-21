#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
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
                 name="", updatefunction=None):
        super(Optimizer, self).__init__(name=name)

        def noupdate(_):
            "No update dummy callback function"
            pass

        self.collector = OptimizableVariableActiveCollector(
            classwithoptvariables)
        self.meritfunction = meritfunction  # function to minimize
        if updatefunction is None:
            updatefunction = noupdate
        self.updatefunction = updatefunction
        # function to be called to update classwithoptvariables
        self.set_backend(backend)
        # eats vector performs optimization, returns vector
        # scipy Nelder-Mead, scipy ..., evolutionary, genetic, ...

        self.meritparameters = {}
        self.updateparameters = {}
        self.number_of_calls = 0
        # how often is the merit function called during one run?

    def setKind(self):
        self.kind = "optimizer"

    def set_backend(self, backend):
        "Setter for backend."
        self.__backend = backend
        self.__backend.init(self.meritfunction_wrapper)

    backend = property(fget=None, fset=set_backend)

    def meritfunction_wrapper(self, vecx):
        """
        Merit function wrapper for backend.
        Notice that vecx and length of active values must have the same size

        :param vecx (np.array): active variable values
        :param meritfunction (function): meritfunction depending on s

        :return value of the merit function
        """
        self.number_of_calls += 1
        self.collector.fromNumpyArrayTransformed(vecx)
        self.updatefunction(self.collector.class_instance,
                            **self.updateparameters)
        res = self.meritfunction(self.collector.class_instance,
                                 **self.meritparameters)
        self.debug("call number " + str(self.number_of_calls) +
                   " meritfunction: " + str(res))
        return res

    def run(self):
        '''
        Funtion to perform a certain number of optimization steps.
        '''
        self.info("optimizer run start")
        vecx0 = self.collector.toNumpyArrayTransformed()

        self.info("initial x: " + str(vecx0))
        self.info("initial merit: " + str(self.meritfunction_wrapper(vecx0)))
        self.debug("calling backend run")
        xfinal = self.__backend.run(vecx0)
        self.debug("finished backend run")
        self.info("final x: " + str(xfinal))
        self.info("final merit: " + str(self.meritfunction_wrapper(xfinal)))
        self.collector.fromNumpyArrayTransformed(xfinal)
        self.info("called merit function " + str(self.number_of_calls) +
                  " times.")
        self.number_of_calls = 0
        self.info("optimizer run finished")
        return self.collector.class_instance
