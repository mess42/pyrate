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

from .log import BaseLogger

'''
# TODO: functions for multi configs

    def resetVariable(self, key, var):
        """
        Resets variable by evaluating deref string from getAllVariables
        dictionary. Maybe this could be solved by using less black magic.
        Although this method is necessary for maintaining multi configs.
        """

        dict_of_vars = self.getAllVariables()
        deref = dict_of_vars["deref"][key]
        exec(compile("self" + deref + " = var", "<string>", "exec"))

    def getVariable(self, key):
        """
        Gets variable from short key.
        """
        dict_of_vars = self.getAllVariables()
        variable = dict_of_vars["vars"][key]
        return variable
'''


class ClassWithOptimizableVariables(BaseLogger):
    """
    Implementation of some class with optimizable variables with the help
    of a dictionary. This class is also able to collect the variables and
    their values from its subclasses per recursion.

    The goal is to provide an abstract class with pretty much functionality
    which gives the user the opportunity to implement a certain logic via
    class inheritance and interface the child class with some type of
    optimizer.
    """
    def __init__(self, name="", kind="classwithoptimizablevariables",
                 **kwargs):
        """
        Initialize with empty dict.
        """
        super(ClassWithOptimizableVariables, self).__init__(
                name=name,
                kind=kind,
                **kwargs)

        self.annotations = {}
        self.list_observers = []
        # for the optimizable variable class it is useful to have some observer
        # links they get informed if variables change their values

    def appendObservers(self, obslist):
        self.list_observers += obslist

    def informObservers(self):
        for obs in self.list_observers:
            obs.informAboutUpdate()
