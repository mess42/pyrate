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
    def __init__(self, name=""):
        """
        Initialize with empty dict.
        """
        super(ClassWithOptimizableVariables, self).__init__(name=name)

        self.annotations = {}

    def setKind(self):
        self.kind = "classwithoptimizablevariables"

    @classmethod
    def cast(cls, some_classwithoptvariables):
        """Cast an ClassWithOptimizableVariables into a child class.

        According to:
            https://stackoverflow.com/questions/15404256/changing-the-class-of-a-python-object-casting

        This is some kind of workaround and we are not satisfied about it.
        A better variant would be to remove the initialization in the
        child classes from the constructor such that a constructor creates
        an empty, but valid type. And the incode creation is done by some
        ....create() class function
        """
        assert isinstance(some_classwithoptvariables,
                          ClassWithOptimizableVariables)
        some_classwithoptvariables.__class__ = cls  # now type cast
        assert isinstance(some_classwithoptvariables, cls)
        return some_classwithoptvariables
