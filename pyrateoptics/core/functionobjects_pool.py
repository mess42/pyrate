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
from .functionobject import FunctionObject


class FunctionObjectsPool(BaseLogger):
    """
    For serializing a collection of function objects
    """

    def __init__(self, functionobjectsdictionary, name=""):
        super(FunctionObjectsPool, self).__init__(name=name)

        self.functionobjects_dictionary = functionobjectsdictionary

    def setKind(self):
        self.kind = "functionobjectspool"

    def toDictionary(self):
        result_dictionary = {}
        for (key, fo) in self.functionobjects_dictionary.items():
            result_dictionary[key] = fo.to_dictionary()
        return result_dictionary

    @staticmethod
    def fromDictionary(dictionary, source_checked, variables_checked, name=""):
        result_dictionary = {}
        for (key, fodict) in dictionary.items():
            result_dictionary[key] = FunctionObject.from_dictionary(
                fodict, source_checked, variables_checked)
        return FunctionObjectsPool(result_dictionary, name=name)
