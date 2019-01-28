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

from .base import OptimizableVariable


class UIInterfaceClassWithOptimizableVariables:  # maybe derive from BaseLogger
    """
    This class is intended to provide an easy to use interface to some
    GUI. It should provide a dict on query with the most important class
    data an it should accept a dict in the same format to be transfered to
    the class. Optimizable variables should be given as (value, bool) pair
    in a dict with names as keys where bool denotes whether those values
    can be changed. bool is read only. (True for variable status "fixed",
                                        and "variable", False for "pickup",
                                        and external)
    """
    def __init__(self, some_class_with_optimizable_variables):
        self.myclass = some_class_with_optimizable_variables

    def queryForDictionary(self):

        def get_value_and_modstate(variable):
            return (float(variable()), variable.var_type == "fixed" or
                    variable.var_type == "variable")

        dict_to_ui = self.myclass.getDictionary()
        dict_to_ui.pop("variables")
        dict_to_ui.pop("classes")
        dict_to_ui["variables"] =\
            self.myclass.getTypesForDict(typ=OptimizableVariable,
                                         func=get_value_and_modstate)
        # maybe some flat list is sufficient due to unique name within
        # the class (structural reconstruction is not necessary)
        # format: [("name", value, bool), ("name2", value2, bool), ...]
        return dict_to_ui

    def modifyFromDictionary(self, dict_from_gui):
        pass
