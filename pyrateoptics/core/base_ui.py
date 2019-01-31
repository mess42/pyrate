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

from copy import copy

from .base import OptimizableVariable

# TODO: provide transformation dict in the form
# {"kind1": {"var_name": ("shown_varname", transform_func, inverse_transform_func, ...}, ...}
# where transform_func(x) transforms the value into the written value in the table
# and inverse_transform_func(x) transforms the value back
# e.g.: {"shape_conic": {"curv": ("radius", lambda x: 1./x, lambda x: 1./x)}}
# or: {"localcoordinates": {"tiltx": ("tiltx_deg", lambda x: x/degree, lambda x: x*degree)}}

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

        myvarlist = []

        def get_value_and_modstate(variable):
            value = float(variable())
            can_be_modified = variable.var_type == "fixed" or\
                variable.var_type == "variable"
            variable_triple = (variable.name, value, can_be_modified)
            myvarlist.append(variable_triple)
            return (value, can_be_modified)

        dict_to_ui = self.myclass.getDictionary()
        dict_to_ui.pop("variables")
        dict_to_ui.pop("classes")
        # the following statement can be used to obtain
        # a structural list of variable triples
        self.myclass.getTypesForDict(typ=OptimizableVariable,
                                     func=get_value_and_modstate)
        # although we only want to use a flat variable list
        # in a first step
        dict_to_ui["variables_list"] = myvarlist
        return dict_to_ui

    def modifyFromDictionary(self, dict_from_gui, override_unique_id=False):
        # make copy from dict to prevent modification
        dict_from_gui_copy = copy(dict_from_gui)
        # check later if protocol version is changed
        protocol_version = dict_from_gui_copy.pop("protocol_version")
        variables_list = dict_from_gui_copy.pop("variables_list")
        if not override_unique_id:
            dict_from_gui_copy.pop("unique_id")
        # we removed all variables which are not necessary
        # this may depend on protocol_number
        for (key_dict, value_dict) in dict_from_gui_copy.items():
            setattr(self.myclass, key_dict, value_dict)

        def set_value_from_modstate(variable):
            for variable_triple in variables_list:
                (variable_name, variable_value, can_be_modified) =\
                    variable_triple
                if variable_name == variable.name and can_be_modified:
                    variable.setvalue(variable_value)
        # update values of modifyable variables
        # this algorithm is of O(N^2) but since N is of order 10 this is not
        # critical
        self.myclass.getTypesForDict(typ=OptimizableVariable,
                                     func=set_value_from_modstate)
