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

from .base_ui_transform import (transformation_dictionary_to_ui,
                                transformation_dictionary_from_ui,
                                transformation_dictionary_ui_string)
from .iterators import (OptimizableVariableCollector,
                        SerializationIterator)
from .log import BaseLogger


class UIInterfaceClassWithOptimizableVariables(BaseLogger):
    """
    This class is intended to provide an easy to use interface to some
    GUI. It should provide a dict on query with the most important class
    data and it should accept a dict in the same format to be transfered to
    the class. Optimizable variables should be given as (value, bool) pair
    in a dict with names as keys where bool denotes whether those values
    can be changed. bool is read only. (True for variable status "fixed",
                                        and "variable", False for "pickup",
                                        and external)
    """
    def __init__(self, some_class_with_optimizable_variables,
                 name="", kind="uiinterface"):
        super(UIInterfaceClassWithOptimizableVariables,
              self).__init__(name=name,
                             kind=kind)

        self.myclass = some_class_with_optimizable_variables

    def queryForDictionary(self):

        dict_to_ui = SerializationIterator(self.myclass).dictionary

        myvarlist = [(ov.name, ov(),
                      ov.var_type() == "fixed" or ov.var_type() == "variable")
                     for ov in OptimizableVariableCollector(self.myclass).
                     variables_list]

        dict_to_ui["variables_list"] = myvarlist
        return dict_to_ui

    def modifyFromDictionary(self, dict_from_ui, override_unique_id=False):
        # make copy from dict to prevent modification
        dict_from_ui_copy = copy(dict_from_ui)
        # check later if protocol version is changed

        protocol_version = dict_from_ui_copy.pop("protocol_version")
        variables_list = dict_from_ui_copy.pop("variables_list")
        if not override_unique_id:
            dict_from_ui_copy.pop("unique_id")
        # we removed all variables which are not necessary
        # this may depend on protocol_number
        for (key_dict, value_dict) in dict_from_ui_copy.items():
            setattr(self.myclass, key_dict, value_dict)

        myvarlist = OptimizableVariableCollector(self.myclass).variables_list

        for ov in myvarlist:
            for (variable_name, variable_value, can_be_modified)\
             in variables_list:
                    if ov.name == variable_name and can_be_modified:
                        ov.set_value(variable_value)

        # update values of modifyable variables
        # this algorithm is of O(N^2) but since N is of order 10 this is not
        # critical

    def transformDictionaryForUI(
            self, dict_to_ui,
            transform_dictionary_value=transformation_dictionary_to_ui,
            transform_dictionary_strings=transformation_dictionary_ui_string):
        # prevent modification
        string_dict_to_ui = copy(dict_to_ui)
        string_dict_to_ui["variables_list"] = []
        transform_kind_to_use = transform_dictionary_value.get(dict_to_ui["kind"], None)

        for (var_name, var_value, var_modifiable) in dict_to_ui["variables_list"]:
            # if no transform is found just convert to string
            final_value = var_value
            final_var_name = var_name

            if transform_kind_to_use is not None:
                var_transform = transform_kind_to_use.get(var_name, None)
                if var_transform is not None:
                    # if transform is found and variable is found
                    # transform value and convert to string
                    (transformed_var_name,
                     transformed_value) = var_transform
                    final_value = transformed_value(var_value)
                    final_var_name = transformed_var_name

            string_final_value = str(final_value)
            string_transform = transform_dictionary_strings.get(
                    final_var_name, None)
            if string_transform is not None:
                string_final_value = string_transform(string_final_value)
            final_type = float if isinstance(final_value, (int, float)) else\
                type(final_value)

            string_dict_to_ui["variables_list"].append(
                    (final_var_name,
                     string_final_value,
                     var_modifiable,
                     final_type))

        return string_dict_to_ui

    def transformDictionaryFromUI(
            self, string_dict_from_ui,
            transform_dictionary_value=transformation_dictionary_from_ui,
            transform_dictionary_strings=transformation_dictionary_ui_string):
        dict_from_ui = copy(string_dict_from_ui)
        dict_from_ui["variables_list"] = []
        transform_kind_to_use = transform_dictionary_value.get(
                dict_from_ui["kind"], None)
        for (transformed_var_name,
             string_transformed_var_value,
             var_modifiable, var_type) in string_dict_from_ui["variables_list"]:

            string_transform = transform_dictionary_strings.get(
                    transformed_var_name, None)
            string_var_value = string_transformed_var_value

            # reverse string transform
            if string_transform is not None:
                string_var_value = string_transform(string_var_value)

            # restore variable value and variable name
            final_value = var_type(string_var_value)
            final_var_name = transformed_var_name

            if transform_kind_to_use is not None:
                var_transform =\
                    transform_kind_to_use.get(transformed_var_name, None)
                if var_transform is not None:
                    # if transform is found and variable is found
                    # transform value and convert to string
                    (final_var_name,
                     invtransformed_value) = var_transform
                    final_value = invtransformed_value(final_value)

            dict_from_ui["variables_list"].append(
                    (final_var_name, final_value, var_modifiable))
        return dict_from_ui

    def queryUI(self):
        dict_to_ui = self.queryForDictionary()
        transformed_string_dict_to_ui =\
            self.transformDictionaryForUI(dict_to_ui)
        return transformed_string_dict_to_ui

    def modifyUI(self, string_dict_from_ui):
        invtransformed_dict_from_ui =\
            self.transformDictionaryFromUI(string_dict_from_ui)
        self.modifyFromDictionary(invtransformed_dict_from_ui)
