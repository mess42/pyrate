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

import copy
import logging


from pyrateoptics import listOptimizableVariables
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.core.base import ClassWithOptimizableVariables
from pyrateoptics.core.iterators import OptimizableVariableKeyIterator
from pyrateoptics.core.optimizable_variable import (FloatOptimizableVariable,
                                                    OptimizableVariable,
                                                    FixedState, VariableState,
                                                    PickupState)

from pyrateoptics.core.log import BaseLogger


class ConfigContainer(ClassWithOptimizableVariables):
    """
    Contains all instances from config manager to be compatible with the optimizer
    frontend.
    """

    def __init__(self, instance_list=None, name=""):
        super(ConfigContainer, self).__init__(name=name)
        self.instance_list = instance_list

    def setKind(self):
        self.kind = "configcontainer"

    # TODO: to be tested


class ConfigManager(BaseLogger):
    """
    Purpose is to make a deep copy of all variables which are subject to change
    and to perform a swallow copy of the variables which are shared between the
    different instances.
    """

    def __init__(self, base_instance=None, name=""):
        super(ConfigManager, self).__init__(name=name)
        self.base_instance = base_instance

    def setKind(self):
        self.kind = "configmanager"

    def set_optimizable_variables(self, names_tuple,
                                  dict_of_keys_and_value_tuples):
        """
        Set values in different instance of base instance.
        (deep copy). What about pickups?

        @param names_tuple (tuple of strings): names of the new systems
                                                (changes also keys)
        @param dict_of_keys_and_tuples (dict of tuples): consists of values
                                                        for fixed or
                                                        variable optimizable
                                                        variables.

        @return instance_list (list of instances)

        """
        length_value_tuples = 0
        if self.base_instance is None:
            instance_list = None
        else:
            instance_list = []
            if names_tuple is None:
                names_tuple = tuple([self.base_instance.name
                                     for r in dict_of_keys_and_value_tuples.
                                     values()[0]])
            self.info("Constructing multiple configs with names: %s" %
                      (str(names_tuple)))
            if len(names_tuple) > 0:
                length_value_tuples = len(names_tuple)
                # use names_tuple to obtain no. of copies
                # notice that first flat copying instances and only deep copying
                # of variables does not work. Therefore we did it the other way
                # around.

                self.debug("Deep copying base instance")
                for index in range(length_value_tuples):
                    self.debug(str(index))
                    instance_list.append(copy.deepcopy(self.base_instance))
                    # deep copy instance

                self.debug("Constructing variables")
                for (index, instance) in enumerate(instance_list):
                    self.debug(str(index))
                    # reset all non-changing variables to instances of base_instance
                    # components.
                    for (key, variable) in OptimizableVariableKeyIterator(
                            instance).variables_dictionary.items():
                        print(key)
                        if key in dict_of_keys_and_value_tuples:
                            val_tuple = dict_of_keys_and_value_tuples[key]
                            self.debug("%s to %s" % (key, val_tuple[index]))
                            #variable.setvalue(val_tuple[index])
                            if not isinstance(val_tuple, tuple):
                                self.warning("Incorrect format for multi config values!")
                                self.debug("Values must be of type (string, contents):")
                                self.debug("For fixed values: (\"fixed\", fixed_val)")
                                self.debug("For pickup values: (\"pickup\", function)")
                            else:
                                (mcv_type, mcv_contents) = val_tuple[index]
                                if mcv_type.lower() == "fixed":
                                    instance.resetVariable(
                                        key,
                                        OptimizableVariable(
                                            FixedState(mcv_contents),
                                            name=variable.name))
                                elif mcv_type.lower() == "pickup":
                                    instance.resetVariable(
                                        key,
                                        OptimizableVariable(
                                            PickupState(
                                                mcv_contents,
                                                variable),
                                            name=variable.name))
                                else:
                                    self.warning("Unknown type for multi config values")
                        else:
                            self.debug("Reseting %s" % (key,))
                            instance.resetVariable(key, self.base_instance.getVariable(key))
                    # set names of the new instances accordingly
                    self.debug("Setting name of instance %s" % (names_tuple[index]))
                    instance.set_name(names_tuple[index])

        return instance_list

# TODO: functions for multi configs
#
#    def resetVariable(self, key, var):
#        """
#        Resets variable by evaluating deref string from getAllVariables
#        dictionary. Maybe this could be solved by using less black magic.
#        Although this method is necessary for maintaining multi configs.
#        """
#
#        dict_of_vars = self.getAllVariables()
#        deref = dict_of_vars["deref"][key]
#        exec(compile("self" + deref + " = var", "<string>", "exec"))
#
#    def getVariable(self, key):
#        """
#        Gets variable from short key.
#        """
#        dict_of_vars = self.getAllVariables()
#        variable = dict_of_vars["vars"][key]
#        return variable
#


if __name__ == "__main__":

    def main():
        "Main code for demo purposes"
        logging.basicConfig(level=logging.DEBUG)

        mysys = OpticalSystem(name="s")
        mysys.lst = []
        mysys.lst.append({})
        mysys.lst.append({})
        mysys.lst[0]["a"] = FloatOptimizableVariable(
            FixedState(3.0), name="v1")
        mysys.lst[1]["b"] = FloatOptimizableVariable(
            VariableState(7.0), name="v2")

        mysys.rootcoordinatesystem.decz = FloatOptimizableVariable(
            FixedState(-99.0),
            name="decz")

        listOptimizableVariables(mysys)

        confmanager = ConfigManager(mysys, name="mc")

        [mysys2, mysys3] = confmanager.set_optimizable_variables(
            ("s2", "s3"),
            {"s.global.decz": (("pickup", lambda x: x + 2.0),
                               ("pickup", lambda x: x + 3.0)),
             "s.global.decy": (("fixed", -2.), ("fixed", -3.))})
        mysys.rootcoordinatesystem.decx.setvalue(-98.0)
        for syscopy in (mysys2, mysys3):
            mydict = listOptimizableVariables(syscopy)
            print(mydict)

main()
