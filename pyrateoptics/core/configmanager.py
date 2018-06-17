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


# TODO: YAML for serialization
# TODO: Can ZMX converted into YAML easily?
# TODO: serialization via JSONEncode derived class decorator (thomas annotator concept)

from pyrateoptics import listOptimizableVariables
from pyrateoptics.raytracer.optical_system import OpticalSystem
from pyrateoptics.core.base import OptimizableVariable

from log import BaseLogger
import copy
import logging


class ConfigManager(BaseLogger):
    """
    Purpose is to make a deep copy of all variables which are subject to change
    and to perform a swallow copy of the variables which are shared between the
    different instances.
    """

    
    def __init__(self, base_instance=None, **kwargs):
        super(ConfigManager, self).__init__(**kwargs)
        self.base_instance = base_instance
        
    def setOptimizableVariables(self, names_tuple, dict_of_keys_and_value_tuples):
        """
        Set values in different instance of base instance.
        (deep copy). What about pickups?
        
        @param names_tuple (tuple of strings): names of the new systems (changes also keys)
        @param dict_of_keys_and_tuples (dict of tuples): consists of values for fixed or
                variable optimizable variables.
        
        """
        length_value_tuples = 0
        if self.base_instance is None:
            instance_list = None
        else:
            instance_list = []
            if names_tuple is None:
                names_tuple = tuple([self.base_instance.name for r in dict_of_keys_and_value_tuples.values()[0]])             
            self.info("Constructing multiple configs with names: %s" % (str(names_tuple)))
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
                    for (key, variable) in instance.getAllVariables()["vars"].iteritems():
                        
                        if key in dict_of_keys_and_value_tuples:
                            val_tuple = dict_of_keys_and_value_tuples[key]
                            self.debug("%s to %s" % (key, val_tuple[index]))
                            #variable.setvalue(val_tuple[index])
                            # TODO: how to get pickup functions with multiple arguments?
                            # TODO: how to improve this code?                            
                            if type(val_tuple) is not tuple:
                                self.warning("Incorrect format for multi config values!")
                                self.debug("Values must be of type (string, contents):")
                                self.debug("For fixed values: (\"fixed\", fixed_val)")
                                self.debug("For pickup values: (\"pickup\", function)")
                            else:
                                (mcv_type, mcv_contents) = val_tuple[index]
                                if mcv_type.lower() == "fixed":
                                    instance.resetVariable(key, OptimizableVariable(variable_type="fixed", value=mcv_contents, name=variable.name))
                                elif mcv_type.lower() == "pickup":
                                    instance.reserVariable(key, OptimizableVariable(variable_type="pickup", function=mcv_contents, args=(variable,), name=variable.name))
                        else:
                            self.debug("Reseting %s" % (key,))
                            instance.resetVariable(key, self.base_instance.getVariable(key))
                    # set names of the new instances accordingly
                    self.debug("Setting name of instance %s" % (names_tuple[index]))
                    instance.setName(names_tuple[index])
                
        return instance_list

if __name__ == "__main__":
    
    logging.basicConfig(level=logging.DEBUG)
    
    s = OpticalSystem(name="s")
    s.lst = []
    s.lst.append({})
    s.lst.append({})
    s.lst[0]["a"] = OptimizableVariable(variable_type="fixed", name="v1", value=3.0)
    s.lst[1]["b"] = OptimizableVariable(variable_type="variable", name="v2", value=7.0)

    s.rootcoordinatesystem.decz = OptimizableVariable(name="decz", value=-99.0)

    listOptimizableVariables(s)
        
    m = ConfigManager(s, name="mc")

    [s2, s3, s4] = m.setOptimizableVariables(("s2", "s3", "s4"), 
                {"s.global.decz": (("fixed", 2.0), ("fixed", 3.0), ("fixed", 4.0)), 
                "s.global.decy": (("fixed", -2.), ("fixed", -3.), ("fixed", -4.))})
    s.rootcoordinatesystem.decx.setvalue(-98.0)    
    for ss in (s2, s3, s4):    
        mydict = listOptimizableVariables(ss)
        
