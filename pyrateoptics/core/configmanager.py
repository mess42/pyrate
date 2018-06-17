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
import re



class MulticonfigManager(BaseLogger):
    """
    Purpose is to make a deep copy of all variables which are subject to change
    and to perform a swallow copy of the variables which are shared between the
    different instances.
    """

    
    def __init__(self, base_instance=None, **kwargs):
        super(MulticonfigManager, self).__init__(**kwargs)
        self.base_instance = base_instance
        
    def setOptimizableVariables(self, dict_of_keys_and_value_tuples):
        """
        Set values in different instance of base instance.
        (deep copy). What about pickups?
        """
        length_value_tuples = 0
        if self.base_instance is None:
            instance_tuple = None
        else:
            instance_tuple = []
            if len(dict_of_keys_and_value_tuples) > 0:
                length_value_tuples = len(dict_of_keys_and_value_tuples.values()[0])
                # use "first" element in dict to obtain no. of copies
 
                for index in range(length_value_tuples):
                    new_instance = copy.copy(self.base_instance) # shallow copy

                    for (key, val_tuples) in dict_of_keys_and_value_tuples.iteritems():
                        variable = self.base_instance.getVariable(key)
                        new_variable = copy.deepcopy(variable)
                        new_variable.setvalue(val_tuples[index])
                        new_instance.resetVariable(key, new_variable)

                    instance_tuple.append(new_instance)                
            
        return instance_tuple

if __name__ == "__main__":
    
    s = OpticalSystem(name="s")
    s.lst = []
    s.lst.append({})
    s.lst.append({})
    s.lst[0]["a"] = OptimizableVariable(variable_type="fixed", name="v1", value=3.0)
    s.lst[1]["b"] = OptimizableVariable(variable_type="fixed", name="v2", value=7.0)

    s.rootcoordinatesystem.decz = OptimizableVariable(name="decz")
        
    m = MulticonfigManager(s)

    [s2, s3, s4] = m.setOptimizableVariables({"s.global.decz": (2.0, 3.0, 4.0)})
    for ss in (s2, s3, s4):    
        mydict = listOptimizableVariables(ss)
        var = mydict["vars"]["s.global.decz"]
        longkey = mydict["longkeystrings"]["s.global.decz"]
        print("sid: ", id(ss))
        print("varval: ", var())        
        print("varid: ", id(var))
        
