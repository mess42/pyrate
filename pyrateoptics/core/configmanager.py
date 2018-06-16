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

from log import BaseLogger
import copy



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
            """            
            for index in range(length_value_tuples):
                print(index)
                new_instance = copy.deepcopy(self.base_instance)
                dictoptvariables = new_instance.getAllOptimizableVariables()    
                for (key, val_tuples) in dict_of_keys_and_value_tuples.items():
                    dictoptvariables[key].setvalue(val_tuples[index])
                instance_tuple.append(new_instance)                
            """
            
        return (copy.deepcopy(self.base_instance))

if __name__ == "__main__":

    s = OpticalSystem()
    m = MulticonfigManager(s)

    [s2, s3, s4] = m.setOptimizableVariables({"s.global.object.decz": (2.0, 3.0, 4.0)})
    listOptimizableVariables(s2)
