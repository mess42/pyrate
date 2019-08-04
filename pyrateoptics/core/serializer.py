#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

import uuid

from .log import BaseLogger
from .iterators import SerializationIterator
from .optimizable_variables_pool import OptimizableVariablesPool


class Serializer(BaseLogger):
    def __init__(self, class_instance, name=""):
        super(Serializer, self).__init__(name=name)
        self.class_instance = class_instance
        self.serialize()

    def setKind(self):
        self.kind = "serializer"

    def serialize(self):
        serialization = SerializationIterator(self.class_instance,
                                              remove=["annotations",
                                                      "list_observers"])
        serialization.collectStructure(remove=["annotations",
                                               "list_observers"])
        optimizable_variables_pool = OptimizableVariablesPool(
            serialization.variables_dictionary)
        functionobjects_pool = optimizable_variables_pool.generateFunctionObjectsPool()

        self.serialization = [
                serialization.dictionary,
                dict([(k,
                       SerializationIterator(v,
                                             remove=["annotations",
                                                     "list_observers"]
                                             ).dictionary
                       ) for (k, v) in serialization.classes_dictionary.items()]
                     ),
                optimizable_variables_pool.toDictionary(),
                functionobjects_pool.toDictionary()
        ]


class Deserializer(BaseLogger):
    def __init__(self, serialization_list, name=""):
        super(Deserializer, self).__init__(name=name)
        self.serialization_list = serialization_list
        self.deserialize()

    def isUUID(self, uuidstr):
        if isinstance(uuidstr, str):
            try:
                uuid.UUID(uuidstr)
            except ValueError:
                return False
            else:
                return True
        else:
            return False


    def deserialize(self):
        (class_to_be_deserialized, subclasses_dict,
         optimizable_variables_pool_dict, functionobjects_pool_dict) =\
             self.serialization_list

        print(class_to_be_deserialized)
        print(subclasses_dict)
        print(optimizable_variables_pool_dict)
        print(functionobjects_pool_dict)

        # reconstruct function objects and optimizable variables
        # go through structure of subclasses
        # search for uuids in structure; if they are optimizable variables,
        # set them; if not let them alone
        # if structure dict does not contain any uuids anymore construct class
        # (DUCKTYPING!)

        self.class_instance = None


