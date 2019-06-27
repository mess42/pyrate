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

from .log import BaseLogger
from .iterators import SerializationIterator
from .optimizable_variables_pool import OptimizableVariablesPool


class Serializer(BaseLogger):
    def __init__(self, class_instance, name=""):
        super(Serializer, self).__init__(name=name, kind="serializer")
        self.class_instance = class_instance
        self.serialize()

    def serialize(self):
        serialization = SerializationIterator(self.class_instance,
                                              remove=["annotations",
                                                      "list_observers"])
        serialization.collectStructure(remove=["annotations",
                                               "list_observers"])
        self.serialization = [
                serialization.dictionary,
                dict([(k,
                       SerializationIterator(v,
                                             remove=["annotations",
                                                     "list_observers"]
                                             ).dictionary
                       ) for (k, v) in serialization.classes_dictionary.items()]
                     ),
                OptimizableVariablesPool(serialization.variables_dictionary).toDictionary()
        ]
