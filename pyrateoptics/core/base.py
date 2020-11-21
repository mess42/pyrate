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

from .log import BaseLogger


class ClassWithOptimizableVariables(BaseLogger):
    """
    Implementation of some class with optimizable variables with the help
    of a dictionary. This class is also able to collect the variables and
    their values from its subclasses per recursion.

    The goal is to provide an abstract class with pretty much functionality
    which gives the user the opportunity to implement a certain logic via
    class inheritance and interface the child class with some type of
    optimizer.
    """
    def __init__(self, annotations_dict=None,
                 structure_dict=None, name="",
                 serializationfilter=None):
        """
        Initialize with empty dict.

        :param: annotations_dict ... contains data which has to be provided
        :param: structure_dict ... contains structural elements (variables,
                                                                 subclasses)

        """
        super(ClassWithOptimizableVariables, self).__init__(name=name)
        if annotations_dict is None:
            annotations_dict = {}
        if structure_dict is None:
            structure_dict = {}

        self.annotations = annotations_dict
        self.serializationfilter =\
            [] if serializationfilter is None else serializationfilter

        for (key, value) in structure_dict.items():
            if not hasattr(self, key):
                self.__setattr__(key, value)
        self.initialize_from_annotations()

    @classmethod
    def p(cls, **kwargs):
        """
        Class initialization via parameters and filling up
        annotations and structure dictionary. To be implemented
        by every child class.
        """
        name = kwargs.get("name", "")
        return cls({}, {}, name)

    def initialize_from_annotations(self):
        """
        Further initialization stages from annotations which need to be
        done to get a valid object.
        """
        pass

    def setKind(self):
        """
        Set kind of object. Should be overridden in child classes.
        """
        self.kind = "classwithoptimizablevariables"
