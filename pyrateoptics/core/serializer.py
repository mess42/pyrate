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
from pprint import pformat

from .log import BaseLogger
from .iterators import SerializationIterator
from .optimizable_variables_pool import OptimizableVariablesPool
from .base import ClassWithOptimizableVariables

from ..raytracer.localcoordinates import LocalCoordinates
from ..raytracer.surface import Surface
from ..raytracer.surface_shape import ZernikeFringe, Asphere, Conic
from ..raytracer.optical_element import OpticalElement
from ..raytracer.optical_system import OpticalSystem
from ..raytracer.material.material_glasscat import CatalogMaterial
from ..raytracer.material.material_isotropic import ConstantIndexGlass


class Serializer(BaseLogger):
    def __init__(self, class_instance, name=""):
        super(Serializer, self).__init__(name=name)
        self.class_instance = class_instance
        self.serialize()

    def setKind(self):
        self.kind = "serializer"

    def serialize(self):
        default_to_be_removed = ["annotations",
                                 "list_observers",
                                 "serializationfilter"]
        serialization = SerializationIterator(self.class_instance,
                                              remove=default_to_be_removed)
        serialization.collectStructure(remove=default_to_be_removed)
        optimizable_variables_pool = OptimizableVariablesPool(
            serialization.variables_dictionary)
        functionobjects_pool = optimizable_variables_pool.generateFunctionObjectsPool()

        self.serialization = [
                serialization.dictionary,
                dict([(k,
                       SerializationIterator(v,
                                             remove=default_to_be_removed
                                             ).dictionary
                       ) for (k, v) in serialization.classes_dictionary.items()]
                     ),
                optimizable_variables_pool.toDictionary(),
                functionobjects_pool.toDictionary()
        ]


class Deserializer(BaseLogger):
    def __init__(self, serialization_list,
                 source_checked, variables_checked, name="",
                 register_classes=None):
        super(Deserializer, self).__init__(name=name)
        self.serialization_list = serialization_list
        self.classes_dictionary = {
            "shape_ZernikeFringe": ZernikeFringe,
            "shape_Conic": Conic,
            "shape_Asphere": Asphere,
            "localcoordinates": LocalCoordinates,
            "constantindexglass": ConstantIndexGlass,
            "opticalelement": OpticalElement,
            "opticalsystem": OpticalSystem,
            "surface": Surface,
            "material_from_catalog": CatalogMaterial,
            }
        if register_classes is not None:
            for (kind_name, class_name) in register_classes:
                self.classes_dictionary[kind_name] = class_name

        self.deserialize(source_checked, variables_checked)

    def isUUID(self, uuidstr):
        """
        Check whether a variable is a valid uuid. This could be improved,
        once in the classes appear regular uuids which have nothing to do
        with the serialization process.
        """
        if isinstance(uuidstr, str):
            try:
                uuid.UUID(uuidstr)
            except ValueError:
                return False
            else:
                return True
        else:
            return False

    def deserialize(self, source_checked, variables_checked):
        """
        Convert the list obtained from a file or another source back
        into a class with optimizable variables, via a recursive
        reconstruction of subclasses. Source checked and variables checked
        are to be set to True by the user.
        """
        (class_to_be_deserialized, subclasses_dict,
         optimizable_variables_pool_dict, functionobjects_pool_dict) =\
            self.serialization_list

        self.debug("Deserializing class")

        def reconstruct_class(class_to_be_reconstructed, subclasses_dict,
                              reconstructed_variables_dict):
            """
            Main reconstruction function. Is to be called recursively in
            the tree. Reconstructs first in-class variables by accessing
            the reconstructed_variables_dict. In a second step, it
            reconstructs the whole class by setting attributes from the
            structure dictionary, sets the annotations from the annotations
            dictionary and finally constructs the ClassWithOptimizableVariables
            object.
            """

            def is_structure_free_of_uuids(structure_dict):
                """
                Checks whether structure dict is free of uuids
                If this is the case the reconstruction process is finished.
                Else there are either unreconstructed variables or classes
                in this structure dict.
                """

                def free_of_uuid(var):
                    """
                    This function is called recursively to verify that a
                    given nested structure is free of UUIDs as checked by
                    isUUID.
                    """
                    if self.isUUID(var):
                        return False
                    elif isinstance(var, list):
                        return all([free_of_uuid(v) for v in var])
                    elif isinstance(var, dict):
                        return all([free_of_uuid(v) for v in var.values()])
                    else:
                        return True
                return free_of_uuid(structure_dict)

            def reconstruct_variables(structure_dict,
                                      reconstructed_variables_dict):
                """
                This function is called to reconstruct variables from a pool
                with in a necessary sub class of a class which is to be
                reconstructed. It uses structure_dict of the class to be
                reconstructed and returns a modified version where the
                variables UUIDs are substituted by the appropriate objects.
                """

                def reconstructRecursively(variable,
                                           reconstructed_variables_dict):
                    if self.isUUID(variable):
                        if variable in reconstructed_variables_dict:
                            return reconstructed_variables_dict[variable]
                        else:
                            return variable
                    elif isinstance(variable, list):
                        return [reconstructRecursively(part,
                                                       reconstructed_variables_dict)
                                for part in variable]
                    elif isinstance(variable, dict):
                        return dict([(key, reconstructRecursively(part,
                                                                  reconstructed_variables_dict))
                                    for (key, part) in variable.items()])
                    else:
                        return variable

                new_structure_dict = reconstructRecursively(
                        structure_dict,
                        reconstructed_variables_dict)

                return new_structure_dict

            def reconstruct_subclasses(structure_dict, subclasses_dict,
                                       reconstructed_variables_dict):
                """
                This function reconstructs necessary subclasses of a given
                class to be reconstructed. This is done by first reconstructing
                their variables via the function reconstruct_variables and
                afterwards parsing the structure for other classes which are
                also reconstructed by using the sub function
                reconstructRecursively. Modified structure dict, annotations
                and all other things are to be used to construct a
                ClassWithOptimizableVariables.
                """

                def reconstructRecursively(variable, subclasses_dict,
                                           reconstructed_variables_dict):
                    if self.isUUID(variable) and variable in subclasses_dict:
                        return reconstruct_class(
                                subclasses_dict[variable],
                                subclasses_dict,
                                reconstructed_variables_dict)
                    elif isinstance(variable, list):
                        return [reconstructRecursively(part,
                                                  subclasses_dict,
                                                  reconstructed_variables_dict)
                                for part in variable]
                    elif isinstance(variable, dict):
                        return dict([(key,
                                      reconstructRecursively(part,
                                                        subclasses_dict,
                                                        reconstructed_variables_dict))
                                    for (key, part) in variable.items()])
                    else:
                        return variable

                new_structure_dict = reconstructRecursively(structure_dict,
                                                            subclasses_dict,
                                                            reconstructed_variables_dict)
                return new_structure_dict

            self.debug("Reconstructing structure dictionary")
            structure_dict = class_to_be_reconstructed["structure"]
            self.debug(class_to_be_reconstructed["name"])
            self.debug(pformat(structure_dict))
            self.debug("Reconstructing optimizable variables")
            structure_dict = reconstruct_variables(structure_dict,
                                                   reconstructed_variables_dict)
            self.debug("Reconstructing sub classes")
            # self.debug(pformat(structure_dict))
            # self.debug(pformat(subclasses_dict))
            # self.debug(pformat(reconstructed_variables_dict))
            structure_dict = reconstruct_subclasses(structure_dict,
                                                    subclasses_dict,
                                                    reconstructed_variables_dict)
            self.debug("No UUIDs in structure dict left? (Shouldn\'t be after reconstruction!)" +
                  str(is_structure_free_of_uuids(structure_dict)))

            self.debug("Generating final object (constructor)")
            kind = class_to_be_reconstructed["kind"]
            name = class_to_be_reconstructed["name"]
            anno = class_to_be_reconstructed["annotations"]
            print("Reconstructing " + kind + " " + name)

            self.debug("Name: " + name)
            self.debug("Kind: " + kind)

            self.debug("Annotations")

            self.debug(pformat(anno))

            self.debug("Creating attributes")

            self.debug(pformat(structure_dict))

            myclass = self.classes_dictionary[kind](
                    anno,
                    structure_dict, name=name)
            return myclass

        """
        Reconstruct the variables pool by its own reconstruction functions.
        Then reconstruct the class to be reconstructed by calling
        reconstruct_class in a recursive manner. Return the final object.
        """

        optimizable_variables_pool = OptimizableVariablesPool.fromDictionary(
                optimizable_variables_pool_dict,
                functionobjects_pool_dict, source_checked, variables_checked)

        mynewobject = reconstruct_class(class_to_be_deserialized,
                                        subclasses_dict,
                                        optimizable_variables_pool.variables_pool)

        self.class_instance = mynewobject


