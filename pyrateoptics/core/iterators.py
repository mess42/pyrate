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

import numpy as np

from .log import BaseLogger
from .base import ClassWithOptimizableVariables
from .optimizable_variable import OptimizableVariable


class AbstractIterator(BaseLogger):
    """
    This class traverses through a class by investigating the __dict__
    and collects several information of the underlying recursive structure.
    """
    def __init__(self, class_instance, run=True, *args, **kwargs):
        self.class_instance = class_instance
        self.sub_instance = class_instance
        self.initVariables()
        if run:
            self.run(*args, **kwargs)

    def initVariables(self):
        raise NotImplementedError()

    def isCollectableElement(self, variable, *args, **kwargs):
        raise NotImplementedError()

    def isTraversableElement(self, variable, *args, **kwargs):
        raise NotImplementedError()

    def collectElement(self, variable, keystring, *args, **kwargs):
        raise NotImplementedError()

    def collectRest(self, variable, keystring, *args, **kwargs):
        pass

    def collectAccessSpecifier(self, var, accessspecifier_string,
                               *args, **kwargs):
        return ""

    def collectParentChild(self, typparent, parent, child, *args, **kwargs):
        pass

    def getCollectableVariableIdentifier(self, variable):
        raise NotImplementedError()

    def getSubInstanceVariableIdentifier(self, variable):
        raise NotImplementedError

    def isSubInstance(self, variable, *args, **kwargs):
        raise NotImplementedError()

    def preRun(self, *args, **kwargs):
        self.sub_instance = self.class_instance
        self.initVariables()

    def postRun(self, *args, **kwargs):
        pass

    def traverse(self, variable, keystring,
                 *args, **kwargs):
        if self.isTraversableElement(variable, *args, **kwargs):
            if self.isCollectableElement(variable, *args, **kwargs):
                newkeystring = keystring + self.collectAccessSpecifier(
                    variable,
                    self.getCollectableVariableIdentifier(variable),
                    *args, **kwargs)
                self.collectParentChild(type(variable), variable, None)
                self.collectElement(variable, newkeystring, *args, **kwargs)
                return
            elif isinstance(variable, list)\
                or isinstance(variable, tuple)\
                    or isinstance(variable, set):
                for (ind, part) in enumerate(variable):
                    newkeystring = keystring +\
                        self.collectAccessSpecifier(variable, str(ind),
                                                    *args, **kwargs)
                    self.collectParentChild(type(variable), variable, part)
                    self.traverse(part, newkeystring, *args, **kwargs)
            elif isinstance(variable, dict):
                for (key, value) in variable.items():
                    newkeystring = keystring + self.collectAccessSpecifier(
                        variable, key, *args, **kwargs)
                    self.collectParentChild(type(variable), variable, value)
                    self.traverse(value, newkeystring, *args, **kwargs)
            elif self.isSubInstance(variable):
                sub_instance_backup = self.sub_instance
                self.sub_instance = variable
                newkeystring = keystring +\
                    self.collectAccessSpecifier(
                        variable,
                        self.getSubInstanceVariableIdentifier(
                            variable
                        ),
                        *args, **kwargs)
                self.collectParentChild(type(variable), variable,
                                        variable.__dict__)
                self.traverse(variable.__dict__, newkeystring, *args, **kwargs)
                self.sub_instance = sub_instance_backup
            else:
                # all others
                self.collectRest(variable, keystring, *args, **kwargs)

    def run(self, *args, **kwargs):
        self.preRun(*args, **kwargs)
        self.traverse(self.class_instance, "", *args, **kwargs)
        self.postRun(*args, **kwargs)


class VariableReference(object):
    def __init__(self, value):
        self.value = value


class AbstractModifyingIterator(AbstractIterator):

    def modifyElement(self, variable_reference, newkeystring, *args, **kwargs):
        raise NotImplementedError()

    def traverse_modify(self, variable_reference, keystring, *args, **kwargs):
        value = variable_reference.value
        if self.isTraversableElement(value, *args, **kwargs):
            if self.isCollectableElement(value, *args, **kwargs):
                newkeystring = keystring + self.collectAccessSpecifier(
                    value,
                    self.getCollectableVariableIdentifier(value),
                    *args, **kwargs)
                self.collectParentChild(type(value), value, None)
                self.modifyElement(variable_reference, newkeystring,
                                   *args, **kwargs)
                return
            elif isinstance(value, list)\
                or isinstance(value, tuple)\
                    or isinstance(value, set):
                new_list = []
                for (ind, part) in enumerate(value):
                    newkeystring = keystring +\
                        self.collectAccessSpecifier(value, str(ind),
                                                    *args, **kwargs)
                    self.collectParentChild(type(value), value, part)
                    part_reference = VariableReference(part)
                    self.traverse_modify(part_reference, newkeystring,
                                         *args, **kwargs)
                    new_list.append(part_reference.value)
                variable_reference.value = type(value)(new_list)
            elif isinstance(value, dict):
                new_dict = {}
                for (key, dict_val) in value.items():
                    newkeystring = keystring + self.collectAccessSpecifier(
                        value, key, *args, **kwargs)
                    self.collectParentChild(type(value), value, dict_val)
                    dict_val_reference = VariableReference(dict_val)
                    self.traverse_modify(dict_val_reference, newkeystring,
                                         *args, **kwargs)
                    new_dict[key] = dict_val_reference.value
                variable_reference.value = new_dict
            elif self.isSubInstance(value):
                sub_instance_backup = self.sub_instance
                self.sub_instance = value
                newkeystring = keystring +\
                    self.collectAccessSpecifier(
                        value,
                        self.getSubInstanceVariableIdentifier(
                            value
                        ),
                        *args, **kwargs)
                self.collectParentChild(type(value), value,
                                        value.__dict__)
                class_dict_reference = VariableReference(value.__dict__)
                self.traverse_modify(class_dict_reference,
                                     newkeystring, *args, **kwargs)
                value.__dict__ = class_dict_reference.value
                self.sub_instance = sub_instance_backup
            else:
                # all others
                self.collectRest(value, keystring, *args, **kwargs)

    def run(self, *args, **kwargs):
        self.preRun(*args, **kwargs)
        class_instance_reference = VariableReference(self.class_instance)
        self.traverse_modify(class_instance_reference, "", *args, **kwargs)
        self.class_instance = class_instance_reference.value
        self.postRun(*args, **kwargs)


class OptimizableVariableIterator(AbstractIterator):
    """
    Traverse ClassWithOptimizableVariables and their subclasses
    and substructures. It does not work for OptimizableVariables which
    are arguments to a pickup and were defined outside of the actual
    ClassWithOptimizableVariables. This is done in the
    OptimizableVariablesPool.

    Accumulates optimizable variables and its linked objects.
    Ignores ring-links and double links.
    """

    def initVariables(self):
        self.idlist = []

    def isTraversableElement(self, variable, *args, **kwargs):
        variable_id = id(variable)
        if variable_id not in self.idlist:
            self.idlist.append(variable_id)
            return True
        else:
            return False

    def getCollectableVariableIdentifier(self, variable):
        return variable.name

    def getSubInstanceVariableIdentifier(self, variable):
        return variable.name

    def isCollectableElement(self, variable, *args, **kwargs):
        return isinstance(variable, OptimizableVariable)

    def isSubInstance(self, variable):
        return isinstance(variable, ClassWithOptimizableVariables)


class AbstractOptimizableVariableModifyingCollector(AbstractModifyingIterator,
                                                    OptimizableVariableIterator):
    """
    AbstractOptimizableVariableModifyingCollector takes the description
    of special objects (isSubInstance, isCollectableElement) from
    OptimizableVariableIterator, but the run method uses traverse_modify
    from AbstractModifyingIterator.

    The method modifyElement still has to be implemented by subclass.
    """
    def run(self, *args, **kwargs):
        AbstractModifyingIterator.run(*args, **kwargs)


class OptimizableVariableCollector(OptimizableVariableIterator):
    """
    For fast evaluation of value vector. Collects variables in list.
    Provides interface for transformation into numpy arrays and back,
    even with transformation.
    """

    def initVariables(self):
        super(OptimizableVariableCollector, self).initVariables()
        self.variables_list = []

    def collectElement(self, variable, keystring, *args, **kwargs):
        self.variables_list.append(variable)

    def toNumpyArray(self):
        """
        Function to get all values into one large np.array.
        Supports only float at the moment.
        """
        return np.fromiter([v() for v in self.variables_list],
                           dtype=float, count=len(self.variables_list))

    def toNumpyArrayTransformed(self):
        """
        Function to get all transformed values into one large np.array.
        Supports only float at the moment.
        """
        return np.fromiter([v.evaluate_transformed()
                            for v in self.variables_list],
                           dtype=float, count=len(self.variables_list))

    def fromNumpyArray(self, x):
        """
        Function to set all values of active variables to the values in the
        large np.array x. Supports only float at the moment.
        """
        [variable.set_value(value)
         for (variable, value) in zip(self.variables_list, x.tolist())]

    def fromNumpyArrayTransformed(self, x):
        """
        Function to set all values of active variables to the transformed
        values in the large np.array x. Supports only float at the moment.
        """
        [variable.set_value_transformed(value)
         for (variable, value) in zip(self.variables_list, x.tolist())]


class OptimizableVariableKeyIterator(OptimizableVariableCollector):

    """
    Creates a dictionary where the keys are easy to memoize names of
    the which can be used to refer to the variables.
    """

    def initVariables(self):
        super(OptimizableVariableKeyIterator, self).initVariables()
        self.variables_dictionary = {}

    def collectAccessSpecifier(self, variable,
                               accessspecifier, *args, **kwargs):
        (shortkeys,) = args

        if isinstance(variable, dict) or isinstance(variable, list):
            return "" if shortkeys else "[\"" + accessspecifier + "\"]"
        elif self.isSubInstance(variable):
            return accessspecifier + "."\
                    if shortkeys else "(" + accessspecifier + ")."
        elif self.isCollectableElement(variable):
            return accessspecifier\
                    if shortkeys else "(" + accessspecifier + ")"
        else:
            return ""

    def collectElement(self, variable, keystring, *args, **kwargs):
        self.variables_dictionary[keystring] = variable

    def collectRest(self, variable, keystring, *args, **kwargs):
        pass

    def run(self, shortkeys=True):
        super(OptimizableVariableKeyIterator, self).run(shortkeys)


class OptimizableVariableSetKeyIterator(
        OptimizableVariableKeyIterator,
        AbstractOptimizableVariableModifyingCollector):

    def __init__(self, class_instance, key_assignment_dictionary, run=True,
                 *args, **kwargs):
        super().__init__(
            class_instance,
            run=run, *args, **kwargs)
        self.key_assignment_dictionary = key_assignment_dictionary

    def run(self, shortkeys=True, *args, **kwargs):
        AbstractOptimizableVariableModifyingCollector.run(self, shortkeys)

    def modifyElement(self, variable_reference, newkeystring, *args, **kwargs):
        if self.key_assignment_dictionary.get(newkeystring, None) is not None:
            variable_reference.value =\
                self.key_assignment_dictionary[newkeystring]


class OptimizableVariableActiveCollector(OptimizableVariableCollector):
    """
    Returns a list of active variables the names are lost
    but it does not matter since the variable references are still in the
    original dictionary in the class.
    """

    def isCollectableElement(self, variable, *args, **kwargs):
        return super(OptimizableVariableCollector,
                     self).isCollectableElement(variable,
                                                *args,
                                                **kwargs) and\
                variable.var_type() == "variable"


class SerializationIterator(OptimizableVariableCollector):
    """
    Generates a list which gives enough information to reconstruct a
    given ClassWithOptimizableVariables.
    """

    def initVariables(self):
        super(SerializationIterator, self).initVariables()
        self.classes_dictionary = {}
        self.variables_dictionary = {}
        self.structure = None

    def collectStructure(self, remove=[]):

        def traverse_structure(var, remove=[]):
            # this code is in some sense doubled traverse code
            # but its intention is another: while traverse only
            # traverses, traverse_structure modifies a data object with several
            # return statements
            if isinstance(var, dict):
                newdict = {}
                for (k, v) in var.items():
                    newitem = traverse_structure(v, remove=remove)
                    if k not in remove and newitem is not None:
                        newdict[k] = newitem
                return newdict
            elif isinstance(var, list):
                newlist = []
                for (ind, v) in enumerate(var):
                    newitem = traverse_structure(v, remove=remove)
                    if newitem is not None:
                        newlist.append(newitem)
                return newlist
            elif self.isCollectableElement(var):
                return var.unique_id
            elif self.isSubInstance(var):
                return var.unique_id
            else:
                # ignore all other variables
                return None

        #mymodifieddict = {}
        #for (k, v) in self.class_instance.__dict__.items():
        #    if k not in self.class_instance.serializationfilter:
        #        mymodifieddict[k] = v
        self.structure = traverse_structure(self.class_instance.__dict__,
                                            remove=remove)

        # build structure via getTypeFromDict in classwithoptvars
        # structure should contain unique ids which identify classes
        # and variables
        # This structure should be created by the constructors such
        # that it gives a valid empty object
        # all other stuff should be in annotations

    def isSubInstance(self, variable):
        result = super(SerializationIterator, self).isSubInstance(variable)
        if result and variable.unique_id not in self.classes_dictionary:
            self.classes_dictionary[variable.unique_id] =\
                variable
        return result

    def postRun(self, remove=[]):
        self.variables_dictionary = dict([(v.unique_id, v)
                                          for v in self.variables_list])
        self.classes_dictionary.pop(self.class_instance.unique_id)
        # remove class_instance from classes_dictionary
        self.collectStructure(remove=remove)
        self.dictionary = {}
        self.dictionary.update(self.class_instance.getBasicInfo())
        self.dictionary["annotations"] = self.class_instance.annotations
        self.dictionary["structure"] = self.structure
        # standard things kind, name, version, id can be given by
        # function from class, structure and rest comes from iterator
